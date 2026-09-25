#include "gh_prs.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Round half-to-even at `digits` decimals, matching Python's round().
 *
 * Formatting to a decimal string and parsing back reproduces Python's
 * behaviour: both round the exact binary value correctly to that many
 * decimals, then take the nearest double. Doing the arithmetic with
 * pow(10, digits) instead would disagree on ties. */
double gh_round(double value, int digits)
{
    if (!isfinite(value))
        return value;

    char buf[64];
    snprintf(buf, sizeof(buf), "%.*f", digits, value);
    return strtod(buf, NULL);
}

double gh_prs_z_to_percentile(double z)
{
    double raw = 0.5 * (1.0 + erf(z / sqrt(2.0))) * 100.0;
    if (raw < 0.1)
        return 0.1;
    if (raw > 99.9)
        return 99.9;
    return raw;
}

const char *gh_prs_categorize(double percentile)
{
    if (percentile < 20.0)
        return "low";
    if (percentile < 80.0)
        return "average";
    if (percentile < 95.0)
        return "elevated";
    return "high";
}

static size_t count_allele(const char *genotype, char allele)
{
    size_t n = 0;
    for (const char *p = genotype; *p; p++)
        if (*p == allele)
            n++;
    return n;
}

/* Contributing SNPs sort by descending contribution, ties keeping model
 * order -- Python sorts on -contribution with a stable sort. */
typedef struct {
    gh_prs_contribution item;
    size_t order;
} ranked_contribution;

static int compare_contribution(const void *a, const void *b)
{
    const ranked_contribution *x = a, *y = b;
    if (x->item.contribution != y->item.contribution)
        return x->item.contribution > y->item.contribution ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

/* Format a percentage the way Python's "{:.0%}" does. */
static void format_percent(char *buf, size_t size, double fraction)
{
    snprintf(buf, size, "%.0f%%", fraction * 100.0);
}

static void score_model(gh_prs_result *out, gh_arena *arena,
                        const gh_prs_model *model, const gh_genome *g,
                        double eur_prop)
{
    memset(out, 0, sizeof(*out));
    out->id = model->id;
    out->name = model->name;
    out->reference = model->reference;
    out->ancestry_applicable = true;
    out->ancestry_warning = "";

    double raw_score = 0.0, pop_mean = 0.0, pop_var = 0.0;
    size_t snps_found = 0, unique_total = 0;

    ranked_contribution *ranked =
        gh_calloc(arena, model->nsnps, sizeof(*ranked));
    size_t ncontrib = 0;

    /* Models may list an rsID more than once; only the first occurrence
     * counts, as in the Python's seen-set. */
    const char **seen = gh_calloc(arena, model->nsnps, sizeof(*seen));
    size_t nseen = 0;

    for (size_t i = 0; i < model->nsnps; i++) {
        const gh_prs_snp *snp = &model->snps[i];

        bool duplicate = false;
        for (size_t s = 0; s < nseen && !duplicate; s++)
            if (strcmp(seen[s], snp->rsid) == 0)
                duplicate = true;
        if (duplicate)
            continue;
        seen[nseen++] = snp->rsid;
        unique_total++;

        const char *genotype = gh_genome_genotype(g, snp->rsid);
        if (!genotype)
            continue;

        size_t copies = count_allele(genotype, snp->risk_allele);
        double contribution = snp->log_or * (double)copies;
        raw_score += contribution;

        /* Expected mean and variance under the model's European frequencies. */
        pop_mean += snp->log_or * 2.0 * snp->eur_freq;
        pop_var += snp->log_or * snp->log_or * 2.0 * snp->eur_freq
                 * (1.0 - snp->eur_freq);

        snps_found++;
        if (copies > 0) {
            ranked[ncontrib].item.rsid = snp->rsid;
            ranked[ncontrib].item.gene = snp->gene;
            ranked[ncontrib].item.risk_allele = snp->risk_allele;
            ranked[ncontrib].item.copies = (int)copies;
            ranked[ncontrib].item.log_or = snp->log_or;
            ranked[ncontrib].item.contribution = contribution;
            ranked[ncontrib].order = ncontrib;
            ncontrib++;
        }
    }

    double z_score, percentile, ci_lower, ci_upper;
    if (snps_found > 0 && pop_var > 0.0) {
        double sd = sqrt(pop_var);
        z_score = (raw_score - pop_mean) / sd;
        /* Standard error assuming independent SNPs, expressed in z units. */
        double se = snps_found > 1 ? sd / sqrt((double)snps_found) : sd;
        double z_se = se / sd;
        percentile = gh_prs_z_to_percentile(z_score);
        ci_lower = gh_prs_z_to_percentile(z_score - 1.96 * z_se);
        ci_upper = gh_prs_z_to_percentile(z_score + 1.96 * z_se);
    } else {
        z_score = 0.0;
        percentile = 50.0;
        ci_lower = 50.0;
        ci_upper = 50.0;
    }

    const char *category = gh_prs_categorize(percentile);

    /* These models are calibrated on European cohorts, so the z-score is
     * scaled by estimated European ancestry and flagged when it drops far. */
    if (eur_prop < 0.95) {
        double adjusted = gh_prs_z_to_percentile(z_score * eur_prop);

        char pct[32], unadjusted[32];
        format_percent(pct, sizeof(pct), eur_prop);
        snprintf(unadjusted, sizeof(unadjusted), "%.0f", percentile);

        char *warning = gh_alloc(arena, 256);
        if (eur_prop < 0.5) {
            out->ancestry_applicable = false;
            snprintf(warning, 256,
                "PRS calibrated on European populations. Your EUR ancestry "
                "is ~%s, so these scores have reduced accuracy. "
                "Unadjusted percentile: %sth.", pct, unadjusted);
        } else {
            snprintf(warning, 256,
                "PRS adjusted for %s European ancestry "
                "(unadjusted: %sth percentile).", pct, unadjusted);
        }
        out->ancestry_warning = warning;

        percentile = gh_round(adjusted, 1);
        category = gh_prs_categorize(percentile);
    }

    qsort(ranked, ncontrib, sizeof(*ranked), compare_contribution);

    out->raw_score = gh_round(raw_score, 4);
    out->z_score = gh_round(z_score, 3);
    out->percentile = gh_round(percentile, 1);
    out->ci_95_lower = gh_round(ci_lower, 1);
    out->ci_95_upper = gh_round(ci_upper, 1);
    out->risk_category = category;
    out->snps_found = snps_found;
    out->snps_total = unique_total;

    out->ncontributing = ncontrib < GH_PRS_MAX_CONTRIBUTIONS
        ? ncontrib : GH_PRS_MAX_CONTRIBUTIONS;
    for (size_t i = 0; i < out->ncontributing; i++)
        out->contributing[i] = ranked[i].item;
}

void gh_calculate_prs(gh_prs_result *out, gh_arena *arena, const gh_genome *g,
                      const gh_ancestry_result *ancestry)
{
    double eur_prop = 1.0;
    if (ancestry) {
        for (size_t p = 0; p < GH_POPULATION_COUNT; p++)
            if (strcmp(GH_POPULATIONS[p], "EUR") == 0) {
                eur_prop = ancestry->proportions[p];
                break;
            }
    }

    for (size_t i = 0; i < GH_PRS_MODEL_COUNT; i++)
        score_model(&out[i], arena, &GH_PRS_MODELS[i], g, eur_prop);
}

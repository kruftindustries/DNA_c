/* Unit tests for ancestry estimation and polygenic risk scores. */

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gh_ancestry.h"
#include "gh_genome.h"
#include "gh_mem.h"
#include "gh_prs.h"

static int failures;
static int checks;

static void check(int cond, const char *what)
{
    checks++;
    if (!cond) {
        failures++;
        fprintf(stderr, "  FAIL: %s\n", what);
    }
}

static void check_str(const char *got, const char *want, const char *what)
{
    checks++;
    if (!got || strcmp(got, want) != 0) {
        failures++;
        fprintf(stderr, "  FAIL: %s\n    got  %s\n    want %s\n",
                what, got ? got : "(null)", want);
    }
}

static void check_close(double got, double want, double tol, const char *what)
{
    checks++;
    if (!(fabs(got - want) <= tol)) {
        failures++;
        fprintf(stderr, "  FAIL: %s\n    got  %.10g\n    want %.10g\n",
                what, got, want);
    }
}

/* ------------------------------------------------------------------ */
/* Fixtures                                                            */
/* ------------------------------------------------------------------ */

static void write_rows(FILE *f, int *line, const char *rsid,
                       const char *genotype)
{
    fprintf(f, "%s\t1\t%d\t%s\n", rsid, 1000 + (*line)++, genotype);
}

static void load_fixture(gh_genome *g, gh_arena *arena, ...)
{
    static int counter;
    char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_anc_fixture_%d.txt", counter++);

    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "  FAIL: cannot write %s\n", path);
        failures++;
        return;
    }
    fputs("# fixture\n", f);

    va_list ap;
    va_start(ap, arena);
    int line = 0;
    for (;;) {
        const char *rsid = va_arg(ap, const char *);
        if (!rsid)
            break;
        write_rows(f, &line, rsid, va_arg(ap, const char *));
    }
    va_end(ap);
    fclose(f);

    if (!gh_genome_load(g, arena, path, NULL)) {
        fprintf(stderr, "  FAIL: cannot load %s\n", path);
        failures++;
    }
}

/* Build a genome from every AIM, each set to the dosage that the named
 * population most (or least) expects. */
static void load_ancestry_fixture(gh_genome *g, gh_arena *arena,
                                  size_t population, bool favour)
{
    static int counter;
    char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_anc_pop_%d.txt", counter++);

    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "  FAIL: cannot write %s\n", path);
        failures++;
        return;
    }
    fputs("# fixture\n", f);

    int line = 0;
    for (size_t i = 0; i < GH_AIM_COUNT; i++) {
        const gh_aim *aim = &GH_AIMS[i];
        char allele = aim->allele;
        char other = allele == 'A' ? 'C' : 'A';
        bool high = aim->frequencies[population] >= 0.5;
        bool want = favour ? high : !high;
        char genotype[3] = {want ? allele : other, want ? allele : other, '\0'};
        write_rows(f, &line, aim->rsid, genotype);
    }
    fclose(f);

    if (!gh_genome_load(g, arena, path, NULL)) {
        fprintf(stderr, "  FAIL: cannot load %s\n", path);
        failures++;
    }
}

/* ------------------------------------------------------------------ */
/* Ancestry                                                            */
/* ------------------------------------------------------------------ */

static void test_ancestry_no_markers(void)
{
    gh_arena *a = gh_arena_new(1 << 18);
    gh_genome g;
    gh_ancestry_result r;

    load_fixture(&g, a, "rs00000000", "AA", NULL);
    gh_estimate_ancestry(&r, a, &g);

    check(r.markers_found == 0, "no AIMs found");
    check_str(r.confidence, "none", "confidence is none");
    check_str(r.top_ancestry, "Unknown", "top ancestry is unknown");
    check(r.top_index == -1, "no top index");
    check(r.ndetails == 0, "no details");
    check(!r.sub.present, "no sub-ancestry");

    /* With nothing to go on, the prior is uniform and still sums to one. */
    double total = 0.0;
    for (size_t p = 0; p < GH_POPULATION_COUNT; p++) {
        check_close(r.proportions[p], 1.0 / (double)GH_POPULATION_COUNT, 1e-12,
                    "uniform prior");
        total += r.proportions[p];
    }
    check_close(total, 1.0, 1e-12, "proportions sum to one");

    gh_arena_free(a);
}

static void test_ancestry_recovers_each_population(void)
{
    gh_arena *a = gh_arena_new(1 << 20);
    gh_genome g;
    gh_ancestry_result r;

    /* A genome built to match a population's expected dosage should be
     * assigned to that population -- except AMR, which is admixed: its
     * allele frequencies sit between the others across this panel, so
     * dosage matching resolves it to SAS. The Python does the same. */
    for (size_t p = 0; p < GH_POPULATION_COUNT; p++) {
        bool identifiable = strcmp(GH_POPULATIONS[p], "AMR") != 0;

        load_ancestry_fixture(&g, a, p, true);
        gh_estimate_ancestry(&r, a, &g);

        check(r.markers_found == GH_AIM_COUNT, "every AIM found");
        if (identifiable) {
            check(r.top_index == (int)p, "top population recovered");
            check_str(r.top_ancestry, GH_POPULATION_LABELS[p], "label matches");
            check(r.proportions[p] > 0.5, "winning proportion dominates");
        } else {
            check(r.top_index >= 0, "an admixed genome still gets a call");
            check(r.proportions[r.top_index] > 0.5,
                  "the call it does make is confident");
        }

        double total = 0.0;
        for (size_t q = 0; q < GH_POPULATION_COUNT; q++) {
            check(r.proportions[q] >= 0.0 && r.proportions[q] <= 1.0,
                  "proportion is a probability");
            total += r.proportions[q];
        }
        check_close(total, 1.0, 1e-9, "proportions sum to one");
    }

    /* All 54 markers, highly discriminating, is a high-confidence call. */
    load_ancestry_fixture(&g, a, 0, true);
    gh_estimate_ancestry(&r, a, &g);
    check_str(r.confidence, "high", "full marker set is high confidence");
    check(r.ndetails == GH_AIM_COUNT, "a detail per marker");
    check(r.details[0].rsid != NULL, "detail carries its rsID");
    check(r.details[0].allele_count >= 0 && r.details[0].allele_count <= 2,
          "allele count is 0-2");

    gh_arena_free(a);
}

static void test_ancestry_confidence_bands(void)
{
    gh_arena *a = gh_arena_new(1 << 18);
    gh_genome g;
    gh_ancestry_result r;

    /* A handful of markers cannot support a confident call. */
    load_fixture(&g, a, GH_AIMS[0].rsid, "AA", GH_AIMS[1].rsid, "AA", NULL);
    gh_estimate_ancestry(&r, a, &g);
    check(r.markers_found == 2, "two markers found");
    check_str(r.confidence, "low", "two markers is low confidence");
    check(r.top_index >= 0, "a top population is still chosen");

    gh_arena_free(a);
}

static void test_sub_ancestry(void)
{
    gh_arena *a = gh_arena_new(1 << 20);
    gh_genome g;
    gh_ancestry_result r;

    /* EUR has sub-population markers; a EUR-leaning genome should get one. */
    size_t eur = 0;
    for (size_t p = 0; p < GH_POPULATION_COUNT; p++)
        if (strcmp(GH_POPULATIONS[p], "EUR") == 0)
            eur = p;

    load_ancestry_fixture(&g, a, eur, true);
    gh_estimate_ancestry(&r, a, &g);
    check(r.top_index == (int)eur, "EUR genome resolves to EUR");
    check(r.sub.present, "EUR has a sub-ancestry estimate");
    check(r.sub.nsubs == 3, "three European sub-populations");
    check(r.sub.top_label != NULL, "a top sub-population is named");
    check(r.sub.markers_used >= 3, "at least three sub-markers used");

    double total = 0.0;
    for (size_t i = 0; i < r.sub.nsubs; i++) {
        check(r.sub.labels[i] != NULL, "sub-population labelled");
        total += r.sub.proportions[i];
    }
    check_close(total, 1.0, 1e-9, "sub-proportions sum to one");

    /* Fewer than three sub-markers means no estimate rather than a guess. */
    load_fixture(&g, a, "rs12913832", "GG", "rs16891982", "GG", NULL);
    gh_estimate_ancestry(&r, a, &g);
    check(!r.sub.present, "two sub-markers is not enough");

    gh_arena_free(a);
}

/* ------------------------------------------------------------------ */
/* PRS helpers                                                         */
/* ------------------------------------------------------------------ */

static void test_round_half_even(void)
{
    /* Python's round() breaks ties to even, which ordinary
     * floor(x*10^n + 0.5) arithmetic would get wrong. */
    check_close(gh_round(0.5, 0), 0.0, 0.0, "0.5 rounds to 0");
    check_close(gh_round(1.5, 0), 2.0, 0.0, "1.5 rounds to 2");
    check_close(gh_round(2.5, 0), 2.0, 0.0, "2.5 rounds to 2");
    check_close(gh_round(-0.5, 0), 0.0, 0.0, "-0.5 rounds to 0");
    check_close(gh_round(-1.5, 0), -2.0, 0.0, "-1.5 rounds to -2");
    check_close(gh_round(1.25, 1), 1.2, 1e-12, "1.25 rounds to 1.2");
    check_close(gh_round(1.35, 1), 1.4, 1e-12, "1.35 rounds to 1.4");
    check_close(gh_round(3.14159, 3), 3.142, 1e-12, "pi to three places");
    check_close(gh_round(-3.14159, 2), -3.14, 1e-12, "negative rounding");
    check_close(gh_round(0.0, 4), 0.0, 0.0, "zero");
}

static void test_z_to_percentile(void)
{
    check_close(gh_prs_z_to_percentile(0.0), 50.0, 1e-9, "z=0 is the median");
    check_close(gh_prs_z_to_percentile(1.0), 84.134, 0.01, "z=1");
    check_close(gh_prs_z_to_percentile(-1.0), 15.866, 0.01, "z=-1");
    check_close(gh_prs_z_to_percentile(1.96), 97.5, 0.01, "z=1.96");

    /* Extremes clamp rather than saturating at 0 or 100. */
    check_close(gh_prs_z_to_percentile(10.0), 99.9, 1e-9, "clamped above");
    check_close(gh_prs_z_to_percentile(-10.0), 0.1, 1e-9, "clamped below");
}

static void test_categorize(void)
{
    check_str(gh_prs_categorize(0.0), "low", "0th percentile is low");
    check_str(gh_prs_categorize(19.9), "low", "just under 20 is low");
    check_str(gh_prs_categorize(20.0), "average", "20 is average");
    check_str(gh_prs_categorize(79.9), "average", "just under 80 is average");
    check_str(gh_prs_categorize(80.0), "elevated", "80 is elevated");
    check_str(gh_prs_categorize(94.9), "elevated", "just under 95 is elevated");
    check_str(gh_prs_categorize(95.0), "high", "95 is high");
    check_str(gh_prs_categorize(99.9), "high", "the top is high");
}

/* ------------------------------------------------------------------ */
/* PRS scoring                                                         */
/* ------------------------------------------------------------------ */

static const gh_prs_result *find_model(const gh_prs_result *results,
                                       const char *id)
{
    for (size_t i = 0; i < GH_PRS_MODEL_COUNT; i++)
        if (strcmp(results[i].id, id) == 0)
            return &results[i];
    return NULL;
}

/* Write a genome carrying `copies` of the risk allele for every SNP of the
 * named model. */
static void load_prs_fixture(gh_genome *g, gh_arena *arena, const char *id,
                             int copies)
{
    static int counter;
    char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_prs_fixture_%d.txt", counter++);

    const gh_prs_model *model = NULL;
    for (size_t i = 0; i < GH_PRS_MODEL_COUNT; i++)
        if (strcmp(GH_PRS_MODELS[i].id, id) == 0)
            model = &GH_PRS_MODELS[i];
    if (!model) {
        fprintf(stderr, "  FAIL: no PRS model %s\n", id);
        failures++;
        return;
    }

    FILE *f = fopen(path, "wb");
    if (!f) {
        failures++;
        return;
    }
    fputs("# fixture\n", f);
    int line = 0;
    for (size_t i = 0; i < model->nsnps; i++) {
        char risk = model->snps[i].risk_allele;
        char other = risk == 'A' ? 'C' : 'A';
        char genotype[3];
        genotype[0] = copies >= 1 ? risk : other;
        genotype[1] = copies >= 2 ? risk : other;
        genotype[2] = '\0';
        write_rows(f, &line, model->snps[i].rsid, genotype);
    }
    fclose(f);

    if (!gh_genome_load(g, arena, path, NULL)) {
        failures++;
    }
}

static void test_prs_no_data(void)
{
    gh_arena *a = gh_arena_new(1 << 20);
    gh_genome g;
    gh_prs_result *results = gh_calloc(a, GH_PRS_MODEL_COUNT, sizeof(*results));

    load_fixture(&g, a, "rs00000000", "AA", NULL);
    gh_calculate_prs(results, a, &g, NULL);

    for (size_t i = 0; i < GH_PRS_MODEL_COUNT; i++) {
        const gh_prs_result *r = &results[i];
        check(r->snps_found == 0, "no SNPs found");
        check(r->snps_total > 0, "model still reports its size");
        check_close(r->percentile, 50.0, 1e-9, "no data means median");
        check_close(r->z_score, 0.0, 1e-9, "no data means z of zero");
        check_str(r->risk_category, "average", "median is average risk");
        check(r->ncontributing == 0, "nothing contributed");
        check(r->name != NULL && r->reference != NULL, "metadata present");
    }

    gh_arena_free(a);
}

static void test_prs_scoring_direction(void)
{
    gh_arena *a = gh_arena_new(1 << 20);
    gh_genome g;
    gh_prs_result *results = gh_calloc(a, GH_PRS_MODEL_COUNT, sizeof(*results));

    /* Two copies of every risk allele must score above the population mean;
     * zero copies must score below it. */
    /* `results` is overwritten by each call, so keep copies rather than
     * pointers into it. */
    load_prs_fixture(&g, a, "type2_diabetes", 2);
    gh_calculate_prs(results, a, &g, NULL);
    gh_prs_result high = *find_model(results, "type2_diabetes");
    check(high.snps_found > 0, "SNPs found at full risk dosage");
    check(high.z_score > 0.0, "all risk alleles scores above the mean");
    check(high.percentile > 50.0, "and above the median percentile");
    check(high.ncontributing > 0, "contributing SNPs recorded");
    check(high.ncontributing <= GH_PRS_MAX_CONTRIBUTIONS,
          "at most ten contributors reported");

    /* Contributions must come back in descending order. */
    for (size_t i = 1; i < high.ncontributing; i++)
        check(high.contributing[i - 1].contribution
                  >= high.contributing[i].contribution,
              "contributions sorted descending");
    for (size_t i = 0; i < high.ncontributing; i++) {
        check(high.contributing[i].copies > 0, "only carriers contribute");
        check_close(high.contributing[i].contribution,
                    high.contributing[i].log_or
                        * (double)high.contributing[i].copies,
                    1e-12, "contribution is log-OR times copies");
    }

    load_prs_fixture(&g, a, "type2_diabetes", 0);
    gh_calculate_prs(results, a, &g, NULL);
    gh_prs_result low = *find_model(results, "type2_diabetes");
    check(low.z_score < 0.0, "no risk alleles scores below the mean");
    check(low.percentile < 50.0, "and below the median percentile");
    check(low.ncontributing == 0, "no contributors without risk alleles");
    check_close(low.raw_score, 0.0, 1e-12, "raw score is zero");

    /* One copy sits between the two. */
    load_prs_fixture(&g, a, "type2_diabetes", 1);
    gh_calculate_prs(results, a, &g, NULL);
    gh_prs_result mid = *find_model(results, "type2_diabetes");
    check(mid.z_score > low.z_score && mid.z_score < high.z_score,
          "one copy sits between none and two");

    /* The confidence interval brackets the point estimate. */
    check(high.ci_95_lower <= high.percentile, "CI lower bound");
    check(high.ci_95_upper >= high.percentile, "CI upper bound");

    gh_arena_free(a);
}

static void test_prs_ancestry_adjustment(void)
{
    gh_arena *a = gh_arena_new(1 << 20);
    gh_genome g;
    gh_prs_result *results = gh_calloc(a, GH_PRS_MODEL_COUNT, sizeof(*results));

    /* One risk copy per SNP, so the percentile lands mid-range. Two copies
     * clamps at 99.9, where scaling the z-score cannot move it and the
     * adjustment would be invisible. */
    load_prs_fixture(&g, a, "type2_diabetes", 1);

    /* No ancestry supplied: European assumed, no adjustment, no warning. */
    gh_calculate_prs(results, a, &g, NULL);
    gh_prs_result plain = *find_model(results, "type2_diabetes");
    double unadjusted = plain.percentile;
    check(unadjusted > 50.0 && unadjusted < 99.0,
          "fixture percentile is mid-range, not clamped");
    check(plain.ancestry_applicable, "applicable without ancestry data");
    check_str(plain.ancestry_warning, "", "no warning without ancestry data");

    /* Fully European: still no adjustment. */
    gh_ancestry_result ancestry;
    memset(&ancestry, 0, sizeof(ancestry));
    for (size_t p = 0; p < GH_POPULATION_COUNT; p++)
        ancestry.proportions[p] = strcmp(GH_POPULATIONS[p], "EUR") == 0 ? 1.0 : 0.0;
    gh_calculate_prs(results, a, &g, &ancestry);
    gh_prs_result eur = *find_model(results, "type2_diabetes");
    check_close(eur.percentile, unadjusted, 1e-9, "no adjustment at EUR 1.0");
    check_str(eur.ancestry_warning, "", "no warning at EUR 1.0");

    /* Partly European: the score is pulled toward the median and explained. */
    for (size_t p = 0; p < GH_POPULATION_COUNT; p++)
        ancestry.proportions[p] = strcmp(GH_POPULATIONS[p], "EUR") == 0 ? 0.7 : 0.075;
    gh_calculate_prs(results, a, &g, &ancestry);
    gh_prs_result mixed = *find_model(results, "type2_diabetes");
    check(mixed.ancestry_applicable, "still applicable at EUR 0.7");
    check(mixed.percentile < unadjusted, "adjusted toward the median");
    check(strstr(mixed.ancestry_warning, "adjusted for") != NULL,
          "warning explains the adjustment");
    check(strstr(mixed.ancestry_warning, "70%") != NULL,
          "warning states the ancestry proportion");

    /* Mostly non-European: flagged as not applicable. */
    for (size_t p = 0; p < GH_POPULATION_COUNT; p++)
        ancestry.proportions[p] = strcmp(GH_POPULATIONS[p], "EUR") == 0 ? 0.2 : 0.2;
    gh_calculate_prs(results, a, &g, &ancestry);
    gh_prs_result low_eur = *find_model(results, "type2_diabetes");
    check(!low_eur.ancestry_applicable, "not applicable below EUR 0.5");
    check(strstr(low_eur.ancestry_warning, "reduced accuracy") != NULL,
          "warning states reduced accuracy");
    check(strstr(low_eur.ancestry_warning, "Unadjusted percentile") != NULL,
          "warning reports the unadjusted value");

    gh_arena_free(a);
}

static void test_prs_model_integrity(void)
{
    /* The Python validates these at import time; the same must hold here. */
    for (size_t i = 0; i < GH_PRS_MODEL_COUNT; i++) {
        const gh_prs_model *m = &GH_PRS_MODELS[i];
        check(m->id && m->name && m->reference, "model metadata present");
        check(m->nsnps > 0, "model has SNPs");

        for (size_t j = 0; j < m->nsnps; j++) {
            const gh_prs_snp *snp = &m->snps[j];
            check(snp->log_or > -1.0 && snp->log_or < 1.0,
                  "effect size within the expected range");
            check(snp->eur_freq > 0.0 && snp->eur_freq < 1.0,
                  "frequency is a proper probability");
            check(strchr("ACGT", snp->risk_allele) != NULL,
                  "risk allele is a base");

            /* No duplicate rsIDs within a model. */
            for (size_t k = j + 1; k < m->nsnps; k++)
                check(strcmp(snp->rsid, m->snps[k].rsid) != 0,
                      "no duplicate rsID in a model");
        }
    }
}

int main(void)
{
    test_ancestry_no_markers();
    test_ancestry_recovers_each_population();
    test_ancestry_confidence_bands();
    test_sub_ancestry();
    test_round_half_even();
    test_z_to_percentile();
    test_categorize();
    test_prs_no_data();
    test_prs_scoring_direction();
    test_prs_ancestry_adjustment();
    test_prs_model_integrity();

    printf("%s: %d checks, %d failures\n",
           failures ? "FAILED" : "ok", checks, failures);
    return failures ? 1 : 0;
}

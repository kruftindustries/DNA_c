#include "gh_dependent_profiles.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gh_prs.h"   /* gh_round */

/* ------------------------------------------------------------------ */
/* Shared helpers                                                      */
/* ------------------------------------------------------------------ */

static size_t count_allele(const char *genotype, char allele)
{
    size_t n = 0;
    for (const char *p = genotype; p && *p; p++)
        if (*p == allele)
            n++;
    return n;
}

static bool list_has(const char *const *list, size_t n, const char *value)
{
    if (!value)
        return false;
    for (size_t i = 0; i < n; i++)
        if (strcmp(list[i], value) == 0)
            return true;
    return false;
}

/* The status these modules attribute to a gene.
 *
 * They build `gene_status[gene] = status` while walking the findings, so the
 * LAST finding for a gene wins -- unlike insights and recommendations, which
 * keep the strongest. Findings arrive sorted by descending magnitude, so this
 * is the weakest call, not the strongest. Reproduced deliberately. */
static const char *last_status(const gh_analysis *a, const char *gene)
{
    const char *status = NULL;
    for (size_t i = 0; a && i < a->nfindings; i++) {
        const gh_finding *f = &a->findings[i];
        if (f->gene && f->status && *f->gene && *f->status &&
            strcmp(f->gene, gene) == 0)
            status = f->status;
    }
    return status;
}

/* Python's round() on a positive value, half to even. */
static long round_half_even(double value)
{
    double floored = floor(value);
    double frac = value - floored;
    if (frac > 0.5)
        return (long)floored + 1;
    if (frac < 0.5)
        return (long)floored;
    return ((long)floored % 2 == 0) ? (long)floored : (long)floored + 1;
}

/* ------------------------------------------------------------------ */
/* Sleep                                                               */
/* ------------------------------------------------------------------ */

void gh_profile_sleep(gh_sleep_result *out, gh_arena *arena,
                      const gh_genome *g, const gh_analysis *analysis)
{
    memset(out, 0, sizeof(*out));

    double evening_score = 0.0, total_weight = 0.0;

    for (size_t i = 0; i < GH_CHRONOTYPE_SNP_COUNT; i++) {
        const gh_chronotype_snp *snp = &GH_CHRONOTYPE_SNPS[i];
        const char *genotype = gh_genome_genotype(g, snp->rsid);
        if (!genotype || !*genotype)
            continue;

        size_t copies = count_allele(genotype, snp->evening_allele);
        evening_score += (double)copies * snp->weight;
        total_weight += 2.0 * snp->weight;   /* the most this SNP could add */
        out->markers_found++;

        if (out->nmarkers < GH_DP_MAX_ITEMS) {
            gh_sleep_marker *m = &out->markers[out->nmarkers++];
            m->rsid = snp->rsid;
            m->gene = snp->gene;
            m->trait = snp->trait;
            m->genotype = genotype;
            m->evening_allele_copies = (int)copies;
            m->description = snp->description;
        }
    }

    out->chronotype_score = total_weight > 0.0
        ? gh_round(evening_score / total_weight * 100.0, 1) : 50.0;

    double score = out->chronotype_score;
    if (score >= 70) {
        out->chronotype = "Definite Evening (Night Owl)";
        out->optimal_sleep_window = "12:00 AM \342\200\223 8:00 AM";
        out->caffeine_cutoff = "12:00 PM (noon)";
        out->peak_alertness = "10:00 AM \342\200\223 2:00 PM and 7:00 PM \342\200\223 11:00 PM";
    } else if (score >= 55) {
        out->chronotype = "Moderate Evening";
        out->optimal_sleep_window = "11:30 PM \342\200\223 7:30 AM";
        out->caffeine_cutoff = "1:00 PM";
        out->peak_alertness = "10:00 AM \342\200\223 1:00 PM and 6:00 PM \342\200\223 10:00 PM";
    } else if (score >= 45) {
        out->chronotype = "Intermediate (Neither)";
        out->optimal_sleep_window = "11:00 PM \342\200\223 7:00 AM";
        out->caffeine_cutoff = "2:00 PM";
        out->peak_alertness = "9:00 AM \342\200\223 12:00 PM and 3:00 PM \342\200\223 7:00 PM";
    } else if (score >= 30) {
        out->chronotype = "Moderate Morning";
        out->optimal_sleep_window = "10:00 PM \342\200\223 6:00 AM";
        out->caffeine_cutoff = "2:00 PM";
        out->peak_alertness = "8:00 AM \342\200\223 12:00 PM";
    } else {
        out->chronotype = "Definite Morning (Early Bird)";
        out->optimal_sleep_window = "9:30 PM \342\200\223 5:30 AM";
        out->caffeine_cutoff = "12:00 PM (noon)";
        out->peak_alertness = "6:00 AM \342\200\223 11:00 AM";
    }

    /* Slow caffeine clearance overrides the chronotype's cutoff. */
    for (size_t i = 0; analysis && i < analysis->nfindings; i++) {
        const gh_finding *f = &analysis->findings[i];
        if (!f->gene || !f->status)
            continue;
        if (strcmp(f->gene, "CYP1A2") == 0 &&
            (strcmp(f->status, "slow") == 0 ||
             strcmp(f->status, "intermediate") == 0))
            out->caffeine_sensitive = true;
        if (strcmp(f->gene, "ADORA2A") == 0 &&
            strcmp(f->status, "anxiety_prone") == 0)
            out->caffeine_sensitive = true;
    }
    if (out->caffeine_sensitive)
        out->caffeine_cutoff = "10:00 AM or avoid entirely";

    out->confidence = out->markers_found >= 6 ? "high"
                    : out->markers_found >= 3 ? "moderate" : "low";

    char *line = gh_alloc(arena, 128);
    snprintf(line, 128, "Target sleep window: %s", out->optimal_sleep_window);
    out->recommendations[out->nrecommendations++] = line;

    line = gh_alloc(arena, 128);
    snprintf(line, 128, "Last caffeine by: %s", out->caffeine_cutoff);
    out->recommendations[out->nrecommendations++] = line;

    line = gh_alloc(arena, 160);
    snprintf(line, 160, "Peak alertness hours: %s", out->peak_alertness);
    out->recommendations[out->nrecommendations++] = line;

    if (score >= 60) {
        out->recommendations[out->nrecommendations++] =
            "Bright light exposure within 30 min of waking helps shift circadian rhythm earlier";
        out->recommendations[out->nrecommendations++] =
            "Blue light filter (f.lux/Night Shift) starting 2 hours before target bedtime";
        out->recommendations[out->nrecommendations++] =
            "Melatonin (0.3-1mg) 2-3 hours before target bedtime may help if shifting schedule";
    } else if (score <= 40) {
        out->recommendations[out->nrecommendations++] =
            "Avoid bright light in the evening to preserve early chronotype";
        out->recommendations[out->nrecommendations++] =
            "Early morning exercise reinforces morning circadian rhythm";
    }
    if (out->caffeine_sensitive)
        out->recommendations[out->nrecommendations++] =
            "Your caffeine metabolism genes suggest high sensitivity \342\200\224 "
            "consider switching to green tea or eliminating caffeine";

    out->deep_sleep_note = "";
    for (size_t i = 0; i < out->nmarkers; i++)
        if (strcmp(out->markers[i].gene, "ADA") == 0 &&
            out->markers[i].evening_allele_copies > 0)
            out->deep_sleep_note =
                "ADA variant detected: you likely build sleep pressure faster, "
                "meaning better deep sleep but higher sensitivity to sleep deprivation.";
}

/* ------------------------------------------------------------------ */
/* Nutrigenomics                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    gh_nutrient_need item;
    size_t order;
} ranked_need;

/* Largest absolute severity first; Python sorts on -abs(severity) with a
 * stable sort, so ties keep profile order. */
static int compare_need(const void *a, const void *b)
{
    const ranked_need *x = a, *y = b;
    double ax = fabs(x->item.severity), ay = fabs(y->item.severity);
    if (ax != ay)
        return ax > ay ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

void gh_profile_nutrigenomics(gh_nutrigenomics_result *out, gh_arena *arena,
                              const gh_analysis *analysis)
{
    memset(out, 0, sizeof(*out));

    ranked_need *ranked =
        gh_calloc(arena, GH_NUTRIENT_PROFILE_COUNT, sizeof(*ranked));

    for (size_t i = 0; i < GH_NUTRIENT_PROFILE_COUNT; i++) {
        const gh_nutrient_profile *profile = &GH_NUTRIENT_PROFILES[i];
        gh_nutrient_need *need = &ranked[i].item;
        need->profile = profile;
        ranked[i].order = i;

        double total = 0.0;
        bool has_signal = false;

        for (size_t gi = 0; gi < profile->ngenes; gi++) {
            const gh_nutrient_gene *gene = &profile->genes[gi];
            const char *status = last_status(analysis, gene->gene);
            if (!status)
                continue;

            for (size_t si = 0; si < gene->nstatuses; si++) {
                if (strcmp(gene->statuses[si].status, status) != 0)
                    continue;
                double severity = gene->statuses[si].weight;
                total += severity;
                has_signal = true;
                if (need->nimpacts < GH_DP_MAX_ITEMS) {
                    gh_nutrient_impact *impact = &need->impacts[need->nimpacts++];
                    impact->gene = gene->gene;
                    impact->status = status;
                    impact->impact = gene->impact;
                    impact->severity = severity;
                }
                break;
            }
        }

        char *recommendation = gh_alloc(arena, 160);
        if (!has_signal) {
            need->need_level = "normal";
            need->recommendation = "Standard dietary intake likely sufficient";
        } else if (total >= 0.8) {
            need->need_level = "high";
            snprintf(recommendation, 160, "Supplementation recommended: %s",
                     profile->supplement_form);
            need->recommendation = recommendation;
        } else if (total >= 0.4) {
            need->need_level = "moderate";
            snprintf(recommendation, 160,
                     "Dietary focus + consider supplementation: %s",
                     profile->supplement_form);
            need->recommendation = recommendation;
        } else if (total < 0) {
            /* Some genes push toward excess rather than deficiency. */
            need->need_level = "caution_excess";
            need->recommendation =
                "Caution: genetic tendency toward excess. Avoid supplementation "
                "unless deficient.";
        } else {
            need->need_level = "low";
            need->recommendation =
                "Minor genetic signal \342\200\224 dietary optimization sufficient";
        }

        need->severity = gh_round(total, 2);
    }

    qsort(ranked, GH_NUTRIENT_PROFILE_COUNT, sizeof(*ranked), compare_need);

    out->needs = gh_calloc(arena, GH_NUTRIENT_PROFILE_COUNT, sizeof(*out->needs));
    for (size_t i = 0; i < GH_NUTRIENT_PROFILE_COUNT; i++)
        out->needs[i] = ranked[i].item;
    out->nneeds = GH_NUTRIENT_PROFILE_COUNT;

    size_t high = 0, moderate = 0;
    for (size_t i = 0; i < out->nneeds; i++) {
        if (strcmp(out->needs[i].need_level, "high") == 0)
            high++;
        else if (strcmp(out->needs[i].need_level, "moderate") == 0)
            moderate++;
    }

    char *summary = gh_alloc(arena, 128);
    if (high >= 2)
        snprintf(summary, 128,
                 "Your genetics suggest increased needs for %zu key nutrients.",
                 high);
    else if (high == 1)
        snprintf(summary, 128,
                 "One nutrient shows high genetic need; %zu moderate.", moderate);
    else
        snprintf(summary, 128,
                 "No major nutrient deficiency risks detected from genetics.");
    out->summary = summary;
}

/* ------------------------------------------------------------------ */
/* Mental health                                                       */
/* ------------------------------------------------------------------ */

void gh_profile_mental_health(gh_mental_result *out, gh_arena *arena,
                              const gh_genome *g, const gh_analysis *analysis,
                              const gh_star_result *stars, size_t nstars)
{
    memset(out, 0, sizeof(*out));

    for (size_t i = 0; i < GH_MENTAL_SNP_COUNT; i++) {
        const gh_mental_snp *snp = &GH_MENTAL_SNPS[i];
        const char *genotype = gh_genome_genotype(g, snp->rsid);
        if (!genotype || !*genotype)
            continue;
        if (out->nmarkers >= GH_DP_MAX_ITEMS)
            break;
        gh_mental_marker *m = &out->markers[out->nmarkers++];
        m->snp = snp;
        m->genotype = genotype;
        m->risk_copies = (int)count_allele(genotype, snp->risk_allele);
    }

    for (size_t d = 0; d < GH_MENTAL_DOMAIN_COUNT; d++) {
        const gh_mental_domain *domain = &GH_MENTAL_DOMAINS[d];
        gh_mental_domain_result *result = &out->domains[out->ndomains++];
        result->name = domain->name;

        double risk = 0.0, protective = 0.0;

        for (size_t gi = 0; gi < domain->ngenes; gi++) {
            const gh_mental_gene *gene = &domain->genes[gi];
            const char *status = last_status(analysis, gene->gene);

            if (status && list_has(gene->risk_statuses, gene->nrisk_statuses,
                                   status)) {
                risk += 1.0;
                if (result->nsignals < 5) {
                    size_t len = strlen(gene->gene) + strlen(gene->impact) + 4;
                    char *signal = gh_alloc(arena, len);
                    snprintf(signal, len, "%s: %s", gene->gene, gene->impact);
                    result->signals[result->nsignals++] = signal;
                }
            } else if (status && *status && strcmp(status, "normal") != 0 &&
                       strcmp(status, "reference") != 0 &&
                       strcmp(status, "typical") != 0) {
                /* A non-risk, non-baseline status counts neither way. */
            } else {
                protective += 0.3;
            }
        }

        for (size_t mi = 0; mi < out->nmarkers; mi++) {
            const gh_mental_marker *m = &out->markers[mi];
            if (strcmp(m->snp->domain, domain->name) != 0 || m->risk_copies <= 0)
                continue;
            risk += m->risk_copies * 0.5;
            if (result->nsignals < 5) {
                size_t len = strlen(m->snp->gene) + strlen(m->snp->effect) + 4;
                char *signal = gh_alloc(arena, len);
                snprintf(signal, len, "%s: %s", m->snp->gene, m->snp->effect);
                result->signals[result->nsignals++] = signal;
            }
        }

        /* The denominator floors the protective term at 1, so a domain with
         * no protective signal still divides by something. */
        if (risk + protective > 0.0) {
            double denominator = risk + (protective > 1.0 ? protective : 1.0);
            result->risk_score = (int)round_half_even(risk / denominator * 100.0);
        } else {
            result->risk_score = 30;
        }

        result->risk_level = result->risk_score >= 60 ? "elevated"
                           : result->risk_score >= 40 ? "moderate" : "low";
    }

    const char *comt = last_status(analysis, "COMT");
    const char *bdnf = last_status(analysis, "BDNF");
    const char *aldh2 = last_status(analysis, "ALDH2");
    const char *slc6a4 = last_status(analysis, "SLC6A4");

    /* Several of these cut both ways, appearing as risk and resilience. */
    if (comt && strcmp(comt, "slow") == 0) {
        out->risk_factors[out->nrisk_factors++] =
            "Slow COMT: higher anxiety proneness, but better working memory";
        out->resilience_factors[out->nresilience_factors++] =
            "Slow COMT: enhanced cognitive performance under low stress";
    }
    if (bdnf && (strcmp(bdnf, "reduced") == 0 ||
                 strcmp(bdnf, "met_carrier") == 0)) {
        out->risk_factors[out->nrisk_factors++] =
            "BDNF Met carrier: reduced neuroplasticity under stress";
        out->resilience_factors[out->nresilience_factors++] =
            "BDNF Met: responds strongly to exercise-induced BDNF boost";
    }
    if (aldh2 && (strcmp(aldh2, "reduced") == 0 ||
                  strcmp(aldh2, "non_functional") == 0))
        out->resilience_factors[out->nresilience_factors++] =
            "ALDH2 variant: natural alcohol aversion (protective)";
    if (slc6a4 && (strcmp(slc6a4, "short") == 0 ||
                   strcmp(slc6a4, "reduced") == 0)) {
        out->risk_factors[out->nrisk_factors++] =
            "Short serotonin transporter: increased environmental sensitivity";
        out->resilience_factors[out->nresilience_factors++] =
            "Short 5-HTTLPR: also more responsive to positive environments";
    }

    for (size_t mi = 0; mi < out->nmarkers && out->nrisk_factors < 8; mi++) {
        const gh_mental_marker *m = &out->markers[mi];
        if (m->risk_copies >= 2) {
            size_t len = strlen(m->snp->gene) + strlen(m->snp->effect) + 32;
            char *line = gh_alloc(arena, len);
            snprintf(line, len, "%s: Homozygous risk \342\200\224 %s",
                     m->snp->gene, m->snp->effect);
            out->risk_factors[out->nrisk_factors++] = line;
        } else if (m->risk_copies == 1 &&
                   strcmp(m->snp->domain, "addiction") != 0) {
            size_t len = strlen(m->snp->gene) + strlen(m->snp->domain) + 48;
            char *line = gh_alloc(arena, len);
            snprintf(line, len, "%s: Heterozygous \342\200\224 moderate %s signal",
                     m->snp->gene, m->snp->domain);
            out->risk_factors[out->nrisk_factors++] = line;
        }
    }

    const char *cyp2c19 = "normal";
    for (size_t i = 0; i < nstars; i++)
        if (strcmp(stars[i].gene, "CYP2C19") == 0 && stars[i].phenotype)
            cyp2c19 = stars[i].phenotype;

    if (strcmp(cyp2c19, "poor") == 0 || strcmp(cyp2c19, "intermediate") == 0)
        out->treatment_notes[out->ntreatment_notes++] =
            "CYP2C19 reduced: Start SSRIs (citalopram, escitalopram, sertraline) "
            "at 50% dose. Consider alternatives less dependent on CYP2C19.";
    else if (strcmp(cyp2c19, "ultrarapid") == 0 || strcmp(cyp2c19, "rapid") == 0)
        out->treatment_notes[out->ntreatment_notes++] =
            "CYP2C19 ultrarapid: Standard SSRI doses may be insufficient. "
            "Consider higher doses or non-CYP2C19 alternatives.";

    if (comt && strcmp(comt, "slow") == 0)
        out->treatment_notes[out->ntreatment_notes++] =
            "Slow COMT: May be sensitive to stimulant medications and methyl donors. "
            "Mindfulness-based therapy particularly effective.";
    if (bdnf && (strcmp(bdnf, "reduced") == 0 ||
                 strcmp(bdnf, "met_carrier") == 0))
        out->treatment_notes[out->ntreatment_notes++] =
            "BDNF Met carrier: Exercise prescription is evidence-based \342\200\224 "
            "aerobic exercise 150+ min/week is as effective as SSRIs for "
            "mild-moderate depression.";

    out->recommendations[out->nrecommendations++] =
        "Regular aerobic exercise (strongest natural antidepressant \342\200\224 "
        "boosts BDNF, serotonin, dopamine)";
    out->recommendations[out->nrecommendations++] =
        "Sleep optimization (circadian disruption worsens all mental health domains)";
    out->recommendations[out->nrecommendations++] =
        "Social connection (protective across all genetic risk profiles)";

    for (size_t d = 0; d < out->ndomains; d++) {
        if (strcmp(out->domains[d].name, "anxiety") == 0 &&
            strcmp(out->domains[d].risk_level, "elevated") == 0)
            out->recommendations[out->nrecommendations++] =
                "Stress management: meditation, breathwork, or CBT-based anxiety protocol";
    }
    for (size_t d = 0; d < out->ndomains; d++) {
        if (strcmp(out->domains[d].name, "addiction") == 0 &&
            strcmp(out->domains[d].risk_level, "elevated") == 0)
            out->recommendations[out->nrecommendations++] =
                "Awareness of addiction vulnerability \342\200\224 set boundaries "
                "with alcohol/substances";
    }
    for (size_t mi = 0; mi < out->nmarkers; mi++)
        if (strcmp(out->markers[mi].snp->gene, "FKBP5") == 0 &&
            out->markers[mi].risk_copies > 0) {
            out->recommendations[out->nrecommendations++] =
                "Trauma-informed care: FKBP5 variant suggests altered stress recovery "
                "\342\200\224 consider EMDR or trauma-focused CBT if relevant";
            break;
        }

    const char *elevated[GH_DP_MAX_ITEMS];
    size_t nelevated = 0;
    for (size_t d = 0; d < out->ndomains; d++)
        if (strcmp(out->domains[d].risk_level, "elevated") == 0)
            elevated[nelevated++] = out->domains[d].name;

    if (nelevated) {
        size_t len = 96;
        for (size_t i = 0; i < nelevated; i++)
            len += strlen(elevated[i]) + 2;
        char *summary = gh_alloc(arena, len);
        size_t at = (size_t)snprintf(summary, len, "Elevated genetic signals in: ");
        for (size_t i = 0; i < nelevated; i++)
            at += (size_t)snprintf(summary + at, len - at, "%s%s",
                                   i ? ", " : "", elevated[i]);
        snprintf(summary + at, len - at, ". Proactive strategies recommended.");
        out->summary = summary;
    } else {
        out->summary = "No strongly elevated mental health genetic signals. "
                       "Maintain wellness practices.";
    }
}

/* ------------------------------------------------------------------ */
/* Longevity                                                           */
/* ------------------------------------------------------------------ */

void gh_profile_longevity(gh_longevity_result *out, gh_arena *arena,
                          const gh_genome *g, const gh_analysis *analysis,
                          const gh_apoe_result *apoe,
                          const gh_prs_result *prs, size_t nprs)
{
    memset(out, 0, sizeof(*out));

    size_t protective_count = 0;

    for (size_t i = 0; i < GH_LONGEVITY_SNP_COUNT; i++) {
        const gh_longevity_snp *snp = &GH_LONGEVITY_SNPS[i];
        const char *genotype = gh_genome_genotype(g, snp->rsid);
        if (!genotype || !*genotype)
            continue;

        out->alleles_checked++;
        size_t copies = count_allele(genotype, snp->longevity_allele);
        if (copies == 2)
            protective_count += 2;
        else if (copies == 1)
            protective_count += 1;

        if (out->nalleles < GH_DP_MAX_ITEMS) {
            gh_longevity_allele *a = &out->alleles[out->nalleles++];
            a->snp = snp;
            a->genotype = genotype;
            a->copies = (int)copies;
            a->status = copies == 2 ? "homozygous_protective"
                      : copies == 1 ? "heterozygous" : "absent";
        }
    }

    size_t max_possible = out->alleles_checked ? out->alleles_checked * 2 : 1;
    out->longevity_score =
        gh_round((double)protective_count / (double)max_possible * 100.0, 1);

    /* APOE only contributes when it carries an e4 allele, and only when the
     * findings did not already name APOE. */
    const char *apoe_status = NULL;
    if (apoe && apoe->apoe_type && strcmp(apoe->apoe_type, "Unknown") != 0 &&
        strstr(apoe->apoe_type, "e4")) {
        size_t len = strlen(apoe->apoe_type) + 1;
        char *converted = gh_alloc(arena, len);
        for (size_t i = 0; i < len; i++)
            converted[i] = apoe->apoe_type[i] == '/' ? '_' : apoe->apoe_type[i];
        apoe_status = converted;
    }

    for (size_t d = 0; d < GH_HEALTHSPAN_DOMAIN_COUNT; d++) {
        const gh_healthspan_domain *domain = &GH_HEALTHSPAN_DOMAINS[d];
        gh_healthspan_result *result = &out->domains[out->ndomains++];
        result->name = domain->name;

        int score = 50;   /* neutral baseline */
        for (size_t gi = 0; gi < domain->ngenes; gi++) {
            const char *gene = domain->genes[gi];
            const char *status = last_status(analysis, gene);
            if (!status && apoe_status && strcmp(gene, "APOE") == 0)
                status = apoe_status;
            if (!status || !*status)
                continue;
            result->genes_found++;
            if (list_has(domain->protective_statuses,
                         domain->nprotective_statuses, status))
                score += 8;
            else if (list_has(domain->risk_statuses, domain->nrisk_statuses,
                              status))
                score -= 10;
        }
        if (score < 10)
            score = 10;
        if (score > 90)
            score = 90;
        result->score = score;
        result->rating = score >= 60 ? "good" : score >= 40 ? "average"
                                                            : "attention";
    }

    for (size_t i = 0; i < nprs; i++) {
        const gh_prs_result *r = &prs[i];
        bool elevated = strcmp(r->risk_category, "elevated") == 0 ||
                        strcmp(r->risk_category, "high") == 0;
        bool low = strcmp(r->risk_category, "low") == 0;
        if (!elevated && !low)
            continue;

        size_t len = strlen(r->name) + 64;
        char *line = gh_alloc(arena, len);
        snprintf(line, len, "%s PRS for %s (%.0fth percentile)",
                 elevated ? "Elevated" : "Low", r->name, r->percentile);
        if (elevated && out->ntop_risks < 8)
            out->top_risks[out->ntop_risks++] = line;
        else if (low && out->ntop_protective < 8)
            out->top_protective[out->ntop_protective++] = line;
    }

    if (apoe && apoe->risk_level) {
        if (strcmp(apoe->risk_level, "elevated") == 0 ||
            strcmp(apoe->risk_level, "high") == 0) {
            size_t len = strlen(apoe->apoe_type) + strlen(apoe->description) + 16;
            char *line = gh_alloc(arena, len);
            snprintf(line, len, "APOE %s: %s", apoe->apoe_type,
                     apoe->description);
            if (out->ntop_risks < 8)
                out->top_risks[out->ntop_risks++] = line;
        } else if (strcmp(apoe->risk_level, "reduced") == 0) {
            size_t len = strlen(apoe->apoe_type) + 48;
            char *line = gh_alloc(arena, len);
            snprintf(line, len, "APOE %s: Protective for Alzheimer's",
                     apoe->apoe_type);
            if (out->ntop_protective < 8)
                out->top_protective[out->ntop_protective++] = line;
        }
    }

    for (size_t i = 0; i < out->nalleles; i++) {
        const gh_longevity_allele *a = &out->alleles[i];
        size_t len = strlen(a->snp->gene) + 48;
        char *line = gh_alloc(arena, len);
        if (a->copies == 2 && out->ntop_protective < 8) {
            snprintf(line, len, "%s: Homozygous for longevity allele",
                     a->snp->gene);
            out->top_protective[out->ntop_protective++] = line;
        } else if (a->copies == 0 && out->ntop_risks < 8) {
            snprintf(line, len, "%s: No longevity allele copies", a->snp->gene);
            out->top_risks[out->ntop_risks++] = line;
        }
    }

    if (last_status(analysis, "PPARGC1A") || last_status(analysis, "BDNF"))
        out->interventions[out->ninterventions++] = (gh_intervention){
            "Regular aerobic exercise (150+ min/week)", "high",
            "Strongest BDNF booster + mitochondrial biogenesis via PGC-1\316\261"};

    out->interventions[out->ninterventions++] = (gh_intervention){
        "Mediterranean dietary pattern", "high",
        "Reduces inflammation (TNF/IL6), supports APOE-linked cardiovascular health"};

    for (size_t i = 0; i < out->nalleles; i++)
        if (strcmp(out->alleles[i].snp->gene, "FOXO3") == 0 &&
            out->alleles[i].copies > 0) {
            out->interventions[out->ninterventions++] = (gh_intervention){
                "Time-restricted eating (12-16hr overnight fast)", "moderate",
                "FOXO3 protective allele activates autophagy pathways enhanced by fasting"};
            break;
        }
    for (size_t i = 0; i < out->nalleles; i++)
        if (strcmp(out->alleles[i].snp->gene, "SIRT1") == 0) {
            out->interventions[out->ninterventions++] = (gh_intervention){
                "Resveratrol / NAD+ precursors (NMN, NR)", "moderate",
                "SIRT1 pathway responds to NAD+ boosting compounds"};
            break;
        }

    out->interventions[out->ninterventions++] = (gh_intervention){
        "Sleep optimization (7-9 hrs, consistent schedule)", "high",
        "Telomere maintenance (TERC/TERT) requires adequate sleep"};
    out->interventions[out->ninterventions++] = (gh_intervention){
        "Stress management (meditation, social connection)", "moderate",
        "Reduces cortisol-driven telomere shortening and inflammation"};

    out->summary = out->longevity_score >= 65
        ? "Your longevity genetic profile is favorable."
        : out->longevity_score >= 40
        ? "Your longevity profile is average \342\200\224 lifestyle optimization has high impact."
        : "Your longevity profile shows some risk factors \342\200\224 proactive interventions recommended.";
}

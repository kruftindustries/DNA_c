#include "gh_dosing.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Shared helpers                                                      */
/* ------------------------------------------------------------------ */

static bool value_in(const char *const *values, size_t n, const char *value)
{
    if (!value)
        return false;
    for (size_t i = 0; i < n; i++)
        if (strcmp(values[i], value) == 0)
            return true;
    return false;
}

/* The phenotype call_star_alleles reported for a gene, or NULL when the gene
 * is not in the results at all. Python spells this
 * `star_alleles.get(gene, {}).get("phenotype")`, so a missing gene and a
 * missing phenotype both come out as None. */
static const char *star_phenotype(const gh_star_result *stars, size_t nstars,
                                  const char *gene)
{
    for (size_t i = 0; i < nstars; i++)
        if (stars[i].gene && strcmp(stars[i].gene, gene) == 0)
            return stars[i].phenotype;
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Drug dosing                                                         */
/* ------------------------------------------------------------------ */

bool gh_dose_is_critical(const char *action)
{
    static const char *const keywords[] = {
        "CONTRAINDICATED", "FATAL", "AVOID", "LIFE-THREATENING",
    };

    if (!action)
        return false;

    /* Python tests `w in rule["action"].upper()`. Only ASCII letters change
     * under str.upper() for these texts, and the keywords are ASCII, so an
     * ASCII-only fold is equivalent here. */
    size_t len = strlen(action);
    for (size_t k = 0; k < sizeof(keywords) / sizeof(keywords[0]); k++) {
        size_t klen = strlen(keywords[k]);
        if (klen > len)
            continue;
        for (size_t i = 0; i + klen <= len; i++) {
            size_t j = 0;
            while (j < klen) {
                char c = action[i + j];
                if (c >= 'a' && c <= 'z')
                    c = (char)(c - 'a' + 'A');
                if (c != keywords[k][j])
                    break;
                j++;
            }
            if (j == klen)
                return true;
        }
    }
    return false;
}

static bool dose_rule_fires(const gh_dose_rule *rule,
                            const gh_star_result *stars, size_t nstars,
                            const gh_analysis *analysis)
{
    if (rule->source == GH_DOSE_STAR)
        return value_in(rule->values, rule->nvalues,
                        star_phenotype(stars, nstars, rule->gene));

    /* GH_DOSE_FINDING: any lifestyle finding for the gene, not the last or
     * the strongest -- the Python is an `any(...)` over the whole list. */
    for (size_t i = 0; analysis && i < analysis->nfindings; i++) {
        const gh_finding *f = &analysis->findings[i];
        if (f->gene && strcmp(f->gene, rule->gene) == 0 &&
            value_in(rule->values, rule->nvalues, f->status))
            return true;
    }
    return false;
}

void gh_generate_drug_dosing(gh_dosing_result *out, gh_arena *arena,
                             const gh_star_result *stars, size_t nstars,
                             const gh_analysis *analysis)
{
    memset(out, 0, sizeof(*out));

    if (nstars == 0) {
        out->summary = "No pharmacogenomic data available for drug dosing.";
        return;
    }

    size_t cap = 0;
    for (size_t i = 0; i < GH_DOSE_DRUG_COUNT; i++)
        cap += GH_DOSE_DRUGS[i].nrules;
    out->recommendations = gh_calloc(arena, cap, sizeof(*out->recommendations));
    out->warnings = gh_calloc(arena, cap, sizeof(*out->warnings));

    for (size_t i = 0; i < GH_DOSE_DRUG_COUNT; i++) {
        const gh_dose_drug *drug = &GH_DOSE_DRUGS[i];

        /* Skipped outright unless at least one of the drug's genes was
         * called, even when a rule would have read a lifestyle finding. */
        bool has_data = false;
        for (size_t g = 0; g < drug->ngenes && !has_data; g++)
            has_data = star_phenotype(stars, nstars, drug->genes[g]) != NULL;
        if (!has_data)
            continue;

        for (size_t r = 0; r < drug->nrules; r++) {
            const gh_dose_rule *rule = &drug->rules[r];
            if (!dose_rule_fires(rule, stars, nstars, analysis))
                continue;

            gh_dose_rec *rec = &out->recommendations[out->nrecommendations++];
            rec->drug = drug;
            rec->rule = rule;
            rec->critical = gh_dose_is_critical(rule->action);
            if (rec->critical)
                out->warnings[out->nwarnings++] = rec;
        }
    }

    char *summary = gh_alloc(arena, 128);
    if (out->nwarnings)
        snprintf(summary, 128,
                 "%zu drug dosing recommendations, including %zu critical "
                 "warning(s).", out->nrecommendations, out->nwarnings);
    else if (out->nrecommendations)
        snprintf(summary, 128,
                 "%zu drug dosing recommendations based on your "
                 "pharmacogenomic profile.", out->nrecommendations);
    else
        snprintf(summary, 128, "No drug dosing adjustments needed based on "
                               "available pharmacogenomic data.");
    out->summary = summary;
}

/* ------------------------------------------------------------------ */
/* Polypharmacy                                                        */
/* ------------------------------------------------------------------ */

static int poly_severity_rank(const char *severity)
{
    if (strcmp(severity, "high") == 0)
        return 0;
    if (strcmp(severity, "moderate") == 0)
        return 1;
    if (strcmp(severity, "low") == 0)
        return 2;
    return 3;
}

typedef struct {
    gh_poly_warning item;
    size_t order;
} ranked_poly;

static int compare_poly(const void *a, const void *b)
{
    const ranked_poly *x = a, *y = b;
    int rx = poly_severity_rank(x->item.rule->severity);
    int ry = poly_severity_rank(y->item.rule->severity);
    if (rx != ry)
        return rx < ry ? -1 : 1;
    /* Python's sort is stable, so equal severities keep rule order. */
    return (x->order > y->order) - (x->order < y->order);
}

/* The status the rules attribute to a gene.
 *
 * Built from the star alleles first and then overwritten by the lifestyle
 * findings, so a gene named by both is judged on its finding. Within the
 * findings the LAST one wins, and since findings arrive sorted by descending
 * magnitude that is the weakest call -- the same quirk the dependent profiles
 * carry, reproduced deliberately. */
static const char *poly_phenotype(const gh_star_result *stars, size_t nstars,
                                  const gh_analysis *analysis, const char *gene)
{
    const char *phenotype = NULL;

    for (size_t i = 0; i < nstars; i++) {
        if (!stars[i].gene || strcmp(stars[i].gene, gene) != 0)
            continue;
        /* "Unknown" means the gene could not be called; Python leaves it out
         * of the lookup entirely rather than storing it. */
        if (stars[i].phenotype && strcmp(stars[i].phenotype, "Unknown") != 0)
            phenotype = stars[i].phenotype;
    }

    for (size_t i = 0; analysis && i < analysis->nfindings; i++) {
        const gh_finding *f = &analysis->findings[i];
        if (f->gene && f->status && *f->gene && *f->status &&
            strcmp(f->gene, gene) == 0)
            phenotype = f->status;
    }

    return phenotype;
}

void gh_assess_polypharmacy(gh_polypharmacy_result *out, gh_arena *arena,
                            const gh_star_result *stars, size_t nstars,
                            const gh_analysis *analysis)
{
    memset(out, 0, sizeof(*out));

    ranked_poly *ranked = gh_calloc(arena, GH_POLY_RULE_COUNT, sizeof(*ranked));
    size_t n = 0;

    for (size_t i = 0; i < GH_POLY_RULE_COUNT; i++) {
        const gh_poly_rule *rule = &GH_POLY_RULES[i];
        gh_poly_warning w;
        memset(&w, 0, sizeof(w));
        w.rule = rule;

        bool match = true;
        for (size_t g = 0; g < rule->ngenes; g++) {
            const gh_poly_gene *gene = &rule->genes[g];
            const char *phenotype =
                poly_phenotype(stars, nstars, analysis, gene->gene);
            if (phenotype && *phenotype &&
                value_in(gene->phenotypes, gene->nphenotypes, phenotype)) {
                w.matched[g] = phenotype;
            } else {
                match = false;
                break;
            }
        }

        if (match) {
            ranked[n].item = w;
            ranked[n].order = n;
            n++;
        }
    }

    qsort(ranked, n, sizeof(*ranked), compare_poly);

    out->warnings = gh_calloc(arena, n ? n : 1, sizeof(*out->warnings));
    out->nwarnings = n;
    for (size_t i = 0; i < n; i++)
        out->warnings[i] = ranked[i].item;

    /* by_severity: distinct severities in the order the sorted list reaches
     * them, which is Python's setdefault over the sorted warnings. */
    out->severities = gh_calloc(arena, n ? n : 1, sizeof(*out->severities));
    out->severity_counts =
        gh_calloc(arena, n ? n : 1, sizeof(*out->severity_counts));
    for (size_t i = 0; i < n; i++) {
        const char *severity = out->warnings[i].rule->severity;
        size_t slot = out->nseverities;
        for (size_t k = 0; k < out->nseverities; k++) {
            if (strcmp(out->severities[k], severity) == 0) {
                slot = k;
                break;
            }
        }
        if (slot == out->nseverities)
            out->severities[out->nseverities++] = severity;
        out->severity_counts[slot]++;
    }
}

/* ------------------------------------------------------------------ */
/* Preventive care                                                     */
/* ------------------------------------------------------------------ */

int gh_screening_priority_rank(const char *priority)
{
    if (!priority)
        return 5;
    if (strcmp(priority, "urgent") == 0)
        return 0;
    if (strcmp(priority, "high") == 0)
        return 1;
    if (strcmp(priority, "elevated") == 0)
        return 2;
    if (strcmp(priority, "ongoing") == 0)
        return 3;
    if (strcmp(priority, "standard") == 0)
        return 4;
    return 5;
}

typedef struct {
    gh_screening item;
    size_t order;
} ranked_screening;

static int compare_screening(const void *a, const void *b)
{
    const ranked_screening *x = a, *y = b;
    if (x->item.start_age != y->item.start_age)
        return x->item.start_age < y->item.start_age ? -1 : 1;
    int px = gh_screening_priority_rank(x->item.priority);
    int py = gh_screening_priority_rank(y->item.priority);
    if (px != py)
        return px < py ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

static const gh_screening_modifier *find_modifier(
    const gh_screening_modifier *table, size_t n, const char *condition)
{
    for (size_t i = 0; i < n; i++)
        if (strcmp(table[i].condition, condition) == 0)
            return &table[i];
    return NULL;
}

/* "PRS 87th percentile (elevated)". Python formats the percentile with
 * `:.0f`, which rounds half to even exactly as printf does. */
static const char *prs_basis(gh_arena *arena, double percentile,
                             const char *category)
{
    char *text = gh_alloc(arena, 64);
    snprintf(text, 64, "PRS %.0fth percentile (%s)", percentile, category);
    return text;
}

void gh_generate_preventive_timeline(gh_preventive_result *out, gh_arena *arena,
                                     const gh_prs_result *prs, size_t nprs,
                                     const gh_apoe_result *apoe,
                                     const gh_acmg_result *acmg,
                                     const gh_star_result *stars, size_t nstars)
{
    memset(out, 0, sizeof(*out));

    /* The base schedule, keyed by condition and mutated in place by the PRS
     * pass. Kept parallel to GH_BASE_SCREENINGS so insertion order survives,
     * since the entries are appended to the timeline in that order. */
    gh_screening *base =
        gh_calloc(arena, GH_BASE_SCREENING_COUNT, sizeof(*base));
    for (size_t i = 0; i < GH_BASE_SCREENING_COUNT; i++) {
        const gh_base_screening *s = &GH_BASE_SCREENINGS[i];
        char *reason = gh_alloc(arena, 96);
        snprintf(reason, 96, "General population guideline (%s)", s->source);
        base[i].test = s->test;
        base[i].start_age = s->start_age;
        base[i].frequency = s->frequency;
        base[i].reason = reason;
        base[i].priority = "standard";
        base[i].genetic_basis = NULL;
    }

    size_t cap = GH_BASE_SCREENING_COUNT + nprs + 3 + 1;
    if (acmg)
        cap += acmg->nfindings * 2;
    ranked_screening *ranked = gh_calloc(arena, cap, sizeof(*ranked));
    size_t n = 0;

    /* --- polygenic scores --- */
    for (size_t i = 0; i < nprs; i++) {
        const gh_prs_result *result = &prs[i];
        const char *category = result->risk_category ? result->risk_category
                                                     : "average";

        const gh_screening_modifier *mod = NULL;
        if (strcmp(category, "high") == 0)
            mod = find_modifier(GH_PRS_HIGH_MODIFIERS,
                                GH_PRS_HIGH_MODIFIERS_COUNT, result->id);
        else if (strcmp(category, "elevated") == 0)
            mod = find_modifier(GH_PRS_ELEVATED_MODIFIERS,
                                GH_PRS_ELEVATED_MODIFIERS_COUNT, result->id);
        if (!mod)
            continue;

        gh_screening *entry = NULL;
        for (size_t b = 0; b < GH_BASE_SCREENING_COUNT; b++) {
            if (strcmp(GH_BASE_SCREENINGS[b].condition, result->id) == 0) {
                entry = &base[b];
                break;
            }
        }

        if (entry) {
            int start = entry->start_age + mod->start_age_delta;
            entry->start_age = start > 18 ? start : 18;
            entry->frequency = mod->frequency;
            entry->reason = mod->reason;
            entry->priority = strcmp(category, "high") == 0 ? "high"
                                                           : "elevated";
            entry->genetic_basis = prs_basis(arena, result->percentile,
                                             category);
        } else {
            /* A condition with a modifier but no base screening becomes a new
             * entry at 30. Its priority is hard-coded "elevated" in the
             * Python even when the category is "high"; kept as is. */
            char *test = gh_alloc(arena, 96);
            snprintf(test, 96, "%s screening", result->name);
            ranked[n].item.test = test;
            ranked[n].item.start_age = 30;
            ranked[n].item.frequency = mod->frequency;
            ranked[n].item.reason = mod->reason;
            ranked[n].item.priority = "elevated";
            ranked[n].item.genetic_basis =
                prs_basis(arena, result->percentile, category);
            ranked[n].order = n;
            n++;
        }
        out->early_screenings++;
    }

    /* --- APOE --- */
    if (apoe && apoe->risk_level) {
        const gh_apoe_screening *table = NULL;
        size_t ntable = 0;
        if (strcmp(apoe->risk_level, "elevated") == 0) {
            table = GH_APOE_ELEVATED_SCREENINGS;
            ntable = GH_APOE_ELEVATED_SCREENINGS_COUNT;
        } else if (strcmp(apoe->risk_level, "high") == 0) {
            table = GH_APOE_HIGH_SCREENINGS;
            ntable = GH_APOE_HIGH_SCREENINGS_COUNT;
        }

        for (size_t i = 0; i < ntable; i++) {
            char *basis = gh_alloc(arena, 64);
            snprintf(basis, 64, "APOE %s", apoe->apoe_type);
            ranked[n].item.test = table[i].test;
            ranked[n].item.start_age = table[i].start_age;
            ranked[n].item.frequency = table[i].frequency;
            ranked[n].item.reason = table[i].reason;
            ranked[n].item.priority =
                strcmp(apoe->risk_level, "high") == 0 ? "high" : "elevated";
            ranked[n].item.genetic_basis = basis;
            ranked[n].order = n;
            n++;
            out->early_screenings++;
        }
    }

    /* --- ACMG secondary findings --- */
    for (size_t i = 0; acmg && i < acmg->nfindings; i++) {
        const char *gene = acmg->findings[i].finding->gene;
        if (!gene)
            gene = "";

        bool brca = strcmp(gene, "BRCA1") == 0 || strcmp(gene, "BRCA2") == 0;
        bool lynch = strcmp(gene, "MLH1") == 0 || strcmp(gene, "MSH2") == 0 ||
                     strcmp(gene, "MSH6") == 0 || strcmp(gene, "PMS2") == 0;
        if (!brca && !lynch)
            continue;

        char *basis = gh_alloc(arena, 96);
        snprintf(basis, 96, "%s pathogenic variant", gene);

        if (brca) {
            char *breast = gh_alloc(arena, 160);
            snprintf(breast, 160, "Pathogenic %s variant — enhanced breast "
                                  "cancer screening per NCCN", gene);
            ranked[n].item.test = "Breast MRI + mammography";
            ranked[n].item.start_age = 25;
            ranked[n].item.frequency = "Every 6 months (alternating)";
            ranked[n].item.reason = breast;
            ranked[n].item.priority = "urgent";
            ranked[n].item.genetic_basis = basis;
            ranked[n].order = n;
            n++;

            char *ovarian = gh_alloc(arena, 160);
            snprintf(ovarian, 160, "Pathogenic %s variant — discuss "
                                   "risk-reducing surgery", gene);
            ranked[n].item.test =
                "Ovarian cancer screening (CA-125 + ultrasound)";
            ranked[n].item.start_age = 30;
            ranked[n].item.frequency = "Every 6-12 months";
            ranked[n].item.reason = ovarian;
            ranked[n].item.priority = "urgent";
            ranked[n].item.genetic_basis = basis;
            ranked[n].order = n;
            n++;
            out->early_screenings += 2;
        } else {
            char *reason = gh_alloc(arena, 160);
            snprintf(reason, 160, "Lynch syndrome (%s) — early colonoscopy "
                                  "per NCCN", gene);
            ranked[n].item.test = "Colonoscopy";
            ranked[n].item.start_age = 20;
            ranked[n].item.frequency = "Every 1-2 years";
            ranked[n].item.reason = reason;
            ranked[n].item.priority = "urgent";
            ranked[n].item.genetic_basis = basis;
            ranked[n].order = n;
            n++;
            out->early_screenings++;
        }
    }

    /* --- pharmacogenomic card --- */
    if (nstars) {
        size_t len = 40;
        size_t nnon_normal = 0;
        for (size_t i = 0; i < nstars; i++) {
            const char *phenotype = stars[i].phenotype;
            if (!phenotype || strcmp(phenotype, "normal") == 0 ||
                strcmp(phenotype, "Unknown") == 0)
                continue;
            len += strlen(stars[i].gene) + 2;
            nnon_normal++;
        }

        if (nnon_normal) {
            char *reason = gh_alloc(arena, len);
            size_t used = (size_t)snprintf(reason, len,
                                           "Non-standard metabolizer for: ");
            bool first = true;
            for (size_t i = 0; i < nstars; i++) {
                const char *phenotype = stars[i].phenotype;
                if (!phenotype || strcmp(phenotype, "normal") == 0 ||
                    strcmp(phenotype, "Unknown") == 0)
                    continue;
                used += (size_t)snprintf(reason + used, len - used, "%s%s",
                                         first ? "" : ", ", stars[i].gene);
                first = false;
            }

            ranked[n].item.test =
                "Pharmacogenomic card review (bring to every prescriber)";
            ranked[n].item.start_age = 18;
            ranked[n].item.frequency = "Every visit";
            ranked[n].item.reason = reason;
            ranked[n].item.priority = "ongoing";
            ranked[n].item.genetic_basis = "Pharmacogenomic profile";
            ranked[n].order = n;
            n++;
        }
    }

    /* Base screenings go on last, so an equal (age, priority) tie resolves in
     * favour of the genetically-driven entries above. */
    for (size_t i = 0; i < GH_BASE_SCREENING_COUNT; i++) {
        ranked[n].item = base[i];
        ranked[n].order = n;
        n++;
    }

    qsort(ranked, n, sizeof(*ranked), compare_screening);

    out->timeline = gh_calloc(arena, n ? n : 1, sizeof(*out->timeline));
    out->ntimeline = n;
    for (size_t i = 0; i < n; i++)
        out->timeline[i] = ranked[i].item;

    if (out->early_screenings) {
        char *summary = gh_alloc(arena, 160);
        snprintf(summary, 160,
                 "%zu screening(s) recommended earlier or more frequently "
                 "than general population guidelines based on your genetic "
                 "profile.", out->early_screenings);
        out->summary = summary;
    } else {
        out->summary = "Your genetic profile does not indicate need for "
                       "earlier screenings beyond standard guidelines.";
    }
}

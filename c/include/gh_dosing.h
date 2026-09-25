/* Drug dosing, polypharmacy and preventive care.
 *
 * Mirrors genetic_health/drug_dosing.py, polypharmacy.py and
 * preventive_care.py. All three read results the earlier modules produced --
 * star alleles, lifestyle findings, polygenic scores, APOE and ACMG -- rather
 * than the genome directly.
 *
 * The tables are generated from the Python by tools/gen_dosing_data.py.
 */
#ifndef GH_DOSING_H
#define GH_DOSING_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_analysis.h"
#include "gh_clinvar.h"
#include "gh_mem.h"
#include "gh_prs.h"
#include "gh_scorers.h"

/* ------------------------------------------------------------------ */
/* Drug dosing                                                         */
/* ------------------------------------------------------------------ */

/* Where a dosing rule looks for the status it tests. Python writes each rule
 * as a lambda over (star_alleles, lifestyle_findings); every one of them
 * reduces to one of these two lookups, which is what the generator asserts. */
typedef enum {
    GH_DOSE_STAR = 0,     /* star_alleles[gene]["phenotype"] */
    GH_DOSE_FINDING,      /* any lifestyle finding for gene */
} gh_dose_source;

typedef struct {
    gh_dose_source source;
    const char *gene;
    const char *const *values;   /* phenotypes/statuses that satisfy the rule */
    size_t nvalues;
    const char *action;
    const char *dose_guidance;
} gh_dose_rule;

typedef struct {
    const char *drug;        /* table key, "warfarin" */
    const char *display;     /* as reported: drug.replace("_", " ").title() */
    const char *category;
    const char *const *genes;
    size_t ngenes;
    const gh_dose_rule *rules;
    size_t nrules;
    const char *source;      /* guideline citation */
} gh_dose_drug;

extern const gh_dose_drug GH_DOSE_DRUGS[];
extern const size_t GH_DOSE_DRUG_COUNT;

typedef struct {
    const gh_dose_drug *drug;
    const gh_dose_rule *rule;
    bool critical;           /* also listed under warnings */
} gh_dose_rec;

typedef struct {
    gh_dose_rec *recommendations;
    size_t nrecommendations;
    /* The critical subset, in the same order; pointers into the array above
     * rather than copies, matching Python appending the same dict to both. */
    const gh_dose_rec **warnings;
    size_t nwarnings;
    const char *summary;
} gh_dosing_result;

/* `stars` and `analysis` may be NULL/empty; an empty star-allele set produces
 * the "no pharmacogenomic data" summary, as in Python. */
void gh_generate_drug_dosing(gh_dosing_result *out, gh_arena *arena,
                             const gh_star_result *stars, size_t nstars,
                             const gh_analysis *analysis);

/* True when an action's upper-cased text names one of the four keywords that
 * promote a recommendation to a warning. */
bool gh_dose_is_critical(const char *action);

/* ------------------------------------------------------------------ */
/* Polypharmacy                                                        */
/* ------------------------------------------------------------------ */

/* No rule tests more than this many genes; the generator enforces it. */
#define GH_POLY_MAX_GENES 4

typedef struct {
    const char *gene;
    const char *const *phenotypes;   /* qualifying phenotypes/statuses */
    size_t nphenotypes;
} gh_poly_gene;

typedef struct {
    const char *id;
    const char *name;
    const char *severity;            /* high / moderate / low */
    const gh_poly_gene *genes;
    size_t ngenes;
    const char *const *drugs_affected;
    size_t ndrugs;
    const char *warning;
    const char *action;
} gh_poly_rule;

extern const gh_poly_rule GH_POLY_RULES[];
extern const size_t GH_POLY_RULE_COUNT;

typedef struct {
    const gh_poly_rule *rule;
    /* The phenotype that matched each of the rule's genes, in rule order.
     * A rule fires only when every gene matches, so this is always ngenes
     * long and parallel to rule->genes. */
    const char *matched[GH_POLY_MAX_GENES];
} gh_poly_warning;

typedef struct {
    gh_poly_warning *warnings;       /* sorted by severity */
    size_t nwarnings;
    /* Distinct severities in the order the sorted list first reaches them,
     * with a count each -- Python's by_severity dict. */
    const char **severities;
    size_t *severity_counts;
    size_t nseverities;
} gh_polypharmacy_result;

void gh_assess_polypharmacy(gh_polypharmacy_result *out, gh_arena *arena,
                            const gh_star_result *stars, size_t nstars,
                            const gh_analysis *analysis);

/* ------------------------------------------------------------------ */
/* Preventive care                                                     */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *test;
    int start_age;
    const char *frequency;
    const char *condition;
    const char *source;
} gh_base_screening;

extern const gh_base_screening GH_BASE_SCREENINGS[];
extern const size_t GH_BASE_SCREENING_COUNT;

typedef struct {
    const char *condition;           /* PRS condition id */
    int start_age_delta;
    const char *frequency;
    const char *reason;
} gh_screening_modifier;

extern const gh_screening_modifier GH_PRS_ELEVATED_MODIFIERS[];
extern const size_t GH_PRS_ELEVATED_MODIFIERS_COUNT;
extern const gh_screening_modifier GH_PRS_HIGH_MODIFIERS[];
extern const size_t GH_PRS_HIGH_MODIFIERS_COUNT;

typedef struct {
    const char *test;
    int start_age;
    const char *frequency;
    const char *reason;
} gh_apoe_screening;

extern const gh_apoe_screening GH_APOE_ELEVATED_SCREENINGS[];
extern const size_t GH_APOE_ELEVATED_SCREENINGS_COUNT;
extern const gh_apoe_screening GH_APOE_HIGH_SCREENINGS[];
extern const size_t GH_APOE_HIGH_SCREENINGS_COUNT;

typedef struct {
    const char *test;
    int start_age;
    const char *frequency;
    const char *reason;
    const char *priority;            /* urgent/high/elevated/ongoing/standard */
    const char *genetic_basis;       /* NULL for the base guideline entries */
} gh_screening;

typedef struct {
    gh_screening *timeline;          /* sorted by start age, then priority */
    size_t ntimeline;
    const char *summary;
    size_t early_screenings;
} gh_preventive_result;

/* Any argument may be NULL. `carrier_screen` is deliberately absent: the
 * Python accepts it and never reads it. */
void gh_generate_preventive_timeline(gh_preventive_result *out, gh_arena *arena,
                                     const gh_prs_result *prs, size_t nprs,
                                     const gh_apoe_result *apoe,
                                     const gh_acmg_result *acmg,
                                     const gh_star_result *stars, size_t nstars);

/* Sort rank of a priority label; 5 for anything unrecognised, as Python's
 * `priority_order.get(x, 5)` gives. */
int gh_screening_priority_rank(const char *priority);

#endif /* GH_DOSING_H */

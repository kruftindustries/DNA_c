/* Prioritised recommendations built from every other analysis.
 *
 * Mirrors genetic_health/recommendations.py. Each risk group looks for
 * converging signals -- gene variants, polygenic scores, ClinVar findings and
 * gene-gene interactions -- and groups with at least one signal become a
 * prioritised entry with actions, a doctor note and a monitoring schedule.
 *
 * The Python's `ancestry_results` parameter is accepted but never read, so
 * it has no counterpart here.
 */
#ifndef GH_RECOMMENDATIONS_H
#define GH_RECOMMENDATIONS_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_analysis.h"
#include "gh_clinical_context.h"
#include "gh_clinvar.h"
#include "gh_epistasis.h"
#include "gh_mem.h"
#include "gh_prs.h"
#include "gh_scorers.h"

/* ------------------------------------------------------------------ */
/* Generated tables                                                    */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *test;
    const char *frequency;
    const char *reason;
} gh_monitoring_item;

typedef struct {
    const char *id;
    const char *title;
    const char *prs_condition;      /* "" when the group has no PRS model */
    const char *const *genes;       /* sorted, so output order is stable */
    size_t ngenes;
    const char *const *keywords;
    size_t nkeywords;
    const char *const *actions;
    size_t nactions;
    const char *doctor_note;
    const gh_monitoring_item *monitoring;
    size_t nmonitoring;
} gh_risk_group;

extern const gh_risk_group GH_RISK_GROUPS[];
extern const size_t GH_RISK_GROUP_COUNT;

typedef struct {
    const char *source;             /* acmg / pathogenic_count / disease / ... */
    bool match_any;
    int minimum;
    const char *gene;
    const char *group;
    const char *min_priority;
    const char *const *genes;
    size_t ngenes;
    const char *const *keywords;
    size_t nkeywords;
    const char *const *statuses;
    size_t nstatuses;
} gh_referral_trigger;

typedef struct {
    const char *id;
    const char *title;
    const char *reason_template;
    const char *urgency;            /* soon / routine */
    const gh_referral_trigger *triggers;
    size_t ntriggers;
} gh_specialist_referral;

extern const gh_specialist_referral GH_SPECIALIST_REFERRALS[];
extern const size_t GH_SPECIALIST_REFERRAL_COUNT;

/* ------------------------------------------------------------------ */
/* Results                                                             */
/* ------------------------------------------------------------------ */

#define GH_MAX_PRIORITY_ACTIONS 16
#define GH_MAX_CLINICAL_ACTIONS 12

typedef struct {
    const char *id;
    const char *title;
    const char *priority;           /* high / moderate / low */
    const char *why;
    const char *const *actions;
    size_t nactions;
    const char *doctor_note;
    const gh_monitoring_item *monitoring;
    size_t nmonitoring;
    size_t signal_count;

    /* Extra actions drawn from the clinical context table. */
    const char *clinical_actions[GH_MAX_CLINICAL_ACTIONS];
    size_t nclinical_actions;
} gh_priority;

typedef struct {
    const char *rsid;
    const char *genotype;
    const char *status;
    const char *description;
    const char *source;
} gh_drug_entry;

typedef struct {
    const char *gene;
    gh_drug_entry *entries;
    size_t nentries;
} gh_drug_card_gene;

typedef struct {
    const char *gene;
    const char *description;
} gh_good_news;

typedef struct {
    const char *specialist;
    const char *reason;
    const char *urgency;
} gh_referral;

typedef struct {
    const char *gene;
    const char *status;
    int magnitude;
    const gh_clinical_context *context;
    const char *pathways[8];
    size_t npathways;
} gh_clinical_insight;

typedef struct {
    gh_priority *priorities;
    size_t npriorities;

    gh_drug_card_gene *drug_card;
    size_t ndrug_card;

    gh_monitoring_item *monitoring_schedule;
    size_t nmonitoring;

    gh_good_news *good_news;
    size_t ngood_news;

    gh_referral *referrals;
    size_t nreferrals;

    gh_clinical_insight *insights;
    size_t ninsights;
} gh_recommendations;

/* Everything the engine reads. Any pointer may be NULL. */
typedef struct {
    const gh_analysis *analysis;
    const gh_clinvar_result *clinvar;
    const gh_prs_result *prs;
    size_t nprs;
    const gh_epi_result *epistasis;
    size_t nepistasis;
    const gh_star_result *stars;
    size_t nstars;
    const gh_acmg_result *acmg;
} gh_recommendation_inputs;

void gh_generate_recommendations(gh_recommendations *out, gh_arena *arena,
                                 const gh_recommendation_inputs *in);

/* Priority band for a set of signals. Exposed for testing. */
const char *gh_compute_priority(size_t gene_signals, bool pathogenic,
                                bool prs_high, bool prs_elevated,
                                bool high_magnitude, bool epistasis);

#endif /* GH_RECOMMENDATIONS_H */

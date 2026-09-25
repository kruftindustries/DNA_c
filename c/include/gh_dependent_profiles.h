/* Health profiles that read upstream results as well as the genome.
 *
 * Mirrors sleep_profile.py, nutrigenomics.py, mental_health.py and
 * longevity.py. Unlike the six genome-only profiles these also consume
 * lifestyle findings, APOE, star alleles or polygenic scores, so their
 * tables carry gene/status rules rather than per-SNP outcomes.
 */
#ifndef GH_DEPENDENT_PROFILES_H
#define GH_DEPENDENT_PROFILES_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_analysis.h"
#include "gh_genome.h"
#include "gh_mem.h"
#include "gh_prs.h"
#include "gh_scorers.h"

#define GH_DP_MAX_RECOMMENDATIONS 12
#define GH_DP_MAX_ITEMS 16

/* ------------------------------------------------------------------ */
/* Sleep                                                               */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *rsid;
    const char *gene;
    const char *trait;
    char evening_allele;
    double weight;
    const char *description;
} gh_chronotype_snp;

extern const gh_chronotype_snp GH_CHRONOTYPE_SNPS[];
extern const size_t GH_CHRONOTYPE_SNP_COUNT;

typedef struct {
    const char *rsid;
    const char *gene;
    const char *trait;
    const char *genotype;
    int evening_allele_copies;
    const char *description;
} gh_sleep_marker;

typedef struct {
    const char *chronotype;
    double chronotype_score;          /* 0-100, higher is more evening */
    const char *optimal_sleep_window;
    const char *caffeine_cutoff;
    const char *peak_alertness;
    bool caffeine_sensitive;
    const char *deep_sleep_note;      /* "" when it does not apply */
    const char *confidence;
    size_t markers_found;

    gh_sleep_marker markers[GH_DP_MAX_ITEMS];
    size_t nmarkers;
    const char *recommendations[GH_DP_MAX_RECOMMENDATIONS];
    size_t nrecommendations;
} gh_sleep_result;

void gh_profile_sleep(gh_sleep_result *out, gh_arena *arena,
                      const gh_genome *g, const gh_analysis *analysis);

/* ------------------------------------------------------------------ */
/* Nutrigenomics                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *status;
    double weight;
} gh_status_weight;

typedef struct {
    const char *gene;
    const char *impact;
    const gh_status_weight *statuses;
    size_t nstatuses;
} gh_nutrient_gene;

typedef struct {
    const char *id;
    const char *name;
    const char *supplement_form;
    const char *dose_range;
    const char *testing;
    const gh_nutrient_gene *genes;
    size_t ngenes;
    const char *food_sources;       /* one comma-separated string, as in Python */
} gh_nutrient_profile;

extern const gh_nutrient_profile GH_NUTRIENT_PROFILES[];
extern const size_t GH_NUTRIENT_PROFILE_COUNT;

typedef struct {
    const char *gene;
    const char *status;
    const char *impact;
    double severity;
} gh_nutrient_impact;

typedef struct {
    const gh_nutrient_profile *profile;
    const char *need_level;           /* normal / low / moderate / high / caution_excess */
    double severity;                  /* rounded to two decimals */
    const char *recommendation;
    gh_nutrient_impact impacts[GH_DP_MAX_ITEMS];
    size_t nimpacts;
} gh_nutrient_need;

typedef struct {
    gh_nutrient_need *needs;          /* sorted by descending |severity| */
    size_t nneeds;
    const char *summary;
} gh_nutrigenomics_result;

void gh_profile_nutrigenomics(gh_nutrigenomics_result *out, gh_arena *arena,
                              const gh_analysis *analysis);

/* ------------------------------------------------------------------ */
/* Mental health                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *rsid;
    const char *gene;
    const char *domain;
    char risk_allele;
    const char *description;
    const char *effect;
    const char *reference;
} gh_mental_snp;

extern const gh_mental_snp GH_MENTAL_SNPS[];
extern const size_t GH_MENTAL_SNP_COUNT;

typedef struct {
    const char *gene;
    const char *impact;
    const char *const *risk_statuses;
    size_t nrisk_statuses;
} gh_mental_gene;

typedef struct {
    const char *name;
    const gh_mental_gene *genes;
    size_t ngenes;
} gh_mental_domain;

extern const gh_mental_domain GH_MENTAL_DOMAINS[];
extern const size_t GH_MENTAL_DOMAIN_COUNT;

typedef struct {
    const gh_mental_snp *snp;
    const char *genotype;
    int risk_copies;
} gh_mental_marker;

typedef struct {
    const char *name;
    const char *risk_level;           /* low / moderate / elevated */
    int risk_score;                   /* percentage */
    const char *signals[5];
    size_t nsignals;
} gh_mental_domain_result;

typedef struct {
    gh_mental_domain_result domains[GH_DP_MAX_ITEMS];
    size_t ndomains;
    gh_mental_marker markers[GH_DP_MAX_ITEMS];
    size_t nmarkers;
    const char *risk_factors[8];
    size_t nrisk_factors;
    const char *resilience_factors[GH_DP_MAX_ITEMS];
    size_t nresilience_factors;
    const char *treatment_notes[GH_DP_MAX_ITEMS];
    size_t ntreatment_notes;
    const char *recommendations[GH_DP_MAX_RECOMMENDATIONS];
    size_t nrecommendations;
    const char *summary;
} gh_mental_result;

void gh_profile_mental_health(gh_mental_result *out, gh_arena *arena,
                              const gh_genome *g, const gh_analysis *analysis,
                              const gh_star_result *stars, size_t nstars);

/* ------------------------------------------------------------------ */
/* Longevity                                                           */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *rsid;
    const char *gene;
    const char *name;
    char longevity_allele;
    const char *reference;
} gh_longevity_snp;

extern const gh_longevity_snp GH_LONGEVITY_SNPS[];
extern const size_t GH_LONGEVITY_SNP_COUNT;

typedef struct {
    const char *name;
    const char *const *genes;
    size_t ngenes;
    const char *const *protective_statuses;
    size_t nprotective_statuses;
    const char *const *risk_statuses;
    size_t nrisk_statuses;
} gh_healthspan_domain;

extern const gh_healthspan_domain GH_HEALTHSPAN_DOMAINS[];
extern const size_t GH_HEALTHSPAN_DOMAIN_COUNT;

typedef struct {
    const gh_longevity_snp *snp;
    const char *genotype;
    int copies;
    const char *status;               /* homozygous_protective / heterozygous / absent */
} gh_longevity_allele;

typedef struct {
    const char *name;
    int score;                        /* clamped to 10-90 */
    size_t genes_found;
    const char *rating;               /* good / average / attention */
} gh_healthspan_result;

typedef struct {
    const char *intervention;
    const char *genetic_support;
    const char *why;
} gh_intervention;

typedef struct {
    double longevity_score;           /* 0-100, one decimal */
    size_t alleles_checked;
    gh_longevity_allele alleles[GH_DP_MAX_ITEMS];
    size_t nalleles;
    gh_healthspan_result domains[GH_DP_MAX_ITEMS];
    size_t ndomains;
    const char *top_risks[8];
    size_t ntop_risks;
    const char *top_protective[8];
    size_t ntop_protective;
    gh_intervention interventions[GH_DP_MAX_ITEMS];
    size_t ninterventions;
    const char *summary;
} gh_longevity_result;

void gh_profile_longevity(gh_longevity_result *out, gh_arena *arena,
                          const gh_genome *g, const gh_analysis *analysis,
                          const gh_apoe_result *apoe,
                          const gh_prs_result *prs, size_t nprs);

#endif /* GH_DEPENDENT_PROFILES_H */

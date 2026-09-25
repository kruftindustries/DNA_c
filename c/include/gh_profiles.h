/* Health profile modules: histamine, alcohol, pain, thyroid, hormone and eye.
 *
 * Mirrors the six `profile_*` functions in genetic_health that depend only
 * on the genome. Each is built from per-SNP assessors that branch on one
 * allele's copy count; those outcomes are generated from the Python by
 * tools/gen_profile_data.py, while the aggregation is written out here.
 */
#ifndef GH_PROFILES_H
#define GH_PROFILES_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_genome.h"
#include "gh_mem.h"

/* One assessor's verdict for a given copy count.
 *
 * The slots carry whichever module-specific fields that assessor produces.
 * The mapping is fixed per module in tools/gen_profile_data.py:
 *
 *   module     points                    severity  a          b                c
 *   histamine  risk_points               -         -          -                -
 *   alcohol    cancer_modifier           -         speed      flush            -
 *   pain       sensitivity_contribution  -         -          -                -
 *   thyroid    -                         severity  domain     -                -
 *   hormone    -                         -         domain     estrogen_effect  androgen_effect
 *   eye        -                         severity  condition  -                -
 *
 * Unused slots are zero or NULL. */
typedef struct {
    const char *label;          /* "finding" or "label" in the Python */
    const char *description;
    int points;
    double severity;
    const char *category_a;
    const char *category_b;
    const char *category_c;
} gh_profile_outcome;

typedef struct {
    const char *rsid;
    const char *gene;
    char allele;                /* the counted allele */
    /* outcomes[0] = two copies, [1] = one copy, [2] = none */
    gh_profile_outcome outcomes[3];
} gh_profile_marker;

typedef struct {
    const gh_profile_marker *markers;
    size_t count;
} gh_profile_markers;

extern const gh_profile_markers GH_HISTAMINE_MARKERS;
extern const gh_profile_markers GH_ALCOHOL_MARKERS;
extern const gh_profile_markers GH_PAIN_MARKERS;
extern const gh_profile_markers GH_THYROID_MARKERS;
extern const gh_profile_markers GH_HORMONE_MARKERS;
extern const gh_profile_markers GH_EYE_MARKERS;

/* An assessor's result for a particular genome. */
typedef struct {
    bool present;               /* false when the SNP is not genotyped */
    const char *rsid;
    const char *gene;
    const char *genotype;
    const gh_profile_outcome *outcome;
} gh_profile_hit;

#define GH_PROFILE_MAX_MARKERS 5
#define GH_PROFILE_MAX_RECOMMENDATIONS 8

/* Fields shared by every profile result. */
typedef struct {
    gh_profile_hit hits[GH_PROFILE_MAX_MARKERS];
    size_t nhits;               /* markers actually genotyped */
    size_t markers_tested;
    const char *summary;
    const char *recommendations[GH_PROFILE_MAX_RECOMMENDATIONS];
    size_t nrecommendations;
} gh_profile_base;

/* ------------------------------------------------------------------ */

typedef struct {
    gh_profile_base base;
    const char *risk_level;     /* unknown / low / moderate / elevated */
    int total_risk;
    const char *const *foods_to_watch;
    size_t nfoods;
} gh_histamine_result;

void gh_profile_histamine(gh_histamine_result *out, gh_arena *arena,
                          const gh_genome *g);

typedef struct {
    gh_profile_base base;
    const char *metabolism_speed;   /* slow / normal / fast */
    const char *flush_risk;         /* unknown / none / mild / severe */
    const char *cancer_risk;        /* average / elevated / high */
} gh_alcohol_result;

void gh_profile_alcohol(gh_alcohol_result *out, gh_arena *arena,
                        const gh_genome *g);

typedef struct {
    gh_profile_base base;
    int sensitivity_score;          /* 0-100, higher means more sensitive */
} gh_pain_result;

void gh_profile_pain(gh_pain_result *out, gh_arena *arena, const gh_genome *g);

/* Thyroid and eye both group their markers into named domains and classify
 * each from the severities that landed in it. */
typedef struct {
    const char *name;
    const char *level;
} gh_profile_domain;

#define GH_PROFILE_MAX_DOMAINS 3

typedef struct {
    gh_profile_base base;
    gh_profile_domain domains[GH_PROFILE_MAX_DOMAINS];
    size_t ndomains;
} gh_thyroid_result;

void gh_profile_thyroid(gh_thyroid_result *out, gh_arena *arena,
                        const gh_genome *g);

typedef struct {
    gh_profile_base base;
    gh_profile_domain conditions[GH_PROFILE_MAX_DOMAINS];
    size_t nconditions;
} gh_eye_result;

void gh_profile_eye(gh_eye_result *out, gh_arena *arena, const gh_genome *g);

typedef struct {
    gh_profile_base base;
    const char *estrogen_level;
    const char *androgen_level;
    const char *overall;
} gh_hormone_result;

void gh_profile_hormone(gh_hormone_result *out, gh_arena *arena,
                        const gh_genome *g);

/* Level for a named domain, or NULL when the profile has no such domain. */
const char *gh_profile_domain_level(const gh_profile_domain *domains,
                                    size_t n, const char *name);

#endif /* GH_PROFILES_H */

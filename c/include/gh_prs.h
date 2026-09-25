/* Polygenic risk scores.
 *
 * Mirrors genetic_health/prs.py. Each model sums published log-odds effect
 * sizes over the risk alleles present, standardises against the European
 * allele frequencies the model was built on, and converts to a percentile.
 */
#ifndef GH_PRS_H
#define GH_PRS_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_ancestry.h"
#include "gh_genome.h"
#include "gh_mem.h"

typedef struct {
    const char *rsid;
    const char *gene;
    char risk_allele;
    double log_or;
    double eur_freq;
} gh_prs_snp;

typedef struct {
    const char *id;          /* "type2_diabetes" */
    const char *name;
    const char *reference;
    const gh_prs_snp *snps;
    size_t nsnps;
} gh_prs_model;

extern const gh_prs_model GH_PRS_MODELS[];
extern const size_t GH_PRS_MODEL_COUNT;

typedef struct {
    const char *rsid;
    const char *gene;
    char risk_allele;
    int copies;
    double log_or;
    double contribution;
} gh_prs_contribution;

/* Python reports at most ten contributing SNPs per condition. */
#define GH_PRS_MAX_CONTRIBUTIONS 10

typedef struct {
    const char *id;
    const char *name;
    const char *reference;

    double raw_score;      /* all values already rounded as Python reports them */
    double z_score;
    double percentile;
    double ci_95_lower;
    double ci_95_upper;
    const char *risk_category;   /* low / average / elevated / high */

    size_t snps_found;
    size_t snps_total;           /* unique rsIDs in the model */

    bool ancestry_applicable;
    const char *ancestry_warning;  /* "" when none */

    gh_prs_contribution contributing[GH_PRS_MAX_CONTRIBUTIONS];
    size_t ncontributing;
} gh_prs_result;

/* Score every model. `out` must have room for GH_PRS_MODEL_COUNT entries.
 * `ancestry` may be NULL, in which case European ancestry is assumed and no
 * adjustment is applied. */
void gh_calculate_prs(gh_prs_result *out, gh_arena *arena, const gh_genome *g,
                      const gh_ancestry_result *ancestry);

/* Risk band for a percentile: low (<20), average (<80), elevated (<95),
 * otherwise high. */
const char *gh_prs_categorize(double percentile);

/* Percentile from a Z-score via the error function, clamped to [0.1, 99.9]. */
double gh_prs_z_to_percentile(double z);

/* Round half-to-even at `digits` decimal places, matching Python's round().
 * Exposed because several reported PRS fields depend on it. */
double gh_round(double value, int digits);

#endif /* GH_PRS_H */

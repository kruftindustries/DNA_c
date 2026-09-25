/* Gene-gene interaction (epistasis) evaluation.
 *
 * Mirrors genetic_health/epistasis.py. Each model lists conditions that
 * require particular statuses across several genes; a match scores a
 * severity from the matched statuses and their magnitudes, and the
 * interaction's strength is the geometric mean of the per-gene severities.
 */
#ifndef GH_EPISTASIS_H
#define GH_EPISTASIS_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_analysis.h"
#include "gh_mem.h"

typedef struct {
    const char *status;
    double weight;
} gh_status_severity;

extern const gh_status_severity GH_STATUS_SEVERITY[];
extern const size_t GH_STATUS_SEVERITY_COUNT;

typedef struct {
    const char *gene;
    const char *const *statuses;   /* any one of these qualifies */
    size_t nstatuses;
} gh_epi_requirement;

typedef struct {
    const gh_epi_requirement *required;
    size_t nrequired;
    const char *effect;
    const char *risk_level;
    const char *mechanism;
    const char *const *actions;
    size_t nactions;
} gh_epi_condition;

typedef struct {
    const char *id;
    const char *name;
    const gh_epi_condition *conditions;
    size_t nconditions;
} gh_epi_model;

extern const gh_epi_model GH_EPISTASIS_MODELS[];
extern const size_t GH_EPISTASIS_MODEL_COUNT;

#define GH_EPI_MAX_GENES 4

typedef struct {
    const char *gene;
    const char **statuses;   /* the matched statuses, sorted */
    size_t nstatuses;
} gh_epi_match;

typedef struct {
    const char *id;
    const char *name;
    const char *effect;
    const char *risk_level;    /* after severity adjustment */
    const char *mechanism;
    const char *const *actions;
    size_t nactions;
    double severity_score;     /* rounded to two decimals */

    gh_epi_match genes[GH_EPI_MAX_GENES];
    size_t ngenes;
} gh_epi_result;

/* Evaluate every model against the lifestyle findings. Results are sorted
 * high risk first, then by descending severity. `*count` receives how many
 * interactions matched. */
gh_epi_result *gh_evaluate_epistasis(gh_arena *arena, const gh_analysis *a,
                                     size_t *count);

/* Severity weight for a status, defaulting to 0.5 for anything unlisted. */
double gh_status_severity_for(const char *status);

#endif /* GH_EPISTASIS_H */

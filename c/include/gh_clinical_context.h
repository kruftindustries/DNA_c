/* Clinical context for a gene/status combination.
 *
 * Mirrors genetic_health/clinical_context.py. The table is keyed on both the
 * gene and the status, so the same gene can carry different guidance
 * depending on how badly its function is affected.
 */
#ifndef GH_CLINICAL_CONTEXT_H
#define GH_CLINICAL_CONTEXT_H

#include <stddef.h>

typedef struct {
    const char *gene;
    const char *status;
    const char *mechanism;
    const char *const *implications;
    size_t nimplications;
    const char *const *actions;
    size_t nactions;
    const char *const *interactions;
    size_t ninteractions;
} gh_clinical_context;

extern const gh_clinical_context GH_CLINICAL_CONTEXT[];
extern const size_t GH_CLINICAL_CONTEXT_COUNT;

/* Context for a gene/status pair, or NULL when the table has no entry. */
const gh_clinical_context *gh_clinical_context_for(const char *gene,
                                                   const char *status);

#endif /* GH_CLINICAL_CONTEXT_H */

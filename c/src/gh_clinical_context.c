#include "gh_clinical_context.h"

#include <string.h>

const gh_clinical_context *gh_clinical_context_for(const char *gene,
                                                   const char *status)
{
    if (!gene || !status)
        return NULL;
    for (size_t i = 0; i < GH_CLINICAL_CONTEXT_COUNT; i++)
        if (strcmp(GH_CLINICAL_CONTEXT[i].gene, gene) == 0 &&
            strcmp(GH_CLINICAL_CONTEXT[i].status, status) == 0)
            return &GH_CLINICAL_CONTEXT[i];
    return NULL;
}

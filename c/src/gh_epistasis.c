#include "gh_epistasis.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "gh_prs.h"   /* gh_round */

double gh_status_severity_for(const char *status)
{
    for (size_t i = 0; i < GH_STATUS_SEVERITY_COUNT; i++)
        if (strcmp(GH_STATUS_SEVERITY[i].status, status) == 0)
            return GH_STATUS_SEVERITY[i].weight;
    /* Unlisted statuses count as moderate rather than zero, so an
     * interaction is not silently nullified by one unknown status. */
    return 0.5;
}

static int compare_strp(const void *a, const void *b)
{
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

/* Highest risk first, then strongest interaction. Python sorts on
 * (risk_order, -severity) with a stable sort, so ties keep model order. */
typedef struct {
    gh_epi_result result;
    size_t order;
} ranked_result;

static int risk_order(const char *risk)
{
    if (strcmp(risk, "high") == 0)
        return 0;
    if (strcmp(risk, "moderate") == 0)
        return 1;
    if (strcmp(risk, "low") == 0)
        return 2;
    return 3;
}

static int compare_result(const void *a, const void *b)
{
    const ranked_result *x = a, *y = b;
    int rx = risk_order(x->result.risk_level);
    int ry = risk_order(y->result.risk_level);
    if (rx != ry)
        return rx < ry ? -1 : 1;
    if (x->result.severity_score != y->result.severity_score)
        return x->result.severity_score > y->result.severity_score ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

gh_epi_result *gh_evaluate_epistasis(gh_arena *arena, const gh_analysis *a,
                                     size_t *count)
{
    *count = 0;

    /* Every condition of every model could match. */
    size_t capacity = 0;
    for (size_t m = 0; m < GH_EPISTASIS_MODEL_COUNT; m++)
        capacity += GH_EPISTASIS_MODELS[m].nconditions;
    ranked_result *ranked = gh_calloc(arena, capacity, sizeof(*ranked));
    size_t n = 0;

    for (size_t m = 0; m < GH_EPISTASIS_MODEL_COUNT; m++) {
        const gh_epi_model *model = &GH_EPISTASIS_MODELS[m];

        for (size_t c = 0; c < model->nconditions; c++) {
            const gh_epi_condition *cond = &model->conditions[c];

            gh_epi_match matches[GH_EPI_MAX_GENES];
            double severities[GH_EPI_MAX_GENES];
            size_t ngenes = 0;
            bool matched = true;

            for (size_t r = 0; r < cond->nrequired && matched; r++) {
                const gh_epi_requirement *req = &cond->required[r];

                /* Collect the statuses this genome shows for the gene that
                 * also qualify, and the gene's highest magnitude. */
                const char **overlap =
                    gh_calloc(arena, req->nstatuses, sizeof(*overlap));
                size_t noverlap = 0;
                int max_magnitude = -1;
                bool gene_present = false;

                for (size_t f = 0; f < a->nfindings; f++) {
                    const gh_finding *finding = &a->findings[f];
                    if (!finding->gene || !finding->status ||
                        !*finding->gene || !*finding->status)
                        continue;
                    if (strcmp(finding->gene, req->gene) != 0)
                        continue;

                    gene_present = true;
                    if (finding->magnitude > max_magnitude)
                        max_magnitude = finding->magnitude;

                    for (size_t s = 0; s < req->nstatuses; s++) {
                        if (strcmp(finding->status, req->statuses[s]) != 0)
                            continue;
                        /* The Python intersects sets, so each status counts
                         * once however many findings carry it. */
                        bool already = false;
                        for (size_t k = 0; k < noverlap && !already; k++)
                            if (strcmp(overlap[k], req->statuses[s]) == 0)
                                already = true;
                        if (!already)
                            overlap[noverlap++] = req->statuses[s];
                    }
                }

                if (!gene_present || noverlap == 0) {
                    matched = false;
                    break;
                }

                qsort(overlap, noverlap, sizeof(*overlap), compare_strp);

                double severity = 0.0;
                for (size_t k = 0; k < noverlap; k++) {
                    double w = gh_status_severity_for(overlap[k]);
                    if (w > severity)
                        severity = w;
                }

                /* Magnitude runs 0-6; a gene with no finding would default
                 * to 3, though that cannot happen here since a match
                 * requires at least one. */
                double magnitude = max_magnitude >= 0 ? max_magnitude : 3;
                double mag_factor = magnitude / 6.0;
                if (mag_factor > 1.0)
                    mag_factor = 1.0;

                matches[ngenes].gene = req->gene;
                matches[ngenes].statuses = overlap;
                matches[ngenes].nstatuses = noverlap;
                severities[ngenes] = severity * (0.5 + 0.5 * mag_factor);
                ngenes++;
            }

            if (!matched || ngenes == 0)
                continue;

            /* Geometric mean, so one weak gene pulls the whole interaction
             * down rather than being averaged away. */
            double product = 1.0;
            for (size_t i = 0; i < ngenes; i++)
                product *= severities[i];
            double score = pow(product, 1.0 / (double)ngenes);

            const char *risk = cond->risk_level;
            if (score < 0.3 && strcmp(risk, "high") == 0)
                risk = "moderate";
            else if (score >= 0.7 && strcmp(risk, "moderate") == 0)
                risk = "high";

            gh_epi_result *out = &ranked[n].result;
            out->id = model->id;
            out->name = model->name;
            out->effect = cond->effect;
            out->risk_level = risk;
            out->mechanism = cond->mechanism;
            out->actions = cond->actions;
            out->nactions = cond->nactions;
            out->severity_score = gh_round(score, 2);
            out->ngenes = ngenes;
            for (size_t i = 0; i < ngenes; i++)
                out->genes[i] = matches[i];
            ranked[n].order = n;
            n++;
        }
    }

    qsort(ranked, n, sizeof(*ranked), compare_result);

    gh_epi_result *results = gh_calloc(arena, n ? n : 1, sizeof(*results));
    for (size_t i = 0; i < n; i++)
        results[i] = ranked[i].result;

    *count = n;
    return results;
}

/* Ancestry estimation from Ancestry-Informative Markers.
 *
 * Mirrors genetic_health/ancestry.py: per-population log-likelihoods over
 * ~54 AIMs, softmax-normalised to proportions, plus a sub-population
 * estimate within the top superpopulation.
 */
#ifndef GH_ANCESTRY_H
#define GH_ANCESTRY_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_genome.h"
#include "gh_mem.h"

/* Superpopulations, in the order Python iterates them. A tie in the
 * top-ancestry pick resolves to the earliest, matching max() over an
 * insertion-ordered dict. */
extern const char *const GH_POPULATIONS[];
extern const char *const GH_POPULATION_LABELS[];
extern const size_t GH_POPULATION_COUNT;

#define GH_MAX_POPULATIONS 5
#define GH_MAX_SUB_POPULATIONS 4

typedef struct {
    const char *rsid;
    const char *gene;
    const char *description;
    char allele;                          /* the informative allele */
    double frequencies[GH_MAX_POPULATIONS];
} gh_aim;

extern const gh_aim GH_AIMS[];
extern const size_t GH_AIM_COUNT;

typedef struct {
    const char *rsid;
    char allele;
    double frequencies[GH_MAX_SUB_POPULATIONS];
} gh_sub_marker;

typedef struct {
    const char *population;          /* "EUR" */
    const char *const *codes;        /* "EUR_N", ... */
    const char *const *labels;       /* "Northern European", ... */
    size_t nsubs;
    const gh_sub_marker *markers;
    size_t nmarkers;
} gh_sub_population;

extern const gh_sub_population GH_SUB_POPULATIONS[];
extern const size_t GH_SUB_POPULATION_COUNT;

typedef struct {
    const char *rsid;
    const char *gene;
    const char *description;
    const char *genotype;
    int allele_count;
} gh_aim_detail;

typedef struct {
    bool present;                 /* false when fewer than 3 markers were used */
    const char *top_label;        /* "Northern European" */
    const char *confidence;
    size_t markers_used;
    const char *labels[GH_MAX_SUB_POPULATIONS];
    double proportions[GH_MAX_SUB_POPULATIONS];
    size_t nsubs;
} gh_sub_ancestry;

typedef struct {
    double proportions[GH_MAX_POPULATIONS];   /* parallel to GH_POPULATIONS */
    size_t markers_found;
    const char *confidence;       /* high / moderate / low / none */
    const char *top_ancestry;     /* label, or "Unknown" */
    int top_index;                /* -1 when unknown */

    gh_aim_detail *details;
    size_t ndetails;

    gh_sub_ancestry sub;
} gh_ancestry_result;

void gh_estimate_ancestry(gh_ancestry_result *out, gh_arena *arena,
                          const gh_genome *g);

#endif /* GH_ANCESTRY_H */

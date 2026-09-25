/* Lifestyle/health SNP analysis and drug-gene matching.
 *
 * Mirrors genetic_health/analysis.py:analyze_snps.
 */
#ifndef GH_ANALYSIS_H
#define GH_ANALYSIS_H

#include <stddef.h>

#include "gh_genome.h"
#include "gh_mem.h"
#include "gh_pharmgkb.h"
#include "gh_snpdb.h"

typedef struct {
    const char *rsid;
    const char *gene;
    const char *category;
    const char *genotype;
    const char *status;
    const char *description;
    const char *note;        /* "" when the SNP carries no note */
    const double *freq;      /* NULL when unknown */
    int magnitude;
} gh_finding;

typedef struct {
    const char *rsid;
    const char *gene;
    const char *drugs;
    const char *genotype;
    const char *annotation;
    const char *level;
    const char *category;
} gh_drug_finding;

typedef struct {
    gh_finding *findings;
    size_t nfindings;

    gh_drug_finding *drug_findings;
    size_t ndrug_findings;

    size_t total_snps;      /* variants present in the genome */
    size_t analyzed_snps;   /* database SNPs matched */
    size_t high_impact;     /* magnitude >= 3 */
    size_t moderate_impact; /* magnitude == 2 */
    size_t low_impact;      /* magnitude == 1 */
} gh_analysis;

/* Run the analysis. `pgx` may be NULL when the ClinPGx files were not loaded;
 * drug-gene findings are then empty and the rest still runs. */
void gh_analyze(gh_analysis *out, gh_arena *arena, const gh_genome *genome,
                const gh_pharmgkb *pgx);

/* True for the evidence levels the report surfaces: 1A, 1B, 2A, 2B. */
bool gh_pgx_level_reportable(const char *level);

#endif /* GH_ANALYSIS_H */

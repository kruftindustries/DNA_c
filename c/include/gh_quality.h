/* Data quality metrics.
 *
 * Mirrors genetic_health/quality_metrics.py:compute_quality_metrics.
 */
#ifndef GH_QUALITY_H
#define GH_QUALITY_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_genome.h"
#include "gh_mem.h"

typedef struct {
    const char *name;   /* normalised: upper-cased, "CHR" prefix removed */
    size_t count;
} gh_chrom_count;

typedef struct {
    size_t total_snps;
    size_t no_call_count;
    double call_rate;

    /* Sorted for stable output: numeric chromosomes ascending, then the
     * non-numeric ones (MT, X, Y, ...) alphabetically. */
    gh_chrom_count *chromosomes;
    size_t nchromosomes;

    bool has_mt;
    bool has_y;
    size_t mt_snp_count;
    size_t autosomal_count;
    double het_rate;
} gh_quality;

/* Compute metrics from a loaded genome. `genome_path` may be NULL, in which
 * case no-calls are not counted and the call rate is 1.0 (matching the
 * Python, which then divides total by total). */
void gh_quality_compute(gh_quality *q, gh_arena *arena, const gh_genome *g,
                        const char *genome_path);

/* Count for a chromosome, or 0 when absent. */
size_t gh_quality_chrom_count(const gh_quality *q, const char *chrom);

#endif /* GH_QUALITY_H */

/* ClinPGx (formerly PharmGKB) clinical annotations.
 *
 * Mirrors genetic_health/loading.py:load_pharmgkb. The two tables are joined
 * on the annotation ID and keyed by rsID; only annotations whose
 * Variant/Haplotypes field is a plain rsID are kept, matching the reference.
 *
 * ClinPGx renamed the ID column from "Clinical Annotation ID" to "Summary
 * Annotation ID" in 2025, so both spellings are accepted.
 */
#ifndef GH_PHARMGKB_H
#define GH_PHARMGKB_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_map.h"
#include "gh_mem.h"

typedef struct {
    const char *rsid;
    const char *gene;
    const char *drugs;
    const char *phenotype;
    const char *level;      /* level of evidence: 1A, 1B, 2A, 2B, 3, 4 */
    const char *category;   /* phenotype category */
    gh_map genotypes;       /* genotype/allele -> annotation text */
} gh_pgx_entry;

typedef struct {
    gh_map by_rsid;         /* rsID -> gh_pgx_entry* */
    /* Entries in the order their rsID was first seen, which is what the
     * Python's dict preserves. Consumers that sort with a stable sort need
     * this rather than the hash map's arbitrary order. */
    gh_pgx_entry **ordered;
    size_t nordered;
    size_t annotations;     /* rsID-bearing rows read from the annotation table */
    size_t allele_rows;     /* rows read from the allele table */
    gh_arena *arena;
} gh_pharmgkb;

/* Load both tables. Returns false if either file is missing or unreadable;
 * the caller then skips drug-gene analysis, as the Python does. */
bool gh_pharmgkb_load(gh_pharmgkb *p, gh_arena *arena,
                      const char *annotations_path, const char *alleles_path);

const gh_pgx_entry *gh_pharmgkb_find(const gh_pharmgkb *p, const char *rsid);

/* The annotation-ID column under either of its names. */
extern const char *const GH_PGX_ID_COLUMNS[];

#endif /* GH_PHARMGKB_H */

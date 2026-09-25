/* The curated lifestyle/health SNP database.
 *
 * The table itself is generated from genetic_health/snp_database.py by
 * tools/gen_snpdb.py, so the C port and the Python reference cannot drift.
 * Regenerate with `make gen` (or tools/gen_snpdb.py) after editing the
 * Python source; the generated file is checked in so a plain `make` needs
 * no Python.
 */
#ifndef GH_SNPDB_H
#define GH_SNPDB_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_map.h"
#include "gh_mem.h"

/* Population order used by every `freq` array in the generated table. */
enum { GH_POP_EUR, GH_POP_AFR, GH_POP_EAS, GH_POP_SAS, GH_POP_AMR, GH_POP_COUNT };

extern const char *const GH_POP_NAMES[GH_POP_COUNT];

typedef struct {
    const char *genotype;
    const char *status;
    const char *desc;
    int magnitude;   /* 0 = unremarkable, 4 = most clinically significant */
} gh_snp_variant;

typedef struct {
    const char *rsid;
    const char *gene;
    const char *category;
    const char *note;              /* NULL when the entry carries no note */
    const gh_snp_variant *variants;
    size_t nvariants;
    /* Minor (variant) allele frequency per superpopulation, or NULL.
     * These are independent per-population measurements: they do not sum
     * to 1 across populations. */
    const double *freq;
} gh_snp;

extern const gh_snp GH_SNPS[];
extern const size_t GH_SNP_COUNT;

/* Look up a SNP by rsID. Returns NULL when the database has no such entry.
 * Backed by a lazily built index, so repeated lookups are O(1). */
const gh_snp *gh_snpdb_find(const char *rsid);

/* Find the entry for a genotype within a SNP, trying the reversed spelling
 * too ("CT" matches a table entry of "TC"). Returns NULL when unmatched. */
const gh_snp_variant *gh_snp_match(const gh_snp *snp, const char *genotype);

#endif /* GH_SNPDB_H */

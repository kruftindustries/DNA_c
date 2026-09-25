/* Genome loading: 23andMe and AncestryDNA raw exports and WGS-derived files.
 *
 * Mirrors genetic_health/loading.py:load_genome, including its filtering
 * rules — no-calls are dropped silently, anything that is not one or two
 * A/C/G/T bases is counted as skipped, and non-rsID identifiers are resolved
 * through the chrom:pos lookup when one is loaded. The vendor layout is
 * detected per row by gh_genome_parse_row, which the quality metrics share
 * so the two can never disagree about what a no-call is.
 */
#ifndef GH_GENOME_H
#define GH_GENOME_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_map.h"
#include "gh_mem.h"

typedef struct {
    const char *rsid;
    const char *chromosome;
    const char *position;
    const char *genotype;
} gh_variant;

typedef struct {
    gh_map by_rsid;      /* rsID    -> gh_variant* */
    gh_map by_position;  /* chr:pos -> gh_variant* */

    size_t total;        /* variants kept */
    size_t skipped;      /* malformed genotypes */
    size_t no_calls;     /* "--" */
    size_t resolved;     /* positional IDs mapped to an rsID */

    gh_arena *arena;
} gh_genome;

/* Load a genome file. `pos_to_rsid` may be NULL (or empty) when no lookup
 * table is available; positional IDs are then kept as-is.
 * Returns false only if the file cannot be opened. */
bool gh_genome_load(gh_genome *g, gh_arena *arena, const char *path,
                    const gh_map *pos_to_rsid);

/* Genotype for an rsID, or NULL when absent. */
const char *gh_genome_genotype(const gh_genome *g, const char *rsid);

const gh_variant *gh_genome_variant(const gh_genome *g, const char *rsid);

/* Load data/rsid_positions_grch37.json into a chrom:pos -> rsID map.
 * The file is a flat object of {"rs123": {"chrom": "1", "pos": 12345}, ...};
 * this reads that shape directly rather than pulling in a JSON library.
 * Returns false if the file is missing or unreadable. */
bool gh_load_rsid_positions(gh_map *out, gh_arena *arena, const char *path);

/* True if `s` is one or two characters, all of A/C/G/T. */
bool gh_genotype_valid(const char *s);

/* ------------------------------------------------------------------ */
/* Row classification, shared with gh_quality                          */
/* ------------------------------------------------------------------ */

/* Split a line in place on tabs into at most `max` fields, returning how
 * many were found. The line must already have its terminator stripped. */
size_t gh_split_tabs(char *line, char **out, size_t max);

/* Enough fields for either vendor layout. */
#define GH_ROW_MAX_FIELDS 5

typedef enum {
    GH_ROW_OK = 0,     /* a usable genotype */
    GH_ROW_NO_CALL,    /* "--" (23andMe) or a 0 allele (AncestryDNA) */
    GH_ROW_INVALID,    /* genotype is not one or two A/C/G/T bases */
    GH_ROW_SHORT,      /* fewer than four fields */
} gh_row_kind;

typedef struct {
    const char *rsid;
    const char *chromosome;   /* AncestryDNA's 23/24/25/26 already mapped */
    const char *position;
    const char *genotype;     /* NULL unless GH_ROW_OK or GH_ROW_INVALID */
    char joined[3];           /* backing store for a two-column genotype */
} gh_row;

/* Classify one tab-split row. Mirrors loading.py:parse_genome_row:
 *
 *   23andMe / WGS-converted:  rsid chrom pos genotype          (4 fields)
 *   AncestryDNA:              rsid chrom pos allele1 allele2   (5 fields)
 *
 * The AncestryDNA layout is recognised by its two single-character allele
 * columns. Before this existed such a file was read as the 23andMe layout,
 * so allele1 alone passed validation as a haploid call and every
 * heterozygote in the file was silently lost. */
gh_row_kind gh_genome_parse_row(char **parts, size_t nparts, gh_row *out);

#endif /* GH_GENOME_H */

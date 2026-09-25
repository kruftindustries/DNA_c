/* getline is POSIX.1-2008; request it explicitly under -std=c11. */
#define _POSIX_C_SOURCE 200809L

#include "gh_quality.h"
#include "gh_tsv.h"

#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Upper-case `s` and delete every occurrence of "CHR".
 *
 * Python does `str(raw).upper().replace("CHR", "")`, which removes the
 * substring anywhere, not just at the front -- reproduced here so an odd
 * chromosome label normalises the same way in both implementations. */
static char *normalise_chrom(gh_arena *arena, const char *s)
{
    size_t n = strlen(s);
    char *up = gh_alloc(arena, n + 1);
    for (size_t i = 0; i < n; i++)
        up[i] = (char)toupper((unsigned char)s[i]);
    up[n] = '\0';

    char *out = gh_alloc(arena, n + 1);
    size_t w = 0;
    for (size_t i = 0; i < n;) {
        if (n - i >= 3 && strncmp(up + i, "CHR", 3) == 0) {
            i += 3;
            continue;
        }
        out[w++] = up[i++];
    }
    out[w] = '\0';
    return out;
}

/* True when every character is a digit and the string is non-empty --
 * Python's str.isdigit() for the ASCII inputs this sees. */
static bool all_digits(const char *s)
{
    if (!*s)
        return false;
    for (const char *p = s; *p; p++)
        if (!isdigit((unsigned char)*p))
            return false;
    return true;
}

static bool is_autosome(const char *chrom)
{
    if (!all_digits(chrom))
        return false;
    long n = strtol(chrom, NULL, 10);
    return n >= 1 && n <= 22;
}

static int compare_chrom(const void *a, const void *b)
{
    const gh_chrom_count *x = a, *y = b;
    bool xd = all_digits(x->name), yd = all_digits(y->name);

    if (xd && yd) {
        long nx = strtol(x->name, NULL, 10);
        long ny = strtol(y->name, NULL, 10);
        return (nx > ny) - (nx < ny);
    }
    if (xd != yd)
        return xd ? -1 : 1;   /* numeric chromosomes first */
    return strcmp(x->name, y->name);
}

/* Count no-calls straight from the raw file.
 *
 * quality_metrics.py strips each line *before* testing for a leading '#',
 * unlike the genome loader, so a line such as "  # note" is a comment here
 * but not there. The count is recomputed rather than taken from the loader
 * so that difference is preserved. What counts as a no-call, though, comes
 * from the same classifier the loader uses -- "--" in a 23andMe file, a 0
 * allele in an AncestryDNA one -- so the two cannot drift apart. */
static size_t count_no_calls(const char *path)
{
    FILE *f = fopen(path, "rb");
    if (!f)
        return 0;

    size_t no_calls = 0;
    char *line = NULL;
    size_t cap = 0;
    long len;

    while ((len = gh_read_line(&line, &cap, f)) != -1) {
        char *start = line;
        char *end = line + len;
        while (start < end && isspace((unsigned char)*start))
            start++;
        while (end > start && isspace((unsigned char)end[-1]))
            end--;
        *end = '\0';

        if (start == end || *start == '#')
            continue;

        char *parts[GH_ROW_MAX_FIELDS];
        size_t nparts = gh_split_tabs(start, parts, GH_ROW_MAX_FIELDS);
        gh_row row;
        if (gh_genome_parse_row(parts, nparts, &row) == GH_ROW_NO_CALL)
            no_calls++;
    }

    free(line);
    fclose(f);
    return no_calls;
}

void gh_quality_compute(gh_quality *q, gh_arena *arena, const gh_genome *g,
                        const char *genome_path)
{
    memset(q, 0, sizeof(*q));
    q->total_snps = gh_map_len(&g->by_rsid);

    /* chromosome -> index into a growing counts array. */
    gh_map index;
    gh_map_init(&index, arena, 64);

    size_t cap = 64;
    gh_chrom_count *counts = gh_calloc(arena, cap, sizeof(*counts));
    size_t n = 0;

    size_t autosomal_het = 0, autosomal_total = 0;

    size_t iter = 0;
    void *value;
    while (gh_map_next(&g->by_rsid, &iter, NULL, &value)) {
        const gh_variant *v = value;
        if (!v->chromosome || !*v->chromosome)
            continue;

        const char *chrom = normalise_chrom(arena, v->chromosome);
        if (!*chrom)
            continue;

        /* Values are index+1 as a tagged integer, so index 0 is still
         * distinguishable from the map's NULL-for-missing. */
        void *slot = gh_map_get(&index, chrom);
        size_t i;
        if (slot) {
            i = (size_t)(uintptr_t)slot - 1;
        } else {
            if (n == cap) {
                size_t grown_cap = cap * 2;
                gh_chrom_count *grown =
                    gh_calloc(arena, grown_cap, sizeof(*grown));
                memcpy(grown, counts, n * sizeof(*counts));
                counts = grown;
                cap = grown_cap;
            }
            i = n++;
            counts[i].name = chrom;
            counts[i].count = 0;
            gh_map_put(&index, chrom, (void *)(uintptr_t)(i + 1));
        }
        counts[i].count++;

        /* Heterozygosity is measured on autosomes only: MT, and X in males,
         * are haploid and cannot be heterozygous. */
        if (is_autosome(chrom)) {
            autosomal_total++;
            const char *gt = v->genotype;
            if (gt && strlen(gt) == 2 && gt[0] != gt[1])
                autosomal_het++;
        }
    }

    qsort(counts, n, sizeof(*counts), compare_chrom);
    q->chromosomes = counts;
    q->nchromosomes = n;

    for (size_t i = 0; i < n; i++) {
        if (strcmp(counts[i].name, "MT") == 0) {
            q->has_mt = true;
            q->mt_snp_count = counts[i].count;
        } else if (strcmp(counts[i].name, "Y") == 0) {
            q->has_y = true;
        }
        if (is_autosome(counts[i].name))
            q->autosomal_count += counts[i].count;
    }

    q->het_rate = autosomal_total
        ? (double)autosomal_het / (double)autosomal_total : 0.0;

    q->no_call_count = genome_path ? count_no_calls(genome_path) : 0;

    size_t denom = q->total_snps + q->no_call_count;
    q->call_rate = denom ? (double)q->total_snps / (double)denom : 0.0;
}

size_t gh_quality_chrom_count(const gh_quality *q, const char *chrom)
{
    for (size_t i = 0; i < q->nchromosomes; i++)
        if (strcmp(q->chromosomes[i].name, chrom) == 0)
            return q->chromosomes[i].count;
    return 0;
}

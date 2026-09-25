#include "gh_snpdb.h"

#include <stdlib.h>
#include <string.h>

/* Index over the generated table, built once on first lookup.
 *
 * The table is static and read-only, so the index can be too. It uses its own
 * arena rather than a caller-supplied one because its lifetime is the whole
 * process; gh_snpdb_release exists mainly so leak checkers stay quiet. */
static gh_arena *index_arena;
static gh_map index_map;
static bool index_built;

static void build_index(void)
{
    index_arena = gh_arena_new(64 * 1024);
    if (!index_arena)
        gh_oom(64 * 1024);
    gh_map_init(&index_map, index_arena, GH_SNP_COUNT);
    for (size_t i = 0; i < GH_SNP_COUNT; i++)
        gh_map_put(&index_map, GH_SNPS[i].rsid, (void *)&GH_SNPS[i]);
    index_built = true;
}

const gh_snp *gh_snpdb_find(const char *rsid)
{
    if (!rsid)
        return NULL;
    if (!index_built)
        build_index();
    return gh_map_get(&index_map, rsid);
}

const gh_snp_variant *gh_snp_match(const gh_snp *snp, const char *genotype)
{
    if (!snp || !genotype || !*genotype)
        return NULL;

    for (size_t i = 0; i < snp->nvariants; i++)
        if (strcmp(snp->variants[i].genotype, genotype) == 0)
            return &snp->variants[i];

    /* Genotype calls are unordered, so "CT" and "TC" are the same call. The
     * table spells each one only once. */
    size_t n = strlen(genotype);
    if (n == 2) {
        char swapped[3] = {genotype[1], genotype[0], '\0'};
        for (size_t i = 0; i < snp->nvariants; i++)
            if (strcmp(snp->variants[i].genotype, swapped) == 0)
                return &snp->variants[i];
    }
    return NULL;
}

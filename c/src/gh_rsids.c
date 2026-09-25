#include "gh_rsids.h"

#include <stdlib.h>
#include <string.h>

#include "gh_ancestry.h"
#include "gh_dependent_profiles.h"
#include "gh_map.h"
#include "gh_prs.h"
#include "gh_profiles.h"
#include "gh_scorers.h"
#include "gh_snpdb.h"
#include "gh_traits.h"

typedef struct {
    gh_map seen;
    const char **items;
    size_t n, cap;
    gh_arena *arena;
} rsid_set;

static void add(rsid_set *s, const char *rsid)
{
    if (!rsid || strncmp(rsid, "rs", 2) != 0 || gh_map_get(&s->seen, rsid))
        return;
    if (s->n == s->cap) {
        size_t grown = s->cap ? s->cap * 2 : 1024;
        const char **items = gh_calloc(s->arena, grown, sizeof *items);
        if (s->n)
            memcpy(items, s->items, s->n * sizeof *items);
        s->items = items;
        s->cap = grown;
    }
    s->items[s->n++] = rsid;
    gh_map_put(&s->seen, rsid, (void *)rsid);
}

static void add_markers(rsid_set *s, const gh_profile_markers *m)
{
    for (size_t i = 0; i < m->count; i++)
        add(s, m->markers[i].rsid);
}

static int compare_strp(const void *a, const void *b)
{
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

size_t gh_all_rsids(gh_arena *arena, const char ***out)
{
    rsid_set s;
    memset(&s, 0, sizeof s);
    s.arena = arena;
    gh_map_init(&s.seen, arena, 4096);

    for (size_t i = 0; i < GH_SNP_COUNT; i++)
        add(&s, GH_SNPS[i].rsid);
    for (size_t i = 0; i < GH_AIM_COUNT; i++)
        add(&s, GH_AIMS[i].rsid);
    for (size_t p = 0; p < GH_SUB_POPULATION_COUNT; p++)
        for (size_t i = 0; i < GH_SUB_POPULATIONS[p].nmarkers; i++)
            add(&s, GH_SUB_POPULATIONS[p].markers[i].rsid);
    for (size_t m = 0; m < GH_PRS_MODEL_COUNT; m++)
        for (size_t i = 0; i < GH_PRS_MODELS[m].nsnps; i++)
            add(&s, GH_PRS_MODELS[m].snps[i].rsid);

    /* Blood type and APOE read fixed SNPs rather than a table. */
    add(&s, "rs505922"); add(&s, "rs8176746"); add(&s, "rs590787");
    add(&s, "rs429358"); add(&s, "rs7412");

    for (size_t i = 0; i < GH_MT_TREE_COUNT; i++)
        add(&s, GH_MT_TREE[i].rsid);
    for (size_t g = 0; g < GH_STAR_GENE_COUNT; g++)
        for (size_t i = 0; i < GH_STAR_GENES[g].nsnps; i++)
            add(&s, GH_STAR_GENES[g].snps[i]);

    for (size_t i = 0; i < GH_SIMPLE_TRAIT_COUNT; i++)
        add(&s, GH_SIMPLE_TRAITS[i].rsid);
    /* The four multi-SNP traits: eye colour, hair colour, freckling, bitter
     * taste. Their SNPs live in the C selection logic, not a table. */
    add(&s, "rs12913832"); add(&s, "rs1800407");
    add(&s, "rs1805007"); add(&s, "rs1805008");
    add(&s, "rs713598"); add(&s, "rs1726866"); add(&s, "rs10246939");

    add_markers(&s, &GH_HISTAMINE_MARKERS);
    add_markers(&s, &GH_ALCOHOL_MARKERS);
    add_markers(&s, &GH_PAIN_MARKERS);
    add_markers(&s, &GH_THYROID_MARKERS);
    add_markers(&s, &GH_HORMONE_MARKERS);
    add_markers(&s, &GH_EYE_MARKERS);

    for (size_t i = 0; i < GH_CHRONOTYPE_SNP_COUNT; i++)
        add(&s, GH_CHRONOTYPE_SNPS[i].rsid);
    for (size_t i = 0; i < GH_LONGEVITY_SNP_COUNT; i++)
        add(&s, GH_LONGEVITY_SNPS[i].rsid);
    for (size_t i = 0; i < GH_MENTAL_SNP_COUNT; i++)
        add(&s, GH_MENTAL_SNPS[i].rsid);

    qsort(s.items, s.n, sizeof *s.items, compare_strp);
    *out = s.items;
    return s.n;
}

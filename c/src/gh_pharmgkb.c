#include "gh_pharmgkb.h"

#include <string.h>

#include "gh_tsv.h"

const char *const GH_PGX_ID_COLUMNS[] = {
    "Clinical Annotation ID",
    "Summary Annotation ID",
    NULL,
};

/* Row of the annotation table, keyed by annotation ID for the join. */
typedef struct {
    const char *rsid;
    const char *gene;
    const char *drugs;
    const char *phenotype;
    const char *level;
    const char *category;
} pgx_annotation;

bool gh_pharmgkb_load(gh_pharmgkb *p, gh_arena *arena,
                      const char *annotations_path, const char *alleles_path)
{
    memset(p, 0, sizeof(*p));
    p->arena = arena;

    /* Pass 1: the annotation table, keyed by annotation ID. */
    gh_tsv ann;
    if (!gh_tsv_open(&ann, arena, annotations_path))
        return false;
    if (!gh_tsv_read_header(&ann))
        return false;

    int c_id = gh_tsv_col_any(&ann, GH_PGX_ID_COLUMNS);
    int c_variant = gh_tsv_col(&ann, "Variant/Haplotypes");
    if (c_id < 0 || c_variant < 0)
        return false;

    int c_gene = gh_tsv_col(&ann, "Gene");
    int c_drugs = gh_tsv_col(&ann, "Drug(s)");
    int c_phen = gh_tsv_col(&ann, "Phenotype(s)");
    int c_level = gh_tsv_col(&ann, "Level of Evidence");
    int c_cat = gh_tsv_col(&ann, "Phenotype Category");

    gh_map by_ann_id;
    gh_map_init(&by_ann_id, arena, 8192);

    while (gh_tsv_next(&ann)) {
        /* Only annotations on a plain rsID are usable here; haplotype rows
         * (CYP2D6*4 and the like) are handled by the star-allele caller. */
        const char *variant = gh_tsv_at(&ann, c_variant);
        if (strncmp(variant, "rs", 2) != 0)
            continue;

        pgx_annotation *a = gh_alloc(arena, sizeof(*a));
        a->rsid = gh_strdup(arena, variant);
        a->gene = gh_strdup(arena, gh_tsv_at(&ann, c_gene));
        a->drugs = gh_strdup(arena, gh_tsv_at(&ann, c_drugs));
        a->phenotype = gh_strdup(arena, gh_tsv_at(&ann, c_phen));
        a->level = gh_strdup(arena, gh_tsv_at(&ann, c_level));
        a->category = gh_strdup(arena, gh_tsv_at(&ann, c_cat));

        gh_map_put(&by_ann_id, gh_tsv_at(&ann, c_id), a);
        p->annotations++;
    }

    /* Pass 2: the allele table drives which rsIDs exist. An rsID with no
     * allele rows never becomes an entry, and the metadata for an rsID comes
     * from the first annotation that reaches it here -- both behaviours match
     * loading.py, which builds its dict inside this second loop. */
    gh_tsv all;
    if (!gh_tsv_open(&all, arena, alleles_path))
        return false;
    if (!gh_tsv_read_header(&all))
        return false;

    int a_id = gh_tsv_col_any(&all, GH_PGX_ID_COLUMNS);
    int a_geno = gh_tsv_col(&all, "Genotype/Allele");
    int a_text = gh_tsv_col(&all, "Annotation Text");
    if (a_id < 0 || a_geno < 0)
        return false;

    gh_map_init(&p->by_rsid, arena, 4096);
    /* Upper bound: one entry per allele row. */
    size_t ordered_cap = 8192;
    p->ordered = gh_calloc(arena, ordered_cap, sizeof(*p->ordered));

    while (gh_tsv_next(&all)) {
        const pgx_annotation *a = gh_map_get(&by_ann_id, gh_tsv_at(&all, a_id));
        if (!a)
            continue;

        gh_pgx_entry *e = gh_map_get(&p->by_rsid, a->rsid);
        if (!e) {
            e = gh_alloc(arena, sizeof(*e));
            e->rsid = a->rsid;
            e->gene = a->gene;
            e->drugs = a->drugs;
            e->phenotype = a->phenotype;
            e->level = a->level;
            e->category = a->category;
            gh_map_init(&e->genotypes, arena, 8);
            gh_map_put(&p->by_rsid, a->rsid, e);

            if (p->nordered == ordered_cap) {
                size_t grown_cap = ordered_cap * 2;
                gh_pgx_entry **grown =
                    gh_calloc(arena, grown_cap, sizeof(*grown));
                memcpy(grown, p->ordered, p->nordered * sizeof(*grown));
                p->ordered = grown;
                ordered_cap = grown_cap;
            }
            p->ordered[p->nordered++] = e;
        }

        gh_map_put(&e->genotypes, gh_tsv_at(&all, a_geno),
                   gh_strdup(arena, gh_tsv_at(&all, a_text)));
        p->allele_rows++;
    }

    return true;
}

const gh_pgx_entry *gh_pharmgkb_find(const gh_pharmgkb *p, const char *rsid)
{
    return gh_map_get(&p->by_rsid, rsid);
}

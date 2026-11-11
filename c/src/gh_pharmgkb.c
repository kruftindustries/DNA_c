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

/* ------------------------------------------------------------------ */
/* Indel genotypes                                                     */
/* ------------------------------------------------------------------ */

/* ClinPGx spells an indel annotation's genotypes with allele sequences
 * and "del" -- CFTR F508del is CTT/CTT, CTT/del, del/del; a single-base
 * deletion is GG, G/del, del/del -- while an array export writes the same
 * three calls as II, DI, DD. Once an entry's rows are in, the array
 * spelling is derived for each row whose alleles are all either "del" or a
 * run of bases, and the entry answers to both. Entries with no "del" row
 * are left alone: D/I cannot tell (TA)6 from (TA)7, or G from TT. */

static bool allele_letter(const char *s, size_t n, char *out)
{
    if (n == 3 && (s[0] | 0x20) == 'd' && (s[1] | 0x20) == 'e' && (s[2] | 0x20) == 'l') {
        *out = 'D';
        return true;
    }
    if (n == 0)
        return false;
    for (size_t i = 0; i < n; i++)
        if (s[i] != 'A' && s[i] != 'C' && s[i] != 'G' && s[i] != 'T')
            return false;
    *out = 'I';
    return true;
}

/* "X/Y", or a slash-less row that is one allele written twice ("GG"). */
static bool indel_key_for(const char *row, char key[3])
{
    const char *a1 = row, *a2;
    size_t n1, n2;
    const char *slash = strchr(row, '/');
    if (slash) {
        n1 = (size_t)(slash - row);
        a2 = slash + 1;
        n2 = strlen(a2);
    } else {
        size_t n = strlen(row);
        if (n == 0 || n % 2)
            return false;
        n1 = n2 = n / 2;
        a2 = row + n1;
        if (strncmp(a1, a2, n1) != 0)
            return false;
    }
    char c1, c2;
    if (!allele_letter(a1, n1, &c1) || !allele_letter(a2, n2, &c2))
        return false;
    key[0] = c1 < c2 ? c1 : c2;   /* D before I */
    key[1] = c1 < c2 ? c2 : c1;
    key[2] = '\0';
    return true;
}

static bool mentions_del(const char *row)
{
    for (; *row; row++)
        if ((row[0] | 0x20) == 'd' && (row[1] | 0x20) == 'e' && (row[2] | 0x20) == 'l')
            return true;
    return false;
}

static void derive_indel_keys(gh_pgx_entry *e, gh_arena *arena)
{
    size_t it = 0;
    const char *row;
    void *text;
    bool has_del = false;
    while (gh_map_next(&e->genotypes, &it, &row, &text))
        if (mentions_del(row))
            has_del = true;
    if (!has_del)
        return;

    /* Collected first: the map is not grown while it is being walked. */
    struct { char key[3]; const char *text; } add[32];
    size_t n = 0;
    it = 0;
    while (gh_map_next(&e->genotypes, &it, &row, &text) && n < 32) {
        if (indel_key_for(row, add[n].key)) {
            add[n].text = text;
            n++;
        }
    }
    for (size_t i = 0; i < n; i++)
        if (!gh_map_has(&e->genotypes, add[i].key))
            gh_map_put(&e->genotypes, gh_strdup(arena, add[i].key), (void *)add[i].text);
}

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

    /* Pass 2: the allele table drives which annotations exist; one with no
     * allele rows never becomes an entry. One entry per annotation, not per
     * rsID: ClinPGx often carries several annotations on one variant at
     * different evidence levels and for different drugs (CFTR F508del has
     * seven, three of them 1A behind a level-3 and a level-4 one), and the
     * Python's one-dict-per-rsID kept whichever came first, so the report
     * could miss every reportable one. `by_rsid` keeps the first entry for
     * an rsID for lookups; `ordered` has them all. */
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
    gh_map by_ann;   /* annotation ID -> its entry */
    gh_map_init(&by_ann, arena, 8192);
    /* Upper bound: one entry per allele row. */
    size_t ordered_cap = 8192;
    p->ordered = gh_calloc(arena, ordered_cap, sizeof(*p->ordered));

    while (gh_tsv_next(&all)) {
        const char *ann_id = gh_tsv_at(&all, a_id);
        const pgx_annotation *a = gh_map_get(&by_ann_id, ann_id);
        if (!a)
            continue;

        gh_pgx_entry *e = gh_map_get(&by_ann, ann_id);
        if (!e) {
            e = gh_alloc(arena, sizeof(*e));
            e->rsid = a->rsid;
            e->gene = a->gene;
            e->drugs = a->drugs;
            e->phenotype = a->phenotype;
            e->level = a->level;
            e->category = a->category;
            gh_map_init(&e->genotypes, arena, 8);
            gh_map_put(&by_ann, gh_strdup(arena, ann_id), e);
            if (!gh_map_has(&p->by_rsid, a->rsid))
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

    for (size_t i = 0; i < p->nordered; i++)
        derive_indel_keys(p->ordered[i], arena);

    return true;
}

const gh_pgx_entry *gh_pharmgkb_find(const gh_pharmgkb *p, const char *rsid)
{
    return gh_map_get(&p->by_rsid, rsid);
}

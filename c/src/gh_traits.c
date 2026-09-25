#include "gh_traits.h"

#include <stdio.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Helpers                                                             */
/* ------------------------------------------------------------------ */

static const char *genotype_or_empty(const gh_genome *g, const char *rsid)
{
    const char *gt = gh_genome_genotype(g, rsid);
    return gt ? gt : "";
}

static size_t count_allele(const char *genotype, char allele)
{
    size_t n = 0;
    for (const char *p = genotype; *p; p++)
        if (*p == allele)
            n++;
    return n;
}

static const char *snp_label(const char *rsid)
{
    for (size_t i = 0; i < GH_TRAIT_SNP_LABEL_COUNT; i++)
        if (strcmp(GH_TRAIT_SNP_LABELS[i].key, rsid) == 0)
            return GH_TRAIT_SNP_LABELS[i].value;
    return rsid;
}

/* Collects the "rsID (LABEL): GENOTYPE" lines a prediction drew on. */
typedef struct {
    const char **items;
    size_t count;
    gh_arena *arena;
} snp_list;

static void snp_list_init(snp_list *l, gh_arena *arena, size_t capacity)
{
    l->items = gh_calloc(arena, capacity, sizeof(*l->items));
    l->count = 0;
    l->arena = arena;
}

/* Append an entry when the genotype is present. `label` may be NULL to look
 * it up from the shared table. */
static void snp_list_add(snp_list *l, const char *rsid, const char *label,
                         const char *genotype)
{
    if (!genotype || !*genotype)
        return;
    if (!label)
        label = snp_label(rsid);

    size_t n = strlen(rsid) + strlen(label) + strlen(genotype) + 8;
    char *entry = gh_alloc(l->arena, n);
    snprintf(entry, n, "%s (%s): %s", rsid, label, genotype);
    l->items[l->count++] = entry;
}

static void finish(gh_trait_result *out, const char *key,
                   const gh_trait_outcome *outcome, const snp_list *snps)
{
    out->key = key;
    out->prediction = outcome->prediction;
    out->confidence = outcome->confidence;
    out->description = outcome->description;
    out->snps_used = snps->items;
    out->nsnps = snps->count;
}

/* ------------------------------------------------------------------ */
/* Table-driven traits                                                 */
/* ------------------------------------------------------------------ */

/* Outcome indices shared by every simple trait. */
enum { OUT_UNKNOWN = 0, OUT_TWO = 1, OUT_ONE = 2, OUT_ZERO = 3 };

static void predict_simple(gh_trait_result *out, gh_arena *arena,
                           const gh_genome *g, const gh_simple_trait *trait)
{
    const char *genotype = genotype_or_empty(g, trait->rsid);

    snp_list snps;
    snp_list_init(&snps, arena, 1);
    snp_list_add(&snps, trait->rsid, trait->label, genotype);

    if (!*genotype) {
        finish(out, trait->key, &trait->outcomes[OUT_UNKNOWN], &snps);
        return;
    }

    size_t copies = count_allele(genotype, trait->allele);
    int index = copies == 2 ? OUT_TWO : copies == 1 ? OUT_ONE : OUT_ZERO;
    finish(out, trait->key, &trait->outcomes[index], &snps);
}

/* ------------------------------------------------------------------ */
/* Bespoke traits                                                      */
/* ------------------------------------------------------------------ */

/* Eye colour: HERC2 rs12913832 sets the base, OCA2 rs1800407 shifts it. */
enum {
    EYE_UNKNOWN = 0, EYE_BLUE_OCA2 = 1, EYE_BLUE = 2,
    EYE_HET_OCA2 = 3, EYE_HET = 4, EYE_BROWN = 5,
};

static void predict_eye_color(gh_trait_result *out, gh_arena *arena,
                              const gh_genome *g)
{
    const char *herc2 = genotype_or_empty(g, "rs12913832");
    const char *oca2 = genotype_or_empty(g, "rs1800407");

    snp_list snps;
    snp_list_init(&snps, arena, 2);
    snp_list_add(&snps, "rs12913832", NULL, herc2);
    snp_list_add(&snps, "rs1800407", NULL, oca2);

    if (!*herc2) {
        finish(out, "eye_color", &GH_EYE_COLOR_OUTCOMES[EYE_UNKNOWN], &snps);
        return;
    }

    size_t a_count = count_allele(herc2, 'A');
    bool oca2_variant = *oca2 && count_allele(oca2, 'T') >= 1;   /* plus-strand T */

    int index;
    if (a_count == 2)
        index = oca2_variant ? EYE_BLUE_OCA2 : EYE_BLUE;
    else if (a_count == 1)
        index = oca2_variant ? EYE_HET_OCA2 : EYE_HET;
    else
        index = EYE_BROWN;

    finish(out, "eye_color", &GH_EYE_COLOR_OUTCOMES[index], &snps);
}

/* Hair colour and freckling both count MC1R loss-of-function T alleles
 * across the same two SNPs, then branch on two-or-more / one / none. */
static int mc1r_branch(const gh_genome *g, snp_list *snps, bool *have_data)
{
    const char *r151c = genotype_or_empty(g, "rs1805007");
    const char *r160w = genotype_or_empty(g, "rs1805008");

    snp_list_add(snps, "rs1805007", NULL, r151c);
    snp_list_add(snps, "rs1805008", NULL, r160w);

    *have_data = *r151c || *r160w;
    if (!*have_data)
        return 0;

    size_t variants = count_allele(r151c, 'T') + count_allele(r160w, 'T');
    return variants >= 2 ? 1 : variants == 1 ? 2 : 3;
}

static void predict_hair_color(gh_trait_result *out, gh_arena *arena,
                               const gh_genome *g)
{
    snp_list snps;
    snp_list_init(&snps, arena, 2);
    bool have_data;
    int index = mc1r_branch(g, &snps, &have_data);
    finish(out, "hair_color", &GH_HAIR_COLOR_OUTCOMES[index], &snps);
}

static void predict_freckling(gh_trait_result *out, gh_arena *arena,
                              const gh_genome *g)
{
    snp_list snps;
    snp_list_init(&snps, arena, 2);
    bool have_data;
    int index = mc1r_branch(g, &snps, &have_data);
    finish(out, "freckling", &GH_FRECKLING_OUTCOMES[index], &snps);
}

/* Bitter taste: the proportion of taster alleles across three TAS2R38 SNPs.
 * Only SNPs that are present count toward the denominator. */
static void predict_bitter_taste(gh_trait_result *out, gh_arena *arena,
                                 const gh_genome *g)
{
    const char *a49p = genotype_or_empty(g, "rs713598");
    const char *v262a = genotype_or_empty(g, "rs1726866");
    const char *i296v = genotype_or_empty(g, "rs10246939");

    snp_list snps;
    snp_list_init(&snps, arena, 3);
    snp_list_add(&snps, "rs713598", NULL, a49p);
    snp_list_add(&snps, "rs1726866", NULL, v262a);
    snp_list_add(&snps, "rs10246939", NULL, i296v);

    if (!*a49p && !*v262a && !*i296v) {
        finish(out, "bitter_taste", &GH_BITTER_TASTE_OUTCOMES[0], &snps);
        return;
    }

    size_t taster = 0, checked = 0;
    if (*a49p) {
        taster += count_allele(a49p, 'G');
        checked += 2;
    }
    if (*v262a) {
        taster += count_allele(v262a, 'G');   /* plus-strand G */
        checked += 2;
    }
    if (*i296v) {
        taster += count_allele(i296v, 'C');
        checked += 2;
    }

    double ratio = checked ? (double)taster / (double)checked : 0.0;
    int index = ratio >= 0.8 ? 1 : ratio >= 0.4 ? 2 : 3;
    finish(out, "bitter_taste", &GH_BITTER_TASTE_OUTCOMES[index], &snps);
}

/* ------------------------------------------------------------------ */
/* Entry point                                                         */
/* ------------------------------------------------------------------ */

void gh_predict_traits(gh_trait_result *out, gh_arena *arena,
                       const gh_genome *g)
{
    memset(out, 0, GH_TRAIT_COUNT * sizeof(*out));

    /* Fill by key so the output order matches GH_TRAIT_ORDER regardless of
     * how the tables happen to be arranged. */
    gh_trait_result computed[GH_TRAIT_COUNT];
    size_t n = 0;

    predict_eye_color(&computed[n++], arena, g);
    predict_hair_color(&computed[n++], arena, g);
    predict_freckling(&computed[n++], arena, g);
    predict_bitter_taste(&computed[n++], arena, g);
    for (size_t i = 0; i < GH_SIMPLE_TRAIT_COUNT; i++)
        predict_simple(&computed[n++], arena, g, &GH_SIMPLE_TRAITS[i]);

    for (size_t i = 0; GH_TRAIT_ORDER[i]; i++) {
        for (size_t j = 0; j < n; j++) {
            if (strcmp(computed[j].key, GH_TRAIT_ORDER[i]) == 0) {
                out[i] = computed[j];
                break;
            }
        }
    }
}

const gh_trait_result *gh_trait_find(const gh_trait_result *results,
                                     const char *key)
{
    for (size_t i = 0; i < GH_TRAIT_COUNT; i++)
        if (results[i].key && strcmp(results[i].key, key) == 0)
            return &results[i];
    return NULL;
}

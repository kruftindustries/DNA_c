#include "gh_analysis.h"

#include <stdlib.h>
#include <string.h>

/* Both result lists are sorted before returning, matching the tail of
 * analyze_lifestyle_health. Python's sort is stable, so entries that tie
 * keep the order the SNP database gave them. */
typedef struct {
    gh_finding item;
    size_t order;
} ranked_finding;

static int compare_finding(const void *a, const void *b)
{
    const ranked_finding *x = a, *y = b;
    if (x->item.magnitude != y->item.magnitude)
        return x->item.magnitude > y->item.magnitude ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

typedef struct {
    gh_drug_finding item;
    size_t order;
} ranked_drug;

static int compare_drug(const void *a, const void *b)
{
    const ranked_drug *x = a, *y = b;
    int cmp = strcmp(x->item.level, y->item.level);
    if (cmp != 0)
        return cmp;
    return (x->order > y->order) - (x->order < y->order);
}

bool gh_pgx_level_reportable(const char *level)
{
    if (!level)
        return false;
    return strcmp(level, "1A") == 0 || strcmp(level, "1B") == 0 ||
           strcmp(level, "2A") == 0 || strcmp(level, "2B") == 0;
}

/* Genotype lookup in a ClinPGx genotype map, trying the reversed spelling as
 * loading.py does. */
static const char *pgx_lookup(const gh_map *genotypes, const char *genotype)
{
    const char *hit = gh_map_get(genotypes, genotype);
    if (hit)
        return hit;
    if (strlen(genotype) == 2) {
        char swapped[3] = {genotype[1], genotype[0], '\0'};
        return gh_map_get(genotypes, swapped);
    }
    return NULL;
}

void gh_analyze(gh_analysis *out, gh_arena *arena, const gh_genome *genome,
                const gh_pharmgkb *pgx)
{
    memset(out, 0, sizeof(*out));
    out->total_snps = gh_map_len(&genome->by_rsid);

    /* At most one finding per database SNP. */
    out->findings = gh_calloc(arena, GH_SNP_COUNT, sizeof(*out->findings));

    for (size_t i = 0; i < GH_SNP_COUNT; i++) {
        const gh_snp *snp = &GH_SNPS[i];
        const char *genotype = gh_genome_genotype(genome, snp->rsid);
        if (!genotype)
            continue;

        const gh_snp_variant *v = gh_snp_match(snp, genotype);
        if (!v)
            continue;

        gh_finding *f = &out->findings[out->nfindings++];
        f->rsid = snp->rsid;
        f->gene = snp->gene;
        f->category = snp->category;
        f->genotype = genotype;
        f->status = v->status;
        f->description = v->desc;
        f->magnitude = v->magnitude;
        f->note = snp->note ? snp->note : "";
        f->freq = snp->freq;

        out->analyzed_snps++;
        if (v->magnitude >= 3)
            out->high_impact++;
        else if (v->magnitude >= 2)
            out->moderate_impact++;
        else if (v->magnitude >= 1)
            out->low_impact++;
    }

    /* Descending magnitude, ties keeping database order. */
    {
        ranked_finding *ranked =
            gh_calloc(arena, out->nfindings ? out->nfindings : 1,
                      sizeof(*ranked));
        for (size_t i = 0; i < out->nfindings; i++) {
            ranked[i].item = out->findings[i];
            ranked[i].order = i;
        }
        qsort(ranked, out->nfindings, sizeof(*ranked), compare_finding);
        for (size_t i = 0; i < out->nfindings; i++)
            out->findings[i] = ranked[i].item;
    }

    if (!pgx)
        return;

    /* Walk the ClinPGx entries and keep those the genome genotypes and whose
     * evidence level is strong enough to report. */
    size_t cap = gh_map_len(&pgx->by_rsid);
    if (cap == 0)
        return;
    out->drug_findings = gh_calloc(arena, cap, sizeof(*out->drug_findings));

    /* Insertion order, not hash order: the sort below is stable, so this
     * is what decides how equal evidence levels are arranged. */
    for (size_t i = 0; i < pgx->nordered; i++) {
        const gh_pgx_entry *e = pgx->ordered[i];

        const char *genotype = gh_genome_genotype(genome, e->rsid);
        if (!genotype)
            continue;

        const char *annotation = pgx_lookup(&e->genotypes, genotype);
        if (!annotation || !gh_pgx_level_reportable(e->level))
            continue;

        gh_drug_finding *d = &out->drug_findings[out->ndrug_findings++];
        d->rsid = e->rsid;
        d->gene = e->gene;
        d->drugs = e->drugs;
        d->genotype = genotype;
        d->annotation = annotation;
        d->level = e->level;
        d->category = e->category;
    }

    /* Ascending evidence level as a string, so 1A precedes 1B precedes 2A. */
    {
        ranked_drug *ranked =
            gh_calloc(arena, out->ndrug_findings ? out->ndrug_findings : 1,
                      sizeof(*ranked));
        for (size_t i = 0; i < out->ndrug_findings; i++) {
            ranked[i].item = out->drug_findings[i];
            ranked[i].order = i;
        }
        qsort(ranked, out->ndrug_findings, sizeof(*ranked), compare_drug);
        for (size_t i = 0; i < out->ndrug_findings; i++)
            out->drug_findings[i] = ranked[i].item;
    }
}

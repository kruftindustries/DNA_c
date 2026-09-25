#include "gh_scorers.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Shared helpers                                                      */
/* ------------------------------------------------------------------ */

/* Genotype for an rsID, or "" when absent -- the shape the Python gets from
 * `genome_by_rsid.get(rsid, {}).get("genotype", "")`. */
static const char *genotype_or_empty(const gh_genome *g, const char *rsid)
{
    const char *gt = gh_genome_genotype(g, rsid);
    return gt ? gt : "";
}

static const char *kv_lookup(const gh_kv *table, size_t n, const char *key,
                             const char *fallback)
{
    for (size_t i = 0; i < n; i++)
        if (strcmp(table[i].key, key) == 0)
            return table[i].value;
    return fallback;
}

/* ------------------------------------------------------------------ */
/* APOE                                                                */
/* ------------------------------------------------------------------ */

/* Decode one allele pair into an epsilon allele, or NULL if the pair is not
 * a natural APOE haplotype (C at rs429358 with T at rs7412). */
static const char *epsilon_for(char a358, char a7412)
{
    if (a358 == 'T' && a7412 == 'T')
        return "e2";
    if (a358 == 'T' && a7412 == 'C')
        return "e3";
    if (a358 == 'C' && a7412 == 'C')
        return "e4";
    return NULL;
}

/* Decode a specific phasing into a canonical "eX/eY" haplotype. */
static const char *decode_haplotype(gh_arena *arena, const char *g358,
                                    const char *g7412)
{
    if (strlen(g358) != 2 || strlen(g7412) != 2)
        return NULL;

    const char *first = epsilon_for(g358[0], g7412[0]);
    if (!first)
        return NULL;
    const char *second = epsilon_for(g358[1], g7412[1]);
    if (!second)
        return NULL;

    /* Canonical order, matching Python's sort of the two allele names. */
    if (strcmp(first, second) > 0) {
        const char *tmp = first;
        first = second;
        second = tmp;
    }

    char *out = gh_alloc(arena, 8);
    snprintf(out, 8, "%s/%s", first, second);
    return out;
}

void gh_call_apoe(gh_apoe_result *out, gh_arena *arena, const gh_genome *g)
{
    memset(out, 0, sizeof(*out));
    out->rs429358 = genotype_or_empty(g, "rs429358");
    out->rs7412 = genotype_or_empty(g, "rs7412");

    if (!*out->rs429358 || !*out->rs7412) {
        out->apoe_type = "Unknown";
        out->risk_level = "unknown";
        out->description =
            "Insufficient data — both rs429358 and rs7412 are needed.";
        out->confidence = "low";
        return;
    }

    /* 23andMe calls are unphased, so try both orderings of each genotype and
     * take the first combination that decodes. */
    const char *haplotype = NULL;
    if (strlen(out->rs429358) == 2 && strlen(out->rs7412) == 2) {
        char p358[2][3] = {
            {out->rs429358[0], out->rs429358[1], '\0'},
            {out->rs429358[1], out->rs429358[0], '\0'},
        };
        char p7412[2][3] = {
            {out->rs7412[0], out->rs7412[1], '\0'},
            {out->rs7412[1], out->rs7412[0], '\0'},
        };
        for (int i = 0; i < 2 && !haplotype; i++)
            for (int j = 0; j < 2 && !haplotype; j++)
                haplotype = decode_haplotype(arena, p358[i], p7412[j]);
    }

    if (!haplotype) {
        char *desc = gh_alloc(arena, 160);
        snprintf(desc, 160,
                 "Unexpected genotype combination: rs429358=%s, rs7412=%s.",
                 out->rs429358, out->rs7412);
        out->apoe_type = "Unknown";
        out->risk_level = "unknown";
        out->description = desc;
        out->confidence = "low";
        return;
    }

    for (size_t i = 0; i < GH_APOE_RISK_COUNT; i++) {
        if (strcmp(GH_APOE_RISKS[i].haplotype, haplotype) == 0) {
            out->apoe_type = haplotype;
            out->risk_level = GH_APOE_RISKS[i].risk_level;
            out->alzheimer_or = GH_APOE_RISKS[i].alzheimer_or;
            out->has_or = true;
            out->description = GH_APOE_RISKS[i].description;
            out->confidence = "high";
            return;
        }
    }

    /* Every decodable haplotype is in the table; this is unreachable unless
     * the Python table loses an entry. */
    out->apoe_type = haplotype;
    out->risk_level = "unknown";
    out->description = "No risk information for this haplotype.";
    out->confidence = "low";
}

/* ------------------------------------------------------------------ */
/* Blood type                                                          */
/* ------------------------------------------------------------------ */

void gh_predict_blood_type(gh_blood_result *out, gh_arena *arena,
                           const gh_genome *g)
{
    memset(out, 0, sizeof(*out));
    out->rs505922 = genotype_or_empty(g, "rs505922");    /* ABO O proxy */
    out->rs8176746 = genotype_or_empty(g, "rs8176746");  /* B antigen */
    out->rs590787 = genotype_or_empty(g, "rs590787");    /* RhD proxy */

    const char *proxy = out->rs505922;
    const char *b = out->rs8176746;

    int abo_confidence = 0;
    if (*proxy)
        abo_confidence++;
    if (*b)
        abo_confidence++;

    bool has_b = *b && strchr(b, 'T') != NULL;

    const char *abo = NULL;
    if (strcmp(proxy, "TT") == 0) {
        /* Two O alleles; a B signal still wins, as in the Python. */
        abo = has_b ? "B" : "O";
    } else if (strcmp(proxy, "CT") == 0 || strcmp(proxy, "TC") == 0) {
        abo = has_b ? "B" : "A";
    } else if (strcmp(proxy, "CC") == 0) {
        if (strcmp(b, "TT") == 0)
            abo = "B";
        else if (has_b)
            abo = "AB";
        else
            abo = "A";
    } else {
        abo = has_b ? "B" : NULL;
    }

    const char *rh = NULL;
    int rh_confidence = 0;
    if (*out->rs590787) {
        rh_confidence = 1;
        rh = strcmp(out->rs590787, "CC") == 0 ? "-" : "+";   /* plus-strand A/A */
    }

    if (abo && rh) {
        char *bt = gh_alloc(arena, 8);
        snprintf(bt, 8, "%s%s", abo, rh);
        out->blood_type = bt;
    } else if (abo) {
        char *bt = gh_alloc(arena, 8);
        snprintf(bt, 8, "%s?", abo);
        out->blood_type = bt;
    } else if (rh) {
        char *bt = gh_alloc(arena, 8);
        snprintf(bt, 8, "?%s", rh);
        out->blood_type = bt;
    } else {
        out->blood_type = "Unknown";
    }

    out->abo = abo ? abo : "Unknown";
    out->rh = rh ? rh : "Unknown";

    int total = abo_confidence + rh_confidence;
    out->confidence = total >= 3 ? "high" : total >= 1 ? "moderate" : "low";
}

/* ------------------------------------------------------------------ */
/* Mitochondrial haplogroup                                            */
/* ------------------------------------------------------------------ */

void gh_estimate_mt_haplogroup(gh_mt_result *out, gh_arena *arena,
                               const gh_genome *g)
{
    memset(out, 0, sizeof(*out));
    out->markers_tested = GH_MT_TREE_COUNT;
    out->matches = gh_calloc(arena, GH_MT_TREE_COUNT, sizeof(*out->matches));

    for (size_t i = 0; i < GH_MT_TREE_COUNT; i++) {
        const gh_mt_marker *m = &GH_MT_TREE[i];
        const char *genotype = gh_genome_genotype(g, m->rsid);
        if (!genotype || !*genotype)
            continue;

        out->markers_found++;

        /* MT DNA is haploid, so the call is usually one character, but a
         * doubled call means the same thing: look for the allele anywhere. */
        if (strstr(genotype, m->allele))
            out->matches[out->nmatches++] = m;
    }

    if (out->nmatches) {
        /* Later entries in the tree are more specific, so the last match wins. */
        const gh_mt_marker *best = out->matches[out->nmatches - 1];
        out->haplogroup = best->haplogroup;
        out->description = best->description;
    } else {
        out->haplogroup = "Unknown";
        out->description = "No defining mitochondrial SNPs matched";
    }

    if (out->markers_found >= 10)
        out->confidence = "high";
    else if (out->markers_found >= 5)
        out->confidence = "moderate";
    else if (out->markers_found >= 1)
        out->confidence = "low";
    else
        out->confidence = "none";

    const char *region = kv_lookup(GH_MT_LINEAGES, GH_MT_LINEAGE_COUNT,
                                   out->haplogroup, "Unknown");
    size_t n = strlen(region) + strlen(" maternal") + 1;
    char *lineage = gh_alloc(arena, n);
    snprintf(lineage, n, "%s maternal", region);
    out->lineage = lineage;
}

/* ------------------------------------------------------------------ */
/* Star alleles                                                        */
/* ------------------------------------------------------------------ */

const char *gh_star_phenotype(const char *function_a, const char *function_b)
{
    /* The table is keyed on the sorted pair. */
    const char *lo = function_a, *hi = function_b;
    if (strcmp(lo, hi) > 0) {
        const char *tmp = lo;
        lo = hi;
        hi = tmp;
    }
    for (size_t i = 0; i < GH_PHENOTYPE_RULE_COUNT; i++)
        if (strcmp(GH_PHENOTYPE_RULES[i].function_a, lo) == 0 &&
            strcmp(GH_PHENOTYPE_RULES[i].function_b, hi) == 0)
            return GH_PHENOTYPE_RULES[i].phenotype;
    return NULL;
}

/* A candidate star allele in assignment order. */
typedef struct {
    const gh_star_allele *allele;
    size_t initial;   /* copies available before any assignment */
    size_t order;     /* definition order, to keep the sort stable */
} candidate;

static int compare_candidates(const void *a, const void *b)
{
    const candidate *x = a, *y = b;
    /* Python sorts on (-nsnps, -copies) and its sort is stable. */
    if (x->allele->nsnps != y->allele->nsnps)
        return x->allele->nsnps > y->allele->nsnps ? -1 : 1;
    if (x->initial != y->initial)
        return x->initial > y->initial ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

/* Index of a base in the per-SNP copy pool, or -1. */
static int base_index(char c)
{
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default:  return -1;
    }
}

/* Position of `rsid` in the gene's SNP list. */
static size_t snp_index(const gh_star_gene *gene, const char *rsid)
{
    for (size_t i = 0; i < gene->nsnps; i++)
        if (strcmp(gene->snps[i], rsid) == 0)
            return i;
    return gene->nsnps;   /* generated tables never let this happen */
}

/* Copies of a haplotype the pool can still supply: the minimum over its
 * defining SNPs. */
static size_t pool_copies(const gh_star_allele *allele, const gh_star_gene *gene,
                          size_t (*pool)[4])
{
    size_t copies = SIZE_MAX;
    for (size_t j = 0; j < allele->nsnps; j++) {
        int b = base_index(allele->snps[j].allele[0]);
        size_t si = snp_index(gene, allele->snps[j].rsid);
        size_t n = (b < 0 || si >= gene->nsnps) ? 0 : pool[si][b];
        if (n < copies)
            copies = n;
    }
    return copies == SIZE_MAX ? 0 : copies;
}

static void call_gene(gh_star_result *out, gh_arena *arena,
                      const gh_star_gene *gene, const gh_genome *g)
{
    memset(out, 0, sizeof(*out));
    out->gene = gene->gene;
    out->snps_total = gene->nsnps;

    /* Which defining SNPs the genome carries. */
    const char **missing = gh_calloc(arena, gene->nsnps + 1, sizeof(*missing));
    size_t nmissing = 0;
    for (size_t i = 0; i < gene->nsnps; i++) {
        const char *gt = gh_genome_genotype(g, gene->snps[i]);
        if (gt && *gt)
            out->snps_found++;
        else
            missing[nmissing++] = gene->snps[i];
    }

    if (out->snps_found == 0) {
        char *note = gh_alloc(arena, 96);
        snprintf(note, 96, "No defining SNPs found for %s in genome data.",
                 gene->gene);
        out->diplotype = "Unknown";
        out->phenotype = "Unknown";
        out->coverage = 0.0;
        out->confidence = "low";
        out->clinical_note = note;
        return;
    }

    /* Haplotype assignment. The variant copies at each defining SNP form a
     * pool that the star alleles consume, most specific (most SNPs) first:
     * a haplotype defined by several SNPs is present only as often as its
     * least-present SNP, and the copies it uses are no longer available to
     * the alleles it contains. SLCO1B1 rs4149056 T/C + rs2306283 G/G is
     * therefore *15 on one copy and *17 on the other -- not *15 homozygous, which summing the copies produced
     * for anyone homozygous for the common 388G. Two haplotype slots exist;
     * among equally specific alleles the one with more copies goes first,
     * then definition order. */
    size_t (*pool)[4] = gh_calloc(arena, gene->nsnps, sizeof(*pool));
    bool *present = gh_calloc(arena, gene->nsnps, sizeof(*present));
    for (size_t i = 0; i < gene->nsnps; i++) {
        const char *gt = gh_genome_genotype(g, gene->snps[i]);
        if (!gt || !*gt)
            continue;
        present[i] = true;
        for (const char *c = gt; *c; c++) {
            int b = base_index(*c);
            if (b >= 0)
                pool[i][b]++;
        }
    }

    candidate *cands = gh_calloc(arena, gene->nalleles, sizeof(*cands));
    size_t ncands = 0;
    for (size_t i = 0; i < gene->nalleles; i++) {
        const gh_star_allele *allele = &gene->alleles[i];
        bool all_defined = true;
        for (size_t j = 0; j < allele->nsnps && all_defined; j++) {
            size_t si = snp_index(gene, allele->snps[j].rsid);
            all_defined = si < gene->nsnps && present[si];
        }
        if (!all_defined)
            continue;
        cands[ncands].allele = allele;
        cands[ncands].initial = pool_copies(allele, gene, pool);
        cands[ncands].order = ncands;
        ncands++;
    }
    qsort(cands, ncands, sizeof(*cands), compare_candidates);

    const char *haplotypes[2] = {"*1", "*1"};
    size_t nhap = 0;
    for (size_t i = 0; i < ncands && nhap < 2; i++) {
        const gh_star_allele *allele = cands[i].allele;
        size_t copies = pool_copies(allele, gene, pool);
        if (copies > 2 - nhap)
            copies = 2 - nhap;
        if (copies == 0)
            continue;
        for (size_t j = 0; j < allele->nsnps; j++) {
            int b = base_index(allele->snps[j].allele[0]);
            size_t si = snp_index(gene, allele->snps[j].rsid);
            if (b >= 0 && si < gene->nsnps)
                pool[si][b] -= copies;
        }
        for (size_t k = 0; k < copies; k++)
            haplotypes[nhap++] = allele->name;
    }
    const char *allele1 = haplotypes[0];
    const char *allele2 = haplotypes[1];

    /* Python sorts the two names as strings, so "*10" precedes "*2". */
    const char *lo = allele1, *hi = allele2;
    if (strcmp(lo, hi) > 0) {
        const char *tmp = lo;
        lo = hi;
        hi = tmp;
    }
    size_t dn = strlen(lo) + 1 + strlen(hi) + 1;
    char *diplotype = gh_alloc(arena, dn);
    snprintf(diplotype, dn, "%s/%s", lo, hi);
    out->diplotype = diplotype;

    const char *func1 = kv_lookup(gene->functions, gene->nfunctions, lo, "normal");
    const char *func2 = kv_lookup(gene->functions, gene->nfunctions, hi, "normal");

    gh_strbuf notes;
    gh_sb_init(&notes);

    const char *phenotype = gh_star_phenotype(func1, func2);
    if (!phenotype) {
        phenotype = "Indeterminate";
        gh_sb_printf(&notes,
            "Unexpected function combination (%s + %s). "
            "Clinical interpretation needed.", func1, func2);
    }
    out->phenotype = phenotype;

    if (gene->clinical_note) {
        if (notes.len)
            gh_sb_putc(&notes, ' ');
        gh_sb_puts(&notes, gene->clinical_note);
    }
    if (nmissing) {
        if (notes.len)
            gh_sb_putc(&notes, ' ');
        gh_sb_puts(&notes, "Missing SNPs: ");
        for (size_t i = 0; i < nmissing; i++) {
            if (i)
                gh_sb_puts(&notes, ", ");
            gh_sb_puts(&notes, missing[i]);
        }
        gh_sb_puts(&notes, ". Result based on available data only.");
    }

    out->clinical_note = notes.len
        ? gh_strndup(arena, notes.data, notes.len)
        : "All defining SNPs found.";
    gh_sb_free(&notes);

    double coverage = gene->nsnps
        ? (double)out->snps_found / (double)gene->nsnps : 0.0;
    /* Python rounds to 2 decimals before reporting. */
    out->coverage = nearbyint(coverage * 100.0) / 100.0;
    out->confidence = coverage >= 0.8 ? "high"
                    : coverage >= 0.5 ? "moderate" : "low";
}

void gh_call_star_alleles(gh_star_result *out, gh_arena *arena,
                          const gh_genome *g)
{
    for (size_t i = 0; i < GH_STAR_GENE_COUNT; i++)
        call_gene(&out[i], arena, &GH_STAR_GENES[i], g);
}

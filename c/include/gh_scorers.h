/* Self-contained scorers: APOE, blood type, mitochondrial haplogroup and
 * star alleles.
 *
 * Mirrors genetic_health/apoe.py, blood_type.py, mt_haplogroup.py and
 * star_alleles.py. The lookup tables are generated from those modules by
 * tools/gen_scorer_data.py.
 */
#ifndef GH_SCORERS_H
#define GH_SCORERS_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_genome.h"
#include "gh_mem.h"

/* Reused for simple key/value tables. Defined in gh_report.h too, so guard
 * against a duplicate typedef when both headers are included. */
#ifndef GH_KV_DEFINED
#define GH_KV_DEFINED
typedef struct {
    const char *key;
    const char *value;
} gh_kv;
#endif

/* ------------------------------------------------------------------ */
/* APOE                                                                */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *haplotype;    /* "e3/e4" */
    const char *risk_level;
    double alzheimer_or;
    const char *description;
} gh_apoe_risk;

extern const gh_apoe_risk GH_APOE_RISKS[];
extern const size_t GH_APOE_RISK_COUNT;

typedef struct {
    const char *apoe_type;     /* "e3/e4", or "Unknown" */
    const char *risk_level;
    double alzheimer_or;       /* 0 when unknown */
    bool has_or;
    const char *description;
    const char *confidence;
    const char *rs429358;      /* "" when absent */
    const char *rs7412;
} gh_apoe_result;

void gh_call_apoe(gh_apoe_result *out, gh_arena *arena, const gh_genome *g);

/* ------------------------------------------------------------------ */
/* Blood type                                                          */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *blood_type;   /* "A+", "O-", "A?", "?+", "Unknown" */
    const char *abo;          /* "A"/"B"/"AB"/"O"/"Unknown" */
    const char *rh;           /* "+"/"-"/"Unknown" */
    const char *confidence;
    const char *rs505922;     /* "" when absent */
    const char *rs8176746;
    const char *rs590787;
} gh_blood_result;

void gh_predict_blood_type(gh_blood_result *out, gh_arena *arena,
                           const gh_genome *g);

/* ------------------------------------------------------------------ */
/* Mitochondrial haplogroup                                            */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *rsid;
    const char *allele;
    const char *haplogroup;
    const char *description;
} gh_mt_marker;

extern const gh_mt_marker GH_MT_TREE[];
extern const size_t GH_MT_TREE_COUNT;
extern const gh_kv GH_MT_LINEAGES[];
extern const size_t GH_MT_LINEAGE_COUNT;

typedef struct {
    const char *haplogroup;
    const char *description;
    const char *confidence;
    const char *lineage;       /* "<region> maternal" */
    size_t markers_found;
    size_t markers_tested;
    const gh_mt_marker **matches;  /* markers whose allele was present */
    size_t nmatches;
} gh_mt_result;

void gh_estimate_mt_haplogroup(gh_mt_result *out, gh_arena *arena,
                               const gh_genome *g);

/* ------------------------------------------------------------------ */
/* Star alleles                                                        */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *rsid;
    const char *allele;     /* the variant base defining this star allele */
} gh_star_snp;

typedef struct {
    const char *name;       /* "*2" */
    const char *function;   /* normal / increased / decreased / no_function */
    const gh_star_snp *snps;
    size_t nsnps;
} gh_star_allele;

typedef struct {
    const char *gene;
    const gh_star_allele *alleles;
    size_t nalleles;
    const char *const *snps;   /* every SNP the gene is called from */
    size_t nsnps;
    /* Full star-allele -> function map, including *1 and any allele the
     * calling rules cannot currently produce. Looked up by name so the C
     * matches Python's function_map.get(name, "normal") exactly. */
    const gh_kv *functions;
    size_t nfunctions;
    const char *clinical_note;       /* NULL when the gene carries none */
} gh_star_gene;

extern const gh_star_gene GH_STAR_GENES[];
extern const size_t GH_STAR_GENE_COUNT;

typedef struct {
    const char *function_a;   /* sorted pair */
    const char *function_b;
    const char *phenotype;
} gh_phenotype_rule;

extern const gh_phenotype_rule GH_PHENOTYPE_RULES[];
extern const size_t GH_PHENOTYPE_RULE_COUNT;

typedef struct {
    const char *gene;
    const char *diplotype;    /* e.g. star-1 over star-2, written "*1" "/" "*2" */
    const char *phenotype;
    const char *confidence;
    const char *clinical_note;
    size_t snps_found;
    size_t snps_total;
    double coverage;          /* rounded to 2 decimals, as in Python */
} gh_star_result;

/* Call every pharmacogene. `out` must have room for GH_STAR_GENE_COUNT
 * entries; the results are in table order. */
void gh_call_star_alleles(gh_star_result *out, gh_arena *arena,
                          const gh_genome *g);

/* Phenotype for a pair of function names, or NULL when the pair is not in
 * the table. The arguments may be given in either order. */
const char *gh_star_phenotype(const char *function_a, const char *function_b);

#endif /* GH_SCORERS_H */

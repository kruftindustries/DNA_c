/* ClinVar disease-risk analysis, ACMG secondary findings and carrier
 * screening.
 *
 * Mirrors genetic_health/analysis.py:load_clinvar_and_analyze,
 * genetic_health/acmg.py and genetic_health/carrier_screen.py.
 */
#ifndef GH_CLINVAR_H
#define GH_CLINVAR_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_genome.h"
#include "gh_mem.h"

#ifndef GH_KV_DEFINED
#define GH_KV_DEFINED
typedef struct {
    const char *key;
    const char *value;
} gh_kv;
#endif

/* ------------------------------------------------------------------ */
/* Generated tables                                                    */
/* ------------------------------------------------------------------ */

extern const char *const GH_ACMG_GENES[];
extern const size_t GH_ACMG_GENE_COUNT;
extern const gh_kv GH_ACMG_ACTIONABILITY[];
extern const size_t GH_ACMG_ACTIONABILITY_COUNT;
extern const char *const GH_ACMG_DEFAULT_ACTIONABILITY;

extern const gh_kv GH_GENE_SYSTEMS[];
extern const size_t GH_GENE_SYSTEM_COUNT;
extern const gh_kv GH_GENE_INHERITANCE[];
extern const size_t GH_GENE_INHERITANCE_COUNT;
extern const char *const GH_COUPLES_GENES[];
extern const size_t GH_COUPLES_GENE_COUNT;

/* ------------------------------------------------------------------ */
/* ClinVar findings                                                    */
/* ------------------------------------------------------------------ */

/* Categories a matched variant can fall into, tested in this order. */
typedef enum {
    GH_CV_PATHOGENIC = 0,
    GH_CV_LIKELY_PATHOGENIC,
    GH_CV_RISK_FACTOR,
    GH_CV_DRUG_RESPONSE,
    GH_CV_PROTECTIVE,
    GH_CV_CATEGORY_COUNT,
} gh_cv_category;

extern const char *const GH_CV_CATEGORY_NAMES[GH_CV_CATEGORY_COUNT];

typedef struct {
    const char *chromosome;
    const char *position;
    const char *rsid;
    const char *gene;
    const char *ref;
    const char *alt;
    const char *user_genotype;
    bool is_homozygous;
    bool is_heterozygous;
    const char *clinical_significance;
    const char *review_status;
    int gold_stars;
    const char *traits;
    const char *inheritance;
    const char *hgvs_p;
    const char *hgvs_c;
    const char *molecular_consequence;
    const char *xrefs;
    const char *zygosity;         /* human-readable description */
    const char *zygosity_status;  /* AFFECTED / CARRIER / HETEROZYGOUS / UNKNOWN */
    gh_cv_category category;
} gh_cv_finding;

typedef struct {
    bool loaded;                  /* false when the ClinVar file is absent */

    gh_cv_finding **by_category;  /* GH_CV_CATEGORY_COUNT arrays */
    size_t counts[GH_CV_CATEGORY_COUNT];

    size_t total_clinvar;         /* rows scanned */
    size_t matched;               /* rows at a genotyped position */
} gh_clinvar_result;

/* Scan ClinVar for variants the genome carries. Returns false, with
 * `out->loaded` false, when the file cannot be read -- the caller then skips
 * disease-risk analysis as the Python does. */
bool gh_clinvar_analyze(gh_clinvar_result *out, gh_arena *arena,
                        const gh_genome *g, const char *clinvar_path);

/* Zygosity classification. `status` and `description` are set from the
 * variant's zygosity and the gene's inheritance mode. */
void gh_classify_zygosity(bool is_homozygous, bool is_heterozygous,
                          const char *inheritance,
                          const char **status, const char **description);

/* ------------------------------------------------------------------ */
/* ACMG secondary findings                                             */
/* ------------------------------------------------------------------ */

typedef struct {
    const gh_cv_finding *finding;
    const char *acmg_category;     /* "pathogenic" / "likely_pathogenic" */
    const char *actionability;
} gh_acmg_finding;

typedef struct {
    gh_acmg_finding *findings;
    size_t nfindings;
    size_t genes_screened;
    size_t genes_with_variants;
    const char *summary;
} gh_acmg_result;

void gh_flag_acmg(gh_acmg_result *out, gh_arena *arena,
                  const gh_clinvar_result *clinvar);

bool gh_acmg_is_actionable_gene(const char *gene_upper);

/* ------------------------------------------------------------------ */
/* Carrier screening                                                   */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *gene;
    const char *condition;
    const char *inheritance;
    const char *system;
    const char *reproductive_note;   /* "" when none applies */
    bool couples_relevant;
    const char *rsid;
    const char *genotype;
    int gold_stars;
} gh_carrier;

typedef struct {
    gh_carrier *carriers;
    size_t ncarriers;
    size_t ncouples_relevant;

    /* Distinct systems, in first-seen order, with a count each. */
    const char **systems;
    size_t *system_counts;
    size_t nsystems;
} gh_carrier_result;

void gh_organize_carriers(gh_carrier_result *out, gh_arena *arena,
                          const gh_clinvar_result *clinvar);

#endif /* GH_CLINVAR_H */

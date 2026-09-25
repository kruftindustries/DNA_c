/* Visible trait predictions.
 *
 * Mirrors genetic_health/traits.py. Ten of the fourteen traits count one
 * allele in one SNP and branch on two / one / zero copies, so they share a
 * table-driven path; the other four (eye colour, hair colour, freckling and
 * bitter taste) combine several SNPs and have their own logic.
 *
 * The prediction, confidence and description strings are generated from the
 * Python by tools/gen_traits_data.py.
 */
#ifndef GH_TRAITS_H
#define GH_TRAITS_H

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

typedef struct {
    const char *prediction;
    const char *confidence;
    const char *description;
} gh_trait_outcome;

/* A trait decided by counting one allele in one SNP.
 * `outcomes` holds four entries, in order: no data, two copies, one copy,
 * zero copies. */
typedef struct {
    const char *key;        /* "earwax_type" */
    const char *rsid;
    const char *label;      /* gene label used in the SNP list */
    char allele;            /* the counted base */
    const gh_trait_outcome *outcomes;
} gh_simple_trait;

extern const gh_simple_trait GH_SIMPLE_TRAITS[];
extern const size_t GH_SIMPLE_TRAIT_COUNT;

/* Bespoke traits. Branch order matches the C selection logic and is
 * documented in the generated file. */
extern const gh_trait_outcome GH_EYE_COLOR_OUTCOMES[];
extern const gh_trait_outcome GH_HAIR_COLOR_OUTCOMES[];
extern const gh_trait_outcome GH_FRECKLING_OUTCOMES[];
extern const gh_trait_outcome GH_BITTER_TASTE_OUTCOMES[];

/* Trait keys in the order predict_traits() reports them. NULL-terminated. */
extern const char *const GH_TRAIT_ORDER[];

extern const gh_kv GH_TRAIT_SNP_LABELS[];
extern const size_t GH_TRAIT_SNP_LABEL_COUNT;

/* One predicted trait. `snps_used` holds `nsnps` entries already formatted
 * as "rsID (LABEL): GENOTYPE", matching the Python. */
typedef struct {
    const char *key;
    const char *prediction;
    const char *confidence;
    const char *description;
    const char **snps_used;
    size_t nsnps;
} gh_trait_result;

/* Predict every trait. `out` must have room for GH_TRAIT_COUNT entries; the
 * results are in GH_TRAIT_ORDER order. */
#define GH_TRAIT_COUNT 14

void gh_predict_traits(gh_trait_result *out, gh_arena *arena,
                       const gh_genome *g);

/* Result for one trait key, or NULL when the key is unknown. */
const gh_trait_result *gh_trait_find(const gh_trait_result *results,
                                     const char *key);

#endif /* GH_TRAITS_H */

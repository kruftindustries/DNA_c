/* Research-backed genomic insights.
 *
 * Mirrors genetic_health/insights.py: single-gene entries keyed on a
 * gene/status pair, multi-gene narratives that fire when enough of a
 * pattern's genes match, a capped highlight reel, and the protective
 * findings worth calling out.
 *
 * The Python's `ancestry_results` and `epistasis_results` parameters are
 * accepted but never read, so they have no counterpart here.
 */
#ifndef GH_INSIGHTS_H
#define GH_INSIGHTS_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_analysis.h"
#include "gh_clinvar.h"
#include "gh_mem.h"
#include "gh_scorers.h"

typedef struct {
    const char *title;
    const char *finding;
    const char *reference;
    const char *practical;
} gh_insight_entry;

typedef struct {
    const char *gene;
    const char *status;
    const gh_insight_entry *entries;
    size_t nentries;
} gh_single_gene_insight;

extern const gh_single_gene_insight GH_SINGLE_GENE_INSIGHTS[];
extern const size_t GH_SINGLE_GENE_INSIGHT_COUNT;

typedef struct {
    const char *gene;
    const char *const *statuses;   /* any one of these qualifies */
    size_t nstatuses;
} gh_gene_match;

typedef struct {
    const char *id;
    const char *title;
    const gh_gene_match *required_genes;
    size_t nrequired_genes;
    const gh_gene_match *optional_genes;
    size_t noptional_genes;
    int min_matches;
    const char *narrative;         /* may contain {matched_count} */
    const char *practical;
    const char *const *references;
    size_t nreferences;
} gh_narrative_pattern;

extern const gh_narrative_pattern GH_NARRATIVES[];
extern const size_t GH_NARRATIVE_COUNT;

/* ------------------------------------------------------------------ */

#define GH_MAX_HIGHLIGHTS 5
#define GH_MAX_MATCHED_GENES 12

typedef struct {
    const gh_insight_entry *entry;
    const char *gene;
    const char *status;
    int magnitude;
} gh_single_gene_result;

typedef struct {
    const char *id;
    const char *title;
    const char *matched_genes[GH_MAX_MATCHED_GENES];
    size_t nmatched;
    const char *narrative;         /* placeholder already substituted */
    const char *practical;
    const char *const *references;
    size_t nreferences;
} gh_narrative_result;

typedef struct {
    const char *title;
    const char *detail;
    const char *type;              /* protective / clinical / pharmacogenomic / lifestyle */
} gh_highlight;

typedef struct {
    const char *gene;
    const char *status;
    const char *title;
    const char *finding;
    const char *reference;
} gh_protective_finding;

typedef struct {
    gh_single_gene_result *single_gene;
    size_t nsingle_gene;

    gh_narrative_result *narratives;
    size_t nnarratives;

    gh_highlight highlights[GH_MAX_HIGHLIGHTS];
    size_t nhighlights;

    gh_protective_finding *protective;
    size_t nprotective;
} gh_insights;

typedef struct {
    const gh_analysis *analysis;
    const gh_apoe_result *apoe;
    const gh_star_result *stars;
    size_t nstars;
    const gh_clinvar_result *clinvar;
} gh_insight_inputs;

void gh_generate_insights(gh_insights *out, gh_arena *arena,
                          const gh_insight_inputs *in);

#endif /* GH_INSIGHTS_H */

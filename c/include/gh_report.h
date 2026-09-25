/* HTML report generation.
 *
 * Mirrors genetic_health/reports/enhanced_html.py. The palette, paper
 * references, pathway groupings, ELI5 text and the HTML/CSS/JS shell are
 * generated from the Python by tools/gen_report_data.py, so the two reports
 * share one source of truth for their presentation layer.
 */
#ifndef GH_REPORT_H
#define GH_REPORT_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_analysis.h"
#include "gh_ancestry.h"
#include "gh_clinvar.h"
#include "gh_dependent_profiles.h"
#include "gh_dosing.h"
#include "gh_epistasis.h"
#include "gh_insights.h"
#include "gh_mem.h"
#include "gh_profiles.h"
#include "gh_prs.h"
#include "gh_quality.h"
#include "gh_recommendations.h"
#include "gh_scorers.h"
#include "gh_traits.h"

/* ------------------------------------------------------------------ */
/* Generated data                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *red;
    const char *orange;
    const char *amber;
    const char *green;
    const char *blue;
    const char *purple;
    const char *slate;
} gh_palette;

typedef struct {
    const char *rsid;
    const char *title;
    const char *pmid;
    int year;
} gh_paper_ref;

typedef struct {
    const char *name;
    const char *const *genes;
    size_t ngenes;
} gh_pathway;

#ifndef GH_KV_DEFINED
#define GH_KV_DEFINED
typedef struct {
    const char *key;
    const char *value;
} gh_kv;
#endif

extern const gh_palette GH_C;
extern const char *const GH_CHART_COLORS[];
extern const size_t GH_CHART_COLOR_COUNT;

extern const gh_paper_ref GH_PAPER_REFS[];
extern const size_t GH_PAPER_REF_COUNT;

extern const gh_pathway GH_PATHWAYS[];
extern const size_t GH_PATHWAY_COUNT;

extern const gh_kv GH_ELI5_GENES[];
extern const size_t GH_ELI5_GENE_COUNT;
extern const gh_kv GH_ELI5_CONDITIONS[];
extern const size_t GH_ELI5_CONDITION_COUNT;

/* The template is stored as NULL-terminated chunks because C99 only
 * guarantees 4095 characters per string literal. */
extern const char *const GH_HTML_TEMPLATE_PARTS[];
extern const char *const GH_TEMPLATE_SLOTS[];

/* Join the template chunks into one arena-held string. */
const char *gh_report_template(gh_arena *arena);

/* ------------------------------------------------------------------ */
/* Fragments                                                           */
/* ------------------------------------------------------------------ */

/* Links to dbSNP/ClinVar/SNPedia/PharmGKB. Appends nothing when `rsid` is
 * not an rsID, matching db_links_html. */
void gh_report_db_links(gh_strbuf *sb, const char *rsid);

/* Paper references for an rsID; appends nothing when there are none. */
void gh_report_paper_refs(gh_strbuf *sb, const char *rsid);

/* Plain-language gene explanation, or "" when the gene has none. */
const char *gh_report_eli5_gene(const char *gene);

/* Lower-case, underscores to spaces, first letter of each word upper-cased --
 * Python's `status.replace("_", " ").title()`. Result is arena-allocated. */
const char *gh_report_titlecase(gh_arena *arena, const char *status);

/* Magnitude bucket name: "high" (>=3), "mod" (2), "low" (1), "info" (0). */
const char *gh_report_mag_class(int magnitude);

/* ------------------------------------------------------------------ */
/* Text helpers mirroring the Python's                                */
/* ------------------------------------------------------------------ */

/* _clean_condition: the readable name out of ClinVar's pipe-separated
 * trait text. Arena-allocated. */
const char *gh_report_clean_condition(gh_arena *arena, const char *raw);

/* _clean_why: _clean_condition when the text has a "|", then duplicate
 * "; "-separated phrases removed. */
const char *gh_report_clean_why(gh_arena *arena, const char *text);

/* repr(float) for the values the report prints (odds ratios, scores):
 * shortest round-tripping digits, always with a decimal point. */
const char *gh_report_pyfloat(gh_arena *arena, double x);

/* The nutrigenomics supplement_priorities list: high/moderate needs in
 * profile order, stable-sorted high first. Returns how many were written. */
size_t gh_report_supplement_priorities(const gh_nutrigenomics_result *nut,
                                       const gh_nutrient_need **out, size_t cap);

/* ------------------------------------------------------------------ */
/* Charts                                                              */
/* ------------------------------------------------------------------ */

/* svg_impact_bar and svg_category_donut are defined in the Python but no
 * section calls them; kept for parity, not rendered. */
void gh_report_svg_impact_bar(gh_strbuf *sb, const gh_analysis *a);
void gh_report_svg_category_donut(gh_strbuf *sb, gh_arena *arena,
                                  const gh_analysis *a);

void gh_report_svg_metabolism_gauge(gh_strbuf *sb, const char *label,
                                    double level, const char *color);
void gh_report_svg_prs_gauge(gh_strbuf *sb, gh_arena *arena, const char *label,
                             double percentile, const char *category);
void gh_report_svg_risk_heatmap(gh_strbuf *sb, gh_arena *arena,
                                const gh_prs_result *prs, size_t nprs);
void gh_report_svg_ancestry_donut(gh_strbuf *sb, gh_arena *arena,
                                  const gh_ancestry_result *anc);

/* ------------------------------------------------------------------ */
/* Whole report                                                        */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *subject_name;   /* may be NULL */
    const char *generated_date; /* NULL for the current local date */
} gh_report_opts;

/* Everything the sections read. Any pointer may be NULL, in which case the
 * section renders what the Python renders for missing data. */
typedef struct {
    const gh_apoe_result *apoe;
    const gh_blood_result *blood;
    const gh_mt_result *mt;
    const gh_star_result *stars;
    size_t nstars;
    const gh_trait_result *traits;
    size_t ntraits;
    const gh_ancestry_result *ancestry;
    const gh_prs_result *prs;
    size_t nprs;
    const gh_clinvar_result *clinvar;
    const gh_acmg_result *acmg;
    const gh_carrier_result *carriers;
    const gh_epi_result *epistasis;
    size_t nepistasis;
    const gh_recommendations *recs;
    const gh_insights *insights;
    const gh_sleep_result *sleep;
    const gh_nutrigenomics_result *nutrition;
    const gh_mental_result *mental;
    const gh_longevity_result *longevity;
    const gh_dosing_result *dosing;
    const gh_polypharmacy_result *polypharmacy;
    const gh_preventive_result *preventive;
    const gh_pain_result *pain;
    const gh_histamine_result *histamine;
    const gh_thyroid_result *thyroid;
    const gh_hormone_result *hormone;
    const gh_eye_result *eye;
    const gh_alcohol_result *alcohol;
} gh_report_scorers;

/* ------------------------------------------------------------------ */
/* Sections, one per build_* in enhanced_html.py                       */
/* ------------------------------------------------------------------ */

void gh_report_key_findings(gh_strbuf *sb, gh_arena *arena, const gh_analysis *a,
                            const gh_report_scorers *sc);
void gh_report_action_plan(gh_strbuf *sb, gh_arena *arena, const gh_report_scorers *sc);
void gh_report_drug_guide(gh_strbuf *sb, gh_arena *arena, const gh_analysis *a,
                          const gh_report_scorers *sc);
void gh_report_disease_risk(gh_strbuf *sb, gh_arena *arena, const gh_report_scorers *sc);
void gh_report_body_profile(gh_strbuf *sb, gh_arena *arena, const gh_analysis *a,
                            const gh_report_scorers *sc);
void gh_report_mental_health(gh_strbuf *sb, gh_arena *arena, const gh_report_scorers *sc);
void gh_report_clinical_detail(gh_strbuf *sb, gh_arena *arena, const gh_analysis *a,
                               const gh_report_scorers *sc);
void gh_report_ancestry(gh_strbuf *sb, gh_arena *arena, const gh_report_scorers *sc);
void gh_report_nutrigenomics(gh_strbuf *sb, gh_arena *arena, const gh_report_scorers *sc);
void gh_report_quality(gh_strbuf *sb, const gh_quality *q);
void gh_report_doctor_card(gh_strbuf *sb, gh_arena *arena, const gh_report_scorers *sc);
void gh_report_references(gh_strbuf *sb, const gh_analysis *a);

/* Render the full report into `sb`. `q` and `sc` may be NULL. */
void gh_report_render(gh_strbuf *sb, gh_arena *arena, const gh_analysis *a,
                      const gh_quality *q, const gh_report_scorers *sc,
                      const gh_report_opts *opts);

/* Substitute {slot} placeholders in `template` from parallel name/value
 * arrays. Unknown placeholders are left verbatim, as a missing section should
 * be visible rather than silently blank. Exposed for testing. */
void gh_report_fill_template(gh_strbuf *sb, const char *template_text,
                             const char *const *names,
                             const char *const *values, size_t count);

#endif /* GH_REPORT_H */

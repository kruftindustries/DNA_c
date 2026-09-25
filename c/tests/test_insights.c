/* Unit tests for the insights engine.
 *
 * The differential sweep compares output against the Python across crafted
 * genomes; these cover the structural rules and the precedence decisions.
 */

#include <stdio.h>
#include <string.h>

#include "gh_insights.h"
#include "gh_mem.h"

static int failures;
static int checks;

static void check(int cond, const char *what)
{
    checks++;
    if (!cond) {
        failures++;
        fprintf(stderr, "  FAIL: %s\n", what);
    }
}

static void check_str(const char *got, const char *want, const char *what)
{
    checks++;
    if (!got || strcmp(got, want) != 0) {
        failures++;
        fprintf(stderr, "  FAIL: %s\n    got  %s\n    want %s\n",
                what, got ? got : "(null)", want);
    }
}

static void check_has(const char *haystack, const char *needle, const char *what)
{
    checks++;
    if (!haystack || !strstr(haystack, needle)) {
        failures++;
        fprintf(stderr, "  FAIL: %s\n    missing: %s\n", what, needle);
    }
}

static void check_lacks(const char *haystack, const char *needle,
                        const char *what)
{
    checks++;
    if (haystack && strstr(haystack, needle)) {
        failures++;
        fprintf(stderr, "  FAIL: %s\n    unexpected: %s\n", what, needle);
    }
}

static void test_empty(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_analysis analysis = {0};
    gh_insight_inputs in = {&analysis, NULL, NULL, 0, NULL};

    gh_insights out;
    gh_generate_insights(&out, a, &in);
    check(out.nsingle_gene == 0, "no single-gene insights");
    check(out.nnarratives == 0, "no narratives");
    check(out.nhighlights == 0, "no highlights");
    check(out.nprotective == 0, "no protective findings");

    gh_generate_insights(&out, a, &(gh_insight_inputs){NULL, NULL, NULL, 0, NULL});
    check(out.nsingle_gene == 0, "NULL analysis is safe");

    gh_arena_free(a);
}

static void test_single_gene_matching(void)
{
    gh_arena *a = gh_arena_new(1 << 16);

    /* Take a gene/status pair straight from the table so the match is real. */
    const gh_single_gene_insight *target = &GH_SINGLE_GENE_INSIGHTS[0];
    gh_finding findings[] = {
        {"rs1", target->gene, "Cat", "AA", target->status, "d", "", NULL, 4},
    };
    gh_analysis analysis = {0};
    analysis.findings = findings;
    analysis.nfindings = 1;

    gh_insight_inputs in = {&analysis, NULL, NULL, 0, NULL};
    gh_insights out;
    gh_generate_insights(&out, a, &in);

    check(out.nsingle_gene == target->nentries, "entries matched");
    if (out.nsingle_gene) {
        check_str(out.single_gene[0].gene, target->gene, "gene carried");
        check_str(out.single_gene[0].status, target->status, "status carried");
        check(out.single_gene[0].magnitude == 4, "magnitude carried");
        check(out.single_gene[0].entry->title != NULL, "entry has a title");
    }

    /* A different status for the same gene must not match. */
    findings[0].status = "definitely_not_a_real_status";
    gh_generate_insights(&out, a, &in);
    check(out.nsingle_gene == 0, "status must match exactly");

    gh_arena_free(a);
}

static void test_single_gene_sorted_by_magnitude(void)
{
    gh_arena *a = gh_arena_new(1 << 16);

    /* Two matching genes with different magnitudes. */
    const gh_single_gene_insight *first = &GH_SINGLE_GENE_INSIGHTS[0];
    const gh_single_gene_insight *second = NULL;
    for (size_t i = 1; i < GH_SINGLE_GENE_INSIGHT_COUNT; i++)
        if (strcmp(GH_SINGLE_GENE_INSIGHTS[i].gene, first->gene) != 0) {
            second = &GH_SINGLE_GENE_INSIGHTS[i];
            break;
        }
    check(second != NULL, "two distinct genes available in the table");
    if (!second) {
        gh_arena_free(a);
        return;
    }

    gh_finding findings[] = {
        {"rs1", first->gene, "Cat", "AA", first->status, "d", "", NULL, 1},
        {"rs2", second->gene, "Cat", "AA", second->status, "d", "", NULL, 5},
    };
    gh_analysis analysis = {0};
    analysis.findings = findings;
    analysis.nfindings = 2;

    gh_insights out;
    gh_generate_insights(&out, a, &(gh_insight_inputs){&analysis, NULL, NULL, 0, NULL});
    check(out.nsingle_gene >= 2, "both genes matched");
    if (out.nsingle_gene >= 2)
        check(out.single_gene[0].magnitude >= out.single_gene[1].magnitude,
              "sorted by descending magnitude");

    gh_arena_free(a);
}

static void test_narrative_thresholds(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    const gh_narrative_pattern *pattern = &GH_NARRATIVES[0];

    /* Only the first required gene: below min_matches, so nothing fires. */
    gh_finding one[] = {
        {"rs1", pattern->required_genes[0].gene, "Cat", "AA",
         pattern->required_genes[0].statuses[0], "d", "", NULL, 3},
    };
    gh_analysis analysis = {0};
    analysis.findings = one;
    analysis.nfindings = 1;

    gh_insights out;
    gh_generate_insights(&out, a, &(gh_insight_inputs){&analysis, NULL, NULL, 0, NULL});
    bool fired = false;
    for (size_t i = 0; i < out.nnarratives; i++)
        if (strcmp(out.narratives[i].id, pattern->id) == 0)
            fired = true;
    check(pattern->min_matches <= 1 || !fired,
          "a single gene does not reach the threshold");

    /* All required genes present. */
    gh_finding all[GH_MAX_MATCHED_GENES];
    size_t n = 0;
    for (size_t i = 0; i < pattern->nrequired_genes && n < GH_MAX_MATCHED_GENES; i++)
        all[n++] = (gh_finding){"rs1", pattern->required_genes[i].gene, "Cat",
                                "AA", pattern->required_genes[i].statuses[0],
                                "d", "", NULL, 3};
    analysis.findings = all;
    analysis.nfindings = n;
    gh_generate_insights(&out, a, &(gh_insight_inputs){&analysis, NULL, NULL, 0, NULL});

    const gh_narrative_result *result = NULL;
    for (size_t i = 0; i < out.nnarratives; i++)
        if (strcmp(out.narratives[i].id, pattern->id) == 0)
            result = &out.narratives[i];

    if ((int)n >= pattern->min_matches) {
        check(result != NULL, "narrative fires once the threshold is met");
        if (result) {
            check(result->nmatched == n, "matched gene count");
            check_lacks(result->narrative, "{matched_count}",
                        "placeholder is substituted");
            check(result->nreferences > 0, "references carried through");
        }
    }

    gh_arena_free(a);
}

static void test_optional_genes_need_a_required_one(void)
{
    gh_arena *a = gh_arena_new(1 << 16);

    /* Find a pattern with enough optional genes to reach its threshold on
     * their own; it must still not fire without a required gene. */
    const gh_narrative_pattern *pattern = NULL;
    for (size_t i = 0; i < GH_NARRATIVE_COUNT; i++)
        if ((int)GH_NARRATIVES[i].noptional_genes >= GH_NARRATIVES[i].min_matches) {
            pattern = &GH_NARRATIVES[i];
            break;
        }
    if (!pattern) {
        gh_arena_free(a);
        return;
    }

    gh_finding optional_only[GH_MAX_MATCHED_GENES];
    size_t n = 0;
    for (size_t i = 0; i < pattern->noptional_genes && n < GH_MAX_MATCHED_GENES; i++)
        optional_only[n++] = (gh_finding){
            "rs1", pattern->optional_genes[i].gene, "Cat", "AA",
            pattern->optional_genes[i].statuses[0], "d", "", NULL, 3};

    gh_analysis analysis = {0};
    analysis.findings = optional_only;
    analysis.nfindings = n;

    gh_insights out;
    gh_generate_insights(&out, a, &(gh_insight_inputs){&analysis, NULL, NULL, 0, NULL});
    for (size_t i = 0; i < out.nnarratives; i++)
        check(strcmp(out.narratives[i].id, pattern->id) != 0,
              "optional genes alone do not fire a narrative");

    gh_arena_free(a);
}

static void test_highlight_precedence_and_cap(void)
{
    gh_arena *a = gh_arena_new(1 << 16);

    gh_apoe_result apoe = {
        .apoe_type = "e2/e3", .risk_level = "reduced", .alzheimer_or = 0.6,
        .has_or = true, .description = "d", .confidence = "high",
        .rs429358 = "TT", .rs7412 = "TC",
    };
    gh_star_result stars[] = {
        {"CYP2C19", "*1/*2", "poor", "high", "note", 8, 8, 1.0},
        {"CYP2D6", "*1/*4", "intermediate", "high", "note", 6, 6, 1.0},
        {"TPMT", "*1/*1", "normal", "high", "note", 5, 5, 1.0},
    };
    gh_finding findings[] = {
        {"rs1", "GENE1", "Cat", "AA", "reduced", "d1", "", NULL, 5},
        {"rs2", "GENE2", "Cat", "AA", "reduced", "d2", "", NULL, 4},
        {"rs3", "GENE3", "Cat", "AA", "reduced", "d3", "", NULL, 3},
        {"rs4", "GENE4", "Cat", "AA", "reduced", "d4", "", NULL, 3},
    };
    gh_analysis analysis = {0};
    analysis.findings = findings;
    analysis.nfindings = 4;

    gh_insights out;
    gh_generate_insights(&out, a, &(gh_insight_inputs){
        &analysis, &apoe, stars, 3, NULL});

    check(out.nhighlights == GH_MAX_HIGHLIGHTS, "highlights are capped at five");
    check_str(out.highlights[0].type, "protective", "APOE e2 leads");
    check_has(out.highlights[0].title, "e2/e3", "APOE haplotype named");

    /* Only actionable phenotypes appear; a normal metabolizer does not. */
    bool has_normal = false;
    for (size_t i = 0; i < out.nhighlights; i++)
        if (strstr(out.highlights[i].title, "TPMT"))
            has_normal = true;
    check(!has_normal, "a normal metabolizer is not highlighted");

    /* Phenotype is title-cased in the highlight. */
    bool found_poor = false;
    for (size_t i = 0; i < out.nhighlights; i++)
        if (strstr(out.highlights[i].title, "Poor Metabolizer"))
            found_poor = true;
    check(found_poor, "phenotype is title-cased");

    /* An e4-only genome contributes no APOE highlight. */
    gh_apoe_result e4 = apoe;
    e4.apoe_type = "e3/e4";
    gh_generate_insights(&out, a, &(gh_insight_inputs){
        &analysis, &e4, NULL, 0, NULL});
    for (size_t i = 0; i < out.nhighlights; i++)
        check_lacks(out.highlights[i].title, "Alzheimer's protection",
                    "no protective APOE highlight without e2");

    gh_arena_free(a);
}

static void test_star_alleles_override_findings(void)
{
    gh_arena *a = gh_arena_new(1 << 16);

    /* A lifestyle finding and a star-allele call for the same gene: the
     * star-allele phenotype wins. */
    gh_finding findings[] = {
        {"rs1", "CYP2C19", "Drug", "AA", "normal", "d", "", NULL, 1},
    };
    gh_star_result stars[] = {
        {"CYP2C19", "*2/*2", "poor", "high", "note", 8, 8, 1.0},
    };
    gh_analysis analysis = {0};
    analysis.findings = findings;
    analysis.nfindings = 1;

    gh_insights out;
    gh_generate_insights(&out, a, &(gh_insight_inputs){
        &analysis, NULL, stars, 1, NULL});

    bool poor_highlighted = false;
    for (size_t i = 0; i < out.nhighlights; i++)
        if (strstr(out.highlights[i].title, "Poor Metabolizer"))
            poor_highlighted = true;
    check(poor_highlighted, "the star-allele phenotype is what gets reported");

    gh_arena_free(a);
}

static void test_table_integrity(void)
{
    for (size_t i = 0; i < GH_SINGLE_GENE_INSIGHT_COUNT; i++) {
        const gh_single_gene_insight *s = &GH_SINGLE_GENE_INSIGHTS[i];
        check(s->gene && *s->gene, "insight names a gene");
        check(s->status && *s->status, "insight names a status");
        check(s->nentries > 0, "insight has entries");
        for (size_t e = 0; e < s->nentries; e++) {
            check(s->entries[e].title && *s->entries[e].title, "entry title");
            check(s->entries[e].finding != NULL, "entry finding");
            check(s->entries[e].reference != NULL, "entry reference");
            check(s->entries[e].practical != NULL, "entry practical advice");
        }
    }
    for (size_t i = 0; i < GH_NARRATIVE_COUNT; i++) {
        const gh_narrative_pattern *p = &GH_NARRATIVES[i];
        check(p->id && *p->id, "narrative has an id");
        check(p->min_matches > 0, "narrative has a threshold");
        check(p->nrequired_genes > 0, "narrative requires at least one gene");
        check(p->nreferences > 0, "narrative cites references");
        check(p->nrequired_genes + p->noptional_genes <= GH_MAX_MATCHED_GENES,
              "narrative fits the matched-gene limit");
    }
}

int main(void)
{
    test_empty();
    test_single_gene_matching();
    test_single_gene_sorted_by_magnitude();
    test_narrative_thresholds();
    test_optional_genes_need_a_required_one();
    test_highlight_precedence_and_cap();
    test_star_alleles_override_findings();
    test_table_integrity();

    printf("%s: %d checks, %d failures\n",
           failures ? "FAILED" : "ok", checks, failures);
    return failures ? 1 : 0;
}

/* Unit tests for gene-gene interaction evaluation. */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "gh_analysis.h"
#include "gh_epistasis.h"
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

static void check_close(double got, double want, double tol, const char *what)
{
    checks++;
    if (!(fabs(got - want) <= tol)) {
        failures++;
        fprintf(stderr, "  FAIL: %s\n    got %.4f want %.4f\n", what, got, want);
    }
}

static bool risk_is_known(const char *risk)
{
    return strcmp(risk, "high") == 0 || strcmp(risk, "moderate") == 0 ||
           strcmp(risk, "low") == 0;
}

static void test_severity_weights(void)
{
    /* Weights come from the generated table; "slow" and "high" are listed
     * twice in the Python and the surviving value is 1.0. */
    check_close(gh_status_severity_for("slow"), 1.0, 1e-12, "slow weight");
    check_close(gh_status_severity_for("high"), 1.0, 1e-12, "high weight");
    check_close(gh_status_severity_for("reduced"), 0.5, 1e-12, "reduced weight");
    check_close(gh_status_severity_for("severely_reduced"), 1.0, 1e-12,
                "severely reduced weight");
    /* Anything unlisted counts as moderate rather than zero. */
    check_close(gh_status_severity_for("not_a_status"), 0.5, 1e-12,
                "unlisted status defaults to 0.5");
}

/* Build an analysis from alternating gene/status/magnitude triples. */
static void make_analysis(gh_analysis *a, gh_finding *findings, size_t n)
{
    memset(a, 0, sizeof(*a));
    a->findings = findings;
    a->nfindings = n;
}

static const gh_epi_result *find_by_id(const gh_epi_result *r, size_t n,
                                       const char *id)
{
    for (size_t i = 0; i < n; i++)
        if (strcmp(r[i].id, id) == 0)
            return &r[i];
    return NULL;
}

static void test_no_matches(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    gh_analysis a;

    /* Nothing at all. */
    make_analysis(&a, NULL, 0);
    size_t n = 0;
    gh_epi_result *r = gh_evaluate_epistasis(arena, &a, &n);
    check(n == 0, "no findings means no interactions");
    check(r != NULL, "a valid pointer is returned even when empty");

    /* Genes present but with statuses that do not qualify. */
    gh_finding findings[] = {
        {"rs1", "MTHFR", "Methylation", "GG", "normal", "d", "", NULL, 0},
        {"rs2", "COMT", "Neuro", "GG", "fast", "d", "", NULL, 0},
    };
    make_analysis(&a, findings, 2);
    gh_evaluate_epistasis(arena, &a, &n);
    check(n == 0, "non-qualifying statuses do not match");

    /* Only one gene of a pair qualifies. */
    gh_finding half[] = {
        {"rs1", "MTHFR", "Methylation", "AA", "reduced", "d", "", NULL, 3},
    };
    make_analysis(&a, half, 1);
    gh_evaluate_epistasis(arena, &a, &n);
    check(n == 0, "a partial match does not count");

    gh_arena_free(arena);
}

static void test_matching_and_severity(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    gh_analysis a;
    size_t n = 0;

    /* MTHFR reduced + COMT slow is a defined high-risk interaction. */
    gh_finding findings[] = {
        {"rs1801133", "MTHFR", "Methylation", "TT", "reduced", "d", "", NULL, 6},
        {"rs4680", "COMT", "Neuro", "AA", "slow", "d", "", NULL, 6},
    };
    make_analysis(&a, findings, 2);
    gh_epi_result *r = gh_evaluate_epistasis(arena, &a, &n);
    check(n >= 1, "the MTHFR/COMT interaction matches");

    const gh_epi_result *mc = find_by_id(r, n, "mthfr_comt_methylation");
    check(mc != NULL, "interaction found by id");
    if (mc) {
        check(mc->ngenes == 2, "two genes involved");
        check(mc->nactions > 0, "actions are carried through");
        check(mc->effect && *mc->effect, "effect text present");
        check(mc->mechanism && *mc->mechanism, "mechanism text present");

        /* Both genes at magnitude 6: mag_factor 1.0, so each severity is
         * weight * (0.5 + 0.5) = weight. MTHFR reduced is 0.5, COMT slow
         * is 1.0, geometric mean sqrt(0.5) = 0.707. */
        check_close(mc->severity_score, 0.71, 0.005,
                    "severity is the geometric mean, rounded");
        check_str(mc->risk_level, "high", "strong match keeps its high risk");
    }

    /* The same pair at low magnitude scores lower and can be demoted. */
    gh_finding weak[] = {
        {"rs1801133", "MTHFR", "Methylation", "TT", "reduced", "d", "", NULL, 0},
        {"rs4680", "COMT", "Neuro", "AA", "slow", "d", "", NULL, 0},
    };
    make_analysis(&a, weak, 2);
    r = gh_evaluate_epistasis(arena, &a, &n);
    const gh_epi_result *weak_mc = find_by_id(r, n, "mthfr_comt_methylation");
    check(weak_mc != NULL, "weak interaction still matches");
    if (weak_mc) {
        /* mag_factor 0, so each severity halves: 0.25 and 0.5,
         * geometric mean sqrt(0.125) = 0.354. */
        check_close(weak_mc->severity_score, 0.35, 0.005,
                    "low magnitude lowers the severity");
        check(weak_mc->severity_score < 0.5,
              "weak match scores below a strong one");
    }

    gh_arena_free(arena);
}

static void test_ordering(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    gh_analysis a;
    size_t n = 0;

    /* Enough qualifying statuses to trigger several interactions at once. */
    gh_finding findings[] = {
        {"rs1801133", "MTHFR", "Methylation", "TT", "severely_reduced", "d", "", NULL, 6},
        {"rs4680", "COMT", "Neuro", "AA", "slow", "d", "", NULL, 6},
        {"rs762551", "CYP1A2", "Drug", "CC", "slow", "d", "", NULL, 6},
        {"rs5751876", "ADORA2A", "Neuro", "TT", "anxiety_prone", "d", "", NULL, 6},
    };
    make_analysis(&a, findings, 4);
    gh_epi_result *r = gh_evaluate_epistasis(arena, &a, &n);
    check(n >= 2, "several interactions match");

    /* High risk sorts before moderate, and within a level the stronger
     * interaction comes first. */
    for (size_t i = 1; i < n; i++) {
        int prev = strcmp(r[i - 1].risk_level, "high") == 0 ? 0
                 : strcmp(r[i - 1].risk_level, "moderate") == 0 ? 1 : 2;
        int cur = strcmp(r[i].risk_level, "high") == 0 ? 0
                : strcmp(r[i].risk_level, "moderate") == 0 ? 1 : 2;
        check(prev <= cur, "sorted by risk level");
        if (prev == cur)
            check(r[i - 1].severity_score >= r[i].severity_score,
                  "sorted by descending severity within a level");
    }

    gh_arena_free(arena);
}

static void test_model_integrity(void)
{
    for (size_t m = 0; m < GH_EPISTASIS_MODEL_COUNT; m++) {
        const gh_epi_model *model = &GH_EPISTASIS_MODELS[m];
        check(model->id && *model->id, "model has an id");
        check(model->name && *model->name, "model has a name");
        check(model->nconditions > 0, "model has conditions");

        for (size_t c = 0; c < model->nconditions; c++) {
            const gh_epi_condition *cond = &model->conditions[c];
            check(cond->nrequired > 0, "condition requires at least one gene");
            check(cond->nrequired <= GH_EPI_MAX_GENES,
                  "condition fits the per-interaction gene limit");
            check(cond->effect && *cond->effect, "condition has an effect");
            check(cond->mechanism && *cond->mechanism,
                  "condition has a mechanism");
            check(risk_is_known(cond->risk_level), "risk level is recognised");

            for (size_t r = 0; r < cond->nrequired; r++) {
                check(cond->required[r].gene && *cond->required[r].gene,
                      "requirement names a gene");
                check(cond->required[r].nstatuses > 0,
                      "requirement lists qualifying statuses");
            }
        }
    }
}

int main(void)
{
    test_severity_weights();
    test_no_matches();
    test_matching_and_severity();
    test_ordering();
    test_model_integrity();

    printf("%s: %d checks, %d failures\n",
           failures ? "FAILED" : "ok", checks, failures);
    return failures ? 1 : 0;
}

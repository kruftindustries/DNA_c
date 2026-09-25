/* Unit tests for the six genome-only health profile modules.
 *
 * The exhaustive differential sweep already compares every input state
 * against the Python. These cover the structural guarantees and the
 * boundaries that are easy to get wrong when the Python changes.
 */

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "gh_genome.h"
#include "gh_mem.h"
#include "gh_profiles.h"

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

static void load_fixture(gh_genome *g, gh_arena *arena, ...)
{
    static int counter;
    char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_profile_fixture_%d.txt", counter++);

    FILE *f = fopen(path, "wb");
    if (!f) {
        failures++;
        return;
    }
    fputs("# fixture\n", f);

    va_list ap;
    va_start(ap, arena);
    int line = 0;
    for (;;) {
        const char *rsid = va_arg(ap, const char *);
        if (!rsid)
            break;
        fprintf(f, "%s\t1\t%d\t%s\n", rsid, 1000 + line++,
                va_arg(ap, const char *));
    }
    va_end(ap);
    fclose(f);

    if (!gh_genome_load(g, arena, path, NULL))
        failures++;
}

/* Every profile must report cleanly with nothing to go on. */
static void test_empty_genome(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    load_fixture(&g, a, "rs00000000", "AA", NULL);

    gh_histamine_result h;
    gh_profile_histamine(&h, a, &g);
    check(h.base.nhits == 0, "histamine finds nothing");
    check(h.base.markers_tested == 3, "histamine tests three markers");
    check_str(h.risk_level, "unknown", "histamine risk unknown");
    check_has(h.base.summary, "No histamine", "histamine says so");
    check(h.nfoods == 8, "histamine still lists foods to watch");

    gh_alcohol_result al;
    gh_profile_alcohol(&al, a, &g);
    check_str(al.metabolism_speed, "normal", "alcohol defaults to normal speed");
    check_str(al.flush_risk, "unknown", "alcohol flush unknown");
    check_str(al.cancer_risk, "average", "alcohol cancer risk average");
    check(al.base.nrecommendations == 0,
          "no alcohol recommendations without markers");

    gh_pain_result pain;
    gh_profile_pain(&pain, a, &g);
    check(pain.sensitivity_score == 50, "pain defaults to a neutral 50");
    check(pain.base.nrecommendations == 0, "no pain recommendations");

    gh_thyroid_result t;
    gh_profile_thyroid(&t, a, &g);
    check(t.ndomains == 3, "thyroid reports three domains");
    for (size_t i = 0; i < t.ndomains; i++)
        check_str(t.domains[i].level, "unknown", "thyroid domain unknown");

    gh_hormone_result hm;
    gh_profile_hormone(&hm, a, &g);
    check_str(hm.estrogen_level, "unknown", "hormone estrogen unknown");
    check_str(hm.overall, "unknown", "hormone overall unknown");

    gh_eye_result e;
    gh_profile_eye(&e, a, &g);
    check(e.nconditions == 2, "eye reports two conditions");
    check_str(e.conditions[0].level, "unknown", "eye glaucoma unknown");

    gh_arena_free(a);
}

static void test_histamine_bands(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_histamine_result r;

    /* Risk points accumulate across DAO and HNMT. */
    load_fixture(&g, a, "rs10156191", "TT", "rs1049793", "CC",
                 "rs11558538", "TT", NULL);
    gh_profile_histamine(&r, a, &g);
    check(r.base.nhits == 3, "all three markers found");
    check(r.total_risk >= 4, "homozygous risk alleles accumulate");
    check_str(r.risk_level, "elevated", "elevated risk band");
    check(r.base.nrecommendations == 4, "elevated gives four recommendations");

    /* No risk alleles at all. */
    load_fixture(&g, a, "rs10156191", "CC", "rs1049793", "TT",
                 "rs11558538", "CC", NULL);
    gh_profile_histamine(&r, a, &g);
    check(r.total_risk == 0, "no risk points");
    check_str(r.risk_level, "low", "low risk band");
    check(r.base.nrecommendations == 1, "low gives one recommendation");

    gh_arena_free(a);
}

static void test_alcohol_interaction(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_alcohol_result r;

    /* ALDH2 homozygous variant is the severe flush genotype. */
    load_fixture(&g, a, "rs671", "AA", NULL);
    gh_profile_alcohol(&r, a, &g);
    check_str(r.flush_risk, "severe", "severe flush detected");
    check_has(r.base.recommendations[0], "avoiding alcohol entirely",
              "severe flush advises avoidance");

    /* CYP2E1 can only slow a speed that is otherwise normal, never a fast
     * one -- a detail of the Python's ordering. */
    load_fixture(&g, a, "rs1229984", "TT", "rs2031920", "TT", NULL);
    gh_profile_alcohol(&r, a, &g);
    check_str(r.metabolism_speed, "fast",
              "a fast ADH1B result is not downgraded by CYP2E1");

    gh_arena_free(a);
}

static void test_pain_scoring(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_pain_result r;

    /* A single marker means the score is that marker's contribution. */
    load_fixture(&g, a, "rs1799971", "AA", NULL);
    gh_profile_pain(&r, a, &g);
    check(r.base.nhits == 1, "one marker found");
    check(r.sensitivity_score == 50, "single normal marker scores 50");
    check(r.base.nrecommendations >= 1, "a recommendation is still offered");

    /* The score is the mean of the markers present, not of all four. */
    load_fixture(&g, a, "rs1799971", "GG", NULL);
    gh_profile_pain(&r, a, &g);
    check(r.sensitivity_score == 20, "reduced opioid response scores 20");
    check_has(r.base.summary, "notably low", "low score summary");

    gh_arena_free(a);
}

static void test_eye_high_penetrance(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_eye_result r;

    /* MYOC is high-penetrance: one severe marker outranks the average, so a
     * benign marker in the same condition cannot dilute it. */
    load_fixture(&g, a, "rs74315329", "AA", "rs4236601", "GG", NULL);
    gh_profile_eye(&r, a, &g);
    check_str(gh_profile_domain_level(r.conditions, r.nconditions,
                                      "glaucoma_risk"),
              "high", "a severe marker sets the level regardless of average");
    check_has(r.base.summary, "urgent", "summary flags urgency");
    check_has(r.base.recommendations[0], "Urgent",
              "first recommendation is the urgent one");

    gh_arena_free(a);
}

static void test_hormone_pathways(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_hormone_result r;

    /* An estrogen-pathway effect must not leak into the androgen summary,
     * and vice versa. */
    load_fixture(&g, a, "rs4646", "TT", NULL);
    gh_profile_hormone(&r, a, &g);
    check_has(r.estrogen_level, "estrogen", "estrogen pathway summarised");
    check_str(r.androgen_level, "unknown",
              "androgen pathway stays unknown without its markers");

    load_fixture(&g, a, "rs523349", "GG", NULL);
    gh_profile_hormone(&r, a, &g);
    check_str(r.estrogen_level, "unknown",
              "estrogen pathway stays unknown without its markers");
    check(strcmp(r.androgen_level, "unknown") != 0,
          "androgen pathway summarised");

    gh_arena_free(a);
}

static void test_marker_table_integrity(void)
{
    const gh_profile_markers *all[] = {
        &GH_HISTAMINE_MARKERS, &GH_ALCOHOL_MARKERS, &GH_PAIN_MARKERS,
        &GH_THYROID_MARKERS, &GH_HORMONE_MARKERS, &GH_EYE_MARKERS,
    };
    for (size_t m = 0; m < sizeof(all) / sizeof(all[0]); m++) {
        check(all[m]->count > 0, "module has markers");
        check(all[m]->count <= GH_PROFILE_MAX_MARKERS,
              "module fits the marker limit");
        for (size_t i = 0; i < all[m]->count; i++) {
            const gh_profile_marker *marker = &all[m]->markers[i];
            check(strncmp(marker->rsid, "rs", 2) == 0, "marker names an rsID");
            check(marker->gene && *marker->gene, "marker names a gene");
            check(strchr("ACGT", marker->allele) != NULL,
                  "counted allele is a base");
            for (size_t o = 0; o < 3; o++) {
                check(marker->outcomes[o].label != NULL, "outcome has a label");
                check(marker->outcomes[o].description != NULL,
                      "outcome has a description");
            }
        }
    }
}

int main(void)
{
    test_empty_genome();
    test_histamine_bands();
    test_alcohol_interaction();
    test_pain_scoring();
    test_eye_high_penetrance();
    test_hormone_pathways();
    test_marker_table_integrity();

    printf("%s: %d checks, %d failures\n",
           failures ? "FAILED" : "ok", checks, failures);
    return failures ? 1 : 0;
}

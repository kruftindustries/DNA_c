/* Unit tests for the trait predictors. */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gh_genome.h"
#include "gh_mem.h"
#include "gh_traits.h"

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

/* ------------------------------------------------------------------ */
/* Fixtures                                                            */
/* ------------------------------------------------------------------ */

static void load_fixture(gh_genome *g, gh_arena *arena, ...)
{
    static int counter;
    char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_traits_fixture_%d.txt", counter++);

    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "  FAIL: cannot write %s\n", path);
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
        const char *genotype = va_arg(ap, const char *);
        fprintf(f, "%s\t1\t%d\t%s\n", rsid, 1000 + line++, genotype);
    }
    va_end(ap);
    fclose(f);

    if (!gh_genome_load(g, arena, path, NULL)) {
        fprintf(stderr, "  FAIL: cannot load %s\n", path);
        failures++;
    }
}

/* Predict from a fixture and return the named trait. */
static const gh_trait_result *predict_one(gh_trait_result *results,
                                          gh_arena *arena, const gh_genome *g,
                                          const char *key)
{
    gh_predict_traits(results, arena, g);
    return gh_trait_find(results, key);
}

/* ------------------------------------------------------------------ */
/* Structure                                                           */
/* ------------------------------------------------------------------ */

static void test_all_traits_reported(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_trait_result results[GH_TRAIT_COUNT];

    /* With no relevant SNPs at all, every trait must still be reported as
     * Unknown rather than dropped. */
    load_fixture(&g, a, "rs00000000", "AA", NULL);
    gh_predict_traits(results, a, &g);

    size_t reported = 0;
    for (size_t i = 0; i < GH_TRAIT_COUNT; i++) {
        if (!results[i].key)
            continue;
        reported++;
        check(results[i].prediction != NULL, "prediction set");
        check(results[i].confidence != NULL, "confidence set");
        check(results[i].description != NULL, "description set");
        check(results[i].nsnps == 0, "no SNPs listed when none present");
    }
    check(reported == GH_TRAIT_COUNT, "every trait is reported");

    /* The order must match the Python's predict_traits() dict order. */
    for (size_t i = 0; GH_TRAIT_ORDER[i]; i++)
        check_str(results[i].key, GH_TRAIT_ORDER[i], "trait order matches");

    check(gh_trait_find(results, "eye_color") != NULL, "lookup by key works");
    check(gh_trait_find(results, "not_a_trait") == NULL, "unknown key is NULL");

    gh_arena_free(a);
}

/* ------------------------------------------------------------------ */
/* Eye colour                                                          */
/* ------------------------------------------------------------------ */

static void test_eye_color(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_trait_result results[GH_TRAIT_COUNT];
    const gh_trait_result *r;

    /* rs12913832 GG is the strongest brown-eye signal. */
    load_fixture(&g, a, "rs12913832", "GG", NULL);
    r = predict_one(results, a, &g, "eye_color");
    check_str(r->prediction, "Likely brown", "HERC2 GG predicts brown");
    check_str(r->confidence, "high", "brown call is high confidence");
    check(r->nsnps == 1, "one SNP listed");
    check_has(r->snps_used[0], "rs12913832 (HERC2): GG", "SNP line formatted");

    /* AA predicts blue, unless OCA2 shifts it toward green. */
    load_fixture(&g, a, "rs12913832", "AA", NULL);
    r = predict_one(results, a, &g, "eye_color");
    check_str(r->prediction, "Likely blue", "HERC2 AA predicts blue");

    load_fixture(&g, a, "rs12913832", "AA", "rs1800407", "CT", NULL);
    r = predict_one(results, a, &g, "eye_color");
    check_str(r->prediction, "Blue or green", "OCA2 variant shifts blue");
    check(r->nsnps == 2, "both SNPs listed");

    load_fixture(&g, a, "rs12913832", "AA", "rs1800407", "GG", NULL);
    r = predict_one(results, a, &g, "eye_color");
    check_str(r->prediction, "Likely blue", "OCA2 without the A allele");

    /* Heterozygous HERC2 is intermediate. */
    load_fixture(&g, a, "rs12913832", "AG", NULL);
    r = predict_one(results, a, &g, "eye_color");
    check_str(r->prediction, "Green, hazel, or light brown", "heterozygous HERC2");

    load_fixture(&g, a, "rs12913832", "AG", "rs1800407", "TT", NULL);
    r = predict_one(results, a, &g, "eye_color");
    check_str(r->prediction, "Green or hazel", "heterozygous plus OCA2 variant");

    /* Reversed spelling must not change the allele count. */
    load_fixture(&g, a, "rs12913832", "GA", NULL);
    r = predict_one(results, a, &g, "eye_color");
    check_str(r->prediction, "Green, hazel, or light brown", "GA equals AG");

    /* The OCA2 SNP alone is not enough. */
    load_fixture(&g, a, "rs1800407", "AA", NULL);
    r = predict_one(results, a, &g, "eye_color");
    check_str(r->prediction, "Unknown", "OCA2 alone cannot predict");
    check(r->nsnps == 1, "OCA2 still listed among SNPs used");

    gh_arena_free(a);
}

/* ------------------------------------------------------------------ */
/* MC1R traits                                                         */
/* ------------------------------------------------------------------ */

static void test_mc1r_traits(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_trait_result results[GH_TRAIT_COUNT];
    const gh_trait_result *hair, *freckles;

    /* Two loss-of-function alleles on one SNP. */
    load_fixture(&g, a, "rs1805007", "TT", NULL);
    gh_predict_traits(results, a, &g);
    hair = gh_trait_find(results, "hair_color");
    freckles = gh_trait_find(results, "freckling");
    check_has(hair->prediction, "red", "two MC1R variants predict red hair");
    check_str(hair->confidence, "high", "red hair call is high confidence");
    check_has(freckles->prediction, "High freckling", "and high freckling");

    /* One allele on each SNP still totals two. */
    load_fixture(&g, a, "rs1805007", "CT", "rs1805008", "CT", NULL);
    gh_predict_traits(results, a, &g);
    hair = gh_trait_find(results, "hair_color");
    check_has(hair->prediction, "red", "variants are summed across both SNPs");
    check(hair->nsnps == 2, "both MC1R SNPs listed");

    /* Exactly one variant allele. */
    load_fixture(&g, a, "rs1805007", "CT", NULL);
    gh_predict_traits(results, a, &g);
    hair = gh_trait_find(results, "hair_color");
    freckles = gh_trait_find(results, "freckling");
    check_has(hair->prediction, "Possible red tint", "one variant is a maybe");
    check_has(freckles->prediction, "Moderate freckling", "moderate freckling");

    /* No variants. */
    load_fixture(&g, a, "rs1805007", "CC", "rs1805008", "CC", NULL);
    gh_predict_traits(results, a, &g);
    hair = gh_trait_find(results, "hair_color");
    freckles = gh_trait_find(results, "freckling");
    check_has(hair->prediction, "Non-red", "no variants means non-red");
    check_str(freckles->prediction, "Typical sun sensitivity", "typical skin");

    /* Neither SNP present. */
    load_fixture(&g, a, "rs00000000", "AA", NULL);
    gh_predict_traits(results, a, &g);
    check_str(gh_trait_find(results, "hair_color")->prediction, "Unknown",
              "no MC1R data");
    check_str(gh_trait_find(results, "freckling")->prediction, "Unknown",
              "no MC1R data for freckling either");

    gh_arena_free(a);
}

/* ------------------------------------------------------------------ */
/* Bitter taste                                                        */
/* ------------------------------------------------------------------ */

static void test_bitter_taste(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_trait_result results[GH_TRAIT_COUNT];
    const gh_trait_result *r;

    /* All taster alleles: ratio 1.0. */
    load_fixture(&g, a, "rs713598", "GG", "rs1726866", "GG",
                 "rs10246939", "CC", NULL);
    r = predict_one(results, a, &g, "bitter_taste");
    check_has(r->prediction, "supertaster", "all taster alleles");
    check(r->nsnps == 3, "three SNPs listed");

    /* Two of three SNPs fully taster: ratio 4/6 = 0.67, medium. */
    load_fixture(&g, a, "rs713598", "GG", "rs1726866", "GG",
                 "rs10246939", "TT", NULL);
    r = predict_one(results, a, &g, "bitter_taste");
    check_str(r->prediction, "Medium bitter taster", "intermediate ratio");

    /* No taster alleles: ratio 0. */
    load_fixture(&g, a, "rs713598", "AA", "rs1726866", "TT",
                 "rs10246939", "TT", NULL);
    r = predict_one(results, a, &g, "bitter_taste");
    check_has(r->prediction, "Non-taster", "no taster alleles");

    /* Only SNPs that are present count toward the denominator, so a single
     * fully-taster SNP still reads as a supertaster. */
    load_fixture(&g, a, "rs713598", "GG", NULL);
    r = predict_one(results, a, &g, "bitter_taste");
    check_has(r->prediction, "supertaster", "absent SNPs are not counted");
    check(r->nsnps == 1, "only the present SNP is listed");

    load_fixture(&g, a, "rs00000000", "AA", NULL);
    r = predict_one(results, a, &g, "bitter_taste");
    check_str(r->prediction, "Unknown", "no TAS2R38 data");

    gh_arena_free(a);
}

/* ------------------------------------------------------------------ */
/* Table-driven traits                                                 */
/* ------------------------------------------------------------------ */

static void test_simple_traits(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_trait_result results[GH_TRAIT_COUNT];
    const gh_trait_result *r;

    struct {
        const char *rsid;
        const char *genotype;
        const char *key;
        const char *expect;
    } cases[] = {
        /* Earwax: T is the dry allele and C is dominant. */
        {"rs17822931", "TT", "earwax_type", "Dry earwax"},
        {"rs17822931", "CT", "earwax_type", "Wet earwax"},
        {"rs17822931", "CC", "earwax_type", "Wet earwax"},
        /* Lactase persistence. */
        {"rs4988235", "AA", "lactose_tolerance", "Lactose tolerant"},
        {"rs4988235", "AG", "lactose_tolerance", "Likely lactose tolerant"},
        {"rs4988235", "GG", "lactose_tolerance", "Likely lactose intolerant"},
        /* ACTN3: TT is the endurance (XX) genotype. */
        {"rs1815739", "CC", "muscle_fiber_type", "Power/sprint oriented (RR)"},
        {"rs1815739", "CT", "muscle_fiber_type", "Mixed power/endurance (RX)"},
        {"rs1815739", "TT", "muscle_fiber_type", "Endurance oriented (XX)"},
        /* FUT2 secretor status. */
        {"rs601338", "AA", "secretor_status", "Non-secretor"},
        {"rs601338", "GA", "secretor_status", "Secretor (carrier)"},
        {"rs601338", "GG", "secretor_status", "Secretor"},
        /* Cilantro. */
        {"rs72921001", "CC", "cilantro_taste",
         "Likely perceives cilantro as soapy"},
        {"rs72921001", "TT", "cilantro_taste",
         "Normal cilantro taste (no soapy perception)"},
        /* Hair curl. */
        {"rs11803731", "AA", "hair_curl", "Likely straight hair"},
        {"rs11803731", "TT", "hair_curl", "Likely wavy or curly hair"},
        /* Unibrow. */
        {"rs12651896", "CC", "unibrow_tendency", "Higher unibrow tendency"},
        {"rs12651896", "TT", "unibrow_tendency", "Lower unibrow tendency"},
        /* Baldness. */
        {"rs2180439", "CC", "baldness_risk", "Higher baldness risk"},
        {"rs2180439", "TT", "baldness_risk", "Lower baldness risk"},
    };

    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        load_fixture(&g, a, cases[i].rsid, cases[i].genotype, NULL);
        r = predict_one(results, a, &g, cases[i].key);
        check_str(r->prediction, cases[i].expect, "simple trait prediction");
        check(r->nsnps == 1, "one SNP listed for a simple trait");
    }

    /* The photic sneeze reflex shares a prediction across CC and CT but not
     * its description, so both must be distinguishable. */
    load_fixture(&g, a, "rs10427255", "CC", NULL);
    r = predict_one(results, a, &g, "photic_sneeze");
    check_str(r->prediction, "Likely sun sneezer", "CC is a sneezer");
    check_has(r->description, "CC genotype", "description names the genotype");

    load_fixture(&g, a, "rs10427255", "CT", NULL);
    r = predict_one(results, a, &g, "photic_sneeze");
    check_str(r->prediction, "Likely sun sneezer", "CT is also a sneezer");
    check_has(r->description, "CT genotype", "description distinguishes CT");

    load_fixture(&g, a, "rs10427255", "TT", NULL);
    r = predict_one(results, a, &g, "photic_sneeze");
    check_str(r->prediction, "Unlikely sun sneezer", "TT is not a sneezer");

    /* Reversed spellings must count the same. */
    load_fixture(&g, a, "rs4988235", "GA", NULL);
    r = predict_one(results, a, &g, "lactose_tolerance");
    check_str(r->prediction, "Likely lactose tolerant", "GA equals AG");

    /* A haploid call counts one copy, not two. */
    load_fixture(&g, a, "rs17822931", "T", NULL);
    r = predict_one(results, a, &g, "earwax_type");
    check_str(r->prediction, "Wet earwax", "single T counts as one copy");

    gh_arena_free(a);
}

int main(void)
{
    test_all_traits_reported();
    test_eye_color();
    test_mc1r_traits();
    test_bitter_taste();
    test_simple_traits();

    printf("%s: %d checks, %d failures\n",
           failures ? "FAILED" : "ok", checks, failures);
    return failures ? 1 : 0;
}

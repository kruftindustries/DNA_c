/* Unit tests for the APOE, blood type, mitochondrial and star-allele
 * scorers. */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gh_genome.h"
#include "gh_mem.h"
#include "gh_scorers.h"

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
/* Genome fixtures                                                     */
/* ------------------------------------------------------------------ */

/* Build a genome from alternating rsid/genotype arguments, NULL-terminated.
 * Goes through the real loader so the fixtures exercise the same path as
 * production input. */
static void load_fixture(gh_genome *g, gh_arena *arena, ...)
{
    static int counter;
    char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_scorer_fixture_%d.txt", counter++);

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

/* ------------------------------------------------------------------ */
/* APOE                                                                */
/* ------------------------------------------------------------------ */

static void test_apoe(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_apoe_result r;

    /* Each epsilon combination, from the rs429358/rs7412 encoding:
     *   T+T = e2, T+C = e3, C+C = e4. */
    struct {
        const char *g358;
        const char *g7412;
        const char *expect;
        const char *risk;
    } cases[] = {
        {"TT", "TT", "e2/e2", "reduced"},
        {"TT", "TC", "e2/e3", "reduced"},
        {"TT", "CC", "e3/e3", "average"},
        {"TC", "CC", "e3/e4", "elevated"},
        {"CC", "CC", "e4/e4", "high"},
        {"TC", "TC", "e2/e4", "moderate"},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        load_fixture(&g, a, "rs429358", cases[i].g358,
                     "rs7412", cases[i].g7412, NULL);
        gh_call_apoe(&r, a, &g);
        check_str(r.apoe_type, cases[i].expect, "APOE haplotype");
        check_str(r.risk_level, cases[i].risk, "APOE risk level");
        check_str(r.confidence, "high", "APOE confidence when both SNPs present");
        check(r.has_or, "APOE odds ratio present");
    }

    /* Phasing must not matter: CT and TC are the same unphased call. */
    load_fixture(&g, a, "rs429358", "CT", "rs7412", "CC", NULL);
    gh_call_apoe(&r, a, &g);
    check_str(r.apoe_type, "e3/e4", "reversed rs429358 spelling");

    load_fixture(&g, a, "rs429358", "TC", "rs7412", "CT", NULL);
    gh_call_apoe(&r, a, &g);
    check_str(r.apoe_type, "e2/e4", "both genotypes reversed");

    /* C at rs429358 with T at rs7412 is not a natural haplotype. */
    load_fixture(&g, a, "rs429358", "CC", "rs7412", "TT", NULL);
    gh_call_apoe(&r, a, &g);
    check_str(r.apoe_type, "Unknown", "impossible C/T combination");
    check_has(r.description, "Unexpected genotype", "explains the rejection");

    /* Missing either SNP. */
    load_fixture(&g, a, "rs429358", "TT", NULL);
    gh_call_apoe(&r, a, &g);
    check_str(r.apoe_type, "Unknown", "missing rs7412");
    check_has(r.description, "Insufficient data", "explains the missing SNP");
    check(!r.has_or, "no odds ratio when unknown");

    load_fixture(&g, a, "rs00000000", "AA", NULL);
    gh_call_apoe(&r, a, &g);
    check_str(r.apoe_type, "Unknown", "neither SNP present");
    check_str(r.confidence, "low", "confidence low when unknown");

    /* Haploid calls: the loader accepts single-base genotypes, and this
     * crashed the Python with an IndexError before it was fixed. */
    load_fixture(&g, a, "rs429358", "T", "rs7412", "TT", NULL);
    gh_call_apoe(&r, a, &g);
    check_str(r.apoe_type, "Unknown", "haploid rs429358 does not decode");
    load_fixture(&g, a, "rs429358", "TT", "rs7412", "C", NULL);
    gh_call_apoe(&r, a, &g);
    check_str(r.apoe_type, "Unknown", "haploid rs7412 does not decode");

    gh_arena_free(a);
}

/* ------------------------------------------------------------------ */
/* Blood type                                                          */
/* ------------------------------------------------------------------ */

static void test_blood_type(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_blood_result r;

    /* rs505922 is the O proxy, rs8176746 the B antigen, rs590787 the Rh
     * proxy (C/C = Rh negative: the minor-allele homozygote; the site is A/C
     * and the RHD-carrying reference genome is A/A). */
    struct {
        const char *proxy;
        const char *b;
        const char *rh;
        const char *type;
        const char *abo;
        const char *rhsign;
    } cases[] = {
        {"TT", "CC", "AA", "O+", "O", "+"},
        {"TT", "CC", "CC", "O-", "O", "-"},
        {"CT", "CC", "AA", "A+", "A", "+"},
        {"TC", "CC", "AC", "A+", "A", "+"},
        {"CC", "CC", "AA", "A+", "A", "+"},
        {"CC", "CT", "AA", "AB+", "AB", "+"},
        {"CC", "TT", "AA", "B+", "B", "+"},
        {"TT", "CT", "CA", "B+", "B", "+"},
        {"CT", "TT", "CC", "B-", "B", "-"},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        load_fixture(&g, a, "rs505922", cases[i].proxy,
                     "rs8176746", cases[i].b, "rs590787", cases[i].rh, NULL);
        gh_predict_blood_type(&r, a, &g);
        check_str(r.blood_type, cases[i].type, "blood type");
        check_str(r.abo, cases[i].abo, "ABO group");
        check_str(r.rh, cases[i].rhsign, "Rh factor");
        check_str(r.confidence, "high", "all three markers means high confidence");
    }

    /* Partial data yields a partial answer rather than a guess. */
    load_fixture(&g, a, "rs505922", "TT", NULL);
    gh_predict_blood_type(&r, a, &g);
    check_str(r.blood_type, "O?", "ABO known, Rh unknown");
    check_str(r.rh, "Unknown", "Rh unknown");
    check_str(r.confidence, "moderate", "one marker means moderate confidence");

    load_fixture(&g, a, "rs590787", "CC", NULL);
    gh_predict_blood_type(&r, a, &g);
    check_str(r.blood_type, "?-", "Rh known, ABO unknown");
    check_str(r.abo, "Unknown", "ABO unknown");

    load_fixture(&g, a, "rs00000000", "AA", NULL);
    gh_predict_blood_type(&r, a, &g);
    check_str(r.blood_type, "Unknown", "no markers at all");
    check_str(r.confidence, "low", "no markers means low confidence");

    gh_arena_free(a);
}

/* ------------------------------------------------------------------ */
/* Mitochondrial haplogroup                                            */
/* ------------------------------------------------------------------ */

static void test_mt_haplogroup(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_mt_result r;

    /* No MT markers at all. */
    load_fixture(&g, a, "rs00000000", "AA", NULL);
    gh_estimate_mt_haplogroup(&r, a, &g);
    check_str(r.haplogroup, "Unknown", "no markers means unknown haplogroup");
    check_str(r.confidence, "none", "no markers means no confidence");
    check(r.markers_found == 0, "no markers found");
    check(r.markers_tested == GH_MT_TREE_COUNT, "markers tested is the tree size");
    check_has(r.lineage, "maternal", "lineage is still described");

    /* A single defining marker, haploid as MT calls usually are. */
    load_fixture(&g, a, "rs2854128", "A", NULL);    /* H: m.2706A */
    gh_estimate_mt_haplogroup(&r, a, &g);
    check_str(r.haplogroup, "H", "single marker call");
    check_str(r.lineage, "European maternal", "lineage from the map");
    check_str(r.confidence, "low", "one marker means low confidence");
    check(r.markers_found == 1, "one marker found");

    /* A doubled call means the same as a haploid one. */
    load_fixture(&g, a, "rs2854128", "AA", NULL);
    gh_estimate_mt_haplogroup(&r, a, &g);
    check_str(r.haplogroup, "H", "doubled call matches too");

    /* Present but not the defining allele: counted as found, not matched. */
    load_fixture(&g, a, "rs2854128", "GG", NULL);
    gh_estimate_mt_haplogroup(&r, a, &g);
    check_str(r.haplogroup, "Unknown", "wrong allele does not match");
    check(r.markers_found == 1, "marker still counted as present");
    check(r.nmatches == 0, "no matches recorded");

    /* Later tree entries are more specific and win. rs2854128 (H) comes
     * before rs3928306 (H1) in the tree. */
    load_fixture(&g, a, "rs2854128", "A", "rs3928306", "A", NULL);
    gh_estimate_mt_haplogroup(&r, a, &g);
    check_str(r.haplogroup, "H1", "most specific match wins");
    check(r.nmatches == 2, "both matches recorded");

    /* Confidence bands: 5 markers present is "moderate". */
    load_fixture(&g, a,
                 "rs2853499", "G", "rs3088309", "C", "rs28358571", "T",
                 "rs2854128", "A", "rs28358587", "T", NULL);
    gh_estimate_mt_haplogroup(&r, a, &g);
    check(r.markers_found == 5, "five markers found");
    check_str(r.confidence, "moderate", "five markers means moderate");

    /* A sample that is rCRS at every marker is H and nothing else: the old
     * tree labelled reference bases as L1, N, C... and called everyone C. */
    load_fixture(&g, a,
                 "rs28357968", "G", "rs2854128", "A", "rs2015062", "C",
                 "rs3928306", "G", "rs41456348", "T", "rs3088309", "C",
                 "rs28359172", "A", "rs2853499", "G", "rs2853825", "G",
                 "rs28358280", "A", "rs28358571", "T", "rs28358587", "T", NULL);
    gh_estimate_mt_haplogroup(&r, a, &g);
    check_str(r.haplogroup, "H", "reference bases call H only");
    check(r.nmatches == 2, "both H markers match, nothing else");
    check_str(r.confidence, "high", "all markers present means high");

    gh_arena_free(a);
}

/* ------------------------------------------------------------------ */
/* Star alleles                                                        */
/* ------------------------------------------------------------------ */

/* Result for one gene by name. */
static const gh_star_result *find_gene(const gh_star_result *results,
                                       const char *gene)
{
    for (size_t i = 0; i < GH_STAR_GENE_COUNT; i++)
        if (strcmp(results[i].gene, gene) == 0)
            return &results[i];
    return NULL;
}

static void test_star_alleles(void)
{
    gh_arena *a = gh_arena_new(1 << 18);
    gh_genome g;
    gh_star_result *results =
        gh_calloc(a, GH_STAR_GENE_COUNT, sizeof(*results));

    /* No defining SNPs at all. */
    load_fixture(&g, a, "rs00000000", "AA", NULL);
    gh_call_star_alleles(results, a, &g);
    const gh_star_result *r = find_gene(results, "CYP2C19");
    check(r != NULL, "CYP2C19 is in the results");
    check_str(r->diplotype, "Unknown", "no SNPs means unknown diplotype");
    check_str(r->phenotype, "Unknown", "no SNPs means unknown phenotype");
    check_str(r->confidence, "low", "no SNPs means low confidence");
    check(r->coverage == 0.0, "coverage is zero");
    check_has(r->clinical_note, "No defining SNPs found", "explains the gap");

    /* Reference call: every CYP2C19 SNP present, none carrying a variant.
     * Genotypes are the GRCh37 plus-strand reference bases: rs28399504 A,
     * rs72558186 T and rs41291556 T (their variants are *4 G, *7 A, *8 C);
     * rs12248560 is *17 (T), C/C is reference there. */
    load_fixture(&g, a,
                 "rs4244285", "GG", "rs4986893", "GG", "rs28399504", "AA",
                 "rs56337013", "CC", "rs72552267", "GG", "rs72558186", "TT",
                 "rs41291556", "TT", "rs12248560", "CC", NULL);
    gh_call_star_alleles(results, a, &g);
    r = find_gene(results, "CYP2C19");
    check_str(r->diplotype, "*1/*1", "reference diplotype");
    check_str(r->phenotype, "normal", "reference phenotype");
    check_str(r->confidence, "high", "full coverage means high confidence");
    check(r->snps_found == 8 && r->snps_total == 8, "all SNPs found");
    check(r->coverage == 1.0, "full coverage");
    check_has(r->clinical_note, "All defining SNPs found", "no caveats");

    /* Homozygous *2 (rs4244285 A/A): two no-function copies -> poor. */
    load_fixture(&g, a,
                 "rs4244285", "AA", "rs4986893", "GG", "rs28399504", "AA",
                 "rs56337013", "CC", "rs72552267", "GG", "rs72558186", "TT",
                 "rs41291556", "TT", "rs12248560", "CC", NULL);
    gh_call_star_alleles(results, a, &g);
    r = find_gene(results, "CYP2C19");
    check_str(r->diplotype, "*2/*2", "homozygous variant");
    check_str(r->phenotype, "poor", "two no-function alleles means poor");

    /* Heterozygous *2: one no-function copy alongside the reference. */
    load_fixture(&g, a,
                 "rs4244285", "GA", "rs4986893", "GG", "rs28399504", "AA",
                 "rs56337013", "CC", "rs72552267", "GG", "rs72558186", "TT",
                 "rs41291556", "TT", "rs12248560", "CC", NULL);
    gh_call_star_alleles(results, a, &g);
    r = find_gene(results, "CYP2C19");
    check_str(r->diplotype, "*1/*2", "heterozygous variant");
    check_str(r->phenotype, "intermediate", "one no-function allele");

    /* Homozygous *17 (increased function) -> ultrarapid. */
    load_fixture(&g, a,
                 "rs4244285", "GG", "rs4986893", "GG", "rs28399504", "AA",
                 "rs56337013", "CC", "rs72552267", "GG", "rs72558186", "TT",
                 "rs41291556", "TT", "rs12248560", "TT", NULL);
    gh_call_star_alleles(results, a, &g);
    r = find_gene(results, "CYP2C19");
    check_str(r->diplotype, "*17/*17", "homozygous increased-function allele");
    check_str(r->phenotype, "ultrarapid", "two increased-function alleles");

    /* Partial coverage lowers confidence and is called out in the note. */
    load_fixture(&g, a, "rs4244285", "GA", "rs4986893", "GG", NULL);
    gh_call_star_alleles(results, a, &g);
    r = find_gene(results, "CYP2C19");
    check(r->snps_found == 2, "two of eight SNPs found");
    check(r->coverage > 0.24 && r->coverage < 0.26, "coverage is 0.25");
    check_str(r->confidence, "low", "quarter coverage means low confidence");
    check_has(r->clinical_note, "Missing SNPs:", "missing SNPs are listed");
    check_has(r->clinical_note, "Result based on available data only",
              "caveat is spelled out");

    /* Every gene is reported, even with nothing to go on. */
    load_fixture(&g, a, "rs00000000", "AA", NULL);
    gh_call_star_alleles(results, a, &g);
    for (size_t i = 0; i < GH_STAR_GENE_COUNT; i++) {
        check(results[i].gene != NULL, "gene name set");
        check(results[i].diplotype != NULL, "diplotype set");
        check(results[i].phenotype != NULL, "phenotype set");
        check(results[i].clinical_note != NULL, "clinical note set");
    }

    /* Multi-SNP haplotypes. SLCO1B1*15 is 388A>G (rs2306283 G) plus
     * 521T>C (rs4149056 C); 388G alone is *17 and 521C alone is *5. The
     * copies a haplotype uses are consumed, so one 521C on a 388G/G
     * background is *15 on one chromosome and *17 on the other. */
    load_fixture(&g, a, "rs4149056", "TT", "rs2306283", "GG", NULL);
    gh_call_star_alleles(results, a, &g);
    r = find_gene(results, "SLCO1B1");
    check_str(r->diplotype, "*17/*17", "common 388G alone is not *15");

    load_fixture(&g, a, "rs4149056", "TC", "rs2306283", "GG", NULL);
    gh_call_star_alleles(results, a, &g);
    r = find_gene(results, "SLCO1B1");
    check_str(r->diplotype, "*15/*17", "*15 consumes one 388G copy");

    load_fixture(&g, a, "rs4149056", "CC", "rs2306283", "GG", NULL);
    gh_call_star_alleles(results, a, &g);
    r = find_gene(results, "SLCO1B1");
    check_str(r->diplotype, "*15/*15", "homozygous *15");

    load_fixture(&g, a, "rs4149056", "TC", "rs2306283", "AA", NULL);
    gh_call_star_alleles(results, a, &g);
    r = find_gene(results, "SLCO1B1");
    check_str(r->diplotype, "*1/*5", "521C without 388G is *5");

    gh_arena_free(a);
}

static void test_phenotype_map(void)
{
    /* The pair is looked up in sorted order, so argument order must not
     * matter. */
    check_str(gh_star_phenotype("normal", "normal"), "normal", "normal pair");
    check_str(gh_star_phenotype("increased", "increased"), "ultrarapid",
              "two increased");
    check_str(gh_star_phenotype("increased", "normal"), "rapid",
              "increased plus normal");
    check_str(gh_star_phenotype("normal", "increased"), "rapid",
              "argument order does not matter");
    check_str(gh_star_phenotype("no_function", "no_function"), "poor",
              "two no-function");
    check_str(gh_star_phenotype("decreased", "no_function"), "poor",
              "decreased plus no-function");
    check_str(gh_star_phenotype("no_function", "normal"), "intermediate",
              "no-function plus normal");
    check_str(gh_star_phenotype("decreased", "increased"), "normal",
              "one up one down is normal");
    check(gh_star_phenotype("normal", "made_up") == NULL,
          "unknown function has no phenotype");

    /* Every combination of the four function names must be covered, which is
     * what the Python asserts at import time. */
    static const char *const functions[] = {
        "normal", "increased", "decreased", "no_function"
    };
    for (size_t i = 0; i < 4; i++)
        for (size_t j = 0; j < 4; j++)
            check(gh_star_phenotype(functions[i], functions[j]) != NULL,
                  "every function pair maps to a phenotype");
}

int main(void)
{
    test_apoe();
    test_blood_type();
    test_mt_haplogroup();
    test_star_alleles();
    test_phenotype_map();

    printf("%s: %d checks, %d failures\n",
           failures ? "FAILED" : "ok", checks, failures);
    return failures ? 1 : 0;
}

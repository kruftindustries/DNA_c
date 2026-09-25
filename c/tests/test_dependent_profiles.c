/* Unit tests for the profiles that consume upstream results. */

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "gh_dependent_profiles.h"
#include "gh_genome.h"
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

static void load_fixture(gh_genome *g, gh_arena *arena, ...)
{
    static int counter;
    char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_dep_fixture_%d.txt", counter++);

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

/* Build a genome where `carrying` chronotype markers hold the evening
 * allele, homozygous. */
static void load_chronotype(gh_genome *g, gh_arena *arena, size_t carrying)
{
    static int counter;
    char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_chrono_%d.txt", counter++);

    FILE *f = fopen(path, "wb");
    if (!f) {
        failures++;
        return;
    }
    fputs("# fixture\n", f);
    for (size_t i = 0; i < GH_CHRONOTYPE_SNP_COUNT; i++) {
        char allele = GH_CHRONOTYPE_SNPS[i].evening_allele;
        char other = allele == 'A' ? 'C' : 'A';
        char use = i < carrying ? allele : other;
        fprintf(f, "%s\t1\t%zu\t%c%c\n", GH_CHRONOTYPE_SNPS[i].rsid,
                1000 + i, use, use);
    }
    fclose(f);
    if (!gh_genome_load(g, arena, path, NULL))
        failures++;
}

static void test_sleep_bands(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_sleep_result r;
    gh_analysis empty = {0};

    /* Nothing genotyped: neutral score and a low-confidence call. */
    load_fixture(&g, a, "rs00000000", "AA", NULL);
    gh_profile_sleep(&r, a, &g, &empty);
    check(r.markers_found == 0, "no chronotype markers");
    check(r.chronotype_score == 50.0, "neutral score without markers");
    check_str(r.confidence, "low", "low confidence");
    check_str(r.chronotype, "Intermediate (Neither)", "intermediate band");
    check(r.nrecommendations >= 3, "the three standing recommendations");

    /* All markers evening: the top band. */
    load_chronotype(&g, a, GH_CHRONOTYPE_SNP_COUNT);
    gh_profile_sleep(&r, a, &g, &empty);
    check(r.chronotype_score == 100.0, "all evening alleles scores 100");
    check_str(r.chronotype, "Definite Evening (Night Owl)", "evening band");
    check_str(r.confidence, "high", "all markers means high confidence");

    /* None: the bottom band. */
    load_chronotype(&g, a, 0);
    gh_profile_sleep(&r, a, &g, &empty);
    check(r.chronotype_score == 0.0, "no evening alleles scores 0");
    check_str(r.chronotype, "Definite Morning (Early Bird)", "morning band");

    gh_arena_free(a);
}

static void test_sleep_caffeine_override(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_sleep_result r;

    load_chronotype(&g, a, 0);

    /* Without a caffeine signal the cutoff comes from the chronotype. */
    gh_analysis none = {0};
    gh_profile_sleep(&r, a, &g, &none);
    check(!r.caffeine_sensitive, "not caffeine sensitive by default");
    const char *baseline = r.caffeine_cutoff;

    /* Slow CYP1A2 overrides it. */
    gh_finding slow[] = {
        {"rs762551", "CYP1A2", "Drug", "CC", "slow", "d", "", NULL, 3},
    };
    gh_analysis sensitive = {0};
    sensitive.findings = slow;
    sensitive.nfindings = 1;
    gh_profile_sleep(&r, a, &g, &sensitive);
    check(r.caffeine_sensitive, "slow CYP1A2 flags caffeine sensitivity");
    check(strcmp(r.caffeine_cutoff, baseline) != 0, "cutoff is overridden");
    check_has(r.caffeine_cutoff, "10:00 AM", "overridden to the earlier cutoff");

    /* ADORA2A anxiety proneness does the same. */
    gh_finding anxious[] = {
        {"rs5751876", "ADORA2A", "Neuro", "TT", "anxiety_prone", "d", "", NULL, 3},
    };
    sensitive.findings = anxious;
    gh_profile_sleep(&r, a, &g, &sensitive);
    check(r.caffeine_sensitive, "ADORA2A also flags sensitivity");

    gh_arena_free(a);
}

static void test_nutrigenomics_bands(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_nutrigenomics_result r;

    /* No findings: every nutrient reads as normal. */
    gh_analysis empty = {0};
    gh_profile_nutrigenomics(&r, a, &empty);
    check(r.nneeds == GH_NUTRIENT_PROFILE_COUNT, "every nutrient reported");
    for (size_t i = 0; i < r.nneeds; i++) {
        check_str(r.needs[i].need_level, "normal", "normal without signals");
        check(r.needs[i].severity == 0.0, "zero severity");
    }
    check_has(r.summary, "No major nutrient", "summary says nothing found");

    /* Drive one nutrient's first gene to its heaviest risk status. */
    const gh_nutrient_profile *profile = &GH_NUTRIENT_PROFILES[0];
    const gh_nutrient_gene *gene = &profile->genes[0];
    double heaviest = 0.0;
    const char *heaviest_status = NULL;
    for (size_t i = 0; i < gene->nstatuses; i++)
        if (gene->statuses[i].weight > heaviest) {
            heaviest = gene->statuses[i].weight;
            heaviest_status = gene->statuses[i].status;
        }

    gh_finding findings[] = {
        {"rs1", gene->gene, "Cat", "AA", heaviest_status, "d", "", NULL, 3},
    };
    gh_analysis one = {0};
    one.findings = findings;
    one.nfindings = 1;
    gh_profile_nutrigenomics(&r, a, &one);

    const gh_nutrient_need *need = NULL;
    for (size_t i = 0; i < r.nneeds; i++)
        if (strcmp(r.needs[i].profile->id, profile->id) == 0)
            need = &r.needs[i];
    check(need != NULL, "the nutrient is reported");
    if (need) {
        check(need->severity > 0.0, "severity accumulated");
        check(need->nimpacts == 1, "one gene impact recorded");
        check_str(need->impacts[0].gene, gene->gene, "impact names the gene");
        check(strcmp(need->need_level, "normal") != 0,
              "a risk status moves it off normal");
    }

    /* Results are ordered by descending absolute severity. */
    for (size_t i = 1; i < r.nneeds; i++) {
        double prev = r.needs[i - 1].severity;
        double cur = r.needs[i].severity;
        check((prev < 0 ? -prev : prev) >= (cur < 0 ? -cur : cur),
              "sorted by descending absolute severity");
    }

    gh_arena_free(a);
}

static void test_mental_health(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_mental_result r;
    gh_analysis empty = {0};

    load_fixture(&g, a, "rs00000000", "AA", NULL);
    gh_profile_mental_health(&r, a, &g, &empty, NULL, 0);
    check(r.ndomains == GH_MENTAL_DOMAIN_COUNT, "every domain reported");
    check(r.nmarkers == 0, "no markers genotyped");
    check(r.nrecommendations >= 3, "the three standing recommendations");
    check_has(r.summary, "No strongly elevated", "summary when nothing fires");

    /* A homozygous risk marker raises its domain and lands in risk factors. */
    const gh_mental_snp *snp = &GH_MENTAL_SNPS[0];
    char genotype[3] = {snp->risk_allele, snp->risk_allele, '\0'};
    load_fixture(&g, a, snp->rsid, genotype, NULL);
    gh_profile_mental_health(&r, a, &g, &empty, NULL, 0);
    check(r.nmarkers == 1, "marker detected");
    check(r.markers[0].risk_copies == 2, "two risk copies");
    check(r.nrisk_factors >= 1, "a risk factor recorded");
    check_has(r.risk_factors[0], "Homozygous risk", "homozygous wording");

    /* A reduced CYP2C19 metabolizer produces an SSRI dosing note. */
    gh_star_result stars[] = {
        {"CYP2C19", "*2/*2", "poor", "high", "note", 8, 8, 1.0},
    };
    gh_profile_mental_health(&r, a, &g, &empty, stars, 1);
    check(r.ntreatment_notes >= 1, "a treatment note is produced");
    check_has(r.treatment_notes[0], "CYP2C19 reduced", "reduced-metabolizer note");

    stars[0].phenotype = "ultrarapid";
    gh_profile_mental_health(&r, a, &g, &empty, stars, 1);
    check_has(r.treatment_notes[0], "ultrarapid", "ultrarapid note");

    stars[0].phenotype = "normal";
    gh_profile_mental_health(&r, a, &g, &empty, stars, 1);
    for (size_t i = 0; i < r.ntreatment_notes; i++)
        check(strstr(r.treatment_notes[i], "CYP2C19") == NULL,
              "a normal metabolizer produces no CYP2C19 note");

    /* Risk factors are capped at eight. */
    check(r.nrisk_factors <= 8, "risk factors capped");

    gh_arena_free(a);
}

static void test_longevity(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    gh_genome g;
    gh_longevity_result r;
    gh_analysis empty = {0};

    load_fixture(&g, a, "rs00000000", "AA", NULL);
    gh_profile_longevity(&r, a, &g, &empty, NULL, NULL, 0);
    check(r.alleles_checked == 0, "no longevity alleles genotyped");
    check(r.ndomains == GH_HEALTHSPAN_DOMAIN_COUNT, "every domain reported");
    check(r.ninterventions >= 3, "the standing interventions are offered");
    for (size_t i = 0; i < r.ndomains; i++) {
        check(r.domains[i].score >= 10 && r.domains[i].score <= 90,
              "domain score stays within its clamp");
        check(r.domains[i].genes_found == 0, "no genes found");
    }

    /* Every longevity allele homozygous: the top of the range. */
    static int counter;
    char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_longev_%d.txt", counter++);
    FILE *f = fopen(path, "wb");
    if (f) {
        fputs("# fixture\n", f);
        for (size_t i = 0; i < GH_LONGEVITY_SNP_COUNT; i++) {
            char allele = GH_LONGEVITY_SNPS[i].longevity_allele;
            fprintf(f, "%s\t1\t%zu\t%c%c\n", GH_LONGEVITY_SNPS[i].rsid,
                    1000 + i, allele, allele);
        }
        fclose(f);
        gh_genome_load(&g, a, path, NULL);
    }
    gh_profile_longevity(&r, a, &g, &empty, NULL, NULL, 0);
    check(r.alleles_checked == GH_LONGEVITY_SNP_COUNT, "all alleles checked");
    check(r.longevity_score == 100.0, "all protective scores 100");
    check_has(r.summary, "favorable", "favourable summary");
    check(r.ntop_protective > 0, "protective entries recorded");

    /* An elevated APOE result becomes a top risk. */
    gh_apoe_result apoe = {
        .apoe_type = "e3/e4", .risk_level = "elevated", .alzheimer_or = 2.8,
        .has_or = true, .description = "One e4 allele raises risk.",
        .confidence = "high", .rs429358 = "TC", .rs7412 = "CC",
    };
    gh_profile_longevity(&r, a, &g, &empty, &apoe, NULL, 0);
    bool apoe_risk = false;
    for (size_t i = 0; i < r.ntop_risks; i++)
        if (strstr(r.top_risks[i], "APOE"))
            apoe_risk = true;
    check(apoe_risk, "elevated APOE appears among the risks");

    /* A reduced APOE result is protective instead. */
    apoe.risk_level = "reduced";
    apoe.apoe_type = "e2/e3";
    gh_profile_longevity(&r, a, &g, &empty, &apoe, NULL, 0);
    bool apoe_protective = false;
    for (size_t i = 0; i < r.ntop_protective; i++)
        if (strstr(r.top_protective[i], "APOE"))
            apoe_protective = true;
    check(apoe_protective, "reduced APOE appears among the protective entries");

    check(r.ntop_risks <= 8 && r.ntop_protective <= 8, "both lists capped");

    gh_arena_free(a);
}

static void test_tables(void)
{
    for (size_t i = 0; i < GH_CHRONOTYPE_SNP_COUNT; i++) {
        check(strncmp(GH_CHRONOTYPE_SNPS[i].rsid, "rs", 2) == 0, "chronotype rsID");
        check(strchr("ACGT", GH_CHRONOTYPE_SNPS[i].evening_allele) != NULL,
              "chronotype allele is a base");
        check(GH_CHRONOTYPE_SNPS[i].weight > 0, "chronotype weight positive");
    }
    for (size_t i = 0; i < GH_LONGEVITY_SNP_COUNT; i++)
        check(strchr("ACGT", GH_LONGEVITY_SNPS[i].longevity_allele) != NULL,
              "longevity allele is a base");
    for (size_t i = 0; i < GH_MENTAL_SNP_COUNT; i++)
        check(strchr("ACGT", GH_MENTAL_SNPS[i].risk_allele) != NULL,
              "mental risk allele is a base");
    for (size_t i = 0; i < GH_NUTRIENT_PROFILE_COUNT; i++) {
        check(GH_NUTRIENT_PROFILES[i].ngenes > 0, "nutrient has genes");
        check(GH_NUTRIENT_PROFILES[i].supplement_form != NULL, "supplement form");
    }
    check(GH_MENTAL_DOMAIN_COUNT <= GH_DP_MAX_ITEMS, "mental domains fit");
    check(GH_HEALTHSPAN_DOMAIN_COUNT <= GH_DP_MAX_ITEMS, "healthspan domains fit");
    check(GH_LONGEVITY_SNP_COUNT <= GH_DP_MAX_ITEMS, "longevity alleles fit");
    check(GH_CHRONOTYPE_SNP_COUNT <= GH_DP_MAX_ITEMS, "chronotype markers fit");
}

int main(void)
{
    test_sleep_bands();
    test_sleep_caffeine_override();
    test_nutrigenomics_bands();
    test_mental_health();
    test_longevity();
    test_tables();

    printf("%s: %d checks, %d failures\n",
           failures ? "FAILED" : "ok", checks, failures);
    return failures ? 1 : 0;
}

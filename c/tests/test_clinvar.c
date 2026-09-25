/* Unit tests for ClinVar analysis, ACMG flagging and carrier screening.
 *
 * These use a synthetic ClinVar extract rather than the real one, so each
 * case can pin down a specific branch.
 *
 * Note that the real extract's `inheritance_modes` column carries variant
 * origin (germline/somatic) rather than a mode of inheritance. Carrier
 * screening therefore resolves inheritance from the curated gene table
 * first and only falls back to that column, which is what makes it work on
 * real data at all.
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gh_clinvar.h"
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

/* ------------------------------------------------------------------ */
/* Fixtures                                                            */
/* ------------------------------------------------------------------ */

static const char *CLINVAR_HEADER =
    "chrom\tpos\tref\talt\tclinical_significance\treview_status\tgold_stars\t"
    "all_traits\tsymbol\tinheritance_modes\thgvs_p\thgvs_c\t"
    "molecular_consequence\txrefs\n";

/* Write a ClinVar extract from pre-formatted rows. */
static const char *write_clinvar(const char *name, const char *rows)
{
    static char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_clinvar_%s.tsv", name);
    FILE *f = fopen(path, "wb");
    if (!f) {
        failures++;
        return path;
    }
    fputs(CLINVAR_HEADER, f);
    fputs(rows, f);
    fclose(f);
    return path;
}

static const char *write_genome(const char *name, const char *rows)
{
    static char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_cvgenome_%s.txt", name);
    FILE *f = fopen(path, "wb");
    if (!f) {
        failures++;
        return path;
    }
    fputs("# fixture\nrsid\tchromosome\tposition\tgenotype\n", f);
    fputs(rows, f);
    fclose(f);
    return path;
}

/* ------------------------------------------------------------------ */
/* Zygosity                                                            */
/* ------------------------------------------------------------------ */

static void test_classify_zygosity(void)
{
    const char *status, *desc;

    gh_classify_zygosity(true, false, "autosomal recessive", &status, &desc);
    check_str(status, "AFFECTED", "homozygous is affected");
    check_str(desc, "Homozygous for variant", "homozygous description");

    gh_classify_zygosity(false, true, "autosomal recessive", &status, &desc);
    check_str(status, "CARRIER", "heterozygous recessive is a carrier");
    check_str(desc, "Heterozygous carrier (recessive)", "carrier description");

    gh_classify_zygosity(false, true, "autosomal dominant", &status, &desc);
    check_str(status, "AFFECTED", "heterozygous dominant is affected");

    gh_classify_zygosity(false, true, "germline", &status, &desc);
    check_str(status, "HETEROZYGOUS", "unknown inheritance stays heterozygous");
    check_str(desc, "Heterozygous (inheritance unclear)", "unclear description");

    gh_classify_zygosity(false, true, "", &status, &desc);
    check_str(status, "HETEROZYGOUS", "empty inheritance stays heterozygous");

    gh_classify_zygosity(false, true, NULL, &status, &desc);
    check_str(status, "HETEROZYGOUS", "NULL inheritance is safe");

    gh_classify_zygosity(false, false, "autosomal recessive", &status, &desc);
    check_str(status, "UNKNOWN", "neither zygosity is unknown");

    /* Recessive is tested before dominant, so a string with both reads as
     * a carrier. */
    gh_classify_zygosity(false, true, "autosomal dominant / autosomal recessive",
                         &status, &desc);
    check_str(status, "CARRIER", "recessive wins when both are present");
}

/* ------------------------------------------------------------------ */
/* ClinVar scanning                                                    */
/* ------------------------------------------------------------------ */

static void test_clinvar_categories(void)
{
    gh_arena *a = gh_arena_new(1 << 18);

    const char *clinvar = write_clinvar("categories",
        "1\t100\tA\tG\tPathogenic\tcriteria provided\t2\tCondition A\tBRCA1\tautosomal dominant\t\t\t\t\n"
        "1\t200\tC\tT\tLikely_pathogenic\tcriteria provided\t1\tCondition B\tCFTR\tautosomal recessive\t\t\t\t\n"
        "1\t300\tG\tA\trisk factor\tno assertion\t0\tCondition C\tAPOE\t\t\t\t\t\n"
        "1\t400\tT\tC\tdrug response\tcriteria provided\t3\tCondition D\tVKORC1\t\t\t\t\t\n"
        "1\t500\tA\tC\tprotective\tcriteria provided\t1\tCondition E\tCCR5\t\t\t\t\t\n"
        /* Likely pathogenic must not fall into the pathogenic bucket. */
        "1\t600\tA\tG\tConflicting_interpretations_of_pathogenicity\treviewed\t1\tF\tTP53\t\t\t\t\t\n");

    const char *genome = write_genome("categories",
        "rs1\t1\t100\tGG\n"
        "rs2\t1\t200\tCT\n"
        "rs3\t1\t300\tAA\n"
        "rs4\t1\t400\tCC\n"
        "rs5\t1\t500\tAC\n"
        "rs6\t1\t600\tAG\n");

    gh_genome g;
    check(gh_genome_load(&g, a, genome, NULL), "genome loads");

    gh_clinvar_result cv;
    check(gh_clinvar_analyze(&cv, a, &g, clinvar), "ClinVar loads");
    check(cv.loaded, "result marked loaded");
    check(cv.total_clinvar == 6, "all rows scanned");
    check(cv.matched == 6, "all positions matched");

    check(cv.counts[GH_CV_PATHOGENIC] == 1, "one pathogenic");
    check(cv.counts[GH_CV_LIKELY_PATHOGENIC] == 1, "one likely pathogenic");
    check(cv.counts[GH_CV_RISK_FACTOR] == 1, "one risk factor");
    check(cv.counts[GH_CV_DRUG_RESPONSE] == 1, "one drug response");
    check(cv.counts[GH_CV_PROTECTIVE] == 1, "one protective");

    const gh_cv_finding *path = &cv.by_category[GH_CV_PATHOGENIC][0];
    check_str(path->gene, "BRCA1", "pathogenic gene");
    check_str(path->rsid, "rs1", "rsID carried from the genome");
    check_str(path->user_genotype, "GG", "genotype carried from the genome");
    check(path->is_homozygous, "homozygous alt detected");
    check(!path->is_heterozygous, "not also heterozygous");
    check(path->gold_stars == 2, "gold stars parsed");
    check_str(path->zygosity_status, "AFFECTED", "homozygous is affected");

    const gh_cv_finding *likely = &cv.by_category[GH_CV_LIKELY_PATHOGENIC][0];
    check(likely->is_heterozygous, "heterozygous detected");
    check_str(likely->zygosity_status, "CARRIER", "recessive het is a carrier");

    /* Conflicting interpretations are excluded from the pathogenic bucket
     * and match no other category, so they are dropped entirely. */
    size_t total = 0;
    for (size_t i = 0; i < GH_CV_CATEGORY_COUNT; i++)
        total += cv.counts[i];
    check(total == 5, "conflicting interpretation is not categorised");

    gh_arena_free(a);
}

static void test_clinvar_filters(void)
{
    gh_arena *a = gh_arena_new(1 << 18);

    const char *clinvar = write_clinvar("filters",
        /* Indels are skipped: array genotypes cannot represent them. */
        "1\t100\tAT\tG\tPathogenic\treviewed\t2\tX\tGENE1\t\t\t\t\t\n"
        "1\t110\tA\tGC\tPathogenic\treviewed\t2\tX\tGENE2\t\t\t\t\t\n"
        /* Reference-only genotype carries no variant. */
        "1\t200\tA\tG\tPathogenic\treviewed\t2\tX\tGENE3\t\t\t\t\t\n"
        /* Variant present but genome has a different base entirely. */
        "1\t300\tA\tG\tPathogenic\treviewed\t2\tX\tGENE4\t\t\t\t\t\n"
        /* Missing chrom or pos. */
        "\t400\tA\tG\tPathogenic\treviewed\t2\tX\tGENE5\t\t\t\t\t\n"
        "1\t\tA\tG\tPathogenic\treviewed\t2\tX\tGENE6\t\t\t\t\t\n"
        /* Position the genome does not cover. */
        "1\t999\tA\tG\tPathogenic\treviewed\t2\tX\tGENE7\t\t\t\t\t\n"
        /* Unparseable gold stars default to zero. */
        "1\t500\tA\tG\tPathogenic\treviewed\tn/a\tX\tGENE8\t\t\t\t\t\n");

    const char *genome = write_genome("filters",
        "rs1\t1\t100\tGG\n"
        "rs2\t1\t110\tGG\n"
        "rs3\t1\t200\tAA\n"
        "rs4\t1\t300\tCC\n"
        "rs5\t1\t400\tGG\n"
        "rs6\t1\t500\tGG\n");

    gh_genome g;
    gh_genome_load(&g, a, genome, NULL);

    gh_clinvar_result cv;
    check(gh_clinvar_analyze(&cv, a, &g, clinvar), "ClinVar loads");
    check(cv.total_clinvar == 8, "all rows scanned");

    /* Only the last row survives every filter. */
    check(cv.counts[GH_CV_PATHOGENIC] == 1, "only one row survives filtering");
    const gh_cv_finding *f = &cv.by_category[GH_CV_PATHOGENIC][0];
    check_str(f->gene, "GENE8", "the surviving row is the expected one");
    check(f->gold_stars == 0, "unparseable gold stars default to zero");

    gh_arena_free(a);
}

static void test_clinvar_missing_file(void)
{
    gh_arena *a = gh_arena_new(1 << 16);
    const char *genome = write_genome("missing", "rs1\t1\t100\tAA\n");

    gh_genome g;
    gh_genome_load(&g, a, genome, NULL);

    gh_clinvar_result cv;
    check(!gh_clinvar_analyze(&cv, a, &g, "/nonexistent/clinvar.tsv"),
          "missing file reports failure");
    check(!cv.loaded, "result is not marked loaded");

    /* Downstream analysis must cope with an unloaded result. */
    gh_acmg_result acmg;
    gh_flag_acmg(&acmg, a, &cv);
    check(acmg.nfindings == 0, "no ACMG findings without ClinVar");
    check(acmg.genes_screened == GH_ACMG_GENE_COUNT, "gene count still reported");
    check(strstr(acmg.summary, "No ClinVar data") != NULL,
          "ACMG summary explains the gap");

    gh_carrier_result carriers;
    gh_organize_carriers(&carriers, a, &cv);
    check(carriers.ncarriers == 0, "no carriers without ClinVar");

    gh_arena_free(a);
}

/* ------------------------------------------------------------------ */
/* ACMG                                                                */
/* ------------------------------------------------------------------ */

static void test_acmg(void)
{
    gh_arena *a = gh_arena_new(1 << 18);

    const char *clinvar = write_clinvar("acmg",
        "1\t100\tA\tG\tPathogenic\treviewed\t1\tBreast cancer\tBRCA1\t\t\t\t\t\n"
        "1\t200\tA\tG\tPathogenic\treviewed\t4\tLi-Fraumeni\tTP53\t\t\t\t\t\n"
        "1\t300\tA\tG\tLikely_pathogenic\treviewed\t2\tLynch\tMLH1\t\t\t\t\t\n"
        /* Not on the ACMG list. */
        "1\t400\tA\tG\tPathogenic\treviewed\t3\tOther\tNOTAGENE\t\t\t\t\t\n"
        /* Lower-case gene name must still match. */
        "1\t500\tA\tG\tPathogenic\treviewed\t2\tHypercholesterolemia\tldlr\t\t\t\t\t\n"
        /* Risk factors are not screened for ACMG. */
        "1\t600\tA\tG\trisk factor\treviewed\t2\tOther\tBRCA2\t\t\t\t\t\n");

    const char *genome = write_genome("acmg",
        "rs1\t1\t100\tGG\nrs2\t1\t200\tGG\nrs3\t1\t300\tGG\n"
        "rs4\t1\t400\tGG\nrs5\t1\t500\tGG\nrs6\t1\t600\tGG\n");

    gh_genome g;
    gh_genome_load(&g, a, genome, NULL);

    gh_clinvar_result cv;
    gh_clinvar_analyze(&cv, a, &g, clinvar);

    gh_acmg_result acmg;
    gh_flag_acmg(&acmg, a, &cv);

    check(acmg.nfindings == 4, "four ACMG findings");
    check(acmg.genes_with_variants == 4, "four distinct ACMG genes");
    check(acmg.genes_screened == 81, "81 genes screened");

    /* Sorted by descending gold stars, then gene name. */
    check_str(acmg.findings[0].finding->gene, "TP53", "highest stars first");
    check(acmg.findings[0].finding->gold_stars == 4, "four stars");
    check_str(acmg.findings[1].finding->gene, "MLH1",
              "two-star genes sort by name");
    check_str(acmg.findings[2].finding->gene, "ldlr",
              "second two-star gene by name");
    check_str(acmg.findings[3].finding->gene, "BRCA1", "one star last");

    check_str(acmg.findings[0].acmg_category, "pathogenic", "category recorded");
    check_str(acmg.findings[1].acmg_category, "likely_pathogenic",
              "likely category recorded");
    check(strstr(acmg.findings[0].actionability, "Li-Fraumeni") != NULL,
          "curated actionability text used");
    check(strstr(acmg.summary, "4 variant(s)") != NULL, "summary counts variants");
    check(strstr(acmg.summary, "4 ACMG-recommended") != NULL,
          "summary counts genes");

    check(gh_acmg_is_actionable_gene("BRCA1"), "BRCA1 is actionable");
    check(!gh_acmg_is_actionable_gene("NOTAGENE"), "unknown gene is not");

    /* No ACMG hits produces the other summary. */
    const char *none = write_clinvar("acmg_none",
        "1\t100\tA\tG\tPathogenic\treviewed\t1\tOther\tNOTAGENE\t\t\t\t\t\n");
    gh_genome g2;
    gh_genome_load(&g2, a, write_genome("acmg_none", "rs1\t1\t100\tGG\n"), NULL);
    gh_clinvar_result cv2;
    gh_clinvar_analyze(&cv2, a, &g2, none);
    gh_acmg_result acmg2;
    gh_flag_acmg(&acmg2, a, &cv2);
    check(acmg2.nfindings == 0, "no ACMG findings");
    check(strstr(acmg2.summary, "No pathogenic") != NULL,
          "summary states nothing was found");

    gh_arena_free(a);
}

/* ------------------------------------------------------------------ */
/* Carrier screening                                                   */
/* ------------------------------------------------------------------ */

static const gh_carrier *find_carrier(const gh_carrier_result *r,
                                      const char *gene)
{
    for (size_t i = 0; i < r->ncarriers; i++)
        if (strcmp(r->carriers[i].gene, gene) == 0)
            return &r->carriers[i];
    return NULL;
}

static void test_carrier_screening(void)
{
    gh_arena *a = gh_arena_new(1 << 18);

    /* Inheritance strings here are the real patterns, which the production
     * ClinVar extract does not currently supply -- see the file header. */
    const char *clinvar = write_clinvar("carriers",
        "1\t100\tA\tG\tPathogenic\treviewed\t3\tCystic fibrosis;Other\tCFTR\tautosomal recessive\t\t\t\t\n"
        "1\t200\tA\tG\tPathogenic\treviewed\t2\tSickle cell disease\tHBB\tautosomal recessive\t\t\t\t\n"
        /* "X-linked recessive", not bare "X-linked": classify_zygosity only
         * reaches CARRIER via the substrings "recessive" or "dominant", so a
         * bare "X-linked" would be filtered out before the curated
         * gene table could supply the X-linked reproductive note. */
        "1\t300\tA\tG\tLikely_pathogenic\treviewed\t1\tG6PD deficiency\tG6PD\tX-linked recessive\t\t\t\t\n"
        /* Homozygous: affected, not a carrier. */
        "1\t400\tA\tG\tPathogenic\treviewed\t2\tTay-Sachs\tHEXA\tautosomal recessive\t\t\t\t\n"
        /* Dominant: affected, not a carrier. */
        "1\t500\tA\tG\tPathogenic\treviewed\t2\tMarfan\tFBN1\tautosomal dominant\t\t\t\t\n"
        /* Unknown gene falls into the Other system. */
        "1\t600\tA\tG\tPathogenic\treviewed\t1\t\tMYSTERYGENE\tautosomal recessive\t\t\t\t\n");

    const char *genome = write_genome("carriers",
        "rs1\t1\t100\tAG\n"   /* het -> carrier */
        "rs2\t1\t200\tAG\n"   /* het -> carrier */
        "rs3\t1\t300\tAG\n"   /* het X-linked -> carrier */
        "rs4\t1\t400\tGG\n"   /* hom -> affected */
        "rs5\t1\t500\tAG\n"   /* het dominant -> affected */
        "rs6\t1\t600\tAG\n"); /* het, unknown gene */

    gh_genome g;
    gh_genome_load(&g, a, genome, NULL);

    gh_clinvar_result cv;
    gh_clinvar_analyze(&cv, a, &g, clinvar);

    gh_carrier_result carriers;
    gh_organize_carriers(&carriers, a, &cv);

    check(carriers.ncarriers == 4, "four carriers");

    const gh_carrier *cftr = find_carrier(&carriers, "CFTR");
    check(cftr != NULL, "CFTR carrier present");
    check_str(cftr->condition, "Cystic fibrosis",
              "condition is the first trait only");
    check_str(cftr->system, "Metabolic/Pulmonary", "system from the gene table");
    check_str(cftr->inheritance, "autosomal recessive", "curated inheritance");
    check(strstr(cftr->reproductive_note, "25% chance") != NULL,
          "recessive reproductive note");
    check(cftr->couples_relevant, "CFTR is a couples-screened condition");
    check_str(cftr->genotype, "AG", "genotype carried through");
    check(cftr->gold_stars == 3, "gold stars carried through");

    /* Carriers are emitted pathogenic-first, then likely-pathogenic, so
     * look them up by gene rather than by position. */
    const gh_carrier *g6pd = find_carrier(&carriers, "G6PD");
    check(g6pd != NULL, "G6PD carrier present");
    /* The curated table spells this "X-linked" with a capital X, so the
     * note lookup has to be case-insensitive to find it. */
    check_str(g6pd->inheritance, "X-linked", "curated inheritance overrides ClinVar");
    check(strstr(g6pd->reproductive_note, "X-linked") != NULL,
          "X-linked reproductive note is found despite the capital X");
    check(strstr(g6pd->reproductive_note, "sons") != NULL,
          "X-linked note names the inheritance risk");
    check(!g6pd->couples_relevant, "G6PD is not couples-screened");

    /* The unknown gene falls back to Other, and to ClinVar's own field. */
    const gh_carrier *unknown = find_carrier(&carriers, "MYSTERYGENE");
    check(unknown != NULL, "unknown-gene carrier present");
    check_str(unknown->system, "Other", "unknown gene lands in Other");
    check_str(unknown->condition, "Unknown condition",
              "empty traits become Unknown condition");

    /* Systems are grouped in first-seen order. G6PD and HBB share
     * Hematologic, so four carriers span three systems. */
    check(carriers.nsystems == 3, "three distinct systems");
    check_str(carriers.systems[0], "Metabolic/Pulmonary", "first system seen");
    check_str(carriers.systems[1], "Hematologic", "second system seen");
    check_str(carriers.systems[2], "Other", "third system seen");

    check(carriers.ncouples_relevant == 2, "two couples-relevant carriers");

    gh_arena_free(a);
}

/* The production extract supplies variant origin where an inheritance mode
 * is expected. Carrier screening has to survive that by consulting the
 * curated gene table, which is the case this pins down. */
static void test_carrier_screening_with_origin_values(void)
{
    gh_arena *a = gh_arena_new(1 << 18);

    const char *clinvar = write_clinvar("origin",
        "1\t100\tA\tG\tPathogenic\treviewed\t3\tCystic fibrosis\tCFTR\tgermline\t\t\t\t\n"
        "1\t200\tA\tG\tPathogenic\treviewed\t2\tSickle cell\tHBB\tnot applicable\t\t\t\t\n"
        /* Not in the curated table, so the origin value is all there is. */
        "1\t300\tA\tG\tPathogenic\treviewed\t1\tSomething\tUNCURATED\tgermline\t\t\t\t\n");
    const char *genome = write_genome("origin",
        "rs1\t1\t100\tAG\nrs2\t1\t200\tAG\nrs3\t1\t300\tAG\n");

    gh_genome g;
    gh_genome_load(&g, a, genome, NULL);

    gh_clinvar_result cv;
    gh_clinvar_analyze(&cv, a, &g, clinvar);
    check(cv.counts[GH_CV_PATHOGENIC] == 3, "all three variants are found");

    /* The stored zygosity still reflects ClinVar's own field, which cannot
     * express a carrier; carrier screening resolves inheritance separately. */
    check_str(cv.by_category[GH_CV_PATHOGENIC][0].zygosity_status,
              "HETEROZYGOUS", "origin values leave the stored zygosity unclear");

    gh_carrier_result carriers;
    gh_organize_carriers(&carriers, a, &cv);
    check(carriers.ncarriers == 2,
          "curated genes are still found despite the origin values");
    check(find_carrier(&carriers, "CFTR") != NULL, "CFTR carrier found");
    check(find_carrier(&carriers, "HBB") != NULL, "HBB carrier found");
    check(find_carrier(&carriers, "UNCURATED") == NULL,
          "an uncurated gene has nothing usable to go on");

    gh_arena_free(a);
}

/* A gene outside the curated table still uses ClinVar's field when that
 * field happens to carry a real inheritance pattern. */
static void test_carrier_screening_uncurated_fallback(void)
{
    gh_arena *a = gh_arena_new(1 << 18);

    const char *clinvar = write_clinvar("fallback",
        "1\t100\tA\tG\tPathogenic\treviewed\t2\tSomething\tUNCURATED\tautosomal recessive\t\t\t\t\n");
    const char *genome = write_genome("fallback", "rs1\t1\t100\tAG\n");

    gh_genome g;
    gh_genome_load(&g, a, genome, NULL);

    gh_clinvar_result cv;
    gh_clinvar_analyze(&cv, a, &g, clinvar);

    gh_carrier_result carriers;
    gh_organize_carriers(&carriers, a, &cv);
    check(carriers.ncarriers == 1, "ClinVar inheritance is used as a fallback");
    if (carriers.ncarriers) {
        check_str(carriers.carriers[0].inheritance, "autosomal recessive",
                  "fallback inheritance is lower-cased");
        check_str(carriers.carriers[0].system, "Other",
                  "uncurated gene lands in Other");
    }

    gh_arena_free(a);
}

int main(void)
{
    test_classify_zygosity();
    test_clinvar_categories();
    test_clinvar_filters();
    test_clinvar_missing_file();
    test_acmg();
    test_carrier_screening();
    test_carrier_screening_with_origin_values();
    test_carrier_screening_uncurated_fallback();

    printf("%s: %d checks, %d failures\n",
           failures ? "FAILED" : "ok", checks, failures);
    return failures ? 1 : 0;
}

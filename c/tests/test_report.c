/* Unit tests for the HTML report generator and quality metrics. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gh_analysis.h"
#include "gh_genome.h"
#include "gh_mem.h"
#include "gh_quality.h"
#include "gh_report.h"
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

/* Assert that `haystack` contains `needle`. */
static void check_has(const char *haystack, const char *needle, const char *what)
{
    checks++;
    if (!haystack || !strstr(haystack, needle)) {
        failures++;
        fprintf(stderr, "  FAIL: %s\n    missing: %s\n", what, needle);
    }
}

static void check_lacks(const char *haystack, const char *needle, const char *what)
{
    checks++;
    if (haystack && strstr(haystack, needle)) {
        failures++;
        fprintf(stderr, "  FAIL: %s\n    unexpectedly present: %s\n", what, needle);
    }
}

/* ------------------------------------------------------------------ */
/* Fragment helpers                                                    */
/* ------------------------------------------------------------------ */

static void test_mag_class(void)
{
    check_str(gh_report_mag_class(6), "high", "magnitude 6 is high");
    check_str(gh_report_mag_class(3), "high", "magnitude 3 is high");
    check_str(gh_report_mag_class(2), "mod", "magnitude 2 is mod");
    check_str(gh_report_mag_class(1), "low", "magnitude 1 is low");
    check_str(gh_report_mag_class(0), "info", "magnitude 0 is info");
    check_str(gh_report_mag_class(-1), "info", "negative magnitude is info");
}

static void test_titlecase(void)
{
    gh_arena *a = gh_arena_new(4096);

    check_str(gh_report_titlecase(a, "fast"), "Fast", "single word");
    check_str(gh_report_titlecase(a, "slow_metabolizer"), "Slow Metabolizer",
              "underscore becomes space");
    check_str(gh_report_titlecase(a, "POOR"), "Poor", "upper-case is normalised");
    check_str(gh_report_titlecase(a, "ultra_rapid_metabolizer"),
              "Ultra Rapid Metabolizer", "several underscores");
    check_str(gh_report_titlecase(a, ""), "", "empty string");
    check_str(gh_report_titlecase(a, NULL), "", "NULL is empty");
    /* Python's title() treats digits as word separators too. */
    check_str(gh_report_titlecase(a, "type2_diabetes"), "Type2 Diabetes",
              "digit does not start a new word mid-token");
    check_str(gh_report_titlecase(a, "a_b_c"), "A B C", "single letters");

    gh_arena_free(a);
}

static void test_db_links(void)
{
    gh_strbuf sb;

    gh_sb_init(&sb);
    gh_report_db_links(&sb, "rs762551");
    check_has(sb.data, "https://www.ncbi.nlm.nih.gov/snp/rs762551", "dbSNP link");
    check_has(sb.data, "clinvar/?term=rs762551", "ClinVar link");
    check_has(sb.data, "snpedia.com/index.php/rs762551", "SNPedia link");
    check_has(sb.data, "pharmgkb.org/search?query=rs762551", "PharmGKB link");
    check_has(sb.data, "&middot;", "links are separated");
    check_has(sb.data, "rel=\"noopener\"", "links are rel=noopener");
    gh_sb_free(&sb);

    /* Non-rsIDs produce nothing at all, as in db_links_html. */
    gh_sb_init(&sb);
    gh_report_db_links(&sb, "chr1_12345");
    check(sb.len == 0, "positional ID yields no links");
    gh_report_db_links(&sb, "");
    check(sb.len == 0, "empty ID yields no links");
    gh_report_db_links(&sb, NULL);
    check(sb.len == 0, "NULL ID yields no links");
    gh_sb_free(&sb);
}

static void test_paper_refs(void)
{
    gh_strbuf sb;

    /* rs762551 has a reference in the generated table. */
    gh_sb_init(&sb);
    gh_report_paper_refs(&sb, "rs762551");
    check(sb.len > 0, "known rsID has references");
    check_has(sb.data, "pubmed.ncbi.nlm.nih.gov/", "reference links to PubMed");
    check_has(sb.data, "paper-refs", "references carry their CSS class");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_paper_refs(&sb, "rs00000000");
    check(sb.len == 0, "unknown rsID yields no references");
    gh_report_paper_refs(&sb, NULL);
    check(sb.len == 0, "NULL rsID yields no references");
    gh_sb_free(&sb);
}

static void test_eli5(void)
{
    check(*gh_report_eli5_gene("MTHFR") != '\0', "MTHFR has an ELI5 entry");
    check_str(gh_report_eli5_gene("NOT_A_GENE"), "", "unknown gene has none");
    check_str(gh_report_eli5_gene(NULL), "", "NULL gene has none");
}

/* ------------------------------------------------------------------ */
/* Template filling                                                    */
/* ------------------------------------------------------------------ */

static void test_fill_template(void)
{
    gh_strbuf sb;
    static const char *const names[] = {"a", "bb", "empty"};
    static const char *const values[] = {"1", "22", ""};

    gh_sb_init(&sb);
    gh_report_fill_template(&sb, "x{a}y{bb}z", names, values, 3);
    check_str(sb.data, "x1y22z", "basic substitution");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_fill_template(&sb, "{a}{a}{bb}", names, values, 3);
    check_str(sb.data, "1122", "repeated and adjacent placeholders");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_fill_template(&sb, "[{empty}]", names, values, 3);
    check_str(sb.data, "[]", "empty value substitutes to nothing");
    gh_sb_free(&sb);

    /* An unknown placeholder stays visible rather than blanking silently. */
    gh_sb_init(&sb);
    gh_report_fill_template(&sb, "{unknown} {a}", names, values, 3);
    check_str(sb.data, "{unknown} 1", "unknown placeholder is preserved");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_fill_template(&sb, "{a", names, values, 3);
    check_str(sb.data, "{a", "unterminated brace is literal");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_fill_template(&sb, "no placeholders", names, values, 3);
    check_str(sb.data, "no placeholders", "plain text passes through");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_fill_template(&sb, "", names, values, 3);
    check(sb.len == 0, "empty template");
    gh_sb_free(&sb);

    /* CSS braces in the real template must survive untouched. */
    gh_sb_init(&sb);
    gh_report_fill_template(&sb, "body { color: red; }", names, values, 3);
    check_str(sb.data, "body { color: red; }", "CSS braces are not eaten");
    gh_sb_free(&sb);
}

static void test_template_parts(void)
{
    gh_arena *a = gh_arena_new(65536);

    const char *tpl = gh_report_template(a);
    check(tpl != NULL && strlen(tpl) > 10000, "template joins to a full document");
    check_has(tpl, "<!DOCTYPE html>", "template starts as HTML5");
    check_has(tpl, "{clinical_detail_content}", "template exposes its slots");
    check_has(tpl, "prefers-color-scheme: dark", "dark mode CSS survives");

    /* Every advertised slot must actually appear in the template. */
    for (size_t i = 0; GH_TEMPLATE_SLOTS[i]; i++) {
        char needle[128];
        snprintf(needle, sizeof(needle), "{%s}", GH_TEMPLATE_SLOTS[i]);
        check(strstr(tpl, needle) != NULL, "declared slot is present in template");
    }

    /* Joining twice must give identical text. */
    const char *again = gh_report_template(a);
    check(strcmp(tpl, again) == 0, "template join is deterministic");

    gh_arena_free(a);
}

/* ------------------------------------------------------------------ */
/* Fixtures                                                            */
/* ------------------------------------------------------------------ */

static const double FREQ_EXAMPLE[GH_POP_COUNT] = {0.15, 0.18, 0.30, 0.34, 0.0005};

static gh_finding SAMPLE[] = {
    {"rs762551", "CYP1A2", "Drug Metabolism", "AA", "fast",
     "Fast caffeine metabolizer", "", NULL, 1},
    {"rs4244285", "CYP2C19", "Drug Metabolism", "AA", "poor_metabolizer",
     "Poor CYP2C19", "clopidogrel caveat", FREQ_EXAMPLE, 4},
    {"rs1801133", "MTHFR", "Methylation", "TT", "reduced",
     "Reduced MTHFR activity", "", NULL, 3},
    {"rs4988235", "MCM6", "Nutrition", "GG", "intolerant",
     "Likely lactose intolerant", "", NULL, 2},
    {"rs1815739", "ACTN3", "Fitness", "TT", "endurance",
     "Endurance muscle profile", "", NULL, 0},
};

static gh_drug_finding SAMPLE_DRUGS[] = {
    {"rs4244285", "CYP2C19", "clopidogrel", "AA",
     "Reduced response to clopidogrel", "1A", "Efficacy"},
};

static void make_analysis(gh_analysis *a)
{
    memset(a, 0, sizeof(*a));
    a->findings = SAMPLE;
    a->nfindings = sizeof(SAMPLE) / sizeof(SAMPLE[0]);
    a->drug_findings = SAMPLE_DRUGS;
    a->ndrug_findings = sizeof(SAMPLE_DRUGS) / sizeof(SAMPLE_DRUGS[0]);
    a->total_snps = 600000;
    a->analyzed_snps = a->nfindings;
    a->high_impact = 2;     /* magnitudes 4 and 3 */
    a->moderate_impact = 1; /* magnitude 2 */
    a->low_impact = 1;      /* magnitude 1 */
}

/* ------------------------------------------------------------------ */
/* Charts                                                              */
/* ------------------------------------------------------------------ */

static void test_impact_bar(void)
{
    gh_analysis a;
    make_analysis(&a);

    gh_strbuf sb;
    gh_sb_init(&sb);
    gh_report_svg_impact_bar(&sb, &a);

    check_has(sb.data, "<svg viewBox=\"0 0 440 165\"", "impact bar viewBox");
    check_has(sb.data, "Impact Distribution", "impact bar title");
    check_has(sb.data, ">High<", "High label");
    check_has(sb.data, ">Moderate<", "Moderate label");
    check_has(sb.data, ">Low<", "Low label");
    check_has(sb.data, ">Info<", "Info label");
    check_has(sb.data, "role=\"img\"", "chart is labelled for assistive tech");
    /* 2 high, 1 moderate, 1 low, 1 info from the fixture. */
    check_has(sb.data, ">2</text>", "high count rendered");
    gh_sb_free(&sb);

    /* No findings: still a valid chart, with every bar at zero. */
    gh_analysis empty;
    memset(&empty, 0, sizeof(empty));
    gh_sb_init(&sb);
    gh_report_svg_impact_bar(&sb, &empty);
    check_has(sb.data, "<svg", "empty analysis still emits a chart");
    check_has(sb.data, "width=\"0\"", "bars are zero width with no findings");
    check_lacks(sb.data, "nan", "no NaN from dividing by zero");
    gh_sb_free(&sb);
}

static void test_category_donut(void)
{
    gh_arena *arena = gh_arena_new(65536);
    gh_analysis a;
    make_analysis(&a);

    gh_strbuf sb;
    gh_sb_init(&sb);
    gh_report_svg_category_donut(&sb, arena, &a);

    check_has(sb.data, "Findings by Category", "donut title");
    check_has(sb.data, "Drug Metabolism (2)", "largest category first with count");
    check_has(sb.data, "Methylation (1)", "single-finding category");
    check_has(sb.data, ">5</text>", "centre shows the total");
    check_has(sb.data, "<path d=\"M ", "arcs are drawn");
    check_lacks(sb.data, "nan", "no NaN in arc geometry");

    /* Largest slice must be drawn first, matching the Python ordering. */
    const char *drug = strstr(sb.data, "Drug Metabolism");
    const char *meth = strstr(sb.data, "Methylation");
    check(drug && meth && drug < meth, "legend is ordered by descending count");
    gh_sb_free(&sb);

    /* Empty input produces nothing rather than a broken chart. */
    gh_analysis empty;
    memset(&empty, 0, sizeof(empty));
    gh_sb_init(&sb);
    gh_report_svg_category_donut(&sb, arena, &empty);
    check(sb.len == 0, "no findings means no donut");
    gh_sb_free(&sb);

    /* A single category becomes one full sweep; the large-arc flag must be set. */
    gh_finding one[] = {{"rs1", "G", "Solo", "AA", "s", "d", "", NULL, 1}};
    gh_analysis single = {0};
    single.findings = one;
    single.nfindings = 1;
    gh_sb_init(&sb);
    gh_report_svg_category_donut(&sb, arena, &single);
    check_has(sb.data, "Solo (1)", "single category legend");
    check_has(sb.data, " 1 1 ", "full sweep sets the large-arc flag");
    gh_sb_free(&sb);

    gh_arena_free(arena);
}

/* ------------------------------------------------------------------ */
/* Sections                                                            */
/* ------------------------------------------------------------------ */


/* ------------------------------------------------------------------ */
/* Text helpers                                                        */
/* ------------------------------------------------------------------ */

static void test_clean_condition(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    check_str(gh_report_clean_condition(arena,
              "PI S|Alpha-1-antitrypsin deficiency|not provided"),
              "Alpha-1-antitrypsin deficiency", "picks the readable segment");
    check_str(gh_report_clean_condition(arena,
              "HYPERTENSION, DIASTOLIC, RESISTANCE TO|KCNMB1-related disorder"),
              "Hypertension,  Diastolic,  resistance",
              /* The double spaces are real: replace(",", ", ") runs on text
               * that already has a space after each comma. */
              "all-caps is title-cased with the Python's replacements");
    check_str(gh_report_clean_condition(arena, "not provided|see cases"),
              "not provided", "falls back to the first raw segment");
    check_str(gh_report_clean_condition(arena, ""), "Unknown", "empty is Unknown");
    check_str(gh_report_clean_condition(arena, NULL), "Unknown", "null is Unknown");
    check_str(gh_report_clean_condition(arena, "Hemochromatosis type 1"),
              "Hemochromatosis type 1", "plain text passes through");
    gh_arena_free(arena);
}

static void test_clean_why(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    check_str(gh_report_clean_why(arena, "risk A; risk B; risk A"),
              "risk A; risk B", "duplicate phrases dropped, order kept");
    check_str(gh_report_clean_why(arena, "no separators here"),
              "no separators here", "text without ';' is unchanged");
    check_str(gh_report_clean_why(arena, "X|Y-related disorder"),
              "X", "a pipe routes through clean_condition");
    gh_arena_free(arena);
}

static void test_pyfloat(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    check_str(gh_report_pyfloat(arena, 2.8), "2.8", "short decimal");
    check_str(gh_report_pyfloat(arena, 1.0), "1.0", "integral keeps .0");
    check_str(gh_report_pyfloat(arena, 100.0), "100.0", "three digits keep .0");
    check_str(gh_report_pyfloat(arena, 62.5), "62.5", "one decimal");
    check_str(gh_report_pyfloat(arena, 0.05), "0.05", "leading zeros");
    check_str(gh_report_pyfloat(arena, 14.9), "14.9", "no binary noise");
    gh_arena_free(arena);
}

/* ------------------------------------------------------------------ */
/* Sections with nothing to show                                       */
/* ------------------------------------------------------------------ */

static void test_empty_states(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    gh_analysis empty;
    memset(&empty, 0, sizeof(empty));
    gh_strbuf sb;

    gh_sb_init(&sb);
    gh_report_action_plan(&sb, arena, NULL);
    check_str(sb.data, "<p>No personalized recommendations available.</p>", "action plan");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_drug_guide(&sb, arena, &empty, NULL);
    check_str(sb.data, "<p>No drug-gene interaction data available.</p>", "drug guide");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_disease_risk(&sb, arena, NULL);
    check_str(sb.data, "<p>No disease risk data available.</p>", "disease risk");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_mental_health(&sb, arena, NULL);
    check_str(sb.data, "<p>No mental health genetic data available.</p>", "mental health");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_ancestry(&sb, arena, NULL);
    check_str(sb.data, "<p>No ancestry-informative markers found in genome data.</p>",
              "ancestry");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_nutrigenomics(&sb, arena, NULL);
    check_str(sb.data, "<p>No nutrigenomics data available.</p>", "nutrigenomics");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_quality(&sb, NULL);
    check_str(sb.data, "<p>No quality metrics available.</p>", "quality");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_doctor_card(&sb, arena, NULL);
    check_str(sb.data, "<p>No significant clinical findings for doctor review.</p>",
              "doctor card");
    gh_sb_free(&sb);

    /* Key findings and the body profile always open with their intro. */
    gh_sb_init(&sb);
    gh_report_key_findings(&sb, arena, &empty, NULL);
    check_has(sb.data, "Your body has a recipe book called DNA", "key findings intro");
    check_has(sb.data, "Show this report to your doctor!", "key findings callout");
    check_lacks(sb.data, "<h3>", "no headings without data");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_body_profile(&sb, arena, &empty, NULL);
    check_str(sb.data, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "Traits, blood type, chronotype, longevity, and athletic profile.</p>",
              "body profile intro only");
    gh_sb_free(&sb);

    gh_arena_free(arena);
}

/* ------------------------------------------------------------------ */
/* Sections with data                                                  */
/* ------------------------------------------------------------------ */

static void test_clinical_detail(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    gh_analysis a;
    make_analysis(&a);

    gh_strbuf sb;
    gh_sb_init(&sb);
    gh_report_clinical_detail(&sb, arena, &a, NULL);

    check_has(sb.data, "Drug Metabolism <span class=\"badge\">2</span>",
              "category header carries its count");
    check_has(sb.data, "mag-badge mag-high\">4/6", "magnitude badge");
    check_has(sb.data, "<code>rs4244285</code> \xe2\x80\x94 <code>AA</code> \xe2\x80\x94 Poor Metabolizer",
              "header uses literal em dashes and a title-cased status");
    check_has(sb.data, "Note: clopidogrel caveat", "note is rendered");
    check_has(sb.data, "South Asian: 34%", "highest frequency first");
    check_lacks(sb.data, "American: 0%", "frequencies below 0.1% are dropped");
    check_has(sb.data, "Pathway Analysis", "pathway section");
    check_has(sb.data, "Methylation Cycle", "MTHFR pulls in its pathway");
    check_lacks(sb.data, "Carrier Screening", "no carrier block without carriers");
    check_lacks(sb.data, "Epistasis", "no epistasis block without interactions");

    const char *drug = strstr(sb.data, ">Drug Metabolism ");
    const char *fit = strstr(sb.data, ">Fitness ");
    const char *nut = strstr(sb.data, ">Nutrition ");
    check(drug && fit && drug < fit, "categories sorted alphabetically");
    check(fit && nut && fit < nut, "categories sorted alphabetically (2)");
    const char *poor = strstr(sb.data, "rs4244285");
    const char *fast = strstr(sb.data, "rs762551");
    check(poor && fast && poor < fast, "higher magnitude sorts first");
    gh_sb_free(&sb);
    gh_arena_free(arena);
}

static void test_html_escaping(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    gh_finding hostile = {
        "rs1", "<b>G</b>", "Cat & \"Co\"", "A'G", "st<us",
        "desc <script>x</script>", "note & more", NULL, 1,
    };
    gh_analysis a;
    memset(&a, 0, sizeof(a));
    a.findings = &hostile;
    a.nfindings = 1;

    gh_strbuf sb;
    gh_sb_init(&sb);
    gh_report_clinical_detail(&sb, arena, &a, NULL);
    check_has(sb.data, "&lt;b&gt;G&lt;/b&gt;", "gene escaped");
    check_has(sb.data, "A&#x27;G", "apostrophe escaped the way html.escape does");
    check_has(sb.data, "desc &lt;script&gt;x&lt;/script&gt;", "description escaped");
    check_has(sb.data, "Note: note &amp; more", "note escaped");
    check_lacks(sb.data, "<script>", "no raw script tag");
    /* The category heading is not escaped in the Python either. */
    check_has(sb.data, "Cat & \"Co\" <span class=\"badge\">1</span>",
              "category heading is emitted raw, matching the reference");
    gh_sb_free(&sb);
    gh_arena_free(arena);
}

static void test_drug_guide(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    gh_analysis a;
    make_analysis(&a);

    gh_star_result stars[2];
    memset(stars, 0, sizeof(stars));
    stars[0].gene = "CYP2D6"; stars[0].diplotype = "*4/*4";
    stars[0].phenotype = "poor"; stars[0].clinical_note = "x";
    stars[0].snps_found = 3; stars[0].snps_total = 4;
    stars[1].gene = "NAT2"; stars[1].diplotype = "*4/*5";
    stars[1].phenotype = "intermediate"; stars[1].clinical_note = "y";
    gh_report_scorers sc;
    memset(&sc, 0, sizeof(sc));
    sc.stars = stars;
    sc.nstars = 2;

    gh_strbuf sb;
    gh_sb_init(&sb);
    gh_report_drug_guide(&sb, arena, &a, &sc);
    check_has(sb.data, "Your Drug-Processing Enzymes", "enzyme section");
    check_has(sb.data, "<title>CYP2D6: Slow</title>", "poor maps to the Slow gauge label");
    check_has(sb.data, "<title>NAT2: Unknown</title>",
              "half-step levels miss the label table, as in Python");
    check_has(sb.data, "<td>3/4</td>", "SNP coverage cell");
    check_has(sb.data, "Drug-Gene Interactions (PharmGKB)", "PharmGKB table");
    check_has(sb.data, "<code>rs4244285</code> <span class=\"db-links\">", "rsID with links");
    check_has(sb.data, "id=\"drug-checker\"", "interactive checker present");
    check_lacks(sb.data, "Drug Combination Warnings", "no polypharmacy block without warnings");
    gh_sb_free(&sb);
    gh_arena_free(arena);
}

static void test_references(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    gh_analysis a;
    make_analysis(&a);

    gh_strbuf sb;
    gh_sb_init(&sb);
    gh_report_references(&sb, &a);
    check_has(sb.data, "<h3>Key Papers</h3><ul>", "key papers list");
    /* The warfarin guideline is listed for two rsIDs but printed once, and
     * first, because PAPER_REFS is walked in insertion order. */
    const char *first_li = strstr(sb.data, "<li>");
    check(first_li && strstr(first_li, "28198005") == strstr(sb.data, "28198005"),
          "first paper is the warfarin guideline");
    const char *second = strstr(sb.data, "28198005");
    check(second && !strstr(second + 1, "PMID: 28198005, 2017)</li>\n<li><a href=\"https://pubmed.ncbi.nlm.nih.gov/28198005/"),
          "duplicate PMIDs collapse to one entry");
    check_has(sb.data, "Database Links for All Analyzed rsIDs", "rsID grid");
    /* String order: rs1801133 < rs1815739 < rs4244285 < rs4988235 < rs762551. */
    const char *r1 = strstr(sb.data, "<div><code>rs1801133</code>");
    const char *r2 = strstr(sb.data, "<div><code>rs4988235</code>");
    const char *r3 = strstr(sb.data, "<div><code>rs762551</code>");
    check(r1 && r2 && r3 && r1 < r2 && r2 < r3, "rsIDs sorted as strings");
    check_has(sb.data, "Methodology &amp; Disclaimers", "methodology list");
    gh_sb_free(&sb);
    gh_arena_free(arena);
}

static void test_doctor_card(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    gh_star_result star;
    memset(&star, 0, sizeof(star));
    star.gene = "CYP2C19"; star.diplotype = "*1/*2"; star.phenotype = "intermediate";
    gh_apoe_result apoe;
    memset(&apoe, 0, sizeof(apoe));
    apoe.apoe_type = "e3/e4"; apoe.risk_level = "elevated";
    gh_report_scorers sc;
    memset(&sc, 0, sizeof(sc));
    sc.stars = &star; sc.nstars = 1; sc.apoe = &apoe;

    gh_strbuf sb;
    gh_sb_init(&sb);
    gh_report_doctor_card(&sb, arena, &sc);
    check_has(sb.data, "Patient Genetic Summary", "card heading");
    check_has(sb.data, "<td><code>*1/*2</code></td><td>Intermediate</td>", "star row");
    check_has(sb.data, "<strong>e3/e4</strong> \xe2\x80\x94 <span class=\"mag-badge\" "
                       "style=\"background:#c47a2b;color:#fff\">Elevated Risk</span>",
              "APOE line with the orange elevated badge");
    check_lacks(sb.data, "High-Priority Conditions", "no priorities without recommendations");
    gh_sb_free(&sb);
    gh_arena_free(arena);
}

static void test_charts(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);
    gh_strbuf sb;

    gh_sb_init(&sb);
    gh_report_svg_prs_gauge(&sb, arena, "T2D", 87.4, "elevated");
    check_has(sb.data, "<title>T2D: 87th percentile (elevated)</title>", "gauge title");
    check_has(sb.data, ">87%</text>", "gauge percentage");
    check_has(sb.data, ">Elevated</text>", "gauge category title-cased");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_report_svg_metabolism_gauge(&sb, "CYP2D6", 2, "#000");
    check_has(sb.data, "<title>CYP2D6: Fast</title>", "metabolism gauge label");
    gh_sb_free(&sb);

    gh_prs_result prs[2];
    memset(prs, 0, sizeof(prs));
    prs[0].name = "Age-Related Macular Degeneration"; prs[0].percentile = 10;
    prs[0].risk_category = "low";
    prs[1].name = "Systemic Lupus Erythematosus"; prs[1].percentile = 96;
    prs[1].risk_category = "high";
    gh_sb_init(&sb);
    gh_report_svg_risk_heatmap(&sb, arena, prs, 2);
    check_has(sb.data, "viewBox=\"0 0 172 68\"", "two cells, one row");
    check_has(sb.data, ">Lupus (SLE)</text>", "name shortening");
    check_has(sb.data, ">AMD</text>", "AMD shortening");
    /* Highest percentile first. */
    check(strstr(sb.data, "Lupus") < strstr(sb.data, "AMD"), "sorted by percentile desc");
    gh_sb_free(&sb);

    gh_arena_free(arena);
}

static void test_supplement_priorities(void)
{
    gh_nutrient_need needs[3];
    memset(needs, 0, sizeof(needs));
    /* Needs arrive sorted by severity; the priorities list must follow
     * profile-table order within each band instead. */
    needs[0].profile = &GH_NUTRIENT_PROFILES[2]; needs[0].need_level = "moderate";
    needs[1].profile = &GH_NUTRIENT_PROFILES[0]; needs[1].need_level = "high";
    needs[2].profile = &GH_NUTRIENT_PROFILES[1]; needs[2].need_level = "moderate";
    gh_nutrigenomics_result nut;
    memset(&nut, 0, sizeof(nut));
    nut.needs = needs;
    nut.nneeds = 3;

    const gh_nutrient_need *out[8];
    size_t n = gh_report_supplement_priorities(&nut, out, 8);
    check(n == 3, "all three qualify");
    check(n == 3 && out[0]->profile == &GH_NUTRIENT_PROFILES[0], "high first");
    check(n == 3 && out[1]->profile == &GH_NUTRIENT_PROFILES[1],
          "moderates follow in profile order");
    check(n == 3 && out[2]->profile == &GH_NUTRIENT_PROFILES[2], "moderates (2)");
}

static const char *write_temp(const char *name, const char *content)
{
    static char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_report_test_%s", name);
    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "  FAIL: cannot write %s\n", path);
        failures++;
        return path;
    }
    fwrite(content, 1, strlen(content), f);
    fclose(f);
    return path;
}

static void test_quality_metrics(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);

    const char *path = write_temp("genome.txt",
        "# comment\n"
        "rs1\t1\t100\tAA\n"     /* autosomal homozygous */
        "rs2\t1\t200\tAG\n"     /* autosomal heterozygous */
        "rs3\t2\t300\tCT\n"     /* autosomal heterozygous */
        "rs4\tMT\t400\tG\n"     /* mitochondrial, haploid */
        "rs5\tY\t500\tC\n"      /* Y present -> inferred male */
        "rs6\tX\t600\tTT\n"
        "rs7\t1\t700\t--\n"     /* no-call */
        "rs8\t1\t800\t--\n");   /* no-call */

    gh_genome g;
    check(gh_genome_load(&g, arena, path, NULL), "quality fixture loads");

    gh_quality q;
    gh_quality_compute(&q, arena, &g, path);

    check(q.total_snps == 6, "no-calls excluded from total");
    check(q.no_call_count == 2, "no-calls counted from the raw file");
    check(q.autosomal_count == 3, "autosomal count");
    check(q.mt_snp_count == 1, "mitochondrial count");
    check(q.has_mt, "MT detected");
    check(q.has_y, "Y detected");
    /* 2 heterozygous of 3 autosomal. */
    check(q.het_rate > 0.66 && q.het_rate < 0.67, "heterozygosity rate");
    /* 6 of 8 called. */
    check(q.call_rate > 0.749 && q.call_rate < 0.751, "call rate");
    /* chr1 carries four rows but two are no-calls, which never enter the
     * genome, so only two are counted. Cross-checked against the Python. */
    check(gh_quality_chrom_count(&q, "1") == 2, "per-chromosome count");
    check(gh_quality_chrom_count(&q, "ZZ") == 0, "absent chromosome is zero");

    /* Chromosomes sort numerically first, then the rest alphabetically. */
    check(q.nchromosomes == 5, "distinct chromosomes");
    check_str(q.chromosomes[0].name, "1", "chromosome 1 first");
    check_str(q.chromosomes[1].name, "2", "chromosome 2 second");
    check_str(q.chromosomes[2].name, "MT", "MT after the numbers");

    /* Without the raw file there are no no-calls to count. */
    gh_quality q2;
    gh_quality_compute(&q2, arena, &g, NULL);
    check(q2.no_call_count == 0, "no-calls not counted without a path");
    check(q2.call_rate > 0.999, "call rate is 1.0 without a path");

    /* "chr" prefixes normalise away. */
    const char *prefixed = write_temp("prefixed.txt",
        "rs1\tchr1\t100\tAA\nrs2\tCHR1\t200\tAG\nrs3\tChr2\t300\tCC\n");
    gh_genome g2;
    check(gh_genome_load(&g2, arena, prefixed, NULL), "prefixed fixture loads");
    gh_quality q3;
    gh_quality_compute(&q3, arena, &g2, prefixed);
    check(gh_quality_chrom_count(&q3, "1") == 2, "chr prefix normalised");
    check(gh_quality_chrom_count(&q3, "2") == 1, "mixed-case prefix normalised");
    check(q3.autosomal_count == 3, "prefixed chromosomes count as autosomal");

    /* An empty genome must not divide by zero. */
    const char *empty_path = write_temp("empty.txt", "# only a comment\n");
    gh_genome g3;
    check(gh_genome_load(&g3, arena, empty_path, NULL), "empty fixture loads");
    gh_quality q4;
    gh_quality_compute(&q4, arena, &g3, empty_path);
    check(q4.total_snps == 0, "empty genome has no SNPs");
    check(q4.call_rate == 0.0, "empty genome call rate is zero, not NaN");
    check(q4.het_rate == 0.0, "empty genome het rate is zero, not NaN");

    gh_arena_free(arena);
}

static void test_quality_section(void)
{
    gh_quality q;
    memset(&q, 0, sizeof(q));
    q.total_snps = 601234; q.no_call_count = 5000; q.call_rate = 0.9917;
    q.autosomal_count = 580000; q.mt_snp_count = 12; q.het_rate = 0.3456;
    gh_chrom_count chroms[] = {{"1", 50000}, {"2", 40000}, {"X", 20000}};
    q.chromosomes = chroms; q.nchromosomes = 3;

    gh_strbuf sb;
    gh_sb_init(&sb);
    gh_report_quality(&sb, &q);
    check_has(sb.data, "99.2% \xe2\x80\x94 Excellent", "call-rate badge with a literal em dash");
    check_has(sb.data, "<td>601,234</td>", "grouped total");
    check_has(sb.data, "<td>0.346</td>", "het rate at three decimals");
    check_has(sb.data, "Female", "no Y means Female");
    check_has(sb.data, "SNPs per Chromosome", "coverage chart");
    check_lacks(sb.data, ">X</text>", "X is not an autosome");
    gh_sb_free(&sb);
}

static void test_render_full(void)
{
    gh_arena *arena = gh_arena_new(1 << 20);
    gh_analysis a;
    make_analysis(&a);
    gh_quality q;
    memset(&q, 0, sizeof(q));
    q.total_snps = 5; q.call_rate = 1.0;

    gh_report_opts opts = {"Ada <Lovelace>", "2024-01-01 00:00"};
    gh_strbuf sb;
    gh_sb_init(&sb);
    gh_report_render(&sb, arena, &a, &q, NULL, &opts);

    check_has(sb.data, "<title>Genetic Health Report \xe2\x80\x94 Ada &lt;Lovelace&gt;</title>",
              "subject title uses an em dash and is escaped");
    check_has(sb.data, "<span>600,000 SNPs analyzed</span>", "total_snps:, placeholder filled");
    /* The CSS has braces of its own, so look for the slot names instead. */
    check_lacks(sb.data, "_content}", "no unsubstituted section placeholders");
    check_lacks(sb.data, "{total_snps", "total_snps placeholder filled");
    check_lacks(sb.data, "{generated_date}", "date placeholder filled");
    check_lacks(sb.data, "{subject_title}", "title placeholder filled");
    check_lacks(sb.data, "not yet available in the C port", "no placeholder sections remain");
    check_has(sb.data, "2024-01-01 00:00", "date override");
    gh_sb_free(&sb);

    /* Low coverage banner appears under 100k SNPs. */
    a.total_snps = 1200;
    gh_sb_init(&sb);
    gh_report_render(&sb, arena, &a, &q, NULL, &opts);
    check_has(sb.data, "Only 1,200 SNPs loaded (0% of a typical", "low coverage banner");
    gh_sb_free(&sb);

    /* Deterministic. */
    gh_strbuf x, y;
    gh_sb_init(&x); gh_sb_init(&y);
    gh_report_render(&x, arena, &a, &q, NULL, &opts);
    gh_report_render(&y, arena, &a, &q, NULL, &opts);
    check(x.len == y.len && memcmp(x.data, y.data, x.len) == 0, "two renders are identical");
    gh_sb_free(&x); gh_sb_free(&y);

    gh_arena_free(arena);
}

int main(void)
{
    test_mag_class();
    test_titlecase();
    test_db_links();
    test_paper_refs();
    test_eli5();
    test_fill_template();
    test_template_parts();
    test_impact_bar();
    test_category_donut();
    test_clean_condition();
    test_clean_why();
    test_pyfloat();
    test_empty_states();
    test_clinical_detail();
    test_html_escaping();
    test_drug_guide();
    test_references();
    test_doctor_card();
    test_charts();
    test_supplement_priorities();
    test_quality_metrics();
    test_quality_section();
    test_render_full();

    printf("%s: %d checks, %d failures\n",
           failures ? "FAILED" : "ok", checks, failures);
    return failures ? 1 : 0;
}

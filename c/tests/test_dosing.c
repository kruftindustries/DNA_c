/* Unit tests for drug dosing, polypharmacy and the preventive care timeline.
 *
 * These three modules never touch the genome: they read star-allele
 * phenotypes, lifestyle findings, polygenic scores, APOE and ACMG findings.
 * So the fixtures here are those result structures, built directly. */

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "gh_dosing.h"
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
        fprintf(stderr, "  FAIL: %s\n    missing: %s\n    in: %s\n",
                what, needle, haystack ? haystack : "(null)");
    }
}

/* A star-allele result set with every gene at "normal", so a test only has to
 * name the genes it wants to move. */
static size_t normal_stars(gh_star_result *out)
{
    for (size_t i = 0; i < GH_STAR_GENE_COUNT; i++) {
        memset(&out[i], 0, sizeof(out[i]));
        out[i].gene = GH_STAR_GENES[i].gene;
        out[i].diplotype = "*1/*1";
        out[i].phenotype = "normal";
        out[i].confidence = "high";
    }
    return GH_STAR_GENE_COUNT;
}

static void set_phenotype(gh_star_result *stars, size_t n, const char *gene,
                          const char *phenotype)
{
    for (size_t i = 0; i < n; i++)
        if (strcmp(stars[i].gene, gene) == 0)
            stars[i].phenotype = phenotype;
}

static const gh_dose_rec *find_rec(const gh_dosing_result *r, const char *drug)
{
    for (size_t i = 0; i < r->nrecommendations; i++)
        if (strcmp(r->recommendations[i].drug->drug, drug) == 0)
            return &r->recommendations[i];
    return NULL;
}

static const gh_screening *find_screening(const gh_preventive_result *r,
                                          const char *test)
{
    for (size_t i = 0; i < r->ntimeline; i++)
        if (strcmp(r->timeline[i].test, test) == 0)
            return &r->timeline[i];
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Drug dosing                                                         */
/* ------------------------------------------------------------------ */

static void test_critical_keywords(void)
{
    check(gh_dose_is_critical("CONTRAINDICATED at full dose"), "contraindicated");
    check(gh_dose_is_critical("can be FATAL"), "fatal");
    check(gh_dose_is_critical("AVOID codeine"), "avoid upper");
    /* Python upper-cases the action before testing, so the stored casing of
     * the keyword does not matter. */
    check(gh_dose_is_critical("Avoid simvastatin >20mg"), "avoid mixed case");
    check(gh_dose_is_critical("life-threatening myelosuppression"),
          "life-threatening lower case");
    check(!gh_dose_is_critical("Standard dosing is appropriate."), "plain text");
    check(!gh_dose_is_critical(""), "empty");
    check(!gh_dose_is_critical(NULL), "null");
    /* Python tests plain substring containment, so a keyword embedded in a
     * longer word still counts. Reproduced rather than "improved". */
    check(gh_dose_is_critical("FATALITY rates"), "substring, not word, match");
    check(!gh_dose_is_critical("FATA L"), "split keyword does not match");
}

static void test_dosing_no_data(gh_arena *arena)
{
    gh_dosing_result r;
    gh_generate_drug_dosing(&r, arena, NULL, 0, NULL);
    check(r.nrecommendations == 0, "no recommendations without star alleles");
    check_str(r.summary, "No pharmacogenomic data available for drug dosing.",
              "no-data summary");
}

static void test_dosing_all_normal(gh_arena *arena)
{
    gh_star_result stars[GH_STAR_GENE_COUNT];
    size_t n = normal_stars(stars);

    gh_dosing_result r;
    gh_generate_drug_dosing(&r, arena, stars, n, NULL);

    /* CYP1A2 normal still fires the caffeine rule, and CYP3A5 normal the
     * tacrolimus expresser rule, so an all-normal profile is not empty. */
    check(find_rec(&r, "caffeine") != NULL, "caffeine rule fires on normal");
    check(find_rec(&r, "tacrolimus") != NULL, "tacrolimus expresser on normal");
    check(find_rec(&r, "codeine") == NULL, "no codeine rule on normal");
    check(r.nwarnings == 0, "all-normal raises no warnings");
    check_has(r.summary, "based on your pharmacogenomic profile.",
              "recommendations-only summary");
}

static void test_dosing_poor_metabolizer(gh_arena *arena)
{
    gh_star_result stars[GH_STAR_GENE_COUNT];
    size_t n = normal_stars(stars);
    set_phenotype(stars, n, "CYP2D6", "ultrarapid");
    set_phenotype(stars, n, "DPYD", "poor");

    gh_dosing_result r;
    gh_generate_drug_dosing(&r, arena, stars, n, NULL);

    const gh_dose_rec *codeine = find_rec(&r, "codeine");
    check(codeine != NULL, "ultrarapid CYP2D6 fires a codeine rule");
    if (codeine) {
        check_has(codeine->rule->action, "AVOID codeine", "codeine action");
        check(codeine->critical, "codeine ultrarapid is critical");
        check_str(codeine->drug->display, "Codeine", "display name");
    }

    const gh_dose_rec *dpyd = find_rec(&r, "fluoropyrimidines");
    check(dpyd != NULL, "DPYD poor fires a chemotherapy rule");
    if (dpyd)
        check(dpyd->critical, "DPYD poor is critical");

    check(r.nwarnings >= 2, "both criticals reach the warning list");
    check_has(r.summary, "critical warning(s).", "warning summary");

    /* warnings must point into recommendations, not copy them. */
    bool aliased = false;
    for (size_t i = 0; i < r.nwarnings; i++)
        for (size_t k = 0; k < r.nrecommendations; k++)
            if (r.warnings[i] == &r.recommendations[k])
                aliased = true;
    check(aliased, "warnings alias the recommendation entries");
}

static void test_dosing_reads_lifestyle_findings(gh_arena *arena)
{
    gh_star_result stars[GH_STAR_GENE_COUNT];
    size_t n = normal_stars(stars);

    /* The second warfarin rule is the only one that reads a lifestyle
     * finding rather than a star allele. */
    gh_finding findings[2];
    memset(findings, 0, sizeof(findings));
    findings[0].gene = "VKORC1";
    findings[0].status = "sensitive";
    findings[1].gene = "MTHFR";
    findings[1].status = "reduced";
    gh_analysis analysis;
    memset(&analysis, 0, sizeof(analysis));
    analysis.findings = findings;
    analysis.nfindings = 2;

    gh_dosing_result with;
    gh_generate_drug_dosing(&with, arena, stars, n, &analysis);
    const gh_dose_rec *warfarin = find_rec(&with, "warfarin");
    check(warfarin != NULL, "VKORC1 sensitivity fires the warfarin rule");
    if (warfarin)
        check_has(warfarin->rule->action, "VKORC1 sensitivity detected",
                  "warfarin finding action");

    gh_dosing_result without;
    gh_generate_drug_dosing(&without, arena, stars, n, NULL);
    check(find_rec(&without, "warfarin") == NULL,
          "no warfarin rule without the finding");
}

static void test_dosing_requires_gene_data(gh_arena *arena)
{
    /* A drug is skipped entirely unless one of its genes was called -- even
     * when a lifestyle finding would have satisfied one of its rules. */
    gh_star_result stars[1];
    memset(stars, 0, sizeof(stars));
    stars[0].gene = "TPMT";
    stars[0].phenotype = "normal";

    gh_finding finding;
    memset(&finding, 0, sizeof(finding));
    finding.gene = "VKORC1";
    finding.status = "sensitive";
    gh_analysis analysis;
    memset(&analysis, 0, sizeof(analysis));
    analysis.findings = &finding;
    analysis.nfindings = 1;

    gh_dosing_result r;
    gh_generate_drug_dosing(&r, arena, stars, 1, &analysis);
    check(find_rec(&r, "warfarin") == NULL,
          "warfarin skipped when neither CYP2C9 nor VKORC1 was called");
    check_str(r.summary, "No drug dosing adjustments needed based on "
                         "available pharmacogenomic data.", "empty summary");
}

/* ------------------------------------------------------------------ */
/* Polypharmacy                                                        */
/* ------------------------------------------------------------------ */

static bool poly_severity_known(const char *severity)
{
    return strcmp(severity, "high") == 0 || strcmp(severity, "moderate") == 0 ||
           strcmp(severity, "low") == 0;
}

static const gh_poly_warning *find_poly(const gh_polypharmacy_result *r,
                                        const char *id)
{
    for (size_t i = 0; i < r->nwarnings; i++)
        if (strcmp(r->warnings[i].rule->id, id) == 0)
            return &r->warnings[i];
    return NULL;
}

static void test_polypharmacy_single_gene(gh_arena *arena)
{
    gh_star_result stars[GH_STAR_GENE_COUNT];
    size_t n = normal_stars(stars);
    set_phenotype(stars, n, "TPMT", "poor");

    gh_polypharmacy_result r;
    gh_assess_polypharmacy(&r, arena, stars, n, NULL);

    const gh_poly_warning *w = find_poly(&r, "thiopurine_toxicity");
    check(w != NULL, "TPMT poor fires the thiopurine rule");
    if (w) {
        check_str(w->matched[0], "poor", "matched phenotype recorded");
        check_str(w->rule->severity, "high", "severity");
    }
}

static void test_polypharmacy_needs_every_gene(gh_arena *arena)
{
    gh_star_result stars[GH_STAR_GENE_COUNT];
    size_t n = normal_stars(stars);
    set_phenotype(stars, n, "CYP2C9", "poor");

    /* warfarin_compound needs CYP2C9 *and* VKORC1; the star allele alone is
     * not enough. */
    gh_polypharmacy_result alone;
    gh_assess_polypharmacy(&alone, arena, stars, n, NULL);
    check(find_poly(&alone, "warfarin_compound") == NULL,
          "two-gene rule does not fire on one gene");

    gh_finding finding;
    memset(&finding, 0, sizeof(finding));
    finding.gene = "VKORC1";
    finding.status = "sensitive";
    gh_analysis analysis;
    memset(&analysis, 0, sizeof(analysis));
    analysis.findings = &finding;
    analysis.nfindings = 1;

    gh_polypharmacy_result both;
    gh_assess_polypharmacy(&both, arena, stars, n, &analysis);
    const gh_poly_warning *w = find_poly(&both, "warfarin_compound");
    check(w != NULL, "two-gene rule fires when both match");
    if (w) {
        check_str(w->matched[0], "poor", "CYP2C9 phenotype");
        check_str(w->matched[1], "sensitive", "VKORC1 status");
    }
}

static void test_polypharmacy_unknown_is_not_a_phenotype(gh_arena *arena)
{
    gh_star_result stars[GH_STAR_GENE_COUNT];
    size_t n = normal_stars(stars);
    set_phenotype(stars, n, "CYP2C19", "Unknown");

    gh_polypharmacy_result r;
    gh_assess_polypharmacy(&r, arena, stars, n, NULL);
    check(find_poly(&r, "clopidogrel_resistance") == NULL,
          "Unknown does not satisfy a rule");
}

static void test_polypharmacy_finding_overrides_star(gh_arena *arena)
{
    /* The lookup is built from the star alleles and then overwritten by the
     * lifestyle findings, so a gene named by both is judged on the finding --
     * and among findings the last one wins. */
    gh_star_result stars[GH_STAR_GENE_COUNT];
    size_t n = normal_stars(stars);
    set_phenotype(stars, n, "CYP2D6", "poor");

    gh_finding findings[2];
    memset(findings, 0, sizeof(findings));
    findings[0].gene = "CYP2D6";
    findings[0].status = "poor";
    findings[1].gene = "CYP2D6";
    findings[1].status = "normal";
    gh_analysis analysis;
    memset(&analysis, 0, sizeof(analysis));
    analysis.findings = findings;
    analysis.nfindings = 2;

    gh_finding one = findings[0];
    gh_analysis just_poor;
    memset(&just_poor, 0, sizeof(just_poor));
    just_poor.findings = &one;
    just_poor.nfindings = 1;

    gh_polypharmacy_result kept;
    gh_assess_polypharmacy(&kept, arena, stars, n, &just_poor);
    check(find_poly(&kept, "opioid_sensitivity") == NULL,
          "opioid rule still needs OPRM1");

    gh_polypharmacy_result overridden;
    gh_assess_polypharmacy(&overridden, arena, stars, n, &analysis);
    check(find_poly(&overridden, "cyp2d6_ultrarapid_codeine") == NULL,
          "last finding wins over the star-allele call");
}

static void test_polypharmacy_severity_order(gh_arena *arena)
{
    gh_star_result stars[GH_STAR_GENE_COUNT];
    size_t n = normal_stars(stars);
    /* One moderate rule and one high rule, with the moderate gene listed
     * first in the rule table, so the sort has something to do. */
    set_phenotype(stars, n, "CYP2C19", "ultrarapid");   /* ppi: moderate */
    set_phenotype(stars, n, "DPYD", "poor");            /* chemo: high */

    gh_polypharmacy_result r;
    gh_assess_polypharmacy(&r, arena, stars, n, NULL);
    check(r.nwarnings >= 2, "both rules fire");
    if (r.nwarnings >= 2) {
        check_str(r.warnings[0].rule->severity, "high", "high sorts first");
        check_str(r.warnings[r.nwarnings - 1].rule->severity, "moderate",
                  "moderate sorts last");
    }
    check(r.nseverities == 2, "two distinct severities");
    if (r.nseverities == 2) {
        check_str(r.severities[0], "high", "by_severity order follows sort");
        check_str(r.severities[1], "moderate", "second severity");
        size_t total = r.severity_counts[0] + r.severity_counts[1];
        check(total == r.nwarnings, "by_severity counts cover every warning");
    }
}

/* ------------------------------------------------------------------ */
/* Preventive care                                                     */
/* ------------------------------------------------------------------ */

static void test_priority_rank(void)
{
    check(gh_screening_priority_rank("urgent") == 0, "urgent rank");
    check(gh_screening_priority_rank("high") == 1, "high rank");
    check(gh_screening_priority_rank("elevated") == 2, "elevated rank");
    check(gh_screening_priority_rank("ongoing") == 3, "ongoing rank");
    check(gh_screening_priority_rank("standard") == 4, "standard rank");
    check(gh_screening_priority_rank("invented") == 5, "unknown rank");
    check(gh_screening_priority_rank(NULL) == 5, "null rank");
}

static void test_preventive_baseline(gh_arena *arena)
{
    gh_preventive_result r;
    gh_generate_preventive_timeline(&r, arena, NULL, 0, NULL, NULL, NULL, 0);

    check(r.ntimeline == GH_BASE_SCREENING_COUNT, "baseline is the base table");
    check(r.early_screenings == 0, "nothing moved earlier");
    check_str(r.summary, "Your genetic profile does not indicate need for "
                         "earlier screenings beyond standard guidelines.",
              "baseline summary");

    for (size_t i = 1; i < r.ntimeline; i++)
        check(r.timeline[i - 1].start_age <= r.timeline[i].start_age,
              "timeline sorted by start age");
    for (size_t i = 0; i < r.ntimeline; i++) {
        check(r.timeline[i].genetic_basis == NULL, "base entry has no basis");
        check_str(r.timeline[i].priority, "standard", "base entry priority");
        check_has(r.timeline[i].reason, "General population guideline",
                  "base entry reason");
    }
}

static void test_preventive_prs_modifies_base(gh_arena *arena)
{
    gh_prs_result prs[2];
    memset(prs, 0, sizeof(prs));
    prs[0].id = "type2_diabetes";
    prs[0].name = "Type 2 Diabetes";
    prs[0].risk_category = "high";
    prs[0].percentile = 97.4;
    prs[1].id = "breast_cancer";
    prs[1].name = "Breast Cancer";
    prs[1].risk_category = "elevated";
    prs[1].percentile = 88.0;

    gh_preventive_result r;
    gh_generate_preventive_timeline(&r, arena, prs, 2, NULL, NULL, NULL, 0);

    check(r.ntimeline == GH_BASE_SCREENING_COUNT,
          "modified conditions stay in place rather than being added");
    check(r.early_screenings == 2, "two screenings moved");

    const gh_screening *glucose = find_screening(&r, "Fasting glucose / HbA1c");
    check(glucose != NULL, "glucose screening present");
    if (glucose) {
        check(glucose->start_age == 20, "35 - 15 = 20");
        check_str(glucose->frequency, "Every 6 months", "high frequency");
        check_str(glucose->priority, "high", "high priority");
        check_str(glucose->genetic_basis, "PRS 97th percentile (high)",
                  "percentile formatted with no decimals");
    }

    const gh_screening *mammo = find_screening(&r, "Mammography");
    check(mammo != NULL, "mammography present");
    if (mammo) {
        check(mammo->start_age == 40, "50 - 10 = 40");
        check_str(mammo->priority, "elevated", "elevated priority");
        check_str(mammo->genetic_basis, "PRS 88th percentile (elevated)",
                  "elevated basis");
    }
}

static void test_preventive_start_age_floor(gh_arena *arena)
{
    /* Blood pressure starts at 18 and the high modifier is -15, so the
     * max(18, ...) floor is what keeps it at 18 rather than 3. */
    gh_prs_result prs;
    memset(&prs, 0, sizeof(prs));
    prs.id = "hypertension";
    prs.name = "Hypertension";
    prs.risk_category = "high";
    prs.percentile = 99.9;

    gh_preventive_result r;
    gh_generate_preventive_timeline(&r, arena, &prs, 1, NULL, NULL, NULL, 0);
    const gh_screening *bp = find_screening(&r, "Blood pressure");
    check(bp != NULL, "blood pressure present");
    if (bp)
        check(bp->start_age == 18, "start age floored at 18");
}

static void test_preventive_prs_new_entry(gh_arena *arena)
{
    /* coronary_artery_disease has a modifier but no base screening, so it
     * becomes a new entry at 30 -- and its priority is hard-coded "elevated"
     * in the Python even at the "high" band. */
    gh_prs_result prs;
    memset(&prs, 0, sizeof(prs));
    prs.id = "coronary_artery_disease";
    prs.name = "Coronary Artery Disease";
    prs.risk_category = "high";
    prs.percentile = 99.0;

    gh_preventive_result r;
    gh_generate_preventive_timeline(&r, arena, &prs, 1, NULL, NULL, NULL, 0);
    check(r.ntimeline == GH_BASE_SCREENING_COUNT + 1, "one entry added");

    const gh_screening *cad =
        find_screening(&r, "Coronary Artery Disease screening");
    check(cad != NULL, "new screening named after the condition");
    if (cad) {
        check(cad->start_age == 30, "new entries start at 30");
        check_str(cad->priority, "elevated",
                  "new entries are 'elevated' even at the high band");
    }
}

static void test_preventive_average_prs_is_ignored(gh_arena *arena)
{
    gh_prs_result prs;
    memset(&prs, 0, sizeof(prs));
    prs.id = "type2_diabetes";
    prs.name = "Type 2 Diabetes";
    prs.risk_category = "average";
    prs.percentile = 50.0;

    gh_preventive_result r;
    gh_generate_preventive_timeline(&r, arena, &prs, 1, NULL, NULL, NULL, 0);
    check(r.early_screenings == 0, "average risk changes nothing");
    const gh_screening *glucose = find_screening(&r, "Fasting glucose / HbA1c");
    if (glucose)
        check(glucose->start_age == 35, "start age unchanged");
}

static void test_preventive_apoe(gh_arena *arena)
{
    gh_apoe_result apoe;
    memset(&apoe, 0, sizeof(apoe));
    apoe.apoe_type = "e4/e4";
    apoe.risk_level = "high";

    gh_preventive_result r;
    gh_generate_preventive_timeline(&r, arena, NULL, 0, &apoe, NULL, NULL, 0);
    check(r.ntimeline == GH_BASE_SCREENING_COUNT + GH_APOE_HIGH_SCREENINGS_COUNT,
          "APOE high adds its three screenings");
    check(r.early_screenings == GH_APOE_HIGH_SCREENINGS_COUNT,
          "each counts as an early screening");

    const gh_screening *cognitive = find_screening(&r, "Cognitive screening");
    check(cognitive != NULL, "cognitive screening added");
    if (cognitive) {
        check(cognitive->start_age == 50, "e4/e4 cognitive screening at 50");
        check_str(cognitive->priority, "high", "high priority");
        check_str(cognitive->genetic_basis, "APOE e4/e4", "basis names the type");
    }

    /* The elevated table is a different, shorter one. */
    gh_apoe_result elevated;
    memset(&elevated, 0, sizeof(elevated));
    elevated.apoe_type = "e3/e4";
    elevated.risk_level = "elevated";
    gh_preventive_result r2;
    gh_generate_preventive_timeline(&r2, arena, NULL, 0, &elevated, NULL,
                                    NULL, 0);
    check(r2.early_screenings == GH_APOE_ELEVATED_SCREENINGS_COUNT,
          "elevated adds its own count");
    const gh_screening *cog2 = find_screening(&r2, "Cognitive screening");
    if (cog2) {
        check(cog2->start_age == 55, "e3/e4 cognitive screening at 55");
        check_str(cog2->priority, "elevated", "elevated priority");
    }

    /* A risk level with no table adds nothing. */
    gh_apoe_result neutral;
    memset(&neutral, 0, sizeof(neutral));
    neutral.apoe_type = "e3/e3";
    neutral.risk_level = "average";
    gh_preventive_result r3;
    gh_generate_preventive_timeline(&r3, arena, NULL, 0, &neutral, NULL,
                                    NULL, 0);
    check(r3.ntimeline == GH_BASE_SCREENING_COUNT, "average APOE adds nothing");
}

static void test_preventive_acmg(gh_arena *arena)
{
    gh_cv_finding variants[3];
    memset(variants, 0, sizeof(variants));
    variants[0].gene = "BRCA1";
    variants[1].gene = "MSH2";
    variants[2].gene = "TTN";      /* actionable, but not one with a rule */

    gh_acmg_finding acmg_findings[3];
    memset(acmg_findings, 0, sizeof(acmg_findings));
    for (size_t i = 0; i < 3; i++) {
        acmg_findings[i].finding = &variants[i];
        acmg_findings[i].acmg_category = "pathogenic";
    }
    gh_acmg_result acmg;
    memset(&acmg, 0, sizeof(acmg));
    acmg.findings = acmg_findings;
    acmg.nfindings = 3;

    gh_preventive_result r;
    gh_generate_preventive_timeline(&r, arena, NULL, 0, NULL, &acmg, NULL, 0);
    check(r.early_screenings == 3, "BRCA adds two, Lynch one, TTN none");
    check(r.ntimeline == GH_BASE_SCREENING_COUNT + 3, "three entries added");

    const gh_screening *mri = find_screening(&r, "Breast MRI + mammography");
    check(mri != NULL, "BRCA1 adds breast MRI");
    if (mri) {
        check(mri->start_age == 25, "breast MRI at 25");
        check_str(mri->priority, "urgent", "ACMG entries are urgent");
        check_str(mri->genetic_basis, "BRCA1 pathogenic variant", "basis");
        check_has(mri->reason, "Pathogenic BRCA1 variant", "reason names gene");
    }

    const gh_screening *colonoscopy = find_screening(&r, "Colonoscopy");
    check(colonoscopy != NULL, "MSH2 adds a colonoscopy");
    if (colonoscopy) {
        check(colonoscopy->start_age == 20, "Lynch colonoscopy at 20");
        check_has(colonoscopy->reason, "Lynch syndrome (MSH2)", "Lynch reason");
    }

    /* Age dominates priority in the sort, so the age-18 blood pressure check
     * still leads and the urgent Lynch colonoscopy follows it at 20. */
    check_str(r.timeline[0].test, "Blood pressure", "age sorts before priority");
    check_str(r.timeline[1].test, "Colonoscopy", "Lynch colonoscopy at 20");
    check_str(r.timeline[2].test, "Breast MRI + mammography", "BRCA MRI at 25");
}

static void test_preventive_pgx_card(gh_arena *arena)
{
    gh_star_result stars[GH_STAR_GENE_COUNT];
    size_t n = normal_stars(stars);

    gh_preventive_result none;
    gh_generate_preventive_timeline(&none, arena, NULL, 0, NULL, NULL, stars, n);
    check(none.ntimeline == GH_BASE_SCREENING_COUNT,
          "all-normal metabolism adds no card entry");

    set_phenotype(stars, n, "CYP2C19", "poor");
    set_phenotype(stars, n, "TPMT", "intermediate");
    set_phenotype(stars, n, "CYP2D6", "Unknown");

    gh_preventive_result r;
    gh_generate_preventive_timeline(&r, arena, NULL, 0, NULL, NULL, stars, n);
    const gh_screening *card =
        find_screening(&r, "Pharmacogenomic card review (bring to every prescriber)");
    check(card != NULL, "non-normal metabolism adds the card entry");
    if (card) {
        check_str(card->priority, "ongoing", "card priority");
        check_str(card->reason,
                  "Non-standard metabolizer for: CYP2C19, TPMT",
                  "card names the non-normal genes, skipping Unknown");
        check_str(card->genetic_basis, "Pharmacogenomic profile", "card basis");
    }
    /* The card does not count as an early screening. */
    check(r.early_screenings == 0, "the card is not an early screening");
}

static void test_preventive_sort_is_stable(gh_arena *arena)
{
    /* Blood pressure (18, standard) and the pharmacogenomic card (18,
     * ongoing) collide on age; the card must come first, because "ongoing"
     * outranks "standard" -- and the base entries are appended last, so an
     * unstable sort would show up here too. */
    gh_star_result stars[GH_STAR_GENE_COUNT];
    size_t n = normal_stars(stars);
    set_phenotype(stars, n, "CYP2C19", "poor");

    gh_preventive_result r;
    gh_generate_preventive_timeline(&r, arena, NULL, 0, NULL, NULL, stars, n);
    check(r.ntimeline >= 2, "two entries at age 18");
    if (r.ntimeline >= 2) {
        check(r.timeline[0].start_age == 18 && r.timeline[1].start_age == 18,
              "both entries are at 18");
        check_str(r.timeline[0].priority, "ongoing", "ongoing before standard");
        check_str(r.timeline[1].priority, "standard", "standard after ongoing");
    }

    for (size_t i = 1; i < r.ntimeline; i++) {
        const gh_screening *a = &r.timeline[i - 1], *b = &r.timeline[i];
        check(a->start_age < b->start_age ||
              (a->start_age == b->start_age &&
               gh_screening_priority_rank(a->priority) <=
               gh_screening_priority_rank(b->priority)),
              "timeline ordered by (age, priority)");
    }
}

/* ------------------------------------------------------------------ */
/* Generated tables                                                    */
/* ------------------------------------------------------------------ */

static void test_tables(void)
{
    check(GH_DOSE_DRUG_COUNT > 0, "drugs generated");
    for (size_t i = 0; i < GH_DOSE_DRUG_COUNT; i++) {
        const gh_dose_drug *d = &GH_DOSE_DRUGS[i];
        check(d->ngenes > 0, "drug names at least one gene");
        check(d->nrules > 0, "drug has rules");
        check(d->source != NULL && *d->source, "drug cites a source");
        check(d->display != NULL && d->display[0] >= 'A' &&
              d->display[0] <= 'Z', "display name is title-cased");
        for (size_t r = 0; r < d->nrules; r++) {
            check(d->rules[r].nvalues > 0, "rule tests at least one value");
            check(d->rules[r].action != NULL, "rule has an action");
            check(d->rules[r].dose_guidance != NULL, "rule has dose guidance");
            check(d->rules[r].source == GH_DOSE_STAR ||
                  d->rules[r].source == GH_DOSE_FINDING, "rule source known");
        }
    }

    check(GH_POLY_RULE_COUNT > 0, "polypharmacy rules generated");
    for (size_t i = 0; i < GH_POLY_RULE_COUNT; i++) {
        const gh_poly_rule *rule = &GH_POLY_RULES[i];
        check(rule->ngenes > 0 && rule->ngenes <= GH_POLY_MAX_GENES,
              "rule gene count fits the matched array");
        check(rule->ndrugs > 0, "rule names affected drugs");
        check(poly_severity_known(rule->severity), "severity in vocabulary");
    }

    check(GH_BASE_SCREENING_COUNT > 0, "base screenings generated");
    for (size_t i = 0; i < GH_BASE_SCREENING_COUNT; i++) {
        check(GH_BASE_SCREENINGS[i].start_age >= 18, "screening starts at 18+");
        check(GH_BASE_SCREENINGS[i].condition != NULL, "screening has an id");
    }

    /* A modifier for a condition that has no base screening is the new-entry
     * branch; one that does is the in-place branch. Both must exist, or the
     * tests above are not covering what they claim. */
    size_t in_place = 0, new_entry = 0;
    for (size_t i = 0; i < GH_PRS_ELEVATED_MODIFIERS_COUNT; i++) {
        bool found = false;
        for (size_t b = 0; b < GH_BASE_SCREENING_COUNT; b++)
            if (strcmp(GH_BASE_SCREENINGS[b].condition,
                       GH_PRS_ELEVATED_MODIFIERS[i].condition) == 0)
                found = true;
        if (found)
            in_place++;
        else
            new_entry++;
    }
    check(in_place > 0, "some modifiers change a base screening");
    check(new_entry > 0, "some modifiers add a new screening");
}

int main(void)
{
    gh_arena *arena = gh_arena_new(1 << 16);

    test_critical_keywords();
    test_dosing_no_data(arena);
    test_dosing_all_normal(arena);
    test_dosing_poor_metabolizer(arena);
    test_dosing_reads_lifestyle_findings(arena);
    test_dosing_requires_gene_data(arena);

    test_polypharmacy_single_gene(arena);
    test_polypharmacy_needs_every_gene(arena);
    test_polypharmacy_unknown_is_not_a_phenotype(arena);
    test_polypharmacy_finding_overrides_star(arena);
    test_polypharmacy_severity_order(arena);

    test_priority_rank();
    test_preventive_baseline(arena);
    test_preventive_prs_modifies_base(arena);
    test_preventive_start_age_floor(arena);
    test_preventive_prs_new_entry(arena);
    test_preventive_average_prs_is_ignored(arena);
    test_preventive_apoe(arena);
    test_preventive_acmg(arena);
    test_preventive_pgx_card(arena);
    test_preventive_sort_is_stable(arena);

    test_tables();

    gh_arena_free(arena);
    printf("%s: %d checks, %d failures\n",
           failures ? "FAILED" : "ok", checks, failures);
    return failures ? 1 : 0;
}

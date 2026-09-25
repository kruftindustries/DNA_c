#include "gh_profiles.h"

#include <stdio.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Shared helpers                                                      */
/* ------------------------------------------------------------------ */

static size_t count_allele(const char *genotype, char allele)
{
    size_t n = 0;
    for (const char *p = genotype; *p; p++)
        if (*p == allele)
            n++;
    return n;
}

/* Run every marker of a module against the genome. A marker whose SNP is
 * absent stays `present = false`, matching the assessors returning None. */
static void assess(gh_profile_base *base, const gh_genome *g,
                   const gh_profile_markers *markers)
{
    memset(base, 0, sizeof(*base));
    base->markers_tested = markers->count;

    for (size_t i = 0; i < markers->count; i++) {
        const gh_profile_marker *m = &markers->markers[i];
        gh_profile_hit *hit = &base->hits[i];
        hit->rsid = m->rsid;
        hit->gene = m->gene;

        const char *genotype = gh_genome_genotype(g, m->rsid);
        if (!genotype || !*genotype)
            continue;

        size_t copies = count_allele(genotype, m->allele);
        hit->present = true;
        hit->genotype = genotype;
        /* outcomes are ordered two copies, one, none. */
        hit->outcome = &m->outcomes[copies >= 2 ? 0 : copies == 1 ? 1 : 2];
        base->nhits++;
    }
}

/* The hit for a marker by index, or NULL when that SNP was not genotyped. */
static const gh_profile_hit *hit_at(const gh_profile_base *base, size_t index)
{
    const gh_profile_hit *hit = &base->hits[index];
    return hit->present ? hit : NULL;
}

static void recommend(gh_profile_base *base, const char *text)
{
    if (base->nrecommendations < GH_PROFILE_MAX_RECOMMENDATIONS)
        base->recommendations[base->nrecommendations++] = text;
}

/* Join clause fragments as "<prefix><a>; <b>.", the shape several of these
 * summaries use. */
static const char *join_clauses(gh_arena *arena, const char *prefix,
                                const char *const *parts, size_t n)
{
    size_t len = strlen(prefix) + 2;
    for (size_t i = 0; i < n; i++)
        len += strlen(parts[i]) + 2;

    char *out = gh_alloc(arena, len + 1);
    size_t at = (size_t)snprintf(out, len + 1, "%s", prefix);
    for (size_t i = 0; i < n; i++)
        at += (size_t)snprintf(out + at, len + 1 - at, "%s%s",
                               i ? "; " : "", parts[i]);
    snprintf(out + at, len + 1 - at, ".");
    return out;
}

const char *gh_profile_domain_level(const gh_profile_domain *domains,
                                    size_t n, const char *name)
{
    for (size_t i = 0; i < n; i++)
        if (strcmp(domains[i].name, name) == 0)
            return domains[i].level;
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Histamine                                                           */
/* ------------------------------------------------------------------ */

static const char *const HISTAMINE_FOODS[] = {
    "Aged cheeses (parmesan, cheddar, gouda)",
    "Fermented foods (sauerkraut, kimchi, kombucha)",
    "Cured meats (salami, prosciutto, bacon)",
    "Alcoholic beverages (especially red wine and beer)",
    "Smoked fish and canned tuna",
    "Tomatoes, spinach, and eggplant",
    "Vinegar and soy sauce",
    "Leftover/reheated meals (histamine increases with time)",
};

void gh_profile_histamine(gh_histamine_result *out, gh_arena *arena,
                          const gh_genome *g)
{
    memset(out, 0, sizeof(*out));
    assess(&out->base, g, &GH_HISTAMINE_MARKERS);

    out->foods_to_watch = HISTAMINE_FOODS;
    out->nfoods = sizeof(HISTAMINE_FOODS) / sizeof(HISTAMINE_FOODS[0]);

    /* Risk points accumulate across the DAO and HNMT markers. */
    for (size_t i = 0; i < out->base.markers_tested; i++) {
        const gh_profile_hit *hit = hit_at(&out->base, i);
        if (hit)
            out->total_risk += hit->outcome->points;
    }

    if (out->base.nhits == 0)
        out->risk_level = "unknown";
    else if (out->total_risk >= 4)
        out->risk_level = "elevated";
    else if (out->total_risk >= 2)
        out->risk_level = "moderate";
    else
        out->risk_level = "low";

    if (out->base.nhits == 0)
        out->base.summary =
            "No histamine metabolism markers available in genotype data.";
    else if (strcmp(out->risk_level, "elevated") == 0)
        out->base.summary =
            "Your genetic profile suggests elevated histamine intolerance risk. "
            "Both DAO and/or HNMT pathways show reduced activity, which may lead to "
            "symptoms like headaches, flushing, nasal congestion, or GI issues after "
            "histamine-rich foods.";
    else if (strcmp(out->risk_level, "moderate") == 0)
        out->base.summary =
            "Your genetic profile suggests moderate histamine intolerance risk. "
            "Some reduction in histamine-degrading enzyme activity detected. "
            "You may notice symptoms with large amounts of histamine-rich foods.";
    else
        out->base.summary =
            "Your genetic profile suggests low histamine intolerance risk. "
            "Histamine-degrading enzyme activity appears normal based on tested markers.";

    if (strcmp(out->risk_level, "elevated") == 0) {
        recommend(&out->base, "Consider a low-histamine diet trial (2-4 weeks) to assess symptom improvement.");
        recommend(&out->base, "DAO supplementation before meals may help with food-triggered symptoms.");
        recommend(&out->base, "Track symptoms with a food diary to identify personal triggers.");
        recommend(&out->base, "Discuss with your physician -- symptoms overlap with allergies and mast cell disorders.");
    } else if (strcmp(out->risk_level, "moderate") == 0) {
        recommend(&out->base, "Be aware of cumulative histamine load -- spacing high-histamine meals may help.");
        recommend(&out->base, "Freshly prepared foods have lower histamine than leftovers or fermented foods.");
        recommend(&out->base, "Consider DAO supplementation if you notice post-meal flushing or headaches.");
    } else {
        /* Covers both "low" and the no-data case, as in the Python. */
        recommend(&out->base, "No specific dietary modifications indicated based on histamine genetics alone.");
    }
}

/* ------------------------------------------------------------------ */
/* Alcohol                                                             */
/* ------------------------------------------------------------------ */

enum { ALC_ADH1B = 0, ALC_ALDH2 = 1, ALC_CYP2E1 = 2 };

void gh_profile_alcohol(gh_alcohol_result *out, gh_arena *arena,
                        const gh_genome *g)
{
    memset(out, 0, sizeof(*out));
    assess(&out->base, g, &GH_ALCOHOL_MARKERS);

    const gh_profile_hit *adh1b = hit_at(&out->base, ALC_ADH1B);
    const gh_profile_hit *aldh2 = hit_at(&out->base, ALC_ALDH2);
    const gh_profile_hit *cyp2e1 = hit_at(&out->base, ALC_CYP2E1);

    /* ADH1B sets the baseline speed; both "fast" and "very_fast" report as
     * fast, and CYP2E1 can only slow a normal result. */
    out->metabolism_speed = "normal";
    if (adh1b) {
        const char *speed = adh1b->outcome->category_a;
        if (speed && (strcmp(speed, "very_fast") == 0 ||
                      strcmp(speed, "fast") == 0))
            out->metabolism_speed = "fast";
    }
    if (cyp2e1 && cyp2e1->outcome->category_a &&
        strcmp(cyp2e1->outcome->category_a, "reduced") == 0 &&
        strcmp(out->metabolism_speed, "normal") == 0)
        out->metabolism_speed = "slow";

    out->flush_risk = aldh2 && aldh2->outcome->category_b
        ? aldh2->outcome->category_b : "unknown";

    int cancer_score = 0;
    for (size_t i = 0; i < out->base.markers_tested; i++) {
        const gh_profile_hit *hit = hit_at(&out->base, i);
        if (hit)
            cancer_score += hit->outcome->points;
    }
    out->cancer_risk = cancer_score >= 3 ? "high"
                     : cancer_score >= 1 ? "elevated" : "average";

    if (out->base.nhits == 0) {
        out->base.summary =
            "No alcohol metabolism markers available in genotype data.";
    } else {
        const char *parts[3];
        size_t nparts = 0;
        if (strcmp(out->metabolism_speed, "fast") == 0)
            parts[nparts++] = "fast ethanol-to-acetaldehyde conversion (ADH1B)";
        else if (strcmp(out->metabolism_speed, "slow") == 0)
            parts[nparts++] = "slower ethanol metabolism";
        if (strcmp(out->flush_risk, "severe") == 0)
            parts[nparts++] = "severe alcohol flush reaction (ALDH2 deficiency)";
        else if (strcmp(out->flush_risk, "mild") == 0)
            parts[nparts++] = "mild alcohol flush risk (ALDH2 partial deficiency)";
        if (strcmp(out->cancer_risk, "high") == 0)
            parts[nparts++] = "high alcohol-related cancer risk";
        else if (strcmp(out->cancer_risk, "elevated") == 0)
            parts[nparts++] = "elevated alcohol-related cancer risk";

        out->base.summary = nparts
            ? join_clauses(arena, "Your alcohol metabolism profile shows: ",
                           parts, nparts)
            : "Your alcohol metabolism profile is typical. No notable genetic "
              "variants affecting alcohol processing detected.";
    }

    if (strcmp(out->flush_risk, "severe") == 0) {
        recommend(&out->base, "ALDH2 deficiency detected -- strongly consider avoiding alcohol entirely.");
        recommend(&out->base, "Even small amounts of alcohol cause toxic acetaldehyde accumulation, increasing esophageal cancer risk up to 10x.");
        recommend(&out->base, "Do not take medications that inhibit ALDH (e.g., disulfiram) without medical supervision.");
    } else if (strcmp(out->flush_risk, "mild") == 0) {
        recommend(&out->base, "Partial ALDH2 deficiency detected -- limit alcohol consumption significantly.");
        recommend(&out->base, "Drinking despite flush reaction substantially increases esophageal and head/neck cancer risk.");
        recommend(&out->base, "If you choose to drink, limit to very small quantities and never on an empty stomach.");
    }
    if (strcmp(out->metabolism_speed, "fast") == 0 &&
        (strcmp(out->flush_risk, "none") == 0 ||
         strcmp(out->flush_risk, "unknown") == 0))
        recommend(&out->base, "Fast ADH1B metabolism provides some natural protection against alcohol dependence, but does not eliminate health risks from heavy drinking.");
    if (strcmp(out->cancer_risk, "elevated") == 0 ||
        strcmp(out->cancer_risk, "high") == 0)
        recommend(&out->base, "Discuss alcohol-related cancer screening (esophageal, head/neck) with your physician, especially if you have a history of alcohol consumption.");
    if (out->base.nrecommendations == 0 && out->base.nhits > 0)
        recommend(&out->base, "Standard alcohol guidelines apply. No genetic variants requiring special precautions detected.");
}

/* ------------------------------------------------------------------ */
/* Pain sensitivity                                                    */
/* ------------------------------------------------------------------ */

enum { PAIN_OPRM1 = 0, PAIN_COMT = 1, PAIN_SCN9A = 2, PAIN_TRPV1 = 3 };

void gh_profile_pain(gh_pain_result *out, gh_arena *arena, const gh_genome *g)
{
    memset(out, 0, sizeof(*out));
    assess(&out->base, g, &GH_PAIN_MARKERS);

    /* Mean of the per-marker contributions, or a neutral 50 with no data. */
    int total = 0;
    size_t n = 0;
    for (size_t i = 0; i < out->base.markers_tested; i++) {
        const gh_profile_hit *hit = hit_at(&out->base, i);
        if (hit) {
            total += hit->outcome->points;
            n++;
        }
    }
    if (n) {
        /* Python rounds the mean half-to-even. */
        double mean = (double)total / (double)n;
        double floored = (double)(long)mean;
        double frac = mean - floored;
        long rounded;
        if (frac > 0.5)
            rounded = (long)floored + 1;
        else if (frac < 0.5)
            rounded = (long)floored;
        else
            rounded = ((long)floored % 2 == 0) ? (long)floored
                                               : (long)floored + 1;
        out->sensitivity_score = (int)rounded;
    } else {
        out->sensitivity_score = 50;
    }

    if (out->base.nhits == 0)
        out->base.summary =
            "No pain sensitivity markers available in genotype data.";
    else if (out->sensitivity_score >= 70)
        out->base.summary =
            "Your genetic profile suggests higher-than-average pain sensitivity. "
            "You may benefit from proactive pain management strategies and should "
            "discuss opioid dosing with your physician.";
    else if (out->sensitivity_score >= 50)
        out->base.summary =
            "Your genetic profile suggests average pain sensitivity across the "
            "markers tested.";
    else if (out->sensitivity_score >= 30)
        out->base.summary =
            "Your genetic profile suggests lower-than-average pain sensitivity. "
            "You may have higher tolerance to pain stimuli but could need adjusted "
            "opioid dosing.";
    else
        out->base.summary =
            "Your genetic profile suggests notably low pain sensitivity. Opioid "
            "medications may be less effective; discuss alternative analgesics with "
            "your physician.";

    if (out->base.nhits == 0)
        return;

    const gh_profile_hit *oprm1 = hit_at(&out->base, PAIN_OPRM1);
    const gh_profile_hit *comt = hit_at(&out->base, PAIN_COMT);
    const gh_profile_hit *trpv1 = hit_at(&out->base, PAIN_TRPV1);

    if (oprm1 && oprm1->outcome->points <= 35)
        recommend(&out->base, "Inform your anesthesiologist about OPRM1 AG/GG status -- you may require higher opioid doses or alternative analgesics.");
    if (comt && comt->outcome->points >= 70)
        recommend(&out->base, "COMT Met/Met carriers may benefit from stress-reduction techniques (meditation, yoga) to manage heightened pain perception.");
    if (trpv1 && trpv1->outcome->points >= 70)
        recommend(&out->base, "Higher capsaicin sensitivity detected -- consider gradual spice exposure and topical capsaicin patches at lower concentrations.");
    recommend(&out->base, "Share your pain sensitivity profile with your healthcare provider for personalized pain management planning.");
}

/* ------------------------------------------------------------------ */
/* Thyroid                                                             */
/* ------------------------------------------------------------------ */

static const char *const THYROID_DOMAINS[] = {
    "autoimmune_risk", "conversion_efficiency", "cancer_risk",
};

/* Average severity decides the level; an empty domain is unknown. */
static const char *classify_mean(const double *scores, size_t n)
{
    if (!n)
        return "unknown";
    double total = 0.0;
    for (size_t i = 0; i < n; i++)
        total += scores[i];
    double avg = total / (double)n;
    if (avg >= 1.5)
        return "elevated";
    if (avg >= 0.8)
        return "moderate";
    return "average";
}

void gh_profile_thyroid(gh_thyroid_result *out, gh_arena *arena,
                        const gh_genome *g)
{
    memset(out, 0, sizeof(*out));
    assess(&out->base, g, &GH_THYROID_MARKERS);

    out->ndomains = sizeof(THYROID_DOMAINS) / sizeof(THYROID_DOMAINS[0]);
    for (size_t d = 0; d < out->ndomains; d++) {
        double scores[GH_PROFILE_MAX_MARKERS];
        size_t n = 0;
        for (size_t i = 0; i < out->base.markers_tested; i++) {
            const gh_profile_hit *hit = hit_at(&out->base, i);
            if (hit && hit->outcome->category_a &&
                strcmp(hit->outcome->category_a, THYROID_DOMAINS[d]) == 0)
                scores[n++] = hit->outcome->severity;
        }
        out->domains[d].name = THYROID_DOMAINS[d];
        out->domains[d].level = classify_mean(scores, n);
    }

    const char *ai = out->domains[0].level;
    const char *conv = out->domains[1].level;
    const char *cancer = out->domains[2].level;

    if (out->base.nhits == 0) {
        out->base.summary =
            "No thyroid-related markers available in genotype data.";
    } else {
        const char *parts[3];
        size_t nparts = 0;
        if (strcmp(ai, "elevated") == 0)
            parts[nparts++] = "elevated autoimmune thyroid risk (Hashimoto's/Graves')";
        else if (strcmp(ai, "moderate") == 0)
            parts[nparts++] = "moderate autoimmune thyroid risk";
        if (strcmp(conv, "elevated") == 0)
            parts[nparts++] = "impaired T4-to-T3 conversion in both DIO1 and DIO2";
        else if (strcmp(conv, "moderate") == 0)
            parts[nparts++] = "mildly reduced T4-to-T3 conversion";
        if (strcmp(cancer, "elevated") == 0)
            parts[nparts++] = "elevated thyroid cancer susceptibility";
        else if (strcmp(cancer, "moderate") == 0)
            parts[nparts++] = "moderate thyroid cancer susceptibility";

        out->base.summary = nparts
            ? join_clauses(arena, "Your thyroid genetic profile shows: ",
                           parts, nparts)
            : "Your thyroid genetic profile appears average across all tested markers.";
    }

    bool ai_flag = strcmp(ai, "elevated") == 0 || strcmp(ai, "moderate") == 0;
    bool conv_flag = strcmp(conv, "elevated") == 0 || strcmp(conv, "moderate") == 0;
    bool cancer_flag = strcmp(cancer, "elevated") == 0 ||
                       strcmp(cancer, "moderate") == 0;

    if (ai_flag) {
        recommend(&out->base, "Monitor thyroid function annually (TSH, free T4, anti-TPO antibodies).");
        recommend(&out->base, "Ensure adequate selenium intake (200 mcg/day) -- supports both deiodinase function and may reduce anti-TPO antibodies.");
        recommend(&out->base, "Maintain adequate iodine intake but avoid excess supplementation, which can trigger autoimmune thyroid disease.");
    }
    if (conv_flag) {
        recommend(&out->base, "If hypothyroid, discuss combination T4/T3 therapy with your endocrinologist -- DIO variants may impair T4-only response.");
        recommend(&out->base, "Selenium and zinc support deiodinase enzyme function.");
        recommend(&out->base, "Monitor both free T4 and free T3 levels, not just TSH.");
    }
    if (cancer_flag)
        recommend(&out->base, "Discuss thyroid ultrasound screening with your physician, especially if there is a family history of thyroid cancer.");
    if (out->base.nrecommendations == 0 && out->base.nhits > 0)
        recommend(&out->base, "No specific thyroid interventions indicated based on genetic markers tested.");
}

/* ------------------------------------------------------------------ */
/* Eye health                                                          */
/* ------------------------------------------------------------------ */

static const char *const EYE_CONDITIONS[] = {"glaucoma_risk", "myopia_risk"};

/* A single high-severity marker outranks the average, so a high-penetrance
 * variant is not diluted by benign ones. */
static const char *classify_eye(const double *scores, size_t n)
{
    if (!n)
        return "unknown";
    double total = 0.0, max = scores[0];
    for (size_t i = 0; i < n; i++) {
        total += scores[i];
        if (scores[i] > max)
            max = scores[i];
    }
    if (max >= 3.0)
        return "high";
    double avg = total / (double)n;
    if (avg >= 1.5)
        return "elevated";
    if (avg >= 0.8)
        return "moderate";
    return "average";
}

void gh_profile_eye(gh_eye_result *out, gh_arena *arena, const gh_genome *g)
{
    memset(out, 0, sizeof(*out));
    assess(&out->base, g, &GH_EYE_MARKERS);

    out->nconditions = sizeof(EYE_CONDITIONS) / sizeof(EYE_CONDITIONS[0]);
    for (size_t d = 0; d < out->nconditions; d++) {
        double scores[GH_PROFILE_MAX_MARKERS];
        size_t n = 0;
        for (size_t i = 0; i < out->base.markers_tested; i++) {
            const gh_profile_hit *hit = hit_at(&out->base, i);
            if (hit && hit->outcome->category_a &&
                strcmp(hit->outcome->category_a, EYE_CONDITIONS[d]) == 0)
                scores[n++] = hit->outcome->severity;
        }
        out->conditions[d].name = EYE_CONDITIONS[d];
        out->conditions[d].level = classify_eye(scores, n);
    }

    const char *glaucoma = out->conditions[0].level;
    const char *myopia = out->conditions[1].level;

    if (out->base.nhits == 0) {
        out->base.summary = "No eye health markers available in genotype data.";
    } else {
        const char *parts[2];
        size_t nparts = 0;
        if (strcmp(glaucoma, "high") == 0)
            parts[nparts++] = "high-penetrance glaucoma variant detected (MYOC) -- urgent ophthalmology referral recommended";
        else if (strcmp(glaucoma, "elevated") == 0)
            parts[nparts++] = "elevated glaucoma risk";
        else if (strcmp(glaucoma, "moderate") == 0)
            parts[nparts++] = "moderate glaucoma risk";
        if (strcmp(myopia, "elevated") == 0 || strcmp(myopia, "high") == 0)
            parts[nparts++] = "elevated myopia susceptibility";
        else if (strcmp(myopia, "moderate") == 0)
            parts[nparts++] = "moderate myopia susceptibility";

        out->base.summary = nparts
            ? join_clauses(arena, "Your eye health profile shows: ", parts, nparts)
            : "Your eye health genetic profile appears average across all tested markers.";
    }

    if (strcmp(glaucoma, "high") == 0) {
        recommend(&out->base, "Urgent: MYOC glaucoma variant detected. Schedule an ophthalmology evaluation for intraocular pressure and optic nerve assessment.");
        recommend(&out->base, "Annual comprehensive eye exams with IOP measurement and visual field testing.");
        recommend(&out->base, "First-degree relatives should also be screened for this variant.");
    } else if (strcmp(glaucoma, "elevated") == 0 ||
               strcmp(glaucoma, "moderate") == 0) {
        recommend(&out->base, "Regular comprehensive eye exams including intraocular pressure measurement, especially after age 40.");
        recommend(&out->base, "Report any vision changes (halos, peripheral vision loss) to your ophthalmologist promptly.");
    }
    if (strcmp(myopia, "elevated") == 0 || strcmp(myopia, "high") == 0) {
        recommend(&out->base, "For children: encourage outdoor time (at least 2 hours/day) -- the strongest modifiable factor for myopia prevention.");
        recommend(&out->base, "Follow the 20-20-20 rule: every 20 minutes of near work, look at something 20 feet away for 20 seconds.");
        recommend(&out->base, "Regular eye exams to monitor refractive changes and screen for myopia complications (retinal detachment, macular degeneration).");
    } else if (strcmp(myopia, "moderate") == 0) {
        recommend(&out->base, "Moderate myopia genetic risk -- regular vision screening and outdoor time recommended.");
    }
    if (out->base.nrecommendations == 0 && out->base.nhits > 0)
        recommend(&out->base, "No specific eye health interventions indicated based on tested genetic markers. Routine eye exams still recommended.");
}

/* ------------------------------------------------------------------ */
/* Hormone metabolism                                                  */
/* ------------------------------------------------------------------ */

void gh_profile_hormone(gh_hormone_result *out, gh_arena *arena,
                        const gh_genome *g)
{
    memset(out, 0, sizeof(*out));
    assess(&out->base, g, &GH_HORMONE_MARKERS);

    /* category_a names the pathway a marker belongs to; category_b and
     * category_c carry its estrogen and androgen effects. A marker only
     * fills the one that applies to it. */
    bool e_increased = false, e_enhanced = false;
    bool e_decreased = false, e_reduced = false, e_any = false;
    bool a_increased = false, a_decreased = false, a_any = false;

    for (size_t i = 0; i < out->base.markers_tested; i++) {
        const gh_profile_hit *hit = hit_at(&out->base, i);
        if (!hit)
            continue;
        const char *domain = hit->outcome->category_a;
        const char *estrogen = hit->outcome->category_b;
        const char *androgen = hit->outcome->category_c;

        /* Python partitions the markers by domain before summarising each
         * side, so an effect only counts toward its own pathway. */
        bool is_estrogen = domain &&
                           strcmp(domain, "estrogen_metabolism") == 0;
        bool is_androgen = domain &&
                           strcmp(domain, "androgen_metabolism") == 0;

        if (is_estrogen && estrogen) {
            e_any = true;
            if (strcmp(estrogen, "increased") == 0) e_increased = true;
            if (strcmp(estrogen, "enhanced_signaling") == 0) e_enhanced = true;
            if (strcmp(estrogen, "decreased") == 0) e_decreased = true;
            if (strcmp(estrogen, "reduced_signaling") == 0) e_reduced = true;
        }
        if (is_androgen && androgen) {
            a_any = true;
            if (strcmp(androgen, "increased_dht") == 0) a_increased = true;
            if (strcmp(androgen, "decreased_dht") == 0) a_decreased = true;
        }
    }

    if (!e_any)
        out->estrogen_level = "unknown";
    else if (e_increased && e_enhanced)
        out->estrogen_level = "high estrogen activity";
    else if (e_increased || e_enhanced)
        out->estrogen_level = "moderately increased estrogen activity";
    else if (e_decreased || e_reduced)
        out->estrogen_level = "reduced estrogen activity";
    else
        out->estrogen_level = "normal estrogen activity";

    if (!a_any)
        out->androgen_level = "unknown";
    else if (a_increased)
        out->androgen_level = "elevated DHT pathway";
    else if (a_decreased)
        out->androgen_level = "reduced DHT pathway";
    else
        out->androgen_level = "normal androgen metabolism";

    if (out->base.nhits == 0) {
        out->overall = "unknown";
        out->base.summary =
            "No hormone metabolism markers available in genotype data.";
        return;
    }

    const char *parts[2];
    size_t nparts = 0;
    if (strstr(out->estrogen_level, "increased") ||
        strstr(out->estrogen_level, "high"))
        parts[nparts++] = "higher estrogen activity";
    else if (strstr(out->estrogen_level, "reduced"))
        parts[nparts++] = "lower estrogen activity";
    if (strstr(out->androgen_level, "elevated"))
        parts[nparts++] = "elevated DHT conversion";
    else if (strstr(out->androgen_level, "reduced"))
        parts[nparts++] = "reduced DHT conversion";

    if (nparts) {
        size_t len = 1;
        for (size_t i = 0; i < nparts; i++)
            len += strlen(parts[i]) + 2;
        char *joined = gh_alloc(arena, len + 1);
        size_t at = 0;
        for (size_t i = 0; i < nparts; i++)
            at += (size_t)snprintf(joined + at, len + 1 - at, "%s%s",
                                   i ? ", " : "", parts[i]);
        out->overall = joined;
    } else {
        out->overall = "balanced hormone metabolism";
    }

    size_t slen = strlen(out->overall) + 128;
    char *summary = gh_alloc(arena, slen);
    if (strcmp(out->overall, "balanced hormone metabolism") == 0)
        snprintf(summary, slen,
                 "Your hormone metabolism profile indicates %s."
                 " No notable deviations detected in the tested markers.",
                 out->overall);
    else
        snprintf(summary, slen,
                 "Your hormone metabolism profile indicates %s.", out->overall);
    out->base.summary = summary;

    if (strstr(out->estrogen_level, "increased") ||
        strstr(out->estrogen_level, "high")) {
        recommend(&out->base, "Higher aromatase activity may increase estrogen levels. Cruciferous vegetables (broccoli, cauliflower) support healthy estrogen metabolism via DIM/I3C.");
        recommend(&out->base, "Maintain a healthy body fat percentage -- adipose tissue is a significant source of aromatase activity.");
        recommend(&out->base, "Discuss estrogen-related cancer screening frequency with your physician.");
    }
    if (strstr(out->estrogen_level, "reduced")) {
        recommend(&out->base, "Lower estrogen activity may affect bone density. Ensure adequate calcium, vitamin D, and weight-bearing exercise.");
        recommend(&out->base, "Consider bone density screening, especially post-menopause or if other osteoporosis risk factors are present.");
    }
    if (strstr(out->androgen_level, "elevated")) {
        recommend(&out->base, "Higher DHT conversion may increase androgenetic alopecia and prostate growth risk.");
        recommend(&out->base, "Saw palmetto and green tea (EGCG) are natural 5-alpha reductase inhibitors -- discuss with your physician.");
    }
    if (out->base.nrecommendations == 0 && out->base.nhits > 0)
        recommend(&out->base, "No specific hormone interventions indicated based on tested genetic markers.");
}

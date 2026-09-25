#include "gh_recommendations.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gh_report.h"   /* GH_PATHWAYS, for the insight pathway list */

/* ------------------------------------------------------------------ */
/* Helpers                                                             */
/* ------------------------------------------------------------------ */

static char *lower_copy(gh_arena *arena, const char *s)
{
    size_t n = s ? strlen(s) : 0;
    char *out = gh_alloc(arena, n + 1);
    for (size_t i = 0; i < n; i++)
        out[i] = (char)tolower((unsigned char)s[i]);
    out[n] = '\0';
    return out;
}

static char *upper_copy(gh_arena *arena, const char *s)
{
    size_t n = s ? strlen(s) : 0;
    char *out = gh_alloc(arena, n + 1);
    for (size_t i = 0; i < n; i++)
        out[i] = (char)toupper((unsigned char)s[i]);
    out[n] = '\0';
    return out;
}

static bool list_has(const char *const *list, size_t n, const char *value)
{
    for (size_t i = 0; i < n; i++)
        if (strcmp(list[i], value) == 0)
            return true;
    return false;
}

/* The highest-magnitude finding for a gene, or NULL when absent. Mirrors the
 * Python's gene_findings lookup, which keeps the strongest per gene. */
static const gh_finding *gene_finding(const gh_analysis *a, const char *gene)
{
    const gh_finding *best = NULL;
    if (!a)
        return NULL;
    for (size_t i = 0; i < a->nfindings; i++) {
        const gh_finding *f = &a->findings[i];
        if (!f->gene || strcmp(f->gene, gene) != 0)
            continue;
        if (!best || f->magnitude > best->magnitude)
            best = f;
    }
    return best;
}

const char *gh_compute_priority(size_t gene_signals, bool pathogenic,
                                bool prs_high, bool prs_elevated,
                                bool high_magnitude, bool epistasis)
{
    if (pathogenic || prs_high || gene_signals >= 3)
        return "high";
    if (prs_elevated || gene_signals >= 2 || high_magnitude || epistasis)
        return "moderate";
    return "low";
}

/* Join detail strings with "; ", falling back to the Python's placeholder. */
static const char *build_why(gh_arena *arena, const char *const *details,
                             size_t n)
{
    if (!n)
        return "Single genetic signal detected.";
    size_t len = 1;
    for (size_t i = 0; i < n; i++)
        len += strlen(details[i]) + 2;
    char *out = gh_alloc(arena, len);
    size_t at = 0;
    for (size_t i = 0; i < n; i++)
        at += (size_t)snprintf(out + at, len - at, "%s%s",
                               i ? "; " : "", details[i]);
    return out;
}

/* ------------------------------------------------------------------ */
/* Priorities                                                          */
/* ------------------------------------------------------------------ */

static int priority_rank(const char *priority)
{
    if (strcmp(priority, "high") == 0)
        return 0;
    if (strcmp(priority, "moderate") == 0)
        return 1;
    if (strcmp(priority, "low") == 0)
        return 2;
    return 3;
}

typedef struct {
    gh_priority item;
    size_t order;
} ranked_priority;

/* Highest priority first, then the most corroborated. Python sorts on
 * (rank, -signal_count) with a stable sort, so ties keep group order. */
static int compare_priority(const void *a, const void *b)
{
    const ranked_priority *x = a, *y = b;
    int rx = priority_rank(x->item.priority);
    int ry = priority_rank(y->item.priority);
    if (rx != ry)
        return rx < ry ? -1 : 1;
    if (x->item.signal_count != y->item.signal_count)
        return x->item.signal_count > y->item.signal_count ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

/* The first three pipe-separated traits. ClinVar appends loosely associated
 * somatic conditions later in the field, which would otherwise match a
 * group's keywords spuriously. */
static const char *primary_traits(gh_arena *arena, const char *traits)
{
    const char *lowered = lower_copy(arena, traits ? traits : "");
    const char *p = lowered;
    int bars = 0;
    while (*p && bars < 3) {
        if (*p == '|')
            bars++;
        if (bars == 3)
            break;
        p++;
    }
    size_t len = (size_t)(p - lowered);
    /* Python strips the joined prefix. */
    while (len && isspace((unsigned char)lowered[len - 1]))
        len--;
    const char *start = lowered;
    while (len && isspace((unsigned char)*start)) {
        start++;
        len--;
    }
    return gh_strndup(arena, start, len);
}

static bool any_keyword(const char *haystack, const char *const *keywords,
                        size_t n)
{
    for (size_t i = 0; i < n; i++)
        if (strstr(haystack, keywords[i]))
            return true;
    return false;
}

static void build_priorities(gh_recommendations *out, gh_arena *arena,
                             const gh_recommendation_inputs *in)
{
    ranked_priority *ranked =
        gh_calloc(arena, GH_RISK_GROUP_COUNT, sizeof(*ranked));
    size_t n = 0;

    for (size_t g = 0; g < GH_RISK_GROUP_COUNT; g++) {
        const gh_risk_group *group = &GH_RISK_GROUPS[g];

        size_t signal_count = 0, gene_signals = 0;
        bool pathogenic = false, prs_high = false, prs_elevated = false;
        bool high_magnitude = false, epistasis = false;

        const char *details[32];
        size_t ndetails = 0;

        /* 1. Gene variants. Normal or low-magnitude calls are not signals. */
        const char *matched[16];
        size_t nmatched = 0;
        for (size_t i = 0; i < group->ngenes; i++) {
            const gh_finding *f = gene_finding(in->analysis, group->genes[i]);
            if (!f)
                continue;
            const char *status = f->status ? f->status : "";
            if (!*status || strcmp(status, "normal") == 0 ||
                strcmp(status, "reference") == 0 ||
                strcmp(status, "typical") == 0 || f->magnitude < 2)
                continue;

            signal_count++;
            gene_signals++;
            if (f->magnitude >= 3) {
                signal_count++;
                high_magnitude = true;
            }
            if (nmatched < 16) {
                size_t len = strlen(group->genes[i]) + strlen(status) + 8;
                char *entry = gh_alloc(arena, len);
                snprintf(entry, len, "%s (%s)", group->genes[i], status);
                matched[nmatched++] = entry;
            }
        }
        if (nmatched) {
            size_t len = 32;
            for (size_t i = 0; i < nmatched; i++)
                len += strlen(matched[i]) + 2;
            char *line = gh_alloc(arena, len);
            size_t at = (size_t)snprintf(line, len, "Gene variants: ");
            for (size_t i = 0; i < nmatched; i++)
                at += (size_t)snprintf(line + at, len - at, "%s%s",
                                       i ? ", " : "", matched[i]);
            details[ndetails++] = line;
        }

        /* 2. Polygenic score for the group's condition. */
        if (in->prs && group->prs_condition && *group->prs_condition) {
            for (size_t i = 0; i < in->nprs; i++) {
                const gh_prs_result *r = &in->prs[i];
                if (strcmp(r->id, group->prs_condition) != 0)
                    continue;
                bool high = strcmp(r->risk_category, "high") == 0;
                bool elevated = strcmp(r->risk_category, "elevated") == 0;
                if (!high && !elevated)
                    continue;
                signal_count++;
                if (high)
                    prs_high = true;
                else
                    prs_elevated = true;
                if (ndetails < 32) {
                    size_t len = strlen(r->name) + 64;
                    char *line = gh_alloc(arena, len);
                    snprintf(line, len, "PRS %s: %.0fth percentile (%s)",
                             r->name, r->percentile,
                             high ? "high" : "elevated");
                    details[ndetails++] = line;
                }
            }
        }

        /* 3. ClinVar findings, by gene or by keyword on the primary traits. */
        if (in->clinvar && in->clinvar->loaded) {
            static const gh_cv_category cats[] = {
                GH_CV_PATHOGENIC, GH_CV_LIKELY_PATHOGENIC,
            };
            for (size_t c = 0; c < 2; c++) {
                for (size_t i = 0; i < in->clinvar->counts[cats[c]]; i++) {
                    const gh_cv_finding *f =
                        &in->clinvar->by_category[cats[c]][i];
                    const char *gene_upper = f->gene ? f->gene : "";

                    /* Python upper-cases the ClinVar gene before comparing
                     * against the group's genes. */
                    const char *folded = upper_copy(arena, gene_upper);
                    bool by_gene = false;
                    for (size_t k = 0; k < group->ngenes && !by_gene; k++)
                        if (strcmp(folded, group->genes[k]) == 0)
                            by_gene = true;

                    const char *primary = primary_traits(arena, f->traits);
                    bool by_keyword = !by_gene && f->gold_stars >= 2 &&
                        any_keyword(primary, group->keywords, group->nkeywords);

                    if (!by_gene && !by_keyword)
                        continue;

                    signal_count++;
                    pathogenic = true;
                    if (ndetails < 32) {
                        const char *tail = by_gene
                            ? (f->traits && *f->traits ? f->traits : "unknown")
                            : primary;
                        size_t len = strlen(f->gene) + strlen(tail) + 40;
                        char *line = gh_alloc(arena, len);
                        if (by_gene)
                            snprintf(line, len, "Pathogenic variant: %s (%s)",
                                     f->gene, tail);
                        else
                            snprintf(line, len, "Pathogenic variant: %s (%.60s)",
                                     f->gene, tail);
                        details[ndetails++] = line;
                    }
                }
            }

            /* Risk factors contribute a signal, but only one detail line. */
            bool noted = false;
            for (size_t i = 0; i < in->clinvar->counts[GH_CV_RISK_FACTOR]; i++) {
                const gh_cv_finding *f =
                    &in->clinvar->by_category[GH_CV_RISK_FACTOR][i];
                const char *traits = lower_copy(arena, f->traits);
                if (!any_keyword(traits, group->keywords, group->nkeywords))
                    continue;
                signal_count++;
                if (!noted && ndetails < 32) {
                    const char *title_lower = lower_copy(arena, group->title);
                    size_t len = strlen(title_lower) + 48;
                    char *line = gh_alloc(arena, len);
                    snprintf(line, len,
                             "ClinVar risk factors related to %s", title_lower);
                    details[ndetails++] = line;
                    noted = true;
                }
            }
        }

        /* 4. Gene-gene interactions touching this group. */
        for (size_t i = 0; i < in->nepistasis; i++) {
            const gh_epi_result *e = &in->epistasis[i];
            bool overlaps = false;
            for (size_t j = 0; j < e->ngenes && !overlaps; j++)
                if (list_has(group->genes, group->ngenes, e->genes[j].gene))
                    overlaps = true;
            if (!overlaps)
                continue;
            signal_count++;
            epistasis = true;
            if (ndetails < 32) {
                size_t len = strlen(e->name) + strlen(e->risk_level) + 24;
                char *line = gh_alloc(arena, len);
                snprintf(line, len, "Epistasis: %s (%s)", e->name,
                         e->risk_level);
                details[ndetails++] = line;
            }
        }

        if (signal_count == 0)
            continue;

        gh_priority *p = &ranked[n].item;
        p->id = group->id;
        p->title = group->title;
        p->priority = gh_compute_priority(gene_signals, pathogenic, prs_high,
                                          prs_elevated, high_magnitude,
                                          epistasis);
        p->why = build_why(arena, details, ndetails);
        p->actions = group->actions;
        p->nactions = group->nactions;
        p->doctor_note = group->doctor_note;
        p->monitoring = group->monitoring;
        p->nmonitoring = group->nmonitoring;
        p->signal_count = signal_count;
        ranked[n].order = n;
        n++;
    }

    qsort(ranked, n, sizeof(*ranked), compare_priority);

    out->priorities = gh_calloc(arena, n ? n : 1, sizeof(*out->priorities));
    for (size_t i = 0; i < n; i++)
        out->priorities[i] = ranked[i].item;
    out->npriorities = n;
}

/* Add clinical-context actions that the group's own list does not cover. */
static void overlay_clinical_actions(gh_recommendations *out,
                                     const gh_analysis *analysis)
{
    for (size_t i = 0; i < out->npriorities; i++) {
        gh_priority *p = &out->priorities[i];

        const gh_risk_group *group = NULL;
        for (size_t g = 0; g < GH_RISK_GROUP_COUNT; g++)
            if (strcmp(GH_RISK_GROUPS[g].id, p->id) == 0)
                group = &GH_RISK_GROUPS[g];
        if (!group)
            continue;

        for (size_t k = 0; k < group->ngenes; k++) {
            const gh_finding *f = gene_finding(analysis, group->genes[k]);
            if (!f)
                continue;
            const gh_clinical_context *ctx =
                gh_clinical_context_for(f->gene, f->status);
            if (!ctx || !ctx->nactions)
                continue;

            for (size_t a = 0; a < ctx->nactions; a++) {
                const char *action = ctx->actions[a];
                if (list_has(p->actions, p->nactions, action) ||
                    list_has(p->clinical_actions, p->nclinical_actions, action))
                    continue;
                if (p->nclinical_actions < GH_MAX_CLINICAL_ACTIONS)
                    p->clinical_actions[p->nclinical_actions++] = action;
            }
        }
    }
}

/* ------------------------------------------------------------------ */
/* Drug card                                                           */
/* ------------------------------------------------------------------ */

static const char *const DRUG_GENES[] = {
    "CYP2C9", "CYP2C19", "CYP1A2", "CYP2D6", "CYP3A4",
    "VKORC1", "DPYD", "TPMT", "UGT1A1", "SLCO1B1",
};

static gh_drug_card_gene *card_for(gh_drug_card_gene *card, size_t *ncard,
                                   const char *gene)
{
    for (size_t i = 0; i < *ncard; i++)
        if (strcmp(card[i].gene, gene) == 0)
            return &card[i];
    card[*ncard].gene = gene;
    return &card[(*ncard)++];
}

static void build_drug_card(gh_recommendations *out, gh_arena *arena,
                            const gh_recommendation_inputs *in)
{
    size_t capacity = sizeof(DRUG_GENES) / sizeof(DRUG_GENES[0]) +
                      in->nstars + 32;
    gh_drug_card_gene *card = gh_calloc(arena, capacity, sizeof(*card));
    size_t ncard = 0;

    size_t entry_capacity = capacity * 8;
    gh_drug_entry *pool = gh_calloc(arena, entry_capacity, sizeof(*pool));
    size_t npool = 0;

    /* Entries are appended per gene, so each gene's run is contiguous only
     * if nothing else interleaves; they are collected per gene below. */
    for (size_t pass = 0; pass < 3; pass++) {
        if (pass == 0 && in->stars) {
            /* Actionable star-allele phenotypes come first, as in Python. */
            for (size_t i = 0; i < in->nstars; i++) {
                const gh_star_result *r = &in->stars[i];
                const char *ph = r->phenotype ? r->phenotype : "normal";
                if (strcmp(ph, "poor") && strcmp(ph, "intermediate") &&
                    strcmp(ph, "ultrarapid") && strcmp(ph, "rapid"))
                    continue;
                gh_drug_card_gene *slot = card_for(card, &ncard, r->gene);
                size_t len = strlen(ph) + 16;
                char *status = gh_alloc(arena, len);
                snprintf(status, len, "%s metabolizer", ph);
                pool[npool] = (gh_drug_entry){"", r->diplotype, status,
                                              r->clinical_note,
                                              "Star Allele (CPIC)"};
                if (!slot->entries)
                    slot->entries = &pool[npool];
                slot->nentries++;
                npool++;
            }
        } else if (pass == 1 && in->analysis) {
            for (size_t i = 0; i < in->analysis->nfindings; i++) {
                const gh_finding *f = &in->analysis->findings[i];
                if (!f->gene ||
                    !list_has(DRUG_GENES,
                              sizeof(DRUG_GENES) / sizeof(DRUG_GENES[0]),
                              f->gene))
                    continue;
                gh_drug_card_gene *slot = card_for(card, &ncard, f->gene);
                pool[npool] = (gh_drug_entry){f->rsid, f->genotype, f->status,
                                              f->description,
                                              "Lifestyle/SNP DB"};
                if (!slot->entries)
                    slot->entries = &pool[npool];
                slot->nentries++;
                npool++;
            }
        } else if (pass == 2 && in->clinvar && in->clinvar->loaded) {
            for (size_t i = 0; i < in->clinvar->counts[GH_CV_DRUG_RESPONSE]; i++) {
                const gh_cv_finding *f =
                    &in->clinvar->by_category[GH_CV_DRUG_RESPONSE][i];
                if (!f->gene || !*f->gene)
                    continue;
                gh_drug_card_gene *slot = card_for(card, &ncard, f->gene);
                pool[npool] = (gh_drug_entry){
                    f->rsid, f->user_genotype, "",
                    (f->traits && *f->traits) ? f->traits
                                              : "Drug response variant",
                    "ClinVar"};
                if (!slot->entries)
                    slot->entries = &pool[npool];
                slot->nentries++;
                npool++;
            }
        }
    }

    out->drug_card = card;
    out->ndrug_card = ncard;
}

/* ------------------------------------------------------------------ */
/* Good news, monitoring, referrals, insights                          */
/* ------------------------------------------------------------------ */

static void build_good_news(gh_recommendations *out, gh_arena *arena,
                            const gh_recommendation_inputs *in)
{
    size_t capacity = (in->analysis ? in->analysis->nfindings : 0) + 64;
    gh_good_news *good = gh_calloc(arena, capacity, sizeof(*good));
    size_t n = 0;

    if (in->analysis) {
        for (size_t i = 0; i < in->analysis->nfindings; i++) {
            const gh_finding *f = &in->analysis->findings[i];
            const char *s = f->status ? f->status : "";
            if (strcmp(s, "protective") && strcmp(s, "longevity") &&
                strcmp(s, "optimal") && strcmp(s, "fast"))
                continue;
            good[n].gene = f->gene;
            good[n].description = f->description;
            n++;
        }
    }
    if (in->clinvar && in->clinvar->loaded) {
        for (size_t i = 0; i < in->clinvar->counts[GH_CV_PROTECTIVE]; i++) {
            const gh_cv_finding *f =
                &in->clinvar->by_category[GH_CV_PROTECTIVE][i];
            good[n].gene = (f->gene && *f->gene) ? f->gene : "Unknown";
            good[n].description = (f->traits && *f->traits) ? f->traits
                                                            : "Protective variant";
            n++;
        }
    }
    out->good_news = good;
    out->ngood_news = n;
}

/* Lower is more frequent; anything unrecognised sorts last. */
static int frequency_rank(const char *frequency)
{
    static const char *const order[] = {
        "Weekly (home)", "Monthly", "Quarterly", "Annually",
        "Once (baseline)", "Per guidelines",
    };
    for (size_t i = 0; i < sizeof(order) / sizeof(order[0]); i++)
        if (strcmp(order[i], frequency) == 0)
            return (int)i;
    return 10;
}

typedef struct {
    gh_monitoring_item item;
    size_t order;
} ranked_monitoring;

static int compare_monitoring(const void *a, const void *b)
{
    const ranked_monitoring *x = a, *y = b;
    int rx = frequency_rank(x->item.frequency);
    int ry = frequency_rank(y->item.frequency);
    if (rx != ry)
        return rx < ry ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

static void build_monitoring(gh_recommendations *out, gh_arena *arena)
{
    size_t capacity = 1;
    for (size_t i = 0; i < out->npriorities; i++)
        capacity += out->priorities[i].nmonitoring;

    ranked_monitoring *items = gh_calloc(arena, capacity, sizeof(*items));
    size_t n = 0;

    for (size_t i = 0; i < out->npriorities; i++) {
        const gh_priority *p = &out->priorities[i];
        for (size_t m = 0; m < p->nmonitoring; m++) {
            const gh_monitoring_item *item = &p->monitoring[m];
            /* One entry per test, keeping whichever asks for it most often. */
            size_t existing = n;
            for (size_t k = 0; k < n; k++)
                if (strcmp(items[k].item.test, item->test) == 0) {
                    existing = k;
                    break;
                }
            if (existing == n) {
                items[n].item = *item;
                items[n].order = n;
                n++;
            } else if (frequency_rank(item->frequency) <
                       frequency_rank(items[existing].item.frequency)) {
                items[existing].item = *item;
            }
        }
    }

    qsort(items, n, sizeof(*items), compare_monitoring);

    out->monitoring_schedule = gh_calloc(arena, n ? n : 1,
                                         sizeof(*out->monitoring_schedule));
    for (size_t i = 0; i < n; i++)
        out->monitoring_schedule[i] = items[i].item;
    out->nmonitoring = n;
}

static void build_referrals(gh_recommendations *out, gh_arena *arena,
                            const gh_recommendation_inputs *in)
{
    gh_referral *referrals =
        gh_calloc(arena, GH_SPECIALIST_REFERRAL_COUNT, sizeof(*referrals));
    size_t n = 0;

    size_t pathogenic_count = 0;
    if (in->clinvar && in->clinvar->loaded)
        pathogenic_count = in->clinvar->counts[GH_CV_PATHOGENIC] +
                           in->clinvar->counts[GH_CV_LIKELY_PATHOGENIC];

    size_t nacmg = in->acmg ? in->acmg->nfindings : 0;

    /* Genes with an actionable metabolizer phenotype, in table order. */
    const char *actionable[32];
    size_t nactionable = 0;
    for (size_t i = 0; i < in->nstars && nactionable < 32; i++) {
        const char *ph = in->stars[i].phenotype;
        if (ph && (strcmp(ph, "poor") == 0 || strcmp(ph, "intermediate") == 0 ||
                   strcmp(ph, "rapid") == 0 || strcmp(ph, "ultrarapid") == 0))
            actionable[nactionable++] = in->stars[i].gene;
    }

    for (size_t r = 0; r < GH_SPECIALIST_REFERRAL_COUNT; r++) {
        const gh_specialist_referral *ref = &GH_SPECIALIST_REFERRALS[r];
        bool triggered = false;
        const char *reason = ref->reason_template;

        for (size_t t = 0; t < ref->ntriggers && !triggered; t++) {
            const gh_referral_trigger *trig = &ref->triggers[t];

            if (strcmp(trig->source, "acmg") == 0) {
                if (trig->match_any && nacmg) {
                    triggered = true;
                } else if (trig->ngenes && in->acmg) {
                    /* Sorted order, so a multi-gene trigger names the same
                     * gene every run. */
                    for (size_t gi = 0; gi < trig->ngenes && !triggered; gi++) {
                        for (size_t fi = 0; fi < nacmg; fi++) {
                            const gh_acmg_finding *f = &in->acmg->findings[fi];
                            if (strcmp(f->finding->gene, trig->genes[gi]) != 0)
                                continue;
                            triggered = true;
                            const char *condition =
                                (f->finding->traits && *f->finding->traits)
                                    ? f->finding->traits : "variant detected";
                            size_t len = strlen(ref->reason_template) +
                                         strlen(trig->genes[gi]) +
                                         strlen(condition) + 8;
                            char *filled = gh_alloc(arena, len);
                            /* The template is "{gene}: {condition} ...". */
                            const char *tail = strstr(ref->reason_template,
                                                      "{condition}");
                            snprintf(filled, len, "%s: %s%s",
                                     trig->genes[gi], condition,
                                     tail ? tail + strlen("{condition}") : "");
                            reason = filled;
                            break;
                        }
                    }
                }
            } else if (strcmp(trig->source, "pathogenic_count") == 0) {
                if (pathogenic_count >= (size_t)(trig->minimum ? trig->minimum : 1))
                    triggered = true;
            } else if (strcmp(trig->source, "disease") == 0) {
                if (in->clinvar && in->clinvar->loaded) {
                    static const gh_cv_category cats[] = {
                        GH_CV_PATHOGENIC, GH_CV_LIKELY_PATHOGENIC,
                    };
                    for (size_t c = 0; c < 2 && !triggered; c++)
                        for (size_t i = 0;
                             i < in->clinvar->counts[cats[c]] && !triggered; i++) {
                            const char *traits = lower_copy(
                                arena, in->clinvar->by_category[cats[c]][i].traits);
                            if (any_keyword(traits, trig->keywords,
                                            trig->nkeywords))
                                triggered = true;
                        }
                }
            } else if (strcmp(trig->source, "priority") == 0) {
                for (size_t i = 0; i < out->npriorities; i++) {
                    if (strcmp(out->priorities[i].id, trig->group) != 0)
                        continue;
                    if (priority_rank(out->priorities[i].priority) <=
                        priority_rank(trig->min_priority))
                        triggered = true;
                }
            } else if (strcmp(trig->source, "gene_status") == 0) {
                const gh_finding *f = gene_finding(in->analysis, trig->gene);
                if (f && f->status &&
                    list_has(trig->statuses, trig->nstatuses, f->status))
                    triggered = true;
            } else if (strcmp(trig->source, "star_alleles_actionable") == 0) {
                if (trig->match_any && nactionable) {
                    triggered = true;
                    size_t len = strlen(ref->reason_template) + 8;
                    for (size_t i = 0; i < nactionable; i++)
                        len += strlen(actionable[i]) + 2;
                    char *filled = gh_alloc(arena, len);
                    const char *slot = strstr(ref->reason_template, "{genes}");
                    size_t prefix = slot ? (size_t)(slot - ref->reason_template)
                                         : strlen(ref->reason_template);
                    size_t at = (size_t)snprintf(filled, len, "%.*s",
                                                 (int)prefix,
                                                 ref->reason_template);
                    for (size_t i = 0; i < nactionable; i++)
                        at += (size_t)snprintf(filled + at, len - at, "%s%s",
                                               i ? ", " : "", actionable[i]);
                    if (slot)
                        snprintf(filled + at, len - at, "%s",
                                 slot + strlen("{genes}"));
                    reason = filled;
                }
            }
        }

        if (!triggered)
            continue;
        referrals[n].specialist = ref->title;
        referrals[n].reason = reason;
        referrals[n].urgency = ref->urgency;
        n++;
    }

    out->referrals = referrals;
    out->nreferrals = n;
}

typedef struct {
    gh_clinical_insight item;
    size_t order;
} ranked_insight;

/* Strongest finding first; Python sorts on -magnitude with a stable sort. */
static int compare_insight(const void *a, const void *b)
{
    const ranked_insight *x = a, *y = b;
    if (x->item.magnitude != y->item.magnitude)
        return x->item.magnitude > y->item.magnitude ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

static void build_insights(gh_recommendations *out, gh_arena *arena,
                           const gh_analysis *analysis)
{
    size_t capacity = analysis ? analysis->nfindings + 1 : 1;
    ranked_insight *items = gh_calloc(arena, capacity, sizeof(*items));
    size_t n = 0;

    /* One insight per gene, from that gene's strongest finding. */
    const char *seen[256];
    size_t nseen = 0;

    for (size_t i = 0; analysis && i < analysis->nfindings; i++) {
        const gh_finding *f = &analysis->findings[i];
        if (!f->gene || !*f->gene || f->magnitude < 1)
            continue;
        if (list_has(seen, nseen, f->gene))
            continue;

        const gh_finding *best = gene_finding(analysis, f->gene);
        if (!best || best->magnitude < 1)
            continue;
        if (nseen < 256)
            seen[nseen++] = f->gene;

        const gh_clinical_context *ctx =
            gh_clinical_context_for(best->gene, best->status);
        if (!ctx)
            continue;

        gh_clinical_insight *insight = &items[n].item;
        insight->gene = best->gene;
        insight->status = best->status;
        insight->magnitude = best->magnitude;
        insight->context = ctx;

        for (size_t p = 0; p < GH_PATHWAY_COUNT && insight->npathways < 8; p++)
            for (size_t k = 0; k < GH_PATHWAYS[p].ngenes; k++)
                if (strcmp(GH_PATHWAYS[p].genes[k], best->gene) == 0) {
                    insight->pathways[insight->npathways++] =
                        GH_PATHWAYS[p].name;
                    break;
                }

        items[n].order = n;
        n++;
    }

    qsort(items, n, sizeof(*items), compare_insight);

    out->insights = gh_calloc(arena, n ? n : 1, sizeof(*out->insights));
    for (size_t i = 0; i < n; i++)
        out->insights[i] = items[i].item;
    out->ninsights = n;
}

/* ------------------------------------------------------------------ */

void gh_generate_recommendations(gh_recommendations *out, gh_arena *arena,
                                 const gh_recommendation_inputs *in)
{
    memset(out, 0, sizeof(*out));

    build_priorities(out, arena, in);
    overlay_clinical_actions(out, in->analysis);
    build_drug_card(out, arena, in);
    build_monitoring(out, arena);
    build_good_news(out, arena, in);
    build_referrals(out, arena, in);
    build_insights(out, arena, in->analysis);
}

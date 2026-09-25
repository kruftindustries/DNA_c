#include "gh_insights.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Gene status table                                                   */
/* ------------------------------------------------------------------ */

/* The status the insights engine attributes to a gene, and how strong that
 * call is. Built from the lifestyle findings, then overlaid with APOE and
 * star-allele results, which take precedence. Entries stay in the order
 * their gene was first seen, which is what the Python's dict preserves and
 * what the protective list iterates. */
typedef struct {
    const char *gene;
    const char *status;
    int magnitude;
} gene_state;

typedef struct {
    gene_state *items;
    size_t count;
    size_t capacity;
} gene_states;

static gene_state *state_find(gene_states *s, const char *gene)
{
    for (size_t i = 0; i < s->count; i++)
        if (strcmp(s->items[i].gene, gene) == 0)
            return &s->items[i];
    return NULL;
}

static const char *state_status(const gene_states *s, const char *gene)
{
    for (size_t i = 0; i < s->count; i++)
        if (strcmp(s->items[i].gene, gene) == 0)
            return s->items[i].status;
    return NULL;
}

static int state_magnitude(const gene_states *s, const char *gene)
{
    for (size_t i = 0; i < s->count; i++)
        if (strcmp(s->items[i].gene, gene) == 0)
            return s->items[i].magnitude;
    return 0;
}

static void state_set(gene_states *s, const char *gene, const char *status,
                      int magnitude)
{
    gene_state *existing = state_find(s, gene);
    if (existing) {
        existing->status = status;
        existing->magnitude = magnitude;
        return;
    }
    if (s->count < s->capacity) {
        s->items[s->count].gene = gene;
        s->items[s->count].status = status;
        s->items[s->count].magnitude = magnitude;
        s->count++;
    }
}

/* Title-case a status for display: underscores become spaces and each word
 * is capitalised, matching Python's replace("_", " ").title(). */
static const char *titlecase(gh_arena *arena, const char *s)
{
    size_t n = s ? strlen(s) : 0;
    char *out = gh_alloc(arena, n + 1);
    bool at_word_start = true;
    for (size_t i = 0; i < n; i++) {
        char c = s[i] == '_' ? ' ' : s[i];
        if (isalpha((unsigned char)c)) {
            out[i] = at_word_start ? (char)toupper((unsigned char)c)
                                   : (char)tolower((unsigned char)c);
            at_word_start = false;
        } else {
            out[i] = c;
            at_word_start = true;
        }
    }
    out[n] = '\0';
    return out;
}

/* ------------------------------------------------------------------ */
/* Single-gene insights                                                */
/* ------------------------------------------------------------------ */

typedef struct {
    gh_single_gene_result item;
    size_t order;
} ranked_single;

/* Strongest first; Python sorts on -magnitude with a stable sort, so ties
 * keep the table's order. */
static int compare_single(const void *a, const void *b)
{
    const ranked_single *x = a, *y = b;
    if (x->item.magnitude != y->item.magnitude)
        return x->item.magnitude > y->item.magnitude ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

static void build_single_gene(gh_insights *out, gh_arena *arena,
                              const gene_states *states)
{
    size_t capacity = 1;
    for (size_t i = 0; i < GH_SINGLE_GENE_INSIGHT_COUNT; i++)
        capacity += GH_SINGLE_GENE_INSIGHTS[i].nentries;

    ranked_single *ranked = gh_calloc(arena, capacity, sizeof(*ranked));
    size_t n = 0;

    for (size_t i = 0; i < GH_SINGLE_GENE_INSIGHT_COUNT; i++) {
        const gh_single_gene_insight *insight = &GH_SINGLE_GENE_INSIGHTS[i];
        const char *status = state_status(states, insight->gene);
        if (!status || strcmp(status, insight->status) != 0)
            continue;

        for (size_t e = 0; e < insight->nentries; e++) {
            ranked[n].item.entry = &insight->entries[e];
            ranked[n].item.gene = insight->gene;
            ranked[n].item.status = insight->status;
            ranked[n].item.magnitude = state_magnitude(states, insight->gene);
            ranked[n].order = n;
            n++;
        }
    }

    qsort(ranked, n, sizeof(*ranked), compare_single);

    out->single_gene = gh_calloc(arena, n ? n : 1, sizeof(*out->single_gene));
    for (size_t i = 0; i < n; i++)
        out->single_gene[i] = ranked[i].item;
    out->nsingle_gene = n;
}

/* ------------------------------------------------------------------ */
/* Multi-gene narratives                                               */
/* ------------------------------------------------------------------ */

static bool matches_any(const gh_gene_match *match, const char *status)
{
    if (!status)
        return false;
    for (size_t i = 0; i < match->nstatuses; i++)
        if (strcmp(match->statuses[i], status) == 0)
            return true;
    return false;
}

/* Substitute {matched_count}; anything else is left alone. */
static const char *fill_narrative(gh_arena *arena, const char *text,
                                  size_t matched_count)
{
    const char *slot = strstr(text, "{matched_count}");
    if (!slot)
        return text;

    char number[24];
    int written = snprintf(number, sizeof(number), "%zu", matched_count);
    size_t prefix = (size_t)(slot - text);
    const char *tail = slot + strlen("{matched_count}");

    size_t len = prefix + (size_t)written + strlen(tail) + 1;
    char *out = gh_alloc(arena, len);
    snprintf(out, len, "%.*s%s%s", (int)prefix, text, number, tail);
    return out;
}

static void build_narratives(gh_insights *out, gh_arena *arena,
                             const gene_states *states)
{
    gh_narrative_result *results =
        gh_calloc(arena, GH_NARRATIVE_COUNT, sizeof(*results));
    size_t n = 0;

    for (size_t i = 0; i < GH_NARRATIVE_COUNT; i++) {
        const gh_narrative_pattern *pattern = &GH_NARRATIVES[i];

        const char *matched[GH_MAX_MATCHED_GENES];
        size_t nmatched = 0;
        bool has_required = false;

        for (size_t r = 0; r < pattern->nrequired_genes; r++) {
            const gh_gene_match *m = &pattern->required_genes[r];
            if (!matches_any(m, state_status(states, m->gene)))
                continue;
            has_required = true;
            if (nmatched < GH_MAX_MATCHED_GENES)
                matched[nmatched++] = m->gene;
        }
        for (size_t o = 0; o < pattern->noptional_genes; o++) {
            const gh_gene_match *m = &pattern->optional_genes[o];
            if (!matches_any(m, state_status(states, m->gene)))
                continue;
            if (nmatched < GH_MAX_MATCHED_GENES)
                matched[nmatched++] = m->gene;
        }

        /* Enough genes overall, and at least one of them required. */
        if ((int)nmatched < pattern->min_matches || !has_required)
            continue;

        gh_narrative_result *result = &results[n];
        result->id = pattern->id;
        result->title = pattern->title;
        result->nmatched = nmatched;
        for (size_t k = 0; k < nmatched; k++)
            result->matched_genes[k] = matched[k];
        result->narrative = fill_narrative(arena, pattern->narrative, nmatched);
        result->practical = pattern->practical;
        result->references = pattern->references;
        result->nreferences = pattern->nreferences;
        n++;
    }

    out->narratives = results;
    out->nnarratives = n;
}

/* ------------------------------------------------------------------ */
/* Highlights                                                          */
/* ------------------------------------------------------------------ */

static void add_highlight(gh_insights *out, const char *title,
                          const char *detail, const char *type)
{
    if (out->nhighlights >= GH_MAX_HIGHLIGHTS)
        return;
    out->highlights[out->nhighlights].title = title;
    out->highlights[out->nhighlights].detail = detail;
    out->highlights[out->nhighlights].type = type;
    out->nhighlights++;
}

typedef struct {
    const gh_finding *finding;
    size_t order;
} ranked_finding_ref;

static int compare_high_mag(const void *a, const void *b)
{
    const ranked_finding_ref *x = a, *y = b;
    if (x->finding->magnitude != y->finding->magnitude)
        return x->finding->magnitude > y->finding->magnitude ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

static void build_highlights(gh_insights *out, gh_arena *arena,
                             const gh_insight_inputs *in,
                             const gene_states *states)
{
    /* APOE e2 is protective and therefore interesting. */
    if (in->apoe && in->apoe->apoe_type && strstr(in->apoe->apoe_type, "e2")) {
        size_t len = strlen(in->apoe->apoe_type) + 48;
        char *title = gh_alloc(arena, len);
        snprintf(title, len, "APOE %s \342\200\224 Alzheimer's protection",
                 in->apoe->apoe_type);
        add_highlight(out, title,
                      "APOE e2 carriers have roughly half the Alzheimer's risk.",
                      "protective");
    }

    const char *cetp = state_status(states, "CETP");
    if (cetp && (strcmp(cetp, "favorable") == 0 ||
                 strcmp(cetp, "highly_favorable") == 0))
        add_highlight(out, "CETP longevity variant",
                      "Associated with exceptional longevity and favorable HDL.",
                      "protective");

    /* The two most significant pathogenic ClinVar findings. */
    if (in->clinvar && in->clinvar->loaded) {
        size_t limit = in->clinvar->counts[GH_CV_PATHOGENIC];
        if (limit > 2)
            limit = 2;
        for (size_t i = 0; i < limit; i++) {
            const gh_cv_finding *f = &in->clinvar->by_category[GH_CV_PATHOGENIC][i];
            const char *gene = (f->gene && *f->gene) ? f->gene : "Unknown";

            const char *trait = "Pathogenic";
            if (f->traits && *f->traits) {
                const char *semi = strchr(f->traits, ';');
                size_t n = semi ? (size_t)(semi - f->traits) : strlen(f->traits);
                const char *start = f->traits;
                while (n && isspace((unsigned char)*start)) {
                    start++;
                    n--;
                }
                while (n && isspace((unsigned char)start[n - 1]))
                    n--;
                trait = gh_strndup(arena, start, n);
            }

            size_t tlen = strlen(gene) + strlen(trait) + 8;
            char *title = gh_alloc(arena, tlen);
            snprintf(title, tlen, "%s \342\200\224 %s", gene, trait);

            char *detail = gh_alloc(arena, 64);
            snprintf(detail, 64, "Pathogenic variant detected (%d/4 stars).",
                     f->gold_stars);
            add_highlight(out, title, detail, "clinical");
        }
    }

    /* Actionable metabolizer phenotypes. */
    for (size_t i = 0; i < in->nstars; i++) {
        const gh_star_result *r = &in->stars[i];
        const char *ph = r->phenotype ? r->phenotype : "";
        if (strcmp(ph, "poor") != 0 && strcmp(ph, "intermediate") != 0)
            continue;
        if (out->nhighlights >= GH_MAX_HIGHLIGHTS)
            break;

        const char *pretty = titlecase(arena, ph);
        size_t len = strlen(r->gene) + strlen(r->diplotype) +
                     strlen(pretty) + 32;
        char *title = gh_alloc(arena, len);
        snprintf(title, len, "%s %s \342\200\224 %s Metabolizer",
                 r->gene, r->diplotype, pretty);
        add_highlight(out, title, r->clinical_note, "pharmacogenomic");
    }

    /* Whatever high-magnitude lifestyle findings still fit, strongest first
     * and skipping genes already named above. */
    if (!in->analysis)
        return;

    ranked_finding_ref *high =
        gh_calloc(arena, in->analysis->nfindings + 1, sizeof(*high));
    size_t nhigh = 0;
    for (size_t i = 0; i < in->analysis->nfindings; i++) {
        if (in->analysis->findings[i].magnitude < 3)
            continue;
        high[nhigh].finding = &in->analysis->findings[i];
        high[nhigh].order = nhigh;
        nhigh++;
    }
    qsort(high, nhigh, sizeof(*high), compare_high_mag);

    for (size_t i = 0; i < nhigh && out->nhighlights < GH_MAX_HIGHLIGHTS; i++) {
        const gh_finding *f = high[i].finding;
        const char *gene = f->gene ? f->gene : "";

        bool covered = false;
        for (size_t h = 0; h < out->nhighlights && !covered; h++)
            if (out->highlights[h].title && *gene &&
                strstr(out->highlights[h].title, gene))
                covered = true;
        if (covered)
            continue;

        const char *status = titlecase(arena, f->status ? f->status : "");
        size_t len = strlen(gene) + strlen(status) + 8;
        char *title = gh_alloc(arena, len);
        snprintf(title, len, "%s \342\200\224 %s", gene, status);
        add_highlight(out, title, f->description, "lifestyle");
    }
}

/* ------------------------------------------------------------------ */
/* Protective findings                                                 */
/* ------------------------------------------------------------------ */

static const char *const PROTECTIVE_STATUSES[] = {
    "protective", "longevity", "optimal", "favorable",
    "highly_favorable", "fast", "e2_carrier",
};

static void build_protective(gh_insights *out, gh_arena *arena,
                             const gene_states *states)
{
    gh_protective_finding *found =
        gh_calloc(arena, states->count ? states->count : 1, sizeof(*found));
    size_t n = 0;

    for (size_t i = 0; i < states->count; i++) {
        const gene_state *st = &states->items[i];

        bool protective = false;
        for (size_t k = 0; k < sizeof(PROTECTIVE_STATUSES) /
                               sizeof(PROTECTIVE_STATUSES[0]); k++)
            if (strcmp(PROTECTIVE_STATUSES[k], st->status) == 0)
                protective = true;
        if (!protective)
            continue;

        /* Only genes that have an insight entry are reported. */
        const gh_single_gene_insight *insight = NULL;
        for (size_t k = 0; k < GH_SINGLE_GENE_INSIGHT_COUNT; k++)
            if (strcmp(GH_SINGLE_GENE_INSIGHTS[k].gene, st->gene) == 0 &&
                strcmp(GH_SINGLE_GENE_INSIGHTS[k].status, st->status) == 0) {
                insight = &GH_SINGLE_GENE_INSIGHTS[k];
                break;
            }
        if (!insight || !insight->nentries)
            continue;

        found[n].gene = st->gene;
        found[n].status = st->status;
        found[n].title = insight->entries[0].title;
        found[n].finding = insight->entries[0].finding;
        found[n].reference = insight->entries[0].reference;
        n++;
    }

    out->protective = found;
    out->nprotective = n;
}

/* ------------------------------------------------------------------ */

void gh_generate_insights(gh_insights *out, gh_arena *arena,
                          const gh_insight_inputs *in)
{
    memset(out, 0, sizeof(*out));

    size_t capacity = (in->analysis ? in->analysis->nfindings : 0) +
                      in->nstars + 8;
    gene_states states = {
        gh_calloc(arena, capacity, sizeof(gene_state)), 0, capacity,
    };

    /* Strongest finding wins per gene. */
    for (size_t i = 0; in->analysis && i < in->analysis->nfindings; i++) {
        const gh_finding *f = &in->analysis->findings[i];
        if (!f->gene || !*f->gene || !f->status || !*f->status)
            continue;
        gene_state *existing = state_find(&states, f->gene);
        if (!existing)
            state_set(&states, f->gene, f->status, f->magnitude);
        else if (f->magnitude > existing->magnitude) {
            existing->status = f->status;
            existing->magnitude = f->magnitude;
        }
    }

    /* APOE only contributes when it carries an e2 allele. */
    if (in->apoe && in->apoe->apoe_type &&
        strcmp(in->apoe->apoe_type, "Unknown") != 0 &&
        strstr(in->apoe->apoe_type, "e2"))
        state_set(&states, "APOE", "e2_carrier", 2);

    /* Actionable metabolizer phenotypes override whatever the lifestyle
     * findings said about the same gene. */
    for (size_t i = 0; i < in->nstars; i++) {
        const char *ph = in->stars[i].phenotype;
        if (!ph)
            continue;
        if (strcmp(ph, "poor") == 0 || strcmp(ph, "intermediate") == 0 ||
            strcmp(ph, "rapid") == 0 || strcmp(ph, "ultrarapid") == 0)
            state_set(&states, in->stars[i].gene, ph, 3);
    }

    build_single_gene(out, arena, &states);
    build_narratives(out, arena, &states);
    build_highlights(out, arena, in, &states);
    build_protective(out, arena, &states);
}

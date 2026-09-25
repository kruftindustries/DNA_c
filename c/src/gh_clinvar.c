#include "gh_clinvar.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gh_tsv.h"

const char *const GH_CV_CATEGORY_NAMES[GH_CV_CATEGORY_COUNT] = {
    "pathogenic", "likely_pathogenic", "risk_factor",
    "drug_response", "protective",
};

/* ------------------------------------------------------------------ */
/* Helpers                                                             */
/* ------------------------------------------------------------------ */

static const char *kv_lookup(const gh_kv *table, size_t n, const char *key)
{
    for (size_t i = 0; i < n; i++)
        if (strcmp(table[i].key, key) == 0)
            return table[i].value;
    return NULL;
}

static bool string_list_has(const char *const *list, size_t n, const char *key)
{
    for (size_t i = 0; i < n; i++)
        if (strcmp(list[i], key) == 0)
            return true;
    return false;
}

/* Lower-case copy with underscores turned into spaces, matching the
 * `.lower().replace("_", " ")` the Python applies to clinical significance. */
static char *normalise_significance(gh_arena *arena, const char *s)
{
    size_t n = strlen(s);
    char *out = gh_alloc(arena, n + 1);
    for (size_t i = 0; i < n; i++) {
        char c = s[i] == '_' ? ' ' : s[i];
        out[i] = (char)tolower((unsigned char)c);
    }
    out[n] = '\0';
    return out;
}

static char *upper_copy(gh_arena *arena, const char *s)
{
    size_t n = strlen(s);
    char *out = gh_alloc(arena, n + 1);
    for (size_t i = 0; i < n; i++)
        out[i] = (char)toupper((unsigned char)s[i]);
    out[n] = '\0';
    return out;
}

static char *lower_copy(gh_arena *arena, const char *s)
{
    size_t n = strlen(s);
    char *out = gh_alloc(arena, n + 1);
    for (size_t i = 0; i < n; i++)
        out[i] = (char)tolower((unsigned char)s[i]);
    out[n] = '\0';
    return out;
}

/* Python's int(value) with a default of 0 on anything unparseable. */
static int safe_int(const char *s)
{
    if (!s || !*s)
        return 0;
    char *end = NULL;
    long v = strtol(s, &end, 10);
    if (end == s || (end && *end != '\0'))
        return 0;
    return (int)v;
}

void gh_classify_zygosity(bool is_homozygous, bool is_heterozygous,
                          const char *inheritance,
                          const char **status, const char **description)
{
    if (is_homozygous) {
        *status = "AFFECTED";
        *description = "Homozygous for variant";
        return;
    }
    if (is_heterozygous) {
        /* The inheritance field is free text from ClinVar, so this is a
         * substring test on a lower-cased copy, as in the Python. */
        if (inheritance && strstr(inheritance, "recessive")) {
            *status = "CARRIER";
            *description = "Heterozygous carrier (recessive)";
        } else if (inheritance && strstr(inheritance, "dominant")) {
            *status = "AFFECTED";
            *description = "Heterozygous (dominant)";
        } else {
            *status = "HETEROZYGOUS";
            *description = "Heterozygous (inheritance unclear)";
        }
        return;
    }
    *status = "UNKNOWN";
    *description = "Zygosity unclear";
}

/* ------------------------------------------------------------------ */
/* ClinVar scan                                                        */
/* ------------------------------------------------------------------ */

/* Growable per-category finding lists. */
typedef struct {
    gh_cv_finding **items;
    size_t count;
    size_t cap;
} finding_list;

static void list_push(finding_list *l, gh_arena *arena, gh_cv_finding *f)
{
    if (l->count == l->cap) {
        size_t cap = l->cap ? l->cap * 2 : 32;
        gh_cv_finding **grown = gh_calloc(arena, cap, sizeof(*grown));
        if (l->count)
            memcpy(grown, l->items, l->count * sizeof(*grown));
        l->items = grown;
        l->cap = cap;
    }
    l->items[l->count++] = f;
}

bool gh_clinvar_analyze(gh_clinvar_result *out, gh_arena *arena,
                        const gh_genome *g, const char *clinvar_path)
{
    memset(out, 0, sizeof(*out));

    /* ClinVar runs to hundreds of megabytes, so it is streamed rather than
     * slurped. Fields are only valid until the next record, so everything
     * kept below is copied into the arena. */
    gh_tsv tsv;
    if (!gh_tsv_open_streaming(&tsv, arena, clinvar_path))
        return false;
    if (!gh_tsv_read_header(&tsv)) {
        gh_tsv_close(&tsv);
        return false;
    }

    int c_chrom = gh_tsv_col(&tsv, "chrom");
    int c_pos = gh_tsv_col(&tsv, "pos");
    int c_ref = gh_tsv_col(&tsv, "ref");
    int c_alt = gh_tsv_col(&tsv, "alt");
    int c_sig = gh_tsv_col(&tsv, "clinical_significance");
    int c_review = gh_tsv_col(&tsv, "review_status");
    int c_stars = gh_tsv_col(&tsv, "gold_stars");
    int c_traits = gh_tsv_col(&tsv, "all_traits");
    int c_symbol = gh_tsv_col(&tsv, "symbol");
    int c_inherit = gh_tsv_col(&tsv, "inheritance_modes");
    int c_hgvs_p = gh_tsv_col(&tsv, "hgvs_p");
    int c_hgvs_c = gh_tsv_col(&tsv, "hgvs_c");
    int c_consequence = gh_tsv_col(&tsv, "molecular_consequence");
    int c_xrefs = gh_tsv_col(&tsv, "xrefs");

    finding_list lists[GH_CV_CATEGORY_COUNT];
    memset(lists, 0, sizeof(lists));

    /* Reused for the chrom:pos lookup key. */
    char key[128];

    while (gh_tsv_next(&tsv)) {
        out->total_clinvar++;

        const char *chrom = gh_tsv_at(&tsv, c_chrom);
        const char *pos = gh_tsv_at(&tsv, c_pos);
        if (!*chrom || !*pos)
            continue;

        snprintf(key, sizeof(key), "%s:%s", chrom, pos);
        const gh_variant *v = gh_map_get(&g->by_position, key);
        if (!v)
            continue;

        out->matched++;

        const char *ref = gh_tsv_at(&tsv, c_ref);
        const char *alt = gh_tsv_at(&tsv, c_alt);

        /* Indels cannot be represented reliably by array genotypes, so only
         * single-base substitutions are considered. Without this filter the
         * report fills with false positives in well-known genes. */
        if (strlen(ref) != 1 || strlen(alt) != 1)
            continue;

        const char *genotype = v->genotype;
        bool has_variant = strchr(genotype, alt[0]) != NULL;
        bool is_homozygous = strlen(genotype) == 2 &&
                             genotype[0] == alt[0] && genotype[1] == alt[0];
        bool is_heterozygous = has_variant && !is_homozygous;
        bool has_ref_only = strlen(genotype) == 2 &&
                            genotype[0] == ref[0] && genotype[1] == ref[0];

        if (has_ref_only || !has_variant)
            continue;

        const char *raw_sig = gh_tsv_at(&tsv, c_sig);
        const char *sig = normalise_significance(arena, raw_sig);

        gh_cv_finding *f = gh_alloc(arena, sizeof(*f));
        f->chromosome = gh_strdup(arena, chrom);
        f->position = gh_strdup(arena, pos);
        f->rsid = v->rsid;
        f->gene = gh_strdup(arena, gh_tsv_at(&tsv, c_symbol));
        f->ref = gh_strdup(arena, ref);
        f->alt = gh_strdup(arena, alt);
        f->user_genotype = genotype;   /* owned by the genome */
        f->is_homozygous = is_homozygous;
        f->is_heterozygous = is_heterozygous;
        f->clinical_significance = gh_strdup(arena, raw_sig);
        f->review_status = gh_strdup(arena, gh_tsv_at(&tsv, c_review));
        f->gold_stars = safe_int(gh_tsv_at(&tsv, c_stars));
        f->traits = gh_strdup(arena, gh_tsv_at(&tsv, c_traits));
        f->inheritance = gh_strdup(arena, gh_tsv_at(&tsv, c_inherit));
        f->hgvs_p = gh_strdup(arena, gh_tsv_at(&tsv, c_hgvs_p));
        f->hgvs_c = gh_strdup(arena, gh_tsv_at(&tsv, c_hgvs_c));
        f->molecular_consequence = gh_strdup(arena, gh_tsv_at(&tsv, c_consequence));
        f->xrefs = gh_strdup(arena, gh_tsv_at(&tsv, c_xrefs));

        const char *inherit_lower = lower_copy(arena, f->inheritance);
        gh_classify_zygosity(is_homozygous, is_heterozygous, inherit_lower,
                             &f->zygosity_status, &f->zygosity);

        /* Categories are tested in order; "pathogenic" excludes the likely
         * and conflicting spellings so those fall through. */
        int category = -1;
        if (strstr(sig, "pathogenic") && !strstr(sig, "likely") &&
            !strstr(sig, "conflict"))
            category = GH_CV_PATHOGENIC;
        else if (strstr(sig, "likely pathogenic"))
            category = GH_CV_LIKELY_PATHOGENIC;
        else if (strstr(sig, "risk factor"))
            category = GH_CV_RISK_FACTOR;
        else if (strstr(sig, "drug response"))
            category = GH_CV_DRUG_RESPONSE;
        else if (strstr(sig, "protective"))
            category = GH_CV_PROTECTIVE;

        if (category < 0)
            continue;

        f->category = (gh_cv_category)category;
        list_push(&lists[category], arena, f);
    }

    out->by_category = gh_calloc(arena, GH_CV_CATEGORY_COUNT,
                                 sizeof(*out->by_category));
    for (size_t i = 0; i < GH_CV_CATEGORY_COUNT; i++) {
        out->by_category[i] = NULL;
        out->counts[i] = lists[i].count;
        if (lists[i].count) {
            gh_cv_finding *flat =
                gh_calloc(arena, lists[i].count, sizeof(*flat));
            for (size_t j = 0; j < lists[i].count; j++)
                flat[j] = *lists[i].items[j];
            out->by_category[i] = flat;
        }
    }

    gh_tsv_close(&tsv);
    out->loaded = true;
    return true;
}

/* ------------------------------------------------------------------ */
/* ACMG                                                                */
/* ------------------------------------------------------------------ */

bool gh_acmg_is_actionable_gene(const char *gene_upper)
{
    return string_list_has(GH_ACMG_GENES, GH_ACMG_GENE_COUNT, gene_upper);
}

/* Sorted by descending gold stars, then gene name -- Python's
 * sort(key=lambda x: (-gold_stars, gene)). */
static int compare_acmg(const void *a, const void *b)
{
    const gh_acmg_finding *x = a, *y = b;
    if (x->finding->gold_stars != y->finding->gold_stars)
        return x->finding->gold_stars > y->finding->gold_stars ? -1 : 1;
    return strcmp(x->finding->gene, y->finding->gene);
}

void gh_flag_acmg(gh_acmg_result *out, gh_arena *arena,
                  const gh_clinvar_result *clinvar)
{
    memset(out, 0, sizeof(*out));
    out->genes_screened = GH_ACMG_GENE_COUNT;

    if (!clinvar || !clinvar->loaded) {
        out->summary = "No ClinVar data available for ACMG screening.";
        return;
    }

    size_t capacity = clinvar->counts[GH_CV_PATHOGENIC]
                    + clinvar->counts[GH_CV_LIKELY_PATHOGENIC];
    out->findings = gh_calloc(arena, capacity + 1, sizeof(*out->findings));

    const char **genes_seen = gh_calloc(arena, capacity + 1, sizeof(*genes_seen));
    size_t ngenes = 0;

    static const gh_cv_category categories[] = {
        GH_CV_PATHOGENIC, GH_CV_LIKELY_PATHOGENIC,
    };

    for (size_t c = 0; c < sizeof(categories) / sizeof(categories[0]); c++) {
        gh_cv_category cat = categories[c];
        for (size_t i = 0; i < clinvar->counts[cat]; i++) {
            const gh_cv_finding *f = &clinvar->by_category[cat][i];
            const char *gene = upper_copy(arena, f->gene ? f->gene : "");
            if (!gh_acmg_is_actionable_gene(gene))
                continue;

            if (!string_list_has(genes_seen, ngenes, gene))
                genes_seen[ngenes++] = gene;

            const char *actionability = kv_lookup(
                GH_ACMG_ACTIONABILITY, GH_ACMG_ACTIONABILITY_COUNT, gene);

            gh_acmg_finding *entry = &out->findings[out->nfindings++];
            entry->finding = f;
            entry->acmg_category = GH_CV_CATEGORY_NAMES[cat];
            entry->actionability = actionability ? actionability
                                                 : GH_ACMG_DEFAULT_ACTIONABILITY;
        }
    }

    qsort(out->findings, out->nfindings, sizeof(*out->findings), compare_acmg);
    out->genes_with_variants = ngenes;

    if (out->nfindings) {
        char *summary = gh_alloc(arena, 160);
        snprintf(summary, 160,
                 "%zu variant(s) found in %zu ACMG-recommended gene(s). "
                 "Genetic counseling recommended.",
                 out->nfindings, ngenes);
        out->summary = summary;
    } else {
        out->summary =
            "No pathogenic/likely pathogenic variants in ACMG SF v3.2 genes.";
    }
}

/* ------------------------------------------------------------------ */
/* Carrier screening                                                   */
/* ------------------------------------------------------------------ */

/* Inheritance pattern for a gene, preferring the curated table.
 *
 * ClinVar's inheritance field in this extract holds variant origin
 * (germline/somatic) rather than a mode of inheritance, so the curated table
 * is the only usable source for the genes it covers. Mirrors
 * carrier_screen.py:_resolve_inheritance. */
static const char *resolve_inheritance(gh_arena *arena, const char *gene_upper,
                                       const char *clinvar_inheritance)
{
    const char *curated = kv_lookup(GH_GENE_INHERITANCE,
                                    GH_GENE_INHERITANCE_COUNT, gene_upper);
    if (curated)
        return curated;
    if (clinvar_inheritance && *clinvar_inheritance)
        return lower_copy(arena, clinvar_inheritance);
    return "unknown";
}

/* True when the person carries one copy of a recessively acting variant.
 *
 * Homozygotes are affected rather than carriers, and a heterozygote for a
 * dominant condition is affected too. X-linked counts: a heterozygous female
 * carries an X-linked recessive condition without being affected. Matching is
 * case-insensitive because the curated table capitalises "X-linked".
 * Mirrors carrier_screen.py:_is_carrier. */
static bool is_carrier(gh_arena *arena, const gh_cv_finding *f,
                       const char *inheritance)
{
    if (f->is_homozygous || !f->is_heterozygous)
        return false;
    const char *lowered = lower_copy(arena, inheritance);
    return strstr(lowered, "recessive") != NULL ||
           strstr(lowered, "x-linked") != NULL;
}

/* First semicolon-separated trait, trimmed; "Unknown condition" when empty. */
static const char *first_trait(gh_arena *arena, const char *traits)
{
    if (!traits)
        return "Unknown condition";

    const char *end = strchr(traits, ';');
    size_t n = end ? (size_t)(end - traits) : strlen(traits);

    const char *start = traits;
    while (n && isspace((unsigned char)*start)) {
        start++;
        n--;
    }
    while (n && isspace((unsigned char)start[n - 1]))
        n--;

    if (n == 0)
        return "Unknown condition";
    return gh_strndup(arena, start, n);
}

void gh_organize_carriers(gh_carrier_result *out, gh_arena *arena,
                          const gh_clinvar_result *clinvar)
{
    memset(out, 0, sizeof(*out));

    if (!clinvar || !clinvar->loaded)
        return;

    size_t capacity = clinvar->counts[GH_CV_PATHOGENIC]
                    + clinvar->counts[GH_CV_LIKELY_PATHOGENIC];
    if (capacity == 0)
        return;

    out->carriers = gh_calloc(arena, capacity, sizeof(*out->carriers));
    out->systems = gh_calloc(arena, capacity, sizeof(*out->systems));
    out->system_counts = gh_calloc(arena, capacity, sizeof(*out->system_counts));

    static const gh_cv_category categories[] = {
        GH_CV_PATHOGENIC, GH_CV_LIKELY_PATHOGENIC,
    };

    for (size_t c = 0; c < sizeof(categories) / sizeof(categories[0]); c++) {
        gh_cv_category cat = categories[c];
        for (size_t i = 0; i < clinvar->counts[cat]; i++) {
            const gh_cv_finding *f = &clinvar->by_category[cat][i];

            const char *gene = f->gene && *f->gene ? f->gene : "Unknown";
            const char *gene_upper = upper_copy(arena, gene);

            /* Resolve inheritance before deciding, so the curated table gets
             * a chance to apply. Doing it the other way round meant the
             * decision was made on ClinVar's unusable origin field. */
            const char *inheritance =
                resolve_inheritance(arena, gene_upper, f->inheritance);
            if (!is_carrier(arena, f, inheritance))
                continue;

            const char *system = kv_lookup(GH_GENE_SYSTEMS,
                                           GH_GENE_SYSTEM_COUNT, gene_upper);
            if (!system)
                system = "Other";

            const char *lowered = lower_copy(arena, inheritance);
            const char *note = "";
            if (strstr(lowered, "recessive"))
                note = "If partner is also a carrier, each child has "
                       "a 25% chance of being affected.";
            else if (strstr(lowered, "x-linked"))
                note = "X-linked: carrier females may pass to sons "
                       "(50% chance affected) and daughters (50% chance carrier).";

            gh_carrier *entry = &out->carriers[out->ncarriers++];
            entry->gene = gene;
            entry->condition = first_trait(arena, f->traits);
            entry->inheritance = inheritance;
            entry->system = system;
            entry->reproductive_note = note;
            entry->couples_relevant = string_list_has(
                GH_COUPLES_GENES, GH_COUPLES_GENE_COUNT, gene_upper);
            entry->rsid = f->rsid ? f->rsid : "";
            entry->genotype = f->user_genotype ? f->user_genotype : "";
            entry->gold_stars = f->gold_stars;

            if (entry->couples_relevant)
                out->ncouples_relevant++;

            /* Systems are grouped in first-seen order, as the Python's
             * setdefault does. */
            size_t s = 0;
            for (; s < out->nsystems; s++)
                if (strcmp(out->systems[s], system) == 0)
                    break;
            if (s == out->nsystems) {
                out->systems[s] = system;
                out->system_counts[s] = 0;
                out->nsystems++;
            }
            out->system_counts[s]++;
        }
    }
}

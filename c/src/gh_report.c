/* localtime_r is POSIX.1-2008; request it explicitly under -std=c11. */
#define _POSIX_C_SOURCE 200809L

#include "gh_report.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* M_PI is not in standard C; define it rather than relying on an extension. */
#define GH_PI 3.14159265358979323846

/* ------------------------------------------------------------------ */
/* Small helpers                                                       */
/* ------------------------------------------------------------------ */

/* Append an integer with thousands separators, like Python's "{:,}". */
static void put_grouped(gh_strbuf *sb, size_t value)
{
    char digits[32];
    int n = snprintf(digits, sizeof(digits), "%zu", value);
    for (int i = 0; i < n; i++) {
        if (i > 0 && (n - i) % 3 == 0)
            gh_sb_putc(sb, ',');
        gh_sb_putc(sb, digits[i]);
    }
}

const char *gh_report_mag_class(int magnitude)
{
    if (magnitude >= 3)
        return "high";
    if (magnitude == 2)
        return "mod";
    if (magnitude == 1)
        return "low";
    return "info";
}

const char *gh_report_titlecase(gh_arena *arena, const char *status)
{
    if (!status)
        return "";
    size_t n = strlen(status);
    char *out = gh_alloc(arena, n + 1);

    /* Python's str.title() upper-cases the first letter of each run of
     * letters and lower-cases the rest, after "_" becomes " ". */
    bool at_word_start = true;
    for (size_t i = 0; i < n; i++) {
        char c = status[i] == '_' ? ' ' : status[i];
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

const char *gh_report_eli5_gene(const char *gene)
{
    if (!gene)
        return "";
    for (size_t i = 0; i < GH_ELI5_GENE_COUNT; i++)
        if (strcmp(GH_ELI5_GENES[i].key, gene) == 0)
            return GH_ELI5_GENES[i].value;
    return "";
}

void gh_report_db_links(gh_strbuf *sb, const char *rsid)
{
    if (!rsid || strncmp(rsid, "rs", 2) != 0)
        return;

    static const struct {
        const char *label;
        const char *url;
    } sources[] = {
        {"dbSNP",    "https://www.ncbi.nlm.nih.gov/snp/"},
        {"ClinVar",  "https://www.ncbi.nlm.nih.gov/clinvar/?term="},
        {"SNPedia",  "https://www.snpedia.com/index.php/"},
        {"PharmGKB", "https://www.pharmgkb.org/search?query="},
    };

    gh_sb_puts(sb, "<span class=\"db-links\">");
    for (size_t i = 0; i < sizeof(sources) / sizeof(sources[0]); i++) {
        if (i)
            gh_sb_puts(sb, " &middot; ");
        gh_sb_printf(sb, "<a href=\"%s%s\" target=\"_blank\" rel=\"noopener\">%s</a>",
                     sources[i].url, rsid, sources[i].label);
    }
    gh_sb_puts(sb, "</span>");
}

void gh_report_paper_refs(gh_strbuf *sb, const char *rsid)
{
    if (!rsid)
        return;

    bool first = true;
    for (size_t i = 0; i < GH_PAPER_REF_COUNT; i++) {
        const gh_paper_ref *r = &GH_PAPER_REFS[i];
        if (strcmp(r->rsid, rsid) != 0)
            continue;
        if (first) {
            gh_sb_puts(sb, "<div class=\"paper-refs\">References: ");
            first = false;
        } else {
            gh_sb_puts(sb, " | ");
        }
        gh_sb_printf(sb,
            "<a href=\"https://pubmed.ncbi.nlm.nih.gov/%s/\" "
            "target=\"_blank\" rel=\"noopener\">%s (%d)</a>",
            r->pmid, r->title, r->year);
    }
    if (!first)
        gh_sb_puts(sb, "</div>");
}

/* ------------------------------------------------------------------ */
/* Charts                                                              */
/* ------------------------------------------------------------------ */

void gh_report_svg_impact_bar(gh_strbuf *sb, const gh_analysis *a)
{
    size_t high = 0, mod = 0, low = 0, info = 0;
    for (size_t i = 0; i < a->nfindings; i++) {
        int m = a->findings[i].magnitude;
        if (m >= 3) high++;
        else if (m == 2) mod++;
        else if (m == 1) low++;
        else info++;
    }

    size_t max_val = high;
    if (mod > max_val) max_val = mod;
    if (low > max_val) max_val = low;
    if (info > max_val) max_val = info;
    if (max_val == 0) max_val = 1;

    const int bar_w = 280;
    const struct {
        int y;
        size_t val;
        const char *color;
        const char *label;
    } bars[] = {
        {10,  high, GH_C.red,   "High"},
        {48,  mod,  GH_C.amber, "Moderate"},
        {86,  low,  GH_C.green, "Low"},
        {124, info, GH_C.slate, "Info"},
    };

    gh_sb_puts(sb,
        "<svg viewBox=\"0 0 440 165\" class=\"chart\" role=\"img\" "
        "aria-label=\"Impact distribution bar chart\">"
        "<title>Impact Distribution</title>");

    for (size_t i = 0; i < sizeof(bars) / sizeof(bars[0]); i++) {
        int w = (int)((double)bars[i].val / (double)max_val * bar_w);
        gh_sb_printf(sb,
            "<rect x=\"110\" y=\"%d\" width=\"%d\" height=\"28\" rx=\"4\" "
            "fill=\"%s\" opacity=\"0.85\"/>"
            "<text x=\"105\" y=\"%d\" text-anchor=\"end\" "
            "fill=\"currentColor\" font-size=\"13\">%s</text>"
            "<text x=\"%d\" y=\"%d\" fill=\"currentColor\" "
            "font-size=\"13\" font-weight=\"bold\">%zu</text>",
            bars[i].y, w, bars[i].color,
            bars[i].y + 19, bars[i].label,
            115 + w, bars[i].y + 19, bars[i].val);
    }
    gh_sb_puts(sb, "</svg>");
}

/* Category tally for the donut, sorted by descending count. */
typedef struct {
    const char *name;
    size_t count;
    size_t first_seen;   /* preserves a stable order among equal counts */
} cat_count;

static int compare_cat(const void *a, const void *b)
{
    const cat_count *x = a, *y = b;
    if (x->count != y->count)
        return x->count > y->count ? -1 : 1;
    /* Python sorts by -count only; its sort is stable, so ties keep
     * dict insertion order. Mirror that with first-seen order. */
    return (x->first_seen > y->first_seen) - (x->first_seen < y->first_seen);
}

void gh_report_svg_category_donut(gh_strbuf *sb, gh_arena *arena,
                                  const gh_analysis *a)
{
    if (a->nfindings == 0)
        return;

    cat_count *cats = gh_calloc(arena, a->nfindings, sizeof(*cats));
    size_t ncats = 0, total = 0;

    for (size_t i = 0; i < a->nfindings; i++) {
        const char *name = a->findings[i].category
            ? a->findings[i].category : "Other";
        size_t j = 0;
        for (; j < ncats; j++)
            if (strcmp(cats[j].name, name) == 0)
                break;
        if (j == ncats) {
            cats[ncats].name = name;
            cats[ncats].count = 0;
            cats[ncats].first_seen = ncats;
            ncats++;
        }
        cats[j].count++;
        total++;
    }
    if (ncats == 0)
        return;

    qsort(cats, ncats, sizeof(*cats), compare_cat);

    const double cx = 120, cy = 120, r = 90, inner_r = 55;
    double angle = -90;

    size_t height = 15 + ncats * 22 + 10;
    if (height < 250)
        height = 250;

    gh_sb_printf(sb,
        "<svg viewBox=\"0 0 440 %zu\" class=\"chart\" role=\"img\" "
        "aria-label=\"Category donut chart\">"
        "<title>Findings by Category</title>", height);

    /* Arc paths first, so the legend and centre text draw on top. */
    for (size_t i = 0; i < ncats; i++) {
        double sweep = (double)cats[i].count / (double)total * 360.0;
        double start_rad = angle * GH_PI / 180.0;
        double end_rad = (angle + sweep) * GH_PI / 180.0;

        double x1 = cx + r * cos(start_rad), y1 = cy + r * sin(start_rad);
        double x2 = cx + r * cos(end_rad),   y2 = cy + r * sin(end_rad);
        double ix1 = cx + inner_r * cos(start_rad);
        double iy1 = cy + inner_r * sin(start_rad);
        double ix2 = cx + inner_r * cos(end_rad);
        double iy2 = cy + inner_r * sin(end_rad);

        int large = sweep > 180 ? 1 : 0;
        gh_sb_printf(sb,
            "<path d=\"M %.1f %.1f L %.1f %.1f A %g %g 0 %d 1 %.1f %.1f "
            "L %.1f %.1f A %g %g 0 %d 0 %.1f %.1f Z\" "
            "fill=\"%s\" opacity=\"0.85\"/>",
            ix1, iy1, x1, y1, r, r, large, x2, y2,
            ix2, iy2, inner_r, inner_r, large, ix1, iy1,
            GH_CHART_COLORS[i % GH_CHART_COLOR_COUNT]);
        angle += sweep;
    }

    gh_sb_printf(sb,
        "<text x=\"%g\" y=\"%g\" text-anchor=\"middle\" "
        "fill=\"currentColor\" font-size=\"24\" font-weight=\"bold\">%zu</text>"
        "<text x=\"%g\" y=\"%g\" text-anchor=\"middle\" "
        "fill=\"currentColor\" font-size=\"11\">findings</text>",
        cx, cy - 5, total, cx, cy + 15);

    for (size_t i = 0; i < ncats; i++) {
        size_t ly = 15 + i * 22;
        gh_sb_printf(sb,
            "<rect x=\"260\" y=\"%zu\" width=\"14\" height=\"14\" rx=\"3\" "
            "fill=\"%s\"/>"
            "<text x=\"280\" y=\"%zu\" fill=\"currentColor\" font-size=\"12\">",
            ly, GH_CHART_COLORS[i % GH_CHART_COLOR_COUNT], ly + 12);
        gh_sb_put_escaped(sb, cats[i].name);
        gh_sb_printf(sb, " (%zu)</text>", cats[i].count);
    }

    gh_sb_puts(sb, "</svg>");
}

/* ------------------------------------------------------------------ */
/* Sections                                                            */
/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */
/* Template filling                                                    */
/* ------------------------------------------------------------------ */

const char *gh_report_template(gh_arena *arena)
{
    size_t total = 0;
    for (size_t i = 0; GH_HTML_TEMPLATE_PARTS[i]; i++)
        total += strlen(GH_HTML_TEMPLATE_PARTS[i]);

    char *out = gh_alloc(arena, total + 1);
    size_t at = 0;
    for (size_t i = 0; GH_HTML_TEMPLATE_PARTS[i]; i++) {
        size_t n = strlen(GH_HTML_TEMPLATE_PARTS[i]);
        memcpy(out + at, GH_HTML_TEMPLATE_PARTS[i], n);
        at += n;
    }
    out[total] = '\0';
    return out;
}

void gh_report_fill_template(gh_strbuf *sb, const char *template_text,
                             const char *const *names,
                             const char *const *values, size_t count)
{
    for (const char *p = template_text; *p;) {
        if (*p != '{') {
            gh_sb_putc(sb, *p++);
            continue;
        }

        const char *close = strchr(p + 1, '}');
        if (!close) {
            gh_sb_putc(sb, *p++);
            continue;
        }

        size_t len = (size_t)(close - (p + 1));
        size_t i = 0;
        for (; i < count; i++)
            if (strlen(names[i]) == len && strncmp(names[i], p + 1, len) == 0)
                break;

        if (i < count) {
            gh_sb_puts(sb, values[i] ? values[i] : "");
            p = close + 1;
        } else {
            /* Unknown placeholder: leave it visible rather than blanking it. */
            gh_sb_putc(sb, *p++);
        }
    }
}

/* ------------------------------------------------------------------ */
/* Whole report                                                        */
/* ------------------------------------------------------------------ */

/* Render one section into an arena-held string. */
#define SECTION(var, call)                                                  \
    const char *var;                                                        \
    do {                                                                    \
        gh_strbuf tmp_;                                                     \
        gh_sb_init(&tmp_);                                                  \
        call;                                                               \
        var = gh_strndup(arena, tmp_.data ? tmp_.data : "", tmp_.len);      \
        gh_sb_free(&tmp_);                                                  \
    } while (0)

void gh_report_render(gh_strbuf *sb, gh_arena *arena, const gh_analysis *a,
                      const gh_quality *q, const gh_report_scorers *sc,
                      const gh_report_opts *opts)
{
    gh_strbuf tmp;

    /* f" — {_esc(subject_name)}" -- an em dash, not a hyphen. */
    char subject_title[256] = "";
    if (opts && opts->subject_name && *opts->subject_name) {
        gh_sb_init(&tmp);
        gh_sb_puts(&tmp, " \xe2\x80\x94 ");
        gh_sb_put_escaped(&tmp, opts->subject_name);
        snprintf(subject_title, sizeof(subject_title), "%s", tmp.data);
        gh_sb_free(&tmp);
    }

    char date_buf[64];
    if (opts && opts->generated_date) {
        snprintf(date_buf, sizeof(date_buf), "%s", opts->generated_date);
    } else {
        time_t now = time(NULL);
        struct tm tm_buf;
#ifdef _WIN32
        /* The Windows C runtime's localtime keeps a per-thread buffer. */
        tm_buf = *localtime(&now);
#else
        localtime_r(&now, &tm_buf);
#endif
        strftime(date_buf, sizeof(date_buf), "%Y-%m-%d %H:%M", &tm_buf);
    }

    char num_findings[32], num_pharmgkb[32];
    snprintf(num_findings, sizeof(num_findings), "%zu", a->nfindings);
    snprintf(num_pharmgkb, sizeof(num_pharmgkb), "%zu", a->ndrug_findings);

    /* The template writes this one as "{total_snps:,}" -- a format spec, so
     * the placeholder text the filler has to match includes the ":,". */
    gh_sb_init(&tmp);
    put_grouped(&tmp, a->total_snps);
    const char *total_snps = gh_strndup(arena, tmp.data ? tmp.data : "0", tmp.len);
    gh_sb_free(&tmp);

    /* Low-coverage warning for genomes with few SNPs, e.g. WGS-derived. */
    const char *low_coverage = "";
    if (a->total_snps < 100000) {
        double pct = a->total_snps ? (double)a->total_snps / 600000.0 * 100.0 : 0.0;
        gh_sb_init(&tmp);
        gh_sb_puts(&tmp,
            "<div class=\"doctor-callout\" style=\"border-color:var(--warn);text-align:left;font-weight:normal\">"
            "<strong style=\"color:var(--warn)\">Low Data Coverage</strong>: Only ");
        put_grouped(&tmp, a->total_snps);
        gh_sb_printf(&tmp,
            " SNPs loaded (%.0f%% of a typical 23andMe/30x WGS dataset). "
            "Many analysis modules may show \"no data available\". "
            "For comprehensive results, use a 30x whole genome sequencing dataset or 23andMe raw data.</div>",
            pct);
        low_coverage = gh_strndup(arena, tmp.data, tmp.len);
        gh_sb_free(&tmp);
    }

    SECTION(key_findings, gh_report_key_findings(&tmp_, arena, a, sc));
    SECTION(action_plan, gh_report_action_plan(&tmp_, arena, sc));
    SECTION(drugs, gh_report_drug_guide(&tmp_, arena, a, sc));
    SECTION(disease_risk, gh_report_disease_risk(&tmp_, arena, sc));
    SECTION(body_profile, gh_report_body_profile(&tmp_, arena, a, sc));
    SECTION(mental, gh_report_mental_health(&tmp_, arena, sc));
    SECTION(clinical, gh_report_clinical_detail(&tmp_, arena, a, sc));
    SECTION(ancestry, gh_report_ancestry(&tmp_, arena, sc));
    SECTION(nutrigenomics, gh_report_nutrigenomics(&tmp_, arena, sc));
    SECTION(quality, gh_report_quality(&tmp_, q));
    SECTION(doctor_card, gh_report_doctor_card(&tmp_, arena, sc));
    SECTION(references, gh_report_references(&tmp_, a));

    const char *names[] = {
        "total_snps:,",
        "subject_title", "generated_date", "num_findings", "num_pharmgkb",
        "low_coverage_warning", "key_findings_content", "action_plan_content",
        "drug_guide_content", "disease_risk_content", "body_profile_content",
        "mental_health_content", "clinical_detail_content", "ancestry_content",
        "nutrigenomics_content", "quality_content", "doctor_card_content",
        "references_content",
    };
    const char *values[] = {
        total_snps,
        subject_title, date_buf, num_findings, num_pharmgkb,
        low_coverage, key_findings, action_plan,
        drugs, disease_risk, body_profile,
        mental, clinical, ancestry,
        nutrigenomics, quality, doctor_card,
        references,
    };

    gh_report_fill_template(sb, gh_report_template(arena), names, values,
                            sizeof(names) / sizeof(names[0]));
}

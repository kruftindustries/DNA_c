/* The twelve report sections.
 *
 * Each function mirrors one build_* in genetic_health/reports/enhanced_html.py
 * closely enough that the rendered HTML is byte-identical, which
 * tools/reportdiff.py checks. That fidelity extends to the reference's
 * inconsistencies -- some builders escape their inputs and some do not, some
 * literals are UTF-8 characters and some are entities -- because the goal is
 * the same report, not a tidier one. Where the Python has a quirk worth
 * knowing about it is called out in a comment at the point it is reproduced.
 */
#define _POSIX_C_SOURCE 200809L

#include "gh_report.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/* math.radians multiplies by a precomputed pi/180; dividing afterwards can
 * differ in the last bit, which the .1f coordinates occasionally expose. */
static const double DEG = M_PI / 180.0;

/* ------------------------------------------------------------------ */
/* Text helpers matching the Python's                                  */
/* ------------------------------------------------------------------ */

static void raw(gh_strbuf *sb, const char *s)
{
    if (s)
        gh_sb_puts(sb, s);
}

/* html.escape(str(text)) if text else "" */
static void esc(gh_strbuf *sb, const char *s)
{
    if (s && *s)
        gh_sb_put_escaped(sb, s);
}

static const char *or_empty(const char *s)
{
    return s ? s : "";
}

/* Python's str.title(): the first letter after any non-letter is
 * upper-cased, every other letter lower-cased. Optionally replaces "_" with
 * " " first, for the `.replace("_", " ").title()` idiom. */
static const char *py_title(gh_arena *arena, const char *s, bool underscores)
{
    if (!s)
        return "";
    size_t n = strlen(s);
    char *out = gh_alloc(arena, n + 1);
    bool word_start = true;
    for (size_t i = 0; i < n; i++) {
        char c = (underscores && s[i] == '_') ? ' ' : s[i];
        if (isalpha((unsigned char)c)) {
            out[i] = word_start ? (char)toupper((unsigned char)c)
                                : (char)tolower((unsigned char)c);
            word_start = false;
        } else {
            out[i] = c;
            word_start = true;
        }
    }
    out[n] = '\0';
    return out;
}

static const char *py_upper(gh_arena *arena, const char *s)
{
    if (!s)
        return "";
    size_t n = strlen(s);
    char *out = gh_alloc(arena, n + 1);
    for (size_t i = 0; i < n; i++)
        out[i] = (char)toupper((unsigned char)s[i]);
    out[n] = '\0';
    return out;
}

/* str.isupper(): at least one cased character, and none of them lower. */
static bool py_isupper(const char *s)
{
    bool cased = false;
    for (const char *p = s; *p; p++) {
        if (islower((unsigned char)*p))
            return false;
        if (isupper((unsigned char)*p))
            cased = true;
    }
    return cased;
}

/* len() counts code points, not bytes. */
static size_t py_len(const char *s)
{
    size_t n = 0;
    for (const unsigned char *p = (const unsigned char *)s; *p; p++)
        if ((*p & 0xC0) != 0x80)
            n++;
    return n;
}

static bool ends_with(const char *s, const char *suffix)
{
    size_t n = strlen(s), m = strlen(suffix);
    return n >= m && strcmp(s + n - m, suffix) == 0;
}

static char *strip_copy(gh_arena *arena, const char *start, const char *end)
{
    while (start < end && isspace((unsigned char)*start))
        start++;
    while (end > start && isspace((unsigned char)end[-1]))
        end--;
    return gh_strndup(arena, start, (size_t)(end - start));
}

static char *replace_all(gh_arena *arena, const char *s, const char *from,
                         const char *to)
{
    size_t flen = strlen(from), tlen = strlen(to), count = 0;
    for (const char *p = s; (p = strstr(p, from)); p += flen)
        count++;
    char *out = gh_alloc(arena, strlen(s) + count * (tlen > flen ? tlen - flen : 0) + 1);
    char *w = out;
    for (const char *p = s;;) {
        const char *hit = strstr(p, from);
        if (!hit) {
            strcpy(w, p);
            break;
        }
        memcpy(w, p, (size_t)(hit - p));
        w += hit - p;
        memcpy(w, to, tlen);
        w += tlen;
        p = hit + flen;
    }
    return out;
}

static const char *lower_copy(gh_arena *arena, const char *s)
{
    size_t n = strlen(s);
    char *out = gh_alloc(arena, n + 1);
    for (size_t i = 0; i < n; i++)
        out[i] = (char)tolower((unsigned char)s[i]);
    out[n] = '\0';
    return out;
}

const char *gh_report_clean_condition(gh_arena *arena, const char *rawtext)
{
    if (!rawtext || !*rawtext)
        return "Unknown";

    /* Split on "|" and ";" and drop the unhelpful segments. */
    const char *text = replace_all(arena, rawtext, "|", ";");
    const char *parts[64];
    size_t nparts = 0;
    for (const char *p = text;;) {
        const char *semi = strchr(p, ';');
        const char *end = semi ? semi : p + strlen(p);
        char *segment = strip_copy(arena, p, end);
        if (*segment && nparts < 64) {
            const char *low = lower_copy(arena, segment);
            bool skip = strcmp(low, "not provided") == 0 ||
                        strcmp(low, "not specified") == 0 ||
                        strcmp(low, "unknown") == 0 ||
                        strcmp(low, "see cases") == 0;
            if (!skip && ends_with(low, "-related disorder") && nparts > 0)
                skip = true;
            if (!skip)
                parts[nparts++] = segment;
        }
        if (!semi)
            break;
        p = semi + 1;
    }

    if (nparts == 0) {
        /* raw.split("|")[0].split(";")[0].strip() */
        const char *end = rawtext + strcspn(rawtext, "|");
        const char *end2 = rawtext + strcspn(rawtext, ";");
        return strip_copy(arena, rawtext, end2 < end ? end2 : end);
    }

    const char *best = parts[0];
    for (size_t i = 0; i < nparts; i++) {
        if (!py_isupper(parts[i]) && py_len(parts[i]) > py_len(best) / 2) {
            best = parts[i];
            break;
        }
    }

    if (py_isupper(best)) {
        best = py_title(arena, replace_all(arena, best, ",", ", "), false);
        best = replace_all(arena, best, " To", "");
        best = replace_all(arena, best, "Resistance", "resistance");
        best = replace_all(arena, best, "Susceptibility", "susceptibility");
    }
    return best;
}

/* _dedup_phrases: unique "; "-separated phrases, order kept. */
static const char *dedup_phrases(gh_arena *arena, const char *text)
{
    if (!text || !strchr(text, ';'))
        return text;

    const char *seen[128];
    size_t nseen = 0;
    gh_strbuf out;
    gh_sb_init(&out);

    for (const char *p = text;;) {
        const char *sep = strstr(p, "; ");
        const char *end = sep ? sep : p + strlen(p);
        char *phrase = strip_copy(arena, p, end);
        if (*phrase) {
            bool dup = false;
            for (size_t i = 0; i < nseen && !dup; i++)
                dup = strcmp(seen[i], phrase) == 0;
            if (!dup) {
                if (nseen < 128)
                    seen[nseen++] = phrase;
                if (out.len)
                    gh_sb_puts(&out, "; ");
                gh_sb_puts(&out, phrase);
            }
        }
        if (!sep)
            break;
        p = sep + 2;
    }
    const char *result = gh_strndup(arena, out.data ? out.data : "", out.len);
    gh_sb_free(&out);
    return result;
}

const char *gh_report_clean_why(gh_arena *arena, const char *text)
{
    if (!text || !*text)
        return or_empty(text);
    if (strchr(text, '|'))
        text = gh_report_clean_condition(arena, text);
    return dedup_phrases(arena, text);
}

/* repr(float) for the values the report prints: shortest round-tripping
 * digits, always with a decimal point. Valid for 1e-4 <= |x| < 1e16, which
 * covers odds ratios and 0-100 scores. */
const char *gh_report_pyfloat(gh_arena *arena, double x)
{
    char buf[64];
    int digits = 1;
    for (; digits <= 17; digits++) {
        snprintf(buf, sizeof buf, "%.*g", digits, x);
        if (strtod(buf, NULL) == x)
            break;
    }
    /* Re-express with that many significant digits as fixed notation. */
    int exponent = (x == 0.0) ? 0 : (int)floor(log10(fabs(x)));
    int decimals = digits - 1 - exponent;
    if (decimals < 0)
        decimals = 0;
    snprintf(buf, sizeof buf, "%.*f", decimals, x);
    if (!strchr(buf, '.'))
        strcat(buf, ".0");
    return gh_strdup(arena, buf);
}

/* "{:,}" */
static void grouped(gh_strbuf *sb, size_t value)
{
    char digits[32];
    int n = snprintf(digits, sizeof digits, "%zu", value);
    for (int i = 0; i < n; i++) {
        if (i && (n - i) % 3 == 0)
            gh_sb_putc(sb, ',');
        gh_sb_putc(sb, digits[i]);
    }
}

static const char *eli5_condition(const char *id)
{
    if (!id)
        return "";
    for (size_t i = 0; i < GH_ELI5_CONDITION_COUNT; i++)
        if (strcmp(GH_ELI5_CONDITIONS[i].key, id) == 0)
            return GH_ELI5_CONDITIONS[i].value;
    return "";
}

static bool contains_ci(const char *haystack, const char *needle)
{
    size_t n = strlen(needle);
    if (n == 0)
        return true;
    for (const char *p = haystack; *p; p++) {
        size_t i = 0;
        while (i < n && p[i] &&
               tolower((unsigned char)p[i]) == tolower((unsigned char)needle[i]))
            i++;
        if (i == n)
            return true;
    }
    return false;
}

static const char *ancestry_color(const char *pop)
{
    if (strcmp(pop, "EUR") == 0) return GH_C.blue;
    if (strcmp(pop, "AFR") == 0) return GH_C.amber;
    if (strcmp(pop, "EAS") == 0) return GH_C.red;
    if (strcmp(pop, "SAS") == 0) return GH_C.purple;
    if (strcmp(pop, "AMR") == 0) return GH_C.green;
    return GH_C.slate;
}

static const char *phenotype_color(const char *phenotype)
{
    if (strcmp(phenotype, "poor") == 0) return GH_C.red;
    if (strcmp(phenotype, "intermediate") == 0) return GH_C.amber;
    if (strcmp(phenotype, "normal") == 0) return GH_C.green;
    if (strcmp(phenotype, "rapid") == 0) return GH_C.blue;
    if (strcmp(phenotype, "ultrarapid") == 0) return GH_C.purple;
    return GH_C.slate;   /* "Unknown" and anything else */
}

static const char *risk_band_color(const char *category)
{
    if (strcmp(category, "low") == 0) return GH_C.green;
    if (strcmp(category, "average") == 0) return GH_C.blue;
    if (strcmp(category, "elevated") == 0) return GH_C.amber;
    if (strcmp(category, "high") == 0) return GH_C.red;
    return GH_C.slate;
}

static const char *priority_color(const char *priority)
{
    if (strcmp(priority, "high") == 0) return GH_C.red;
    if (strcmp(priority, "moderate") == 0) return GH_C.amber;
    if (strcmp(priority, "low") == 0) return GH_C.green;
    return "var(--border)";
}

static const char *frequency_color(const char *frequency, gh_arena *arena)
{
    static const struct { const char *key; int color; } table[] = {
        {"weekly", 0}, {"monthly", 1}, {"quarterly", 2},
        {"semi-annually", 3}, {"annually", 3}, {"baseline", 4},
    };
    const char *low = lower_copy(arena, or_empty(frequency));
    for (size_t i = 0; i < sizeof table / sizeof table[0]; i++)
        if (strstr(low, table[i].key)) {
            switch (table[i].color) {
            case 0: return GH_C.red;
            case 1: return GH_C.orange;
            case 2: return GH_C.amber;
            case 3: return GH_C.blue;
            default: return GH_C.green;
            }
        }
    return "var(--accent2)";
}

static const char *urgency_color(const char *urgency)
{
    if (urgency && strcmp(urgency, "soon") == 0) return GH_C.red;
    if (urgency && strcmp(urgency, "routine") == 0) return GH_C.amber;
    return "var(--border)";
}

static void stars_glyphs(gh_strbuf *sb, int stars)
{
    for (int i = 0; i < stars; i++)
        gh_sb_puts(sb, "&#9733;");
    for (int i = stars; i < 4; i++)
        gh_sb_puts(sb, "&#9734;");
}

static void good_news_grid(gh_strbuf *sb, gh_arena *arena,
                           const gh_recommendations *recs, bool with_eli5)
{
    gh_sb_puts(sb, "<div class=\"good-news-grid\">\n");
    const char *seen[256];
    size_t nseen = 0;
    for (size_t i = 0; i < recs->ngood_news; i++) {
        const gh_good_news *g = &recs->good_news[i];
        bool dup = false;
        for (size_t k = 0; k < nseen && !dup; k++)
            dup = strcmp(seen[k], g->gene) == 0;
        if (dup)
            continue;
        if (nseen < 256)
            seen[nseen++] = g->gene;

        const char *desc = strchr(or_empty(g->description), '|')
            ? gh_report_clean_condition(arena, g->description)
            : or_empty(g->description);
        gh_sb_puts(sb, "<div class=\"good-news-card\"><strong>");
        raw(sb, g->gene);
        gh_sb_puts(sb, "</strong>: ");
        raw(sb, desc);
        if (with_eli5) {
            const char *eli5 = gh_report_eli5_gene(g->gene);
            if (*eli5) {
                gh_sb_puts(sb, " <span class=\"eli5-inline\">(");
                raw(sb, eli5);
                gh_sb_puts(sb, ")</span>");
            }
        }
        gh_sb_puts(sb, "</div>\n");
    }
    gh_sb_puts(sb, "</div>");
}

/* Python's `"\n".join(parts)`: every part but the first is preceded by a
 * newline. Sections build with `nl(sb)` before each appended part. */
static void nl(gh_strbuf *sb, bool *first)
{
    if (!*first)
        gh_sb_putc(sb, '\n');
    *first = false;
}

/* ------------------------------------------------------------------ */
/* Charts used by the sections                                         */
/* ------------------------------------------------------------------ */

void gh_report_svg_metabolism_gauge(gh_strbuf *sb, const char *label,
                                    double level, const char *color)
{
    const double cx = 70, cy = 65, r = 50;
    const double start_math = 150, end_math = 30;
    const double range_math = start_math - end_math;

    double ptr_math = start_math - (level / 2) * range_math;
    double ptr_rad = ptr_math * DEG;
    double ptr_r = r - 10;
    double px = cx + ptr_r * cos(ptr_rad);
    double py = cy - ptr_r * sin(ptr_rad);

    double s_rad = start_math * DEG, e_rad = end_math * DEG;
    double sx = cx + r * cos(s_rad), sy = cy - r * sin(s_rad);
    double ex = cx + r * cos(e_rad), ey = cy - r * sin(e_rad);

    /* speed_labels is keyed 0, 1, 2; the half-step levels miss it. */
    const char *speed = level == 0 ? "Slow" : level == 1 ? "Intermediate"
                      : level == 2 ? "Fast" : "Unknown";

    gh_sb_printf(sb,
        "<svg viewBox=\"0 0 140 100\" class=\"gauge\" role=\"img\" "
        "aria-label=\"%s metabolism gauge\">"
        "<title>%s: %s</title>"
        "<path d=\"M %.1f %.1f A %g %g 0 0 1 %.1f %.1f\" "
        "fill=\"none\" stroke=\"var(--border)\" stroke-width=\"8\" stroke-linecap=\"round\"/>"
        "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"6\" fill=\"%s\"/>"
        "<text x=\"15\" y=\"78\" fill=\"currentColor\" font-size=\"8\">Slow</text>"
        "<text x=\"108\" y=\"78\" fill=\"currentColor\" font-size=\"8\">Fast</text>"
        "<text x=\"%g\" y=\"95\" text-anchor=\"middle\" fill=\"currentColor\" "
        "font-size=\"11\" font-weight=\"bold\">%s: %s</text>"
        "</svg>",
        label, label, speed, sx, sy, r, r, ex, ey, px, py, color, cx, label, speed);
}

void gh_report_svg_prs_gauge(gh_strbuf *sb, gh_arena *arena, const char *label,
                             double percentile, const char *category)
{
    const double cx = 80, cy = 70, r = 55;
    const struct { double a, b; const char *color; } zones[] = {
        {0.0, 0.2, GH_C.green}, {0.2, 0.8, GH_C.blue},
        {0.8, 0.95, GH_C.amber}, {0.95, 1.0, GH_C.red},
    };

    gh_sb_printf(sb,
        "<svg viewBox=\"0 0 160 105\" class=\"prs-gauge\" role=\"img\" "
        "aria-label=\"%s PRS gauge\">"
        "<title>%s: %.0fth percentile (%s)</title>",
        label, label, percentile, category);

    for (size_t i = 0; i < 4; i++) {
        double a1 = (180 - zones[i].a * 180) * DEG;
        double a2 = (180 - zones[i].b * 180) * DEG;
        double x1 = cx + r * cos(a1), y1 = cy - r * sin(a1);
        double x2 = cx + r * cos(a2), y2 = cy - r * sin(a2);
        int large = fabs(zones[i].b - zones[i].a) > 0.5 ? 1 : 0;
        gh_sb_printf(sb,
            "<path d=\"M %.1f %.1f A %g %g 0 %d 1 %.1f %.1f\" "
            "fill=\"none\" stroke=\"%s\" stroke-width=\"10\" stroke-linecap=\"butt\"/>",
            x1, y1, r, r, large, x2, y2, zones[i].color);
    }

    double frac = percentile / 100.0;
    if (frac < 0.0) frac = 0.0;
    if (frac > 1.0) frac = 1.0;
    double ptr_angle = (180 - frac * 180) * DEG;
    double ptr_r = r - 18;
    double px = cx + ptr_r * cos(ptr_angle), py = cy - ptr_r * sin(ptr_angle);

    gh_sb_printf(sb,
        "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"6\" fill=\"%s\"/>"
        "<text x=\"%g\" y=\"%g\" text-anchor=\"middle\" fill=\"currentColor\" "
        "font-size=\"16\" font-weight=\"bold\">%.0f%%</text>"
        "<text x=\"%g\" y=\"%g\" text-anchor=\"middle\" fill=\"currentColor\" "
        "font-size=\"10\">%s</text>"
        "<text x=\"%g\" y=\"100\" text-anchor=\"middle\" fill=\"currentColor\" "
        "font-size=\"11\" font-weight=\"bold\">%s</text>"
        "</svg>",
        px, py, risk_band_color(category), cx, cy + 5, percentile,
        cx, cy + 20, py_title(arena, category, false), cx, label);
}

/* The name shortening shared by the heatmap and the gauge row. */
static const char *short_prs_name(gh_arena *arena, const char *name, bool heatmap)
{
    const char *s = replace_all(arena, name, "Age-Related ", "");
    s = replace_all(arena, s, "Macular Degeneration", "AMD");
    if (!heatmap)
        return s;
    s = replace_all(arena, s, "(Bone Mineral Density)", "BMD");
    s = replace_all(arena, s, "(Cholelithiasis)", "");
    s = replace_all(arena, s, "Systemic Lupus Erythematosus", "Lupus (SLE)");
    s = replace_all(arena, s, "Hashimoto's Thyroiditis", "Hashimoto's");
    if (py_len(s) > 16) {
        /* name[:15] + "…", counting code points. */
        const char *p = s;
        for (size_t n = 0; n < 15 && *p; p++)
            if ((*(const unsigned char *)p & 0xC0) != 0x80 &&
                ++n == 15) {
                p++;
                while ((*(const unsigned char *)p & 0xC0) == 0x80)
                    p++;
                break;
            }
        char *cut = gh_alloc(arena, (size_t)(p - s) + 4);
        memcpy(cut, s, (size_t)(p - s));
        strcpy(cut + (p - s), "\xe2\x80\xa6");
        s = cut;
    }
    return s;
}

typedef struct {
    const gh_prs_result *r;
    size_t order;
} ranked_prs;

static int compare_prs_desc(const void *a, const void *b)
{
    const ranked_prs *x = a, *y = b;
    if (x->r->percentile != y->r->percentile)
        return x->r->percentile > y->r->percentile ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

void gh_report_svg_risk_heatmap(gh_strbuf *sb, gh_arena *arena,
                                const gh_prs_result *prs, size_t nprs)
{
    if (!nprs)
        return;
    ranked_prs *items = gh_calloc(arena, nprs, sizeof *items);
    for (size_t i = 0; i < nprs; i++) {
        items[i].r = &prs[i];
        items[i].order = i;
    }
    qsort(items, nprs, sizeof *items, compare_prs_desc);

    size_t cols = nprs < 5 ? nprs : 5;
    size_t rows = (nprs + cols - 1) / cols;
    const int cell_w = 80, cell_h = 60, pad = 4;
    size_t w = cols * (size_t)(cell_w + pad) + (size_t)pad;
    size_t h = rows * (size_t)(cell_h + pad) + (size_t)pad;

    gh_sb_printf(sb,
        "<svg viewBox=\"0 0 %zu %zu\" class=\"chart\" role=\"img\" "
        "aria-label=\"PRS risk heatmap\">"
        "<title>Polygenic Risk Score Heatmap</title>", w, h);

    for (size_t idx = 0; idx < nprs; idx++) {
        const gh_prs_result *r = items[idx].r;
        size_t col = idx % cols, row = idx / cols;
        size_t x = (size_t)pad + col * (size_t)(cell_w + pad);
        size_t y = (size_t)pad + row * (size_t)(cell_h + pad);
        gh_sb_printf(sb,
            "<rect x=\"%zu\" y=\"%zu\" width=\"%d\" height=\"%d\" rx=\"4\" "
            "fill=\"%s\" opacity=\"0.8\"/>"
            "<text x=\"%zu\" y=\"%zu\" text-anchor=\"middle\" "
            "fill=\"#fff\" font-size=\"9\" font-weight=\"bold\">%s</text>"
            "<text x=\"%zu\" y=\"%zu\" text-anchor=\"middle\" "
            "fill=\"#fff\" font-size=\"14\" font-weight=\"bold\">%.0f%%</text>",
            x, y, cell_w, cell_h, risk_band_color(r->risk_category),
            x + (size_t)cell_w / 2, y + 24, short_prs_name(arena, r->name, true),
            x + (size_t)cell_w / 2, y + 42, r->percentile);
    }
    gh_sb_puts(sb, "</svg>");
}

void gh_report_svg_ancestry_donut(gh_strbuf *sb, gh_arena *arena,
                                  const gh_ancestry_result *anc)
{
    const double cx = 110, cy = 110, r = 85, inner_r = 50;
    double angle = -90;

    size_t order[GH_MAX_POPULATIONS];
    for (size_t i = 0; i < GH_POPULATION_COUNT; i++)
        order[i] = i;
    for (size_t i = 0; i + 1 < GH_POPULATION_COUNT; i++)
        for (size_t j = 0; j + 1 < GH_POPULATION_COUNT - i; j++)
            if (anc->proportions[order[j]] < anc->proportions[order[j + 1]]) {
                size_t t = order[j]; order[j] = order[j + 1]; order[j + 1] = t;
            }

    gh_strbuf paths, legend;
    gh_sb_init(&paths);
    gh_sb_init(&legend);

    for (size_t idx = 0; idx < GH_POPULATION_COUNT; idx++) {
        size_t p = order[idx];
        double prop = anc->proportions[p];
        /* Slices below half a percent are skipped, but the legend index
         * still advances, which is why the legend can have gaps. */
        if (prop < 0.005)
            continue;
        const char *color = ancestry_color(GH_POPULATIONS[p]);
        double sweep = prop * 360;
        double a1 = angle * DEG, a2 = (angle + sweep) * DEG;
        double x1 = cx + r * cos(a1), y1 = cy + r * sin(a1);
        double x2 = cx + r * cos(a2), y2 = cy + r * sin(a2);
        double ix1 = cx + inner_r * cos(a1), iy1 = cy + inner_r * sin(a1);
        double ix2 = cx + inner_r * cos(a2), iy2 = cy + inner_r * sin(a2);
        int large = sweep > 180 ? 1 : 0;
        gh_sb_printf(&paths,
            "<path d=\"M %.1f %.1f L %.1f %.1f A %g %g 0 %d 1 %.1f %.1f "
            "L %.1f %.1f A %g %g 0 %d 0 %.1f %.1f Z\" fill=\"%s\" opacity=\"0.85\"/>",
            ix1, iy1, x1, y1, r, r, large, x2, y2, ix2, iy2,
            inner_r, inner_r, large, ix1, iy1, color);
        angle += sweep;

        size_t ly = 10 + idx * 24;
        gh_sb_printf(&legend,
            "<rect x=\"240\" y=\"%zu\" width=\"14\" height=\"14\" rx=\"3\" fill=\"%s\"/>"
            "<text x=\"260\" y=\"%zu\" fill=\"currentColor\" font-size=\"12\">"
            "%s (%.1f%%)</text>",
            ly, color, ly + 12, GH_POPULATION_LABELS[p], prop * 100.0);
    }

    size_t height = 10 + GH_POPULATION_COUNT * 24 + 10;
    if (height < 230)
        height = 230;
    gh_sb_printf(sb,
        "<svg viewBox=\"0 0 420 %zu\" class=\"chart\" role=\"img\" "
        "aria-label=\"Ancestry proportion donut chart\">"
        "<title>Ancestry Proportions</title>", height);
    gh_sb_write(sb, paths.data ? paths.data : "", paths.len);
    gh_sb_printf(sb,
        "<text x=\"%g\" y=\"%g\" text-anchor=\"middle\" "
        "fill=\"currentColor\" font-size=\"14\" font-weight=\"bold\">%s</text>"
        "<text x=\"%g\" y=\"%g\" text-anchor=\"middle\" "
        "fill=\"currentColor\" font-size=\"11\">%s</text>",
        cx, cy - 5, anc->top_ancestry, cx, cy + 15, anc->confidence);
    gh_sb_write(sb, legend.data ? legend.data : "", legend.len);
    gh_sb_puts(sb, "</svg>");
    gh_sb_free(&paths);
    gh_sb_free(&legend);
    (void)arena;
}

/* ------------------------------------------------------------------ */
/* 1. Key findings                                                     */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *color;
    const char *text;
} bullet;

static void bullets_out(gh_strbuf *sb, bool *first, const bullet *items, size_t n)
{
    nl(sb, first);
    gh_sb_puts(sb, "<ul class=\"key-findings\">");
    for (size_t i = 0; i < n; i++) {
        nl(sb, first);
        gh_sb_printf(sb, "<li class=\"kf-%s\">%s</li>", items[i].color, items[i].text);
    }
    nl(sb, first);
    gh_sb_puts(sb, "</ul>");
}

static const char *sb_take(gh_arena *arena, gh_strbuf *sb)
{
    const char *s = gh_strndup(arena, sb->data ? sb->data : "", sb->len);
    gh_sb_free(sb);
    return s;
}

void gh_report_key_findings(gh_strbuf *sb, gh_arena *arena,
                            const gh_analysis *a, const gh_report_scorers *sc)
{
    bool first = true;
    gh_strbuf t;
    bullet items[512];
    size_t n = 0;

    nl(sb, &first);
    gh_sb_puts(sb,
        "<p class=\"eli5\">Your body has a recipe book called DNA. We read yours and "
        "here is everything important we found, in plain language. "
        "Each topic below links to a detailed section later in the report.</p>");

    const gh_recommendations *recs = sc ? sc->recs : NULL;
    const gh_acmg_result *acmg = sc ? sc->acmg : NULL;

    /* ---- Act On These ---- */
    for (size_t i = 0; acmg && i < acmg->nfindings; i++) {
        const gh_cv_finding *f = acmg->findings[i].finding;
        const char *gene = f->gene ? f->gene : "Unknown";
        const char *condition = gh_report_clean_condition(
            arena, f->traits && *f->traits ? f->traits : "Unknown");
        gh_sb_init(&t);
        gh_sb_printf(&t, "<strong>%s</strong>: You carry a variant linked to "
                         "<em>%s</em> that doctors consider medically actionable. "
                         "You should see a genetic counselor.", gene, condition);
        const char *eli5 = gh_report_eli5_gene(gene);
        if (*eli5)
            gh_sb_printf(&t, "<br><span class=\"eli5-inline\">%s</span>", eli5);
        items[n++] = (bullet){"red", sb_take(arena, &t)};
    }
    for (size_t i = 0; recs && i < recs->npriorities; i++) {
        const gh_priority *p = &recs->priorities[i];
        if (strcmp(p->priority, "high") != 0)
            continue;
        gh_sb_init(&t);
        gh_sb_printf(&t, "<strong>%s</strong>: %s", p->title,
                     gh_report_clean_why(arena, p->why));
        const char *eli5 = eli5_condition(p->id);
        if (*eli5)
            gh_sb_printf(&t, "<br><span class=\"eli5-inline\">%s</span>", eli5);
        items[n++] = (bullet){"red", sb_take(arena, &t)};
    }
    if (n) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Act On These <a href=\"#action-plan\">&rarr; details</a></h3>");
        bullets_out(sb, &first, items, n);
    }

    /* ---- Health Risks ---- */
    n = 0;
    const gh_apoe_result *apoe = sc ? sc->apoe : NULL;
    if (apoe && strcmp(apoe->apoe_type, "Unknown") != 0) {
        const char *risk = apoe->risk_level;
        const char *color = (strcmp(risk, "high") == 0 || strcmp(risk, "elevated") == 0)
            ? "red" : strcmp(risk, "moderate") == 0 ? "yellow" : "green";
        gh_sb_init(&t);
        gh_sb_printf(&t, "<strong>APOE %s</strong>: %s Alzheimer's risk (odds ratio: %sx).",
                     apoe->apoe_type, py_title(arena, risk, false),
                     apoe->has_or ? gh_report_pyfloat(arena, apoe->alzheimer_or) : "N/A");
        const char *eli5 = gh_report_eli5_gene("APOE");
        if (*eli5)
            gh_sb_printf(&t, "<br><span class=\"eli5-inline\">%s</span>", eli5);
        items[n++] = (bullet){color, sb_take(arena, &t)};
    }
    if (sc && sc->prs && sc->nprs) {
        for (size_t i = 0; i < sc->nprs; i++) {
            const gh_prs_result *r = &sc->prs[i];
            if (strcmp(r->risk_category, "elevated") != 0 &&
                strcmp(r->risk_category, "high") != 0)
                continue;
            gh_sb_init(&t);
            gh_sb_printf(&t, "<strong>%s</strong>: %.0fth percentile genetic risk "
                             "\xe2\x80\x94 higher than %.0f%% of people.",
                         r->name, r->percentile, r->percentile);
            items[n++] = (bullet){"yellow", sb_take(arena, &t)};
        }
        for (int pass = 0; pass < 2; pass++) {
            const char *band = pass == 0 ? "average" : "low";
            gh_sb_init(&t);
            bool any = false;
            for (size_t i = 0; i < sc->nprs; i++) {
                if (strcmp(sc->prs[i].risk_category, band) != 0)
                    continue;
                if (any)
                    gh_sb_puts(&t, ", ");
                gh_sb_puts(&t, sc->prs[i].name);
                any = true;
            }
            const char *names = sb_take(arena, &t);
            if (!any)
                continue;
            gh_sb_init(&t);
            if (pass == 0)
                gh_sb_printf(&t, "Average genetic risk for: %s.", names);
            else
                gh_sb_printf(&t, "<em>Lower</em> than average risk for: %s.", names);
            items[n++] = (bullet){"green", sb_take(arena, &t)};
        }
    }
    const gh_clinvar_result *cv = sc ? sc->clinvar : NULL;
    if (cv && cv->loaded) {
        size_t path_count = cv->counts[GH_CV_PATHOGENIC];
        size_t lp_count = cv->counts[GH_CV_LIKELY_PATHOGENIC];
        size_t prot_count = cv->counts[GH_CV_PROTECTIVE];
        if (path_count || lp_count) {
            gh_sb_init(&t);
            gh_sb_printf(&t, "%zu pathogenic and %zu likely pathogenic variants "
                             "found in ClinVar scan.", path_count, lp_count);
            items[n++] = (bullet){"yellow", sb_take(arena, &t)};
        }
        if (prot_count) {
            gh_sb_init(&t);
            gh_sb_printf(&t, "%zu protective variants detected \xe2\x80\x94 these "
                             "<em>lower</em> your risk for certain diseases.", prot_count);
            items[n++] = (bullet){"green", sb_take(arena, &t)};
        }
    }
    if (n) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Health Risks <a href=\"#disease-risk\">&rarr; details</a></h3>");
        bullets_out(sb, &first, items, n);
    }

    /* ---- Drugs & Medications ---- */
    n = 0;
    if (sc && sc->stars && sc->nstars) {
        for (size_t i = 0; i < sc->nstars; i++) {
            const gh_star_result *r = &sc->stars[i];
            if (strcmp(r->phenotype, "normal") == 0 || strcmp(r->phenotype, "Unknown") == 0)
                continue;
            gh_sb_init(&t);
            gh_sb_printf(&t, "<strong>%s (%s)</strong>: %s Metabolizer \xe2\x80\x94 "
                             "some drugs need dose adjustment.",
                         r->gene, r->diplotype, py_title(arena, r->phenotype, true));
            const char *eli5 = gh_report_eli5_gene(r->gene);
            if (*eli5)
                gh_sb_printf(&t, "<br><span class=\"eli5-inline\">%s</span>", eli5);
            items[n++] = (bullet){"yellow", sb_take(arena, &t)};
        }
        gh_sb_init(&t);
        bool any = false;
        for (size_t i = 0; i < sc->nstars; i++) {
            if (strcmp(sc->stars[i].phenotype, "normal") != 0)
                continue;
            if (any)
                gh_sb_puts(&t, ", ");
            gh_sb_puts(&t, sc->stars[i].gene);
            any = true;
        }
        const char *names = sb_take(arena, &t);
        if (any) {
            gh_sb_init(&t);
            gh_sb_printf(&t, "Normal metabolism for: %s \xe2\x80\x94 standard drug "
                             "doses should work.", names);
            items[n++] = (bullet){"green", sb_take(arena, &t)};
        }
    }
    if (n) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Drugs &amp; Medications <a href=\"#drug-guide\">&rarr; details</a></h3>");
        bullets_out(sb, &first, items, n);
    }

    /* ---- Nutrition & Supplements ---- */
    n = 0;
    static const char *const nutrition_groups[] = {
        "methylation", "iron", "caffeine", "nutrition", "vitamin"};
    for (size_t i = 0; recs && i < recs->npriorities; i++) {
        const gh_priority *p = &recs->priorities[i];
        if (strcmp(p->priority, "high") != 0 && strcmp(p->priority, "moderate") != 0)
            continue;
        /* "not in urgent titles": every high-priority title is urgent. */
        bool urgent = false;
        for (size_t k = 0; k < recs->npriorities && !urgent; k++)
            urgent = strcmp(recs->priorities[k].priority, "high") == 0 &&
                     strcmp(recs->priorities[k].title, p->title) == 0;
        if (urgent)
            continue;
        const char *low = lower_copy(arena, or_empty(p->id));
        bool hit = false;
        for (size_t g = 0; g < 5 && !hit; g++)
            hit = strstr(low, nutrition_groups[g]) != NULL;
        if (!hit)
            continue;
        gh_sb_init(&t);
        gh_sb_printf(&t, "<strong>%s</strong>: %s", p->title,
                     gh_report_clean_why(arena, p->why));
        const char *eli5 = eli5_condition(p->id);
        if (*eli5)
            gh_sb_printf(&t, "<br><span class=\"eli5-inline\">%s</span>", eli5);
        items[n++] = (bullet){"yellow", sb_take(arena, &t)};
    }
    const gh_nutrigenomics_result *nut = sc ? sc->nutrition : NULL;
    if (nut) {
        /* supplement_priorities: high/moderate needs in profile order, then
         * a stable sort by priority -- see supplement_priorities() below. */
        const gh_nutrient_need *supps[64];
        size_t nsupps = gh_report_supplement_priorities(nut, supps, 64);
        if (nsupps) {
            gh_sb_init(&t);
            for (size_t i = 0; i < nsupps && i < 4; i++) {
                if (i)
                    gh_sb_puts(&t, ", ");
                gh_sb_printf(&t, "%s (%s)", supps[i]->profile->name,
                             supps[i]->profile->supplement_form);
            }
            const char *top = sb_take(arena, &t);
            gh_sb_init(&t);
            gh_sb_printf(&t, "<strong>Top supplement priorities</strong>: %s.", top);
            items[n++] = (bullet){"yellow", sb_take(arena, &t)};
        }
    }
    if (n) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Nutrition &amp; Supplements <a href=\"#nutrigenomics\">&rarr; details</a></h3>");
        bullets_out(sb, &first, items, n);
    }

    /* ---- Mental Health ---- */
    const gh_mental_result *mh = sc ? sc->mental : NULL;
    if (mh && mh->ndomains) {
        n = 0;
        bool elevated = false;
        for (int pass = 0; pass < 2; pass++) {
            const char *band = pass == 0 ? "elevated" : "low";
            gh_sb_init(&t);
            bool any = false;
            for (size_t i = 0; i < mh->ndomains; i++) {
                if (strcmp(mh->domains[i].risk_level, band) != 0)
                    continue;
                if (any)
                    gh_sb_puts(&t, ", ");
                gh_sb_puts(&t, py_title(arena, mh->domains[i].name, true));
                any = true;
            }
            const char *names = sb_take(arena, &t);
            if (!any)
                continue;
            if (pass == 0)
                elevated = true;
            gh_sb_init(&t);
            if (pass == 0)
                gh_sb_printf(&t, "Elevated genetic susceptibility for: %s. Lifestyle "
                                 "and support make a huge difference.", names);
            else
                gh_sb_printf(&t, "Low genetic susceptibility for: %s.", names);
            items[n++] = (bullet){pass == 0 ? "yellow" : "green", sb_take(arena, &t)};
        }
        if (mh->summary && *mh->summary)
            items[n++] = (bullet){elevated ? "yellow" : "green", mh->summary};

        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Mental Health <a href=\"#mental-health\">&rarr; details</a></h3>");
        nl(sb, &first);
        gh_sb_puts(sb,
            "<p class=\"eli5-inline\">Your genes influence how your brain handles "
            "stress, mood, and anxiety \xe2\x80\x94 but environment, relationships, and exercise "
            "matter just as much. Genes are the cards; lifestyle is how you play them.</p>");
        bullets_out(sb, &first, items, n);
    }

    /* ---- Your Body ---- */
    n = 0;
    if (sc && sc->traits && sc->ntraits) {
        gh_sb_init(&t);
        gh_sb_puts(&t, "<strong>Traits</strong>: ");
        for (size_t i = 0; i < sc->ntraits; i++) {
            if (i)
                gh_sb_puts(&t, " &middot; ");
            gh_sb_printf(&t, "%s: %s", py_title(arena, sc->traits[i].key, true),
                         sc->traits[i].prediction);
        }
        items[n++] = (bullet){"green", sb_take(arena, &t)};
    }
    if (sc && sc->blood && strcmp(sc->blood->blood_type, "Unknown") != 0) {
        gh_sb_init(&t);
        gh_sb_printf(&t, "<strong>Blood type</strong>: %s (%s confidence)",
                     sc->blood->blood_type, sc->blood->confidence);
        items[n++] = (bullet){"green", sb_take(arena, &t)};
    }
    if (sc && sc->sleep) {
        gh_sb_init(&t);
        gh_sb_printf(&t, "<strong>Chronotype</strong>: %s (score %s). Optimal sleep: %s.",
                     sc->sleep->chronotype,
                     gh_report_pyfloat(arena, sc->sleep->chronotype_score),
                     or_empty(sc->sleep->optimal_sleep_window));
        items[n++] = (bullet){"green", sb_take(arena, &t)};
    }
    if (sc && sc->longevity) {
        double score = sc->longevity->longevity_score;
        const char *color = score >= 60 ? "green" : score >= 40 ? "yellow" : "red";
        gh_sb_init(&t);
        gh_sb_printf(&t, "<strong>Longevity Score</strong>: %s/100. %s",
                     gh_report_pyfloat(arena, score), or_empty(sc->longevity->summary));
        items[n++] = (bullet){color, sb_take(arena, &t)};
    }
    for (size_t i = 0; a && i < a->nfindings; i++) {
        if (!a->findings[i].gene || strcmp(a->findings[i].gene, "ACTN3") != 0)
            continue;
        const char *eli5 = gh_report_eli5_gene("ACTN3");
        gh_sb_init(&t);
        gh_sb_printf(&t, "<strong>Athletic profile (ACTN3)</strong>: %s.%s%s",
                     py_title(arena, a->findings[i].status, true),
                     *eli5 ? " " : "", eli5);
        items[n++] = (bullet){"green", sb_take(arena, &t)};
        break;
    }
    if (sc && sc->ancestry && sc->ancestry->top_ancestry && *sc->ancestry->top_ancestry) {
        gh_sb_init(&t);
        gh_sb_printf(&t, "<strong>Ancestry</strong>: Primarily %s (%s confidence, %zu markers).",
                     sc->ancestry->top_ancestry, or_empty(sc->ancestry->confidence),
                     sc->ancestry->markers_found);
        items[n++] = (bullet){"green", sb_take(arena, &t)};
    }
    if (n) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Your Body <a href=\"#body-profile\">&rarr; details</a></h3>");
        bullets_out(sb, &first, items, n);
    }

    /* ---- Good News ---- */
    if (recs && recs->ngood_news) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Good News</h3>");
        nl(sb, &first);
        good_news_grid(sb, arena, recs, true);
    }

    nl(sb, &first);
    gh_sb_puts(sb,
        "<div class=\"doctor-callout\">Show this report to your doctor! "
        "They can help you use this information wisely. "
        "Each section below has the detailed data behind these findings.</div>");
}

/* ------------------------------------------------------------------ */
/* 2. Action plan                                                      */
/* ------------------------------------------------------------------ */

/* Python builds {gene: insight} over the list, so a gene keeps its first
 * position but its last entry. */
static size_t insights_by_gene(const gh_recommendations *recs,
                               const gh_clinical_insight **out, size_t cap)
{
    size_t n = 0;
    for (size_t i = 0; i < recs->ninsights; i++) {
        const char *gene = recs->insights[i].gene;
        size_t k = 0;
        for (; k < n; k++)
            if (strcmp(out[k]->gene, gene) == 0)
                break;
        if (k < n)
            out[k] = &recs->insights[i];
        else if (n < cap)
            out[n++] = &recs->insights[i];
    }
    return n;
}

static void monitoring_table(gh_strbuf *sb, bool *first, gh_arena *arena,
                             const gh_monitoring_item *items, size_t n,
                             bool colored)
{
    nl(sb, first);
    gh_sb_puts(sb, "<table><tr><th>Test</th><th>Frequency</th><th>Reason</th></tr>");
    for (size_t i = 0; i < n; i++) {
        nl(sb, first);
        if (colored)
            gh_sb_printf(sb,
                "<tr><td><strong>%s</strong></td>"
                "<td><span class=\"mag-badge\" style=\"background:%s;color:#fff\">%s</span></td>"
                "<td>%s</td></tr>",
                items[i].test, frequency_color(items[i].frequency, arena),
                items[i].frequency, items[i].reason);
        else
            gh_sb_printf(sb, "<tr><td>%s</td><td>%s</td><td>%s</td></tr>",
                         items[i].test, items[i].frequency, items[i].reason);
    }
    nl(sb, first);
    gh_sb_puts(sb, "</table>");
}

void gh_report_action_plan(gh_strbuf *sb, gh_arena *arena,
                           const gh_report_scorers *sc)
{
    const gh_recommendations *recs = sc ? sc->recs : NULL;
    if (!recs) {
        gh_sb_puts(sb, "<p>No personalized recommendations available.</p>");
        return;
    }
    bool first = true;

    nl(sb, &first);
    gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                   "Actions grouped by urgency, with genetic rationale for each.</p>");

    const gh_clinical_insight *by_gene[128];
    size_t nby_gene = insights_by_gene(recs, by_gene, 128);

    if (!recs->npriorities) {
        nl(sb, &first);
        gh_sb_puts(sb, "<p>No convergent risk patterns detected.</p>");
    }
    for (size_t i = 0; i < recs->npriorities; i++) {
        const gh_priority *p = &recs->priorities[i];
        const char *color = priority_color(p->priority);
        const char *eli5 = eli5_condition(p->id);

        nl(sb, &first);
        gh_sb_printf(sb,
            "<details open class=\"rec-card\" style=\"border-left:4px solid %s;\">"
            "<summary><span class=\"mag-badge\" style=\"background:%s;color:#fff\">"
            "%s</span> <strong>%s</strong></summary>",
            color, color, py_upper(arena, p->priority), p->title);
        if (*eli5) {
            nl(sb, &first);
            gh_sb_printf(sb, "<p class=\"eli5\">%s</p>", eli5);
        }
        nl(sb, &first);
        gh_sb_printf(sb, "<p><strong>Why:</strong> %s</p>",
                     gh_report_clean_why(arena, p->why));
        nl(sb, &first);
        gh_sb_puts(sb, "<p><strong>Actions:</strong></p><ol>");
        for (size_t k = 0; k < p->nactions; k++) {
            nl(sb, &first);
            gh_sb_printf(sb, "<li>%s</li>", p->actions[k]);
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</ol>");

        if (p->nclinical_actions) {
            nl(sb, &first);
            gh_sb_puts(sb, "<p><strong>Clinical actions:</strong></p><ul>");
            for (size_t k = 0; k < p->nclinical_actions; k++) {
                nl(sb, &first);
                gh_sb_printf(sb, "<li>%s</li>", p->clinical_actions[k]);
            }
            nl(sb, &first);
            gh_sb_puts(sb, "</ul>");
        }
        if (p->doctor_note && *p->doctor_note) {
            nl(sb, &first);
            gh_sb_printf(sb,
                "<div class=\"doctor-callout\" style=\"text-align:left;font-weight:normal\">"
                "<strong>Doctor note:</strong> %s</div>", p->doctor_note);
        }
        if (p->nmonitoring)
            monitoring_table(sb, &first, arena, p->monitoring, p->nmonitoring, false);

        /* Clinical context whose gene appears in the title or the why. */
        bool any = false;
        for (size_t k = 0; k < nby_gene; k++) {
            const gh_clinical_insight *ci = by_gene[k];
            if (!contains_ci(p->title, ci->gene) &&
                !contains_ci(or_empty(p->why), ci->gene))
                continue;
            if (!any) {
                nl(sb, &first);
                gh_sb_puts(sb, "<details style=\"margin-top:.5em\">"
                               "<summary style=\"font-size:.85em;color:var(--accent2)\">"
                               "Clinical Context</summary>");
                any = true;
            }
            nl(sb, &first);
            gh_sb_printf(sb, "<p><strong>%s:</strong> %s</p>", ci->gene,
                         ci->context ? or_empty(ci->context->mechanism) : "");
            if (ci->context && ci->context->nactions) {
                nl(sb, &first);
                gh_sb_puts(sb, "<ul>");
                for (size_t x = 0; x < ci->context->nactions; x++) {
                    nl(sb, &first);
                    gh_sb_printf(sb, "<li>%s</li>", ci->context->actions[x]);
                }
                nl(sb, &first);
                gh_sb_puts(sb, "</ul>");
            }
        }
        if (any) {
            nl(sb, &first);
            gh_sb_puts(sb, "</details>");
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</details>");
    }

    if (recs->nmonitoring) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Monitoring Schedule</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "Print this and bring it to your next doctor visit.</p>");
        monitoring_table(sb, &first, arena, recs->monitoring_schedule,
                         recs->nmonitoring, true);
    }

    if (recs->nreferrals) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Doctors to See</h3>");
        for (size_t i = 0; i < recs->nreferrals; i++) {
            const gh_referral *ref = &recs->referrals[i];
            const char *color = urgency_color(ref->urgency);
            nl(sb, &first);
            gh_sb_printf(sb,
                "<div class=\"finding-card\" style=\"border-left-color:%s\">"
                "<strong>%s</strong>: %s <span class=\"badge\" style=\"background:%s;color:#fff\">"
                "%s</span></div>",
                color, ref->specialist, gh_report_clean_why(arena, ref->reason),
                color, ref->urgency ? ref->urgency : "routine");
        }
    }

    const gh_preventive_result *pc = sc->preventive;
    if (pc && pc->ntimeline) {
        bool any = false;
        for (size_t i = 0; i < pc->ntimeline && !any; i++)
            any = strcmp(pc->timeline[i].priority, "standard") != 0 &&
                  strcmp(pc->timeline[i].priority, "ongoing") != 0;
        if (any) {
            nl(sb, &first);
            gh_sb_puts(sb, "<h3>Preventive Screening Timeline</h3>");
            nl(sb, &first);
            gh_sb_printf(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">%s</p>",
                         or_empty(pc->summary));
            nl(sb, &first);
            gh_sb_puts(sb, "<table><tr><th>Test</th><th>Start Age</th>"
                           "<th>Frequency</th><th>Why</th></tr>");
            for (size_t i = 0; i < pc->ntimeline; i++) {
                const gh_screening *s = &pc->timeline[i];
                if (strcmp(s->priority, "standard") == 0 || strcmp(s->priority, "ongoing") == 0)
                    continue;
                const char *color =
                    strcmp(s->priority, "urgent") == 0 ? GH_C.red :
                    strcmp(s->priority, "high") == 0 ? GH_C.orange :
                    strcmp(s->priority, "elevated") == 0 ? GH_C.amber : "var(--border)";
                nl(sb, &first);
                gh_sb_puts(sb, "<tr><td><strong>");
                esc(sb, s->test);
                gh_sb_printf(sb, "</strong></td><td>%d</td><td>", s->start_age);
                esc(sb, s->frequency);
                gh_sb_puts(sb, "</td><td>");
                esc(sb, s->reason);
                if (s->genetic_basis && *s->genetic_basis) {
                    gh_sb_printf(sb, " <span class=\"badge\" style=\"background:%s;color:#fff\">", color);
                    esc(sb, s->genetic_basis);
                    gh_sb_puts(sb, "</span>");
                }
                gh_sb_puts(sb, "</td></tr>");
            }
            nl(sb, &first);
            gh_sb_puts(sb, "</table>");
        }
    }

    if (recs->ngood_news) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Good News</h3>");
        nl(sb, &first);
        good_news_grid(sb, arena, recs, false);
    }
}

/* ------------------------------------------------------------------ */
/* 3. Drug guide                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    const gh_dose_rec *rec;
    size_t order;
} ranked_dose;

static int dose_category_rank(const char *category)
{
    static const char *const order[] = {
        "chemotherapy", "anticoagulant", "antiplatelet", "analgesic",
        "immunosuppressant", "statin", "antibiotic", "dietary"};
    for (int i = 0; i < 8; i++)
        if (strcmp(category, order[i]) == 0)
            return i;
    return 99;
}

static int compare_dose(const void *a, const void *b)
{
    const ranked_dose *x = a, *y = b;
    int rx = dose_category_rank(x->rec->drug->category);
    int ry = dose_category_rank(y->rec->drug->category);
    if (rx != ry)
        return rx < ry ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

void gh_report_drug_guide(gh_strbuf *sb, gh_arena *arena, const gh_analysis *a,
                          const gh_report_scorers *sc)
{
    const gh_star_result *stars = sc ? sc->stars : NULL;
    size_t nstars = sc ? sc->nstars : 0;
    if ((!stars || !nstars) && (!a || !a->ndrug_findings)) {
        gh_sb_puts(sb, "<p>No drug-gene interaction data available.</p>");
        return;
    }
    bool first = true;

    nl(sb, &first);
    gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                   "How your genes affect drug processing. "
                   "Show this section to every prescribing physician.</p>");

    if (stars && nstars) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Your Drug-Processing Enzymes</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "Enzyme speed determines how fast your liver clears each drug.</p>");
        nl(sb, &first);
        gh_sb_puts(sb, "<div class=\"gauge-row\" style=\"justify-content:center\">");
        for (size_t i = 0; i < nstars; i++) {
            const char *ph = stars[i].phenotype;
            double level = strcmp(ph, "poor") == 0 ? 0 : strcmp(ph, "intermediate") == 0 ? 0.5
                         : strcmp(ph, "normal") == 0 ? 1 : strcmp(ph, "rapid") == 0 ? 1.5
                         : strcmp(ph, "ultrarapid") == 0 ? 2 : 1;
            nl(sb, &first);
            gh_report_svg_metabolism_gauge(sb, stars[i].gene, level, phenotype_color(ph));
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</div>");

        nl(sb, &first);
        gh_sb_puts(sb, "<table><tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th>"
                       "<th>SNPs</th><th>What This Means</th></tr>");
        for (size_t i = 0; i < nstars; i++) {
            const gh_star_result *r = &stars[i];
            const char *eli5 = gh_report_eli5_gene(r->gene);
            /* clinical_note[:120] -- code points, not bytes. */
            const char *note = eli5;
            if (!*eli5) {
                const char *cn = or_empty(r->clinical_note);
                const char *p = cn;
                for (size_t n = 0; *p && n < 120; p++)
                    if ((*(const unsigned char *)p & 0xC0) != 0x80 && ++n == 120) {
                        p++;
                        while ((*(const unsigned char *)p & 0xC0) == 0x80)
                            p++;
                        break;
                    }
                note = gh_strndup(arena, cn, (size_t)(p - cn));
            }
            nl(sb, &first);
            gh_sb_printf(sb,
                "<tr><td><strong>%s</strong></td><td><code>%s</code></td><td>%s</td>"
                "<td>%zu/%zu</td><td style=\"font-size:.85em\">%s</td></tr>",
                r->gene, r->diplotype, py_title(arena, r->phenotype, true),
                r->snps_found, r->snps_total, note);
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</table>");
    }

    if (a && a->ndrug_findings) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Drug-Gene Interactions (PharmGKB)</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<table><tr><th>Gene</th><th>RSID</th><th>Level</th>"
                       "<th>Drugs</th><th>Genotype</th></tr>");
        for (size_t i = 0; i < a->ndrug_findings; i++) {
            const gh_drug_finding *d = &a->drug_findings[i];
            nl(sb, &first);
            gh_sb_puts(sb, "<tr><td><strong>");
            esc(sb, d->gene);
            gh_sb_puts(sb, "</strong></td><td><code>");
            esc(sb, d->rsid);
            gh_sb_puts(sb, "</code> ");
            gh_report_db_links(sb, d->rsid);
            gh_sb_puts(sb, "</td><td>");
            esc(sb, d->level);
            gh_sb_puts(sb, "</td><td>");
            esc(sb, d->drugs);
            gh_sb_puts(sb, "</td><td><code>");
            esc(sb, d->genotype);
            gh_sb_puts(sb, "</code></td></tr>");
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</table>");
    }

    const gh_polypharmacy_result *pp = sc ? sc->polypharmacy : NULL;
    if (pp && pp->nwarnings) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Drug Combination Warnings</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "Drug combinations that may interfere based on your enzyme profile.</p>");
        for (size_t i = 0; i < pp->nwarnings; i++) {
            const gh_poly_warning *w = &pp->warnings[i];
            const char *color = priority_color(w->rule->severity);
            nl(sb, &first);
            gh_sb_printf(sb,
                "<details open class=\"rec-card\" style=\"border-left:4px solid %s\">"
                "<summary><span class=\"mag-badge\" style=\"background:%s;color:#fff\">"
                "%s</span> <strong>", color, color, py_upper(arena, w->rule->severity));
            esc(sb, w->rule->name);
            gh_sb_puts(sb, "</strong></summary>");
            if (w->rule->ngenes) {
                nl(sb, &first);
                gh_sb_puts(sb, "<p><strong>Your genotype:</strong> ");
                gh_strbuf g;
                gh_sb_init(&g);
                for (size_t k = 0; k < w->rule->ngenes; k++)
                    gh_sb_printf(&g, "%s%s: %s", k ? ", " : "",
                                 w->rule->genes[k].gene, w->matched[k]);
                esc(sb, sb_take(arena, &g));
                gh_sb_puts(sb, "</p>");
            }
            nl(sb, &first);
            gh_sb_puts(sb, "<p><strong>Drugs affected:</strong> ");
            gh_strbuf d;
            gh_sb_init(&d);
            for (size_t k = 0; k < w->rule->ndrugs; k++)
                gh_sb_printf(&d, "%s%s", k ? ", " : "", w->rule->drugs_affected[k]);
            esc(sb, sb_take(arena, &d));
            gh_sb_puts(sb, "</p>");
            nl(sb, &first);
            gh_sb_puts(sb, "<p>");
            esc(sb, w->rule->warning);
            gh_sb_puts(sb, "</p>");
            nl(sb, &first);
            gh_sb_puts(sb, "<p><strong>Clinical action:</strong> ");
            esc(sb, w->rule->action);
            gh_sb_puts(sb, "</p>");
            nl(sb, &first);
            gh_sb_puts(sb, "</details>");
        }
    }

    const gh_dosing_result *dd = sc ? sc->dosing : NULL;
    if (dd && dd->nrecommendations) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Drug-Specific Dosing Guide</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "CPIC/DPWG guideline-based dosing adjustments for your genotype.</p>");
        ranked_dose *ranked = gh_calloc(arena, dd->nrecommendations, sizeof *ranked);
        for (size_t i = 0; i < dd->nrecommendations; i++) {
            ranked[i].rec = &dd->recommendations[i];
            ranked[i].order = i;
        }
        qsort(ranked, dd->nrecommendations, sizeof *ranked, compare_dose);
        for (size_t i = 0; i < dd->nrecommendations; i++) {
            const gh_dose_rec *rec = ranked[i].rec;
            const char *color = rec->critical ? GH_C.red : GH_C.amber;
            const char *badge = rec->critical ? "CRITICAL"
                                              : py_upper(arena, rec->drug->category);
            nl(sb, &first);
            gh_sb_printf(sb,
                "<details open class=\"rec-card\" style=\"border-left:4px solid %s\">"
                "<summary><span class=\"mag-badge\" style=\"background:%s;color:#fff\">"
                "%s</span> <strong>", color, color, badge);
            esc(sb, rec->drug->display);
            gh_sb_puts(sb, "</strong></summary><p><strong>Genes:</strong> ");
            for (size_t k = 0; k < rec->drug->ngenes; k++)
                gh_sb_printf(sb, "%s%s", k ? ", " : "", rec->drug->genes[k]);
            gh_sb_puts(sb, "</p><p>");
            esc(sb, rec->rule->action);
            gh_sb_puts(sb, "</p><p><strong>Dose guidance:</strong> ");
            esc(sb, rec->rule->dose_guidance);
            gh_sb_puts(sb, "</p><p class=\"paper-refs\">Source: ");
            esc(sb, rec->drug->source);
            gh_sb_puts(sb, "</p></details>");
        }
    }

    nl(sb, &first);
    gh_sb_puts(sb, "<h3>Drug Interaction Checker</h3>");
    nl(sb, &first);
    gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                   "Type a drug name to see if your genotype affects it.</p>");
    nl(sb, &first);
    gh_sb_puts(sb,
        "<input type=\"text\" id=\"drug-checker\" placeholder=\"Search for a drug (e.g., warfarin, codeine, simvastatin)...\" "
        "style=\"width:100%;padding:.5em .75em;border:2px solid var(--border);border-radius:6px;"
        "background:var(--bg);color:var(--fg);font-size:.92em;font-family:var(--body-font);margin:.5em 0\">"
        "<div id=\"drug-checker-results\"></div>");
    nl(sb, &first);
    gh_sb_puts(sb, "<div class=\"doctor-callout\">Share this entire drug section with "
                   "every prescribing physician. Print it or save as PDF.</div>");
}

/* ------------------------------------------------------------------ */
/* 4. Disease risk                                                     */
/* ------------------------------------------------------------------ */

static void variant_table(gh_strbuf *sb, bool *first, gh_arena *arena,
                          const gh_cv_finding *variants, size_t n,
                          const char *label, const char *color)
{
    /* Python appends "" for an empty table, which still contributes its
     * "\n" to the join at this position. */
    nl(sb, first);
    if (!n)
        return;
    gh_sb_printf(sb, "<h3>%s <span class=\"badge\" style=\"background:%s\">%zu</span></h3>",
                 label, color, n);
    nl(sb, first);
    gh_sb_puts(sb, "<table><tr><th>Gene</th><th>Condition</th><th>Genotype</th>"
                   "<th>Stars</th><th>Zygosity</th></tr>");
    for (size_t i = 0; i < n && i < 50; i++) {
        const gh_cv_finding *v = &variants[i];
        nl(sb, first);
        gh_sb_puts(sb, "<tr><td><strong>");
        esc(sb, v->gene ? v->gene : "Unknown");
        gh_sb_puts(sb, "</strong></td>"
                       "<td style=\"max-width:400px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap\" "
                       "title=\"");
        esc(sb, v->traits);
        gh_sb_puts(sb, "\">");
        esc(sb, gh_report_clean_condition(arena, v->traits && *v->traits ? v->traits : "Unknown"));
        gh_sb_puts(sb, "</td><td><code>");
        esc(sb, v->user_genotype);
        gh_sb_puts(sb, "</code></td><td>");
        stars_glyphs(sb, v->gold_stars);
        gh_sb_puts(sb, "</td><td>");
        esc(sb, py_title(arena, v->zygosity, true));
        gh_sb_puts(sb, "</td></tr>");
    }
    nl(sb, first);
    gh_sb_puts(sb, "</table>");
    if (n > 50) {
        nl(sb, first);
        gh_sb_printf(sb, "<p style=\"font-size:.85em;color:var(--accent2)\">"
                         "Showing 50 of %zu variants.</p>", n);
    }
}

void gh_report_disease_risk(gh_strbuf *sb, gh_arena *arena,
                            const gh_report_scorers *sc)
{
    const gh_prs_result *prs = sc ? sc->prs : NULL;
    size_t nprs = sc ? sc->nprs : 0;
    const gh_clinvar_result *cv = sc && sc->clinvar && sc->clinvar->loaded ? sc->clinvar : NULL;
    const gh_acmg_result *acmg = sc ? sc->acmg : NULL;

    if ((!prs || !nprs) && !cv && !acmg) {
        gh_sb_puts(sb, "<p>No disease risk data available.</p>");
        return;
    }
    bool first = true;

    nl(sb, &first);
    gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                   "Genetic risk estimates from polygenic scores and ClinVar variants. "
                   "Lifestyle and environment also affect risk significantly.</p>");

    if (prs && nprs) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Risk Heatmap</h3>");
        nl(sb, &first);
        gh_report_svg_risk_heatmap(sb, arena, prs, nprs);

        bool any_na = false;
        for (size_t i = 0; i < nprs; i++)
            if (!prs[i].ancestry_applicable)
                any_na = true;
        if (any_na) {
            nl(sb, &first);
            gh_sb_puts(sb,
                "<div class=\"doctor-callout\" style=\"border-color:var(--warn)\">"
                "PRS models are calibrated on European-ancestry populations. "
                "Your ancestry profile is substantially non-European \xe2\x80\x94 interpret with caution."
                "</div>");
        }
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Polygenic Risk Scores</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "Your combined variant score vs. the general population.</p>");
        nl(sb, &first);
        gh_sb_puts(sb, "<div class=\"gauge-row\" style=\"justify-content:center\">");
        for (size_t i = 0; i < nprs; i++) {
            nl(sb, &first);
            gh_report_svg_prs_gauge(sb, arena, short_prs_name(arena, prs[i].name, false),
                                    prs[i].percentile, prs[i].risk_category);
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</div>");

        nl(sb, &first);
        gh_sb_puts(sb, "<table class=\"sortable\"><tr><th>Condition</th><th>Percentile</th>"
                       "<th>Category</th><th>SNPs</th><th>Reference</th></tr>");
        for (size_t i = 0; i < nprs; i++) {
            const gh_prs_result *r = &prs[i];
            nl(sb, &first);
            gh_sb_printf(sb,
                "<tr><td><strong>%s</strong></td><td>%.0fth</td><td>%s</td>"
                "<td>%zu/%zu</td><td style=\"font-size:.8em\">%s</td></tr>",
                r->name, r->percentile, py_title(arena, r->risk_category, false),
                r->snps_found, r->snps_total, r->reference);
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</table>");

        bool any_el = false;
        for (size_t i = 0; i < nprs; i++)
            if (strcmp(prs[i].risk_category, "elevated") == 0 ||
                strcmp(prs[i].risk_category, "high") == 0)
                any_el = true;
        if (any_el) {
            nl(sb, &first);
            gh_sb_puts(sb, "<h3>Elevated Risk Details</h3>");
            for (size_t i = 0; i < nprs; i++) {
                const gh_prs_result *r = &prs[i];
                if (strcmp(r->risk_category, "elevated") != 0 &&
                    strcmp(r->risk_category, "high") != 0)
                    continue;
                nl(sb, &first);
                gh_sb_printf(sb, "<details open><summary><strong>%s</strong> \xe2\x80\x94 "
                                 "%.0fth percentile (%s)</summary>",
                             r->name, r->percentile, r->risk_category);
                if (r->ncontributing) {
                    nl(sb, &first);
                    gh_sb_puts(sb, "<table><tr><th>Gene</th><th>rsID</th><th>Copies</th>"
                                   "<th>Effect</th></tr>");
                    for (size_t k = 0; k < r->ncontributing && k < 5; k++) {
                        const gh_prs_contribution *c = &r->contributing[k];
                        nl(sb, &first);
                        gh_sb_printf(sb, "<tr><td>%s</td><td><code>%s</code> ", c->gene, c->rsid);
                        gh_report_db_links(sb, c->rsid);
                        gh_sb_printf(sb, "</td><td>%d</td><td>%.3f</td></tr>",
                                     c->copies, c->contribution);
                    }
                    nl(sb, &first);
                    gh_sb_puts(sb, "</table>");
                }
                nl(sb, &first);
                gh_sb_puts(sb, "</details>");
            }
        }
    }

    if (acmg) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Medically Actionable Genes (ACMG SF v3.2)</h3>");
        nl(sb, &first);
        gh_sb_printf(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                         "Screened %zu medically actionable genes (ACMG SF v3.2).</p>",
                     acmg->genes_screened);
        if (!acmg->nfindings) {
            nl(sb, &first);
            gh_sb_puts(sb, "<div class=\"doctor-callout\" style=\"border-color:var(--green)\">"
                           "No pathogenic/likely pathogenic variants found in ACMG genes.</div>");
        } else {
            nl(sb, &first);
            gh_sb_printf(sb, "<div class=\"doctor-callout\" style=\"border-color:var(--warn)\">"
                             "<strong>%zu variant(s)</strong> found in %zu ACMG gene(s). "
                             "Genetic counseling recommended.</div>",
                         acmg->nfindings, acmg->genes_with_variants);
            nl(sb, &first);
            gh_sb_puts(sb, "<table><tr><th>Gene</th><th>Condition</th>"
                           "<th>Genotype</th><th>Stars</th><th>Actionability</th></tr>");
            for (size_t i = 0; i < acmg->nfindings; i++) {
                const gh_cv_finding *f = acmg->findings[i].finding;
                const char *gene = f->gene ? f->gene : "Unknown";
                const char *eli5 = gh_report_eli5_gene(gene);
                nl(sb, &first);
                gh_sb_printf(sb, "<tr><td><strong>%s</strong>", gene);
                if (*eli5)
                    gh_sb_printf(sb, "<br><span class=eli5-inline>%s</span>", eli5);
                gh_sb_printf(sb, "</td><td>%s</td><td><code>%s</code></td><td>",
                             gh_report_clean_condition(arena, f->traits && *f->traits ? f->traits : "Unknown"),
                             or_empty(f->user_genotype));
                stars_glyphs(sb, f->gold_stars);
                gh_sb_printf(sb, "</td><td style=\"font-size:.85em\">%s</td></tr>",
                             or_empty(acmg->findings[i].actionability));
            }
            nl(sb, &first);
            gh_sb_puts(sb, "</table>");
        }
    }

    if (cv) {
        variant_table(sb, &first, arena, cv->by_category[GH_CV_PATHOGENIC],
                      cv->counts[GH_CV_PATHOGENIC], "Pathogenic Variants", GH_C.red);
        variant_table(sb, &first, arena, cv->by_category[GH_CV_LIKELY_PATHOGENIC],
                      cv->counts[GH_CV_LIKELY_PATHOGENIC], "Likely Pathogenic Variants", GH_C.orange);
        variant_table(sb, &first, arena, cv->by_category[GH_CV_RISK_FACTOR],
                      cv->counts[GH_CV_RISK_FACTOR], "Risk Factor Variants", GH_C.amber);
        variant_table(sb, &first, arena, cv->by_category[GH_CV_DRUG_RESPONSE],
                      cv->counts[GH_CV_DRUG_RESPONSE], "Drug Response Variants", GH_C.blue);

        size_t nprot = cv->counts[GH_CV_PROTECTIVE];
        if (nprot) {
            nl(sb, &first);
            gh_sb_printf(sb, "<h3>Protective Variants <span class=\"badge\" "
                             "style=\"background:var(--green)\">%zu</span></h3>", nprot);
            nl(sb, &first);
            gh_sb_puts(sb, "<div class=\"good-news-grid\">");
            for (size_t i = 0; i < nprot; i++) {
                const gh_cv_finding *v = &cv->by_category[GH_CV_PROTECTIVE][i];
                nl(sb, &first);
                gh_sb_puts(sb, "<div class=\"good-news-card\"><strong>");
                esc(sb, v->gene ? v->gene : "Unknown");
                gh_sb_puts(sb, "</strong>: ");
                esc(sb, gh_report_clean_condition(arena, or_empty(v->traits)));
                gh_sb_puts(sb, "</div>");
            }
            nl(sb, &first);
            gh_sb_puts(sb, "</div>");
        }
    }

    nl(sb, &first);
    gh_sb_puts(sb, "<p style=\"font-size:.85em;color:var(--accent2)\">"
                   "PRS estimates relative genetic risk from common variants. "
                   "Lifestyle, environment, and rare variants also affect risk. "
                   "Not a clinical diagnosis.</p>");
}

/* ------------------------------------------------------------------ */
/* 5. Body profile                                                     */
/* ------------------------------------------------------------------ */

static const char *trait_label(const char *key)
{
    static const gh_kv labels[] = {
        {"eye_color", "Eye Color"}, {"hair_color", "Hair Color"},
        {"earwax_type", "Earwax Type"}, {"freckling", "Freckling / Sun Sensitivity"},
        {"lactose_tolerance", "Lactose Tolerance"}, {"bitter_taste", "Bitter Taste (PTC)"},
        {"cilantro_taste", "Cilantro Taste"}, {"asparagus_smell", "Asparagus Smell Detection"},
        {"muscle_fiber_type", "Muscle Fiber Type (ACTN3)"},
        {"secretor_status", "FUT2 Secretor Status"}, {"photic_sneeze", "Photic Sneeze Reflex"},
        {"hair_curl", "Hair Curl / Texture"}, {"baldness_risk", "Baldness Susceptibility"},
        {"unibrow_tendency", "Unibrow Tendency (PAX3)"},
    };
    for (size_t i = 0; i < sizeof labels / sizeof labels[0]; i++)
        if (strcmp(labels[i].key, key) == 0)
            return labels[i].value;
    return NULL;
}

static const char *level_color(const char *level)
{
    /* {"elevated": red, "moderate": amber, "low": green} */
    if (strcmp(level, "elevated") == 0) return GH_C.red;
    if (strcmp(level, "moderate") == 0) return GH_C.amber;
    if (strcmp(level, "low") == 0) return GH_C.green;
    return "var(--border)";
}

static void list_items(gh_strbuf *sb, bool *first, const char *const *items,
                       size_t n, bool escape)
{
    for (size_t i = 0; i < n; i++) {
        nl(sb, first);
        gh_sb_puts(sb, "<li>");
        if (escape)
            esc(sb, items[i]);
        else
            raw(sb, items[i]);
        gh_sb_puts(sb, "</li>");
    }
}

void gh_report_body_profile(gh_strbuf *sb, gh_arena *arena, const gh_analysis *a,
                            const gh_report_scorers *sc)
{
    bool first = true;
    nl(sb, &first);
    gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                   "Traits, blood type, chronotype, longevity, and athletic profile.</p>");

    if (sc && sc->traits && sc->ntraits) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Predicted Traits</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<table><tr><th>Trait</th><th>Prediction</th>"
                       "<th>Confidence</th><th>What It Means</th></tr>");
        for (size_t i = 0; i < sc->ntraits; i++) {
            const gh_trait_result *t = &sc->traits[i];
            const char *label = trait_label(t->key);
            if (!label)
                label = py_title(arena, t->key, true);
            const char *conf_color =
                strcmp(t->confidence, "high") == 0 ? "var(--green)" :
                strcmp(t->confidence, "moderate") == 0 ? "var(--accent)" :
                strcmp(t->confidence, "low") == 0 ? "var(--warn)" : "inherit";
            nl(sb, &first);
            gh_sb_printf(sb,
                "<tr><td><strong>%s</strong></td><td>%s</td>"
                "<td style=\"color:%s\">%s</td><td style=\"font-size:.85em\">%s</td></tr>",
                label, t->prediction, conf_color, py_title(arena, t->confidence, false),
                t->description);
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</table>");
    }

    if (sc && sc->blood && strcmp(sc->blood->blood_type, "Unknown") != 0) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Blood Type</h3>");
        nl(sb, &first);
        gh_sb_printf(sb,
            "<div style=\"text-align:center;margin:1em 0\">"
            "<span style=\"font-size:3em;font-weight:bold;color:var(--warn)\">%s</span>"
            "<br><span style=\"font-size:.9em;color:var(--accent2)\">"
            "Confidence: %s</span></div>",
            sc->blood->blood_type, py_title(arena, sc->blood->confidence, false));
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.85em;color:var(--accent2)\">Blood type from SNP data has limitations. "
                       "Confirm with clinical blood typing.</p>");
    }

    const gh_sleep_result *sl = sc ? sc->sleep : NULL;
    if (sl) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Sleep &amp; Chronotype</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "Genetic influence on circadian rhythm and sleep architecture.</p>");
        nl(sb, &first);
        gh_sb_puts(sb, "<div style=\"text-align:center;margin:1em 0\">"
                       "<p style=\"font-size:1.5em;font-weight:bold\">");
        esc(sb, sl->chronotype ? sl->chronotype : "Unknown");
        gh_sb_printf(sb, "</p><p style=\"color:var(--accent2)\">Chronotype Score: %s/100 "
                         "(0=extreme morning, 100=extreme evening)</p></div>",
                     gh_report_pyfloat(arena, sl->chronotype_score));
        nl(sb, &first);
        gh_sb_puts(sb, "<table><tr><td><strong>Optimal Sleep Window</strong></td><td>");
        esc(sb, sl->optimal_sleep_window);
        gh_sb_puts(sb, "</td></tr><tr><td><strong>Peak Alertness</strong></td><td>");
        esc(sb, sl->peak_alertness);
        gh_sb_puts(sb, "</td></tr><tr><td><strong>Last Caffeine By</strong></td><td>");
        esc(sb, sl->caffeine_cutoff);
        gh_sb_puts(sb, "</td></tr></table>");
        if (sl->nrecommendations) {
            nl(sb, &first);
            gh_sb_puts(sb, "<p><strong>Sleep Recommendations:</strong></p><ul>");
            list_items(sb, &first, sl->recommendations, sl->nrecommendations, true);
            nl(sb, &first);
            gh_sb_puts(sb, "</ul>");
        }
    }

    const gh_longevity_result *lg = sc ? sc->longevity : NULL;
    if (lg) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Longevity &amp; Healthspan</h3>");
        double score = lg->longevity_score;
        const char *score_color = score >= 60 ? GH_C.green : score >= 40 ? GH_C.amber : GH_C.red;
        nl(sb, &first);
        gh_sb_printf(sb,
            "<div style=\"text-align:center;margin:1em 0\">"
            "<span style=\"font-size:2.5em;font-weight:bold;color:%s\">%s</span>"
            "<span style=\"font-size:1.2em;color:var(--accent2)\">/100</span>"
            "<p style=\"color:var(--accent2)\">Longevity Genetic Score</p>"
            "<p style=\"font-size:.9em\">", score_color, gh_report_pyfloat(arena, score));
        esc(sb, lg->summary);
        gh_sb_puts(sb, "</p></div>");
        if (lg->ndomains) {
            nl(sb, &first);
            gh_sb_puts(sb, "<div style=\"display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:1em\">");
            for (size_t i = 0; i < lg->ndomains; i++) {
                const gh_healthspan_result *d = &lg->domains[i];
                const char *c = d->score >= 60 ? GH_C.green : d->score >= 40 ? GH_C.amber : GH_C.red;
                nl(sb, &first);
                gh_sb_puts(sb, "<div style=\"border:1px solid var(--border);border-radius:8px;padding:1em;text-align:center\"><strong>");
                esc(sb, py_title(arena, d->name, true));
                gh_sb_printf(sb, "</strong><br><span style=\"font-size:1.5em;color:%s\">%d</span>/90<br>"
                                 "<span style=\"font-size:.85em;color:var(--accent2)\">%s</span></div>",
                             c, d->score, or_empty(d->rating));
            }
            nl(sb, &first);
            gh_sb_puts(sb, "</div>");
        }
        if (lg->ninterventions) {
            nl(sb, &first);
            gh_sb_puts(sb, "<p><strong>Genetically-Supported Interventions:</strong></p><ul>");
            for (size_t i = 0; i < lg->ninterventions; i++) {
                nl(sb, &first);
                gh_sb_puts(sb, "<li><strong>");
                esc(sb, lg->interventions[i].intervention);
                gh_sb_puts(sb, "</strong> \xe2\x80\x94 ");
                esc(sb, lg->interventions[i].why);
                gh_sb_puts(sb, "</li>");
            }
            nl(sb, &first);
            gh_sb_puts(sb, "</ul>");
        }
    }

    for (size_t i = 0; a && i < a->nfindings; i++) {
        const gh_finding *f = &a->findings[i];
        if (!f->gene || strcmp(f->gene, "ACTN3") != 0)
            continue;
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Athletic Profile</h3>");
        nl(sb, &first);
        gh_sb_printf(sb, "<p class=\"eli5\">%s</p>", gh_report_eli5_gene("ACTN3"));
        nl(sb, &first);
        gh_sb_printf(sb, "<p><strong>ACTN3</strong>: %s \xe2\x80\x94 %s</p>",
                     py_title(arena, f->status, true), or_empty(f->description));
        nl(sb, &first);
        gh_report_paper_refs(sb, "rs1815739");
        break;
    }

    const gh_pain_result *pain = sc ? sc->pain : NULL;
    if (pain && pain->base.nhits > 0) {
        int score = pain->sensitivity_score;
        const char *color = score >= 70 ? GH_C.red : score >= 40 ? GH_C.amber : GH_C.green;
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Pain Sensitivity Profile</h3>");
        nl(sb, &first);
        gh_sb_printf(sb,
            "<div style=\"text-align:center;margin:1em 0\">"
            "<span style=\"font-size:2em;font-weight:bold;color:%s\">%d</span>"
            "<span style=\"color:var(--accent2)\">/100</span>"
            "<p style=\"color:var(--accent2)\">Pain Sensitivity Score (higher = more sensitive)</p></div>",
            color, score);
        nl(sb, &first);
        gh_sb_puts(sb, "<p>");
        esc(sb, pain->base.summary);
        gh_sb_puts(sb, "</p>");
        if (pain->base.nrecommendations) {
            nl(sb, &first);
            gh_sb_puts(sb, "<ul>");
            list_items(sb, &first, pain->base.recommendations, pain->base.nrecommendations, true);
            nl(sb, &first);
            gh_sb_puts(sb, "</ul>");
        }
    }

    const gh_histamine_result *hi = sc ? sc->histamine : NULL;
    if (hi && hi->base.nhits > 0) {
        const char *level = hi->risk_level ? hi->risk_level : "low";
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Histamine Intolerance Risk</h3>");
        nl(sb, &first);
        gh_sb_printf(sb, "<p><span class=\"mag-badge\" style=\"background:%s;color:#fff\">%s</span> ",
                     level_color(level), py_upper(arena, level));
        esc(sb, hi->base.summary);
        gh_sb_puts(sb, "</p>");
        if (hi->nfoods && strcmp(level, "low") != 0) {
            nl(sb, &first);
            gh_sb_puts(sb, "<p><strong>Foods to watch:</strong> ");
            for (size_t i = 0; i < hi->nfoods && i < 8; i++) {
                if (i)
                    gh_sb_puts(sb, ", ");
                esc(sb, hi->foods_to_watch[i]);
            }
            gh_sb_puts(sb, "</p>");
        }
    }

    const gh_alcohol_result *al = sc ? sc->alcohol : NULL;
    if (al && al->base.nhits > 0) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Alcohol Metabolism</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p>");
        esc(sb, al->base.summary);
        gh_sb_puts(sb, "</p>");
        nl(sb, &first);
        gh_sb_puts(sb, "<table><tr><td><strong>Metabolism Speed</strong></td><td>");
        esc(sb, al->metabolism_speed ? al->metabolism_speed : "Unknown");
        gh_sb_puts(sb, "</td></tr><tr><td><strong>Flush Risk</strong></td><td>");
        esc(sb, al->flush_risk ? al->flush_risk : "Unknown");
        gh_sb_puts(sb, "</td></tr><tr><td><strong>Cancer Risk from Alcohol</strong></td><td>");
        esc(sb, al->cancer_risk ? al->cancer_risk : "Unknown");
        gh_sb_puts(sb, "</td></tr></table>");
    }

    /* Eye and thyroid: the Python builder reads info["risk_level"] and
     * info["detail"], but the profile modules write "level" and "markers".
     * So every badge comes out as the default "AVERAGE" in green with an
     * empty detail, whatever the genotype says. Reproduced as is. */
    const gh_eye_result *eye = sc ? sc->eye : NULL;
    if (eye && eye->base.nhits > 0) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Eye Health Genetics</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p>");
        esc(sb, eye->base.summary);
        gh_sb_puts(sb, "</p>");
        for (size_t i = 0; i < eye->nconditions; i++) {
            nl(sb, &first);
            gh_sb_printf(sb, "<p><span class=\"mag-badge\" style=\"background:%s;color:#fff\">"
                             "AVERAGE</span> <strong>%s</strong>: </p>",
                         GH_C.green, py_title(arena, eye->conditions[i].name, true));
        }
    }

    const gh_thyroid_result *th = sc ? sc->thyroid : NULL;
    if (th && th->base.nhits > 0) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Thyroid Genetics</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p>");
        esc(sb, th->base.summary);
        gh_sb_puts(sb, "</p>");
        for (size_t i = 0; i < th->ndomains; i++) {
            nl(sb, &first);
            gh_sb_printf(sb, "<p><span class=\"mag-badge\" style=\"background:%s;color:#fff\">"
                             "AVERAGE</span> <strong>%s</strong></p>",
                         GH_C.green, py_title(arena, th->domains[i].name, true));
        }
    }

    const gh_hormone_result *ho = sc ? sc->hormone : NULL;
    if (ho && ho->base.nhits > 0) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Hormone Metabolism</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p>");
        esc(sb, ho->base.summary);
        gh_sb_puts(sb, "</p>");
    }
}

/* ------------------------------------------------------------------ */
/* 6. Mental health                                                    */
/* ------------------------------------------------------------------ */

void gh_report_mental_health(gh_strbuf *sb, gh_arena *arena,
                             const gh_report_scorers *sc)
{
    const gh_mental_result *mh = sc ? sc->mental : NULL;
    if (!mh) {
        gh_sb_puts(sb, "<p>No mental health genetic data available.</p>");
        return;
    }
    bool first = true;
    nl(sb, &first);
    gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                   "Genetic susceptibility markers. Environment and lifestyle are equally important.</p>");
    nl(sb, &first);
    gh_sb_puts(sb, "<p><strong>");
    esc(sb, mh->summary);
    gh_sb_puts(sb, "</strong></p>");

    if (mh->ndomains) {
        nl(sb, &first);
        gh_sb_puts(sb, "<div style=\"display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:1em;margin:1em 0\">");
        for (size_t i = 0; i < mh->ndomains; i++) {
            const gh_mental_domain_result *d = &mh->domains[i];
            nl(sb, &first);
            gh_sb_puts(sb, "<div style=\"border:1px solid var(--border);border-radius:8px;padding:1em;text-align:center\"><strong>");
            esc(sb, py_title(arena, d->name, true));
            gh_sb_printf(sb, "</strong><br><span style=\"font-size:1.3em;color:%s\">",
                         level_color(d->risk_level));
            esc(sb, py_title(arena, d->risk_level, false));
            gh_sb_puts(sb, "</span></div>");
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</div>");
    }

    if (mh->nrisk_factors) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Risk Factors</h3><ul>");
        list_items(sb, &first, mh->risk_factors, mh->nrisk_factors, true);
        nl(sb, &first);
        gh_sb_puts(sb, "</ul>");
    }
    if (mh->nresilience_factors) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Resilience Factors</h3><ul>");
        list_items(sb, &first, mh->resilience_factors, mh->nresilience_factors, true);
        nl(sb, &first);
        gh_sb_puts(sb, "</ul>");
    }
    if (mh->ntreatment_notes) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Treatment Matching Notes</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "Share these with your therapist or psychiatrist.</p>");
        nl(sb, &first);
        gh_sb_puts(sb, "<ul>");
        list_items(sb, &first, mh->treatment_notes, mh->ntreatment_notes, true);
        nl(sb, &first);
        gh_sb_puts(sb, "</ul>");
    }
    if (mh->nrecommendations) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Recommendations</h3><ol>");
        list_items(sb, &first, mh->recommendations, mh->nrecommendations, true);
        nl(sb, &first);
        gh_sb_puts(sb, "</ol>");
    }
}

/* ------------------------------------------------------------------ */
/* 7. Clinical detail                                                  */
/* ------------------------------------------------------------------ */

typedef struct {
    const gh_finding *f;
    size_t order;
} ranked_finding;

static int compare_ranked(const void *a, const void *b)
{
    const ranked_finding *x = a, *y = b;
    if (x->f->magnitude != y->f->magnitude)
        return x->f->magnitude > y->f->magnitude ? -1 : 1;
    return (x->order > y->order) - (x->order < y->order);
}

static int compare_strp(const void *a, const void *b)
{
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

void gh_report_clinical_detail(gh_strbuf *sb, gh_arena *arena, const gh_analysis *a,
                               const gh_report_scorers *sc)
{
    bool first = true;
    nl(sb, &first);
    gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                   "Complete technical findings for clinical review. "
                   "Searchable, sortable, and exportable.</p>");

    size_t nf = a ? a->nfindings : 0;
    const char **cats = gh_calloc(arena, nf + 1, sizeof *cats);
    size_t ncats = 0;
    for (size_t i = 0; i < nf; i++) {
        const char *name = a->findings[i].category ? a->findings[i].category : "Other";
        size_t j = 0;
        for (; j < ncats; j++)
            if (strcmp(cats[j], name) == 0)
                break;
        if (j == ncats)
            cats[ncats++] = name;
    }
    qsort(cats, ncats, sizeof *cats, compare_strp);
    ranked_finding *ranked = gh_calloc(arena, nf + 1, sizeof *ranked);

    for (size_t c = 0; c < ncats; c++) {
        size_t n = 0;
        for (size_t i = 0; i < nf; i++) {
            const char *name = a->findings[i].category ? a->findings[i].category : "Other";
            if (strcmp(name, cats[c]) == 0) {
                ranked[n].f = &a->findings[i];
                ranked[n].order = n;
                n++;
            }
        }
        qsort(ranked, n, sizeof *ranked, compare_ranked);

        nl(sb, &first);
        gh_sb_printf(sb, "<details class=\"category-section\" open>"
                         "<summary><h3 class=\"inline-h3\">%s <span class=\"badge\">%zu</span></h3></summary>",
                     cats[c], n);
        for (size_t i = 0; i < n; i++) {
            const gh_finding *f = ranked[i].f;
            const char *cls = gh_report_mag_class(f->magnitude);
            nl(sb, &first);
            gh_sb_printf(sb, "<div class=\"finding-card %s\">", cls);
            nl(sb, &first);
            gh_sb_printf(sb, "<div class=\"finding-header\"><span class=\"mag-badge mag-%s\">%d/6</span> <strong>",
                         cls, f->magnitude);
            esc(sb, f->gene ? f->gene : "Unknown");
            gh_sb_puts(sb, "</strong> <code>");
            esc(sb, f->rsid);
            gh_sb_puts(sb, "</code> \xe2\x80\x94 <code>");
            esc(sb, f->genotype);
            gh_sb_puts(sb, "</code> \xe2\x80\x94 ");
            esc(sb, py_title(arena, f->status, true));
            gh_sb_puts(sb, "</div>");
            const char *eli5 = gh_report_eli5_gene(f->gene);
            if (*eli5) {
                nl(sb, &first);
                gh_sb_printf(sb, "<p class=\"eli5-inline\" style=\"font-size:.85em;margin:.2em 0\">%s</p>", eli5);
            }
            nl(sb, &first);
            gh_sb_puts(sb, "<p class=\"finding-desc\">");
            esc(sb, f->description);
            gh_sb_puts(sb, "</p>");
            if (f->note && *f->note) {
                nl(sb, &first);
                gh_sb_puts(sb, "<p class=\"finding-note\">Note: ");
                esc(sb, f->note);
                gh_sb_puts(sb, "</p>");
            }
            if (f->freq) {
                static const char *const labels[GH_MAX_POPULATIONS] = {
                    "European", "African", "East Asian", "South Asian", "American"};
                int order[GH_MAX_POPULATIONS];
                size_t count = 0;
                for (int p = 0; p < (int)GH_POPULATION_COUNT; p++)
                    if (f->freq[p] > 0.001)
                        order[count++] = p;
                for (size_t x = 0; x + 1 < count; x++)
                    for (size_t y = 0; y + 1 < count - x; y++)
                        if (f->freq[order[y]] < f->freq[order[y + 1]]) {
                            int t = order[y]; order[y] = order[y + 1]; order[y + 1] = t;
                        }
                if (count) {
                    nl(sb, &first);
                    gh_sb_puts(sb, "<p class=\"finding-note\" style=\"font-size:.8em\">Population freq: ");
                    for (size_t x = 0; x < count; x++) {
                        if (x)
                            gh_sb_puts(sb, " &middot; ");
                        gh_sb_printf(sb, "%s: %.0f%%", labels[order[x]],
                                     nearbyint(f->freq[order[x]] * 100.0));
                    }
                    gh_sb_puts(sb, "</p>");
                }
            }
            nl(sb, &first);
            gh_report_db_links(sb, f->rsid);
            nl(sb, &first);
            gh_report_paper_refs(sb, f->rsid);
            nl(sb, &first);
            gh_sb_puts(sb, "</div>");
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</details>");
    }

    /* Pathway analysis: {f["gene"]: f} keeps the LAST finding for a gene. */
    nl(sb, &first);
    gh_sb_puts(sb, "<h3>Pathway Analysis</h3>");
    for (size_t p = 0; p < GH_PATHWAY_COUNT; p++) {
        const gh_pathway *pw = &GH_PATHWAYS[p];
        const gh_finding *hits[64];
        size_t nhits = 0;
        for (size_t gi = 0; gi < pw->ngenes && nhits < 64; gi++) {
            const gh_finding *hit = NULL;
            for (size_t i = 0; i < nf; i++)
                if (a->findings[i].gene && strcmp(a->findings[i].gene, pw->genes[gi]) == 0)
                    hit = &a->findings[i];
            if (hit)
                hits[nhits++] = hit;
        }
        if (!nhits)
            continue;
        nl(sb, &first);
        gh_sb_puts(sb, "<details><summary><strong>");
        esc(sb, pw->name);
        gh_sb_puts(sb, "</strong></summary><ul>");
        for (size_t i = 0; i < nhits; i++) {
            nl(sb, &first);
            gh_sb_printf(sb, "<li><span class=\"mag-dot mag-%s\"></span> <strong>",
                         gh_report_mag_class(hits[i]->magnitude));
            esc(sb, hits[i]->gene);
            gh_sb_puts(sb, "</strong>: ");
            esc(sb, py_title(arena, hits[i]->status, true));
            gh_sb_puts(sb, "</li>");
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</ul></details>");
    }

    /* Epistasis. pipeline.py strips genes_involved before the report reads
     * the JSON, so the "Genes" line falls back to the escaped name. */
    if (sc && sc->epistasis && sc->nepistasis) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Gene-Gene Interactions (Epistasis)</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "Combined effects when multiple gene variants interact.</p>");
        for (size_t i = 0; i < sc->nepistasis; i++) {
            const gh_epi_result *e = &sc->epistasis[i];
            const char *color = strcmp(e->risk_level, "high") == 0 ? "var(--warn)" :
                                strcmp(e->risk_level, "moderate") == 0 ? "var(--accent)" :
                                strcmp(e->risk_level, "low") == 0 ? "var(--green)" : "var(--accent)";
            nl(sb, &first);
            gh_sb_printf(sb, "<details open><summary><span class=\"mag-badge\" "
                             "style=\"background:%s;color:#fff\">%s</span> <strong>%s</strong></summary>",
                         color, py_upper(arena, e->risk_level), e->name);
            nl(sb, &first);
            gh_sb_puts(sb, "<p><strong>Genes:</strong> ");
            esc(sb, e->name);
            gh_sb_puts(sb, "</p>");
            nl(sb, &first);
            gh_sb_printf(sb, "<p><strong>Effect:</strong> %s</p>", or_empty(e->effect));
            nl(sb, &first);
            gh_sb_printf(sb, "<p><strong>Mechanism:</strong> %s</p>", or_empty(e->mechanism));
            nl(sb, &first);
            gh_sb_puts(sb, "<p><strong>Recommended Actions:</strong></p><ul>");
            list_items(sb, &first, e->actions, e->nactions, false);
            nl(sb, &first);
            gh_sb_puts(sb, "</ul></details>");
        }
    }

    const gh_carrier_result *cr = sc ? sc->carriers : NULL;
    if (cr && cr->ncarriers) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Carrier Screening</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "Carrier status: you have one copy of a recessive variant. "
                       "Relevant for reproductive planning.</p>");
        nl(sb, &first);
        gh_sb_printf(sb, "<p>%zu carrier finding(s) organized by disease system.</p>", cr->ncarriers);

        const char **systems = gh_calloc(arena, cr->nsystems + 1, sizeof *systems);
        for (size_t s = 0; s < cr->nsystems; s++)
            systems[s] = cr->systems[s];
        qsort(systems, cr->nsystems, sizeof *systems, compare_strp);
        for (size_t s = 0; s < cr->nsystems; s++) {
            nl(sb, &first);
            gh_sb_printf(sb, "<h4>%s</h4><ul>", systems[s]);
            for (size_t i = 0; i < cr->ncarriers; i++) {
                const gh_carrier *c = &cr->carriers[i];
                if (strcmp(c->system, systems[s]) != 0)
                    continue;
                nl(sb, &first);
                gh_sb_printf(sb, "<li><strong>%s</strong> (%s): %s <span class=\"badge\">%s</span>",
                             c->gene, or_empty(c->rsid), c->condition, c->inheritance);
                if (c->reproductive_note && *c->reproductive_note)
                    gh_sb_printf(sb, " &mdash; <em>%s</em>", c->reproductive_note);
                gh_sb_puts(sb, "</li>");
            }
            nl(sb, &first);
            gh_sb_puts(sb, "</ul>");
        }
        if (cr->ncouples_relevant) {
            nl(sb, &first);
            gh_sb_puts(sb, "<h4>Couples-Relevant Conditions</h4><ul>");
            for (size_t i = 0; i < cr->ncarriers; i++) {
                if (!cr->carriers[i].couples_relevant)
                    continue;
                nl(sb, &first);
                gh_sb_printf(sb, "<li><strong>%s</strong>: %s</li>",
                             cr->carriers[i].gene, cr->carriers[i].condition);
            }
            nl(sb, &first);
            gh_sb_puts(sb, "</ul>");
        }
    }
}

/* ------------------------------------------------------------------ */
/* 8. Ancestry                                                         */
/* ------------------------------------------------------------------ */

void gh_report_ancestry(gh_strbuf *sb, gh_arena *arena, const gh_report_scorers *sc)
{
    const gh_ancestry_result *anc = sc ? sc->ancestry : NULL;
    const gh_mt_result *mt = sc ? sc->mt : NULL;
    bool mt_known = mt && strcmp(mt->haplogroup, "Unknown") != 0;
    bool have_aims = anc && anc->markers_found > 0;

    if (!have_aims && !mt_known) {
        gh_sb_puts(sb, "<p>No ancestry-informative markers found in genome data.</p>");
        return;
    }
    bool first = true;

    if (have_aims) {
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "Superpopulation estimates from ~55 ancestry-informative markers.</p>");
        nl(sb, &first);
        gh_sb_puts(sb, "<div class=\"chart-grid\">");
        nl(sb, &first);
        gh_sb_puts(sb, "<div>");
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Ancestry Proportions</h3>");
        nl(sb, &first);
        gh_report_svg_ancestry_donut(sb, arena, anc);
        nl(sb, &first);
        gh_sb_puts(sb, "</div>");
        nl(sb, &first);
        gh_sb_puts(sb, "<div>");
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Details</h3>");
        nl(sb, &first);
        gh_sb_printf(sb, "<p><strong>Confidence:</strong> %s (%zu markers)</p>",
                     py_title(arena, anc->confidence, false), anc->markers_found);
        nl(sb, &first);
        gh_sb_puts(sb, "<table><tr><th>Population</th><th>Proportion</th></tr>");

        size_t order[GH_MAX_POPULATIONS];
        for (size_t i = 0; i < GH_POPULATION_COUNT; i++)
            order[i] = i;
        for (size_t i = 0; i + 1 < GH_POPULATION_COUNT; i++)
            for (size_t j = 0; j + 1 < GH_POPULATION_COUNT - i; j++)
                if (anc->proportions[order[j]] < anc->proportions[order[j + 1]]) {
                    size_t t = order[j]; order[j] = order[j + 1]; order[j + 1] = t;
                }
        for (size_t i = 0; i < GH_POPULATION_COUNT; i++) {
            size_t p = order[i];
            nl(sb, &first);
            gh_sb_printf(sb, "<tr><td>%s (%s)</td><td>%.1f%%</td></tr>",
                         GH_POPULATION_LABELS[p], GH_POPULATIONS[p],
                         anc->proportions[p] * 100.0);
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</table>");
        nl(sb, &first);
        gh_sb_puts(sb, "</div>");
        nl(sb, &first);
        gh_sb_puts(sb, "</div>");
    }

    if (anc && anc->sub.present && anc->sub.nsubs) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Sub-Population Estimate</h3>");
        nl(sb, &first);
        gh_sb_printf(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                         "Finer-grained ancestry estimate within %s. %s confidence (%zu markers).</p>",
                     anc->top_ancestry && *anc->top_ancestry ? anc->top_ancestry : "top population",
                     py_title(arena, anc->sub.confidence, false), anc->sub.markers_used);
        nl(sb, &first);
        gh_sb_puts(sb, "<table><tr><th>Sub-Population</th><th>Proportion</th></tr>");
        size_t order[GH_MAX_SUB_POPULATIONS];
        for (size_t i = 0; i < anc->sub.nsubs; i++)
            order[i] = i;
        for (size_t i = 0; i + 1 < anc->sub.nsubs; i++)
            for (size_t j = 0; j + 1 < anc->sub.nsubs - i; j++)
                if (anc->sub.proportions[order[j]] < anc->sub.proportions[order[j + 1]]) {
                    size_t t = order[j]; order[j] = order[j + 1]; order[j + 1] = t;
                }
        for (size_t i = 0; i < anc->sub.nsubs; i++) {
            nl(sb, &first);
            gh_sb_printf(sb, "<tr><td>%s</td><td>%.1f%%</td></tr>",
                         anc->sub.labels[order[i]], anc->sub.proportions[order[i]] * 100.0);
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</table>");
    }

    if (mt_known) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Maternal Haplogroup (Mitochondrial DNA)</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "Maternal lineage from mitochondrial DNA.</p>");
        nl(sb, &first);
        gh_sb_printf(sb,
            "<div style=\"text-align:center;margin:1em 0\">"
            "<span style=\"font-size:2.5em;font-weight:bold;color:var(--accent)\">%s</span>"
            "<br><span style=\"font-size:.95em\">%s</span>"
            "<br><span style=\"font-size:.85em;color:var(--accent2)\">"
            "%s \xe2\x80\x94 %s confidence (%zu/%zu markers)</span></div>",
            mt->haplogroup, or_empty(mt->description), or_empty(mt->lineage),
            py_title(arena, mt->confidence, false), mt->markers_found, mt->markers_tested);
    }
}

/* ------------------------------------------------------------------ */
/* 9. Nutrigenomics                                                    */
/* ------------------------------------------------------------------ */

size_t gh_report_supplement_priorities(const gh_nutrigenomics_result *nut,
                                       const gh_nutrient_need **out, size_t cap)
{
    /* Python appends these while walking the profiles in table order, then
     * stable-sorts by priority; the needs list itself is sorted by severity
     * afterwards, so it cannot be used as the source order. */
    size_t n = 0;
    for (size_t p = 0; p < GH_NUTRIENT_PROFILE_COUNT && n < cap; p++) {
        for (size_t i = 0; i < nut->nneeds; i++) {
            const gh_nutrient_need *need = &nut->needs[i];
            if (need->profile != &GH_NUTRIENT_PROFILES[p])
                continue;
            if (strcmp(need->need_level, "high") == 0 || strcmp(need->need_level, "moderate") == 0)
                out[n++] = need;
            break;
        }
    }
    /* Stable partition: high before moderate, original order within each. */
    size_t w = 0;
    const gh_nutrient_need *tmp[64];
    size_t k = 0;
    for (size_t i = 0; i < n; i++)
        if (strcmp(out[i]->need_level, "high") == 0)
            tmp[k++] = out[i];
    for (size_t i = 0; i < n; i++)
        if (strcmp(out[i]->need_level, "high") != 0)
            tmp[k++] = out[i];
    for (; w < n; w++)
        out[w] = tmp[w];
    return n;
}

static bool nutrition_narrative(const char *id)
{
    return strcmp(id, "methylation_choline") == 0 || strcmp(id, "vitamin_d_profile") == 0 ||
           strcmp(id, "iron_profile") == 0 || strcmp(id, "caffeine_metabolism") == 0;
}

void gh_report_nutrigenomics(gh_strbuf *sb, gh_arena *arena, const gh_report_scorers *sc)
{
    const gh_nutrigenomics_result *nut = sc ? sc->nutrition : NULL;
    const gh_insights *ins = sc ? sc->insights : NULL;
    bool first = true;

    if (nut) {
        nl(sb, &first);
        gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                       "How your genes affect nutrient needs. Personalized supplement and dietary guidance.</p>");
        nl(sb, &first);
        gh_sb_puts(sb, "<p>");
        esc(sb, nut->summary);
        gh_sb_puts(sb, "</p>");

        for (size_t i = 0; i < nut->nneeds; i++) {
            const gh_nutrient_need *n = &nut->needs[i];
            if (strcmp(n->need_level, "normal") == 0 && !n->nimpacts)
                continue;
            const char *color =
                strcmp(n->need_level, "high") == 0 ? GH_C.red :
                strcmp(n->need_level, "moderate") == 0 ? GH_C.amber :
                strcmp(n->need_level, "low") == 0 ? GH_C.blue :
                strcmp(n->need_level, "normal") == 0 ? GH_C.green :
                strcmp(n->need_level, "caution_excess") == 0 ? GH_C.purple : "var(--border)";
            nl(sb, &first);
            gh_sb_printf(sb, "<details class=\"rec-card\" style=\"border-left:4px solid %s\">"
                             "<summary><span class=\"mag-badge\" style=\"background:%s;color:#fff\">",
                         color, color);
            esc(sb, replace_all(arena, py_upper(arena, n->need_level), "_", " "));
            gh_sb_puts(sb, "</span> <strong>");
            esc(sb, n->profile->name);
            gh_sb_puts(sb, "</strong></summary>");
            nl(sb, &first);
            gh_sb_puts(sb, "<p>");
            esc(sb, n->recommendation);
            gh_sb_puts(sb, "</p>");
            nl(sb, &first);
            gh_sb_puts(sb, "<p><strong>Food sources:</strong> ");
            esc(sb, n->profile->food_sources);
            gh_sb_puts(sb, "</p>");
            for (size_t k = 0; k < n->nimpacts; k++) {
                nl(sb, &first);
                gh_sb_puts(sb, "<p style=\"font-size:.9em\">\xe2\x80\x94 ");
                esc(sb, n->impacts[k].gene);
                gh_sb_puts(sb, ": ");
                esc(sb, n->impacts[k].impact);
                gh_sb_puts(sb, "</p>");
            }
            nl(sb, &first);
            gh_sb_puts(sb, "</details>");
        }

        const gh_nutrient_need *supps[64];
        size_t nsupps = gh_report_supplement_priorities(nut, supps, 64);
        if (nsupps) {
            nl(sb, &first);
            gh_sb_puts(sb, "<h3>Supplement Priority List</h3>");
            nl(sb, &first);
            gh_sb_puts(sb, "<table><tr><th>Nutrient</th><th>Form</th><th>Dose</th><th>Priority</th></tr>");
            for (size_t i = 0; i < nsupps; i++) {
                nl(sb, &first);
                gh_sb_puts(sb, "<tr><td><strong>");
                esc(sb, supps[i]->profile->name);
                gh_sb_puts(sb, "</strong></td><td>");
                esc(sb, supps[i]->profile->supplement_form);
                gh_sb_puts(sb, "</td><td>");
                esc(sb, supps[i]->profile->dose_range);
                gh_sb_puts(sb, "</td><td>");
                esc(sb, supps[i]->need_level);
                gh_sb_puts(sb, "</td></tr>");
            }
            nl(sb, &first);
            gh_sb_puts(sb, "</table>");
        }

        /* testing_recommendations: one per profile, in table order. */
        if (GH_NUTRIENT_PROFILE_COUNT) {
            nl(sb, &first);
            gh_sb_puts(sb, "<h3>Recommended Lab Tests</h3><ul>");
            for (size_t p = 0; p < GH_NUTRIENT_PROFILE_COUNT; p++) {
                nl(sb, &first);
                gh_sb_puts(sb, "<li><strong>");
                esc(sb, GH_NUTRIENT_PROFILES[p].name);
                gh_sb_puts(sb, ":</strong> ");
                esc(sb, GH_NUTRIENT_PROFILES[p].testing);
                gh_sb_puts(sb, "</li>");
            }
            nl(sb, &first);
            gh_sb_puts(sb, "</ul>");
        }
    }

    bool any_narrative = false;
    for (size_t i = 0; ins && i < ins->nnarratives; i++)
        if (nutrition_narrative(ins->narratives[i].id))
            any_narrative = true;
    if (any_narrative) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Nutritional Gene Stories</h3>");
        for (size_t i = 0; i < ins->nnarratives; i++) {
            const gh_narrative_result *n = &ins->narratives[i];
            if (!nutrition_narrative(n->id))
                continue;
            nl(sb, &first);
            gh_sb_printf(sb, "<details open class=\"rec-card\" style=\"border-left:4px solid var(--accent2)\">"
                             "<summary><strong>%s</strong> <span class=\"badge\" style=\"background:var(--accent2)\">"
                             "%zu genes</span></summary><p><strong>Genes:</strong> ", n->title, n->nmatched);
            for (size_t k = 0; k < n->nmatched; k++)
                gh_sb_printf(sb, "%s%s", k ? ", " : "", n->matched_genes[k]);
            gh_sb_printf(sb, "</p><p>%s</p><p><strong>What this means for you:</strong> %s</p>",
                         or_empty(n->narrative), or_empty(n->practical));
            if (n->nreferences) {
                /* Two parts in the Python: the opening text, then the joined
                 * links with the closing tag -- so a newline sits between. */
                nl(sb, &first);
                gh_sb_puts(sb, "<div class=\"paper-refs\">References: ");
                nl(sb, &first);
                for (size_t k = 0; k < n->nreferences; k++) {
                    const char *ref = n->references[k];
                    if (k)
                        gh_sb_puts(sb, " | ");
                    const char *tag = strstr(ref, "PMID:");
                    if (tag) {
                        /* ref.split("PMID:")[1].strip().rstrip(",.") */
                        const char *p = tag + 5;
                        while (*p && isspace((unsigned char)*p))
                            p++;
                        const char *e = p + strlen(p);
                        while (e > p && isspace((unsigned char)e[-1]))
                            e--;
                        while (e > p && (e[-1] == ',' || e[-1] == '.'))
                            e--;
                        gh_sb_printf(sb, "<a href=\"https://pubmed.ncbi.nlm.nih.gov/%.*s/\" "
                                         "target=\"_blank\" rel=\"noopener\">%s</a>",
                                     (int)(e - p), p, ref);
                    } else {
                        raw(sb, ref);
                    }
                }
                gh_sb_puts(sb, "</div>");
            }
            nl(sb, &first);
            gh_sb_puts(sb, "</details>");
        }
    }

    if (first)
        gh_sb_puts(sb, "<p>No nutrigenomics data available.</p>");
}

/* ------------------------------------------------------------------ */
/* 10. Quality                                                         */
/* ------------------------------------------------------------------ */

void gh_report_quality(gh_strbuf *sb, const gh_quality *q)
{
    if (!q) {
        gh_sb_puts(sb, "<p>No quality metrics available.</p>");
        return;
    }
    bool first = true;
    nl(sb, &first);
    gh_sb_puts(sb, "<p style=\"font-size:.9em;color:var(--accent2)\">"
                   "Raw data quality assessment \xe2\x80\x94 call rate, heterozygosity, and coverage.</p>");

    double call_pct = q->call_rate * 100.0;
    const char *badge_color, *badge_text;
    if (call_pct >= 99)      { badge_color = "var(--green)"; badge_text = "Excellent"; }
    else if (call_pct >= 97) { badge_color = GH_C.blue;      badge_text = "Good"; }
    else if (call_pct >= 95) { badge_color = GH_C.amber;     badge_text = "Fair"; }
    else                     { badge_color = "var(--warn)";  badge_text = "Low"; }

    nl(sb, &first);
    gh_sb_printf(sb, "<p><strong>Call Rate:</strong> <span class=\"mag-badge\" "
                     "style=\"background:%s;color:#fff\">%.1f%% \xe2\x80\x94 %s</span></p>",
                 badge_color, call_pct, badge_text);

    nl(sb, &first);
    gh_sb_puts(sb, "<table><tr><th>Metric</th><th>Value</th></tr><tr><td>Total SNPs</td><td>");
    grouped(sb, q->total_snps);
    gh_sb_puts(sb, "</td></tr><tr><td>No-call positions</td><td>");
    grouped(sb, q->no_call_count);
    gh_sb_puts(sb, "</td></tr><tr><td>Autosomal SNPs</td><td>");
    grouped(sb, q->autosomal_count);
    gh_sb_printf(sb, "</td></tr><tr><td>Mitochondrial SNPs</td><td>%zu</td></tr>"
                     "<tr><td>Heterozygosity rate</td><td>%.3f</td></tr>"
                     "<tr><td>Inferred sex</td><td>%s</td></tr></table>",
                 q->mt_snp_count, q->het_rate, q->has_y ? "Male (Y detected)" : "Female");

    size_t n = 0, max_cnt = 0;
    for (size_t i = 0; i < q->nchromosomes; i++) {
        const char *name = q->chromosomes[i].name;
        bool digits = *name != '\0';
        for (const char *p = name; *p && digits; p++)
            if (!isdigit((unsigned char)*p))
                digits = false;
        if (!digits)
            continue;
        n++;
        if (q->chromosomes[i].count > max_cnt)
            max_cnt = q->chromosomes[i].count;
    }
    if (!n)
        return;

    const int bar_w = 300, bar_h = 14;
    size_t height = n * (size_t)(bar_h + 4) + 20;
    nl(sb, &first);
    gh_sb_printf(sb, "<h3>Chromosome Coverage</h3><svg viewBox=\"0 0 440 %zu\" class=\"chart\" "
                     "role=\"img\" aria-label=\"Chromosome SNP coverage\">"
                     "<title>SNPs per Chromosome</title>", height);
    size_t index = 0;
    for (size_t i = 0; i < q->nchromosomes; i++) {
        const char *name = q->chromosomes[i].name;
        bool digits = *name != '\0';
        for (const char *p = name; *p && digits; p++)
            if (!isdigit((unsigned char)*p))
                digits = false;
        if (!digits)
            continue;
        size_t cnt = q->chromosomes[i].count;
        size_t y = 10 + index * (size_t)(bar_h + 4);
        int w = max_cnt ? (int)((double)cnt / (double)max_cnt * bar_w) : 0;
        gh_sb_printf(sb,
            "<text x=\"30\" y=\"%zu\" text-anchor=\"end\" fill=\"currentColor\" font-size=\"10\">%s</text>"
            "<rect x=\"35\" y=\"%zu\" width=\"%d\" height=\"%d\" rx=\"2\" fill=\"var(--accent)\" opacity=\"0.7\"/>"
            "<text x=\"%d\" y=\"%zu\" fill=\"currentColor\" font-size=\"9\">",
            y + 11, name, y, w, bar_h, 40 + w, y + 11);
        grouped(sb, cnt);
        gh_sb_puts(sb, "</text>");
        index++;
    }
    gh_sb_puts(sb, "</svg>");
}

/* ------------------------------------------------------------------ */
/* 11. Doctor card                                                     */
/* ------------------------------------------------------------------ */

void gh_report_doctor_card(gh_strbuf *sb, gh_arena *arena, const gh_report_scorers *sc)
{
    const gh_recommendations *recs = sc ? sc->recs : NULL;
    const gh_star_result *stars = sc ? sc->stars : NULL;
    size_t nstars = sc ? sc->nstars : 0;
    const gh_apoe_result *apoe = sc ? sc->apoe : NULL;
    const gh_acmg_result *acmg = sc ? sc->acmg : NULL;

    size_t nhigh = 0;
    for (size_t i = 0; recs && i < recs->npriorities; i++)
        if (strcmp(recs->priorities[i].priority, "high") == 0)
            nhigh++;
    size_t nreferrals = recs ? recs->nreferrals : 0;
    size_t nacmg = acmg ? acmg->nfindings : 0;
    if (!nhigh && !nreferrals && !(stars && nstars) && !nacmg) {
        gh_sb_puts(sb, "<p>No significant clinical findings for doctor review.</p>");
        return;
    }
    bool first = true;

    nl(sb, &first);
    gh_sb_puts(sb, "<div style=\"border:2px solid var(--accent);border-radius:8px;padding:1.5em\">");
    nl(sb, &first);
    gh_sb_puts(sb, "<h3 style=\"margin-top:0;color:var(--accent)\">Patient Genetic Summary</h3>");

    if (nhigh) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h4>High-Priority Conditions</h4>");
        nl(sb, &first);
        gh_sb_puts(sb, "<table><tr><th>Condition</th><th>Doctor Note</th></tr>");
        for (size_t i = 0; i < recs->npriorities; i++) {
            const gh_priority *p = &recs->priorities[i];
            if (strcmp(p->priority, "high") != 0)
                continue;
            const char *note = p->doctor_note && *p->doctor_note ? p->doctor_note : or_empty(p->why);
            nl(sb, &first);
            gh_sb_printf(sb, "<tr><td><strong>%s</strong></td><td>%s</td></tr>",
                         p->title, gh_report_clean_why(arena, note));
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</table>");
    }

    if (nreferrals) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h4>Specialist Referrals</h4>");
        for (size_t i = 0; i < nreferrals; i++) {
            const gh_referral *ref = &recs->referrals[i];
            const char *color = urgency_color(ref->urgency);
            nl(sb, &first);
            gh_sb_printf(sb, "<p><strong>%s</strong>: %s <span class=\"mag-badge\" "
                             "style=\"background:%s;color:#fff\">%s</span></p>",
                         ref->specialist, gh_report_clean_why(arena, ref->reason), color,
                         ref->urgency ? ref->urgency : "routine");
        }
    }

    if (stars && nstars) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h4>Pharmacogenomic Profile</h4>");
        nl(sb, &first);
        gh_sb_puts(sb, "<table><tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th></tr>");
        for (size_t i = 0; i < nstars; i++) {
            nl(sb, &first);
            gh_sb_printf(sb, "<tr><td><strong>%s</strong></td><td><code>%s</code></td><td>%s</td></tr>",
                         stars[i].gene, stars[i].diplotype,
                         py_title(arena, stars[i].phenotype, true));
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</table>");
    }

    if (apoe && strcmp(apoe->apoe_type, "Unknown") != 0) {
        const char *rl = or_empty(apoe->risk_level);
        const char *color =
            strcmp(rl, "reduced") == 0 ? "var(--green)" :
            strcmp(rl, "average") == 0 ? GH_C.blue :
            strcmp(rl, "moderate") == 0 ? GH_C.amber :
            strcmp(rl, "elevated") == 0 ? GH_C.orange :
            strcmp(rl, "high") == 0 ? GH_C.red : "var(--accent2)";
        nl(sb, &first);
        gh_sb_printf(sb, "<h4>APOE Status</h4><p><strong>%s</strong> \xe2\x80\x94 "
                         "<span class=\"mag-badge\" style=\"background:%s;color:#fff\">%s Risk</span></p>",
                     apoe->apoe_type, color, py_title(arena, rl, false));
    }

    if (nacmg) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h4>ACMG Medically Actionable Findings</h4>");
        nl(sb, &first);
        gh_sb_puts(sb, "<table><tr><th>Gene</th><th>Condition</th><th>Actionability</th></tr>");
        for (size_t i = 0; i < nacmg; i++) {
            const gh_cv_finding *f = acmg->findings[i].finding;
            nl(sb, &first);
            gh_sb_printf(sb, "<tr><td><strong>%s</strong></td><td>%s</td><td>%s</td></tr>",
                         f->gene ? f->gene : "Unknown",
                         gh_report_clean_condition(arena, f->traits && *f->traits ? f->traits : "Unknown"),
                         or_empty(acmg->findings[i].actionability));
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</table>");
        nl(sb, &first);
        gh_sb_puts(sb, "<div class=\"doctor-callout\" style=\"border-color:var(--warn)\">"
                       "ACMG actionable findings detected. Refer to genetic counselor.</div>");
    }

    if (recs && recs->nmonitoring) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h4>Recommended Monitoring</h4>");
        nl(sb, &first);
        gh_sb_puts(sb, "<table><tr><th>Test</th><th>Frequency</th></tr>");
        for (size_t i = 0; i < recs->nmonitoring; i++) {
            nl(sb, &first);
            gh_sb_printf(sb, "<tr><td>%s</td><td>%s</td></tr>",
                         recs->monitoring_schedule[i].test, recs->monitoring_schedule[i].frequency);
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</table>");
    }

    nl(sb, &first);
    gh_sb_puts(sb, "</div>");
}

/* ------------------------------------------------------------------ */
/* 12. References                                                      */
/* ------------------------------------------------------------------ */

void gh_report_references(gh_strbuf *sb, const gh_analysis *a)
{
    bool first = true;
    nl(sb, &first);
    gh_sb_puts(sb, "<h3>Key Papers</h3><ul>");
    const char *seen[128];
    size_t nseen = 0;
    for (size_t i = 0; i < GH_PAPER_REF_COUNT; i++) {
        const gh_paper_ref *r = &GH_PAPER_REFS[i];
        bool dup = false;
        for (size_t k = 0; k < nseen && !dup; k++)
            dup = strcmp(seen[k], r->pmid) == 0;
        if (dup)
            continue;
        if (nseen < 128)
            seen[nseen++] = r->pmid;
        nl(sb, &first);
        gh_sb_printf(sb, "<li><a href=\"https://pubmed.ncbi.nlm.nih.gov/%s/\" "
                         "target=\"_blank\" rel=\"noopener\">%s</a> (PMID: %s, %d)</li>",
                     r->pmid, r->title, r->pmid, r->year);
    }
    nl(sb, &first);
    gh_sb_puts(sb, "</ul>");

    /* sorted(set(rsid ...)): distinct, string order. */
    size_t nf = a ? a->nfindings : 0;
    const char **rsids = calloc(nf + 1, sizeof *rsids);
    size_t n = 0;
    for (size_t i = 0; i < nf; i++) {
        const char *r = a->findings[i].rsid;
        if (!r || strncmp(r, "rs", 2) != 0)
            continue;
        bool dup = false;
        for (size_t k = 0; k < n && !dup; k++)
            dup = strcmp(rsids[k], r) == 0;
        if (!dup)
            rsids[n++] = r;
    }
    qsort(rsids, n, sizeof *rsids, compare_strp);
    if (n) {
        nl(sb, &first);
        gh_sb_puts(sb, "<h3>Database Links for All Analyzed rsIDs</h3>");
        nl(sb, &first);
        gh_sb_puts(sb, "<div class=\"rsid-grid\">");
        for (size_t i = 0; i < n; i++) {
            nl(sb, &first);
            gh_sb_printf(sb, "<div><code>%s</code> ", rsids[i]);
            gh_report_db_links(sb, rsids[i]);
            gh_report_paper_refs(sb, rsids[i]);
            gh_sb_puts(sb, "</div>");
        }
        nl(sb, &first);
        gh_sb_puts(sb, "</div>");
    }
    free(rsids);

    nl(sb, &first);
    gh_sb_puts(sb, "<h3>Methodology &amp; Disclaimers</h3>");
    nl(sb, &first);
    gh_sb_puts(sb,
        "<ul>"
        "<li>Lifestyle findings from curated SNP database (~260 variants)</li>"
        "<li>Disease risk from ClinVar (~341K variants scanned)</li>"
        "<li>Drug interactions from PharmGKB clinical annotations</li>"
        "<li>Polygenic risk scores from published GWAS (8 conditions)</li>"
        "<li>Star allele calling for 6 pharmacogenes (CPIC-style)</li>"
        "<li>Only true SNPs analyzed (indels filtered to prevent false positives)</li>"
        "<li>This report is for <strong>informational purposes only</strong> \xe2\x80\x94 not a clinical diagnosis</li>"
        "<li>Genetic associations are probabilistic, not deterministic</li>"
        "<li>Always consult healthcare providers before making medical decisions</li>"
        "</ul>");
}

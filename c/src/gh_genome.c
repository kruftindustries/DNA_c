/* getline is POSIX.1-2008; request it explicitly under -std=c11. */
#define _POSIX_C_SOURCE 200809L

#include "gh_genome.h"
#include "gh_tsv.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Genotype validation                                                 */
/* ------------------------------------------------------------------ */

static bool is_base(char c)
{
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

bool gh_genotype_valid(const char *s)
{
    if (!s)
        return false;
    size_t n = strlen(s);
    if (n < 1 || n > 2)
        return false;
    for (size_t i = 0; i < n; i++)
        if (!is_base(s[i]))
            return false;
    return true;
}

/* ------------------------------------------------------------------ */
/* Genome loading                                                      */
/* ------------------------------------------------------------------ */

size_t gh_split_tabs(char *line, char **out, size_t max)
{
    size_t n = 0;
    char *p = line;
    while (n < max) {
        out[n++] = p;
        char *tab = strchr(p, '\t');
        if (!tab)
            break;
        *tab = '\0';
        p = tab + 1;
    }
    return n;
}

/* AncestryDNA numbers the non-autosomal chromosomes: 23 = X, 24 = Y, 25 = the
 * pseudoautosomal region (which 23andMe reports as X), 26 = MT. Mapped so
 * ClinVar matching and the quality metrics see one set of labels. */
static const char *ancestry_chromosome(const char *chrom)
{
    if (strcmp(chrom, "23") == 0 || strcmp(chrom, "25") == 0)
        return "X";
    if (strcmp(chrom, "24") == 0)
        return "Y";
    if (strcmp(chrom, "26") == 0)
        return "MT";
    return chrom;
}

gh_row_kind gh_genome_parse_row(char **parts, size_t nparts, gh_row *out)
{
    memset(out, 0, sizeof(*out));

    if (nparts >= 5 && parts[3][0] && !parts[3][1] &&
        parts[4][0] && !parts[4][1]) {
        out->rsid = parts[0];
        out->chromosome = ancestry_chromosome(parts[1]);
        out->position = parts[2];
        if (parts[3][0] == '0' || parts[4][0] == '0')
            return GH_ROW_NO_CALL;
        out->joined[0] = parts[3][0];
        out->joined[1] = parts[4][0];
        out->joined[2] = '\0';
        out->genotype = out->joined;
    } else if (nparts >= 4) {
        out->rsid = parts[0];
        out->chromosome = parts[1];
        out->position = parts[2];
        if (strcmp(parts[3], "--") == 0)
            return GH_ROW_NO_CALL;
        out->genotype = parts[3];
    } else {
        return GH_ROW_SHORT;
    }

    return gh_genotype_valid(out->genotype) ? GH_ROW_OK : GH_ROW_INVALID;
}

bool gh_genome_load(gh_genome *g, gh_arena *arena, const char *path,
                    const gh_map *pos_to_rsid)
{
    memset(g, 0, sizeof(*g));
    g->arena = arena;

    FILE *f = fopen(path, "rb");
    if (!f)
        return false;

    /* 23andMe exports run to ~600k rows; size the maps to avoid rehashing. */
    gh_map_init(&g->by_rsid, arena, 1 << 16);
    gh_map_init(&g->by_position, arena, 1 << 16);

    char *line = NULL;
    size_t cap = 0;
    long len;

    while ((len = gh_read_line(&line, &cap, f)) != -1) {
        if (len > 0 && line[len - 1] == '\n')
            line[--len] = '\0';
        if (len > 0 && line[len - 1] == '\r')
            line[--len] = '\0';
        if (len == 0 || line[0] == '#')
            continue;

        char *parts[GH_ROW_MAX_FIELDS];
        size_t nparts = gh_split_tabs(line, parts, GH_ROW_MAX_FIELDS);

        gh_row row;
        switch (gh_genome_parse_row(parts, nparts, &row)) {
        case GH_ROW_SHORT:
            continue;
        case GH_ROW_NO_CALL:
            g->no_calls++;
            continue;
        case GH_ROW_INVALID:
            g->skipped++;
            continue;
        case GH_ROW_OK:
            break;
        }

        const char *rsid = row.rsid;
        const char *chrom = row.chromosome;
        const char *pos = row.position;
        const char *genotype = row.genotype;

        /* chrom:pos key, built once and reused for both maps. */
        size_t klen = strlen(chrom) + 1 + strlen(pos);
        char *pos_key = gh_alloc(arena, klen + 1);
        snprintf(pos_key, klen + 1, "%s:%s", chrom, pos);

        /* WGS-derived files carry positional IDs; map them back to rsIDs. */
        if (strncmp(rsid, "rs", 2) != 0 && pos_to_rsid) {
            const char *mapped = gh_map_get(pos_to_rsid, pos_key);
            if (mapped) {
                rsid = mapped;
                g->resolved++;
            }
        }

        gh_variant *v = gh_alloc(arena, sizeof(*v));
        v->rsid = gh_strdup(arena, rsid);
        v->chromosome = gh_strdup(arena, chrom);
        v->position = gh_strdup(arena, pos);
        v->genotype = gh_strdup(arena, genotype);

        /* Later rows win, matching the Python dict assignment. */
        gh_map_put(&g->by_rsid, v->rsid, v);
        gh_map_put(&g->by_position, pos_key, v);
        g->total++;
    }

    free(line);
    fclose(f);
    return true;
}

const gh_variant *gh_genome_variant(const gh_genome *g, const char *rsid)
{
    return gh_map_get(&g->by_rsid, rsid);
}

const char *gh_genome_genotype(const gh_genome *g, const char *rsid)
{
    const gh_variant *v = gh_genome_variant(g, rsid);
    return v ? v->genotype : NULL;
}

/* ------------------------------------------------------------------ */
/* rsid_positions_grch37.json                                          */
/* ------------------------------------------------------------------ */

/* Minimal reader for the flat {"rsID": {"chrom": "..", "pos": ".."}} shape
 * this file always has. Not a general JSON parser: it understands strings
 * (with the standard escapes), objects and whitespace, which is all the file
 * contains. Anything unexpected aborts the load rather than guessing. */

typedef struct {
    const char *p;
    const char *end;
    gh_arena *arena;
} json_scan;

static void js_ws(json_scan *s)
{
    while (s->p < s->end &&
           (*s->p == ' ' || *s->p == '\t' || *s->p == '\n' || *s->p == '\r'))
        s->p++;
}

static bool js_eat(json_scan *s, char c)
{
    js_ws(s);
    if (s->p < s->end && *s->p == c) {
        s->p++;
        return true;
    }
    return false;
}

/* Parse a JSON string into the arena. Returns NULL on malformed input. */
static char *js_string(json_scan *s)
{
    js_ws(s);
    if (s->p >= s->end || *s->p != '"')
        return NULL;
    s->p++;

    const char *start = s->p;
    size_t extra = 0;
    while (s->p < s->end && *s->p != '"') {
        if (*s->p == '\\') {
            s->p++;
            extra++;
            if (s->p >= s->end)
                return NULL;
        }
        s->p++;
    }
    if (s->p >= s->end)
        return NULL;

    size_t raw = (size_t)(s->p - start);
    char *out = gh_alloc(s->arena, raw - extra + 1);
    char *w = out;
    for (const char *r = start; r < s->p; r++) {
        if (*r != '\\') {
            *w++ = *r;
            continue;
        }
        r++;
        switch (*r) {
        case 'n':  *w++ = '\n'; break;
        case 't':  *w++ = '\t'; break;
        case 'r':  *w++ = '\r'; break;
        case 'b':  *w++ = '\b'; break;
        case 'f':  *w++ = '\f'; break;
        case '/':  *w++ = '/';  break;
        case '\\': *w++ = '\\'; break;
        case '"':  *w++ = '"';  break;
        default:
            /* \uXXXX and anything else: keep the escape literally rather
             * than silently mangling it. These identifiers are ASCII. */
            *w++ = '\\';
            *w++ = *r;
            break;
        }
    }
    *w = '\0';
    s->p++; /* closing quote */
    return out;
}

/* Skip a value we do not care about (string, object, array, or literal). */
static bool js_skip_value(json_scan *s)
{
    js_ws(s);
    if (s->p >= s->end)
        return false;
    if (*s->p == '"')
        return js_string(s) != NULL;
    if (*s->p == '{' || *s->p == '[') {
        char open = *s->p, close = (open == '{') ? '}' : ']';
        int depth = 0;
        while (s->p < s->end) {
            if (*s->p == '"') {
                if (!js_string(s))
                    return false;
                continue;
            }
            if (*s->p == open)
                depth++;
            else if (*s->p == close && --depth == 0) {
                s->p++;
                return true;
            }
            s->p++;
        }
        return false;
    }
    while (s->p < s->end && *s->p != ',' && *s->p != '}' && *s->p != ']')
        s->p++;
    return true;
}

bool gh_load_rsid_positions(gh_map *out, gh_arena *arena, const char *path)
{
    gh_map_init(out, arena, 1024);

    FILE *f = fopen(path, "rb");
    if (!f)
        return false;

    if (fseek(f, 0, SEEK_END) != 0) {
        fclose(f);
        return false;
    }
    long end = ftell(f);
    if (end < 0) {
        fclose(f);
        return false;
    }
    rewind(f);

    size_t size = (size_t)end;
    char *buf = gh_alloc(arena, size + 1);
    size_t got = size ? fread(buf, 1, size, f) : 0;
    fclose(f);
    if (got != size)
        return false;
    buf[size] = '\0';

    json_scan s = {buf, buf + size, arena};
    if (!js_eat(&s, '{'))
        return false;

    if (js_eat(&s, '}'))
        return true; /* empty object */

    for (;;) {
        char *rsid = js_string(&s);
        if (!rsid || !js_eat(&s, ':'))
            return false;

        if (!js_eat(&s, '{'))
            return false;

        const char *chrom = NULL, *pos = NULL;
        if (!js_eat(&s, '}')) {
            for (;;) {
                char *key = js_string(&s);
                if (!key || !js_eat(&s, ':'))
                    return false;

                if (strcmp(key, "chrom") == 0) {
                    chrom = js_string(&s);
                    if (!chrom)
                        return false;
                } else if (strcmp(key, "pos") == 0) {
                    pos = js_string(&s);
                    if (!pos)
                        return false;
                } else if (!js_skip_value(&s)) {
                    return false;
                }

                if (js_eat(&s, ','))
                    continue;
                if (js_eat(&s, '}'))
                    break;
                return false;
            }
        }

        if (chrom && pos) {
            size_t klen = strlen(chrom) + 1 + strlen(pos);
            char *key = gh_alloc(arena, klen + 1);
            snprintf(key, klen + 1, "%s:%s", chrom, pos);
            gh_map_put(out, key, rsid);
        }

        if (js_eat(&s, ','))
            continue;
        if (js_eat(&s, '}'))
            break;
        return false;
    }

    return true;
}

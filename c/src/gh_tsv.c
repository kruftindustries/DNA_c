/* getline is POSIX.1-2008; request it explicitly under -std=c11. */
#define _POSIX_C_SOURCE 200809L

#include "gh_tsv.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GH_TSV_INIT_FIELDS 32

static void push_field(gh_tsv *t, char *field)
{
    if (t->nfields == t->fields_cap) {
        size_t cap = t->fields_cap ? t->fields_cap * 2 : GH_TSV_INIT_FIELDS;
        char **grown = gh_alloc(t->arena, cap * sizeof(*grown));
        if (t->nfields)
            memcpy(grown, t->fields, t->nfields * sizeof(*grown));
        t->fields = grown;
        t->fields_cap = cap;
    }
    t->fields[t->nfields++] = field;
}

bool gh_tsv_open(gh_tsv *t, gh_arena *arena, const char *path)
{
    memset(t, 0, sizeof(*t));
    t->arena = arena;
    t->path = gh_strdup(arena, path);

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

    t->buf = buf;
    t->size = size;
    t->pos = 0;
    return true;
}

bool gh_tsv_open_streaming(gh_tsv *t, gh_arena *arena, const char *path)
{
    memset(t, 0, sizeof(*t));
    t->arena = arena;
    t->path = gh_strdup(arena, path);

    FILE *f = fopen(path, "rb");
    if (!f)
        return false;
    t->stream = f;
    return true;
}

void gh_tsv_close(gh_tsv *t)
{
    if (t->stream) {
        fclose((FILE *)t->stream);
        t->stream = NULL;
    }
    free(t->line);
    t->line = NULL;
    t->line_cap = 0;
    free(t->record);
    t->record = NULL;
    t->record_cap = 0;
}

/* Append `n` bytes to the record buffer, growing it as needed. */
static bool record_append(gh_tsv *t, const char *data, size_t n, size_t *len)
{
    if (*len + n + 1 > t->record_cap) {
        size_t cap = t->record_cap ? t->record_cap : 4096;
        while (cap < *len + n + 1)
            cap *= 2;
        char *grown = realloc(t->record, cap);
        if (!grown)
            return false;
        t->record = grown;
        t->record_cap = cap;
    }
    memcpy(t->record + *len, data, n);
    *len += n;
    t->record[*len] = '\0';
    return true;
}

/* Read one logical record into t->buf, joining continuation lines when a
 * quoted field spans them. Returns false at end of input. */
static bool read_record(gh_tsv *t)
{
    FILE *f = (FILE *)t->stream;
    size_t len = 0;
    bool in_quotes = false;

    for (;;) {
        long got = gh_read_line(&t->line, &t->line_cap, f);
        if (got < 0)
            return len > 0;   /* trailing record with no newline */

        if (!record_append(t, t->line, (size_t)got, &len))
            gh_oom(len);

        /* A field stays open while the quote count on the record is odd. */
        for (long i = 0; i < got; i++)
            if (t->line[i] == '"')
                in_quotes = !in_quotes;

        if (!in_quotes)
            break;
    }

    /* Strip the record terminator; the parser handles the rest. */
    while (len && (t->record[len - 1] == '\n' || t->record[len - 1] == '\r'))
        t->record[--len] = '\0';

    t->buf = t->record;
    t->size = len;
    t->pos = 0;
    return true;
}

bool gh_tsv_next(gh_tsv *t)
{
    if (t->stream) {
        if (!read_record(t))
            return false;
        /* An empty final line is not a record. */
        if (t->size == 0 && feof((FILE *)t->stream))
            return false;
    } else if (t->pos >= t->size) {
        return false;
    }

    t->nfields = 0;
    t->row++;

    /* Fields are unescaped in place. That is safe because unescaping only
     * ever shrinks a field: "" collapses to " and the surrounding quotes go
     * away, so the write cursor never passes the read cursor. */
    for (;;) {
        char *start = t->buf + t->pos;
        char *w = start;
        bool quoted = (t->pos < t->size && *start == '"');

        if (quoted) {
            t->pos++; /* consume the opening quote */
            while (t->pos < t->size) {
                char c = t->buf[t->pos];
                if (c == '"') {
                    if (t->pos + 1 < t->size && t->buf[t->pos + 1] == '"') {
                        *w++ = '"';
                        t->pos += 2;
                        continue;
                    }
                    t->pos++; /* closing quote */
                    break;
                }
                *w++ = c;
                t->pos++;
            }
            /* Trailing text after the closing quote is literal, matching
             * Python's csv reader. */
        }

        while (t->pos < t->size) {
            char c = t->buf[t->pos];
            if (c == '\t' || c == '\n' || c == '\r')
                break;
            *w++ = c;
            t->pos++;
        }

        char terminator = (t->pos < t->size) ? t->buf[t->pos] : '\0';
        *w = '\0';
        push_field(t, start);

        if (terminator == '\t') {
            t->pos++;
            continue;
        }
        if (terminator == '\r') {
            t->pos++;
            if (t->pos < t->size && t->buf[t->pos] == '\n')
                t->pos++;
        } else if (terminator == '\n') {
            t->pos++;
        }
        break;
    }

    return true;
}

bool gh_tsv_read_header(gh_tsv *t)
{
    if (!gh_tsv_next(t))
        return false;

    t->nheaders = t->nfields;
    t->headers = gh_alloc(t->arena, t->nheaders * sizeof(*t->headers));
    for (size_t i = 0; i < t->nheaders; i++)
        t->headers[i] = gh_strdup(t->arena, t->fields[i]);

    gh_map_init(&t->col, t->arena, t->nheaders);
    for (size_t i = 0; i < t->nheaders; i++) {
        /* Store index+1 so that "absent" and "column 0" stay distinguishable
         * through the map's NULL-for-missing convention. */
        gh_map_put(&t->col, t->headers[i], (void *)(uintptr_t)(i + 1));
    }
    return true;
}

int gh_tsv_col(const gh_tsv *t, const char *name)
{
    void *v = gh_map_get(&t->col, name);
    return v ? (int)((uintptr_t)v - 1) : -1;
}

int gh_tsv_col_any(const gh_tsv *t, const char *const *names)
{
    for (size_t i = 0; names[i]; i++) {
        int c = gh_tsv_col(t, names[i]);
        if (c >= 0)
            return c;
    }
    return -1;
}

const char *gh_tsv_at(const gh_tsv *t, int index)
{
    if (index < 0 || (size_t)index >= t->nfields)
        return "";
    return t->fields[index];
}

const char *gh_tsv_get(const gh_tsv *t, const char *name)
{
    return gh_tsv_at(t, gh_tsv_col(t, name));
}

long gh_read_line(char **buf, size_t *cap, FILE *f)
{
    if (!*buf || *cap < 128) {
        char *nb = realloc(*buf, 4096);
        if (!nb)
            return -1;
        *buf = nb;
        *cap = 4096;
    }
    size_t len = 0;
    for (;;) {
        if (!fgets(*buf + len, (int)(*cap - len), f))
            return len ? (long)len : -1;
        len += strlen(*buf + len);
        if (len && (*buf)[len - 1] == '\n')
            return (long)len;
        if (len + 1 < *cap)          /* short read without a newline: EOF */
            return (long)len;
        char *nb = realloc(*buf, *cap * 2);
        if (!nb)
            return -1;
        *buf = nb;
        *cap *= 2;
    }
}

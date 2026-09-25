#include "gh_mem.h"

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GH_ALIGN 16
#define GH_DEFAULT_CHUNK (256 * 1024)

typedef struct gh_chunk {
    struct gh_chunk *next;
    size_t used;
    size_t cap;
    /* Payload follows the header; `data` marks its aligned start. */
    _Alignas(GH_ALIGN) unsigned char data[];
} gh_chunk;

struct gh_arena {
    gh_chunk *head;
    size_t chunk_size;
    size_t used;
};

void gh_oom(size_t requested)
{
    fprintf(stderr, "fatal: out of memory requesting %zu bytes\n", requested);
    exit(1);
}

static size_t align_up(size_t n)
{
    return (n + (GH_ALIGN - 1)) & ~(size_t)(GH_ALIGN - 1);
}

static gh_chunk *chunk_new(size_t payload)
{
    gh_chunk *c = malloc(sizeof(*c) + payload);
    if (!c)
        gh_oom(payload);
    c->next = NULL;
    c->used = 0;
    c->cap = payload;
    return c;
}

gh_arena *gh_arena_new(size_t chunk_size)
{
    gh_arena *a = calloc(1, sizeof(*a));
    if (!a)
        return NULL;
    a->chunk_size = chunk_size ? chunk_size : GH_DEFAULT_CHUNK;
    a->head = chunk_new(a->chunk_size);
    return a;
}

void gh_arena_free(gh_arena *a)
{
    if (!a)
        return;
    gh_chunk *c = a->head;
    while (c) {
        gh_chunk *next = c->next;
        free(c);
        c = next;
    }
    free(a);
}

void *gh_alloc(gh_arena *a, size_t n)
{
    if (n == 0)
        n = 1;
    size_t need = align_up(n);

    if (a->head->used + need > a->head->cap) {
        /* Oversized requests get a dedicated chunk so they do not strand the
         * remaining space in the current one. */
        size_t payload = need > a->chunk_size ? need : a->chunk_size;
        gh_chunk *c = chunk_new(payload);
        c->next = a->head;
        a->head = c;
    }

    void *p = a->head->data + a->head->used;
    a->head->used += need;
    a->used += need;
    return p;
}

void *gh_calloc(gh_arena *a, size_t count, size_t size)
{
    if (count != 0 && size > SIZE_MAX / count)
        gh_oom(SIZE_MAX);
    size_t total = count * size;
    void *p = gh_alloc(a, total);
    memset(p, 0, total);
    return p;
}

char *gh_strndup(gh_arena *a, const char *s, size_t n)
{
    char *p = gh_alloc(a, n + 1);
    if (n)
        memcpy(p, s, n);
    p[n] = '\0';
    return p;
}

char *gh_strdup(gh_arena *a, const char *s)
{
    if (!s)
        return NULL;
    return gh_strndup(a, s, strlen(s));
}

size_t gh_arena_used(const gh_arena *a)
{
    return a ? a->used : 0;
}

/* ------------------------------------------------------------------ */
/* gh_strbuf                                                           */
/* ------------------------------------------------------------------ */

void gh_sb_init(gh_strbuf *sb)
{
    sb->data = NULL;
    sb->len = 0;
    sb->cap = 0;
}

void gh_sb_free(gh_strbuf *sb)
{
    free(sb->data);
    gh_sb_init(sb);
}

void gh_sb_reserve(gh_strbuf *sb, size_t extra)
{
    /* +1 keeps room for the terminator. */
    if (sb->len + extra + 1 <= sb->cap)
        return;
    size_t cap = sb->cap ? sb->cap : 256;
    while (cap < sb->len + extra + 1) {
        if (cap > SIZE_MAX / 2)
            gh_oom(extra);
        cap *= 2;
    }
    char *p = realloc(sb->data, cap);
    if (!p)
        gh_oom(cap);
    sb->data = p;
    sb->cap = cap;
}

void gh_sb_write(gh_strbuf *sb, const char *s, size_t n)
{
    if (!n)
        return;
    gh_sb_reserve(sb, n);
    memcpy(sb->data + sb->len, s, n);
    sb->len += n;
    sb->data[sb->len] = '\0';
}

void gh_sb_putc(gh_strbuf *sb, char c)
{
    gh_sb_write(sb, &c, 1);
}

void gh_sb_puts(gh_strbuf *sb, const char *s)
{
    if (s)
        gh_sb_write(sb, s, strlen(s));
}

void gh_sb_printf(gh_strbuf *sb, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    va_list probe;
    va_copy(probe, ap);
    int n = vsnprintf(NULL, 0, fmt, probe);
    va_end(probe);

    if (n < 0) {
        va_end(ap);
        return;
    }
    gh_sb_reserve(sb, (size_t)n);
    vsnprintf(sb->data + sb->len, (size_t)n + 1, fmt, ap);
    sb->len += (size_t)n;
    va_end(ap);
}

void gh_sb_put_escaped(gh_strbuf *sb, const char *s)
{
    if (!s)
        return;
    for (const char *p = s; *p; p++) {
        switch (*p) {
        case '&':  gh_sb_puts(sb, "&amp;");  break;
        case '<':  gh_sb_puts(sb, "&lt;");   break;
        case '>':  gh_sb_puts(sb, "&gt;");   break;
        case '"':  gh_sb_puts(sb, "&quot;"); break;
        /* html.escape's spelling, so escaped text is byte-identical. */
        case '\'': gh_sb_puts(sb, "&#x27;"); break;
        default:   gh_sb_putc(sb, *p);       break;
        }
    }
}

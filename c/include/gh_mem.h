/* Arena allocator and growable buffers.
 *
 * Everything the pipeline parses (genome rows, TSV fields, report fragments)
 * lives for the whole run, so allocations come from an arena and are released
 * in one go at exit. That removes almost all of the lifetime bookkeeping a
 * port of this shape would otherwise need.
 */
#ifndef GH_MEM_H
#define GH_MEM_H

#include <stddef.h>

typedef struct gh_arena gh_arena;

/* Create an arena. `chunk_size` is a hint; larger requests get their own
 * chunk. Returns NULL only if the initial chunk cannot be allocated. */
gh_arena *gh_arena_new(size_t chunk_size);

/* Release the arena and every pointer handed out by it. Safe on NULL. */
void gh_arena_free(gh_arena *a);

/* Allocate `n` bytes aligned for any type. Never returns NULL: allocation
 * failure is fatal (see gh_oom), because there is no useful recovery in a
 * batch analysis run. */
void *gh_alloc(gh_arena *a, size_t n);

/* Allocate zeroed memory for `count` items of `size` bytes, checking the
 * multiplication for overflow. */
void *gh_calloc(gh_arena *a, size_t count, size_t size);

/* Copy a NUL-terminated string into the arena. NULL in, NULL out. */
char *gh_strdup(gh_arena *a, const char *s);

/* Copy `n` bytes into the arena and NUL-terminate. */
char *gh_strndup(gh_arena *a, const char *s, size_t n);

/* Bytes currently handed out, for diagnostics. */
size_t gh_arena_used(const gh_arena *a);

/* Report a fatal allocation failure and exit. Does not return. */
void gh_oom(size_t requested);

/* ------------------------------------------------------------------ */
/* Growable string buffer, used to assemble the HTML report.           */
/* ------------------------------------------------------------------ */

typedef struct {
    char *data;    /* NUL-terminated; owned by malloc, not the arena */
    size_t len;
    size_t cap;
} gh_strbuf;

void gh_sb_init(gh_strbuf *sb);
void gh_sb_free(gh_strbuf *sb);
void gh_sb_reserve(gh_strbuf *sb, size_t extra);
void gh_sb_putc(gh_strbuf *sb, char c);
void gh_sb_write(gh_strbuf *sb, const char *s, size_t n);
void gh_sb_puts(gh_strbuf *sb, const char *s);

/* printf-style append. */
void gh_sb_printf(gh_strbuf *sb, const char *fmt, ...)
#if defined(__GNUC__)
    __attribute__((format(printf, 2, 3)))
#endif
    ;

/* Append `s` with &, <, >, " and ' escaped, for interpolation into HTML.
 * NULL is treated as the empty string. */
void gh_sb_put_escaped(gh_strbuf *sb, const char *s);

#endif /* GH_MEM_H */

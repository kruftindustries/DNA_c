/* Tab-separated value reader.
 *
 * Matches the quoting rules of Python's `csv` module with delimiter='\t' and
 * the default QUOTE_MINIMAL dialect, because the reference implementation
 * reads these files with csv.DictReader:
 *
 *   - A double quote opens a quoted field only at the start of that field.
 *   - Inside a quoted field, "" is a literal quote, and tabs and newlines are
 *     ordinary characters.
 *   - Text after the closing quote is appended literally.
 *   - \r\n and \n both end a record; \r\n inside quotes is preserved as-is.
 *
 * The ClinPGx allele table relies on this: 21 of its annotation texts are
 * quoted fields with embedded tabs, which a naive split on '\t' would truncate.
 */
#ifndef GH_TSV_H
#define GH_TSV_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#include "gh_map.h"
#include "gh_mem.h"

typedef struct {
    /* Whole-file mode: the file is slurped into the arena and parsed in
     * place. Streaming mode leaves `buf` holding one record at a time. */
    char *buf;
    size_t size;
    size_t pos;

    void *stream;        /* FILE* in streaming mode, NULL otherwise */
    char *line;          /* getline buffer, malloc-owned */
    size_t line_cap;
    char *record;        /* assembled record, malloc-owned */
    size_t record_cap;

    char **fields;       /* current record */
    size_t nfields;
    size_t fields_cap;

    char **headers;      /* header record, if gh_tsv_read_header ran */
    size_t nheaders;
    gh_map col;          /* header name -> (index + 1), as a tagged pointer */

    gh_arena *arena;
    const char *path;
    size_t row;          /* 1-based record number, header included */
} gh_tsv;

/* Open and slurp a TSV. Returns false if the file cannot be read.
 * Fields stay valid for the life of the arena. */
bool gh_tsv_open(gh_tsv *t, gh_arena *arena, const char *path);

/* Open a TSV for streaming, holding one record in memory at a time.
 * Use for files too large to slurp -- ClinVar runs to hundreds of megabytes.
 * Fields are only valid until the next gh_tsv_next call, so anything kept
 * must be copied. Call gh_tsv_close when finished. */
bool gh_tsv_open_streaming(gh_tsv *t, gh_arena *arena, const char *path);

/* Release the streaming reader's file handle and buffers. Safe to call on a
 * whole-file reader, and safe to call twice. */
void gh_tsv_close(gh_tsv *t);

/* Read the next record. Returns false at end of input. Fields are valid
 * until the following call. */
bool gh_tsv_next(gh_tsv *t);

/* Read the first record as a header and build the column index. */
bool gh_tsv_read_header(gh_tsv *t);

/* Column index for a header name, or -1 if the file has no such column. */
int gh_tsv_col(const gh_tsv *t, const char *name);

/* Field by index; "" when the index is past the end of a short record, so
 * callers never have to bounds-check. Never returns NULL. */
const char *gh_tsv_at(const gh_tsv *t, int index);

/* Field by header name; "" when the column is missing or the record is short. */
const char *gh_tsv_get(const gh_tsv *t, const char *name);

/* First present column among `names` (NULL-terminated), for files whose
 * column names have been renamed upstream. Returns -1 if none match. */
int gh_tsv_col_any(const gh_tsv *t, const char *const *names);

/* Read one line (including its newline, if any) into a malloc'd buffer that
 * grows as needed, the way POSIX getline does. Returns the length, or -1 at
 * end of file or on allocation failure. Portable: getline is not in the
 * Windows C runtime. */
long gh_read_line(char **buf, size_t *cap, FILE *f);

#endif /* GH_TSV_H */

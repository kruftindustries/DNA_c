/* Unit tests for the arena, string buffer, hash map and TSV reader. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gh_genome.h"
#include "gh_map.h"
#include "gh_mem.h"
#include "gh_tsv.h"

static int failures;
static int checks;

static void check(int cond, const char *what)
{
    checks++;
    if (!cond) {
        failures++;
        fprintf(stderr, "  FAIL: %s\n", what);
    }
}

static void check_str(const char *got, const char *want, const char *what)
{
    checks++;
    if (!got || strcmp(got, want) != 0) {
        failures++;
        fprintf(stderr, "  FAIL: %s\n    got  %s\n    want %s\n",
                what, got ? got : "(null)", want);
    }
}

static void test_arena(void)
{
    gh_arena *a = gh_arena_new(64);
    check(a != NULL, "arena_new");

    char *s = gh_strdup(a, "hello");
    check_str(s, "hello", "strdup");

    /* Force several oversized allocations through their own chunks. */
    for (int i = 0; i < 10; i++) {
        char *big = gh_alloc(a, 4096);
        memset(big, 'x', 4096);
        check(big[0] == 'x' && big[4095] == 'x', "oversized alloc writable");
    }
    check_str(s, "hello", "earlier pointer survives new chunks");

    int *zeros = gh_calloc(a, 100, sizeof(int));
    int all_zero = 1;
    for (int i = 0; i < 100; i++)
        if (zeros[i] != 0)
            all_zero = 0;
    check(all_zero, "calloc zeroes");

    char *part = gh_strndup(a, "abcdef", 3);
    check_str(part, "abc", "strndup truncates");
    check(gh_strdup(a, NULL) == NULL, "strdup(NULL)");

    gh_arena_free(a);
}

static void test_strbuf(void)
{
    gh_strbuf sb;
    gh_sb_init(&sb);

    gh_sb_puts(&sb, "a");
    gh_sb_putc(&sb, 'b');
    gh_sb_printf(&sb, "%d-%s", 42, "c");
    check_str(sb.data, "ab42-c", "strbuf append");

    /* Grow well past the initial capacity. */
    for (int i = 0; i < 5000; i++)
        gh_sb_puts(&sb, "0123456789");
    check(sb.len == 6 + 50000, "strbuf length after growth");
    check(sb.data[sb.len] == '\0', "strbuf stays terminated");

    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_sb_put_escaped(&sb, "<a href=\"x\">R&D 'q'</a>");
    check_str(sb.data,
              "&lt;a href=&quot;x&quot;&gt;R&amp;D &#x27;q&#x27;&lt;/a&gt;",
              "html escaping");
    gh_sb_free(&sb);

    gh_sb_init(&sb);
    gh_sb_put_escaped(&sb, NULL);
    check(sb.len == 0, "escaping NULL is empty");
    gh_sb_free(&sb);
}

static void test_map(void)
{
    gh_arena *a = gh_arena_new(1024);
    gh_map m;
    gh_map_init(&m, a, 0);

    check(gh_map_get(&m, "missing") == NULL, "get on empty map");
    check(gh_map_put(&m, "k", (void *)1), "put returns true when new");
    check(!gh_map_put(&m, "k", (void *)2), "put returns false when replacing");
    check(gh_map_get(&m, "k") == (void *)2, "value overwritten");
    check(gh_map_len(&m) == 1, "len after overwrite");

    /* Enough entries to force several rehashes. */
    enum { N = 5000 };
    for (long i = 0; i < N; i++) {
        char key[32];
        snprintf(key, sizeof(key), "rs%ld", i);
        gh_map_put(&m, key, (void *)(i + 100));
    }
    check(gh_map_len(&m) == N + 1, "len after bulk insert");

    int all_found = 1;
    for (long i = 0; i < N; i++) {
        char key[32];
        snprintf(key, sizeof(key), "rs%ld", i);
        if (gh_map_get(&m, key) != (void *)(i + 100))
            all_found = 0;
    }
    check(all_found, "all keys retrievable after growth");
    check(gh_map_get(&m, "rs999999") == NULL, "absent key after growth");

    /* Keys are copied, so a mutated source buffer must not matter. */
    char temp[16];
    snprintf(temp, sizeof(temp), "borrowed");
    gh_map_put(&m, temp, (void *)7);
    memset(temp, 'z', sizeof(temp) - 1);
    temp[sizeof(temp) - 1] = '\0';
    check(gh_map_get(&m, "borrowed") == (void *)7, "keys are copied not borrowed");

    /* getn respects the given length rather than reading to NUL. */
    check(gh_map_getn(&m, "borrowedXXX", 8) == (void *)7, "getn honours length");

    size_t iter = 0, counted = 0;
    const char *k;
    void *v;
    while (gh_map_next(&m, &iter, &k, &v))
        counted++;
    check(counted == gh_map_len(&m), "iteration visits every entry");

    gh_arena_free(a);
}

/* Write a temporary file and return its path in static storage. */
static const char *write_temp(const char *name, const char *content)
{
    static char path[256];
    snprintf(path, sizeof(path), "/tmp/gh_test_%s", name);
    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "  FAIL: cannot write %s\n", path);
        failures++;
        return path;
    }
    fwrite(content, 1, strlen(content), f);
    fclose(f);
    return path;
}

static void test_tsv_basic(void)
{
    gh_arena *a = gh_arena_new(4096);
    const char *path = write_temp("basic.tsv",
        "rsid\tchrom\tpos\tgenotype\n"
        "rs1\t1\t100\tAA\n"
        "rs2\t2\t200\tCT\n");

    gh_tsv t;
    check(gh_tsv_open(&t, a, path), "tsv open");
    check(gh_tsv_read_header(&t), "read header");
    check(gh_tsv_col(&t, "rsid") == 0, "column index rsid");
    check(gh_tsv_col(&t, "genotype") == 3, "column index genotype");
    check(gh_tsv_col(&t, "nope") == -1, "missing column is -1");

    check(gh_tsv_next(&t), "first row");
    check_str(gh_tsv_get(&t, "rsid"), "rs1", "row1 rsid");
    check_str(gh_tsv_get(&t, "genotype"), "AA", "row1 genotype");
    check_str(gh_tsv_get(&t, "nope"), "", "missing column reads empty");

    check(gh_tsv_next(&t), "second row");
    check_str(gh_tsv_get(&t, "genotype"), "CT", "row2 genotype");
    check(!gh_tsv_next(&t), "end of file");

    gh_arena_free(a);
}

static void test_tsv_quoting(void)
{
    gh_arena *a = gh_arena_new(4096);
    /* Mirrors the shapes found in the ClinPGx allele table. */
    const char *path = write_temp("quoted.tsv",
        "id\ttext\ttail\n"
        "a\t\"has\ttab\"\tz\n"
        "b\t\"has\nnewline\"\tz\n"
        "c\t\"doubled \"\"quote\"\"\"\tz\n"
        "d\t\"closed\" then literal\tz\n"
        "e\tplain \"mid\" quote\tz\n");

    gh_tsv t;
    check(gh_tsv_open(&t, a, path), "open quoted");
    check(gh_tsv_read_header(&t), "header quoted");

    check(gh_tsv_next(&t), "row a");
    check_str(gh_tsv_get(&t, "text"), "has\ttab", "embedded tab kept");
    check_str(gh_tsv_get(&t, "tail"), "z", "field after embedded tab");

    check(gh_tsv_next(&t), "row b");
    check_str(gh_tsv_get(&t, "text"), "has\nnewline", "embedded newline kept");
    check_str(gh_tsv_get(&t, "tail"), "z", "field after embedded newline");

    check(gh_tsv_next(&t), "row c");
    check_str(gh_tsv_get(&t, "text"), "doubled \"quote\"", "doubled quotes unescaped");

    check(gh_tsv_next(&t), "row d");
    check_str(gh_tsv_get(&t, "text"), "closed then literal", "text after closing quote");

    check(gh_tsv_next(&t), "row e");
    check_str(gh_tsv_get(&t, "text"), "plain \"mid\" quote", "mid-field quote is literal");

    check(!gh_tsv_next(&t), "eof quoted");
    gh_arena_free(a);
}

static void test_tsv_edges(void)
{
    gh_arena *a = gh_arena_new(4096);

    gh_tsv t;
    check(!gh_tsv_open(&t, a, "/nonexistent/path.tsv"), "open missing file fails");

    const char *empty = write_temp("empty.tsv", "");
    check(gh_tsv_open(&t, a, empty), "open empty file");
    check(!gh_tsv_read_header(&t), "no header in empty file");

    /* Short rows, a trailing empty field, and no final newline. */
    const char *ragged = write_temp("ragged.tsv", "a\tb\tc\nx\ny\t\nz\t1\t2");
    check(gh_tsv_open(&t, a, ragged), "open ragged");
    check(gh_tsv_read_header(&t), "header ragged");

    check(gh_tsv_next(&t), "ragged row 1");
    check_str(gh_tsv_get(&t, "a"), "x", "short row field present");
    check_str(gh_tsv_get(&t, "c"), "", "short row missing field is empty");

    check(gh_tsv_next(&t), "ragged row 2");
    check_str(gh_tsv_get(&t, "b"), "", "explicit empty field");

    check(gh_tsv_next(&t), "ragged row 3");
    check_str(gh_tsv_get(&t, "c"), "2", "final row without trailing newline");
    check(!gh_tsv_next(&t), "eof ragged");

    /* CRLF line endings. */
    const char *crlf = write_temp("crlf.tsv", "a\tb\r\n1\t2\r\n");
    check(gh_tsv_open(&t, a, crlf), "open crlf");
    check(gh_tsv_read_header(&t), "header crlf");
    check_str(t.headers[1], "b", "crlf header not polluted by \\r");
    check(gh_tsv_next(&t), "crlf row");
    check_str(gh_tsv_get(&t, "b"), "2", "crlf field clean");

    /* col_any picks the first header that exists. */
    const char *renamed = write_temp("renamed.tsv", "Summary Annotation ID\tx\nv\t1\n");
    check(gh_tsv_open(&t, a, renamed), "open renamed");
    check(gh_tsv_read_header(&t), "header renamed");
    static const char *const ids[] = {"Clinical Annotation ID", "Summary Annotation ID", NULL};
    check(gh_tsv_col_any(&t, ids) == 0, "col_any finds the second name");
    static const char *const none[] = {"Nope", NULL};
    check(gh_tsv_col_any(&t, none) == -1, "col_any returns -1 when none match");

    gh_arena_free(a);
}

/* ------------------------------------------------------------------ */
/* Genome row classification                                           */
/* ------------------------------------------------------------------ */

static gh_row_kind classify(const char *text, gh_row *row)
{
    char line[128];
    snprintf(line, sizeof(line), "%s", text);
    char *parts[GH_ROW_MAX_FIELDS];
    size_t n = gh_split_tabs(line, parts, GH_ROW_MAX_FIELDS);
    gh_row_kind kind = gh_genome_parse_row(parts, n, row);
    /* The row points into `line`, so copy out what the checks need. */
    static char chrom[16], geno[8];
    if (row->chromosome) { snprintf(chrom, sizeof chrom, "%s", row->chromosome); row->chromosome = chrom; }
    if (row->genotype && row->genotype != row->joined) { snprintf(geno, sizeof geno, "%s", row->genotype); row->genotype = geno; }
    return kind;
}

static void test_genome_rows(void)
{
    gh_row row;

    /* 23andMe layout, unchanged. */
    check(classify("rs1\t1\t100\tAG", &row) == GH_ROW_OK, "23andMe diploid");
    check_str(row.genotype, "AG", "23andMe genotype");
    check(classify("rs1\t1\t100\t--", &row) == GH_ROW_NO_CALL, "23andMe no-call");
    check(classify("rs1\t1\t100\tXY", &row) == GH_ROW_INVALID, "23andMe malformed");
    check(classify("rs1\t1\t100", &row) == GH_ROW_SHORT, "three fields");

    /* AncestryDNA layout: the fix. A heterozygote used to come back as its
     * first allele alone. */
    check(classify("rs429358\t19\t45411941\tC\tT", &row) == GH_ROW_OK,
          "AncestryDNA het accepted");
    check_str(row.genotype, "CT", "AncestryDNA het joined from both columns");
    check(classify("rs1\t1\t100\t0\t0", &row) == GH_ROW_NO_CALL, "0 0 no-call");
    check(classify("rs1\t1\t100\tA\t0", &row) == GH_ROW_NO_CALL, "A 0 no-call");
    check(classify("rs1\t1\t100\tN\tN", &row) == GH_ROW_INVALID, "N N malformed");

    /* Chromosome codes. */
    classify("rs1\t23\t1\tA\tA", &row); check_str(row.chromosome, "X", "23 -> X");
    classify("rs1\t24\t1\tA\tA", &row); check_str(row.chromosome, "Y", "24 -> Y");
    classify("rs1\t25\t1\tA\tA", &row); check_str(row.chromosome, "X", "25 (PAR) -> X");
    classify("rs1\t26\t1\tA\tA", &row); check_str(row.chromosome, "MT", "26 -> MT");
    classify("rs1\t22\t1\tA\tA", &row); check_str(row.chromosome, "22", "autosome untouched");
    /* The mapping is a property of the two-column layout only. */
    classify("rs1\t23\t1\tAA", &row); check_str(row.chromosome, "23", "4-column 23 stays 23");

    /* The vendor header is not '#'-prefixed and must be malformed, not a
     * variant -- and must not be mistaken for the two-column layout. */
    check(classify("rsid\tchromosome\tposition\tallele1\tallele2", &row)
          == GH_ROW_INVALID, "AncestryDNA header is skipped");
}

int main(void)
{
    test_arena();
    test_strbuf();
    test_map();
    test_tsv_basic();
    test_tsv_quoting();
    test_tsv_edges();
    test_genome_rows();

    printf("%s: %d checks, %d failures\n",
           failures ? "FAILED" : "ok", checks, failures);
    return failures ? 1 : 0;
}

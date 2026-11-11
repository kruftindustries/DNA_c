/* Unit tests for the ClinPGx loader and the drug-gene matching it feeds.
 *
 * The real tables spell an indel annotation's genotypes with the allele
 * sequences and "del" (CFTR F508del: CTT/CTT, CTT/del, del/del) while an
 * array export writes the same calls as II, DI, DD. The loader derives the
 * array spelling so the two meet; these tests pin that mapping and the
 * cases it must leave alone.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gh_analysis.h"
#include "gh_genome.h"
#include "gh_mem.h"
#include "gh_pharmgkb.h"

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

static const char *write_file(const char *name, const char *text)
{
    static char paths[8][256];
    static int slot;
    char *path = paths[slot++ & 7];
    snprintf(path, 256, "/tmp/gh_pgx_%s", name);
    FILE *f = fopen(path, "wb");
    if (!f) {
        failures++;
        return path;
    }
    fputs(text, f);
    fclose(f);
    return path;
}

/* Annotation IDs 1-5 cover the forms the real tables use. */
static const char *ANNOTATIONS =
    "Clinical Annotation ID\tVariant/Haplotypes\tGene\tLevel of Evidence\tLevel Override\t"
    "Level Modifiers\tScore\tPMID Count\tEvidence Count\tDrug(s)\tPhenotype(s)\tPhenotype Category\n"
    "7\trs113993960\tCFTR\t3\t\t\t\t\t\tcysteamine\tCystic Fibrosis\tEfficacy\n"
    "1\trs113993960\tCFTR\t1A\t\t\t\t\t\tivacaftor\tCystic Fibrosis\tEfficacy\n"
    "2\trs72549303\tDPYD\t1A\t\t\t\t\t\tfluorouracil\t\tToxicity\n"
    "3\trs11322783\tIFNL4\t1B\t\t\t\t\t\tpeginterferon alfa-2a\t\tEfficacy\n"
    "4\trs8175347\tUGT1A1\t1A\t\t\t\t\t\tirinotecan\t\tToxicity\n"
    "5\trs4244285\tCYP2C19\t1A\t\t\t\t\t\tclopidogrel\t\tEfficacy\n"
    "6\tCYP2D6*4\tCYP2D6\t1A\t\t\t\t\t\tcodeine\t\tEfficacy\n";

static const char *ALLELES =
    "Clinical Annotation ID\tGenotype/Allele\tAnnotation Text\tAllele Function\n"
    /* a level-3 annotation on the same variant comes first in the file */
    "7\tCTT/CTT\tcysteamine: no F508del\t\n"
    "7\tCTT/del\tcysteamine: F508del carrier\t\n"
    /* sequence/del spelling, the common indel form */
    "1\tCTT/CTT\tno F508del\t\n"
    "1\tCTT/del\tone copy of F508del: ivacaftor indicated\t\n"
    "1\tdel/del\ttwo copies of F508del\t\n"
    /* a single-base deletion: the homozygous sequence form has no slash */
    "2\tGG\tnormal DPYD\t\n"
    "2\tG/del\tDPYD c.1898delC carrier\t\n"
    "2\tdel/del\tDPYD deficient\t\n"
    /* two sequences and no del: D/I cannot say which is which */
    "3\tTT/TT\tfavourable\t\n"
    "3\tG/TT\tintermediate\t\n"
    "3\tGG\tunfavourable\t\n"
    /* repeat alleles: likewise undecidable from D/I */
    "4\t(TA)6/(TA)6\tnormal\t\n"
    "4\t(TA)6/(TA)7\tintermediate\t\n"
    "4\t(TA)7/(TA)7\t*28/*28\t\n"
    /* an ordinary SNV */
    "5\tAA\tpoor metabolizer\t\n"
    "5\tAG\tintermediate metabolizer\t\n"
    "5\tGG\tnormal metabolizer\t\n"
    "6\t*4/*4\tpoor\t\n";

/* The entry for `rsid` at evidence `level`, from the per-annotation list. */
static const gh_pgx_entry *find_level(const gh_pharmgkb *p, const char *rsid, const char *level)
{
    for (size_t i = 0; i < p->nordered; i++)
        if (strcmp(p->ordered[i]->rsid, rsid) == 0 && strcmp(p->ordered[i]->level, level) == 0)
            return p->ordered[i];
    return NULL;
}

static void test_indel_keys_derived(void)
{
    gh_arena *a = gh_arena_new(1 << 18);
    gh_pharmgkb p;
    check(gh_pharmgkb_load(&p, a, write_file("ann.tsv", ANNOTATIONS),
                           write_file("all.tsv", ALLELES)),
          "tables load");
    check(p.annotations == 6, "haplotype annotation skipped, six rsID ones kept");

    /* One entry per annotation: the level-3 one that comes first in the
     * file does not swallow the 1A one. */
    check(p.nordered == 6, "six entries, one per annotation with allele rows");
    const gh_pgx_entry *first = gh_pharmgkb_find(&p, "rs113993960");
    check(first && strcmp(first->level, "3") == 0, "by_rsid gives the first annotation");
    check_str(gh_map_get(&first->genotypes, "DI"), "cysteamine: F508del carrier", "level-3 entry has its own texts and D/I keys");
    const gh_pgx_entry *cftr = find_level(&p, "rs113993960", "1A");
    check(cftr != NULL, "CFTR 1A entry exists");
    if (cftr) {
        check_str(gh_map_get(&cftr->genotypes, "CTT/del"),
                  "one copy of F508del: ivacaftor indicated", "original spelling kept");
        check_str(gh_map_get(&cftr->genotypes, "DI"),
                  "one copy of F508del: ivacaftor indicated", "CTT/del -> DI");
        check_str(gh_map_get(&cftr->genotypes, "DD"), "two copies of F508del", "del/del -> DD");
        check_str(gh_map_get(&cftr->genotypes, "II"), "no F508del", "CTT/CTT -> II");
    }

    const gh_pgx_entry *dpyd = gh_pharmgkb_find(&p, "rs72549303");
    check(dpyd != NULL, "DPYD entry exists");
    if (dpyd) {
        check_str(gh_map_get(&dpyd->genotypes, "II"), "normal DPYD", "GG (no slash) -> II when the entry has a del form");
        check_str(gh_map_get(&dpyd->genotypes, "DI"), "DPYD c.1898delC carrier", "G/del -> DI");
        check_str(gh_map_get(&dpyd->genotypes, "DD"), "DPYD deficient", "del/del -> DD");
        check_str(gh_map_get(&dpyd->genotypes, "GG"), "normal DPYD", "GG itself still there");
    }

    const gh_pgx_entry *ifnl4 = gh_pharmgkb_find(&p, "rs11322783");
    check(ifnl4 && !gh_map_has(&ifnl4->genotypes, "DI") && !gh_map_has(&ifnl4->genotypes, "II"),
          "two sequence alleles and no del: nothing derived");
    const gh_pgx_entry *ugt = gh_pharmgkb_find(&p, "rs8175347");
    check(ugt && !gh_map_has(&ugt->genotypes, "DI") && !gh_map_has(&ugt->genotypes, "II")
              && gh_map_len(&ugt->genotypes) == 3,
          "repeat alleles: nothing derived");
    const gh_pgx_entry *cyp = gh_pharmgkb_find(&p, "rs4244285");
    check(cyp && gh_map_len(&cyp->genotypes) == 3 && !gh_map_has(&cyp->genotypes, "II"),
          "an SNV entry gains nothing");

    gh_arena_free(a);
}

static void test_indel_call_reaches_drug_guide(void)
{
    gh_arena *a = gh_arena_new(1 << 18);
    gh_pharmgkb p;
    check(gh_pharmgkb_load(&p, a, write_file("ann2.tsv", ANNOTATIONS),
                           write_file("all2.tsv", ALLELES)),
          "tables load");

    const char *genome = write_file("genome.txt",
        "# 23andMe-style export with indel calls\n"
        "rsid\tchromosome\tposition\tgenotype\n"
        "rs113993960\t7\t117199645\tDI\n"     /* F508del carrier */
        "rs72549303\t1\t97915621\tII\n"       /* DPYD normal: reported, level 1A */
        "rs11322783\t19\t39739154\tDI\n"      /* undecidable: not reported */
        "rs4244285\t10\t96541616\tAG\n");     /* ordinary SNV, as before */
    gh_genome g;
    check(gh_genome_load(&g, a, genome, NULL), "genome loads");
    check(g.total == 4 && g.skipped == 1, "indel calls are kept (the one skipped row is the vendor header)");

    gh_analysis out;
    gh_analyze(&out, a, &g, &p);
    check(out.ndrug_findings == 3, "three drug findings (CFTR 1A, DPYD, CYP2C19); the level-3 CFTR one is not shown");
    int cftr = 0, dpyd = 0, ifnl4 = 0;
    for (size_t i = 0; i < out.ndrug_findings; i++) {
        const gh_drug_finding *d = &out.drug_findings[i];
        if (strcmp(d->rsid, "rs113993960") == 0) {
            cftr++;
            check_str(d->genotype, "DI", "CFTR finding carries the array genotype");
            check_str(d->annotation, "one copy of F508del: ivacaftor indicated", "CFTR annotation text");
            check_str(d->level, "1A", "CFTR level");
        } else if (strcmp(d->rsid, "rs72549303") == 0) {
            dpyd++;
            check_str(d->annotation, "normal DPYD", "DPYD II matched the GG row");
        } else if (strcmp(d->rsid, "rs11322783") == 0) {
            ifnl4++;
        }
    }
    check(cftr == 1 && dpyd == 1 && ifnl4 == 0, "F508del and DPYD reported, IFNL4 not");

    /* Allele order in the export does not matter: ID is DI. */
    const char *genome2 = write_file("genome2.txt",
        "rsid\tchromosome\tposition\tgenotype\n"
        "rs113993960\t7\t117199645\tID\n");
    gh_genome g2;
    gh_genome_load(&g2, a, genome2, NULL);
    gh_analysis out2;
    gh_analyze(&out2, a, &g2, &p);
    check(out2.ndrug_findings == 1, "ID matches the DI row");

    gh_arena_free(a);
}

int main(void)
{
    test_indel_keys_derived();
    test_indel_call_reaches_drug_guide();
    printf("%s: %d checks, %d failures\n",
           failures ? "FAILED" : "ok", checks, failures);
    return failures ? 1 : 0;
}

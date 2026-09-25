# The analysis engine

A C11 implementation of the genetic health analysis. No external libraries —
libc and libm only.

```bash
make            # build bin/genetic-health (Unix; the top-level qmake build also builds it, to ../)
make test       # unit tests (11 suites, ~6,000 checks)
make debug      # rebuild under ASan + UBSan and run the tests
```

`engine.pro` is the same engine as a qmake project, which is how the
top-level `genetic-health.pro` builds it alongside the Qt tools on every
platform (Qt's MinGW kit on Windows). The code is C11 with no POSIX
dependencies: line reading and local time go through small portable
wrappers in `gh_tsv.c` and `gh_report.c`.

```bash
./bin/genetic-health --data ../data --name "Your Name" ../data/genome.txt
./bin/genetic-health --data ../data --json ../data/genome.txt
./bin/genetic-health --data ../data --html report.html ../data/genome.txt
./bin/genetic-health --data ../data --html - --date "2024-01-01 00:00" genome.txt  # reproducible
./bin/genetic-health --list-rsids     # every rsID the analysis reads
```

The genome may be a 23andMe export, an AncestryDNA export or a WGS-converted
file; the layout is detected per row.

## Origin

This engine was ported, module by module, from a Python implementation of
the same analysis (about 17,000 lines), and the port was held to a strict
standard: for any genome file, the C produces the same findings, the same
numbers and the same HTML report as the Python, byte for byte. Everything
that was analysis moved across; what did not were the Python's process
wrappers — its CLI, its subprocess-driven WGS pipeline, its HTTP downloader
and a PDF exporter — whose roles `main.c` and the Qt tools (`qt/`) fill
natively.

The curated data tables (`src/gh_*_data.c`) were generated from the Python's
and are now maintained here directly. Genotype keys and alleles are GRCh37
plus-strand bases, exactly as export files report them.

The Python reference, its test suite, the table generators and the
differential-testing tools stay in the upstream repository. The section
below records how the port's correctness was established with them, and what
the audits they made possible found in the data; the *Design notes* record
the behaviours of the original that the C reproduces on purpose.

## How correctness was established

Two layers, because a rewrite of numeric and string-processing code fails
quietly otherwise. The first, the unit tests, is in this tree; the second,
the differential tests against the Python, ran during the port and lives
with the Python upstream.

**Unit tests** (`make test`, 6,069 checks) cover four layers. The plumbing —
arena, string buffer, hash map, TSV reader — including the awkward cases:
embedded tabs and newlines inside quoted fields, doubled quotes, text after a
closing quote, short rows, CRLF, and files with no trailing newline. Then the
report — title-casing, template substitution (unknown placeholders, CSS
braces, unterminated braces), chart geometry with zero and single categories,
section ordering, empty states, quality-metric bands, and HTML escaping of
hostile input in every interpolated field. Then the scorers — every APOE
epsilon combination and both phasings, the whole blood-type decision table,
haplogroup specificity and confidence bands, star-allele calling from
reference through homozygous variant to partial coverage, every trait branch
including the multi-SNP ones, ancestry recovery for each superpopulation, and
PRS scoring direction, confidence intervals and all three ancestry-adjustment
branches. Then the modules that consume those results rather than the genome —
drug dosing across every rule shape, polypharmacy's multi-gene matching and
severity ordering, and the preventive-care timeline: the start-age floor, the
in-place and new-entry modifier branches, both APOE tables, the BRCA and Lynch
rules, the pharmacogenomic card, and the (age, priority) sort.

**Differential tests** (`difftest.py`, upstream) generate random genomes —
including no-calls, malformed genotypes and haploid calls — run both
implementations over each, and compare *every field of every finding* plus
every quality metric, not just counts. The port matched the Python exactly
on every genome tested.

Random genomes are sampled by seed; `--exhaustive` instead enumerates the *entire* input domain of the small-domain
scorers — every rs429358 x rs7412 genotype pair, every ABO/B-antigen/Rh
combination, every genotype combination of each trait's defining SNPs, and a
sweep of genomes constructed to move estimated European ancestry across its
whole range, and the complete input space of each health profile module --
those read a handful of SNPs and branch only on copy counts, so absent/0/1/2
per marker covers every state they can distinguish (5,978 in total).

The last three modules needed a different approach again, because they read no
SNPs at all: their inputs are star-allele phenotypes, lifestyle statuses,
polygenic scores, APOE and ACMG findings, none of which a random genome moves
off its default. So those cases are constructed — every copy-count combination
of each pharmacogene's defining SNPs (which is also how the genotype producing
a given phenotype is *found*, by probing `call_star_alleles` rather than by
reading the allele definitions), each polypharmacy rule with every gene set to
every value it accepts, each polygenic condition driven to its band one at a
time, all four APOE risk levels, and pathogenic BRCA and Lynch variants sampled
out of ClinVar, since the general-purpose ClinVar sample happens to contain
none. Both were necessary: random genomes almost never produce a
decodable APOE haplotype, and they land almost exclusively below 0.5 European
ancestry, so two of the three PRS adjustment branches would never have run.

Both modes spread cases across worker processes (`--jobs`, defaulting to all
cores) and reduce the ClinVar extract once to the positions the genomes can
actually match. That matters more than it sounds: without it a single case
spends about nine seconds in Python and two in C re-reading the same 621MB
file, and the exhaustive sweep took hours. With both, it runs in about ten
seconds. A dropped row sits at a position no genome touches, so it cannot
change a result -- `--verify-reduction` demonstrates that by running cases
against both extracts and comparing findings, and `--full-clinvar` forces the
whole file.

**The rendered report was compared too**, not just the numbers behind it.
`reportdiff.py` renders the same genome through both implementations
with the same subject name and a fixed timestamp (`--date`, mirroring
`render_report(generated_date=...)`) and requires the two HTML documents to be
byte-identical. `difftest.py` runs that comparison on every case unless told
`--no-report`, so the 5,989 exhaustive cases and the random genomes cover
presentation as well as analysis. On a mismatch it names the section and
shows a unified diff of it, which is how the twelve builders were brought
into line one by one. Byte identity is a deliberately high bar: it is the
only standard under which "the C report is the Python report" is a fact
rather than an opinion.

**Curated alleles are audited against the genome build**, because parity
with the Python is worthless if both read a table that says the wrong
letter. Running a real genome (1000 Genomes NA12878, built by `gh-data
test-genome`) produced a European sample called CYP2C19 *8/*8, DPYD *13/*13
and CYP3A4 *20/*20 with 1.00 East Asian ancestry, all from alleles that had
been transcribed from the minus strand or were in fact the GRCh37 reference
base. `check_alleles.py` (upstream, `make check-alleles` there) fetches every
curated rsID's GRCh37 record from Ensembl — allele string, variant class,
1000 Genomes population frequencies, cached in `data/difftest_cache/` — and
checks the 800-odd allele uses across the star-allele, SNP-database,
ancestry, PRS, trait, profile and MT tables: a "variant" that is the
reference, an allele that only exists as its complement, an allele with ~0
frequency whose complement is common, a single-base definition on an indel,
a curated frequency that matches `1 - f`, an MT marker on a nuclear
chromosome, an MT marker equal to the rCRS base. `--fix` rewrites the Python
tables from that evidence, after which the C tables were regenerated. The first pass
changed 27 SNP-database entries, 9 star-allele definitions (and removed four
alleles whose defining SNPs do not exist as SNPs), 57 PRS risk alleles, the
ancestry, sleep, longevity and mental-health alleles, two trait allele counts
and the Rh proxy genotype. Two consequences are worth knowing:

- *The MT haplogroup tree was rebuilt.* Of the original 25 markers, three
  were on chromosomes Y, 8 and 1, nine were the rCRS base labelled as another
  lineage (so an ordinary European sample was called haplogroup C), and
  several more were near-universal or unverifiable positions. The tree is now
  twelve PhyloTree-defining positions with the rCRS coordinate in each
  description; H is recognised by the rCRS bases at m.2706 and m.7028, since
  rCRS is itself an H2a2a sequence.
- *Multi-SNP star alleles are assigned from a pool.* The old caller summed the
  variant copies over a haplotype's defining SNPs, so SLCO1B1*15 (388A>G +
  521T>C) fired homozygous on anyone carrying only the common 388G. Both
  implementations now take the minimum over the defining SNPs, assign the
  most specific haplotype first and consume the copies it uses, giving
  *15/*17 for one 521C on a 388G/G background.

The remaining flags are informational and listed in the script: sites with no
GRCh37 SNP record (indels such as ACE I/D, UGT1A1*28 and the ABO O deletion),
thirteen SNP-database entries keyed in an orientation that matches neither
strand, and PRS entries whose risk allele could not be confirmed against the
source GWAS and whose frequency was therefore taken from 1000 Genomes. The
NA12878 genome also shows the limit of that test data at CYP2D6: the 1000
Genomes call set represents the *4 and *10 positions as copy-number records,
so its CYP2D6 diplotype is not reliable.

**The reference assemblies are run as samples.** `gh-data test-genome
--sample GRCh37` (and `GRCh38`) writes a genome that is the build's own base
at every position the analysis reads, homozygous; both are checked in under
`samples/`, and upstream a test pins what they must
report. This is the one test whose right answers come from outside the
project: GRCh37's chrM *is* rCRS, so the haplogroup has to be H; the
reference carries no variant at any pharmacogene except the two whose
"variant" is the reference haplotype (CYP3A5\*3, CYP2D6\*2), so every other
gene must be \*1/\*1; and every remaining claim the report makes is a
statement about the reference genome that must be either a documented GRCh37
quirk — the reference carries the minor risk allele at Factor V Leiden
(rs6025), PTPN22 R620W, CFH Y402H, NOS3 E298D and a handful more — or a
table error. The first run found the second kind: the reference genome was
reported homozygous for GBA N370S (magnitude 4, so 99% of real users were
too), Rh-negative, and ADRB3 Trp64Arg, and eight SNP-database entries had
the variant's meaning on the common homozygote. That added three checks to
the audit: `REF_GENOTYPE_IS_RISK` (ClinVar's pathogenic allele is the
alternate but the reference homozygote carries the top magnitude),
`RARE_KEY_ALLELE` (the entry is keyed on a base with ~0 frequency at the
site — ten entries, VKORC1 −1639G>A among them, had never matched a real
genotype), and `COMMON_HOMOZYGOTE_IS_RISK`, which is informational and
resolved per rsID by the `INVERTED_ENTRIES` / `COMMON_RISK_ALLELE` lists in
the script, each with its source. The GRCh38 sample differs from GRCh37
at nine of the 601 sites, each a real change of reference base: GRCh38 has
no Factor V Leiden, is a CYP3A5 and CYP2D6 \*1 reference rather than \*3 and
\*2, and carries an A-group ABO haplotype where GRCh37 carried O. The bases
come from the assemblies themselves — `samtools faidx` over HTTPS on the
bgzipped, indexed FASTAs Ensembl publishes for each build, with GRCh37
positions carried to GRCh38 through Ensembl's `GRCh37_to_GRCh38.chain.gz`,
whose strand flag is what writes a locus GRCh38 inverted (MSMB rs10993994)
in GRCh37 orientation. The Ensembl REST API is only a fallback, for the
sites the chain does not place and for machines without samtools; the two
sources were built independently and agree base for base. The fixtures
also rode along with every differential-test run, so the C and Python were
compared on genomes with known answers as well as random ones.

**Gene names are checked against the annotation.** Every table entry names
a gene; the audit now looks up what Ensembl's GRCh37 gene set (the
`genes.gff3.gz` of the flat-file release, fetched once into the cache) has at
the rsID's position and compares. Compound names, locus names (8q24),
descriptive tags and older HGNC symbols are handled; an intergenic SNP that
names the nearest gene passes. The result was two distinct findings. In the
descriptive tables — PRS, ancestry, sleep — 44 labels named a gene on
another chromosome; several PRS models had their gene column shifted by one
row (the atrial-fibrillation model read KCNN3/NEURL1/SYNE2/SCN10A where the
SNPs sit in HCN4/WNT8A/NEURL1/SYNE2). Those are relabelled to the gene that
spans the site when one does (`--fix`), and left as informational
`NEAREST_GENE` / `GENE_NEARBY` entries when the site is intergenic, because
the nearest gene is not necessarily the locus's gene. In the SNP database,
27 entries — `WRONG_LOCUS` in the report — interpret a gene the rsID is
nowhere near: "CETP I405V" at an APOC3 promoter SNP, "MTNR1B" at a BST1
intron, "FGF21 sweet taste" at the OR2M7 asparagus-smell SNP, "LOXL1" on the
wrong chromosome. There the rsID is what is wrong, the correct one is not
derivable from the data, and the entries are listed for re-curation rather
than relabelled or removed.

The TSV reader was additionally checked byte-for-byte against Python's `csv`
module over the full 6.4 MB of real ClinPGx data. That matters: 21 annotation
texts in the allele table are quoted fields containing literal tabs, which a
naive `split('\t')` silently truncates.

## Design notes

**Arena allocation.** Everything parsed lives for the whole run, so
allocations come from an arena (`gh_mem.h`) freed in one call at exit. This
removes nearly all lifetime bookkeeping, which is where a port like this would
otherwise leak or double-free.

**Probing rather than transcribing.** The Python's `traits.py` held its
predictions and long descriptions inside branches rather than in a table, so
the generator called each predictor with a genome crafted to reach each
branch and recorded what came back. The C only has to pick the right branch
index; the wording could not drift during the port.

**Generated data tables.** The SNP database was ~2,700 lines of literals and
the report carries another ~19 KB of HTML, CSS and JS plus its palette,
pathway and reference tables. Transcribing any of that by hand would have
been error-prone, so `src/gh_snpdb_data.c`, `src/gh_report_data.c` and the
other `*_data.c` files were emitted by generators from the Python tables.
They are ordinary source now, edited directly; the generators stay upstream
with the Python.

The report template is stored as an array of chunks rather than one literal,
because C99 only guarantees 4095 characters per string literal.

**Matching Python's semantics deliberately.** Several behaviours are quirks of
the reference rather than obvious choices, and are reproduced on purpose:

- Genotype lookups retry the reversed spelling, so `CT` matches a table entry
  of `TC`.
- ClinPGx entries are keyed off the *allele* table, not the annotation table,
  so an rsID with annotations but no allele rows produces no entry — and an
  rsID's metadata comes from the first annotation to reach it in that pass.
- Only evidence levels 1A, 1B, 2A and 2B reach the report.
- The header row of a 23andMe export is not special-cased; it fails genotype
  validation and lands in the "malformed" count, exactly as in Python.
- The same goes for AncestryDNA's header row. Its two-column allele layout
  is recognised per row by two single-character allele fields; `0` is a
  no-call and chromosome codes 23/24/25/26 map to X/Y/X/MT. One classifier
  (`gh_genome_parse_row`) serves both the loader and the quality metrics'
  no-call count, so the two cannot disagree about what a no-call is.
- `quality_metrics.py` strips each line *before* testing for a leading `#`,
  unlike the genome loader, so its no-call count is recomputed from the file
  rather than reused from the loader.
- Chromosome labels have every occurrence of `CHR` removed after
  upper-casing, not just a leading prefix — Python uses `str.replace`.
- Category legends sort by descending count with ties keeping first-seen
  order, matching Python's stable sort.
- Star-allele names are compared as strings, so a diplotype sorts `*10`
  before `*2` exactly as Python's `sorted()` does.
- Star alleles are assigned from a pool of variant copies, most specific
  haplotype first, and consume the copies they use (see *Curated alleles*
  below); among equally specific alleles the one with more copies goes
  first, then definition order, matching Python's stable sort.
- A blood-type B signal overrides an otherwise-O proxy call, which is the
  Python's behaviour rather than an obvious biological rule.
- The top-ancestry pick resolves ties to the earliest population in
  `GH_POPULATIONS`, matching `max()` over an insertion-ordered dict.
- A sub-population marker that omits a sub-population falls back to a
  frequency of 0.5, as Python's `.get(sp, 0.5)` does.
- PRS models may list an rsID twice; only the first occurrence counts.
- The four dependent profiles key off the *last* finding for a gene, not the
  strongest, unlike insights and recommendations. Since findings arrive
  sorted by descending magnitude that is the weakest call, which is what
  `last_status` reproduces.
- `polypharmacy.py` builds its gene lookup from the star alleles and then
  overwrites it from the findings, so a gene named by both is judged on its
  finding -- and, again, on the last one rather than the strongest.
- `drug_dosing.py` skips a drug entirely unless one of its genes was called,
  even for a rule that would have read a lifestyle finding instead. Warfarin
  therefore reports nothing on a genome with no CYP2C9 coverage, however
  clear the VKORC1 finding is.
- A recommendation is promoted to a warning by a plain substring search for
  CONTRAINDICATED, FATAL, AVOID or LIFE-THREATENING in the upper-cased
  action, so a keyword inside a longer word counts too.
- `preventive_care.py` hard-codes the priority of a newly added screening as
  "elevated" even when the polygenic score put the condition in the "high"
  band; only screenings that modify a base guideline get "high".
- `generate_preventive_timeline` takes a `carrier_screen` argument and never
  reads it. The C signature leaves it out rather than carrying a parameter
  that does nothing.

- The report reproduces the reference's presentation exactly, including the
  parts of it that are inconsistent or wrong:
  - Some builders HTML-escape their inputs and some do not. The C escapes
    where the Python escapes, field by field.
  - Some dashes are literal UTF-8 em dashes and some are `&mdash;`; the
    apostrophe is written `&#x27;`, as `html.escape` does.
  - `svg_impact_bar` and `svg_category_donut` are defined in the Python but
    called by no section, so they render nowhere. The C keeps them for parity.
  - The eye-health and thyroid builders read `info["risk_level"]` and
    `info["detail"]`, but the profile modules write `level` and `markers`, so
    every one of those badges is "AVERAGE" in green with an empty detail.
  - `pipeline.py` strips `genes_involved` from each epistasis entry before
    saving, so the report's "Genes" line for an interaction shows its name.
  - `PAPER_REFS` is listed in insertion order, not sorted; the generator now
    preserves that order.

One quirk that used to be on this list has been removed from the reference
instead. `snp_database.py` held three entries keyed `rs4680_pain`,
`rs1800795_longevity` and `rs28929474_serpina` -- suffixed so that a second
curation of the same variant could sit alongside the first in a plain dict.
The key is used directly as the genome lookup key, and no genome carries an
rsID shaped like that, so none of the three ever produced a finding in either
implementation. They were identically dead in C and Python, which is exactly
why the differential test could not see them. All three were redundant
(`pain_sensitivity.py` reads rs4680 itself, the Inflammation entry for
rs1800795 already drives longevity's inflammation domain, and the SERPINA1
pair were near-identical) so they were deleted rather than made reachable.
The only output change is the SERPINA1 note, which moved to the live entry.

**Streaming the large input.** ClinVar's extract is 621 MB and 4.5 million
rows. Slurping it the way the other tables are read cost 658 MB resident, so
`gh_tsv_open_streaming` reads one record at a time, joining continuation
lines when a quoted field spans them. That brought peak memory to 22 MB with
byte-identical output, and the whole run still takes under two seconds
against nine for the Python.

**A missing sort the port hid from itself.** `analyze_lifestyle_health`
sorts its two result lists before returning -- findings by descending
magnitude, drug findings by evidence level -- and the C did not. It went
unnoticed for a long time because the differential test compared findings
keyed by rsID, which throws order away. It only surfaced once
`recommendations.py` was ported, since that consumes the list order. The
sorts are now in `gh_analyze`, the ClinPGx loader keeps insertion order so
stable-sort ties break the same way, and the differential test compares both
lists in order.

**Two bugs the port surfaced in carrier screening, since fixed in the
Python.** `update_data.py` fills the `inheritance_modes` column from
ClinVar's `OriginSimple` field, which holds germline/somatic origin rather
than a mode of inheritance — across all 4.5 million rows, not one contains
"recessive" or "dominant". Carrier screening decided status from that field,
so it reported zero carriers for every genome. Separately, the curated gene
table spells its X-linked entries with a capital X while the note lookup
tested for lowercase "x-linked", so that reproductive note was unreachable.

`carrier_screen.py` now resolves inheritance from the curated table before
deciding, and matches case-insensitively; `gh_clinvar.c` mirrors both. The
practical limit is that carrier screening only covers the 27 genes the
curated table gives a recessive or X-linked pattern for — everything else
still falls back to a column that cannot answer the question.

**Two determinism bugs in `recommendations.py`, since fixed.** Risk groups
store their genes as a `set`, and Python randomises string hashing per
process, so iterating them directly leaked hash order into the user-facing
"why" text and the clinical action list: eight runs of the same genome gave
eight different reports. Both iterations are now sorted. A third site picked
an arbitrary element from a set intersection; only one trigger uses it and it
lists a single gene, so it was latent rather than live, but it is sorted too.

**A duplicated sentence.** The guard meant to emit the ClinVar risk-factor
line once counted occurrences in the wrong list, so the same sentence
repeated once per matching variant. Fixed with a regression test.

**Floating point and rounding.** The PRS module reports values Python has
already passed through `round()`, which breaks ties to even. Computing
`floor(x * 10^n + 0.5)` would disagree on exact halves, so `gh_round` formats
to a decimal string and parses it back — the same correctly-rounded path
Python takes. Ancestry proportions come out of a softmax and match to within
1e-9 across every case tested.

**A bug the port surfaced.** The exhaustive APOE sweep crashed the *Python*
with an `IndexError`: `_try_all_phasings` indexed `genotype[1]` on both SNPs,
but `load_genome` accepts single-base (haploid) genotypes. Its own
`_decode_haplotype` already guarded `len != 2`, so the intended behaviour was
clearly "Unknown" — the guard simply ran too late. Fixed in `apoe.py` with
regression tests upstream; the C had the guard from the start.

**The annotation ID column has two names.** ClinPGx renamed
`Clinical Annotation ID` to `Summary Annotation ID`; both are accepted, in
`GH_PGX_ID_COLUMNS`.

## Layout

```
c/
├── include/         public headers, one per module
├── src/             implementations; *_data.c files are the curated tables
├── tests/           C unit tests
└── Makefile
```

Built with `-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion` and
warning-free.

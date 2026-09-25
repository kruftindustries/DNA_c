# Qt front end

A desktop shell for the analysis, plus the two pieces the C engine does not
do itself: the annotation downloader and the FASTQ pipeline, implemented
natively in `qt/core`. The analysis itself is the C binary, driven through
`QProcess`; the GUI adds no analysis logic of its own.

```bash
make build                # top-level qmake: ./genetic-health-qt, ./gh-data and the engine ./genetic-health
make qt-test              # Qt unit tests (offscreen)
make validate-wgs         # native chr22 toolchain validation
./genetic-health-qt
```

The build is one qmake project (`../genetic-health.pro`) that also compiles
the C engine (`c/engine.pro`), so a plain Qt installation — the MinGW kit on
Windows — builds everything; the three executables land in the repository
root and the app looks for the engine beside itself first.

## gh-data

The core without the GUI, as a console tool:

```bash
gh-data clinvar                 # download variant_summary.txt.gz, reduce to GRCh37
gh-data pharmgkb                # download ClinPGx annotations
gh-data status                  # installed releases
gh-data reference               # fetch + index the GRCh37 reference for the WGS pipeline (~12 GB)
gh-data wgs reads.fastq.gz --name "Subject"   # FASTQ -> genotypes -> report
gh-data wgs-validate            # simulate reads on chr22, run, compare to truth
gh-data tools                   # where minimap2/samtools/bcftools/fastp were found
gh-data test-genome             # a public genome to test with (below)
```

### A genome to test with

`gh-data test-genome` (or the button on the Input page) builds a real,
public genome: 1000 Genomes Project sample NA12878 -- the most-studied human
genome, openly consented -- extracted at the ~620 positions this analysis
reads. The phase 3 VCFs are indexed, so `bcftools` fetches just those rows
over HTTPS; positions with no variant in the 2,504-sample cohort are written
as homozygous reference, with the base from Ensembl, so the file has a call
at every position the way an array export would. It lands in
`data/genome_1000g_NA12878.txt` and is gitignored like the rest of `data/`.
`--sample` picks any other 1000 Genomes sample. The build can be stopped
(the Stop button, or Ctrl-C on `gh-data`) and resumed: every chromosome's
calls and every fetched reference base are kept under
`data/test_genome_cache/` as they arrive, and the next run reads them
instead of asking again. A chromosome's cache records the positions it was
queried for, so it is refetched if the rsID list changes.

`--sample GRCh37` or `--sample GRCh38` writes the reference assembly itself as
a sample instead: its own base at every position, homozygous. That is a null
test of the curated tables -- GRCh37's mtDNA is rCRS, so the haplogroup must
be H, and every other variant the report claims is either a documented
reference quirk (GRCh37 carries Factor V Leiden) or a table error. The bases
are read with `samtools faidx` straight from the bgzipped, indexed FASTA
Ensembl publishes for each build (a few hundred HTTPS range requests, no API);
GRCh37 positions are carried to GRCh38 through Ensembl's chain file, cached
in `data/ensembl/`, whose strand flag keeps inverted loci in GRCh37
orientation. Without samtools the builder falls back to the Ensembl REST API.
The two files are checked in under `samples/`; `make reference-samples`
regenerates them.

### Theme

The app follows the desktop's light/dark preference. Qt 5 takes fonts and
dialogs from the platform theme but, on most Linux desktops, not the colours,
so a dark GTK theme would otherwise leave the window light. `SystemTheme`
asks the desktop -- the XDG settings portal's `color-scheme`, then gsettings,
then the GTK theme name; the registry on Windows -- and when it prefers dark
switches to the Fusion style with a dark palette, then keeps listening so the
app flips with the desktop. Qt 6.5+ propagates the scheme itself and is left
to it, as is anyone who configures Qt theming explicitly
(`QT_QPA_PLATFORMTHEME`, `QT_STYLE_OVERRIDE`).

### Requirements

The FASTQ pipeline and the test-genome builder shell out to minimap2,
samtools, bcftools and fastp. They are looked up on `PATH`; when a required
one is missing, the WGS page, the Input page and `gh-data tools` say so and
print the install command for the platform:

```bash
sudo apt install minimap2 samtools bcftools fastp
```

Every step was written against the original Python pipeline and verified to
match it: the ClinVar filter writes byte-identical output, the target regions
come from the same positions, the VCF conversion applies the same rules, and
the pipeline runs the same commands. Two things are done natively: the rsID
list the target regions are built from comes from the C binary
(`genetic-health --list-rsids`) rather than a hand-copied list, which had
missed the sleep, longevity, mental-health and sub-population markers; and
gzip and zip are decoded in-process by the vendored
[miniz](core/third_party/miniz) (MIT), so ClinVar's `variant_summary.txt.gz`,
the ClinPGx zip, Ensembl's chain file and the reference genome need no
`gzip` or `unzip` binary on any platform. The reader accepts multi-member
files (bgzip) and trailing non-gzip data: the 1000 Genomes reference is a
samtools RAZF file, a gzip stream with a block index appended, which `gzip`
reports as "trailing garbage ignored" (`tst_compression` covers all three
layouts and zip members).

Needs Qt 5.15 or Qt 6 with `qmake` (Widgets, Network; Test for the tests).
Qt WebEngine is optional: with it the interactive report renders in-app,
without it the Report page shows Qt's basic HTML rendering and opens the
full report in the system browser.

## Pages

| Page | What it does |
|---|---|
| **Input** | Pick a genome file. Its layout is detected from the contents, not the extension: 23andMe, AncestryDNA, VCF or FASTQ. A FASTQ is routed to the WGS page. |
| **Run** | Live log, progress and cancel for the analysis (`--json`, then `--html`). |
| **Report** | The generated report, with *Open in browser* and *Save copy*. |
| **Findings** | The findings table: search, category and magnitude filters, TSV export of the visible rows. |
| **Data Sources** | ClinVar and ClinPGx: release, last update, files present. Update and validate natively, with progress, on a worker thread. |
| **WGS Pipeline** | Tool discovery, reference status with *Set up GRCh37 reference…* (download, decompress, samtools and minimap2 indexes; resumable), the chr22 toolchain validation, and the full FASTQ run -- all native, on a worker thread. Choosing a FASTQ on the Input page and pressing Run comes here, offers to fetch the reference if it is missing, and runs. A finished run opens its report. |

*Settings* (File menu) overrides the data directory and the analysis binary;
by default they are found relative to the checkout that contains the
executable.

## What is deliberately not here yet

PDF export (Qt WebEngine's `printToPdf` once WebEngine is installed). It
does not change what the analysis computes.

## Layout

```
qt/
  qt.pro           subdirs: app, cli, tests (the top-level ../genetic-health.pro adds the engine)
  common.pri       sources shared by the app and the tests
  core/            no GUI: TextTable (csv-compatible TSV), GzipStream and
                   ZipArchive (miniz), Downloader, ClinVarUpdater,
                   PharmgkbUpdater, EnsemblLookup, TargetRegions, VcfConverter,
                   WgsPipeline, WgsValidator, ReferenceGenome, TestGenome,
                   ChainLift, RemoteFasta, ToolLocator
  cli/             gh-data
  app/             GenomeFormat, RepoPaths, ProcessRunner, AnalysisRunner,
                   FindingsModel, JobRunner, SystemTheme, the six pages, MainWindow
  tests/           QtTest suites: format sniffer, findings model, table I/O,
                   VCF conversion, target regions, Ensembl merge, ClinVar
                   processing, test-genome calls, chain lift, compression, and
                   an end-to-end run of the analysis binary (test.pri gives
                   each its own object directory so make -j is safe)
```

# Genetic Health Analysis (C engine + Qt desktop app)

Turn raw DNA data into a health report you can read and a doctor can use.
Takes **23andMe and AncestryDNA raw-data exports** and **whole-genome
sequencing FASTQ files**, cross-references the genome against ClinVar
(~340K clinical variants), ClinPGx/PharmGKB (drug–gene interactions) and a
curated database of ~260 lifestyle and health SNPs, and writes a single
self-contained HTML report: key findings in plain language, drug guide,
disease risk, pharmacogenomic star alleles, polygenic risk scores, ancestry,
traits, and a print-optimised doctor card.

Everything runs on your machine. The only network use is fetching the public
annotation data, the reference genome for the sequencing path, and,
optionally, public test genomes — all from inside the app.

## Credit

This is a C and Qt implementation of the Python pipeline in
**[unbalancedparentheses/DNA](https://github.com/unbalancedparentheses/DNA)**,
which designed the analysis, curated the SNP tables and wrote the report. The
port was verified against that code down to byte-identical HTML output, and
the curated tables were audited against Ensembl GRCh37 and corrected along
the way (the original repository carries the Python reference, the test
suite and the audit tooling). The original is not for clinical use; neither
is this.



## What is in the box

| Executable | What it is |
|---|---|
| `genetic-health-qt` | The desktop app: choose a file, run, browse findings, open the report, manage the data. |
| `gh-data` | The same machinery on the command line: data downloads, reference setup, the WGS pipeline, test genomes. |
| `genetic-health` | The analysis engine (C11, libc and libm only). The app drives it; it also runs on its own. |

All three are built into the repository root by one build. The trees:
`c/` (engine, unit tests), `qt/` (app, `gh-data`, QtTest suites),
`samples/` (the GRCh37 and GRCh38 reference assemblies as genomes, for null
tests), `data/`, `reference/`, `reports/` (all gitignored; filled in by the
app).

## Downloads

Release builds for all three platforms are attached to each
[GitHub release](https://github.com/kruftindustries/DNA_c/releases), built by
the workflow in `.github/workflows/release.yml` from a version tag:

| Platform | Package | Run |
|---|---|---|
| Linux x86_64 | `Genetic_Health-x86_64.AppImage` | `chmod +x` and run (needs FUSE 2: `sudo apt install libfuse2`, or `--appimage-extract` and run `squashfs-root/usr/bin/genetic-health-qt`). `genetic-health-linux-x86_64.tar.gz` is the same tree as a folder, with `usr/bin/gh-data` and the engine beside the app. |
| Windows x64 | `genetic-health-windows-x64.zip` | Unzip anywhere; run `genetic-health-qt.exe`. `gh-data.exe` and `genetic-health.exe` are in the same folder. The build is unsigned, so SmartScreen asks once (*More info → Run anyway*). |
| macOS x86_64 | `genetic-health-macos-x86_64.dmg` | Drag *Genetic Health* to Applications. Unsigned: first launch is right-click → *Open*. `gh-data` and the engine are inside the bundle at `Genetic Health.app/Contents/MacOS/`. Runs under Rosetta on Apple silicon. |

A packaged build keeps its data — annotation downloads, the reference genome,
reports — under the per-user application data directory
(`~/.local/share/genetic-health`, `%LOCALAPPDATA%\genetic-health`,
`~/Library/Application Support/genetic-health`); a build inside a source
checkout uses the checkout's `data/`, `reference/` and `reports/`. Settings
in the app can point either elsewhere.

The packages are produced by the scripts in `packaging/`
(`linux-appimage.sh`, `windows-package.bat`, `macos-dmg.sh`), which can be
run locally after a build; the workflow runs exactly those.

## Building

### Windows

Prerequisites: [Qt](https://www.qt.io/download-qt-installer) **6.8** with
the **MinGW 64-bit** kit (the online installer offers it under the Qt
version; it brings the compiler and `mingw32-make`, so nothing else is
needed), and [Git for Windows](https://git-scm.com/download/win) or a zip of
this repository. Qt 6 matters on Windows: its built-in Schannel TLS backend
gives the app HTTPS for the data downloads, whereas Qt 5.15 needs OpenSSL 1.1
DLLs that Qt no longer distributes.

Either open `genetic-health.pro` in **Qt Creator**, pick the MinGW kit when
asked, and press *Build* — or from the **Qt MinGW command prompt** (Start
menu → the Qt folder → *Qt 6.8.x (MinGW 64-bit)*):

```bat
cd \path\to\DNA_c
mkdir build
cd build
qmake ..\genetic-health.pro
mingw32-make -j
```

`genetic-health-qt.exe`, `gh-data.exe` and `genetic-health.exe` appear in the
repository root. To run the app from Explorer or a shortcut, the Qt DLLs must
be findable: run `windeployqt genetic-health-qt.exe` once from the same Qt
prompt (it copies them next to the executable), or start the app from that
prompt. The MSVC kit is untested; the engine is standard C11 and the Qt code
is standard Qt, so it should work, but the MinGW kit is the supported path.

Everything in the 23andMe/AncestryDNA workflow — the app, the data downloads,
the report — is native and needs nothing beyond Qt. The FASTQ workflow needs
the sequencing tools (below), which do not have Windows builds; use WSL2
(Ubuntu) and follow the Linux steps there.

### Linux

Tested on Debian/Ubuntu with Qt 5.15:

```bash
sudo apt install build-essential qtbase5-dev qt5-qmake   # compiler + Qt
git clone https://github.com/kruftindustries/DNA_c && cd DNA_c
make build          # = mkdir build && cd build && qmake ../genetic-health.pro && make
./genetic-health-qt
```

`make test` runs the engine's unit tests (via `c/Makefile`) and the Qt
suites. Other distributions: the packages are the Qt 5 (or 6) base
development files and qmake; on Fedora `qt5-qtbase-devel`, on Arch
`qt5-base`.

### macOS

TODO. Expected to be `brew install qt@5` and the same `qmake` + `make`
sequence as Linux (`app_bundle` is left on, so the app builds as a `.app`);
untested.

## How to use

### 23andMe or AncestryDNA export

Prerequisites: the built app, an internet connection for the first-run data
download, and your raw-data file — the `.txt` (or `.zip`, unzipped) that
23andMe's *Browse Raw Data → Download* or AncestryDNA's *Download DNA data*
gives you. No other tools.

1. Start `genetic-health-qt`. On **Data Sources**, press *Update ClinVar* and
   *Update ClinPGx* (a few minutes; ClinVar is a ~430 MB download reduced to
   the GRCh37 variants on the fly). Without them the report still runs, minus
   the disease-risk and drug-interaction sections.
2. On **Input**, browse to your file. The format is detected from the
   contents (23andMe's four columns, AncestryDNA's five); a subject name is
   optional.
3. Press **Run analysis**. The report opens on the **Report** page and the
   **Findings** page lists every finding with filters; *Open in browser*
   gives the full interactive report, and *Save copy…* keeps it.

No export to hand? *Get a public test genome…* builds 1000 Genomes sample
NA12878 at the analysed positions (needs `bcftools`), or writes the GRCh37 or
GRCh38 reference assembly itself as a sample (needs `samtools`) — the latter
are also checked in under `samples/`.

Command line equivalent:

```bash
./gh-data clinvar && ./gh-data pharmgkb
./genetic-health ~/Downloads/genome_Your_Name.txt --name "Your Name" --html reports/report.html
```

### Sequencing reads (FASTQ)

Prerequisites, on Linux or WSL2: **minimap2**, **samtools**, **bcftools**,
and optionally **fastp** for read QC (`sudo apt install minimap2 samtools
bcftools fastp`; the app and `gh-data tools` print this command when
something is missing), about **17 GB of disk** (the GRCh37 reference and its
indexes take ~12 GB, the alignment a few more), **16 GB of RAM** (minimap2
holds the whole-genome index in memory, ~10.5 GB), and time: alignment of a
30x genome takes hours.

1. On **Input**, browse to the FASTQ (`.fastq`, `.fq`, or gzipped). The
   button becomes **Run WGS pipeline + analysis**.
2. Press it. If the GRCh37 reference is not set up yet the app offers to
   fetch and index it (~900 MB download, 10–20 minutes, resumable) and then
   continues; otherwise the pipeline starts at once. The **WGS Pipeline**
   page shows progress and the log, and lets you cancel; the same page has
   *Set up GRCh37 reference…* and *Validate toolchain (chr22)*, a
   twenty-second end-to-end check of the tools against simulated reads with
   known genotypes.
3. The pipeline runs fastp → minimap2 → bcftools (targeted calling at the
   positions the analysis reads; *Genome-wide calling* is an option) →
   conversion to the array layout → the analysis. When it finishes the report
   opens as above.

Command line equivalent:

```bash
./gh-data reference                 # once: fetch + index GRCh37
./gh-data wgs reads.fastq.gz --name "Your Name" --threads 12
```

### Managing the data

The **Data Sources** page lists everything the app has fetched, with sizes
and dates, and manages it: ClinVar and ClinPGx (update, validate), the rsID →
GRCh37 position lookup (update), the reference genome (set up), and the
caches behind the test genomes (clear). Refresh ClinVar and ClinPGx roughly
quarterly; the files live in `data/`, the reference in `reference/`.

#### Reference data and toolchains are not included with this software, these works are licensed by their respective owners

## The report

One self-contained HTML file in twelve sections. **Part 1 — what you should
do:** Key Findings (colour-coded, plain-language), Action Plan, Drug Guide,
Disease Risk, Body Profile, Mental Health. **Part 2 — details for your
doctor:** Clinical Findings (ClinVar, ACMG secondary findings, carrier
screening), Ancestry, Nutrigenomics, Data Quality, Doctor Card
(print-optimised), References. SVG charts, search and filter, sortable
tables, CSV export, dark mode, and a database link for every rsID. For a
PDF, print the page from a browser.

## What is analysed

Curated SNPs (~260 variants, 16 categories, per-genotype interpretation and
impact magnitude 0–6) · ClinVar pathogenic SNVs with zygosity and gold
stars, ACMG SF v3.2 secondary findings (81 genes), carrier screening ·
CPIC-style star alleles for eleven pharmacogenes with metabolizer
phenotypes, dosing for nine drugs, polypharmacy warnings · polygenic risk
scores for 25 conditions · ancestry from ~55 AIMs with sub-population
estimates · blood type, APOE, mitochondrial haplogroup (twelve PhyloTree
markers), 14 traits, ten gene–gene interactions, body profiles (pain,
histamine, thyroid, hormones, eye health, alcohol, sleep, nutrigenomics,
mental health, longevity), a preventive-care timeline, and data quality.
`DEVELOPING.md` has the module-by-module reference; `c/README.md` the design
notes and what the audits of the tables found.

## Testing

`make test` runs the engine's ~6,000 unit checks (11 suites) and the 11
QtTest suites. `samples/genome_reference_GRCh37.txt` is a null test with
answers known from outside the project: run through the engine it must
report haplogroup H (rCRS), every pharmacogene \*1/\*1 except CYP3A5\*3 and
CYP2D6\*2 (the GRCh37 reference haplotypes), and Factor V Leiden (which
GRCh37 carries and GRCh38 corrected). `make validate-wgs` checks the
sequencing toolchain on chr22.

## Interpretation and limitations

Magnitude 0 informational, 1 low, 2 moderate, 3 actionable, 4–6 clinical
attention. ClinVar gold stars 4 (expert panel) down to 0 (no criteria).
Heterozygous findings in recessive genes are carrier status.

Not a clinical diagnosis — see a genetic counsellor or physician for
decisions. Most evidence is European-biased and the report flags where that
matters. A pathogenic variant is not a diagnosis. Indels are not analysed.
Refresh the annotation data quarterly. Low-coverage sequencing misses
positions; the Data Quality section shows how many were read.

## License

Developed for personal use. Not for clinical or diagnostic purposes. It may produce erroneous results.
### This software is not intended to diagnose, treat, cure or prevent any medical condition

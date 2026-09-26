# Developer reference

A C analysis engine (`c/`) and a Qt desktop app plus the `gh-data` console
tool (`qt/`) that turn a genome file — a 23andMe or AncestryDNA export, or a
WGS-converted file — into a twelve-section HTML health report. This tree was
ported from a Python implementation and verified against it down to
byte-identical HTML; that reference, its test suite, the table generators and
the differential tests live in the upstream repository. Nothing here needs
Python or Nix.

## Quick start

```bash
make build                       # one qmake build -> ./genetic-health ./genetic-health-qt ./gh-data
make update-data                 # ClinVar + ClinPGx into data/ (gh-data)
make analysis NAME="John Doe"    # data/genome.txt -> reports/GENETIC_HEALTH_REPORT.html
./genetic-health /path/to/genome.txt --name "Jane Doe" --html out.html
make gui                         # the desktop app
make test                        # C unit tests + Qt unit tests
```

The build is the top-level `genetic-health.pro` (subdirs: `c/engine.pro`,
`qt/app`, `qt/cli`, `qt/tests`); `c/Makefile` is the engine's own Unix
build with its unit tests. On Windows the same project builds with Qt's
MinGW kit (`qmake` + `mingw32-make`); see README.md.

## Layout

```
DNA_c/
├── genetic-health.pro  top-level qmake project (engine + app + gh-data + tests)
├── c/
│   ├── engine.pro      the engine as a qmake project (used by the top-level build)
│   ├── include/        one public header per module
│   ├── src/            the engine; gh_*_data.c are the curated tables
│   ├── tests/          C unit tests (make -C c test)
│   ├── bin/            c/Makefile's build output (gitignored)
│   └── README.md       design notes, reproduced behaviours, what the audits found
├── qt/
│   ├── core/           QtCore+Network: downloader, gzip/zip (miniz), WGS pipeline, reference setup, test genomes
│   ├── app/            the desktop app (pages, format sniffer, process runners, theme)
│   ├── cli/            gh-data
│   ├── tests/          QtTest suites (make qt-test)
│   └── README.md
├── samples/            GRCh37 and GRCh38 reference assemblies as genomes
├── data/               annotation data, genomes, caches (gitignored except rsid_positions_grch37.json)
├── reference/          GRCh37 reference + indexes for the WGS pipeline (gitignored)
└── reports/            generated reports (gitignored)
```

## The engine (`c/`)

C11, libc and libm only, warning-free under `-Wall -Wextra -Wpedantic
-Wconversion -Wsign-conversion`. `main.c` orchestrates: load the genome
(layout detected per row: 23andMe's four columns, AncestryDNA's five with
`0` no-calls and 23/24/25/26 chromosome codes), load ClinVar and ClinPGx,
then run every module and render the report.

| Module | Files | What it does |
|--------|-------|--------------|
| Genome + data loading | `gh_genome`, `gh_tsv`, `gh_pgx` | vendor layouts, streaming TSV (ClinVar is 4.5M rows read in 22 MB), csv-compatible quoting |
| Curated SNP analysis | `gh_analyze`, `gh_snpdb_data` | ~260 variants, 16 categories, per-genotype status/magnitude |
| ClinVar | `gh_clinvar`, `gh_clinvar_data` | pathogenic SNVs with zygosity and gold stars; ACMG SF v3.2 (81 genes); carrier screening by system |
| Scorers | `gh_scorers`, `gh_scorer_data` | APOE, blood type (rs505922/rs8176746/rs590787), MT haplogroup (12 PhyloTree markers), star alleles (11 genes, pool-based haplotype assignment) |
| Traits, ancestry, PRS | `gh_traits`, `gh_ancestry`, `gh_prs`, `gh_ancestry_prs_data` | 14 traits; ~55 AIMs by maximum likelihood; 25 GWAS models |
| Profiles | `gh_profiles`, `gh_dependent_profiles`, `gh_profile_data` | pain, histamine, thyroid, hormones, eye, alcohol, sleep, nutrigenomics, mental health, longevity |
| Consumers | `gh_dosing`, `gh_epistasis`, `gh_recommendations`, `gh_insights` | CPIC/DPWG dosing (9 drugs), polypharmacy, preventive care, 10 interaction models, recommendations, narratives |
| Report | `gh_report`, `gh_report_sections`, `gh_report_data` | the twelve sections, SVG charts, escaping |

Flags: `--html PATH` (`-` for stdout), `--json`, `--name`, `--data DIR`,
`--date TEXT` (fixed timestamp for reproducible output), `--list-rsids`
(every rsID the analysis reads — the single source of truth the WGS target
regions and the test genomes are built from).

**Editing the curated tables.** They are the `gh_*_data.c` files. Genotype
keys and alleles are GRCh37 plus-strand bases, exactly as export files report
them. After a change, rebuild and run `samples/genome_reference_GRCh37.txt`:
the reference genome must not be reported as carrying a variant unless the
GRCh37 base really is the variant allele (CYP3A5\*3, Factor V Leiden and a
handful of others are; `c/README.md` lists them). That single check catches
strand flips and reference-as-variant errors, which were the two most common
mistakes the audit found in the original tables.

## The Qt tools (`qt/`)

Qt 5.15+ or Qt 6 with qmake; Widgets, Network (D-Bus on Linux for the
desktop theme), Test for the tests; WebEngine optional for in-app report
rendering.

- **Desktop app** (`./genetic-health-qt`): Input (file chooser with format
  detection; Run on an array export runs the engine, Run on a FASTQ runs the
  WGS pipeline first, offering to fetch the reference; public test genomes
  with stop/resume), Run, Report, Findings (filterable table), Data Sources
  (ClinVar, ClinPGx, rsID positions, the GRCh37 reference, caches — status,
  update, set up, clear), WGS Pipeline (tools, reference setup, chr22
  validation, the run). Follows the desktop's light/dark preference.
- **`gh-data`**: `clinvar`, `pharmgkb`, `validate-pharmgkb`, `status`,
  `reference [--status]`, `wgs FASTQ`, `wgs-validate`, `build-targets`,
  `convert-vcf`, `lookup-rsids`, `test-genome [--sample NA12878|GRCh37|GRCh38]`,
  `tools`.
- **Compressed inputs** are decoded in-process by the vendored miniz
  (`qt/core/third_party/miniz`, MIT): no `gzip`/`unzip` binaries.
- **Test genomes**: NA12878 from the 1000 Genomes VCFs by `bcftools` range
  requests (no full download); the reference assemblies from Ensembl's
  indexed FASTAs by `samtools faidx`, with GRCh37 positions carried to GRCh38
  through Ensembl's chain file (`data/ensembl/`). Everything fetched is cached
  under `data/test_genome_cache/`, so a stopped build resumes.
- **WGS**: fastp → minimap2 → bcftools targeted calling at the analysed
  positions → conversion to the array layout → the engine. The GRCh37
  reference (`reference/`, ~12 GB) is fetched, decompressed and indexed by
  `ReferenceGenome::setup` (`gh-data reference`, or the buttons on the WGS
  and Data Sources pages); each step resumes if interrupted. `make
  validate-wgs` checks the toolchain on chr22 against simulated reads with
  known genotypes. Tools are found on `PATH`; `make tools` prints the install
  command. Linux or WSL2 only — the tools have no Windows builds.

## Data files (`data/`)

| File | Source | Used for |
|------|--------|----------|
| `clinvar_alleles.tsv` | `gh-data clinvar` (ClinVar FTP, GRCh37 SNVs) | disease risk, ACMG, carrier screening |
| `clinical_annotations.tsv`, `clinical_ann_alleles.tsv` | `gh-data pharmgkb` (ClinPGx) | drug–gene interactions |
| `rsid_positions_grch37.json` | `gh-data lookup-rsids` (Ensembl); tracked, and compiled into the executables (`qt/core/core.qrc`) so a packaged build seeds its data directory with it | WGS conversion, test genomes |
| `data_versions.json` | written by `gh-data` | release bookkeeping |

Missing annotation files skip their sections; the report is still produced.

External tools (minimap2, samtools, bcftools; fastp optional) are found by
`ToolLocator`: the `tools/` folder beside the executable first, then PATH,
then the Homebrew/MacPorts/`~/.local` prefixes a desktop-launched app does
not have on PATH. The Windows release ships them in `tools\`, built by
`packaging/windows-tools.sh` in MSYS2 MinGW64 (the `windows-tools` CI job;
fastp needs `packaging/fastp-mingw.{h,cpp}` for pwrite and 64-bit ftell);
Linux and macOS builds tell the user the package-manager command.

## Interpretation

Magnitude 0 informational, 1 low, 2 moderate, 3 actionable, 4–6 clinical
attention. ClinVar gold stars 4 (expert panel) down to 0 (no criteria).
Heterozygous findings in recessive genes are carrier status. Not a clinical
diagnosis; most evidence is European-biased; indels are not analysed.

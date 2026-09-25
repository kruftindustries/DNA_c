SHELL := /bin/bash

# The C analysis engine (c/) and the Qt desktop app plus the gh-data console
# tool (qt/). This tree was ported from a Python implementation, which is
# kept in the upstream repository as the reference the port was verified
# against; nothing here needs Python.

# FASTQ input for `make pipeline`: make pipeline FASTQ=/path/to/reads.fastq
FASTQ ?=
NAME ?=
THREADS ?= 12

# Directories
REF_DIR := reference
WORK_DIR := wgs_work
DATA_DIR := data
REPORTS_DIR := reports

# Reference files (WGS path only)
REF_FASTA := $(REF_DIR)/human_g1k_v37.fasta
REF_MMI := $(REF_DIR)/human_g1k_v37.mmi
REF_FAI := $(REF_FASTA).fai

# Inputs and outputs
GENOME := $(DATA_DIR)/genome.txt
REPORT := $(REPORTS_DIR)/GENETIC_HEALTH_REPORT.html
RSID_LOOKUP := $(DATA_DIR)/rsid_positions_grch37.json

# Binaries: the top-level qmake build puts all three here.
BUILD_DIR := build
ANALYSIS := genetic-health
GH_DATA := gh-data
GUI := genetic-health-qt

# Tools needed only by the WGS (FASTQ) pipeline
WGS_TOOLS := minimap2 samtools bcftools

ifdef NAME
  NAME_FLAG := --name "$(NAME)"
else
  NAME_FLAG :=
endif

.PHONY: help build c qt gui analysis pipeline test test-c qt-test \
        reference-samples test-genome update-data update-clinvar \
        update-pharmgkb validate-pharmgkb data-status validate-wgs validate-wgs-setup \
        setup tools check-tools clean clean-all

## Default target: show help
help:
	@echo "Genetic Health Analysis"
	@echo "======================="
	@echo ""
	@echo "Build:"
	@echo "  make build          Build the C engine (c/bin/genetic-health) and the Qt tools"
	@echo "  make gui            Build and start the desktop app"
	@echo ""
	@echo "Run:"
	@echo "  make analysis       Analyse data/genome.txt (23andMe or AncestryDNA export)"
	@echo "  make pipeline       Full WGS pipeline: FASTQ -> report (gh-data wgs)"
	@echo "  make test-genome    Build a public test genome (1000 Genomes NA12878)"
	@echo ""
	@echo "Data:"
	@echo "  make update-data    Download ClinVar + ClinPGx annotation data"
	@echo "  make data-status    Show installed data releases"
	@echo "  make setup          Download the GRCh37 reference + build indexes (WGS only)"
	@echo "  make tools          Print the install command for the WGS tools"
	@echo ""
	@echo "Test:"
	@echo "  make test           C unit tests and Qt unit tests"
	@echo "  make validate-wgs   chr22 truth-set validation of the WGS toolchain"
	@echo "  make reference-samples  Rebuild the GRCh37/GRCh38 reference-as-sample genomes"
	@echo ""
	@echo "Options:"
	@echo "  FASTQ=path          Input FASTQ file for 'make pipeline'"
	@echo "  NAME=\"John Doe\"     Subject name for the report"
	@echo "  THREADS=12          CPU threads for the WGS pipeline"
	@echo ""
	@echo "Cleanup:"
	@echo "  make clean          Remove WGS intermediate files"
	@echo "  make clean-all      Remove all generated files"

## Build everything through the top-level qmake project (genetic-health.pro):
## the engine, the desktop app, gh-data and the Qt tests. Needs a C/C++
## compiler and Qt 5.15+ or Qt 6 with qmake.
build:
	mkdir -p $(BUILD_DIR) && cd $(BUILD_DIR) && qmake ../genetic-health.pro && $(MAKE)

c: build
qt: build

gui: build
	./$(GUI)

## Analyse an array-style export with the C engine
analysis: build
	@test -f "$(GENOME)" || { echo "genome file not found: $(GENOME). Copy your export there, or run 'make test-genome' for a public one."; exit 1; }
	@mkdir -p $(REPORTS_DIR)
	./$(ANALYSIS) "$(GENOME)" $(NAME_FLAG) --html $(REPORT)
	@echo "Report: $(REPORT)"

## Full pipeline: FASTQ -> QC -> align -> call -> convert -> analyse
pipeline: check-tools build
	@test -n "$(FASTQ)" -a -f "$(FASTQ)" || { echo "Usage: make pipeline FASTQ=/path/to/reads.fastq [NAME=\"Subject\"]"; exit 1; }
	./$(GH_DATA) wgs "$(FASTQ)" --threads $(THREADS) $(NAME_FLAG)

## A public genome to test with (needs bcftools)
test-genome: build
	./$(GH_DATA) test-genome

## Annotation data
update-data: build
	./$(GH_DATA) clinvar
	./$(GH_DATA) pharmgkb

update-clinvar: build
	./$(GH_DATA) clinvar

update-pharmgkb: build
	./$(GH_DATA) pharmgkb

validate-pharmgkb: build
	./$(GH_DATA) validate-pharmgkb

data-status: build
	./$(GH_DATA) status

## Tests
test: test-c qt-test

test-c:
	$(MAKE) -C c test   # the engine's own build + unit tests

qt-test: build
	cd $(BUILD_DIR)/qt/tests && for t in tst_genomeformat tst_findings tst_texttable tst_vcfconvert \
		tst_targets tst_ensembl tst_clinvarprocess tst_analysisrunner tst_testgenome tst_chainlift tst_compression; do \
		QT_QPA_PLATFORM=offscreen ./$$t || exit 1; done

## Rebuild the reference-assembly samples (samples/): each is the build's own
## base at every position the analysis reads, a null test of the tables.
reference-samples: build
	./$(GH_DATA) test-genome --sample GRCh37 --output samples/genome_reference_GRCh37.txt
	./$(GH_DATA) test-genome --sample GRCh38 --output samples/genome_reference_GRCh38.txt

## chr22 truth-set validation of the WGS toolchain: simulate reads carrying
## known genotypes, run the pipeline, compare (~25s after a ~10MB setup).
validate-wgs-setup: build
	./$(GH_DATA) wgs-validate --setup

validate-wgs: build
	./$(GH_DATA) wgs-validate

## Print the install command for the WGS tools on this platform
tools:
	@echo "The WGS (FASTQ) pipeline needs: $(WGS_TOOLS) [+ fastp for QC]"
	@echo ""
	@if [ -f /etc/debian_version ]; then \
		echo "  sudo apt install minimap2 samtools bcftools fastp"; \
	elif [ "$$(uname)" = "Darwin" ]; then \
		echo "  brew install minimap2 samtools bcftools fastp"; \
	else \
		echo "  Install via your package manager"; \
	fi
	@echo ""
	@echo "The array-export path ('make analysis') needs none of these."

## Verify the WGS tools are on PATH
check-tools:
	@missing=""; \
	for t in $(WGS_TOOLS); do \
		command -v $$t >/dev/null 2>&1 || missing="$$missing $$t"; \
	done; \
	if [ -n "$$missing" ]; then \
		echo "Missing WGS tools:$$missing"; \
		echo "Run 'make tools' for the install command."; \
		exit 1; \
	fi; \
	echo "All WGS tools found: $(WGS_TOOLS)"
	@command -v fastp >/dev/null 2>&1 || echo "NOTE: fastp not found - the pipeline runs with --skip-qc."

## One-time setup for the FASTQ path: the GRCh37 reference + indexes (~12 GB),
## fetched and indexed by gh-data (also a button on the WGS and Data pages).
setup: build
	./$(GH_DATA) reference

## Remove intermediate files (keep BAM, VCF, and reports)
clean:
	rm -f $(WORK_DIR)/trimmed.fastq.gz
	rm -f $(WORK_DIR)/fastp.*
	@echo "Cleaned intermediate files"

## Remove all generated files (keep reference and data)
clean-all: clean
	rm -rf $(WORK_DIR)
	rm -f $(DATA_DIR)/genome.txt
	rm -f $(DATA_DIR)/genome.txt.bak
	rm -f $(RSID_LOOKUP)
	rm -rf $(REPORTS_DIR)/*.md $(REPORTS_DIR)/*.json
	rm -rf $(BUILD_DIR) $(ANALYSIS) $(GH_DATA) $(GUI)
	$(MAKE) -C c clean
	@echo "Cleaned all generated files"

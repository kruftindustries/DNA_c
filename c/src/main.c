/* genetic-health: command line entry point for the C port. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gh_analysis.h"
#include "gh_ancestry.h"
#include "gh_clinvar.h"
#include "gh_epistasis.h"
#include "gh_genome.h"
#include "gh_mem.h"
#include "gh_pharmgkb.h"
#include "gh_prs.h"
#include "gh_profiles.h"
#include "gh_dependent_profiles.h"
#include "gh_dosing.h"
#include "gh_insights.h"
#include "gh_recommendations.h"
#include "gh_quality.h"
#include "gh_report.h"
#include "gh_rsids.h"
#include "gh_scorers.h"
#include "gh_snpdb.h"
#include "gh_traits.h"

#ifndef GH_VERSION            /* set by engine.pro / c/Makefile from the release version */
#define GH_VERSION "0.2.5"
#endif

/* Everything the scorers produce, gathered for output. */
typedef struct {
    gh_apoe_result apoe;
    gh_blood_result blood;
    gh_mt_result mt;
    gh_star_result *stars;   /* GH_STAR_GENE_COUNT entries */
    gh_trait_result *traits; /* GH_TRAIT_COUNT entries */
    gh_ancestry_result ancestry;
    gh_prs_result *prs;      /* GH_PRS_MODEL_COUNT entries */
    gh_clinvar_result clinvar;
    gh_acmg_result acmg;
    gh_carrier_result carriers;
    gh_epi_result *epistasis;
    size_t nepistasis;
    gh_histamine_result histamine;
    gh_alcohol_result alcohol;
    gh_pain_result pain;
    gh_thyroid_result thyroid;
    gh_hormone_result hormone;
    gh_eye_result eye;
    gh_recommendations recs;
    gh_insights insights;
    gh_sleep_result sleep;
    gh_nutrigenomics_result nutrition;
    gh_mental_result mental;
    gh_longevity_result longevity;
    gh_dosing_result dosing;
    gh_polypharmacy_result polypharmacy;
    gh_preventive_result preventive;
} gh_scorer_results;

typedef struct {
    const char *genome;
    const char *data_dir;
    const char *name;
    const char *html;        /* output path, or NULL */
    const char *date;        /* report timestamp override, or NULL */
    bool json;
    bool quiet;
} options;

static void usage(FILE *out, const char *argv0)
{
    fprintf(out,
        "genetic-health " GH_VERSION " - genetic health analysis\n"
        "\n"
        "Usage: %s [options] [genome.txt]\n"
        "\n"
        "Options:\n"
        "  --data DIR     Directory holding the annotation data (default: data)\n"
        "  --name NAME    Subject name to label the output with\n"
        "  --html PATH    Write the HTML report to PATH (- for stdout)\n"
        "  --date TEXT    Timestamp to print in the report instead of now\n"
        "  --json         Emit findings as JSON instead of text\n"
        "  --list-rsids   Print every rsID the analysis reads, one per line\n"
        "  --quiet        Suppress progress output\n"
        "  -h, --help     Show this help\n"
        "  --version      Show the version\n"
        "\n"
        "The genome file defaults to DATA/genome.txt.\n",
        argv0);
}

/* Join a directory and filename into arena memory. */
static const char *path_join(gh_arena *a, const char *dir, const char *file)
{
    size_t n = strlen(dir) + 1 + strlen(file);
    char *p = gh_alloc(a, n + 1);
    snprintf(p, n + 1, "%s/%s", dir, file);
    return p;
}

static bool parse_args(options *o, int argc, char **argv)
{
    o->genome = NULL;
    o->data_dir = "data";
    o->name = NULL;
    o->html = NULL;
    o->date = NULL;
    o->json = false;
    o->quiet = false;

    for (int i = 1; i < argc; i++) {
        const char *arg = argv[i];
        if (strcmp(arg, "-h") == 0 || strcmp(arg, "--help") == 0) {
            usage(stdout, argv[0]);
            exit(0);
        } else if (strcmp(arg, "--version") == 0) {
            puts("genetic-health " GH_VERSION);
            exit(0);
        } else if (strcmp(arg, "--list-rsids") == 0) {
            gh_arena *arena = gh_arena_new(1 << 20);
            const char **rsids;
            size_t n = gh_all_rsids(arena, &rsids);
            for (size_t k = 0; k < n; k++)
                puts(rsids[k]);
            gh_arena_free(arena);
            exit(0);
        } else if (strcmp(arg, "--json") == 0) {
            o->json = true;
        } else if (strcmp(arg, "--quiet") == 0) {
            o->quiet = true;
        } else if (strcmp(arg, "--data") == 0 || strcmp(arg, "--name") == 0 ||
                   strcmp(arg, "--html") == 0 || strcmp(arg, "--date") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "error: %s needs a value\n", arg);
                return false;
            }
            if (strcmp(arg, "--data") == 0)
                o->data_dir = argv[++i];
            else if (strcmp(arg, "--name") == 0)
                o->name = argv[++i];
            else if (strcmp(arg, "--date") == 0)
                o->date = argv[++i];
            else
                o->html = argv[++i];
        } else if (arg[0] == '-' && arg[1] != '\0') {
            fprintf(stderr, "error: unknown option %s\n", arg);
            return false;
        } else if (!o->genome) {
            o->genome = arg;
        } else {
            fprintf(stderr, "error: unexpected argument %s\n", arg);
            return false;
        }
    }
    return true;
}

/* Escape a string as a JSON string value, quotes included. */
static void json_string(gh_strbuf *sb, const char *s)
{
    gh_sb_putc(sb, '"');
    for (const char *p = s ? s : ""; *p; p++) {
        unsigned char c = (unsigned char)*p;
        switch (c) {
        case '"':  gh_sb_puts(sb, "\\\""); break;
        case '\\': gh_sb_puts(sb, "\\\\"); break;
        case '\n': gh_sb_puts(sb, "\\n");  break;
        case '\r': gh_sb_puts(sb, "\\r");  break;
        case '\t': gh_sb_puts(sb, "\\t");  break;
        default:
            if (c < 0x20)
                gh_sb_printf(sb, "\\u%04x", c);
            else
                gh_sb_putc(sb, (char)c);
        }
    }
    gh_sb_putc(sb, '"');
}

/* Fields every profile shares. Leaves the object open so each caller can
 * append its own. */
static void emit_profile_base(gh_strbuf *sb, const char *key,
                              const gh_profile_base *base)
{
    gh_sb_puts(sb, "  ");
    json_string(sb, key);
    gh_sb_printf(sb, ": {\"markers_found\": %zu, \"markers_tested\": %zu, "
                     "\"summary\": ", base->nhits, base->markers_tested);
    json_string(sb, base->summary);
    gh_sb_puts(sb, ", \"recommendations\": [");
    for (size_t i = 0; i < base->nrecommendations; i++) {
        if (i)
            gh_sb_puts(sb, ", ");
        json_string(sb, base->recommendations[i]);
    }
    gh_sb_puts(sb, "]");
}

static void emit_json(const gh_analysis *a, const gh_quality *q,
                      const gh_scorer_results *sc, const options *o)
{
    gh_strbuf sb;
    gh_sb_init(&sb);

    gh_sb_puts(&sb, "{\n  \"name\": ");
    json_string(&sb, o->name ? o->name : "");
    gh_sb_printf(&sb,
        ",\n  \"summary\": {\n"
        "    \"total_snps\": %zu,\n"
        "    \"analyzed_snps\": %zu,\n"
        "    \"high_impact\": %zu,\n"
        "    \"moderate_impact\": %zu,\n"
        "    \"low_impact\": %zu\n  },\n",
        a->total_snps, a->analyzed_snps, a->high_impact,
        a->moderate_impact, a->low_impact);

    gh_sb_puts(&sb, "  \"findings\": [\n");
    for (size_t i = 0; i < a->nfindings; i++) {
        const gh_finding *f = &a->findings[i];
        gh_sb_puts(&sb, "    {\"rsid\": ");
        json_string(&sb, f->rsid);
        gh_sb_puts(&sb, ", \"gene\": ");
        json_string(&sb, f->gene);
        gh_sb_puts(&sb, ", \"category\": ");
        json_string(&sb, f->category);
        gh_sb_puts(&sb, ", \"genotype\": ");
        json_string(&sb, f->genotype);
        gh_sb_puts(&sb, ", \"status\": ");
        json_string(&sb, f->status);
        gh_sb_puts(&sb, ", \"description\": ");
        json_string(&sb, f->description);
        /* note and freq both reach the rendered report, so both have to be
         * comparable against the Python -- they were emitted nowhere and
         * compared nowhere, which left two whole finding fields untested. */
        gh_sb_puts(&sb, ", \"note\": ");
        json_string(&sb, f->note);
        gh_sb_puts(&sb, ", \"freq\": ");
        if (f->freq) {
            gh_sb_puts(&sb, "{");
            for (size_t k = 0; k < GH_POPULATION_COUNT; k++) {
                if (k)
                    gh_sb_puts(&sb, ", ");
                json_string(&sb, GH_POPULATIONS[k]);
                gh_sb_printf(&sb, ": %.12g", f->freq[k]);
            }
            gh_sb_puts(&sb, "}");
        } else {
            gh_sb_puts(&sb, "null");
        }
        gh_sb_printf(&sb, ", \"magnitude\": %d}", f->magnitude);
        gh_sb_puts(&sb, i + 1 < a->nfindings ? ",\n" : "\n");
    }
    gh_sb_puts(&sb, "  ],\n  \"drug_findings\": [\n");
    for (size_t i = 0; i < a->ndrug_findings; i++) {
        const gh_drug_finding *d = &a->drug_findings[i];
        gh_sb_puts(&sb, "    {\"rsid\": ");
        json_string(&sb, d->rsid);
        gh_sb_puts(&sb, ", \"gene\": ");
        json_string(&sb, d->gene);
        gh_sb_puts(&sb, ", \"drugs\": ");
        json_string(&sb, d->drugs);
        gh_sb_puts(&sb, ", \"genotype\": ");
        json_string(&sb, d->genotype);
        gh_sb_puts(&sb, ", \"level\": ");
        json_string(&sb, d->level);
        gh_sb_puts(&sb, ", \"annotation\": ");
        json_string(&sb, d->annotation);
        gh_sb_putc(&sb, '}');
        gh_sb_puts(&sb, i + 1 < a->ndrug_findings ? ",\n" : "\n");
    }
    gh_sb_puts(&sb, "  ],\n  \"quality\": {\n");
    gh_sb_printf(&sb,
        "    \"total_snps\": %zu,\n"
        "    \"no_call_count\": %zu,\n"
        "    \"call_rate\": %.12g,\n"
        "    \"autosomal_count\": %zu,\n"
        "    \"mt_snp_count\": %zu,\n"
        "    \"has_mt\": %s,\n"
        "    \"has_y\": %s,\n"
        "    \"het_rate\": %.12g,\n",
        q->total_snps, q->no_call_count, q->call_rate, q->autosomal_count,
        q->mt_snp_count, q->has_mt ? "true" : "false",
        q->has_y ? "true" : "false", q->het_rate);
    gh_sb_puts(&sb, "    \"chromosomes\": {");
    for (size_t i = 0; i < q->nchromosomes; i++) {
        if (i)
            gh_sb_putc(&sb, ',');
        gh_sb_puts(&sb, "\n      ");
        json_string(&sb, q->chromosomes[i].name);
        gh_sb_printf(&sb, ": %zu", q->chromosomes[i].count);
    }
    gh_sb_puts(&sb, "\n    }\n  },\n");

    /* APOE */
    gh_sb_puts(&sb, "  \"apoe\": {\"apoe_type\": ");
    json_string(&sb, sc->apoe.apoe_type);
    gh_sb_puts(&sb, ", \"risk_level\": ");
    json_string(&sb, sc->apoe.risk_level);
    gh_sb_puts(&sb, ", \"confidence\": ");
    json_string(&sb, sc->apoe.confidence);
    gh_sb_puts(&sb, ", \"description\": ");
    json_string(&sb, sc->apoe.description);
    if (sc->apoe.has_or)
        gh_sb_printf(&sb, ", \"alzheimer_or\": %.12g},\n", sc->apoe.alzheimer_or);
    else
        gh_sb_puts(&sb, ", \"alzheimer_or\": null},\n");

    /* Blood type */
    gh_sb_puts(&sb, "  \"blood_type\": {\"blood_type\": ");
    json_string(&sb, sc->blood.blood_type);
    gh_sb_puts(&sb, ", \"abo\": ");
    json_string(&sb, sc->blood.abo);
    gh_sb_puts(&sb, ", \"rh\": ");
    json_string(&sb, sc->blood.rh);
    gh_sb_puts(&sb, ", \"confidence\": ");
    json_string(&sb, sc->blood.confidence);
    gh_sb_puts(&sb, "},\n");

    /* Mitochondrial haplogroup */
    gh_sb_puts(&sb, "  \"mt_haplogroup\": {\"haplogroup\": ");
    json_string(&sb, sc->mt.haplogroup);
    gh_sb_puts(&sb, ", \"description\": ");
    json_string(&sb, sc->mt.description);
    gh_sb_puts(&sb, ", \"confidence\": ");
    json_string(&sb, sc->mt.confidence);
    gh_sb_puts(&sb, ", \"lineage\": ");
    json_string(&sb, sc->mt.lineage);
    gh_sb_printf(&sb, ", \"markers_found\": %zu, \"markers_tested\": %zu},\n",
                 sc->mt.markers_found, sc->mt.markers_tested);

    /* Star alleles */
    gh_sb_puts(&sb, "  \"star_alleles\": {");
    for (size_t i = 0; i < GH_STAR_GENE_COUNT; i++) {
        const gh_star_result *r = &sc->stars[i];
        if (i)
            gh_sb_putc(&sb, ',');
        gh_sb_puts(&sb, "\n    ");
        json_string(&sb, r->gene);
        gh_sb_puts(&sb, ": {\"gene\": ");
        json_string(&sb, r->gene);
        gh_sb_puts(&sb, ", \"diplotype\": ");
        json_string(&sb, r->diplotype);
        gh_sb_puts(&sb, ", \"phenotype\": ");
        json_string(&sb, r->phenotype);
        gh_sb_puts(&sb, ", \"confidence\": ");
        json_string(&sb, r->confidence);
        gh_sb_puts(&sb, ", \"clinical_note\": ");
        json_string(&sb, r->clinical_note);
        gh_sb_printf(&sb,
            ", \"snps_found\": %zu, \"snps_total\": %zu, \"coverage\": %.12g}",
            r->snps_found, r->snps_total, r->coverage);
    }
    gh_sb_puts(&sb, "\n  },\n");

    /* Traits */
    gh_sb_puts(&sb, "  \"traits\": {");
    for (size_t i = 0; i < GH_TRAIT_COUNT; i++) {
        const gh_trait_result *t = &sc->traits[i];
        if (!t->key)
            continue;
        if (i)
            gh_sb_putc(&sb, ',');
        gh_sb_puts(&sb, "\n    ");
        json_string(&sb, t->key);
        gh_sb_puts(&sb, ": {\"prediction\": ");
        json_string(&sb, t->prediction);
        gh_sb_puts(&sb, ", \"confidence\": ");
        json_string(&sb, t->confidence);
        gh_sb_puts(&sb, ", \"description\": ");
        json_string(&sb, t->description);
        gh_sb_puts(&sb, ", \"snps_used\": [");
        for (size_t j = 0; j < t->nsnps; j++) {
            if (j)
                gh_sb_puts(&sb, ", ");
            json_string(&sb, t->snps_used[j]);
        }
        gh_sb_puts(&sb, "]}");
    }
    gh_sb_puts(&sb, "\n  },\n");

    /* Ancestry */
    gh_sb_puts(&sb, "  \"ancestry\": {\"proportions\": {");
    for (size_t i = 0; i < GH_POPULATION_COUNT; i++) {
        if (i)
            gh_sb_putc(&sb, ',');
        gh_sb_puts(&sb, "\n      ");
        json_string(&sb, GH_POPULATIONS[i]);
        gh_sb_printf(&sb, ": %.12g", sc->ancestry.proportions[i]);
    }
    gh_sb_printf(&sb, "\n    }, \"markers_found\": %zu, \"confidence\": ",
                 sc->ancestry.markers_found);
    json_string(&sb, sc->ancestry.confidence);
    gh_sb_puts(&sb, ", \"top_ancestry\": ");
    json_string(&sb, sc->ancestry.top_ancestry);
    if (sc->ancestry.sub.present) {
        gh_sb_puts(&sb, ", \"sub_ancestry\": {\"top_sub_ancestry\": ");
        json_string(&sb, sc->ancestry.sub.top_label);
        gh_sb_puts(&sb, ", \"confidence\": ");
        json_string(&sb, sc->ancestry.sub.confidence);
        gh_sb_printf(&sb, ", \"markers_used\": %zu, \"sub_proportions\": {",
                     sc->ancestry.sub.markers_used);
        for (size_t i = 0; i < sc->ancestry.sub.nsubs; i++) {
            if (i)
                gh_sb_putc(&sb, ',');
            gh_sb_puts(&sb, "\n        ");
            json_string(&sb, sc->ancestry.sub.labels[i]);
            gh_sb_printf(&sb, ": %.12g", sc->ancestry.sub.proportions[i]);
        }
        gh_sb_puts(&sb, "\n      }}");
    } else {
        gh_sb_puts(&sb, ", \"sub_ancestry\": null");
    }
    gh_sb_puts(&sb, "},\n");

    /* Polygenic risk scores */
    gh_sb_puts(&sb, "  \"prs\": {");
    for (size_t i = 0; i < GH_PRS_MODEL_COUNT; i++) {
        const gh_prs_result *r = &sc->prs[i];
        if (i)
            gh_sb_putc(&sb, ',');
        gh_sb_puts(&sb, "\n    ");
        json_string(&sb, r->id);
        gh_sb_puts(&sb, ": {\"name\": ");
        json_string(&sb, r->name);
        gh_sb_printf(&sb,
            ", \"raw_score\": %.12g, \"z_score\": %.12g, "
            "\"percentile\": %.12g, \"ci_95_lower\": %.12g, "
            "\"ci_95_upper\": %.12g, \"snps_found\": %zu, "
            "\"snps_total\": %zu, \"ancestry_applicable\": %s, "
            "\"risk_category\": ",
            r->raw_score, r->z_score, r->percentile, r->ci_95_lower,
            r->ci_95_upper, r->snps_found, r->snps_total,
            r->ancestry_applicable ? "true" : "false");
        json_string(&sb, r->risk_category);
        gh_sb_puts(&sb, ", \"ancestry_warning\": ");
        json_string(&sb, r->ancestry_warning);
        gh_sb_puts(&sb, ", \"reference\": ");
        json_string(&sb, r->reference);
        gh_sb_puts(&sb, ", \"contributing_snps\": [");
        for (size_t j = 0; j < r->ncontributing; j++) {
            const gh_prs_contribution *c = &r->contributing[j];
            if (j)
                gh_sb_puts(&sb, ", ");
            gh_sb_puts(&sb, "{\"rsid\": ");
            json_string(&sb, c->rsid);
            gh_sb_puts(&sb, ", \"gene\": ");
            json_string(&sb, c->gene);
            gh_sb_printf(&sb, ", \"copies\": %d, \"contribution\": %.12g}",
                         c->copies, c->contribution);
        }
        gh_sb_puts(&sb, "]}");
    }
    gh_sb_puts(&sb, "\n  },\n");

    /* ClinVar disease findings */
    gh_sb_puts(&sb, "  \"disease_findings\": ");
    if (!sc->clinvar.loaded) {
        gh_sb_puts(&sb, "null,\n");
    } else {
        gh_sb_puts(&sb, "{");
        for (size_t c = 0; c < GH_CV_CATEGORY_COUNT; c++) {
            if (c)
                gh_sb_putc(&sb, ',');
            gh_sb_puts(&sb, "\n    ");
            json_string(&sb, GH_CV_CATEGORY_NAMES[c]);
            gh_sb_puts(&sb, ": [");
            for (size_t i = 0; i < sc->clinvar.counts[c]; i++) {
                const gh_cv_finding *f = &sc->clinvar.by_category[c][i];
                if (i)
                    gh_sb_puts(&sb, ", ");
                gh_sb_puts(&sb, "{\"rsid\": ");
                json_string(&sb, f->rsid);
                gh_sb_puts(&sb, ", \"gene\": ");
                json_string(&sb, f->gene);
                gh_sb_puts(&sb, ", \"chromosome\": ");
                json_string(&sb, f->chromosome);
                gh_sb_puts(&sb, ", \"position\": ");
                json_string(&sb, f->position);
                gh_sb_puts(&sb, ", \"user_genotype\": ");
                json_string(&sb, f->user_genotype);
                gh_sb_puts(&sb, ", \"zygosity_status\": ");
                json_string(&sb, f->zygosity_status);
                gh_sb_puts(&sb, ", \"zygosity\": ");
                json_string(&sb, f->zygosity);
                gh_sb_puts(&sb, ", \"clinical_significance\": ");
                json_string(&sb, f->clinical_significance);
                gh_sb_printf(&sb, ", \"gold_stars\": %d, \"is_homozygous\": %s}",
                             f->gold_stars, f->is_homozygous ? "true" : "false");
            }
            gh_sb_puts(&sb, "]");
        }
        gh_sb_printf(&sb, ",\n    \"_stats\": {\"total_clinvar\": %zu, "
                          "\"matched\": %zu}\n  },\n",
                     sc->clinvar.total_clinvar, sc->clinvar.matched);
    }

    /* ACMG secondary findings */
    gh_sb_printf(&sb,
        "  \"acmg\": {\"genes_screened\": %zu, \"genes_with_variants\": %zu, "
        "\"summary\": ", sc->acmg.genes_screened, sc->acmg.genes_with_variants);
    json_string(&sb, sc->acmg.summary);
    gh_sb_puts(&sb, ", \"acmg_findings\": [");
    for (size_t i = 0; i < sc->acmg.nfindings; i++) {
        const gh_acmg_finding *e = &sc->acmg.findings[i];
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"gene\": ");
        json_string(&sb, e->finding->gene);
        gh_sb_puts(&sb, ", \"acmg_category\": ");
        json_string(&sb, e->acmg_category);
        gh_sb_puts(&sb, ", \"acmg_actionability\": ");
        json_string(&sb, e->actionability);
        gh_sb_printf(&sb, ", \"gold_stars\": %d}", e->finding->gold_stars);
    }
    gh_sb_puts(&sb, "]},\n");

    /* Carrier screening */
    gh_sb_printf(&sb, "  \"carrier_screen\": {\"total_carriers\": %zu, "
                      "\"carriers\": [", sc->carriers.ncarriers);
    for (size_t i = 0; i < sc->carriers.ncarriers; i++) {
        const gh_carrier *c = &sc->carriers.carriers[i];
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"gene\": ");
        json_string(&sb, c->gene);
        gh_sb_puts(&sb, ", \"condition\": ");
        json_string(&sb, c->condition);
        gh_sb_puts(&sb, ", \"inheritance\": ");
        json_string(&sb, c->inheritance);
        gh_sb_puts(&sb, ", \"system\": ");
        json_string(&sb, c->system);
        gh_sb_puts(&sb, ", \"reproductive_note\": ");
        json_string(&sb, c->reproductive_note);
        gh_sb_puts(&sb, ", \"rsid\": ");
        json_string(&sb, c->rsid);
        gh_sb_puts(&sb, ", \"genotype\": ");
        json_string(&sb, c->genotype);
        gh_sb_printf(&sb, ", \"couples_relevant\": %s, \"gold_stars\": %d}",
                     c->couples_relevant ? "true" : "false", c->gold_stars);
    }
    gh_sb_puts(&sb, "], \"by_system\": {");
    for (size_t i = 0; i < sc->carriers.nsystems; i++) {
        if (i)
            gh_sb_putc(&sb, ',');
        gh_sb_puts(&sb, "\n      ");
        json_string(&sb, sc->carriers.systems[i]);
        gh_sb_printf(&sb, ": %zu", sc->carriers.system_counts[i]);
    }
    gh_sb_puts(&sb, "\n    }},\n");

    /* Gene-gene interactions */
    gh_sb_puts(&sb, "  \"epistasis\": [");
    for (size_t i = 0; i < sc->nepistasis; i++) {
        const gh_epi_result *e = &sc->epistasis[i];
        if (i)
            gh_sb_putc(&sb, ',');
        gh_sb_puts(&sb, "\n    {\"id\": ");
        json_string(&sb, e->id);
        gh_sb_puts(&sb, ", \"name\": ");
        json_string(&sb, e->name);
        gh_sb_puts(&sb, ", \"risk_level\": ");
        json_string(&sb, e->risk_level);
        gh_sb_printf(&sb, ", \"severity_score\": %.12g", e->severity_score);
        gh_sb_puts(&sb, ", \"genes_involved\": {");
        for (size_t j = 0; j < e->ngenes; j++) {
            if (j)
                gh_sb_puts(&sb, ", ");
            json_string(&sb, e->genes[j].gene);
            gh_sb_puts(&sb, ": [");
            for (size_t k = 0; k < e->genes[j].nstatuses; k++) {
                if (k)
                    gh_sb_puts(&sb, ", ");
                json_string(&sb, e->genes[j].statuses[k]);
            }
            gh_sb_puts(&sb, "]");
        }
        gh_sb_puts(&sb, "}}");
    }
    gh_sb_puts(&sb, "\n  ],\n");

    /* Health profiles */
    emit_profile_base(&sb, "histamine", &sc->histamine.base);
    gh_sb_puts(&sb, ", \"risk_level\": ");
    json_string(&sb, sc->histamine.risk_level);
    gh_sb_puts(&sb, "},\n");

    emit_profile_base(&sb, "alcohol_profile", &sc->alcohol.base);
    gh_sb_puts(&sb, ", \"metabolism_speed\": ");
    json_string(&sb, sc->alcohol.metabolism_speed);
    gh_sb_puts(&sb, ", \"flush_risk\": ");
    json_string(&sb, sc->alcohol.flush_risk);
    gh_sb_puts(&sb, ", \"cancer_risk\": ");
    json_string(&sb, sc->alcohol.cancer_risk);
    gh_sb_puts(&sb, "},\n");

    emit_profile_base(&sb, "pain_sensitivity", &sc->pain.base);
    gh_sb_printf(&sb, ", \"pain_sensitivity_score\": %d},\n",
                 sc->pain.sensitivity_score);

    emit_profile_base(&sb, "thyroid", &sc->thyroid.base);
    gh_sb_puts(&sb, ", \"risk_profile\": {");
    for (size_t i = 0; i < sc->thyroid.ndomains; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, sc->thyroid.domains[i].name);
        gh_sb_puts(&sb, ": ");
        json_string(&sb, sc->thyroid.domains[i].level);
    }
    gh_sb_puts(&sb, "}},\n");

    emit_profile_base(&sb, "hormone_metabolism", &sc->hormone.base);
    gh_sb_puts(&sb, ", \"estrogen_level\": ");
    json_string(&sb, sc->hormone.estrogen_level);
    gh_sb_puts(&sb, ", \"androgen_level\": ");
    json_string(&sb, sc->hormone.androgen_level);
    gh_sb_puts(&sb, ", \"overall\": ");
    json_string(&sb, sc->hormone.overall);
    gh_sb_puts(&sb, "},\n");

    emit_profile_base(&sb, "eye_health", &sc->eye.base);
    gh_sb_puts(&sb, ", \"conditions\": {");
    for (size_t i = 0; i < sc->eye.nconditions; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, sc->eye.conditions[i].name);
        gh_sb_puts(&sb, ": ");
        json_string(&sb, sc->eye.conditions[i].level);
    }
    gh_sb_puts(&sb, "}},\n");

    /* Recommendations */
    const gh_recommendations *r = &sc->recs;
    gh_sb_puts(&sb, "  \"recommendations\": {\"priorities\": [");
    for (size_t i = 0; i < r->npriorities; i++) {
        const gh_priority *p = &r->priorities[i];
        if (i)
            gh_sb_putc(&sb, ',');
        gh_sb_puts(&sb, "\n    {\"id\": ");
        json_string(&sb, p->id);
        gh_sb_puts(&sb, ", \"title\": ");
        json_string(&sb, p->title);
        gh_sb_puts(&sb, ", \"priority\": ");
        json_string(&sb, p->priority);
        gh_sb_puts(&sb, ", \"why\": ");
        json_string(&sb, p->why);
        gh_sb_printf(&sb, ", \"signal_count\": %zu, \"clinical_actions\": [",
                     p->signal_count);
        for (size_t k = 0; k < p->nclinical_actions; k++) {
            if (k)
                gh_sb_puts(&sb, ", ");
            json_string(&sb, p->clinical_actions[k]);
        }
        gh_sb_puts(&sb, "]}");
    }
    gh_sb_puts(&sb, "],\n    \"monitoring_schedule\": [");
    for (size_t i = 0; i < r->nmonitoring; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"test\": ");
        json_string(&sb, r->monitoring_schedule[i].test);
        gh_sb_puts(&sb, ", \"frequency\": ");
        json_string(&sb, r->monitoring_schedule[i].frequency);
        gh_sb_puts(&sb, "}");
    }
    gh_sb_puts(&sb, "],\n    \"good_news\": [");
    for (size_t i = 0; i < r->ngood_news; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"gene\": ");
        json_string(&sb, r->good_news[i].gene);
        gh_sb_puts(&sb, ", \"description\": ");
        json_string(&sb, r->good_news[i].description);
        gh_sb_puts(&sb, "}");
    }
    gh_sb_puts(&sb, "],\n    \"specialist_referrals\": [");
    for (size_t i = 0; i < r->nreferrals; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"specialist\": ");
        json_string(&sb, r->referrals[i].specialist);
        gh_sb_puts(&sb, ", \"reason\": ");
        json_string(&sb, r->referrals[i].reason);
        gh_sb_puts(&sb, ", \"urgency\": ");
        json_string(&sb, r->referrals[i].urgency);
        gh_sb_puts(&sb, "}");
    }
    gh_sb_puts(&sb, "],\n    \"drug_card\": [");
    for (size_t i = 0; i < r->ndrug_card; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"gene\": ");
        json_string(&sb, r->drug_card[i].gene);
        gh_sb_printf(&sb, ", \"entries\": %zu}", r->drug_card[i].nentries);
    }
    gh_sb_puts(&sb, "],\n    \"clinical_insights\": [");
    for (size_t i = 0; i < r->ninsights; i++) {
        const gh_clinical_insight *ins = &r->insights[i];
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"gene\": ");
        json_string(&sb, ins->gene);
        gh_sb_puts(&sb, ", \"status\": ");
        json_string(&sb, ins->status);
        gh_sb_printf(&sb, ", \"magnitude\": %d, \"pathways\": [",
                     ins->magnitude);
        for (size_t k = 0; k < ins->npathways; k++) {
            if (k)
                gh_sb_puts(&sb, ", ");
            json_string(&sb, ins->pathways[k]);
        }
        gh_sb_puts(&sb, "]}");
    }
    gh_sb_puts(&sb, "]\n  },\n");

    /* Insights */
    const gh_insights *ins = &sc->insights;
    gh_sb_puts(&sb, "  \"insights\": {\"single_gene\": [");
    for (size_t i = 0; i < ins->nsingle_gene; i++) {
        const gh_single_gene_result *sg = &ins->single_gene[i];
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"gene\": ");
        json_string(&sb, sg->gene);
        gh_sb_puts(&sb, ", \"status\": ");
        json_string(&sb, sg->status);
        gh_sb_printf(&sb, ", \"magnitude\": %d, \"title\": ", sg->magnitude);
        json_string(&sb, sg->entry->title);
        gh_sb_puts(&sb, "}");
    }
    gh_sb_puts(&sb, "],\n    \"narratives\": [");
    for (size_t i = 0; i < ins->nnarratives; i++) {
        const gh_narrative_result *nr = &ins->narratives[i];
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"id\": ");
        json_string(&sb, nr->id);
        gh_sb_puts(&sb, ", \"title\": ");
        json_string(&sb, nr->title);
        gh_sb_puts(&sb, ", \"narrative\": ");
        json_string(&sb, nr->narrative);
        gh_sb_puts(&sb, ", \"matched_genes\": [");
        for (size_t k = 0; k < nr->nmatched; k++) {
            if (k)
                gh_sb_puts(&sb, ", ");
            json_string(&sb, nr->matched_genes[k]);
        }
        gh_sb_puts(&sb, "]}");
    }
    gh_sb_puts(&sb, "],\n    \"genome_highlights\": [");
    for (size_t i = 0; i < ins->nhighlights; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"title\": ");
        json_string(&sb, ins->highlights[i].title);
        gh_sb_puts(&sb, ", \"detail\": ");
        json_string(&sb, ins->highlights[i].detail);
        gh_sb_puts(&sb, ", \"type\": ");
        json_string(&sb, ins->highlights[i].type);
        gh_sb_puts(&sb, "}");
    }
    gh_sb_puts(&sb, "],\n    \"protective_findings\": [");
    for (size_t i = 0; i < ins->nprotective; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"gene\": ");
        json_string(&sb, ins->protective[i].gene);
        gh_sb_puts(&sb, ", \"status\": ");
        json_string(&sb, ins->protective[i].status);
        gh_sb_puts(&sb, ", \"title\": ");
        json_string(&sb, ins->protective[i].title);
        gh_sb_puts(&sb, "}");
    }
    gh_sb_puts(&sb, "]\n  },\n");

    /* Sleep */
    gh_sb_puts(&sb, "  \"sleep_profile\": {\"chronotype\": ");
    json_string(&sb, sc->sleep.chronotype);
    gh_sb_printf(&sb, ", \"chronotype_score\": %.12g, \"markers_found\": %zu, "
                      "\"confidence\": ", sc->sleep.chronotype_score,
                 sc->sleep.markers_found);
    json_string(&sb, sc->sleep.confidence);
    gh_sb_puts(&sb, ", \"optimal_sleep_window\": ");
    json_string(&sb, sc->sleep.optimal_sleep_window);
    gh_sb_puts(&sb, ", \"caffeine_cutoff\": ");
    json_string(&sb, sc->sleep.caffeine_cutoff);
    gh_sb_puts(&sb, ", \"peak_alertness\": ");
    json_string(&sb, sc->sleep.peak_alertness);
    gh_sb_puts(&sb, ", \"deep_sleep_note\": ");
    json_string(&sb, sc->sleep.deep_sleep_note);
    gh_sb_printf(&sb, ", \"caffeine_sensitive\": %s, \"recommendations\": [",
                 sc->sleep.caffeine_sensitive ? "true" : "false");
    for (size_t i = 0; i < sc->sleep.nrecommendations; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, sc->sleep.recommendations[i]);
    }
    gh_sb_puts(&sb, "]},\n");

    /* Nutrigenomics */
    gh_sb_puts(&sb, "  \"nutrigenomics\": {\"summary\": ");
    json_string(&sb, sc->nutrition.summary);
    gh_sb_puts(&sb, ", \"nutrient_needs\": [");
    for (size_t i = 0; i < sc->nutrition.nneeds; i++) {
        const gh_nutrient_need *n = &sc->nutrition.needs[i];
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"id\": ");
        json_string(&sb, n->profile->id);
        gh_sb_puts(&sb, ", \"need_level\": ");
        json_string(&sb, n->need_level);
        gh_sb_printf(&sb, ", \"severity\": %.12g, \"recommendation\": ",
                     n->severity);
        json_string(&sb, n->recommendation);
        gh_sb_puts(&sb, "}");
    }
    gh_sb_puts(&sb, "]},\n");

    /* Mental health */
    gh_sb_puts(&sb, "  \"mental_health\": {\"summary\": ");
    json_string(&sb, sc->mental.summary);
    gh_sb_puts(&sb, ", \"domains\": {");
    for (size_t i = 0; i < sc->mental.ndomains; i++) {
        const gh_mental_domain_result *d = &sc->mental.domains[i];
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, d->name);
        gh_sb_puts(&sb, ": {\"risk_level\": ");
        json_string(&sb, d->risk_level);
        gh_sb_printf(&sb, ", \"risk_score\": %d, \"signals\": [", d->risk_score);
        for (size_t k = 0; k < d->nsignals; k++) {
            if (k)
                gh_sb_puts(&sb, ", ");
            json_string(&sb, d->signals[k]);
        }
        gh_sb_puts(&sb, "]}");
    }
    gh_sb_puts(&sb, "}, \"risk_factors\": [");
    for (size_t i = 0; i < sc->mental.nrisk_factors; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, sc->mental.risk_factors[i]);
    }
    gh_sb_puts(&sb, "], \"resilience_factors\": [");
    for (size_t i = 0; i < sc->mental.nresilience_factors; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, sc->mental.resilience_factors[i]);
    }
    gh_sb_puts(&sb, "], \"treatment_notes\": [");
    for (size_t i = 0; i < sc->mental.ntreatment_notes; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, sc->mental.treatment_notes[i]);
    }
    gh_sb_puts(&sb, "], \"recommendations\": [");
    for (size_t i = 0; i < sc->mental.nrecommendations; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, sc->mental.recommendations[i]);
    }
    gh_sb_puts(&sb, "]},\n");

    /* Longevity */
    gh_sb_printf(&sb, "  \"longevity\": {\"longevity_score\": %.12g, "
                      "\"alleles_checked\": %zu, \"summary\": ",
                 sc->longevity.longevity_score, sc->longevity.alleles_checked);
    json_string(&sb, sc->longevity.summary);
    gh_sb_puts(&sb, ", \"healthspan_domains\": {");
    for (size_t i = 0; i < sc->longevity.ndomains; i++) {
        const gh_healthspan_result *d = &sc->longevity.domains[i];
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, d->name);
        gh_sb_printf(&sb, ": {\"score\": %d, \"genes_found\": %zu, "
                          "\"rating\": ", d->score, d->genes_found);
        json_string(&sb, d->rating);
        gh_sb_puts(&sb, "}");
    }
    gh_sb_puts(&sb, "}, \"top_risks\": [");
    for (size_t i = 0; i < sc->longevity.ntop_risks; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, sc->longevity.top_risks[i]);
    }
    gh_sb_puts(&sb, "], \"top_protective\": [");
    for (size_t i = 0; i < sc->longevity.ntop_protective; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, sc->longevity.top_protective[i]);
    }
    gh_sb_puts(&sb, "], \"interventions\": [");
    for (size_t i = 0; i < sc->longevity.ninterventions; i++) {
        const gh_intervention *iv = &sc->longevity.interventions[i];
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"intervention\": ");
        json_string(&sb, iv->intervention);
        gh_sb_puts(&sb, ", \"genetic_support\": ");
        json_string(&sb, iv->genetic_support);
        gh_sb_puts(&sb, "}");
    }
    gh_sb_puts(&sb, "]},\n");

    /* Drug dosing */
    gh_sb_puts(&sb, "  \"drug_dosing\": {\"summary\": ");
    json_string(&sb, sc->dosing.summary);
    gh_sb_puts(&sb, ", \"recommendations\": [");
    for (size_t i = 0; i < sc->dosing.nrecommendations; i++) {
        const gh_dose_rec *rec = &sc->dosing.recommendations[i];
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"drug\": ");
        json_string(&sb, rec->drug->display);
        gh_sb_puts(&sb, ", \"category\": ");
        json_string(&sb, rec->drug->category);
        gh_sb_puts(&sb, ", \"genes\": [");
        for (size_t g = 0; g < rec->drug->ngenes; g++) {
            if (g)
                gh_sb_puts(&sb, ", ");
            json_string(&sb, rec->drug->genes[g]);
        }
        gh_sb_puts(&sb, "], \"action\": ");
        json_string(&sb, rec->rule->action);
        gh_sb_puts(&sb, ", \"dose_guidance\": ");
        json_string(&sb, rec->rule->dose_guidance);
        gh_sb_puts(&sb, ", \"source\": ");
        json_string(&sb, rec->drug->source);
        gh_sb_puts(&sb, "}");
    }
    gh_sb_puts(&sb, "], \"warnings\": [");
    for (size_t i = 0; i < sc->dosing.nwarnings; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, sc->dosing.warnings[i]->rule->action);
    }
    gh_sb_puts(&sb, "]},\n");

    /* Polypharmacy */
    gh_sb_printf(&sb, "  \"polypharmacy\": {\"total_warnings\": %zu, "
                      "\"warnings\": [", sc->polypharmacy.nwarnings);
    for (size_t i = 0; i < sc->polypharmacy.nwarnings; i++) {
        const gh_poly_warning *w = &sc->polypharmacy.warnings[i];
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"id\": ");
        json_string(&sb, w->rule->id);
        gh_sb_puts(&sb, ", \"name\": ");
        json_string(&sb, w->rule->name);
        gh_sb_puts(&sb, ", \"severity\": ");
        json_string(&sb, w->rule->severity);
        gh_sb_puts(&sb, ", \"matched_genes\": {");
        for (size_t g = 0; g < w->rule->ngenes; g++) {
            if (g)
                gh_sb_puts(&sb, ", ");
            json_string(&sb, w->rule->genes[g].gene);
            gh_sb_puts(&sb, ": ");
            json_string(&sb, w->matched[g]);
        }
        gh_sb_puts(&sb, "}, \"drugs_affected\": [");
        for (size_t d = 0; d < w->rule->ndrugs; d++) {
            if (d)
                gh_sb_puts(&sb, ", ");
            json_string(&sb, w->rule->drugs_affected[d]);
        }
        gh_sb_puts(&sb, "], \"warning\": ");
        json_string(&sb, w->rule->warning);
        gh_sb_puts(&sb, ", \"action\": ");
        json_string(&sb, w->rule->action);
        gh_sb_puts(&sb, "}");
    }
    gh_sb_puts(&sb, "], \"by_severity\": {");
    for (size_t i = 0; i < sc->polypharmacy.nseverities; i++) {
        if (i)
            gh_sb_puts(&sb, ", ");
        json_string(&sb, sc->polypharmacy.severities[i]);
        gh_sb_printf(&sb, ": %zu", sc->polypharmacy.severity_counts[i]);
    }
    gh_sb_puts(&sb, "}},\n");

    /* Preventive care */
    gh_sb_printf(&sb, "  \"preventive_care\": {\"early_screenings\": %zu, "
                      "\"summary\": ", sc->preventive.early_screenings);
    json_string(&sb, sc->preventive.summary);
    gh_sb_puts(&sb, ", \"timeline\": [");
    for (size_t i = 0; i < sc->preventive.ntimeline; i++) {
        const gh_screening *s = &sc->preventive.timeline[i];
        if (i)
            gh_sb_puts(&sb, ", ");
        gh_sb_puts(&sb, "{\"test\": ");
        json_string(&sb, s->test);
        gh_sb_printf(&sb, ", \"start_age\": %d, \"frequency\": ",
                     s->start_age);
        json_string(&sb, s->frequency);
        gh_sb_puts(&sb, ", \"reason\": ");
        json_string(&sb, s->reason);
        gh_sb_puts(&sb, ", \"priority\": ");
        json_string(&sb, s->priority);
        gh_sb_puts(&sb, ", \"genetic_basis\": ");
        if (s->genetic_basis)
            json_string(&sb, s->genetic_basis);
        else
            gh_sb_puts(&sb, "null");
        gh_sb_puts(&sb, "}");
    }
    gh_sb_puts(&sb, "]}\n}\n");

    fwrite(sb.data, 1, sb.len, stdout);
    gh_sb_free(&sb);
}

static void emit_text(const gh_analysis *a, const options *o)
{
    printf("\n");
    printf("======================================================================\n");
    printf("GENETIC HEALTH ANALYSIS%s%s\n",
           o->name ? " - " : "", o->name ? o->name : "");
    printf("======================================================================\n\n");
    printf("  Variants in genome:   %zu\n", a->total_snps);
    printf("  Database SNPs found:  %zu of %zu\n", a->analyzed_snps, GH_SNP_COUNT);
    printf("  High impact:          %zu\n", a->high_impact);
    printf("  Moderate impact:      %zu\n", a->moderate_impact);
    printf("  Low impact:           %zu\n", a->low_impact);
    printf("  Drug-gene findings:   %zu\n\n", a->ndrug_findings);

    if (a->high_impact) {
        printf("High-impact findings\n");
        printf("--------------------\n");
        for (size_t i = 0; i < a->nfindings; i++) {
            const gh_finding *f = &a->findings[i];
            if (f->magnitude < 3)
                continue;
            printf("  [%d] %-12s %-10s %-4s  %s\n",
                   f->magnitude, f->rsid, f->gene, f->genotype, f->description);
        }
        printf("\n");
    }

    if (a->ndrug_findings) {
        printf("Drug-gene interactions (ClinPGx)\n");
        printf("--------------------------------\n");
        for (size_t i = 0; i < a->ndrug_findings; i++) {
            const gh_drug_finding *d = &a->drug_findings[i];
            printf("  %-12s %-10s %-4s  level %-3s  %s\n",
                   d->rsid, d->gene, d->genotype, d->level, d->drugs);
        }
        printf("\n");
    }
}

int main(int argc, char **argv)
{
    options o;
    if (!parse_args(&o, argc, argv)) {
        usage(stderr, argv[0]);
        return 2;
    }

    gh_arena *arena = gh_arena_new(1 << 20);
    if (!arena) {
        fprintf(stderr, "error: cannot allocate arena\n");
        return 1;
    }

    const char *genome_path = o.genome ? o.genome
                                       : path_join(arena, o.data_dir, "genome.txt");

    /* The positional lookup is optional: without it, WGS-derived files keep
     * their chr_pos identifiers and simply match fewer database SNPs. */
    gh_map pos_to_rsid;
    bool have_lookup = gh_load_rsid_positions(
        &pos_to_rsid, arena, path_join(arena, o.data_dir, "rsid_positions_grch37.json"));

    gh_genome genome;
    if (!gh_genome_load(&genome, arena, genome_path,
                        have_lookup ? &pos_to_rsid : NULL)) {
        fprintf(stderr, "error: cannot read genome file: %s\n", genome_path);
        gh_arena_free(arena);
        return 1;
    }
    if (!o.quiet) {
        fprintf(stderr, ">>> Loaded %zu variants from %s\n", genome.total, genome_path);
        if (genome.skipped || genome.no_calls)
            fprintf(stderr, "    %zu no-calls, %zu malformed genotypes skipped\n",
                    genome.no_calls, genome.skipped);
        if (genome.resolved)
            fprintf(stderr, "    %zu positional IDs resolved to rsIDs\n", genome.resolved);
    }

    gh_pharmgkb pgx;
    bool have_pgx = gh_pharmgkb_load(
        &pgx, arena,
        path_join(arena, o.data_dir, "clinical_annotations.tsv"),
        path_join(arena, o.data_dir, "clinical_ann_alleles.tsv"));
    if (!o.quiet) {
        if (have_pgx)
            fprintf(stderr, ">>> Loaded %zu ClinPGx drug-gene interactions\n",
                    gh_map_len(&pgx.by_rsid));
        else
            fprintf(stderr, ">>> ClinPGx files not found, skipping drug interactions\n");
    }

    gh_analysis analysis;
    gh_analyze(&analysis, arena, &genome, have_pgx ? &pgx : NULL);


    gh_quality quality;
    gh_quality_compute(&quality, arena, &genome, genome_path);

    gh_scorer_results scorers;
    scorers.stars = gh_calloc(arena, GH_STAR_GENE_COUNT, sizeof(*scorers.stars));
    gh_call_apoe(&scorers.apoe, arena, &genome);
    gh_predict_blood_type(&scorers.blood, arena, &genome);
    gh_estimate_mt_haplogroup(&scorers.mt, arena, &genome);
    gh_call_star_alleles(scorers.stars, arena, &genome);
    scorers.traits = gh_calloc(arena, GH_TRAIT_COUNT, sizeof(*scorers.traits));
    gh_predict_traits(scorers.traits, arena, &genome);
    gh_estimate_ancestry(&scorers.ancestry, arena, &genome);
    scorers.prs = gh_calloc(arena, GH_PRS_MODEL_COUNT, sizeof(*scorers.prs));
    gh_calculate_prs(scorers.prs, arena, &genome, &scorers.ancestry);

    bool have_clinvar = gh_clinvar_analyze(
        &scorers.clinvar, arena, &genome,
        path_join(arena, o.data_dir, "clinvar_alleles.tsv"));
    if (!o.quiet) {
        if (have_clinvar)
            fprintf(stderr,
                    ">>> Scanned %zu ClinVar entries, %zu at genotyped positions\n",
                    scorers.clinvar.total_clinvar, scorers.clinvar.matched);
        else
            fprintf(stderr, ">>> ClinVar file not found, skipping disease risk\n");
    }
    gh_flag_acmg(&scorers.acmg, arena, &scorers.clinvar);
    gh_organize_carriers(&scorers.carriers, arena, &scorers.clinvar);
    scorers.epistasis = gh_evaluate_epistasis(arena, &analysis,
                                              &scorers.nepistasis);
    gh_profile_histamine(&scorers.histamine, arena, &genome);
    gh_profile_alcohol(&scorers.alcohol, arena, &genome);
    gh_profile_pain(&scorers.pain, arena, &genome);
    gh_profile_thyroid(&scorers.thyroid, arena, &genome);
    gh_profile_hormone(&scorers.hormone, arena, &genome);
    gh_profile_eye(&scorers.eye, arena, &genome);

    gh_recommendation_inputs rec_in = {
        &analysis, &scorers.clinvar, scorers.prs, GH_PRS_MODEL_COUNT,
        scorers.epistasis, scorers.nepistasis,
        scorers.stars, GH_STAR_GENE_COUNT, &scorers.acmg,
    };
    gh_generate_recommendations(&scorers.recs, arena, &rec_in);

    gh_insight_inputs ins_in = {
        &analysis, &scorers.apoe, scorers.stars, GH_STAR_GENE_COUNT,
        &scorers.clinvar,
    };
    gh_generate_insights(&scorers.insights, arena, &ins_in);

    gh_profile_sleep(&scorers.sleep, arena, &genome, &analysis);
    gh_profile_nutrigenomics(&scorers.nutrition, arena, &analysis);
    gh_profile_mental_health(&scorers.mental, arena, &genome, &analysis,
                             scorers.stars, GH_STAR_GENE_COUNT);
    gh_profile_longevity(&scorers.longevity, arena, &genome, &analysis,
                         &scorers.apoe, scorers.prs, GH_PRS_MODEL_COUNT);

    gh_generate_drug_dosing(&scorers.dosing, arena, scorers.stars,
                            GH_STAR_GENE_COUNT, &analysis);
    gh_assess_polypharmacy(&scorers.polypharmacy, arena, scorers.stars,
                           GH_STAR_GENE_COUNT, &analysis);
    gh_generate_preventive_timeline(&scorers.preventive, arena, scorers.prs,
                                    GH_PRS_MODEL_COUNT, &scorers.apoe,
                                    &scorers.acmg, scorers.stars,
                                    GH_STAR_GENE_COUNT);

    if (o.html) {
        gh_strbuf html;
        gh_sb_init(&html);
        gh_report_opts ropts = {o.name, o.date};
        gh_report_scorers rsc = {
            .apoe = &scorers.apoe, .blood = &scorers.blood, .mt = &scorers.mt,
            .stars = scorers.stars, .nstars = GH_STAR_GENE_COUNT,
            .traits = scorers.traits, .ntraits = GH_TRAIT_COUNT,
            .ancestry = &scorers.ancestry,
            .prs = scorers.prs, .nprs = GH_PRS_MODEL_COUNT,
            .clinvar = &scorers.clinvar, .acmg = &scorers.acmg,
            .carriers = &scorers.carriers,
            .epistasis = scorers.epistasis, .nepistasis = scorers.nepistasis,
            .recs = &scorers.recs, .insights = &scorers.insights,
            .sleep = &scorers.sleep, .nutrition = &scorers.nutrition,
            .mental = &scorers.mental, .longevity = &scorers.longevity,
            .dosing = &scorers.dosing, .polypharmacy = &scorers.polypharmacy,
            .preventive = &scorers.preventive,
            .pain = &scorers.pain, .histamine = &scorers.histamine,
            .thyroid = &scorers.thyroid, .hormone = &scorers.hormone,
            .eye = &scorers.eye, .alcohol = &scorers.alcohol,
        };
        gh_report_render(&html, arena, &analysis, &quality, &rsc, &ropts);

        if (strcmp(o.html, "-") == 0) {
            fwrite(html.data, 1, html.len, stdout);
        } else {
            FILE *out = fopen(o.html, "wb");
            if (!out) {
                fprintf(stderr, "error: cannot write %s\n", o.html);
                gh_sb_free(&html);
                gh_arena_free(arena);
                return 1;
            }
            fwrite(html.data, 1, html.len, out);
            fclose(out);
            if (!o.quiet)
                fprintf(stderr, ">>> Report written to %s (%zu bytes)\n",
                        o.html, html.len);
        }
        gh_sb_free(&html);
    } else if (o.json) {
        emit_json(&analysis, &quality, &scorers, &o);
    } else {
        emit_text(&analysis, &o);
    }

    if (!o.quiet)
        fprintf(stderr, ">>> Arena high-water mark: %.1f MB\n",
                (double)gh_arena_used(arena) / (1024.0 * 1024.0));

    gh_arena_free(arena);
    return 0;
}

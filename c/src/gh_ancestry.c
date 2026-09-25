#include "gh_ancestry.h"

#include <math.h>
#include <string.h>

static size_t count_allele(const char *genotype, char allele)
{
    size_t n = 0;
    for (const char *p = genotype; *p; p++)
        if (*p == allele)
            n++;
    return n;
}

/* Clamp a frequency away from 0 and 1 so the log terms stay finite. */
static double clamp_freq(double f)
{
    if (f < 0.001)
        return 0.001;
    if (f > 0.999)
        return 0.999;
    return f;
}

/* Softmax in place, shifting by the maximum first for numerical stability --
 * the same guard the Python uses. */
static void softmax(double *values, size_t n)
{
    double max = values[0];
    for (size_t i = 1; i < n; i++)
        if (values[i] > max)
            max = values[i];

    double total = 0.0;
    for (size_t i = 0; i < n; i++) {
        values[i] = exp(values[i] - max);
        total += values[i];
    }
    for (size_t i = 0; i < n; i++)
        values[i] /= total;
}

/* Index of the largest value. Ties go to the earliest, matching Python's
 * max() over an insertion-ordered dict. */
static size_t argmax(const double *values, size_t n)
{
    size_t best = 0;
    for (size_t i = 1; i < n; i++)
        if (values[i] > values[best])
            best = i;
    return best;
}

/* Sub-population estimate within the winning superpopulation. */
static void estimate_sub(gh_sub_ancestry *out, const gh_genome *g,
                         const char *population)
{
    memset(out, 0, sizeof(*out));

    const gh_sub_population *set = NULL;
    for (size_t i = 0; i < GH_SUB_POPULATION_COUNT; i++)
        if (strcmp(GH_SUB_POPULATIONS[i].population, population) == 0) {
            set = &GH_SUB_POPULATIONS[i];
            break;
        }
    if (!set)
        return;   /* this superpopulation has no sub-population markers */

    double log_lls[GH_MAX_SUB_POPULATIONS] = {0};
    size_t markers_used = 0;

    for (size_t m = 0; m < set->nmarkers; m++) {
        const gh_sub_marker *marker = &set->markers[m];
        const char *genotype = gh_genome_genotype(g, marker->rsid);
        if (!genotype)
            continue;

        size_t n = count_allele(genotype, marker->allele);
        markers_used++;

        for (size_t s = 0; s < set->nsubs; s++) {
            double freq = clamp_freq(marker->frequencies[s]);
            log_lls[s] += (double)n * log(freq)
                        + (2.0 - (double)n) * log(1.0 - freq);
        }
    }

    /* Python returns None below three markers rather than guessing. */
    if (markers_used < 3)
        return;

    softmax(log_lls, set->nsubs);

    out->present = true;
    out->markers_used = markers_used;
    out->nsubs = set->nsubs;
    for (size_t s = 0; s < set->nsubs; s++) {
        out->labels[s] = set->labels[s];
        out->proportions[s] = log_lls[s];
    }
    out->top_label = set->labels[argmax(log_lls, set->nsubs)];
    out->confidence = markers_used >= 6 ? "moderate" : "low";
}

void gh_estimate_ancestry(gh_ancestry_result *out, gh_arena *arena,
                          const gh_genome *g)
{
    memset(out, 0, sizeof(*out));
    out->top_index = -1;

    double log_lls[GH_MAX_POPULATIONS] = {0};
    double total_informativeness = 0.0;

    out->details = gh_calloc(arena, GH_AIM_COUNT, sizeof(*out->details));

    for (size_t i = 0; i < GH_AIM_COUNT; i++) {
        const gh_aim *aim = &GH_AIMS[i];
        const char *genotype = gh_genome_genotype(g, aim->rsid);
        if (!genotype)
            continue;

        size_t n = count_allele(genotype, aim->allele);
        out->markers_found++;

        /* Informativeness: the spread between the most and least common
         * frequency for this marker, an Fst-like discriminating power. */
        double lo = aim->frequencies[0], hi = aim->frequencies[0];
        for (size_t p = 1; p < GH_POPULATION_COUNT; p++) {
            if (aim->frequencies[p] < lo)
                lo = aim->frequencies[p];
            if (aim->frequencies[p] > hi)
                hi = aim->frequencies[p];
        }
        total_informativeness += hi - lo;

        for (size_t p = 0; p < GH_POPULATION_COUNT; p++) {
            double freq = clamp_freq(aim->frequencies[p]);
            log_lls[p] += (double)n * log(freq)
                        + (2.0 - (double)n) * log(1.0 - freq);
        }

        gh_aim_detail *d = &out->details[out->ndetails++];
        d->rsid = aim->rsid;
        d->gene = aim->gene;
        d->description = aim->description;
        d->genotype = genotype;
        d->allele_count = (int)n;
    }

    if (out->markers_found == 0) {
        for (size_t p = 0; p < GH_POPULATION_COUNT; p++)
            out->proportions[p] = 1.0 / (double)GH_POPULATION_COUNT;
        out->confidence = "none";
        out->top_ancestry = "Unknown";
        out->ndetails = 0;
        return;
    }

    for (size_t p = 0; p < GH_POPULATION_COUNT; p++)
        out->proportions[p] = log_lls[p];
    softmax(out->proportions, GH_POPULATION_COUNT);

    /* Confidence weighs how many markers were found against how
     * discriminating they were on average. */
    double avg_info = total_informativeness / (double)out->markers_found;
    double effective = (double)out->markers_found * avg_info;

    if (out->markers_found >= 40 && effective >= 15.0)
        out->confidence = "high";
    else if (out->markers_found >= 20 && effective >= 8.0)
        out->confidence = "moderate";
    else
        out->confidence = "low";

    size_t top = argmax(out->proportions, GH_POPULATION_COUNT);
    out->top_index = (int)top;
    out->top_ancestry = GH_POPULATION_LABELS[top];

    estimate_sub(&out->sub, g, GH_POPULATIONS[top]);
}

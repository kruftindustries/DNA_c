/* Every rsID the analysis reads, across all modules.
 *
 * The WGS pipeline calls variants only at positions the analysis will look
 * at, and needs the list to build them. It used to come from a Python
 * function that hand-copied each module's markers and forgot some (sleep,
 * longevity, mental health, sub-population ancestry). Taking it from the C
 * tables gives one source of truth; tests/test_qt_parity.py checks the old
 * list is a subset.
 */
#ifndef GH_RSIDS_H
#define GH_RSIDS_H

#include <stddef.h>

#include "gh_mem.h"

/* Distinct rsIDs, sorted as strings. Returns the count; *out is arena-held. */
size_t gh_all_rsids(gh_arena *arena, const char ***out);

#endif /* GH_RSIDS_H */

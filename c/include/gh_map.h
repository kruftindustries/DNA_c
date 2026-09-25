/* String-keyed hash map with open addressing.
 *
 * Used for rsID -> genotype, annotation ID -> record, and similar lookups.
 * Keys are copied into the caller-supplied arena; values are borrowed
 * pointers the map does not own.
 */
#ifndef GH_MAP_H
#define GH_MAP_H

#include <stdbool.h>
#include <stddef.h>

#include "gh_mem.h"

typedef struct {
    const char *key;   /* NULL marks an empty slot */
    void *value;
    unsigned hash;
} gh_map_entry;

typedef struct {
    gh_map_entry *slots;
    size_t cap;        /* always a power of two */
    size_t len;
    gh_arena *arena;
} gh_map;

/* Initialise a map. `hint` is an expected element count (0 for a default). */
void gh_map_init(gh_map *m, gh_arena *arena, size_t hint);

/* Insert or overwrite. The key is copied; the value pointer is stored as-is.
 * Returns true if the key was newly added, false if it replaced a value. */
bool gh_map_put(gh_map *m, const char *key, void *value);

/* Look up a key. Returns NULL when absent. */
void *gh_map_get(const gh_map *m, const char *key);

/* Look up a key given an explicit length (key need not be NUL-terminated). */
void *gh_map_getn(const gh_map *m, const char *key, size_t len);

bool gh_map_has(const gh_map *m, const char *key);

size_t gh_map_len(const gh_map *m);

/* Iterate. Start with *iter = 0; returns false when exhausted.
 * `key` and/or `value` may be NULL if not wanted. */
bool gh_map_next(const gh_map *m, size_t *iter, const char **key, void **value);

/* FNV-1a, exposed so callers can pre-hash when it helps. */
unsigned gh_hash_bytes(const char *s, size_t n);

#endif /* GH_MAP_H */

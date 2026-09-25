#include "gh_map.h"

#include <string.h>

#define GH_MAP_MIN_CAP 16
/* Grow at 70% load: linear probing degrades sharply past that. */
#define GH_MAP_LOAD_NUM 7
#define GH_MAP_LOAD_DEN 10

unsigned gh_hash_bytes(const char *s, size_t n)
{
    unsigned h = 2166136261u;
    for (size_t i = 0; i < n; i++) {
        h ^= (unsigned char)s[i];
        h *= 16777619u;
    }
    /* Ensure a non-zero hash so 0 can stay a sentinel if ever needed. */
    return h ? h : 1u;
}

static size_t round_pow2(size_t n)
{
    size_t c = GH_MAP_MIN_CAP;
    while (c < n)
        c *= 2;
    return c;
}

void gh_map_init(gh_map *m, gh_arena *arena, size_t hint)
{
    /* Size for the hint at the target load factor so a map built from a known
     * row count never rehashes. */
    size_t want = hint ? (hint * GH_MAP_LOAD_DEN) / GH_MAP_LOAD_NUM + 1 : GH_MAP_MIN_CAP;
    m->cap = round_pow2(want);
    m->len = 0;
    m->arena = arena;
    m->slots = gh_calloc(arena, m->cap, sizeof(*m->slots));
}

/* Find the slot for `key`: either the one holding it, or the first free slot
 * where it belongs. */
static gh_map_entry *find_slot(gh_map_entry *slots, size_t cap,
                               const char *key, size_t len, unsigned hash)
{
    size_t mask = cap - 1;
    size_t i = hash & mask;
    for (;;) {
        gh_map_entry *e = &slots[i];
        if (!e->key)
            return e;
        if (e->hash == hash && strncmp(e->key, key, len) == 0 && e->key[len] == '\0')
            return e;
        i = (i + 1) & mask;
    }
}

static void grow(gh_map *m)
{
    size_t new_cap = m->cap * 2;
    gh_map_entry *new_slots = gh_calloc(m->arena, new_cap, sizeof(*new_slots));

    for (size_t i = 0; i < m->cap; i++) {
        gh_map_entry *e = &m->slots[i];
        if (!e->key)
            continue;
        gh_map_entry *dst = find_slot(new_slots, new_cap, e->key,
                                      strlen(e->key), e->hash);
        *dst = *e;
    }
    /* The old slot array is arena memory; it is reclaimed with the arena. */
    m->slots = new_slots;
    m->cap = new_cap;
}

bool gh_map_put(gh_map *m, const char *key, void *value)
{
    if ((m->len + 1) * GH_MAP_LOAD_DEN > m->cap * GH_MAP_LOAD_NUM)
        grow(m);

    size_t len = strlen(key);
    unsigned hash = gh_hash_bytes(key, len);
    gh_map_entry *e = find_slot(m->slots, m->cap, key, len, hash);

    if (e->key) {
        e->value = value;
        return false;
    }
    e->key = gh_strndup(m->arena, key, len);
    e->hash = hash;
    e->value = value;
    m->len++;
    return true;
}

void *gh_map_getn(const gh_map *m, const char *key, size_t len)
{
    if (!m->slots)
        return NULL;
    unsigned hash = gh_hash_bytes(key, len);
    gh_map_entry *e = find_slot(m->slots, m->cap, key, len, hash);
    return e->key ? e->value : NULL;
}

void *gh_map_get(const gh_map *m, const char *key)
{
    return key ? gh_map_getn(m, key, strlen(key)) : NULL;
}

bool gh_map_has(const gh_map *m, const char *key)
{
    return gh_map_get(m, key) != NULL;
}

size_t gh_map_len(const gh_map *m)
{
    return m->len;
}

bool gh_map_next(const gh_map *m, size_t *iter, const char **key, void **value)
{
    for (size_t i = *iter; i < m->cap; i++) {
        if (m->slots[i].key) {
            if (key)
                *key = m->slots[i].key;
            if (value)
                *value = m->slots[i].value;
            *iter = i + 1;
            return true;
        }
    }
    *iter = m->cap;
    return false;
}

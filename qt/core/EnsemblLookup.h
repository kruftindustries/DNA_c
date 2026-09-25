// rsID -> GRCh37 position, from the Ensembl REST API, cached in
// data/rsid_positions_grch37.json. Mirrors wgs_pipeline.step_rsid_lookup,
// including the fix for merged rsIDs: Ensembl answers under a variant's
// current ID and the requested one comes back in `synonyms`.
#pragma once

#include <QJsonObject>
#include <QNetworkAccessManager>
#include <QSet>
#include <QString>
#include <QStringList>

#include "Reporter.h"

namespace EnsemblLookup {

extern const char *const kEndpoint;        // GRCh37
extern const char *const kEndpointGRCh38;

// The first mapping of a variation record that lies on a real chromosome
// (patch scaffolds sit alongside the assembly), with "chrom" added as the
// bare label. Empty when there is none.
QJsonObject primaryMapping(const QJsonObject &info);

// POST up to 200 ids to a variation endpoint; `results` is keyed by the
// variant's current name, with requested ids that were merged in `synonyms`.
bool postVariationBatch(QNetworkAccessManager &nam, const char *endpoint,
                        const QStringList &ids, QJsonObject *results, QString *error);

// GET one id from a variation endpoint. The GRCh38 host in particular goes
// through spells of answering every batch POST with a 500 while single GETs
// still (slowly) work.
bool getVariation(QNetworkAccessManager &nam, const char *endpoint, const QString &id,
                  QJsonObject *result, QString *error);

// Loads the lookup. A path that does not exist yet is first seeded from
// the copy compiled into the executable (data/rsid_positions_grch37.json at
// build time), which is what a packaged build starts from.
QJsonObject loadLookup(const QString &path);
bool seedLookup(const QString &path);   // true if the file exists afterwards
bool saveLookup(const QString &path, const QJsonObject &lookup);

// Fold one batch response into `lookup`, keyed by the ids that were asked
// for. Pure; returns how many requested ids gained a position.
int mergeBatch(const QJsonObject &results, const QSet<QString> &requested,
               QJsonObject &lookup);

// Query for whatever in `rsids` is not yet in the cache and save it.
// `unresolved` receives the ids that still have no position. False when
// Ensembl could not be asked about some of them (with `error` set); the
// ids that were answered are saved regardless.
bool update(const QString &lookupPath, const QStringList &rsids,
            const Reporter &reporter, QStringList *unresolved, QString *error = nullptr);

} // namespace EnsemblLookup

// HTTP download to a file with progress, following redirects and sending a
// real User-Agent (ClinPGx answers 403 to library defaults). Synchronous:
// spins a local event loop, so it is called from the console tool's thread
// or from a worker thread, never from the GUI thread.
#pragma once

#include <QString>
#include <QUrl>

#include "Reporter.h"

namespace Downloader {

extern const char *const kUserAgent;

// Returns false and sets *error on failure. Honours reporter.cancelled().
bool download(const QUrl &url, const QString &dest, const Reporter &reporter,
              QString *error);

} // namespace Downloader

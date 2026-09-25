#include "Reporter.h"

#include <QPair>

#include <cstdio>
#include <memory>

Reporter stdoutReporter()
{
    Reporter r;
    r.log = [](const QString &s) {
        std::fputs(qPrintable(s), stdout);
        std::fputc('\n', stdout);
        std::fflush(stdout);
    };
    // Stage changes are chatty during downloads; print on a name change or
    // every ten percent.
    auto last = std::make_shared<QPair<QString, int>>(QString(), -100);
    r.stage = [last](const QString &name, int pct) {
        if (name != last->first || pct >= last->second + 10 || pct == 100) {
            std::fprintf(stdout, "[%3d%%] %s\n", pct, qPrintable(name));
            std::fflush(stdout);
            *last = {name, pct};
        }
    };
    return r;
}

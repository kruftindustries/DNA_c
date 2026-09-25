// How long-running core work talks back: log lines, a stage with a percent,
// and a cancellation check. Plain callables rather than signals, so the same
// code runs under the console tool (print) and the GUI (queued signals from
// a worker thread) without QObject thread-affinity rules getting involved.
#pragma once

#include <QString>
#include <functional>

struct Reporter {
    std::function<void(const QString &)> log = [](const QString &) {};
    std::function<void(const QString &, int)> stage = [](const QString &, int) {};
    std::function<bool()> cancelled = [] { return false; };

    void info(const QString &s) const { log(s); }
};

// A reporter that prints to stdout, for the console tool and tests.
Reporter stdoutReporter();

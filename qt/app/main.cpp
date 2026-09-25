#include <QApplication>

#include "MainWindow.h"
#include "SystemTheme.h"

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    QCoreApplication::setOrganizationName("genetic-health");
    QCoreApplication::setApplicationName("genetic-health-qt");
    SystemTheme theme(app);   // follow the desktop's light/dark preference
    MainWindow w;
    w.show();
    return app.exec();
}

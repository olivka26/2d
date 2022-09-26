#include <qapplication.h>
#include <QMainWindow>
#include <QMenuBar>
#include "glwidget.hpp"

int main(int argc, char *argv[]){
    QApplication app(argc, argv);
    MainWindow w;
    app.setActiveWindow(&w);
    if(w.parse_command_line(argc, argv)){
        qWarning("Wrong input arguments!");
        return -1;
    }
    w.show();
    return app.exec();
}




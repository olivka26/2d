QT += opengl
win32:LIBS += -lOpenGL32
TEMPLATE = app
TARGET = 2d
INCLUDEPATH += .

DEFINES += QT_DEPRECATED_WARNINGS

HEADERS += glwidget.hpp chebyshev.hpp help.hpp
SOURCES += glwidget.cpp main.cpp chebyshev.cpp help.cpp
FORMS += \
        mainwindow.ui


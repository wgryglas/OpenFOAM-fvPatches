TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    ../iteracje.cpp \
    ../funkcje.cpp \
    ../main.cpp \
    ../f_struktury.cpp \
    ../obliczenia.cpp \
    ../src/solvers.cpp \
    ../src/interpolation.cpp \
    ../src/dataanalysis.cpp \
    ../src/ap.cpp \
    ../src/alglibmisc.cpp \
    ../src/diffequations.cpp \
    ../src/linalg.cpp \
    ../src/statistics.cpp \
    ../src/specialfunctions.cpp \
    ../src/optimization.cpp \
    ../src/integration.cpp \
    ../src/fasttransforms.cpp \
    ../src/alglibinternal.cpp \
    ../postproc.cpp

HEADERS += \
    ../iteracje.h \
    ../struktury.h \
    ../f_struktury.h \
    ../funkcje.h \
    ../obliczenia.h \
    ../src/solvers.h \
    ../src/ap.h \
    ../src/stdafx.h \
    ../src/statistics.h \
    ../src/specialfunctions.h \
    ../src/optimization.h \
    ../src/linalg.h \
    ../src/interpolation.h \
    ../src/integration.h \
    ../src/fasttransforms.h \
    ../src/diffequations.h \
    ../src/dataanalysis.h \
    ../src/alglibmisc.h \
    ../src/alglibinternal.h \
    ../postproc.h


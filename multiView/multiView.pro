# A multi-view application, with constrained camera displacements.

# Four viewers are created, each displaying the same scene. The camera displacements are constrained
# for three of the viewers to create the classical top, front, side views. The last viewer is a
# classical 3D viewer.

# Note that the four viewers share the same OpenGL context.

TEMPLATE = app
TARGET   = multiView

HEADERS  = \
    camerapathplayer.h \
    controlpoint.h \
    curvepoint.h \
    curve.h \
    mainwindow.h \
    mesh.h \
    meshreader.h \
    plane.h \
    standardcamera.h \
    viewer.h \
    Triangle.h \
    Vec3D.h \
    viewerfibula.h
SOURCES  = main.cpp \
    camerapathplayer.cpp \
    controlpoint.cpp \
    curvepoint.cpp \
    curve.cpp \
    mainwindow.cpp \
    mesh.cpp \
    plane.cpp \
    standardcamera.cpp \
    viewer.cpp \
    viewerfibula.cpp

include( ../baseInclude.pri )

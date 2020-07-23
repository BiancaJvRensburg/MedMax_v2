#ifndef PANDA_H
#define PANDA_H

#include <QGLViewer/qglviewer.h>
#include "simplemesh.h"
#include "QGLViewer/manipulatedFrame.h"
#include "Manipulator/simplemanipulator.h"

using namespace qglviewer;

class Panda
{
public:
    Panda();

    void draw();

    void setLocation(const Vec &v);
    void setOrientation(const Quaternion &q);
    void setToPlane(const Vec &v, const Quaternion &q);
    void rotate(Vec axis, double angle){ f.rotate(Quaternion(axis, angle)); }
    void checkPanda();
    void setFrames();

private:
    void drawFrame(Frame &frame);
    void openOFF(QString f, SimpleMesh &m);


    SimpleMesh effector, navex, marker, wrj;
    Frame f, fsw, ftcp;     // f : panda_link8
    Vec tcpPanda[4];
};

#endif // PANDA_H

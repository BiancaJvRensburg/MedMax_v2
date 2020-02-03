#ifndef CONTROLPOINT_H
#define CONTROLPOINT_H

#include "camerapathplayer.h"
#include <QGLViewer/qglviewer.h>
#include <QGLViewer/manipulatedFrame.h>

using namespace qglviewer;

class ControlPoint : public QObject
{

    Q_OBJECT

public:
    ControlPoint(Vec p);
    ControlPoint(double x, double y, double z);

    const Vec& getPoint(){ return p; }
    double getX(){ return p.x; }
    double getY(){ return p.y; }
    double getZ(){ return p.z; }
    void setPoint(double x, double y, double z){
        p.x = x;
        p.y = y;
        p.z = z;
    }
    void setPoint(Vec p){ this->p = p; }
    Frame* getFrame(){  return mf; }

    void moveControlPoint(Vec newPos);

    void initialise();
    virtual void draw();

public Q_SLOTS:
    virtual void cntrlMoved();

Q_SIGNALS:
    void cntrlPointTranslated();

protected:
    Frame* mf;
    Vec p;
};

#endif // CONTROLPOINT_H

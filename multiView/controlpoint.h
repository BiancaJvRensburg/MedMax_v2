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
    ControlPoint();
    ControlPoint(const Vec& p);
    ControlPoint(double x, double y, double z);

    Vec& getPoint(){ return p; }
    double& getX(){ return p.x; }
    double& getY(){ return p.y; }
    double& getZ(){ return p.z; }
    void setPosition(double& x, double& y, double& z){
        p.x = x;
        p.y = y;
        p.z = z;
    }
    void setPosition(Vec& p){ this->p = p; }
    ManipulatedFrame& getFrame(){ return mf; }
    const Frame* getReferenceFrame(){ return  mf.referenceFrame(); }
    const Quaternion& getOrientation();

    void initialise();
    virtual void draw();

    void toggleSwitchFrames(){ isSwitchFrames = !isSwitchFrames; }

public Q_SLOTS:
    virtual void cntrlMoved();

Q_SIGNALS:
    void cntrlPointTranslated();

protected:
    ManipulatedFrame mf;
    Vec p;
    bool isSwitchFrames;
};

#endif // CONTROLPOINT_H

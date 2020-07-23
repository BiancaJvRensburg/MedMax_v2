#ifndef PANDAMANIPULATOR_H
#define PANDAMANIPULATOR_H

#include <vector>
using std::vector;
#include <map>
using std::map;
using std::pair;

#include <QGLViewer/qglviewer.h>
#include "Tools/GLUtilityMethods.h"
class GLVWidget;
#include <QtGui>
#include <QGLViewer/mouseGrabber.h>

class PandaManipulator: public QObject , public qglviewer::MouseGrabber
{
    Q_OBJECT

public:
    PandaManipulator();

    void checkIfGrabsMouse(int x, int y,const qglviewer::Camera* const cam);
    int getState( ){ return this->state; }

    void activate(){ this->setState(1); }
    void deactivate(){ this->setState(0); }
    void switchStates(){ if(getState()==0) activate(); else deactivate(); }
    void setState(const bool &b){ if(b) activate(); else deactivate(); }
    void setID(unsigned int i){ id = i; }

    void setOrigin( qglviewer::Vec const & p ){ Origin = p; }
    void setRepX( qglviewer::Vec const & p ){ RepX = p; }
    void setRepY( qglviewer::Vec const & p ){ RepY = p; }
    void setRepZ( qglviewer::Vec const & p ){ RepZ = p; }

    void setRotationActivated(bool b){ isRotationActivated = b; }

    const qglviewer::Vec& getPosition(){ return Origin; }
    void getOrientation(qglviewer::Vec &x, qglviewer::Vec &y, qglviewer::Vec &z){
        x = RepX;
        y = RepY;
        z = RepZ;
    }

    void setDisplayScale(float ds){ display_scale = ds; }
    int getModification(){ return this->mode_modification; }
    void resetScales(){ Xscale = Yscale = Zscale = 1.f; }

    void draw();
    void clear();

    void manipulatedCallback();
    void fakeMouseDoubleClickEvent( QMouseEvent* const );
    void mousePressEvent( QMouseEvent* const event  , qglviewer::Camera* const cam );
    void mouseReleaseEvent( QMouseEvent* const , qglviewer::Camera* const  );
    void mouseMoveEvent(QMouseEvent* const event, qglviewer::Camera* const cam);

Q_SIGNALS:
    void mouseReleased();
    void moved(qglviewer::Vec);

private:
    void setState( int e ){ this->state = e; }

    int state;
    //  0:  desactive.
    //  1:  il vient d'etre cree car l'user a clique sur une zone selectionnee quand aucun mode de selection n'etait actif

    int mode_modification;
    //  1:  tx
    //  2:  ty
    //  3:  tz
    //  4:  rx
    //  5:  ry
    //  6:  rz
    //  7:  sx
    //  8:  sy
    //  9:  sz

    int mode_grabbing;
    // 0 : grabs it !
    // 1 : axis + rotations + scales

    bool mouse_released;


    qglviewer::Vec Origin;              // (1)
    qglviewer::Vec RepX, RepY, RepZ;    // (2)

    double uTeta , vTeta ;
    qglviewer::Vec PrevPos ;

    double display_scale;
    double Xscale , Yscale , Zscale;     // (3)

    vector< pair< int , qglviewer::Vec > > coordinates;
    map< int,int > idpoints;


    float m_xx_default , m_yy_default;

    unsigned int id;        // same as the plane's id
    bool isRotationActivated = true;
};

#endif // PANDAMANIPULATOR_H




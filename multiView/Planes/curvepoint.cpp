#include "curvepoint.h"

// A sub class of Control point. This is a point placed on the curve
CurvePoint::CurvePoint(const unsigned int &id) : ControlPoint ()
{
    this->mf = ManipulatedFrame();
    this->id = id;
    //connect(&mf, &ManipulatedFrame::manipulated, this, &ControlPoint::cntrlMoved);
}

void CurvePoint::draw(){
    if(isSwitchFrames){
        glPushMatrix();
        glMultMatrixd(mf.matrix());
    }

    /*if(mf.grabsMouse()) glColor3f(0, 1, 1);
    else glColor3f(0.6f, 0, 0.4f);*/

    glPointSize(10.0);
    glBegin(GL_POINTS);
        glVertex3d(0, 0, 0);
    glEnd();

    glPointSize(1.0);
    glColor3f(1.0,1.0,1.0);

    if(isSwitchFrames) glPopMatrix();
}

// If a point is moved by hand, update it and send the id of the point and its new position
void CurvePoint::cntrlMoved(){
    double x,y,z;

    mf.getPosition(x,y,z);

    Vec position(x,y,z);
    p.x = x;
    p.y = y;
    p.z = z;

    Q_EMIT CurvePoint::curvePointTranslated(id, position);
}

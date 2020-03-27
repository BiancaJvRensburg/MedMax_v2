#ifndef CURVEPOINT_H
#define CURVEPOINT_H

#include "controlpoint.h"

class CurvePoint : public ControlPoint
{
    Q_OBJECT

public:
    CurvePoint(Vec& p);
    //CurvePoint(CurvePoint &cp);
    /*~CurvePoint(){
        disconnect((ManipulatedFrame*)mf, &ManipulatedFrame::manipulated, this, &ControlPoint::cntrlMoved);
    }*/

    void setPosition(Vec& p){ this->p = p; mf.setPosition(getX(), getY(), getZ()); }

    void matchCurvepoint(CurvePoint &c);
    Quaternion getOrientation(){ return mf.orientation(); }

    void draw();

public Q_SLOTS:
    void cntrlMoved();

Q_SIGNALS:
    void curvePointTranslated(Vec offset);
};

#endif // CURVEPOINT_H

#ifndef CURVEPOINT_H
#define CURVEPOINT_H

#include "Curve/controlpoint.h"

class CurvePoint : public ControlPoint
{
    Q_OBJECT

public:
    CurvePoint(const unsigned int &id);

    void setPosition(Vec& p){ mf.setPosition(p.x, p.y, p.z); } //this->p = p;  }
    void setID(unsigned int id){ this->id = id; }
    Vec getPosition(){ return  mf.position(); }

    Quaternion getOrientation(){ return mf.orientation(); }

    void draw();

public Q_SLOTS:
    void cntrlMoved();

Q_SIGNALS:
    void curvePointTranslated(unsigned int pointIndex, Vec offset);

private:
    unsigned int id;

};

#endif // CURVEPOINT_H

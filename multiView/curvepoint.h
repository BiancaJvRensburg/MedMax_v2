#ifndef CURVEPOINT_H
#define CURVEPOINT_H

#include "controlpoint.h"

class CurvePoint : public ControlPoint
{
    Q_OBJECT

public:
    CurvePoint(Vec p);
    ~CurvePoint(){
        disconnect((ManipulatedFrame*)mf, &ManipulatedFrame::manipulated, this, &ControlPoint::cntrlMoved);
    }

    void setPosition(Vec p){ this->p.x = p.x;
                                this->p.y = p.y;
                                this->p.z = p.z; }

    void draw();

public Q_SLOTS:
    void cntrlMoved();

Q_SIGNALS:
    void curvePointTranslated(Vec offset);
};

#endif // CURVEPOINT_H

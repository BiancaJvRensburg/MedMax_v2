#ifndef PLANE_H
#define PLANE_H

#include <QGLViewer/qglviewer.h>
#include <QGLViewer/manipulatedFrame.h>

#include "curvepoint.h"

enum Movable {STATIC, DYNAMIC};

using namespace qglviewer;

class Plane
{
public:
    Plane(double s, Movable status);
    ~Plane(){
        //std::cout << "destroying plane " << this << std::endl;
        // delete cp;
    }

    void setSize(double s){ size = s; }
    void setPosition(Vec pos);
    void setOrientation(Quaternion q){ cp->getFrame().setOrientation(q); }
    Quaternion fromRotatedBasis(Vec x, Vec y, Vec z);

    Vec getPolylineVector(Vec v){ return cp->getFrame().localCoordinatesOf(v); }  // Return the vector v in the coordinates of this plane (could be done w/ another function)

    void rotatePlaneXY(double percentage);   // rotate around the z axis    // NOTE could be useless in the near future
    void rotatePlane(Vec axis, double angle);
    void constrainZRotation(){ cp->getFrame().setConstraint(&constraint); }
    void freeZRotation(){ cp->getFrame().setConstraint(&constraintFree); }
    void draw();

    // Mesh calculations
    bool isIntersection(Vec v0, Vec v1, Vec v2);
    double getSign(Vec v);

    Vec getNormal(){ return normal; }
    //const Frame& getFrame(){ return *cp->getFrame(); }
    Vec getProjection(Vec p);
    const Vec& getPosition(){ return cp->getPoint(); }
    CurvePoint* getCurvePoint(){ return cp; }

    Vec getLocalCoordinates(Vec v) { return cp->getFrame().localCoordinatesOf(v); }    // same as get polyline
    Vec getMeshCoordinatesFromLocal(Vec v){ return cp->getFrame().localInverseCoordinatesOf(v); }
    Vec getLocalVector(Vec v) { return cp->getFrame().localTransformOf(v); }    // same as get polyline
    Vec getMeshVectorFromLocal(Vec v){ return cp->getFrame().localInverseTransformOf(v); }

    Movable status;

private:
    void initBasePlane();
    void rotate(Quaternion q) { cp->getFrame().rotate(q); }

    AxisPlaneConstraint constraint;
    AxisPlaneConstraint constraintFree;
    Vec points[4];
    double size;
    double rotationPercentage;
    Vec normal;
    CurvePoint* cp;
};

#endif // PLANE_H

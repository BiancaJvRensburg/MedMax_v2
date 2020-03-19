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
    Plane(double s, Movable status, Vec& pos, float alpha);

    void toggleIsVisible(){ isVisible = !isVisible; }
    void setAlpha(float alpha){ this->alpha = alpha; }
    float getAlpha(){ return alpha; }

    void setSize(double s){ size = s; }
    void setPosition(Vec pos);
    void setOrientation(Quaternion q){ cp.getFrame().setOrientation(q); }
    Quaternion fromRotatedBasis(Vec x, Vec y, Vec z);
    void setFrameFromBasis(Vec x, Vec y, Vec z);

    Vec getPolylineVector(Vec v){ return cp.getFrame().localCoordinatesOf(v); }  // Return the vector v in the coordinates of this plane (could be done w/ another function)

    void rotatePlaneXY(double percentage);   // rotate around the z axis    // NOTE could be useless in the near future
    void rotatePlane(Vec axis, double angle);
    void setPlaneRotation(Vec axis, double angle);
    void constrainZRotation(){ cp.getFrame().setConstraint(&constraint); }
    void freeZRotation(){ cp.getFrame().setConstraint(&constraintFree); }
    void draw();
    const Frame* getReferenceFrame(){ return cp.getReferenceFrame(); }

    // Mesh calculations
    bool isIntersection(Vec v0, Vec v1, Vec v2);
    double getSign(Vec v);

    Vec getNormal(){ return normal; }
    //const Frame& getFrame(){ return *cp->getFrame(); }
    Vec getProjection(Vec p);
    Vec getLocalProjection(Vec p);      // for vectors already in local coordinates
    Vec& getPosition(){ return cp.getPoint(); }
    CurvePoint& getCurvePoint(){ return cp; }

    Vec getLocalCoordinates(Vec v) { return cp.getFrame().localCoordinatesOf(v); }    // same as get polyline
    Vec getMeshCoordinatesFromLocal(Vec v){ return cp.getFrame().localInverseCoordinatesOf(v); }
    Vec getLocalVector(Vec v) { return cp.getFrame().localTransformOf(v); }    // same as get polyline
    Vec getMeshVectorFromLocal(Vec v){ return cp.getFrame().localInverseTransformOf(v); }

    Frame getFrameCopy();

    void setOrientationFromOtherReference(std::vector<Vec> &frame, unsigned int startIndex, Plane* reference);
    void setRotation(Quaternion q) { cp.getFrame().setRotation(q); }
    void rotate(Quaternion q) { cp.getFrame().rotate(q); }

    bool isIntersectionPlane(Vec &v0, Vec &v1, Vec &v2, Vec &v3);
    void getCorners(Vec &v0, Vec &v1, Vec &v2, Vec &v3);

    void matchPlane(Plane *p);

    Movable status;

private:
    void initBasePlane();
    AxisPlaneConstraint constraint;
    AxisPlaneConstraint constraintFree;
    Vec points[4];
    double size;
    double rotationPercentage;
    Vec normal;
    CurvePoint cp;
    bool isVisible;
    float alpha;
};

#endif // PLANE_H

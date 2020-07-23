#ifndef BOX_H
#define BOX_H

#include <QGLViewer/qglviewer.h>
#include <QGLViewer/frame.h>
#include "Manipulator/simplemanipulator.h"

using namespace qglviewer;

class Box
{
public:
    Box();
    void draw(double offset);
    void init(const Frame *ref);
    void setPosition(const Vec& v){ f.setPosition(v); }
    //void setEndPosition(const Vec &v);
    void setFrameFromBasis(Vec x, Vec y, Vec z);
    void setLength(double l){ dimensions.x = l; }
    void setWidth(double l){ dimensions.y = l; }
    void setHeight(double l){ dimensions.z = l; }
    const double& getLength(){ return dimensions.x; }
    const double& getWidth(){ return dimensions.y; }
    const double& getHeight(){ return dimensions.z; }
    void rotateOnAxis(double angle);
    Vec localTransform(Vec v){ return f.localTransformOf(v); }
    Vec localCoordinates(Vec v){ return f.localCoordinatesOf(v); }
    Vec worldCoordinates(Vec v){ return f.localInverseCoordinatesOf(v); }
    Vec worldTransform(Vec v){ return f.localInverseTransformOf(v); }
    Vec worldTangent(){ return worldTransform(tangent); }
    Vec worldBinormal(){ return worldTransform(binormal); }
    Vec worldNormal(){ return worldTransform(normal); }
    Vec getLocation();
    Vec getEnd();
    Vec getMidPoint();
    Vec getHighPoint();
    Vec getHighEnd();
    const Vec& getNormal(){ return normal; }
    const Vec& getBinormal(){ return binormal; }
    const Vec& getTangent(){ return tangent; }
    void restoreRotation();
    void getOrientation(Vec &x, Vec &y, Vec &z);

private:
    Frame f;
    Vec normal, binormal, tangent;
    Vec dimensions;
    double prevRotation;
    //SimpleManipulator manipulator;
};

#endif // BOX_H

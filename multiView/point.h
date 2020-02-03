#ifndef POINT_H
#define POINT_H

#include <QGLViewer/qglviewer.h>

class Point
{
public:
    Point();
    Point(double x, double y, double z);

    Point operator-(const Point& p){
        Point r;
        r.setX(this->x - p.x);
        r.setY(this->y - p.y);
        r.setZ(this->z - p.z);
        return r;
    }

    Point operator*(const long& n){
        Point r;
        r.setX(this->x * n);
        r.setY(this->y * n);
        r.setZ(this->z * n);
        return r;
    }

    double getX(){ return x;}
    double getY(){ return y;}
    double getZ(){ return z;}

    void setX(double x){this->x = x;}
    void setY(double y){this->y = y;}
    void setZ(double z){this->z = z;}

    qglviewer::Vec toVec(){ return qglviewer::Vec(x, y, z); }

private:
    double x;
    double y;
    double z;
};

#endif // POINT_H

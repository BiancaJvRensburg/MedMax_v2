#ifndef CURVE_H
#define CURVE_H

#include <QGLViewer/qglviewer.h>
#include "controlpoint.h"

using namespace qglviewer;

class Curve : public QObject
{
    Q_OBJECT

public:
    Curve(unsigned int nbCP, std::vector<Vec>& cntrlPoints);

    void generateBSpline(unsigned int& nbU, unsigned int degree);
    void generateCatmull(unsigned int& nbU);

    Vec* getCurve(){ return curve; }
    Vec& getPoint(int index){ return curve[index]; }
    unsigned int& getNbU(){ return nbU; }

    Vec tangent(int index);
    Vec normal(int index);
    Vec binormal(int index);
    void getFrame(int index, Vec& t, Vec& n, Vec& b);

    void draw();
    void drawControl();
    void drawTangent(int index);

    double discreteLength(unsigned int indexS, unsigned int indexE);      // Returns the discrete length between 2 points (Straight line distance)
    double discreteChordLength(unsigned int indexS, unsigned int indexE); // To use for the initial visualisation
    unsigned int indexForLength(unsigned int indexS, double length);   // Returns the end index which will create a segment of a certain length

public Q_SLOTS:
    void reintialiseCurve();

Q_SIGNALS:
    void curveReinitialised();

private:
    std::vector<ControlPoint*> TabControlPoint;
    unsigned int nbControlPoint;
    Vec *curve;
    unsigned int nbU;

    unsigned int degree;
    std::vector<double> knotVector;
    unsigned int knotIndex;

    // BSpline
    std::vector<double> generateUniformKnotVector(unsigned int k);
    Vec deBoor(double u, unsigned int i, unsigned int r);
    Vec* splineDerivative(unsigned int k);
    Vec deBoorDerivative(double u, unsigned int i, unsigned int r, unsigned int k);

    // Catmull rom
    void catmullrom();  // calculate the spline and the first derivative
    void calculateCatmullPoints(Vec& c, Vec& cp, Vec& cpp, double t);

    std::vector<double> generateCatmullKnotVector(double alpha);

    void initConnections();

    // Frenet frame
    Vec* dt;
    Vec* d2t;
};

#endif // CURVE_H

#ifndef CURVE_H
#define CURVE_H

#include <QGLViewer/qglviewer.h>
#include "controlpoint.h"

using namespace qglviewer;

class Curve : public QObject
{
    Q_OBJECT

public:
    Curve();
    void init(unsigned int nbCP, std::vector<Vec>& cntrlPoints);

    void generateCatmull(unsigned int& nbU);

    std::vector<Vec>& getCurve(){ return curve; }
    Vec& getPoint(unsigned int index){ return curve[index]; }
    unsigned int& getNbU(){ return nbU; }

    Vec tangent(unsigned int index);
    Vec normal(unsigned int index);
    Vec binormal(unsigned int index);
    void getFrame(unsigned int index, Vec& t, Vec& n, Vec& b);

    void draw();
    void drawControl();
    void drawFrame(unsigned int index);

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
    std::vector<Vec> curve;
    unsigned int nbU;

    unsigned int degree;
    std::vector<double> knotVector;
    unsigned int knotIndex;

    // Catmull rom
    void catmullrom();  // calculate the spline and the first derivative
    void calculateCatmullPoints(Vec& c, Vec& cp, Vec& cpp, double t);

    void generateCatmullKnotVector(double alpha, std::vector<double>& knotV);

    void initConnections();

    unsigned int getClosestDistance(double target, unsigned int indexS, unsigned int a, unsigned int b);       // get the index which is closest to the target distance

    // Frenet frame
    std::vector<Vec> dt;
    std::vector<Vec> d2t;
};

#endif // CURVE_H

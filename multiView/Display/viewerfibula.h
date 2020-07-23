#ifndef VIEWERFIBULA_H
#define VIEWERFIBULA_H

#include "viewer.h"

class ViewerFibula : public Viewer
{
    Q_OBJECT

public:
    ViewerFibula(QWidget *parent, StandardCamera *camera, int sliderMax, int fibulaOffsetMax);
    void initGhostPlanes(Movable s);
    void updateFibPolyline(const Vec& firstPoint, std::vector<double>& distances);
    void constructCurve();
    void toggleIsPolyline();
    void repositionPlanesOnPolyline();  
    void initSignals();
   // void draw();

public Q_SLOTS:
    void bendPolylineNormals(std::vector<Vec>&, std::vector<double>&);
    void bendPolyline(unsigned int pointIndex, Vec v);
    void constructPolyline(std::vector<double>&, const std::vector<Vec>&);
    void updateDistances(std::vector<double>&);
    void movePlanes(double);
    void updatePlaneOrientations(std::vector<Vec>&);
    void rotatePolylineOnAxisFibula(double);
    void cut();
    void uncut();
    void recieveFromFibulaMesh(std::vector<int>&, std::vector<Vec>&, std::vector<std::vector<int>>&, std::vector<int>&, std::vector<Vec>&, int);
    void tryOffsetAngle();
    void reinitBox(unsigned int, std::vector<double>&);
    void reinitPoly(unsigned int);
    void reprojectToMesh();
    void slidePolyline(int);
    void setPanda();
    void handlePandaManipulated(Vec);

Q_SIGNALS:
    void okToPlacePlanes(const std::vector<Vec>&);
    void sendToManible(std::vector<int>&, std::vector<Vec>, std::vector<std::vector<int>>&, std::vector<int>&, std::vector<Vec>, int);
    void requestNewNorms();
    void requestFakeBend();

private:
    void rotatePolyline();
    void rotatePolyToCurve();
    void positionBoxes();
    Vec getOffsetDistanceToMeshBorder(std::vector<Vec>& projections, Plane &p);
    Vec getGreatestOffsetDistance(Plane &p1, Plane &p2, unsigned int p1I, unsigned int p2I);
    Vec getMinOfMin(Vec &a, Vec &b);
    void bendWithRelativeVector(Plane &p, Vec v);
    void setPlanesInPolyline(std::vector<Vec> &normals);
    void setPlaneOrientations(std::vector<Vec> &normals);
    void setDistances(std::vector<double> &distances);
    void constructSegmentPoints(unsigned int nbU);
    void projectToMesh(const std::vector<double>& distances);
    void matchDistances(const std::vector<double>& distances, std::vector<unsigned int> &segIndexes, std::vector<Vec> &outputPoints, double epsilon, const unsigned int& searchRadius);
    unsigned int getClosestDistance(unsigned int index, const double &targetDistance, std::vector<unsigned int> &segIndexes, std::vector<Vec> &outputPoints, unsigned int searchRadius);
    double getOffsetDistance(double angle);
    double getBoxPlaneAngle(Plane &p);
    void modifyDistances(std::vector<double> &distances){ if(distances.size() > 0) distances[0] = polylineOffset; }

    void positionPanda();

    double prevRotation = 0.;
    std::vector<double> saveDistances;
    double polylineOffset = 0.001;

};

#endif // VIEWERFIBULA_H

#ifndef VIEWERFIBULA_H
#define VIEWERFIBULA_H

#include "viewer.h"

class ViewerFibula : public Viewer
{
    Q_OBJECT

public:
    ViewerFibula(QWidget *parent, StandardCamera *camera, int sliderMax, int fibulaOffsetMax);
    void addGhostPlanes(unsigned int nb);
    void handleCut();
    std::vector<Vec> getPolyline();

public Q_SLOTS:
    void movePlanes(int);
    void planesMoved();
    void movePlaneDistance(double, std::vector<Vec>, std::vector<Vec>);
    void ghostPlanesRecieved(unsigned int, double[], std::vector<Vec>, std::vector<Vec>);
    void middlePlaneMoved(unsigned int, double[], std::vector<Vec>, std::vector<Vec>);

    void initCurve();
    void constructCurve();
    void cutMesh();
    void uncutMesh();

    void noGhostPlanesToRecieve(std::vector<Vec>, std::vector<Vec>, double);

    void recieveAxes(std::vector<Vec>);
    void recieveFromFibulaMesh(const std::vector<int>&, const std::vector<Vec>&, const std::vector<std::vector<int>>&, const std::vector<int>&, const std::vector<Vec>&, int);

Q_SIGNALS:
    void setPlaneSliderValue(int);
    void sendToManible(const std::vector<int>&, std::vector<Vec>, const std::vector<std::vector<int>>&, const std::vector<int>&, std::vector<Vec>, int);
    void requestAxes();

private:
    void findGhostLocations(unsigned int nb, double distance[]); // finds the location of the ghost planes + the right plane
    void initSignals();
    void createPolyline(std::vector<Vec> &polyline);
    void repositionPlanes(std::vector<Vec>& polyline, std::vector<Vec>& axes);
    void setPlaneOrientations();
    void setPlanePositions();
    void resetMandibleInfo(std::vector<Vec>& polyline, std::vector<Vec>& axes);
    void swivelToPolyline(std::vector<Vec>& fibulaPolyline);
    void swivelPlane(Plane *p, const Vec& fibPoly, const Vec& mandPoly);
    void findIndexesFromDistances();

    void findClosestPoint(unsigned int pNb, Vec &a, Vec &b);
    Vec findMinZ(const std::vector<unsigned int> &tIndex, Plane &tempPlane);
    Vec findMaxZ(const std::vector<unsigned int> &tIndex, Plane &tempPlane);
    void approachPlanes(unsigned int pStart);
    double euclideanDistance(Vec &a, Vec &b);

    bool isCutSignal;
    bool isPlanesRecieved;

    int indexOffset;
    int prevOffset;
    int maxOffset;
    std::vector<Vec> mandiblePolyline;      // the last mandible polyline we recieved
    std::vector<Vec> mandibleAxes;          // the last mandible axes we recieved
    std::vector<double> distances;
    bool isNeedToFlip;

    const double securityMargin = 30;       // this is temporary
};

#endif // VIEWERFIBULA_H

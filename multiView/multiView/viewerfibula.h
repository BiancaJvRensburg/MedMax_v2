#ifndef VIEWERFIBULA_H
#define VIEWERFIBULA_H

#include "viewer.h"

class ViewerFibula : public Viewer
{
    Q_OBJECT

public:
    ViewerFibula(QWidget *parent, StandardCamera *camera, int sliderMax, int fibulaOffsetMax);
    void addGhostPlanes(int nb);
    void handleCut();
    std::vector<Vec> getPolyline();


public Q_SLOTS:
    void movePlanes(int);
    void planesMoved();
    void movePlaneDistance(double, std::vector<Vec>, std::vector<Vec>);
    void ghostPlanesRecieved(int, double[], std::vector<Vec>, std::vector<Vec>);
    void middlePlaneMoved(int, double[], std::vector<Vec>, std::vector<Vec>);

    void initCurve();
    void cutMesh();
    void uncutMesh();

    void noGhostPlanesToRecieve();

    void handleMovementStart();
    void handleMovementEnd();

    void recieveTest(std::vector<Vec>);

    void recieveFromFibulaMesh(std::vector<int>, std::vector<Vec>, std::vector<std::vector<int>>, std::vector<int>, std::vector<Vec>, int);

Q_SIGNALS:
    void setPlaneSliderValue(int);
    void sendToManible(std::vector<int>, std::vector<Vec>, std::vector<std::vector<int>>, std::vector<Vec>, std::vector<int>, std::vector<Vec>, int);

private:
    void findGhostLocations(int nb, double distance[]); // finds the location of the ghost planes + the right plane
    void setPlaneOrientations(std::vector<Vec> mandPolyline, std::vector<Vec> axes);
    void reinitialisePlanes(unsigned int nbToInit);      // Reinitialises the position and orientation of the planes
    void initSignals();
    void createPolyline();
    void matchToMandibleFrame(Plane* p1, Plane* p2, Vec a, Vec b, Vec c, Vec x, Vec y, Vec z);

    bool isCutSignal;
    bool isPlanesRecieved;

    int indexOffset;
    int maxOffset;
    std::vector<Vec> mandiblePolyline;      // the last mandible polyline we recieved
    std::vector<Vec> mandibleAxes;          // the last mandible axes we recieved
};

#endif // VIEWERFIBULA_H

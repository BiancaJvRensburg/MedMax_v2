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
    void ghostPlanesRecieved(unsigned int, double[], std::vector<Vec>, std::vector<Vec>);
    void middlePlaneMoved(unsigned int, double[], std::vector<Vec>, std::vector<Vec>);

    void initCurve();
    void cutMesh();
    void uncutMesh();

    void noGhostPlanesToRecieve();

    void handleMovementStart();
    void handleMovementEnd();

    void recieveAxes(std::vector<Vec>);

    void recieveFromFibulaMesh(std::vector<int>, std::vector<Vec>, std::vector<std::vector<int>>, std::vector<int>, std::vector<Vec>, int);

Q_SIGNALS:
    void setPlaneSliderValue(int);
    void sendToManible(std::vector<int>, std::vector<Vec>, std::vector<std::vector<int>>, std::vector<Vec>, std::vector<int>, std::vector<Vec>, int);
    void requestAxes();

private:
    void findGhostLocations(unsigned int nb, double distance[]); // finds the location of the ghost planes + the right plane
    void initSignals();
    void createPolyline();
    void matchToMandibleFrame(Plane* p1, Plane* p2, Vec a, Vec b, Vec c, Vec x, Vec y, Vec z);
    void repositionPlanes(std::vector<Vec> polyline, std::vector<Vec> axes);
    void setPlaneOrientations();
    void setPlanePositions();
    void resetMandibleInfo(std::vector<Vec> polyline, std::vector<Vec> axes);
    void swivelToPolyline();

    bool isCutSignal;
    bool isPlanesRecieved;

    int indexOffset;
    int maxOffset;
    std::vector<Vec> mandiblePolyline;      // the last mandible polyline we recieved
    std::vector<Vec> mandibleAxes;          // the last mandible axes we recieved
};

#endif // VIEWERFIBULA_H

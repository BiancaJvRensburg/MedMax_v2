#ifndef VIEWER_H
#define VIEWER_H

#include <QGLViewer/qglviewer.h>
#include "meshreader.h"
#include "mesh.h"
#include "standardcamera.h"
#include "plane.h"
#include "curve.h"

using namespace qglviewer;

class Viewer : public QGLViewer
{
    Q_OBJECT

public :
    Viewer(QWidget *parent, StandardCamera *camera, int sliderMax);
    void openOFF(QString f);
    Mesh mesh;

public Q_SLOTS:
    void moveLeftPlane(int);
    void moveRightPlane(int);
    void rotateLeftPlane(int);
    void rotateRightPlane(int);
    void updatePlanes();
    virtual void cutMesh();
    virtual void uncutMesh();
    void ghostPlaneMoved();
    void drawMesh();
    void onLeftSliderReleased();
    void onRightSliderReleased();
    void recieveFromFibulaMesh(std::vector<int>, std::vector<Vec>, std::vector<std::vector<int>>, std::vector<Vec>, std::vector<int>, std::vector<Vec>, int);
    void toUpdate();

Q_SIGNALS:
    void leftPosChanged(double, std::vector<Vec>);
    void rightPosChanged(double, std::vector<Vec>);
    void ghostPlanesAdded(int, double[], std::vector<Vec>);
    void ghostPlanesTranslated(int, double[], std::vector<Vec>);
    void okToCut();
    // set the slider to the value
    void setLRSliderValue(int);   // Left rotation
    void setRRSliderValue(int);   // Right rotation
    void setLMSliderValue(int);   // Left movement
    void setRMSliderValue(int);   // Right movement
    void sendFibulaToMesh(std::vector<Vec>, std::vector<std::vector<int>>, std::vector<int>, std::vector<Vec>, int);

    void noGhostPlanesToSend();     // tells the fibula not to wait for ghost planes before cutting
    void preparingToCut();          // tells the fibula to reset its planes

protected:
    void draw();
    std::vector<Vec> updatePolyline();   // returns the new angles between the polyline and the planes
    void drawPolyline();
    void init();
    QString helpString() const;
    void updateCamera(const Vec3Df & center, float radius);

    virtual void initSignals();
    virtual void initCurve();
    void initPlanes(Movable status);
    virtual void addGhostPlanes(int nb);
    void initGhostPlanes();
    Quaternion getNewOrientation(int index);
    Quaternion updateOrientation(int index);

    double angle(Vec a, Vec b);
    double segmentLength(const Vec a, const Vec b);
    std::vector<Vec> getPolylinePlaneAngles();
    int partition(int sorted[], int start, int end);
    void quicksort(int sorted[], int start, int end);

    ManipulatedFrame* viewerFrame;

    Plane *leftPlane;
    Plane *rightPlane;
    Curve *curve;
    std::vector<Plane> ghostPlanes;
    std::vector<int> ghostLocation;
    std::vector<Vec> polyline;
    int nbGhostPlanes;
    int currentNbGhostPlanes;
    bool isGhostPlanes;

    int curveIndexL;
    int curveIndexR;
    long* nbU;
    int sliderMax;
    bool isDrawMesh;

    const double constraint = 25;

};

#endif // VIEWER_H

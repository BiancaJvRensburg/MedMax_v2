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
    void readJSON(const QJsonArray &json);
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
    void recieveFromFibulaMesh(const std::vector<int>&, std::vector<Vec>, const std::vector<std::vector<int>>&, const std::vector<int>&, std::vector<Vec>, const int);
    void toUpdate();
    void getAxes();
    void toggleIsDrawPlane();
    void setAlpha(int);

Q_SIGNALS:
    void leftPosChanged(double, std::vector<Vec>, std::vector<Vec>);
    void rightPosChanged(double, std::vector<Vec>, std::vector<Vec>);
    void ghostPlanesAdded(unsigned int, double[], std::vector<Vec>, std::vector<Vec>);
    void ghostPlanesTranslated(unsigned int, double[], std::vector<Vec>, std::vector<Vec>);
    void okToCut();

    // set the slider to the value
    void setLRSliderValue(int);   // Left rotation
    void setRRSliderValue(int);   // Right rotation
    void setLMSliderValue(int);   // Left movement
    void setRMSliderValue(int);   // Right movement
    void sendFibulaToMesh(std::vector<Vec>, const std::vector<std::vector<int>>&, const std::vector<int>&, std::vector<Vec>, int);

    void noGhostPlanesToSend(std::vector<Vec>, std::vector<Vec>, double);     // tells the fibula not to wait for ghost planes before cutting
    void preparingToCut();          // tells the fibula to reset its planes

    void ghostPlaneMovementStart();      // tells the fibula to "uncut" the mesh while we move the planes

    void sendAxes(std::vector<Vec>);

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
    virtual void addGhostPlanes(unsigned int nb);

    double angle(Vec a, Vec b);
    double segmentLength(const Vec a, const Vec b);

    Quaternion getNewOrientation(unsigned int index);
    Vec getCustomProjection(Vec a, Vec normal);        // project a onto a plane defined by the normal
    void repositionPlane(Plane* p, unsigned int index);

    virtual void constructCurve();

    ManipulatedFrame* viewerFrame;

    Plane *leftPlane;
    Plane *rightPlane;
    Curve *curve;
    std::vector<Plane*> ghostPlanes;
    std::vector<unsigned int> ghostLocation;
    // std::vector<Vec> polyline;  // just a list of vertex coordinates
    int nbGhostPlanes;
    int currentNbGhostPlanes;
    bool isGhostPlanes;
    bool isGhostActive;

    unsigned int curveIndexL;
    unsigned int curveIndexR;
    unsigned int nbU;
    int sliderMax;
    bool isDrawMesh;

    const double constraint = 25;

    std::vector<Vec> control;
    bool isCurve;

    Vec convertToPlane(Plane *base, Plane *p, Vec axis);        // get the z axis of p in relation to base

private:
    void initGhostPlanes();
    Quaternion updateOrientation(unsigned int index);
    void matchPlaneToFrenet(Plane* p, unsigned int index);
    void handlePlaneMoveStart();
    void handlePlaneMoveEnd();
    void updateMeshPolyline(std::vector<Vec> &polyline);
    void rotateFrame(Frame& f, Vec axis, double angle);
    void balanceGhostPlanes();

    std::vector<Vec> getReferenceAxes();        // get all the z axes in terms of their directors
    std::vector<Vec> getPolylinePlaneAngles(std::vector<Vec> &polyline);      // returns the polyline in the coordinates of each plane, one for each side of the plane
    int partition(int sorted[], int start, int end);
    void quicksort(int sorted[], int start, int end);

    void rotatePlane(Plane* p, int position);
    void movePlane(Plane *p, bool isLeft, unsigned int index);
    void addFrameChangeToAxes(std::vector<Vec> &axes, Plane *base, Plane *p);
    void addInverseFrameChangeToAxes(std::vector<Vec> &axes, Plane *base, Plane *p);

    bool isSpaceForGhosts();
};

#endif // VIEWER_H

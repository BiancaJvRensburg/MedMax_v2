#ifndef VIEWER_H
#define VIEWER_H

#include <QGLViewer/qglviewer.h>
#include "Tools/standardcamera.h"
#include "Planes/plane.h"
#include "Curve/curve.h"
#include "Polyline/polyline.h"
#include "Mesh/mesh.h"
#include "Tools/savedstate.h"
#include "Mesh/panda.h"
#include "Manipulator/pandamanipulator.h"

using namespace qglviewer;

class Viewer : public QGLViewer
{
    Q_OBJECT

public :
    Viewer(QWidget *parent, StandardCamera *camera, int sliderMax);
    void openOFF(QString f);
    void readJSON(const QJsonArray &json);

public Q_SLOTS:
    void toUpdate();
    void bendPolylineManually(unsigned int pointIndex, Vec v, unsigned int s);
    void cutMesh();
    void uncutMesh();
    void placePlanes(const std::vector<Vec>&);
    void moveLeftPlane(int);
    void moveRightPlane(int);
    void rotatePolylineOnAxis(int);
    void toggleDrawMesh();
    void toggleDrawPlane();
    void setPlaneAlpha(int);
    void setMeshAlpha(int);
    void setBoxAlpha(int);
    virtual void cut();
    virtual void uncut();
    void toggleWireframe();
    void recieveFromFibulaMesh(std::vector<int>&, std::vector<Vec>, std::vector<std::vector<int>>&, std::vector<int>&, std::vector<Vec>, int);
    void sendNewNorms();
    void setBoxToManipulator(unsigned int, Vec, int s);
    void setBoxToCornerManipulator(unsigned int, Vec);
    void toggleEditPlaneMode(unsigned int id, bool b);
    void toggleEditBoxMode(unsigned int id, bool b);
    void toggleEditFirstCorner(unsigned int id, bool b);
    void toggleEditEndCorner(unsigned int id, bool b);
    void toggleAllPlanes(bool b);
    void toggleAllBoxes(bool b);
    void toggleAllFirstCorners(bool b);
    void toggleAllEndCorners(bool b);
    void toggleDrawPolyline();
    void toggleDrawBoxes();
    void toggleDrawCurve();
    void saveCurrentState(Modification m);
    void manipulatorReleasedPlane();
    void manipulatorReleasedBox();
    void keyPressEvent(QKeyEvent *e);
    void fakeBend();

Q_SIGNALS:
    void polylineUpdate(const std::vector<Vec>&);
    void polylineBent(std::vector<Vec>&, std::vector<double>&);
    void constructPoly(std::vector<double>&, const std::vector<Vec>&);
    void toUpdateDistances(std::vector<double>&);
    void toUpdatePlaneOrientations(std::vector<Vec>&);
    void toRotatePolylineOnAxis(double);
    void planeMoved(double);
    void sendFibulaToMesh(std::vector<int>&, std::vector<Vec>, std::vector<std::vector<int>>&, std::vector<int>&, std::vector<Vec>, int);
    void cutFibula();
    void uncutFibula();
    void toReinitBox(unsigned int, std::vector<double>&);
    void toReinitPoly(unsigned int);
    void enableFragmentEditing();
    void disableFragmentEditing();
    void editPlane(unsigned int);
    void editBoxCentre(unsigned int);
    void editBoxStart(unsigned int);
    void editBoxEnd(unsigned int);

protected:
    void draw();
    void drawWithNames();
    void postSelection(const QPoint &point);
    void init();
    virtual void initSignals();
    virtual void constructCurve();
    virtual void initGhostPlanes(Movable s);
    virtual void bendPolyline(unsigned int pointIndex, Vec v);
    void initPolyPlanes(Movable s);
    void initCurvePlanes(Movable s);
    void updateCamera(const Vec& center, float radius);
    void updatePolyline(const std::vector<Vec> &newPoints);
    void deleteGhostPlanes();
    void repositionPlane(Plane* p, unsigned int index);
    virtual void repositionPlanesOnPolyline();
    void matchPlaneToFrenet(Plane *p, unsigned int index);
    virtual void toggleIsPolyline();
    unsigned int partition(std::vector<unsigned int>& sorted, unsigned int start, unsigned int end);
    void quicksort(std::vector<unsigned int>& sorted, int start, int end);
    void findGhostLocations(unsigned int nbGhostPlanes, std::vector<unsigned int>& ghostLocations);
    void deconstructPolyline();
    void getPlaneBoxOrientations(std::vector<Vec>& norms);
    void simpleBend(const unsigned int &pointIndex, Vec v, std::vector<Vec>& planeNormals, std::vector<Vec>& planeBinormals);
    void lowerPoints(double size, Vec localDirection);
    void changePlaneDisplaySize(double width, double height);
    Vec projectBoxToPlane(Plane &p, Plane &endP, double& distShift);
    double euclideanDistance(const Vec &a, const Vec &b);
    bool isViolatesContraint();
    void restoreLastState();


    double angle(Vec a, Vec b);
    double segmentLength(const Vec a, const Vec b);

    ManipulatedFrame* viewerFrame;
    class Polyline poly;
    std::vector<Plane*> ghostPlanes;
    Plane* leftPlane;
    Plane* rightPlane;
    Plane* endLeft;
    Plane* endRight;

    std::vector<Vec> control;
    Curve curve;
    Mesh mesh;

    unsigned int curveIndexL;
    unsigned int curveIndexR;
    unsigned int nbU;
    bool isCurve;
    const double constraint = 25.;
    int sliderMax;
    bool isCut;
    bool isPoly = false;
    bool isDrawMesh, isDrawCurve;
    bool isFibula;

    int polyRotation;

    std::vector<Vec> segmentPoints;

    std::deque< SavedState > Q;     // cntrl z

    // temporary for testing
   // std::vector<Vec> testPoints;
    Vec camCentre;
   // Vec projPoint;


    Panda panda;
    PandaManipulator pandaManipulator;

private:
    void constructPolyline(const std::vector<Vec>& polyPoints);
    void movePlane(Plane *p, unsigned int curveIndex);
    void setPlaneOrientation(Plane& p, std::vector<Vec>& norms, std::vector<Vec>& binorms);
    void setPlaneOrientations(std::vector<Vec>& norms, std::vector<Vec>& binorms);
    Plane& getPlaneFromID(unsigned int id);
    Plane& getOppositePlaneFromID(unsigned int id);
    double maxDouble(double a, double b);
    SavedState saveState(Modification m);
    void resetState(SavedState s);
    void resetUndoQueue();

};

#endif // VIEWER_H

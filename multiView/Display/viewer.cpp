#include "viewer.h"
#include "Mesh/meshreader.h"
#include <QGLViewer/manipulatedFrame.h>

static unsigned int max_operation_saved = 10;

Viewer::Viewer(QWidget *parent, StandardCamera *cam, int sliderMax) : QGLViewer(parent) {
    Camera *c = camera();       // switch the cameras
    setCamera(c);
    isCurve = false;
    this->sliderMax = sliderMax;
    this->isCut = false;
    this->isDrawMesh = false;
    this->isDrawCurve = false;
    this->isPoly = false;
    this->isFibula = false;
}

void Viewer::draw() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    Vec base = viewerFrame->localInverseCoordinatesOf(camCentre);
    camera()->setSceneCenter(base);

    glPushMatrix();
    glMultMatrixd(viewerFrame->matrix());

     mesh.draw();

     if(isCurve){

         if(isCut && !isFibula){
             for(unsigned int i=0; i<ghostPlanes.size()+1; i++){
                 mesh.drawFragment(i);
             }
         }

         glColor4f(0., 1., 0., leftPlane->getAlpha());
         leftPlane->draw();
         glColor4f(1., 0, 0., leftPlane->getAlpha());
         rightPlane->draw();

         for(unsigned int i=0; i<ghostPlanes.size(); i++){
             glColor4f(0., 1., 1., ghostPlanes[i]->getAlpha());
             ghostPlanes[i]->draw();
         }

         if(isDrawCurve) curve.draw();
     }

     if(isDrawMesh) mesh.drawCutMand();

    if(isPoly) poly.draw();

    if(isFibula){ panda.draw(); pandaManipulator.draw();}

    glPopMatrix();
}

void Viewer::drawWithNames(){
    if(isCurve){
        glPushName(0);
        glColor4f(0., 1., 0., leftPlane->getAlpha());
        leftPlane->draw();
        glPopName();

        glPushName(1);
        glColor4f(1., 0, 0., leftPlane->getAlpha());
        rightPlane->draw();
        glPopName();


        for(unsigned int i=0; i<ghostPlanes.size(); i++){
            glPushName(i+2);
            glColor4f(0., 1., 1., ghostPlanes[i]->getAlpha());
            ghostPlanes[i]->draw();
            glPopName();
        }

        for(unsigned int i=0; i<ghostPlanes.size()+1; i++){
            glPushName(ghostPlanes.size()+2+i);
            mesh.drawFragment(i);
            glPopName();
        }

    }
}

void Viewer::postSelection(const QPoint &point) {
  Vec orig, dir;
  camera()->convertClickToLine(point, orig, dir);

  bool found;
  Vec selectedPoint = camera()->pointUnderPixel(point, found);
  selectedPoint -= 0.01f * dir;

  int name = selectedName();

  //if(name == 0) leftPlane->toggleEditMode(true);
  //else if(name == 1) rightPlane->toggleEditMode(true);
  if(name < ghostPlanes.size() + 2) Q_EMIT editPlane(name-2); //ghostPlanes[name-2]->toggleEditMode(true);
  else if (name != -1) Q_EMIT editBoxCentre(name-ghostPlanes.size()-1);// poly.toggleBoxManipulators(name-ghostPlanes.size()-1, true);
}

void Viewer::toUpdate(){
    update();
}

void Viewer::init() {
  setMouseTracking(true);
  restoreStateFromFile();
  setBackgroundColor(QColor("gray"));

  viewerFrame = new ManipulatedFrame();
  setManipulatedFrame(viewerFrame);
  setAxisIsDrawn(false);
  poly.init(viewerFrame, 1);

  // Camera without mesh
  Vec centre(0,0,0);
  float radius(15.);
  camCentre = centre;
  updateCamera(camCentre, radius);

  glEnable(GL_LIGHTING);
  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
  glLineWidth (1.0f);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);

  initSignals();
}

// Initialise the planes on the curve
void Viewer::initCurvePlanes(Movable s){
    curveIndexR = nbU - 1;
    curveIndexL = 0;
    float size = 30.0;

    leftPlane = new Plane(static_cast<double>(size), s, 0.5, 0);
    rightPlane = new Plane(static_cast<double>(size), s, 0.5, 1);

    mesh.addPlane(leftPlane);
    mesh.addPlane(rightPlane);

    repositionPlane(leftPlane, curveIndexL);
    repositionPlane(rightPlane, curveIndexR);

    endLeft = new Plane(static_cast<double>(size), Movable::STATIC, 0.5, 0);
    endRight = new Plane(static_cast<double>(size), Movable::STATIC, 0.5, 2);
    repositionPlane(endLeft, 0);
    repositionPlane(endRight, nbU-1);
}

// Position the plane on the curve (before the mesh is cut)
void Viewer::repositionPlane(Plane *p, unsigned int index){
    p->setPosition(curve.getPoint(index));
    matchPlaneToFrenet(p, index);
}

// Orientate the plane on the curve (before the mesh is cut)
void Viewer::matchPlaneToFrenet(Plane *p, unsigned int index){
    Vec x, y, z;
    curve.getFrame(index,z,x,y);
    p->setFrameFromBasis(x,y,z);
}

// Initialise the position and rotations of each plane on the polyline and set their IDs
void Viewer::initPolyPlanes(Movable s){
    leftPlane->setID(1);
    rightPlane->setID(poly.getNbPoints()-2);
    endRight->setID(poly.getNbPoints()-1);

    leftPlane->setPosition(viewerFrame->localCoordinatesOf(poly.getMeshPoint(leftPlane->getID())));
    leftPlane->setFrameFromBasis(Vec(0,0,1), Vec(0,-1,0), Vec(1,0,0));      // TODO this could be done with an existing function?

    rightPlane->setPosition(viewerFrame->localCoordinatesOf(poly.getMeshPoint(rightPlane->getID())));
    rightPlane->setFrameFromBasis(Vec(0,0,1), Vec(0,-1,0), Vec(1,0,0));

    initGhostPlanes(s);
}

void Viewer::toggleIsPolyline(){
    Vec direction = Vec(1,1,0);
    double displaySize = 30.;

    if(!leftPlane->getIsPoly()){
        direction = -direction;
        displaySize = 15.;
    }

    endLeft->toggleIsPoly();
    endRight->toggleIsPoly();

    leftPlane->toggleIsPoly();
    rightPlane->toggleIsPoly();
    for(unsigned int i=0; i<ghostPlanes.size(); i++) ghostPlanes[i]->toggleIsPoly();

    lowerPoints(poly.getBoxHeight(0)/2., direction);       // all the boxes are the same size

    repositionPlanesOnPolyline();
    changePlaneDisplaySize(displaySize, displaySize);
}

void Viewer::changePlaneDisplaySize(double width, double height){
    for(unsigned int i=0; i<ghostPlanes.size(); i++) ghostPlanes[i]->setDisplayDimensions(width, height);
    leftPlane->setDisplayDimensions(width, height);
    rightPlane->setDisplayDimensions(width, height);
}

void Viewer::lowerPoints(double size, Vec localDirection){
    poly.lowerPoint(endLeft->getID(), endLeft->getMeshVectorFromLocal(localDirection)*size);
    poly.lowerPoint(endRight->getID(), endRight->getMeshVectorFromLocal(localDirection)*size);

    poly.lowerPoint(leftPlane->getID(), leftPlane->getMeshVectorFromLocal(localDirection)*size);
    poly.lowerPoint(rightPlane->getID(), rightPlane->getMeshVectorFromLocal(localDirection)*size);
    for(unsigned int i=0; i<ghostPlanes.size(); i++) poly.lowerPoint(i+2, ghostPlanes[i]->getMeshVectorFromLocal(localDirection)*size);
}

void Viewer::repositionPlanesOnPolyline(){
    leftPlane->setPosition(poly.getMeshPoint(leftPlane->getID()));
    for(unsigned int i=0; i<ghostPlanes.size(); i++) ghostPlanes[i]->setPosition(poly.getMeshPoint(ghostPlanes[i]->getID()));
    rightPlane->setPosition(poly.getMeshPoint(rightPlane->getID()));
    endLeft->setPosition(poly.getMeshPoint(endLeft->getID()));
    endRight->setPosition(poly.getMeshPoint(endRight->getID()));
}

void Viewer::initGhostPlanes(Movable s){
    double size = leftPlane->getSize();     // match everythign to the size of the left plane

    for(unsigned int i=2; i<static_cast<unsigned int>(poly.getNbPoints()-2); i++){
        Plane *p = new Plane(size, s, .5f, i);
        p->setPosition(poly.getPoint(i));
        p->setFrameFromBasis(Vec(0,0,1), Vec(0,-1,0), Vec(1,0,0));
        ghostPlanes.push_back(p);
    }

    //connect(&(leftPlane->getManipulator()), &SimpleManipulator::moved, this, &Viewer::bendPolylineManually);
    //connect(&(rightPlane->getManipulator()), &SimpleManipulator::moved, this, &Viewer::bendPolylineManually);
    for(unsigned int i=0; i<ghostPlanes.size(); i++) {
        connect(&(ghostPlanes[i]->getManipulator()), &SimpleManipulator::moved, this, &Viewer::bendPolylineManually);
        connect(&(ghostPlanes[i]->getManipulator()), &SimpleManipulator::mouseReleased, this, &Viewer::manipulatorReleasedPlane);
    }
    // connnect the ghost planes
}

void Viewer::updateCamera(const Vec& center, float radius){
    camCentre = center;
    camera()->setSceneCenter(camCentre);
    camera()->setSceneRadius(static_cast<double>(radius*1.05f));
    camera()->setZClippingCoefficient(static_cast<double>(10.));
    camera()->showEntireScene();
}

void Viewer::constructPolyline(const std::vector<Vec> &polyPoints){
    poly.reinit(polyPoints.size());     // re-initialise the polyline
    std::vector<double> dists;
    poly.getDistances(dists);       // the distances should all be 1 here
    Q_EMIT constructPoly(dists, polyPoints);        // The fibula polyline has to be initialised before we can start bending the polylines
}

void Viewer::deconstructPolyline(){
    poly.reinit(1);
    deleteGhostPlanes();
    toggleIsPolyline();
    repositionPlane(leftPlane, curveIndexL);
    repositionPlane(rightPlane, curveIndexR);
}

/*
 * Physically place planes on the polyline
 * Bend the polyline so it corresponds to the correct point
*/
void Viewer::placePlanes(const std::vector<Vec> &polyPoints){
    initPolyPlanes(Movable::DYNAMIC);       // Create the planes
    std::vector<Vec> norms, binorms;
    for(unsigned int i=0; i<poly.getNbPoints()-1; i++) simpleBend(i, polyPoints[i], norms, binorms);        // Bend the polyline point by point WITHOUT updating rest
    bendPolyline(poly.getNbPoints()-1, polyPoints[poly.getNbPoints()-1]);       // update the fibula and set the planes for all points on the last point
    toggleIsPolyline();     // Change the polyline from going through the centre of the planes to going through the corner

    poly.resetBoxes();      // Set the boxes to the polyline

    std::vector<double> distances;
    poly.getDistances(distances);
    Q_EMIT toUpdateDistances(distances);        // the distances are no longer the same because the polyline has been lowered, so update it in the fibula

    sendNewNorms();
}

double Viewer::segmentLength(const Vec a, const Vec b){
    return sqrt( pow((b.x - a.x), 2) + pow((b.y - a.y), 2) + pow((b.z - a.z), 2));
}

void Viewer::updatePolyline(const std::vector<Vec> &newPoints){
    poly.updatePoints(newPoints);
    poly.resetBoxes();
}


/*
 * Modify the polyline point pointIndex to v
 * Activated when a plane is moved by hand in the mandible
*/
void Viewer::bendPolyline(unsigned int pointIndex, Vec v){
    std::vector<Vec> planeNormals;
    std::vector<Vec> planeBinormals;

    simpleBend(pointIndex, v, planeNormals, planeBinormals);

    // Update the planes
    repositionPlanesOnPolyline();
    setPlaneOrientations(planeNormals, planeBinormals);

    // Get the mandible polyine planes in relation to the boxes so this info can be sent to the fibula
    std::vector<Vec> relativeNorms;
    getPlaneBoxOrientations(relativeNorms);

    // Get the new distances between each point in the mandible
    std::vector<double> distances;
    poly.getDistances(distances);

    Q_EMIT polylineBent(relativeNorms, distances);      // Send the new normals and distances to the fibula TODO only send over the info for the corresponding point
}

void Viewer::fakeBend(){
    std::vector<Vec> relativeNorms;
    getPlaneBoxOrientations(relativeNorms);

    // Get the new distances between each point in the mandible
    std::vector<double> distances;
    poly.getDistances(distances);

    Q_EMIT polylineBent(relativeNorms, distances);
}

void Viewer::bendPolylineManually(unsigned int pointIndex, Vec v, unsigned int s){
    bool isOriginallyCut = isCut;
    if(isOriginallyCut) uncut();
    Q_EMIT toReinitPoly(poly.getNbPoints());
    bendPolyline(pointIndex, v);
    if(isOriginallyCut) cut();
}

// Bend the polyline without updating everything else
void Viewer::simpleBend(const unsigned int &pointIndex, Vec v,std::vector<Vec> &planeNormals, std::vector<Vec> &planeBinormals){
    poly.bend(pointIndex, v, planeNormals, planeBinormals);     // Bend the polyline and get the planeNormals and planeBinormals which will be used to set the planes' orientations
}

/* Plane orientations */

void Viewer::setPlaneOrientation(Plane &p, std::vector<Vec> &norms, std::vector<Vec> &binorms){
    p.setFrameFromBasis(norms[p.getID()], binorms[p.getID()], cross(norms[p.getID()], binorms[p.getID()]));
}

void Viewer::setPlaneOrientations(std::vector<Vec> &norms, std::vector<Vec> &binorms){
    setPlaneOrientation(*leftPlane, norms, binorms);
    setPlaneOrientation(*rightPlane, norms, binorms);
    setPlaneOrientation(*endLeft, norms, binorms);
    setPlaneOrientation(*endRight, norms, binorms);
    for(unsigned int i=0; i<ghostPlanes.size(); i++) setPlaneOrientation(*ghostPlanes[i], norms, binorms);
}

// Get the rotation of the planes in terms of the boxes
void Viewer::getPlaneBoxOrientations(std::vector<Vec> &relativeNorms){
    relativeNorms.clear();
    std::vector<Vec> norms;
    poly.getRelativePlane(*leftPlane, norms);
    relativeNorms.push_back(norms[2]);
    relativeNorms.push_back(norms[3]);

    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        std::vector<Vec> norms;
        poly.getRelativePlane(*ghostPlanes[i], norms);
        for(unsigned int j=0; j<norms.size(); j++){
            relativeNorms.push_back(norms[j]);
        }
    }

    poly.getRelativePlane(*rightPlane, norms);
    relativeNorms.push_back(norms[0]);
    relativeNorms.push_back(norms[1]);
}

// Create the curve from the coordinates given in the JSON
void Viewer::constructCurve(){
    curve.init(control.size(), control);
    nbU = 100;
    curve.generateCatmull(nbU);
    isCurve = true;
    initCurvePlanes(Movable::DYNAMIC);
}

void Viewer::deleteGhostPlanes(){
    for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];
    ghostPlanes.clear();
}

void Viewer::readJSON(const QJsonArray &controlArray){
    control.clear();

    for(int i=0; i<controlArray.size(); i++){
        QJsonArray singleControl = controlArray[i].toArray();
        control.push_back(Vec(singleControl[0].toDouble(), singleControl[1].toDouble(), singleControl[2].toDouble()));
    }

    constructCurve();
}

// Open a .OFF file
void Viewer::openOFF(QString filename) {
    std::vector<Vec3Df> &vertices = mesh.getVertices();
    std::vector<Triangle> &triangles = mesh.getTriangles();

    FileIO::openOFF(filename.toStdString(), vertices, triangles);

    mesh.init();

    // Set the camera
    Vec3Df center;
    float radius;
    mesh.computeBB(center, radius);
    updateCamera(Vec(center),radius);

    update();
}

void Viewer::cutMesh(){
    bool isNumberRecieved;
    unsigned int nbGhostPlanes;
    int nbPieces = QInputDialog::getInt(this, "Cut mesh", "Number of pieces", 0, 1, 10, 1, &isNumberRecieved, Qt::WindowFlags());
    if(isNumberRecieved) nbGhostPlanes = static_cast<unsigned int>(nbPieces)-1;
    else return;

    isCut = true;
    isPoly = true;

    // Construct the polyline : First and last planes are immovable and are at the ends of the meshes
    std::vector<Vec> polylinePoints;
    //polylinePoints.push_back(curve.getPoint(0));        // the start of the curve
    unsigned int boxSize = 10;
    if(boxSize>curveIndexL) polylinePoints.push_back(curve.getPoint(0));
    else polylinePoints.push_back(curve.getPoint(curveIndexL-boxSize));
    polylinePoints.push_back(curve.getPoint(curveIndexL));  // the left plane

    std::vector<unsigned int> ghostLocations;
    findGhostLocations(nbGhostPlanes, ghostLocations);      // find the curve indexes on which the ghost planes must be placed
    for(unsigned int i=0; i<ghostLocations.size(); i++) polylinePoints.push_back(curve.getPoint(ghostLocations[i]));        // get the world location of these indexes

    if(boxSize+curveIndexL >= nbU) polylinePoints.push_back(curve.getPoint(nbU-1));
    else polylinePoints.push_back(curve.getPoint(curveIndexR));      // the right plane
    polylinePoints.push_back(curve.getPoint(curveIndexR+boxSize));
    //polylinePoints.push_back(curve.getPoint(nbU-1));        // the end of the curve

    constructPolyline(polylinePoints);

    cut();

    for(unsigned int i=1; i<poly.getNbPoints()-2; i++) {
        connect(poly.getBoxManipulator(i), &SimpleManipulator::moved, this, &Viewer::setBoxToManipulator);
        connect(poly.getBoxManipulator(i), &SimpleManipulator::mouseReleased, this, &Viewer::manipulatorReleasedBox);
    }
    for(unsigned int i=1; i<(poly.getNbPoints()-2)*2; i++) {
        connect(poly.getCornerManipulator(i), &SimpleManipulator::moved, this, &Viewer::setBoxToCornerManipulator);
        connect(poly.getCornerManipulator(i), &SimpleManipulator::mouseReleased, this, &Viewer::manipulatorReleasedPlane);
    }

    Q_EMIT enableFragmentEditing();
    saveCurrentState(Modification::PLANE);
    update();
}

void Viewer::uncutMesh(){
    uncut();

    isCut = false;
    isPoly = false;
    deconstructPolyline();

    Q_EMIT disableFragmentEditing();
    //update();

    resetUndoQueue();
}

void Viewer::moveLeftPlane(int position){
    if(isPoly) return;       // disable the ability to move the planes once the mesh is cut

    double percentage = static_cast<double>(position) / static_cast<double>(sliderMax);
    unsigned int index = static_cast<unsigned int>(percentage * static_cast<double>(nbU) );

    if(curve.indexForLength(curveIndexR, -constraint) > index){  // Only move if we're going backwards or we haven't met the other plane
        curveIndexL = index;
        if(curveIndexL >= nbU) curveIndexL = nbU-1;     // shouldn't ever happen
    }
    else if( curveIndexL == curve.indexForLength(curveIndexR, -constraint) ){
        return;       // already in the correct position
    }
    else curveIndexL = curve.indexForLength(curveIndexR, -constraint);     // get the new position

    movePlane(leftPlane, curveIndexL);
}

void Viewer::moveRightPlane(int position){
    if(isPoly) return;

    double percentage = static_cast<double>(position) / static_cast<double>(sliderMax);
    unsigned int index = nbU - 1 - static_cast<unsigned int>(percentage * static_cast<double>(nbU) );

    if( index > curve.indexForLength(curveIndexL, constraint)){        // its within the correct boundaries
        curveIndexR = index;
        if(curveIndexR >= nbU) curveIndexR = nbU-1; // shouldn't ever happen
    }
    else if(curveIndexR == curve.indexForLength(curveIndexL, constraint)) return;
    else curveIndexR = curve.indexForLength(curveIndexL, constraint);

    movePlane(rightPlane, curveIndexR);
}

// Move the plane and send the distance between the two planes to the fibula
void Viewer::movePlane(Plane *p, unsigned int curveIndex){
    repositionPlane(p, curveIndex);
    double distance = curve.discreteChordLength(curveIndexL, curveIndexR);
    Q_EMIT planeMoved(distance);
    update();
}

void swap(unsigned int& a, unsigned int& b){
    unsigned int temp = a;
    a = b;
    b = temp;
}

unsigned int Viewer::partition(std::vector<unsigned int>& sorted, unsigned int start, unsigned int end){
    unsigned int p = sorted[end];
    unsigned int index = start - 1;

    for(unsigned int i=start; i<end; i++){
        double tangentAngleA = angle(curve.tangent(sorted[i]-1), curve.tangent(sorted[i]));
        double tangentAngleP = angle(curve.tangent(p-1), curve.tangent(p));

        if(tangentAngleA >= tangentAngleP){
            index++;
            swap(sorted[index], sorted[i]);
        }
    }
    swap(sorted[index+1], sorted[end]);
    return index+1;
}

void Viewer::quicksort(std::vector<unsigned int>& sorted, int start, int end){
    if(start < end){
        unsigned int p = partition(sorted, start, end);
        quicksort(sorted, start, static_cast<int>(p)-1);
        quicksort(sorted, static_cast<int>(p)+1, end);
    }
}

// Find where the ghost planes should be placed
void Viewer::findGhostLocations(unsigned int nbGhostPlanes, std::vector<unsigned int>& ghostLocation){
    if(nbGhostPlanes==0) return;
    unsigned int finalNb = nbGhostPlanes;        // the number we can actually fit in
    std::vector<unsigned int> maxIndicies(nbGhostPlanes);

    const unsigned int startI = curve.indexForLength(curveIndexL, constraint);
    const unsigned int endI = curve.indexForLength(curveIndexR, -constraint);

    if(endI > startI){         // if there's enough space for a plane
        const unsigned int searchArea = endI - startI;       // the space that's left between the left and right planes after the constraint is taken into account
        std::vector<unsigned int> sorted(static_cast<unsigned long long>(searchArea));
        for(unsigned int i=0; i<searchArea; i++) sorted[i] = startI+i;       // the possible indexes for the plane

        quicksort(sorted, 0, searchArea-1);      // Sort the indicies according to their tangent angles

        maxIndicies[0] = sorted[0];
        unsigned int sortedIndex = 1;

        for(unsigned int i=1; i<nbGhostPlanes; i++){
            // the constraint (don't take it if it's too close to another existing plane)
            bool tooClose;
            do{
                tooClose = false;
                for(int j=static_cast<int>(i)-1; j>=0; j--){
                    if(sortedIndex < static_cast<unsigned int>(searchArea) && curve.discreteLength(static_cast<unsigned int>(maxIndicies[static_cast<unsigned int>(j)]),static_cast<unsigned int>(sorted[static_cast<unsigned int>(sortedIndex)]))<constraint){
                        tooClose = true;
                        break;
                    }
                }
                if(tooClose) sortedIndex++;
            }while(tooClose);

            if(sortedIndex >= static_cast<unsigned int>(searchArea)){      // if we leave the search area, stop
                finalNb = i;
                break;
            }
            maxIndicies[i] = sorted[sortedIndex];
            sortedIndex++;  // move with i
        }

        // sort the planes
        for(unsigned int i=0; i<finalNb; i++){
            for(unsigned int j=i+1; j<finalNb; j++){
                if(maxIndicies[i] > maxIndicies[j]) swap(maxIndicies[i], maxIndicies[j]);
            }
        }

        ghostLocation.clear();
        for(unsigned int i=0; i<finalNb; i++) ghostLocation.push_back(maxIndicies[i]);       // get the location for each ghost plane
    }
}

double Viewer::angle(Vec a, Vec b){
    double na = a.norm();
    double nb = b.norm();
    double ab = a*b;

    double val = ab / (na*nb);
    if(val >= static_cast<double>(1)) val = 1;          // protection from floating point errors (comparing it to an epsilon didn't work)
    else if(val < static_cast<double>(-1)) val = -1;
    return acos(val);
}

void Viewer::rotatePolylineOnAxis(int position){
    double r = (position/360.*2.*M_PI);

    // send over the new norms
    sendNewNorms();
    Q_EMIT toRotatePolylineOnAxis(r);
    update();
}

void Viewer::sendNewNorms(){
    std::vector<Vec> norms;
    getPlaneBoxOrientations(norms);
    Q_EMIT toUpdatePlaneOrientations(norms);
}

/* Display settings */

void Viewer::toggleDrawMesh(){
    isDrawMesh = !isDrawMesh;
    update();
}

void Viewer::toggleWireframe(){
    poly.toggleIsWireframe();
    update();
}

void Viewer::toggleDrawPlane(){
    leftPlane->toggleIsVisible();
    rightPlane->toggleIsVisible();
    for(unsigned int i=0; i<ghostPlanes.size(); i++) ghostPlanes[i]->toggleIsVisible();
    update();
}

void Viewer::setMeshAlpha(int position){
    float p = static_cast<float>(position)/100.f;
    mesh.setAlpha(p);
    update();
}

void Viewer::setBoxAlpha(int position){
    float p = static_cast<float>(position)/100.f;
    poly.setAlpha(p);
    update();
}

void Viewer::setPlaneAlpha(int position){
    float p = static_cast<float>(position)/100.f;

    for(unsigned int i=0; i<ghostPlanes.size(); i++) ghostPlanes[i]->setAlpha(p);
    leftPlane->setAlpha(p);
    rightPlane->setAlpha(p);
    update();
}

void Viewer::cut(){
    mesh.setIsCut(Side::INTERIOR, true, true);
    Q_EMIT cutFibula();
    update();
}

void Viewer::uncut(){
    mesh.setIsCut(Side::INTERIOR, false, false);
    Q_EMIT uncutFibula();
    update();
}

void Viewer::initSignals(){
    connect(this, &Viewer::sendFibulaToMesh, &mesh, &Mesh::recieveInfoFromFibula);
    connect(&mesh, &Mesh::updateViewer, this, &Viewer::toUpdate);
}

void Viewer::recieveFromFibulaMesh(std::vector<int> &planes, std::vector<Vec> verticies, std::vector<std::vector<int>> &triangles, std::vector<int> &colours, std::vector<Vec> normals, int nbColours){
    /*
     * 0 : left plane
     * 1 : right plane
     * n : ghostPlane[n-2]
    */

    // For each vertex and normal, convert it from the corresponding plane's coordinates to the mesh coordinates
    for(unsigned int i=0; i<verticies.size(); i++){
        if(planes[i]==0){
            normals[i] = poly.getWorldBoxTransform(leftPlane->getID(), normals[i]);
            verticies[i] = poly.getWorldBoxCoordinates(leftPlane->getID(), verticies[i]);
        }
        else if(planes[i]==1){
            normals[i] = poly.getWorldBoxTransform(rightPlane->getID(), normals[i]);
            verticies[i] = poly.getWorldBoxCoordinates(rightPlane->getID(), verticies[i]);
        }
        else {
            int mandPlane = (planes[i]+2) / 2 - 2;
            planes[i] = (planes[i]+2) / 2 - 1;
            normals[i] = poly.getWorldBoxTransform(ghostPlanes[static_cast<unsigned int>(mandPlane)]->getID(), normals[i]);
            verticies[i] = poly.getWorldBoxCoordinates(ghostPlanes[static_cast<unsigned int>(mandPlane)]->getID(), verticies[i]);
        }
    }

    Q_EMIT sendFibulaToMesh(planes, verticies, triangles, colours, normals, nbColours);
}

void Viewer::toggleEditPlaneMode(unsigned int id, bool b){
    //leftPlane->toggleEditMode(b);
    //rightPlane->toggleEditMode(b);
    //for(unsigned int i=0; i<ghostPlanes.size(); i++) ghostPlanes[i]->toggleEditMode(b);
    ghostPlanes[id]->toggleEditMode(b);
    update();
}

Plane& Viewer::getPlaneFromID(unsigned int id){
    int idOffset = static_cast<int>(id) - static_cast<int>(ghostPlanes[0]->getID());
    if(idOffset < 0) return *leftPlane;
    return *ghostPlanes[static_cast<unsigned int>(idOffset)];
}

Plane& Viewer::getOppositePlaneFromID(unsigned int id){
    int idOffset = static_cast<int>(id) - static_cast<int>(leftPlane->getID());
    if(idOffset >= static_cast<int>(ghostPlanes.size())) return *rightPlane;
    return *ghostPlanes[static_cast<unsigned int>(idOffset)];
}

void Viewer::setBoxToManipulator(unsigned int id, Vec manipulatorPosition, int s){
    SavedState &ss = Q.back();
    poly.setBoxToManipulator(id, manipulatorPosition, s, ss.getBoxX(id), ss.getBoxY(id), ss.getBoxZ(id));

    double distShift;
    Vec projPoint;
    if(ghostPlanes.size() == 0) projPoint = projectBoxToPlane(*leftPlane, *rightPlane, distShift);
    else projPoint = projectBoxToPlane(getPlaneFromID(id), getOppositePlaneFromID(id), distShift);
    poly.setBoxToProjectionPoint(id, projPoint);
    poly.adjustBoxLength(id, distShift); // Shift the box size

    sendNewNorms();     // Need to send new norms and the new rotation of the box - here its the opposite of before -- don't update on this side but send the update to the fibula correspondance
    std::vector<double> distances;
    poly.getDistances(distances);
    Q_EMIT toReinitBox(id, distances);
    update();
}

void Viewer::setBoxToCornerManipulator(unsigned int id, Vec manipulatorPosition){
    poly.setBoxToCornerManipulator(id, manipulatorPosition);

    id = id/2;

    double distShift;
    Vec projPoint;
    if(ghostPlanes.size() == 0) projPoint = projectBoxToPlane(*leftPlane, *rightPlane, distShift);
    else projPoint = projectBoxToPlane(getPlaneFromID(id), getOppositePlaneFromID(id), distShift);
    poly.setBoxToProjectionPoint(id, projPoint);
    poly.adjustBoxLength(id, distShift); // Shift the box size

    sendNewNorms();     // Need to send new norms and the new rotation of the box - here its the opposite of before -- don't update on this side but send the update to the fibula correspondance
    std::vector<double> distances;
    poly.getDistances(distances);
    Q_EMIT toReinitBox(id, distances);
    update();
}

void Viewer::toggleEditBoxMode(unsigned int id, bool b){
    poly.toggleBoxManipulators(id, b);
    update();
}

void Viewer::toggleEditFirstCorner(unsigned int id, bool b){
    poly.toggleFirstCornerManipulators(id, b);
    update();
}

void Viewer::toggleEditEndCorner(unsigned int id, bool b){
    poly.toggleEndCornerManipulators(id, b);
    update();
}

void Viewer::toggleAllPlanes(bool b){
    for(unsigned int i=0; i<ghostPlanes.size(); i++) ghostPlanes[i]->toggleEditMode(b);
    update();
}

void Viewer::toggleAllBoxes(bool b){
    poly.activateBoxManipulators(b);
    update();
}

void Viewer::toggleAllFirstCorners(bool b){
    poly.activateFirstCornerManipulators(b);
    update();
}

void Viewer::toggleAllEndCorners(bool b){
    poly.activateEndCornerManipulators(b);
    update();
}

void Viewer::toggleDrawBoxes(){
    poly.toggleDrawBoxes();
    update();
}

void Viewer::toggleDrawPolyline(){
    poly.toggleDrawLine();
    update();
}

void Viewer::toggleDrawCurve(){
    isDrawCurve = !isDrawCurve;
    update();
}

Vec Viewer::projectBoxToPlane(Plane &p, Plane &endP, double &distShift){
    unsigned int boxIndex = p.getID();

    // Side one projection
    Vec worldBox = poly.getMeshBoxPoint(boxIndex);
    Vec worldProjectionAxis = worldBox - poly.getMeshBoxEnd(boxIndex);
    Vec projPoint = p.getAxisProjection(worldBox, worldProjectionAxis);

    // Side two projection
    Vec worldBoxEnd = poly.getMeshBoxEnd(boxIndex);
    Vec projEnd = endP.getAxisProjection(worldBoxEnd, -worldProjectionAxis);

    // High projection
    /*Vec worldBoxHighPoint = poly.getMeshBoxHighPoint(boxIndex);
    Vec projHighPoint = p.getAxisProjection(worldBoxHighPoint, worldProjectionAxis);

    // High end projection
    Vec worldBoxHighEnd = poly.getMeshBoxHighEnd(boxIndex);
    Vec projHighEnd = endP.getAxisProjection(worldBoxHighEnd, -worldProjectionAxis);

    distShift = maxDouble(euclideanDistance(projPoint, projEnd), euclideanDistance(projHighPoint, projHighEnd));*/

    distShift = euclideanDistance(projPoint, projEnd);

    return projPoint;
}

double Viewer::maxDouble(double a, double b){
    if(a>b) return a;
    else return b;
}

double Viewer::euclideanDistance(const Vec &a, const Vec &b){
    return sqrt(pow(a.x-b.x, 2.) + pow(a.y-b.y, 2.) + pow(a.z-b.z, 2.));
}

bool Viewer::isViolatesContraint(){
    return false;
}

void Viewer::saveCurrentState(Modification m){
    if( Q.size() == max_operation_saved )
        Q.pop_front();

    Q.push_back(saveState(m));
}

SavedState Viewer::saveState(Modification m){
    SavedState s = SavedState(m);
    s.addPlanes(ghostPlanes);
    s.addBoxes(poly.getBoxes());
    return s;
}


void Viewer::restoreLastState(){

    if( Q.size() > 0 ){

        if(Q.size() > 1) Q.pop_back();

        resetState(Q.back());
    }

}

void Viewer::resetState(SavedState s){
    std::vector<Vec> bX, bY, bZ, bPos, pPos;
    std::vector<double> bl;
    std::vector<Quaternion> pOrient;

    if(s.getModification()==Modification::PLANE){
        s.getPlanePositions(pPos);
        s.getPlaneOrientations(pOrient);
        for(unsigned int i=0; i<ghostPlanes.size(); i++){
            ghostPlanes[i]->setPosition(pPos[i]);
            ghostPlanes[i]->setOrientation(pOrient[i]);
        }
        poly.setPointsToGhostPlanes(pPos);
    }

    s.getBoxLengths(bl);
    s.getBoxPositions(bPos);
    s.getBoxXOrientations(bX);
    s.getBoxYOrientations(bY);
    s.getBoxZOrientations(bZ);

    for(unsigned int i=0; i<bl.size(); i++){
        poly.setBoxLocation(i, bPos[i]);
        poly.adjustBoxLength(i, bl[i]);
        poly.setManipulatorOrientations(bX, bY, bZ);
        poly.setBoxToManipulatorOrientation(i);
    }
    poly.setManipulatorsToBoxes();

    sendNewNorms();
    std::vector<double> distances;
    poly.getDistances(distances);
    Q_EMIT toReinitBox(0, distances);       // TODO the id index isn't currently used
    update();
}

void Viewer::resetUndoQueue(){
    Q.clear();
}

void Viewer::manipulatorReleasedPlane(){
    saveCurrentState(Modification::PLANE);
}

void Viewer::manipulatorReleasedBox(){
    saveCurrentState(Modification::BOX);
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
    switch (e->key())
    {
    case Qt::Key_Z :
        if(e->modifiers() & Qt::ControlModifier) restoreLastState(); update(); break;
    default : QGLViewer::keyPressEvent(e);
    }
}

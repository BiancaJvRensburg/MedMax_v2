#include "viewer.h"
#include "Triangle.h"
#include "Vec3D.h"
#include <QGLViewer/manipulatedFrame.h>

Viewer::Viewer(QWidget *parent, StandardCamera *cam, int sliderMax) : QGLViewer(parent) {
    Camera *c = camera();       // switch the cameras
    setCamera(cam);
    delete c;

    this->nbU = 0;
    this->sliderMax = sliderMax;
    this->isDrawMesh = false;
    this->nbGhostPlanes = 0;
    this->isGhostPlanes = false;
    this->isGhostActive = true;
    this->isCurve = false;
}

void Viewer::draw() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();
    glMultMatrixd(manipulatedFrame()->matrix());

    glColor3f(1.,1.,1.);
    mesh.draw();
    if(isDrawMesh) mesh.drawCut();   // draw the cut versions

    if(isCurve){
        //if(isGhostPlanes && isGhostActive) drawPolyline();

        // draw the planes
        glColor4f(1.0, 0, 0, leftPlane->getAlpha());
        leftPlane->draw();

        glColor4f(0, 1.0, 0, rightPlane->getAlpha());
        rightPlane->draw();

        for(unsigned int i=0; i<ghostPlanes.size(); i++){       // draw the ghost planes
            glColor4f(0,0,1.0, ghostPlanes[i]->getAlpha());
            ghostPlanes[i]->draw();
        }

        glDisable(GL_DEPTH);
        glDisable(GL_DEPTH_TEST);

         curve->draw();
         //curve->drawControl();
         //for(unsigned int i=0; i<curve->getNbU(); i++) curve->drawTangent(i);
    }

    glPopMatrix();
}

void Viewer::toggleIsDrawPlane(){
    leftPlane->toggleIsVisible();
    rightPlane->toggleIsVisible();
    for(unsigned int i=0; i<ghostPlanes.size(); i++) ghostPlanes[i]->toggleIsVisible();
    update();
}

// Updates the polyline vector and angles when a plane is moved (sends the angles to the fibula)
std::vector<Vec> Viewer::updatePolyline(){

    // NOTE if this is in comments, the planes always line up on the two sides
    // if its not cut / are no ghost planes
    if(!isGhostPlanes){
        std::vector<Vec> angles;
        return angles;   // return an empty vector
    }

    std::vector<Vec> polyline;
    updateMeshPolyline(polyline);

    return getPolylinePlaneAngles(polyline);
}

void Viewer::updateMeshPolyline(std::vector<Vec> &polyline){
    const Vec &p0 = leftPlane->getPosition();
    polyline.push_back(p0);

    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        const Vec &pGhost = Vec(ghostPlanes[i]->getCurvePoint().getX(), ghostPlanes[i]->getCurvePoint().getY(), ghostPlanes[i]->getCurvePoint().getZ());
        polyline.push_back(pGhost);
    }

    const Vec &p1 = rightPlane->getPosition();
    polyline.push_back(p1);
}

// Get each polyline in the coordinates of its plane
std::vector<Vec> Viewer::getPolylinePlaneAngles(std::vector<Vec> &polyline){
    std::vector<Vec> angles;        // There are always 2 intersections for each polyline segment

    if(polyline.size()==0) return angles;

    angles.push_back(leftPlane->getPolylineVector(polyline[1]));

    for(unsigned int i=1; i<polyline.size()-1; i++){
        angles.push_back(ghostPlanes[i-1]->getPolylineVector(polyline[i-1]));
        angles.push_back(ghostPlanes[i-1]->getPolylineVector(polyline[i+1]));
    }

    angles.push_back(rightPlane->getPolylineVector(polyline[polyline.size()-2]));

    return angles;
}

/*void Viewer::drawPolyline(){
    glEnable(GL_DEPTH);
    glEnable(GL_DEPTH_TEST);

    glLineWidth(5.0);

    glBegin(GL_LINE_STRIP);

        for(unsigned int i=0; i<polyline.size(); i++){      // Looks at the previous
            if(i==0 || segmentLength(polyline[i], polyline[i-1]) > constraint) glColor3f(0.0, 0.0, 1.0);        // if the distance is greater than the constraint : blue
            else glColor3f(1.0, 0.0, 0.0);      // if it violates the constraint, then red
            glVertex3f(static_cast<float>(polyline[i].x), static_cast<float>(polyline[i].y), static_cast<float>(polyline[i].z));
        }

    glEnd();

    glLineWidth(1.0);
    glColor3f(1.,1.,1.);    // reset the colour

    glDisable(GL_DEPTH);
    glDisable(GL_DEPTH_TEST);
}*/

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

  glEnable(GL_LIGHTING);
  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
  glLineWidth (1.0f);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);

  initSignals();
}

void Viewer::initSignals(){
    connect(this, &Viewer::sendFibulaToMesh, &mesh, &Mesh::recieveInfoFromFibula);
    connect(&mesh, &Mesh::updateViewer, this, &Viewer::toUpdate);
}

void Viewer::recieveFromFibulaMesh(const std::vector<int> &planes, std::vector<Vec> verticies, const std::vector<std::vector<int>> &triangles, const std::vector<int> &colours, std::vector<Vec> normals, const int nbColours){
    /*
     * 0 : left plane
     * 1 : right plane
     * n : ghostPlane[n-2]
    */

    // For each vertex and normal, convert it from the corresponding plane's coordinates to the mesh coordinates
    for(unsigned int i=0; i<verticies.size(); i++){
        if(planes[i]==0){
            normals[i] = leftPlane->getMeshVectorFromLocal(normals[i]);
            verticies[i] = leftPlane->getMeshCoordinatesFromLocal(verticies[i]);
        }
        else if(planes[i]==1){
            normals[i] = rightPlane->getMeshVectorFromLocal(normals[i]);
            verticies[i] = rightPlane->getMeshCoordinatesFromLocal(verticies[i]);
        }
        else {
            int mandPlane = (planes[i]+2) / 2 - 2;
            normals[i] = ghostPlanes[static_cast<unsigned int>(mandPlane)]->getMeshVectorFromLocal(normals[i]);
            verticies[i] = ghostPlanes[static_cast<unsigned int>(mandPlane)]->getMeshCoordinatesFromLocal(verticies[i]);
        }
    }

    Q_EMIT sendFibulaToMesh(verticies, triangles, colours, normals, nbColours);
}

QString Viewer::helpString() const {
  QString text("<h2>H e l p  B o x</h2>");
  text += "This is the help text";
  return text;
}

void swap(int& a, int& b){
    int temp = a;
    a = b;
    b = temp;
}

int Viewer::partition(int sorted[], int start, int end){
    int p = sorted[end];
    int index = start - 1;

    for(int i=start; i<end; i++){
        double tangentAngleA = angle(curve->tangent(sorted[i]-1), curve->tangent(sorted[i]));
        double tangentAngleP = angle(curve->tangent(p-1), curve->tangent(p));

        if(tangentAngleA >= tangentAngleP){
            index++;
            swap(sorted[index], sorted[i]);
        }
    }
    swap(sorted[index+1], sorted[end]);
    return index+1;
}

void Viewer::quicksort(int sorted[], int start, int end){
    if(start < end){
        int p = partition(sorted, start, end);
        quicksort(sorted, start, p-1);
        quicksort(sorted, p+1, end);
    }
}

void Viewer::drawMesh(){
    if(isDrawMesh) isDrawMesh = false;      // change states
    else isDrawMesh = true;
    update();
}

bool Viewer::isSpaceForGhosts(){
    int startI = static_cast<int>(curve->indexForLength(curveIndexL, constraint));
    int endI = static_cast<int>(curve->indexForLength(curveIndexR, -constraint));
    int searchArea = endI - startI;
    if(searchArea>0) return true;
    return false;
}

void Viewer::initGhostPlanes(){
    for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];     // get rid of any previous ghost planes
    ghostPlanes.clear();

    int finalNb = nbGhostPlanes;        // the number we can actually fit in
    int maxIndicies[nbGhostPlanes];

    const int startI = static_cast<const int>(curve->indexForLength(curveIndexL, constraint));
    const int endI = static_cast<const int>(curve->indexForLength(curveIndexR, -constraint));
    const int searchArea = endI - startI;       // the space that's left between the left and right planes after the constraint is taken into account


    if(searchArea > 0){         // if there's enough space for a plane
        int sorted[searchArea];

        for(int i=0; i<searchArea; i++) sorted[i] = startI+i;       // the possible indexes for the plane

        quicksort(sorted, 0, searchArea-1);      // Sort the indicies according to their tangent angles

        maxIndicies[0] = sorted[0];
        int sortedIndex = 1;

        for(int i=1; i<nbGhostPlanes; i++){
            // the constraint (don't take it if it's too close to another existing plane)
            bool tooClose;
            do{
                tooClose = false;
                for(int j=i-1; j>=0; j--){
                    if(sortedIndex < searchArea && curve->discreteLength(static_cast<unsigned int>(maxIndicies[static_cast<unsigned int>(j)]),static_cast<unsigned int>(sorted[static_cast<unsigned int>(sortedIndex)]))<constraint){
                        tooClose = true;
                        break;
                    }
                }
                if(tooClose) sortedIndex++;
            }while(tooClose);

            if(sortedIndex >= searchArea){      // if we leave the search area, stop
                finalNb = i;
                break;
            }
            maxIndicies[i] = sorted[sortedIndex];
            sortedIndex++;  // move with i
        }

        // sort the planes
        for(int i=0; i<finalNb; i++){
            for(int j=i+1; j<finalNb; j++){
                if(maxIndicies[i] > maxIndicies[j]) swap(maxIndicies[i], maxIndicies[j]);
            }
        }

        ghostLocation.clear();
        for(int i=0; i<finalNb; i++) ghostLocation.push_back(maxIndicies[i]);       // get the location for each ghost plane

        // the number of ghost planes we can currently fit
        currentNbGhostPlanes = finalNb;

        addGhostPlanes(finalNb);

        for(unsigned int i=0; i<ghostPlanes.size(); i++) connect(&(ghostPlanes[i]->getCurvePoint()), &CurvePoint::curvePointTranslated, this, &Viewer::ghostPlaneMoved);        // connnect the ghost planes

         // Send the info to the fibula
        const std::vector<Vec> &poly = updatePolyline();
        const std::vector<Vec> &axes = getReferenceAxes();
        double distance;        // the distance we moved

        if(finalNb > 0) distance = curve->discreteLength(curveIndexL, static_cast<unsigned int>(ghostLocation[0]));
        else distance = curve->discreteLength(curveIndexL, curveIndexR);
        Q_EMIT leftPosChanged(distance, poly, axes);

        if(finalNb > 0) distance = curve->discreteLength(curveIndexR, ghostLocation[finalNb-1]);
        else distance = curve->discreteLength(curveIndexL, curveIndexR);
        Q_EMIT rightPosChanged(distance, poly, axes);
    }
}

void Viewer::cutMesh(){
    bool isNumberRecieved;
    int nbPieces = QInputDialog::getInt(this, "Cut mesh", "Number of pieces", 0, 1, 10, 1, &isNumberRecieved, Qt::WindowFlags());
    if(isNumberRecieved) nbGhostPlanes = nbPieces-1;
    else return;

    isGhostPlanes = true;

    Q_EMIT preparingToCut();

    if(nbGhostPlanes==0) Q_EMIT noGhostPlanesToSend(updatePolyline(), getReferenceAxes(), curve->discreteLength(curveIndexL, curveIndexR));

    Q_EMIT okToCut();

    mesh.setIsCut(Side::INTERIOR, true, true);      // cut the mandible mesh and initialise the ghost planes
    initGhostPlanes();

    update();
}

void Viewer::uncutMesh(){
    mesh.setIsCut(Side::INTERIOR, false, false);
    isGhostPlanes = false;
    for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];
    ghostPlanes.clear();
    update();
}

// Slide the left plane
void Viewer::moveLeftPlane(int position){
    double percentage = static_cast<double>(position) / static_cast<double>(sliderMax);
    unsigned int index = static_cast<unsigned int>(percentage * static_cast<double>(nbU) );

    if(curve->indexForLength(curveIndexR, -constraint) > index){  // Only move if we're going backwards or we haven't met the other plane
        curveIndexL = index;
        if(curveIndexL >= nbU) curveIndexL = nbU-1;     // shouldn't ever happen
    }
    else if( curveIndexL == curve->indexForLength(curveIndexR, -constraint) ){
        return;       // already in the correct position
    }
    else curveIndexL = curve->indexForLength(curveIndexR, -constraint);     // get the new position

    movePlane(leftPlane, true, curveIndexL);
}

void Viewer::moveRightPlane(int position){
    double percentage = static_cast<double>(position) / static_cast<double>(sliderMax);
    unsigned int index = nbU - 1 - static_cast<unsigned int>(percentage * static_cast<double>(nbU) );

    if( index > curve->indexForLength(curveIndexL, constraint)){        // its within the correct boundaries
        curveIndexR = index;
        if(curveIndexR >= nbU) curveIndexR = nbU-1; // shouldn't ever happen
    }
    else if(curveIndexR == curve->indexForLength(curveIndexL, constraint)) return;
    else curveIndexR = curve->indexForLength(curveIndexL, constraint);

    movePlane(rightPlane, false, curveIndexR);
}

void Viewer::movePlane(Plane *p, bool isLeft, unsigned int curveIndex){
    repositionPlane(p, curveIndex);

    mesh.updatePlaneIntersections(p);       // update the mesh
    update();

    // Communication with the fibula
    if(isGhostPlanes && isGhostActive) handlePlaneMoveStart();

    double distance;

    /*if(!(isGhostPlanes && isGhostActive)) distance = curve->discreteChordLength(curveIndexL, curveIndexR);       // get the distance of the mandible curve
    else*/ if(ghostPlanes.size()==0) distance = curve->discreteLength(curveIndexL, curveIndexR);  // is cut but no ghost planes
    else{
        if(isLeft) distance = curve->discreteLength(curveIndexL, ghostLocation[0]);
        else distance = curve->discreteLength(ghostLocation[ghostPlanes.size()-1], curveIndexR);
    }

    const std::vector<Vec> &poly = updatePolyline();
    const std::vector<Vec> &axes = getReferenceAxes();

    if(isLeft){
        Q_EMIT leftPosChanged(distance, poly, axes);
        Q_EMIT setLRSliderValue(0);     // Reset the rotation slider
    }
    else{
        Q_EMIT rightPosChanged(distance, poly, axes);
        Q_EMIT setRRSliderValue(0);
    }
}

void Viewer::onLeftSliderReleased(){
    if(isGhostPlanes) handlePlaneMoveEnd();
    // Q_EMIT setLMSliderValue( static_cast<int>( (static_cast<double>(curveIndexL)/static_cast<double>(nbU)) * static_cast<double>(sliderMax) ) );
}

void Viewer::onRightSliderReleased(){
    if(isGhostPlanes) handlePlaneMoveEnd();
    // Q_EMIT setRMSliderValue( static_cast<int>( sliderMax - (static_cast<double>(curveIndexR)/static_cast<double>(nbU)) * static_cast<double>(sliderMax) ) );
}


// Linked to the sliders
void Viewer::rotateLeftPlane(int position){
    rotatePlane(leftPlane, position);
}

void Viewer::rotateRightPlane(int position){
    rotatePlane(rightPlane, position);
}

void Viewer::rotatePlane(Plane *p, int position){
    double percentage = static_cast<double>(position) / static_cast<double>(sliderMax);

    p->rotatePlaneXY(percentage);
    mesh.updatePlaneIntersections(p);
    update();
}

void Viewer::handlePlaneMoveStart(){
    isGhostActive = false;      // disable drawing the polyline
    mesh.setIsCut(Side::INTERIOR, false, false);        // Uncut the mesh but keep the ghostPlane marker as true
    for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];
    ghostPlanes.clear();

    Q_EMIT ghostPlaneMovementStart();
}

void Viewer::handlePlaneMoveEnd(){
    isGhostActive = true;       // enable the polyline

    // Recut
    Q_EMIT preparingToCut();
    if(!isSpaceForGhosts() || nbGhostPlanes==0) Q_EMIT noGhostPlanesToSend(updatePolyline(), getReferenceAxes(), curve->discreteLength(curveIndexL, curveIndexR));
    Q_EMIT okToCut();

    mesh.setIsCut(Side::INTERIOR, true, true);
    isGhostPlanes = true;
    initGhostPlanes();
}

void Viewer::openOFF(QString filename) {
    std::vector<Vec3Df> &vertices = mesh.getVertices();
    std::vector<Triangle> &triangles = mesh.getTriangles();

    FileIO::openOFF(filename.toStdString(), vertices, triangles);

    mesh.init();

    // Set the camera
    Vec3Df center;
    double radius;
    MeshTools::computeAveragePosAndRadius(vertices, center, radius);
    updateCamera(center, static_cast<float>(radius));

    update();
}

void Viewer::initCurve(){    
    control.push_back(Vec(-56.9335, -13.9973, 8.25454));

    control.push_back(Vec(-50.8191, -20.195, -19.53));
    control.push_back(Vec(-40.155, -34.5957, -50.7005));
    control.push_back(Vec(-27.6007, -69.2743, -67.6769));

    control.push_back(Vec(0, -85.966, -68.3154));

    control.push_back(Vec(26.7572, -69.0705, -65.6261));
    control.push_back(Vec(40.3576, -34.3609, -50.7634));
    control.push_back(Vec(46.2189, -21.3245, -17.9009));

    control.push_back(Vec(52.3669, -15.4613, 8.70223));

    constructCurve();
}

void Viewer::constructCurve(){
    curve = new Curve(control.size(), control);
    nbU = 100;
    curve->generateCatmull(nbU);
    isCurve = true;
    initPlanes(Movable::DYNAMIC);
}

void Viewer::initPlanes(Movable status){
    curveIndexR = nbU - 1;
    curveIndexL = 0;
    Vec pos(0,0,0);
    float size = 40.0;

    leftPlane = new Plane(static_cast<double>(size), status, pos, 0.5);
    rightPlane = new Plane(static_cast<double>(size), status, pos, 0.5);

    repositionPlane(leftPlane, curveIndexL);
    repositionPlane(rightPlane, curveIndexR);

    mesh.addPlane(leftPlane);
    mesh.addPlane(rightPlane);

    connect(&(leftPlane->getCurvePoint()), &CurvePoint::curvePointTranslated, this, &Viewer::ghostPlaneMoved);
    connect(&(rightPlane->getCurvePoint()), &CurvePoint::curvePointTranslated, this, &Viewer::ghostPlaneMoved);
}

void Viewer::repositionPlane(Plane *p, unsigned int index){
    p->setPosition(curve->getPoint(index));
    matchPlaneToFrenet(p, index);
}

void Viewer::addGhostPlanes(unsigned int nb){
    if(nb==0) return;
    double distances[nb+1];     // +1 for the last plane

    for(unsigned int i=0; i<static_cast<unsigned int>(nb); i++){
        Vec pos(0,0,0);
        ghostPlanes.push_back(new Plane(40.0, Movable::DYNAMIC, pos, leftPlane->getAlpha()));
        repositionPlane(ghostPlanes[i], ghostLocation[i]);
        if(i==0) distances[i] = curve->discreteLength(curveIndexL, ghostLocation[i]);
        else distances[i] = curve->discreteLength(ghostLocation[i-1], ghostLocation[i]);
    }

    balanceGhostPlanes();

    distances[nb] = curve->discreteLength(ghostLocation[static_cast<unsigned int>(nb-1)], curveIndexR);

    const std::vector<Vec> &poly = updatePolyline();
    const std::vector<Vec> &axes = getReferenceAxes();

    Q_EMIT ghostPlanesAdded(nb, distances, poly, axes);
}

void Viewer::ghostPlaneMoved(){
    unsigned int nb = static_cast<unsigned int>(ghostPlanes.size());
    double distances[nb+1];     // +1 for the last plane

    for(unsigned int i=0; i<nb; i++){
        if(i==0) distances[i] = segmentLength(leftPlane->getPosition(), ghostPlanes[i]->getCurvePoint().getPoint());
        else distances[i] = segmentLength(ghostPlanes[i-1]->getCurvePoint().getPoint(), ghostPlanes[i]->getCurvePoint().getPoint());
    }

    balanceGhostPlanes();

    distances[nb] = segmentLength(rightPlane->getPosition(), ghostPlanes[nb-1]->getCurvePoint().getPoint());

    const std::vector<Vec> &poly = updatePolyline();
    const std::vector<Vec> &axes = getReferenceAxes();
    Q_EMIT ghostPlanesTranslated(nb, distances, poly, axes);
}

void Viewer::updateCamera(const Vec3Df & center, float radius){
    camera()->setSceneCenter(Vec(static_cast<double>(center[0]), static_cast<double>(center[1]), static_cast<double>(center[2])));
    camera()->setSceneRadius(static_cast<double>(radius*1.05f));
    camera()->setZClippingCoefficient(static_cast<double>(radius/8.0f));
    camera()->showEntireScene();
}

void Viewer::updatePlanes(){
    repositionPlane(leftPlane, curveIndexL);
    repositionPlane(rightPlane, curveIndexR);

    mesh.updatePlaneIntersections();

    update();
}

Quaternion Viewer::getNewOrientation(unsigned int index){
    Quaternion s(Vec(0,0,1.0), curve->tangent(index));
    return s.normalized();
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

double Viewer::segmentLength(const Vec a, const Vec b){
    return sqrt( pow((b.x - a.x), 2) + pow((b.y - a.y), 2) + pow((b.z - a.z), 2));
}

void Viewer::matchPlaneToFrenet(Plane *p, unsigned int index){
    Vec x, y, z;
    curve->getFrame(index,z,x,y);
    p->setFrameFromBasis(x,y,z);
}

Vec Viewer::convertToPlane(Plane *base, Plane *p, Vec axis){
    Vec  a = p->getMeshVectorFromLocal(axis);  // get the z vector of p in the world space
    Vec b = base->getLocalVector(a);
    return b;    // get it in terms of base
}

void Viewer::rotateFrame(Frame &f, Vec axis, double theta){
    f.rotate(Quaternion(cos(theta/2.0)*axis.x, cos(theta/2.0)*axis.y, cos(theta/2.0)*axis.z, sin(theta/2.0)));
}

void Viewer::getAxes(){
    Q_EMIT sendAxes(getReferenceAxes());
}

Vec Viewer::getCustomProjection(Vec a, Vec normal){
    return a - normal * (a * normal);
}

void Viewer::addFrameChangeToAxes(std::vector<Vec> &axes, Plane *base, Plane *p){
    Vec axisX(1,0,0);
    Vec axisY(0,1,0);
    Vec axisZ(0,0,1);

    axes.push_back(convertToPlane(base, p, axisX));
    axes.push_back(convertToPlane(base, p, axisY));
    axes.push_back(convertToPlane(base, p, axisZ));
}

void Viewer::addInverseFrameChangeToAxes(std::vector<Vec> &axes, Plane *base, Plane *p){
    Vec axisX(1,0,0);
    Vec axisY(0,1,0);
    Vec axisZ(0,0,1);

    axes.push_back(convertToPlane(base, p, axisX));
    axes.push_back(convertToPlane(base, p, -axisY));
    axes.push_back(convertToPlane(base, p, -axisZ));
}

std::vector<Vec> Viewer::getReferenceAxes(){
    std::vector<Vec> v;

    if(ghostPlanes.size()==0){
        addFrameChangeToAxes(v, rightPlane, leftPlane);
    }
    else{
        addInverseFrameChangeToAxes(v, leftPlane, ghostPlanes[0]);

       for(unsigned int i=0; i<ghostPlanes.size()-1; i++){      // don't include the last ghost plane, the right plane will set it
            addInverseFrameChangeToAxes(v, ghostPlanes[i], ghostPlanes[i+1]);
        }

        unsigned int lastIndex = static_cast<unsigned int>(ghostPlanes.size())-1;
        addFrameChangeToAxes(v, rightPlane, ghostPlanes[lastIndex]);
    }
    return v;
}

void Viewer::readJSON(const QJsonArray &controlArray){
    control.clear();

    for(int i=0; i<controlArray.size(); i++){
        QJsonArray singleControl = controlArray[i].toArray();
        control.push_back(Vec(singleControl[0].toDouble(), singleControl[1].toDouble(), singleControl[2].toDouble()));
    }
    constructCurve();

    Vec3Df center;
    float radius;
    mesh.computeBB(center, radius);
    updateCamera(center, radius);
}

void Viewer::setAlpha(int position){
    float a = static_cast<float>(position) / 100.f;
    for(unsigned int i=0; i<ghostPlanes.size(); i++) ghostPlanes[i]->setAlpha(a);
    leftPlane->setAlpha(a);
    rightPlane->setAlpha(a);
    update();
}

void Viewer::balanceGhostPlanes(){
    Vec norm(0,0,1);
    const std::vector<Vec> &poly = updatePolyline();

    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        double theta = M_PI - angle(poly[i*2+1], -norm);
        double fi = M_PI - angle(poly[i*2+2], norm);
        double shift = (fi - theta) / 2.0;
        ghostPlanes[i]->rotatePlane(Vec(0,1,0), shift+M_PI);
    }
}

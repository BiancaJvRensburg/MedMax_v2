#include "viewer.h"
#include "Triangle.h"
#include "Vec3D.h"
//#include "controlpoint.h"
#include <QGLViewer/manipulatedFrame.h>

Viewer::Viewer(QWidget *parent, StandardCamera *cam, int sliderMax) : QGLViewer(parent) {
    // Change the camera.
    Camera *c = camera();
    setCamera(cam);
    delete c;

    this->nbU = 0;
    this->sliderMax = sliderMax;
    this->isDrawMesh = false;
    this->nbGhostPlanes = 3;
    this->isGhostPlanes = false;
    this->isGhostActive = true;
}

void Viewer::draw() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();
    glMultMatrixd(manipulatedFrame()->matrix());
    //drawAxis(20.0);

    glColor3f(1.,1.,1.);
    mesh.draw();
    if(isDrawMesh) mesh.drawCut();

    if(isGhostPlanes && isGhostActive) drawPolyline();

    //glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

    glColor3f(1.0, 0, 0);
    leftPlane->draw();

    glColor3f(0, 1.0, 0);
    rightPlane->draw();

    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        glColor3f(0,0,1.0);
        ghostPlanes[i]->draw();
    }

    //glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

    curve->draw();

     // curve->drawControl(); // We want to visualise this at all times
     curve->drawTangent(curveIndexL);

    glPopMatrix();
}

// Updates the polyline vector and angles when a plane is moved (sends the angles to the fibula)
std::vector<Vec> Viewer::updatePolyline(){

    // NOTE if this is in comments, the planes always line up on the two sides
    // if its not cut / are no ghost planes
    if(!isGhostPlanes || !isGhostActive){
        std::vector<Vec> angles;
        return angles;   // return an empty vector
    }

    polyline.clear();

    Vec p = leftPlane->getPosition();
    polyline.push_back(p);

    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        //if(ghostPlanes[i].getCurvePoint()==nullptr) continue;
        p = Vec(ghostPlanes[i]->getCurvePoint().getX(), ghostPlanes[i]->getCurvePoint().getY(), ghostPlanes[i]->getCurvePoint().getZ());
        polyline.push_back(p);
    }

    p = rightPlane->getPosition();
    polyline.push_back(p);

    return getPolylinePlaneAngles();
}

void Viewer::drawPolyline(){
    glEnable(GL_DEPTH);
    glEnable(GL_DEPTH_TEST);

    glLineWidth(5.0);

    glBegin(GL_LINE_STRIP);

        for(unsigned int i=0; i<polyline.size(); i++){
            // Looks at the previous
            // if the distance is greater than the constraint : blue
            if(i==0 || segmentLength(polyline[i], polyline[i-1]) > constraint) glColor3f(0.0, 0.0, 1.0);
            // if it violates the constraint, then red
            else glColor3f(1.0, 0.0, 0.0);
            glVertex3f(static_cast<float>(polyline[i].x), static_cast<float>(polyline[i].y), static_cast<float>(polyline[i].z));
        }

    glEnd();

    glLineWidth(1.0);
    glColor3f(1.,1.,1.);    // reset the colour

    glDisable(GL_DEPTH);
    glDisable(GL_DEPTH_TEST);
}

void Viewer::toUpdate(){
    update();
}

void Viewer::init() {
  setMouseTracking(true);

  restoreStateFromFile();

  viewerFrame = new ManipulatedFrame();

  setManipulatedFrame(viewerFrame);

  setAxisIsDrawn(false);

  initCurve();

  // Set up gl settings
  glEnable(GL_LIGHTING);

  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
  //glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

  glLineWidth (1.0f);

  initSignals();
}

void Viewer::initSignals(){
    connect(this, &Viewer::sendFibulaToMesh, &mesh, &Mesh::recieveInfoFromFibula);
    connect(&mesh, &Mesh::updateViewer, this, &Viewer::toUpdate);
}

void Viewer::recieveFromFibulaMesh(std::vector<int> planes, std::vector<Vec> verticies, std::vector<std::vector<int>> triangles, std::vector<Vec> fibulaPolyline, std::vector<int> colours, std::vector<Vec> normals, int nbColours){
    // TODO : Be careful of this!
    /*
     * 0 : left plane
     * 1 : right plane
     * n : ghostPlane[n-2]
    */

    // Rotate the planes to match the fibula (just the extremities for now)
   /* Vec mandPolylineSegment;
    Vec axis = Vec(0,0,1);
    Quaternion s;
    double anglePoly;

    updatePolyline();

    bool isRotated[planes.size()];
    for(int i=0; i<planes.size(); i++) isRotated[i]=false;*/

    // For each vertex and normal, convert it from the corresponding plane's coordinates to the mesh coordinates
    for(unsigned int i=0; i<verticies.size(); i++){
        if(planes[i]==0){
            //Rotate the plane if it isn't already
            /*if(!isRotated[planes[i]]){
                leftPlane->constrainZRotation();
                leftPlane->setOrientation(Quaternion(fibulaPolyline[0], polyline[1]));
                isRotated[planes[i]] = true;
                leftPlane->freeZRotation();
            }*/
            normals[i] = leftPlane->getMeshVectorFromLocal(normals[i]);
            verticies[i] = leftPlane->getMeshCoordinatesFromLocal(verticies[i]);
        }
        else if(planes[i]==1){
            /*if(!isRotated[planes[i]]){
                rightPlane->constrainZRotation();
                rightPlane->setOrientation(Quaternion(fibulaPolyline[fibulaPolyline.size()-2], polyline[polyline.size()-2]));
                isRotated[planes[i]] = true;
                rightPlane->freeZRotation();
            }*/
            normals[i] = rightPlane->getMeshVectorFromLocal(normals[i]);
            verticies[i] = rightPlane->getMeshCoordinatesFromLocal(verticies[i]);
        }
        else {
            int mandPlane = (planes[i]+2) / 2 - 2;
            normals[i] = ghostPlanes[mandPlane]->getMeshVectorFromLocal(normals[i]);
            verticies[i] = ghostPlanes[mandPlane]->getMeshCoordinatesFromLocal(verticies[i]);
        }
    }

    Q_EMIT sendFibulaToMesh(verticies, triangles, colours, normals, nbColours);
}

QString Viewer::helpString() const {
  QString text("<h2>H e l p  B o x</h2>");
  text += "This is the help text";

  return text;
}

void swap(int* a, int* b){
    int temp = *a;
    *a = *b;
    *b = temp;
}

int Viewer::partition(int sorted[], int start, int end){
    int p = sorted[end];
    int index = start - 1;

    for(int i=start; i<end; i++){
        double tangentAngleA = angle(curve->tangent(sorted[i]-1), curve->tangent(sorted[i]));
        double tangentAngleP = angle(curve->tangent(p-1), curve->tangent(p));

        if(tangentAngleA >= tangentAngleP){
            index++;
            swap(&sorted[index], &sorted[i]);
        }
    }
    swap(&sorted[index+1], &sorted[end]);
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
    if(isDrawMesh) isDrawMesh = false;
    else isDrawMesh = true;
    update();
}

void Viewer::initGhostPlanes(){
    for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];
    ghostPlanes.clear();

    int finalNb = nbGhostPlanes;

    int maxIndicies[nbGhostPlanes];

    const int startI = curve->indexForLength(curveIndexL, constraint);
    const int endI = curve->indexForLength(curveIndexR, -constraint);
    const int searchArea = endI - startI;

    // if there's enough space for a plane
    if(searchArea > 0){

        int sorted[searchArea];

        for(int i=0; i<searchArea; i++){
            sorted[i] = startI+i;
        }

        // Sort the indicies according to their tangent angles
        quicksort(sorted, 0, searchArea-1);

        maxIndicies[0] = sorted[0];
        int sortedIndex = 1;

        for(int i=1; i<nbGhostPlanes; i++){
            // the constraint (don't take it if it's too close to another existing plane)
            bool tooClose;
            do{
                tooClose = false;
                for(int j=i-1; j>=0; j--){
                    if(sortedIndex < searchArea && curve->discreteLength(maxIndicies[j],sorted[sortedIndex])<constraint){
                        tooClose = true;
                        break;
                    }
                }
                if(tooClose) sortedIndex++;
            }while(tooClose);

            if(sortedIndex >= searchArea){
                finalNb = i;
                break;
            }
            maxIndicies[i] = sorted[sortedIndex];
            sortedIndex++;  // move with i
        }

        // sort the planes
        for(int i=0; i<finalNb; i++){
            for(int j=i+1; j<finalNb; j++){
                if(maxIndicies[i] > maxIndicies[j]) swap(&maxIndicies[i], &maxIndicies[j]);
            }
        }

        ghostLocation.clear();
        for(int i=0; i<finalNb; i++) ghostLocation.push_back(maxIndicies[i]);

        // the number of ghost planes we can currently fit
        currentNbGhostPlanes = finalNb;

        addGhostPlanes(finalNb);

        for(unsigned int i=0; i<ghostPlanes.size(); i++) connect(&(ghostPlanes[i]->getCurvePoint()), &CurvePoint::curvePointTranslated, this, &Viewer::ghostPlaneMoved);

        // Update the fibula planes and polyline
        std::vector<Vec> angles = updatePolyline();

       // Q_EMIT haltMeshUpdate();

         // Send the info to the fibula
        double distance;
        if(finalNb > 0) distance = curve->discreteLength(curveIndexL, ghostLocation[0]);
        else distance = curve->discreteLength(curveIndexL, curveIndexR);
        Q_EMIT leftPosChanged(distance, angles);
        if(finalNb > 0) distance = curve->discreteLength(curveIndexR, ghostLocation[finalNb-1]);
        else distance = curve->discreteLength(curveIndexL, curveIndexR);
        Q_EMIT rightPosChanged(distance, angles);

       // Q_EMIT continueMeshUpdate();

    }
    else{
        std::vector<Vec> angles = updatePolyline();
        Q_EMIT ghostPlanesAdded(0,0,angles);
    }
}

void Viewer::cutMesh(){
    // Get the number of ghost planes from the total number of pieces dialog
    bool isNumberRecieved;
    int nbPieces = QInputDialog::getInt(this, "Cut mesh", "Number of pieces", 0, 1, 10, 1, &isNumberRecieved, Qt::WindowFlags());
    if(isNumberRecieved) nbGhostPlanes = nbPieces-1;
    else return;

    Q_EMIT preparingToCut();

    if(nbGhostPlanes==0) Q_EMIT noGhostPlanesToSend();

    // The dialog wasn't cancelled so the fibula can be cut
    Q_EMIT okToCut();

    mesh.setIsCut(Side::INTERIOR, true, true);
    isGhostPlanes = true;
    initGhostPlanes();
    update();
}

void Viewer::uncutMesh(){
    mesh.setIsCut(Side::INTERIOR, false, false);
    isGhostPlanes = false;
    for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];
    ghostPlanes.clear();
    //std::cout << "Ghost plane size : " << ghostPlanes.size() << std::endl;
    update();
}

// Slide the left plane
void Viewer::moveLeftPlane(int position){
    //bool isPassed = false;

    double percentage = static_cast<double>(position) / static_cast<double>(sliderMax);
    int index = static_cast<int>(percentage * static_cast<double>(nbU) );

    if( (curve->indexForLength(curveIndexR, -constraint) > index)){  // Only move if we're going backwards or we haven't met the other plane
        curveIndexL = index;

        if(curveIndexL >= nbU) curveIndexL = nbU-1;
        else if(curveIndexL < 0) curveIndexL = 0;   // shouldn't ever happen
    }
    else if( curveIndexL == curve->indexForLength(curveIndexR, -constraint) ) return;
    else{
        // isPassed = true;
        curveIndexL = curve->indexForLength(curveIndexR, -constraint);
    }

    leftPlane->setPosition(curve->getPoint(curveIndexL));
    //leftPlane->setOrientation(getNewOrientation(curveIndexL));
    matchPlaneToFrenet(leftPlane, curveIndexL);

    mesh.updatePlaneIntersections(leftPlane);

    if(isGhostPlanes && isGhostActive) handlePlaneMoveStart();

    double distance;

    // Get the info and send it to the fibula
    std::vector<Vec> angles = updatePolyline();

    update();

    if(!isGhostPlanes || !isGhostActive) distance = curve->discreteChordLength(curveIndexL, curveIndexR);
    else if(ghostPlanes.size()==0) distance = curve->discreteLength(curveIndexL, curveIndexR);  // is cut but no ghost planes
    else{
        distance = curve->discreteLength(curveIndexL, ghostLocation[0]);
    }

    Q_EMIT leftPosChanged(distance, angles);
    Q_EMIT setLRSliderValue(0);     // Reset the rotation slider
}

void Viewer::onLeftSliderReleased(){
    if(isGhostPlanes) handlePlaneMoveEnd();
    // creates an infinite loop if done when the slider value changes
    // TODO this means moveLeftPlane called twice
    Q_EMIT setLMSliderValue( static_cast<int>( (static_cast<double>(curveIndexL)/static_cast<double>(nbU)) * static_cast<double>(sliderMax) ) );
}

void Viewer::onRightSliderReleased(){
    if(isGhostPlanes) handlePlaneMoveEnd();
    // creates an infinite loop if done when the slider value changes
    // TODO this means moveLeftPlane called twice
    Q_EMIT setRMSliderValue( static_cast<int>( sliderMax - (static_cast<double>(curveIndexR)/static_cast<double>(nbU)) * static_cast<double>(sliderMax) ) );
}

void Viewer::rotateLeftPlane(int position){
    double percentage = static_cast<double>(position) / static_cast<double>(sliderMax);

    leftPlane->rotatePlaneXY(percentage);
    mesh.updatePlaneIntersections(leftPlane);
    update();
}

void Viewer::rotateRightPlane(int position){
    double percentage = static_cast<double>(position) / static_cast<double>(sliderMax);

    rightPlane->rotatePlaneXY(percentage);
    mesh.updatePlaneIntersections(rightPlane);
    update();
}

void Viewer::moveRightPlane(int position){

    double percentage = static_cast<double>(position) / static_cast<double>(sliderMax);
    int index = nbU - 1 - static_cast<int>(percentage * static_cast<double>(nbU) );

    if( index > curve->indexForLength(curveIndexL, constraint)){        // its within the correct boundaries
        curveIndexR = index;

        if(curveIndexR >= nbU) curveIndexR = nbU-1; // shouldn't ever happen either, outside of testing
        else if(curveIndexR < 0) curveIndexR = 0;   // shouldn't ever happen
    }
    else if(curveIndexR == curve->indexForLength(curveIndexL, constraint)) return;
    else curveIndexR = curve->indexForLength(curveIndexL, constraint);

    rightPlane->setPosition(curve->getPoint(curveIndexR));
    //rightPlane->setOrientation(getNewOrientation(curveIndexR));
    matchPlaneToFrenet(rightPlane, curveIndexR);

    mesh.updatePlaneIntersections(rightPlane);

    if(isGhostPlanes && isGhostActive) handlePlaneMoveStart();

    double distance;

    std::vector<Vec> angles = updatePolyline();

    if(!isGhostPlanes || !isGhostActive) distance = curve->discreteChordLength(curveIndexL, curveIndexR);
    else if(ghostPlanes.size()==0) distance = curve->discreteLength(curveIndexL, curveIndexR);  // is cut but no ghost planes
    else{
        distance = curve->discreteLength(ghostLocation[ghostPlanes.size()-1], curveIndexR);
    }

    update();

    Q_EMIT rightPosChanged(distance, angles);
    Q_EMIT setRRSliderValue(0); // Reset the rotation slider
}

void Viewer::handlePlaneMoveStart(){
    isGhostActive = false;      // disable drawing the polyline
    // Uncut the mesh but keep the ghostPlane marker as true
    mesh.setIsCut(Side::INTERIOR, false, false);
    for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];
    ghostPlanes.clear();

    Q_EMIT ghostPlaneMovementStart();
}

void Viewer::handlePlaneMoveEnd(){
    isGhostActive = true;       // enable the polyline
    // Recut
    Q_EMIT preparingToCut();
    if(nbGhostPlanes==0) Q_EMIT noGhostPlanesToSend();
    Q_EMIT okToCut();
    mesh.setIsCut(Side::INTERIOR, true, true);
    isGhostPlanes = true;
    initGhostPlanes();

    Q_EMIT ghostPlaneMovementEnd();
}

void Viewer::openOFF(QString filename) {
    std::vector<Vec3Df> &vertices = mesh.getVertices();
    std::vector<Triangle> &triangles = mesh.getTriangles();
    std::vector< std::vector<unsigned int>> &neighbours = mesh.getVertexNeighbours();
    std::vector< std::vector<unsigned int>> &vertexTriangles = mesh.getVertexTriangles();

    FileIO::openOFF(filename.toStdString(), vertices, triangles, neighbours, vertexTriangles);

    mesh.update();

    // Set the camera
    Vec3Df center;
    double radius;
    MeshTools::computeAveragePosAndRadius(vertices, center, radius);
    updateCamera(center, static_cast<float>(radius));

    update();
}

void Viewer::initCurve(){
    const long nbCP = 9;
    std::vector<Vec> control;

    control.push_back(Vec(-56.9335, -13.9973, 8.25454));

    control.push_back(Vec(-50.8191, -20.195, -19.53));
    control.push_back(Vec(-40.155, -34.5957, -50.7005));
    control.push_back(Vec(-27.6007, -69.2743, -67.6769));

    control.push_back(Vec(0, -85.966, -68.3154));

    control.push_back(Vec(26.7572, -69.0705, -65.6261));
    control.push_back(Vec(40.3576, -34.3609, -50.7634));
    control.push_back(Vec(46.2189, -21.3245, -17.9009));

    control.push_back(Vec(52.3669, -15.4613, 8.70223));


    curve = new Curve(nbCP, control);

    connect(curve, &Curve::curveReinitialised, this, &Viewer::updatePlanes);

    nbU = 200;
    curve->generateCatmull(nbU);

   initPlanes(Movable::STATIC);
}

void Viewer::initPlanes(Movable status){
    curveIndexR = nbU - 1;
    curveIndexL = 0;

    leftPlane = new Plane(40.0, status);
    rightPlane = new Plane(40.0, status);

    leftPlane->setPosition(curve->getPoint(curveIndexL));
    rightPlane->setPosition(curve->getPoint(curveIndexR));

    /*leftPlane->setOrientation(getNewOrientation(curveIndexL));
    rightPlane->setOrientation(getNewOrientation(curveIndexR));*/
    matchPlaneToFrenet(leftPlane, curveIndexL);
    matchPlaneToFrenet(rightPlane, curveIndexL);

    mesh.addPlane(leftPlane);
    mesh.addPlane(rightPlane);
}

void Viewer::addGhostPlanes(int nb){
    if(nb==0) return;
    double distances[nb+1];     // +1 for the last plane

    for(unsigned int i=0; i<static_cast<unsigned int>(nb); i++){
        //std::cout << "adding plane " << std::endl;
        ghostPlanes.push_back(new Plane(40.0, Movable::DYNAMIC));
        //std::cout << "plane added " << std::endl;
        //std::cout << "size after plane added : " << ghostPlanes.size() << std::endl;
        // ! HERE IS WHERE THE DESTRUCTOR CRASHES (in set Orientation)
        ghostPlanes[i]->setOrientation(getNewOrientation(ghostLocation[i]));
        ghostPlanes[i]->setPosition(curve->getPoint(ghostLocation[i]));
        if(i==0) distances[i] = curve->discreteLength(curveIndexL, ghostLocation[i]);
        else distances[i] = curve->discreteLength(ghostLocation[i-1], ghostLocation[i]);
    }

    distances[nb] = curve->discreteLength(ghostLocation[static_cast<unsigned int>(nb-1)], curveIndexR);

    std::vector<Vec> angles = updatePolyline();

    Q_EMIT ghostPlanesAdded(nb, distances, angles);
}

void Viewer::ghostPlaneMoved(){
    int nb = ghostPlanes.size();
    double distances[nb+1];     // +1 for the last plane

    for(unsigned int i=0; i<static_cast<unsigned int>(nb); i++){
        if(i==0) distances[i] = segmentLength(leftPlane->getPosition(), ghostPlanes[i]->getCurvePoint().getPoint());
        else distances[i] = segmentLength(ghostPlanes[i-1]->getCurvePoint().getPoint(), ghostPlanes[i]->getCurvePoint().getPoint());
    }

    distances[nb] = segmentLength(rightPlane->getPosition(), ghostPlanes[nb-1]->getCurvePoint().getPoint());

    std::vector<Vec> angles = updatePolyline();

    Q_EMIT ghostPlanesTranslated(nb, distances, angles);
}

void Viewer::updateCamera(const Vec3Df & center, float radius){
    camera()->setSceneCenter(Vec(static_cast<double>(center[0]), static_cast<double>(center[1]), static_cast<double>(center[2])));
    camera()->setSceneRadius(static_cast<double>(radius*1.05f));
    camera()->showEntireScene();
}

void Viewer::updatePlanes(){
    leftPlane->setPosition(curve->getPoint(curveIndexL));
    rightPlane->setPosition(curve->getPoint(curveIndexR));

    /*leftPlane->setOrientation(getNewOrientation(curveIndexL));
    rightPlane->setOrientation(getNewOrientation(curveIndexR));*/
    matchPlaneToFrenet(leftPlane, curveIndexL);
    matchPlaneToFrenet(rightPlane, curveIndexR);

    mesh.updatePlaneIntersections();

    update();
}

Quaternion Viewer::getNewOrientation(int index){
    Quaternion s = Quaternion(Vec(0,0,1.0), curve->tangent(index));
    return s.normalized();
}

std::vector<Vec> Viewer::getPolylinePlaneAngles(){
    // There are always 2 intersections for each polyline segment
    std::vector<Vec> angles;

    if(polyline.size()==0) return angles;

    angles.push_back(leftPlane->getPolylineVector(polyline[1]));

    for(unsigned int i=1; i<polyline.size()-1; i++){
        angles.push_back(ghostPlanes[i-1]->getPolylineVector(polyline[i-1]));
        angles.push_back(ghostPlanes[i-1]->getPolylineVector(polyline[i+1]));
    }

    angles.push_back(rightPlane->getPolylineVector(polyline[polyline.size()-2]));

    return angles;
}

double Viewer::angle(Vec a, Vec b){

    double na = a.norm();
    double nb = b.norm();

    double ab = a*b;

    return acos(ab / (na*nb));
}

double Viewer::segmentLength(const Vec a, const Vec b){
    return sqrt( pow((b.x - a.x), 2) + pow((b.y - a.y), 2) + pow((b.z - a.z), 2));
}

void Viewer::matchPlaneToFrenet(Plane *p, int index){
    Vec x, y, z;
    curve->getFrame(index,z,x,y);
    p->setFrameFromBasis(x,y,z);
}



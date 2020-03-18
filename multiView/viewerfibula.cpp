#include "viewerfibula.h"

ViewerFibula::ViewerFibula(QWidget *parent, StandardCamera *camera, int sliderMax, int fibulaOffset) : Viewer (parent, camera, sliderMax)
{
    indexOffset = 0;
    maxOffset = fibulaOffset;
    isPlanesRecieved = false;
    isCutSignal = false;
}

void ViewerFibula::initSignals(){
    connect(&mesh, &Mesh::sendInfoToManible, this, &ViewerFibula::recieveFromFibulaMesh);
}

void ViewerFibula::recieveFromFibulaMesh(std::vector<int> planes, std::vector<Vec> verticies, std::vector<std::vector<int>> triangles, std::vector<int> colours, std::vector<Vec> normals, int nbColours){
    std::vector<Vec> polylineInPlanes = getPolyline();

    Q_EMIT sendToManible(planes, verticies, triangles, polylineInPlanes, colours, normals, nbColours);

    polyline.clear();       // to stop it from being drawn
}

std::vector<Vec> ViewerFibula::getPolyline(){
    std::vector<Vec> polylineInPlanes;
    Vec v;

    createPolyline();

    v = leftPlane->getPolylineVector(polyline[1]);
    polylineInPlanes.push_back(v);

    for(unsigned int i=0; i<ghostPlanes.size(); i++){       // +1 offset for the left plane
        if(i%2==0) v = ghostPlanes[i]->getPolylineVector(polyline[i]);     // even: look behind
        else v = ghostPlanes[i]->getPolylineVector(polyline[i+2]);       // odd : look forward
        polylineInPlanes.push_back(v);
    }

    v = rightPlane->getPolylineVector(polyline[polyline.size()-2]);
    polylineInPlanes.push_back(v);

    return polylineInPlanes;
}

void ViewerFibula::createPolyline(){
    polyline.clear();

    polyline.push_back(leftPlane->getPosition());
    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        polyline.push_back(ghostPlanes[i]->getPosition());
    }

    polyline.push_back(rightPlane->getPosition());
}

void ViewerFibula::repositionPlanes(std::vector<Vec> polyline, std::vector<Vec> axes){
    if(isGhostPlanes){
        resetMandibleInfo(polyline, axes);
        repositionPlane(rightPlane, static_cast<unsigned int>(static_cast<int>(curveIndexR)+indexOffset));
        repositionPlane(leftPlane, static_cast<unsigned int>(static_cast<int>(curveIndexL)+indexOffset));
        //setPlanePositions();
        setPlaneOrientations();
   }
    else{
        repositionPlane(rightPlane, static_cast<unsigned int>(static_cast<int>(curveIndexR)+indexOffset));
        repositionPlane(leftPlane, static_cast<unsigned int>(static_cast<int>(curveIndexL)+indexOffset));
    }
    update();
}

void ViewerFibula::resetMandibleInfo(std::vector<Vec> polyline, std::vector<Vec> axes){
    mandiblePolyline.clear();
    mandiblePolyline = polyline;
    mandibleAxes.clear();
    mandibleAxes = axes;
}

void ViewerFibula::setPlanePositions(){
    leftPlane->setPosition(curve->getPoint( static_cast<unsigned int>(static_cast<int>(curveIndexL) + indexOffset)));
    rightPlane->setPosition(curve->getPoint(static_cast<unsigned int>(static_cast<int>(curveIndexR) + indexOffset)));
    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        ghostPlanes[i]->setPosition(curve->getPoint(static_cast<unsigned int>(static_cast<int>(ghostLocation[i]) + indexOffset)));
    }
}

void ViewerFibula::setPlaneOrientations(){
    if(mandiblePolyline.size()==0) return;
    Vec normal = Vec(0,0,1);

    // Initialise the ghost planes' rotation
    for(unsigned int i=0; i<ghostPlanes.size(); i++) repositionPlane(ghostPlanes[i], static_cast<unsigned int>(static_cast<int>(ghostLocation[i])+indexOffset));

    // Reset the rotation to line up with the fibula polyline
    std::vector<Vec> fibulaPolyline = getPolyline();
    Quaternion bLeft = Quaternion(normal, fibulaPolyline[0]);
    Quaternion bRight = Quaternion(-normal, fibulaPolyline[fibulaPolyline.size()-1]);
    leftPlane->rotate(bLeft);
    rightPlane->rotate(bRight);

    // Now we can move the normal to the mandible polyline from the fibula polyline
    Quaternion s = Quaternion(normal, mandiblePolyline[0]);
    leftPlane->rotate(s);
    s = Quaternion(-normal, mandiblePolyline[mandiblePolyline.size()-1]);
    rightPlane->rotate(s);

    // Orientate the ghost planes
    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        // To fibula
        Quaternion b = Quaternion(normal, fibulaPolyline[i+1]);
        ghostPlanes[i]->rotate(b);
        // To mandible
        if(i%2==0) b = Quaternion(-normal, mandiblePolyline[i+1]);       // the mandible polyline is in relation to the forward facing plane
        else b = Quaternion(normal, mandiblePolyline[i+1]);
        ghostPlanes[i]->rotate(b);
    }

    fibulaPolyline = getPolyline();
    swivelToPolyline(fibulaPolyline);

    Q_EMIT requestAxes();

    /*std::cout << "Distances : " << std::endl;
    for(unsigned int i = 0; i<distances.size(); i++){
        std::cout << i << " : " << distances[i] << std::endl;
    }

    std::cout << "Approaching planes :" << std::endl;*/
    for(unsigned int i=0; i<ghostPlanes.size(); i+=2){
        approachPlanes(i);
    }

    /*std::cout << "ACTUAL DISTANCES : " << std::endl;
    std::cout << "0 : " << curve->discreteLength(curveIndexL, ghostLocation[0]) << std::endl;
    for(unsigned int i=0; i<ghostPlanes.size()-1; i++) {
        std::cout << i+1 << " : " <<  curve->discreteLength(ghostLocation[i], ghostLocation[i+1]) << std::endl;
    }
    std::cout << "last : " << curve->discreteLength(ghostLocation[ghostPlanes.size()-1], curveIndexR) << std::endl;*/

    /*std::cout << "Checking intersections" << std::endl;
    for(unsigned int i=0; i<ghostLocation.size()-1; i+=2){
        Vec v0, v1, v2, v3;
        ghostPlanes[i+1]->getCorners(v0, v1, v2, v3);
        std::cout << "Intersection? : " << ghostPlanes[i]->isIntersectionPlane(v0, v1, v2, v3) << std::endl;
    }*/
}

void ViewerFibula::swivelToPolyline(std::vector<Vec>& fibulaPolyline){
    Vec axis = Vec(0,0,1);
    /*leftPlane->rotatePlane(axis, M_PI*2.0);
    rightPlane->rotatePlane(axis, M_PI*2.0);*/

    mandiblePolyline[0].normalize();
    fibulaPolyline[0].normalize();
    Vec mandPoint = leftPlane->getLocalProjection(mandiblePolyline[0]);
    Vec fibPoint = leftPlane->getLocalProjection(fibulaPolyline[0]);
    mandPoint.normalize();
    fibPoint.normalize();
    double alpha = angle(mandPoint, fibPoint) + M_PI;
    leftPlane->rotatePlane(axis, alpha);

    if(ghostPlanes.size()!=0) {     // NOTE : always seems to be a rotation of 180 degrees
        //ghostPlanes[0]->rotatePlane(axis, M_PI*2.0);
        for(unsigned int i=1; i<ghostPlanes.size()-2; i+=2){
            mandiblePolyline[i+1].normalize();
            fibulaPolyline[i+1].normalize();
            Vec mandPoint = ghostPlanes[i]->getLocalProjection(mandiblePolyline[i+1]);
            Vec fibPoint = ghostPlanes[i]->getLocalProjection(fibulaPolyline[i+1]);
            mandPoint.normalize();
            fibPoint.normalize();
            double alpha = angle(mandPoint, fibPoint) + M_PI;
            ghostPlanes[i]->rotatePlane(axis, alpha);
            ghostPlanes[i+1]->rotatePlane(axis, alpha);
        }
        //ghostPlanes[ghostPlanes.size()-1]->rotatePlane(axis, M_PI*2.0);
    }

    unsigned int polySize = mandiblePolyline.size()-1;
    mandiblePolyline[polySize].normalize();
    fibulaPolyline[polySize].normalize();
    mandPoint = rightPlane->getLocalProjection(mandiblePolyline[polySize]);
    fibPoint = rightPlane->getLocalProjection(fibulaPolyline[polySize]);
    mandPoint.normalize();
    fibPoint.normalize();
    alpha = angle(mandPoint, fibPoint) + M_PI;
    rightPlane->rotatePlane(axis, alpha);
}

// Rotate the end plane to match the mandibule
void ViewerFibula::recieveAxes(std::vector<Vec> axes){

    if(ghostPlanes.size()==0){
        leftPlane->setOrientationFromOtherReference(axes, 0, rightPlane);
    }
    else{
        ghostPlanes[0]->setOrientationFromOtherReference(axes, 0, leftPlane);

        unsigned int axesIndex = 3;
        for(unsigned int i=2; i<ghostPlanes.size()-1; i+=2){
            ghostPlanes[i]->setOrientationFromOtherReference(axes, axesIndex, ghostPlanes[i-1]);
            axesIndex+=3;
        }

        unsigned int lastIndex = static_cast<unsigned int>(ghostPlanes.size()-1);
        unsigned int lastAxe = static_cast<unsigned int>(axes.size()-3);
        ghostPlanes[lastIndex]->setOrientationFromOtherReference(axes, lastAxe, rightPlane);
    }
}

// Move all planes by the same offset (right plane INCLUDED) - when the slider is dragged
void ViewerFibula::movePlanes(int position){
    int offset = static_cast<int>(static_cast<double>(position)/ static_cast<double>(maxOffset) * static_cast<double>(nbU));

    // Check that it this offset doesn't exceed the size of the fibula
    if(static_cast<int>(curveIndexL) + offset < static_cast<int>(nbU) && static_cast<int>(curveIndexL) + offset > 0 && static_cast<int>(curveIndexR) + offset < static_cast<int>(nbU) && static_cast<int>(curveIndexR) + offset > 0){
        indexOffset = offset;
        repositionPlanes(mandiblePolyline, mandibleAxes);
    }

    mesh.setTransfer(false);
    mesh.updatePlaneIntersections();
}

void ViewerFibula::planesMoved(){
    mesh.setTransfer(true);
    mesh.sendToManible();
}

// Add the ghost planes (this should only be called once)
void ViewerFibula::addGhostPlanes(unsigned int nb){
    for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];     // remove any ghost planes
    ghostPlanes.clear();
    Vec pos = Vec(0,0,0);

    for(unsigned int i=0; i<static_cast<unsigned int>(nb); i++){
        ghostPlanes.push_back(new Plane(25.0, Movable::STATIC, pos));

        // If we're too far along the fibula, take it all back
        int overload = static_cast<int>(ghostLocation[i]) + indexOffset - static_cast<int>(curve->getNbU()) + 1;   // The amount by which the actual index passes the end of the curve
        if(overload > 0) indexOffset -= overload;
    }
}

// Find the locations of the ghost planes from the distances from the planes in the mandible
void ViewerFibula::findGhostLocations(unsigned int nb, double distance[]){
    distances.clear();
    for(unsigned int i=0; i<nb; i++) {
        distances.push_back(distance[i]);
        distances.push_back(securityMargin);
    }

    distances.push_back(distance[nb]);

    findIndexesFromDistances();
}

void ViewerFibula::findIndexesFromDistances(){
    ghostLocation.clear();
    unsigned int nb = static_cast<unsigned int>(distances.size()-1);

    unsigned int index = curve->indexForLength(curveIndexL+indexOffset, distances[0]);
    ghostLocation.push_back(index);

    for(unsigned int i=1; i<nb; i++){
        index = curve->indexForLength(ghostLocation[i-1], distances[i]);
        ghostLocation.push_back(index);
    }

    curveIndexR = curve->indexForLength(ghostLocation[nb-1], distances[nb]);        // place the right plane after the last ghost plane (left plane doesn't move)
}

// NOT USED
void ViewerFibula::matchToMandibleFrame(Plane* p1, Plane* p2, Vec a, Vec b, Vec c, Vec x, Vec y, Vec z){
    p1->setFrameFromBasis(a,b,c);
    p2->setFrameFromBasis(x,y,z);
}

// Don't wait for ghost planes, go ahead and cut
void ViewerFibula::noGhostPlanesToRecieve(std::vector<Vec> mandPolyline, std::vector<Vec> axes){
    isPlanesRecieved = true;
    for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];
    ghostPlanes.clear();
    mesh.deleteGhostPlanes();
    isGhostPlanes = true;
    repositionPlanes(mandPolyline, axes);
    handleCut();
}

// Add ghost planes that correspond to the ghost planes in the jaw
void ViewerFibula::ghostPlanesRecieved(unsigned int nb, double distance[], std::vector<Vec> mandPolyline, std::vector<Vec> axes){
    if(nb==0) return;

    findGhostLocations(nb, distance);
    addGhostPlanes(2* static_cast<int>(nb));    // 2*nb ghost planes : there are 2 angles for each plane in the manible, so twice the number of ghost planes

    repositionPlanes(mandPolyline, axes);

    isPlanesRecieved = true;
    handleCut();
}

// When we want to move the right plane
void ViewerFibula::movePlaneDistance(double distance, std::vector<Vec> mandPolyline, std::vector<Vec> axes){
    unsigned int newIndex;

    if(ghostPlanes.size()==0) newIndex = curve->indexForLength(curveIndexL, distance);
    else newIndex = curve->indexForLength(ghostLocation[ghostPlanes.size()-1], distance);

    if(static_cast<unsigned int>(static_cast<int>(newIndex) + indexOffset) >= nbU) return;      // This should never happen
    else curveIndexR = newIndex;

    repositionPlanes(mandPolyline, axes);

    mesh.updatePlaneIntersections(rightPlane);
}

// One of the ghost planes is moved in the jaw
// Probably out of date
void ViewerFibula::middlePlaneMoved(unsigned int nb, double distances[], std::vector<Vec> mandPolyline, std::vector<Vec> axes){
    if(nb==0) return;

    findGhostLocations(nb, distances);

    // If we're too far along the fibula, take it all back
    int overload = static_cast<int>(curveIndexR) + indexOffset - static_cast<int>(curve->getNbU()) + 1;   // The amount by which the actual index passes the end of the curve
    if(overload > 0){
        indexOffset -= overload;
        Q_EMIT setPlaneSliderValue(static_cast<int>( (static_cast<double>(indexOffset)/static_cast<double>(nbU)) * static_cast<double>(maxOffset) ));
    }

    repositionPlanes(mandPolyline, axes);

    mesh.updatePlaneIntersections(rightPlane);
}

// Initialise the curve that the planes follow (to eventually be changed to automatically calculate the points)
void ViewerFibula::initCurve(){
    const long nbCP = 6;
    std::vector<Vec> control;

    control.push_back(Vec(108.241, 69.6891, -804.132));
    control.push_back(Vec(97.122, 82.1788, -866.868));
    control.push_back(Vec(93.5364, 90.1045, -956.126));
    control.push_back(Vec(83.3966, 92.5807, -1069.7));
    control.push_back(Vec(80.9, 90.1, -1155));
    control.push_back(Vec(86.4811, 90.9929, -1199.7));

    curve = new Curve(nbCP, control);

    nbU = 2000;

    int nbSeg = nbCP-3;
    nbU -= static_cast<unsigned int>(static_cast<int>(nbU)%nbSeg);

    curve->generateCatmull(nbU);
    connect(curve, &Curve::curveReinitialised, this, &Viewer::updatePlanes);

    initPlanes(Movable::STATIC);
}

void ViewerFibula::cutMesh(){
    isCutSignal = true;
    handleCut();
}

void ViewerFibula::handleCut(){
    if(isCutSignal && isPlanesRecieved){
        for(unsigned int i=0; i<ghostPlanes.size(); i++){
            mesh.addPlane(ghostPlanes[i]);
        }

        if(ghostPlanes.size()==0) mesh.setIsCut(Side::EXTERIOR, true, true);    // call the update if an exterior plane isn't going to
        else mesh.setIsCut(Side::EXTERIOR, true, false);

        isGhostPlanes = true;
        isCutSignal = false;
        isPlanesRecieved = false;
    }
}

void ViewerFibula::uncutMesh(){
    isPlanesRecieved = false;
    mesh.setIsCut(Side::EXTERIOR, false, false);
    isGhostPlanes = false;
    for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];
    ghostPlanes.clear();
    repositionPlanes(mandiblePolyline, mandibleAxes);
    update();
}

void ViewerFibula::findClosestPoint(unsigned int pNb, Vec &a, Vec &b){
    Vec pos(0,0,0);
    Plane tempPlane(40.0, Movable::STATIC, pos);
    repositionPlane(&tempPlane, ghostLocation[pNb]);

    const std::vector<unsigned int> tInd1 = mesh.getIntersectionTriangles(pNb+2);
    a = findMaxZ(tInd1, tempPlane);
    const std::vector<unsigned int> tInd2 = mesh.getIntersectionTriangles(pNb+3);
    b = findMinZ(tInd2, tempPlane);
}

Vec ViewerFibula::findMinZ(const std::vector<unsigned int> &tIndexes, Plane &tempPlane){
    Vec minI(0,0,0);
    double min = DBL_MAX;
    for(unsigned int i=0; i<tIndexes.size(); i++){ // for each triangle aligned with the plane
        for(unsigned int j=0; j<3; j++){
            unsigned int v = mesh.getTriangle(tIndexes[i]).getVertex(j);
            Vec vVec = Vec(mesh.getVertex(v));
            double z = tempPlane.getLocalCoordinates(vVec).z;
            if(z < min){
                min = z;
                minI = vVec;
            }
        }
    }

    return minI;
}

Vec ViewerFibula::findMaxZ(const std::vector<unsigned int> &tIndexes, Plane &tempPlane){
    Vec maxI(0,0,0);
    double max = -DBL_MAX;
    for(unsigned int i=0; i<tIndexes.size(); i++){ // for each triangle aligned with the plane
        for(unsigned int j=0; j<3; j++){
            unsigned int vI = mesh.getTriangle(tIndexes[i]).getVertex(j);
            Vec vVec = Vec(mesh.getVertex(vI));
            double z = tempPlane.getLocalCoordinates(vVec)[0];
            if(z > max){
                max = z;
                maxI = vVec;
            }
        }
    }

    return maxI;
}

void ViewerFibula::approachPlanes(unsigned int pStart){
    Vec p1, p2;
    findClosestPoint(pStart, p1, p2);
    double distZ = p2.z - p1.z;

    Vec pB1, pB2;
    pB1 = curve->getPoint(ghostLocation[pStart]);
    pB2 = curve->getPoint(ghostLocation[pStart+1]);
    double currentDistZ = pB2.z - pB1.z;

    double distPercentage = distZ/currentDistZ;
    double distShift = euclideanDistance(pB1, pB2) * distPercentage;

    /*for(unsigned int i=pStart+1; i<ghostPlanes.size(); i++) ghostLocation[i] -= distShift;
    curveIndexR -= distShift;*/

    /*std::cout << "Points : " << p1.x << " , " << p1.y << " , " << p1.z << " and " << p2.x << " , " << p2.y << " , " << p2.z << std::endl;
    std::cout << "Bases : " << pB1.x << " , " << pB1.y << " , " << pB1.z << " and " << pB2.x << " , " << pB2.y << " , " << pB2.z << std::endl;
    std::cout << "Current distance : " << pB2.z - pB1.z << std::endl;
    std::cout << "Distance " << pStart << " and " << pStart+1 << " : " << dist << std::endl;*/
}

double ViewerFibula::euclideanDistance(Vec &a, Vec &b){
    return sqrt( pow(a.x-b.x, 2.0) + pow(a.y-b.y, 2.0) + pow(a.z-b.z, 2.0));
}

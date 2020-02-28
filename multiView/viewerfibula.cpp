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

    //polyline.clear();       // to stop it from being drawn
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
    if(polyline.size()!=0){
        resetMandibleInfo(polyline, axes);
        setPlanePositions();
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

    for(unsigned int i=0; i<ghostPlanes.size(); i++) repositionPlane(ghostPlanes[i], static_cast<unsigned int>(static_cast<int>(ghostLocation[i])+indexOffset));

    std::vector<Vec> fibulaPolyline = getPolyline();
    Quaternion bLeft = Quaternion(Vec(0,0,1), fibulaPolyline[0]);
    Quaternion bRight = Quaternion(-Vec(0,0,1), fibulaPolyline[fibulaPolyline.size()-1]);
    leftPlane->rotate(bLeft.normalized());
    rightPlane->rotate(bRight.normalized());

    Vec normal = leftPlane->getNormal();        // The normal defines the polyline, so move our polyline to the mandible polyline
    Quaternion s = Quaternion(normal, mandiblePolyline[0]);  // -normal so it doesnt do a 180 flip (a rotation of the normal to the polyline)
    leftPlane->rotate(s.normalized());

    // Orientate the ghost planes
    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        Quaternion b = Quaternion(Vec(0,0,1), fibulaPolyline[i+1]);
        ghostPlanes[i]->rotate(b);
        // To mandible
        if(i%2==0) b = Quaternion(-Vec(0,0,1), mandiblePolyline[i+1]);       // the mandible polyline is in relation to the forward facing plane
        else b = Quaternion(Vec(0,0,1), mandiblePolyline[i+1]);
        ghostPlanes[i]->rotate(b);
    }

    // Orientate the right plane
    normal = rightPlane->getNormal();
    s = Quaternion(-normal, mandiblePolyline[mandiblePolyline.size()-1]);
    rightPlane->rotate(s.normalized());

    swivelToPolyline();

    Q_EMIT requestAxes();
}

double calcNorm(Vec &v){
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

void ViewerFibula::swivelToPolyline(){
    Vec axis = Vec(0,0,1);
    leftPlane->rotatePlane(axis, M_PI*2.0);
    rightPlane->rotatePlane(axis, M_PI*2.0);

    if(ghostPlanes.size()!=0) {
        std::vector<Vec> fibulaPolyline = getPolyline();
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
    }
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
void ViewerFibula::addGhostPlanes(int nb){
    for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];     // remove any ghost planes
    ghostPlanes.clear();

    for(unsigned int i=0; i<static_cast<unsigned int>(nb); i++){
        ghostPlanes.push_back(new Plane(25.0, Movable::STATIC));

        // If we're too far along the fibula, take it all back
        int overload = static_cast<int>(ghostLocation[i]) + indexOffset - static_cast<int>(curve->getNbU()) + 1;   // The amount by which the actual index passes the end of the curve
        if(overload > 0) indexOffset -= overload;
    }
}

// Find the locations of the ghost planes from the distances from the planes in the mandible
void ViewerFibula::findGhostLocations(unsigned int nb, double distance[]){
    distances.clear();
    for(unsigned int i=0; i<=nb; i++) distances.push_back(distance[i]);

    findIndexesFromDistances();
}

void ViewerFibula::findIndexesFromDistances(){
    ghostLocation.clear();
    unsigned int index = curve->indexForLength(curveIndexL+indexOffset, distances[0]);
    ghostLocation.push_back(index);
    std::cout << "Ghost location 0 : " << ghostLocation[0] << std::endl;
    unsigned int nextIndex = curve->indexForLength(index, securityMargin);
    ghostLocation.push_back(nextIndex);     // the mirror plane
    std::cout << "Ghost location 1 : " << ghostLocation[1] << std::endl;
    unsigned int nb = static_cast<unsigned int>(distances.size()-1);

    for(unsigned int i=1; i<nb; i++){
        index = curve->indexForLength(ghostLocation[2*i-1], distances[i]);
        ghostLocation.push_back(index);
        std::cout << "Ghost location " << 2*i << " : " << ghostLocation[2*i] << std::endl;
        unsigned int nbU = curve->getNbU();
        nextIndex = curve->indexForLength(index, securityMargin);
        if((nextIndex)<nbU) ghostLocation.push_back(nextIndex);
        else ghostLocation.push_back(nbU-1);
        std::cout << "Ghost location " << 2*i+1 << " : " << ghostLocation[2*i+1] << std::endl;
    }
    std::cout << "Nb ghost locations found : " << ghostLocation.size() << std::endl;
    curveIndexR = curve->indexForLength(ghostLocation[2*nb-1], distances[nb]);        // place the right plane after the last ghost plane (left plane doesn't move)

}

// NOT USED
void ViewerFibula::matchToMandibleFrame(Plane* p1, Plane* p2, Vec a, Vec b, Vec c, Vec x, Vec y, Vec z){
    p1->setFrameFromBasis(a,b,c);
    p2->setFrameFromBasis(x,y,z);
}

// Don't wait for ghost planes, go ahead and cut
void ViewerFibula::noGhostPlanesToRecieve(){
    isPlanesRecieved = true;
    handleCut();
}

// Add ghost planes that correspond to the ghost planes in the jaw
void ViewerFibula::ghostPlanesRecieved(unsigned int nb, double distance[], std::vector<Vec> mandPolyline, std::vector<Vec> axes){
    std::cout << "Nb ghost planes recieved : " << nb << std::endl;
    if(nb==0){      // if no ghost planes were actually recieved
        for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];
        ghostPlanes.clear();
        mesh.deleteGhostPlanes();
        return;
    }

    findGhostLocations(nb, distance);
    addGhostPlanes(2* static_cast<int>(nb));    // 2*nb ghost planes : there are 2 angles for each plane in the manible, so twice the number of ghost planes

    repositionPlanes(mandPolyline, axes);

    std::cout << "  LOCATIONS : " << std::endl;
    Vec l = leftPlane->getPosition();
    std::cout << "      Left : " << l.x << " , " << l.y << " , " << l.z << std::endl;
    for(unsigned int i=0; i<ghostLocation.size(); i++){
        l = ghostPlanes[i]->getPosition();
        std::cout << "      " << i << " : " << l.x << " , " << l.y << " , " << l.z << std::endl;
    }
    l = rightPlane->getPosition();
    std::cout << "      Right : " << l.x << " , " << l.y << " , " << l.z << std::endl;

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
// NOT USED FOR NOW
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

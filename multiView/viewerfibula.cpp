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
    std::vector<Vec> polylineInPlanes;
    Vec v;

    createPolyline();

    // Get the polyline vector in relation to the planes (in order of the planes)
    v = leftPlane->getPolylineVector(polyline[1]);
    polylineInPlanes.push_back(v);

    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        // +1 offset for the left plane
        if(i%2==0) v = ghostPlanes[i].getPolylineVector(polyline[i]);     // even: look behind
        else v = ghostPlanes[i].getPolylineVector(polyline[i+2]);       // odd : look forward
        polylineInPlanes.push_back(v);
    }

    v = rightPlane->getPolylineVector(polyline[polyline.size()-2]);
    polylineInPlanes.push_back(v);

    Q_EMIT sendToManible(planes, verticies, triangles, polylineInPlanes, colours, normals, nbColours);

    polyline.clear();       // to stop it from being drawn
}

void ViewerFibula::createPolyline(){
    polyline.clear();

    polyline.push_back(leftPlane->getPosition());
    for(int i=0; i<ghostPlanes.size(); i++) polyline.push_back(ghostPlanes[i].getPosition());
    polyline.push_back(rightPlane->getPosition());
}

// Move all planes by the same offset (right plane INCLUDED) - when the slider is dragged
void ViewerFibula::movePlanes(int position){

    int offset = static_cast<int>(static_cast<double>(position)/ static_cast<double>(maxOffset) * static_cast<double>(*nbU));

    // Check that it this offset doesn't exceed the size of the fibula
    if(curveIndexL + offset < *nbU && curveIndexL + offset > 0 && curveIndexR + offset < *nbU && curveIndexR + offset > 0){
        indexOffset = offset;

        // Reset the position and orientations of ALL planes
        leftPlane->setPosition(curve->getPoint(curveIndexL + indexOffset));
        rightPlane->setPosition(curve->getPoint(curveIndexR + indexOffset));
        leftPlane->setOrientation(getNewOrientation(curveIndexL + indexOffset));
        rightPlane->setOrientation(getNewOrientation(curveIndexR + indexOffset));

        for(unsigned int i=0; i<ghostPlanes.size(); i++){
            ghostPlanes[i].setPosition(curve->getPoint(ghostLocation[i] + indexOffset));
            ghostPlanes[i].setOrientation(getNewOrientation(ghostLocation[i] + indexOffset));
        }
    }

    setPlaneOrientations(angleVectors);

    mesh.setTransfer(false);
    mesh.updatePlaneIntersections();

    update();

}

void ViewerFibula::planesMoved(){
    mesh.setTransfer(true);
    mesh.sendToManible();
    update();
}

// Add the ghost planes (this should only be called once)
void ViewerFibula::addGhostPlanes(int nb){
    ghostPlanes.clear();

    for(unsigned int i=0; i<static_cast<unsigned int>(nb); i++){
        ghostPlanes.push_back(Plane(25.0, Movable::STATIC));
        int index = ghostLocation[i];

        // If we're too far along the fibula, take it all back
        int overload = index + indexOffset - curve->getNbU() + 1;   // The amount by which the actual index passes the end of the curve
        if(overload > 0){
            indexOffset -= overload;
            reinitialisePlanes(i+1);
            // Q_EMIT setPlaneSliderValue(static_cast<int>( (static_cast<double>(indexOffset)/static_cast<double>(*nbU)) * static_cast<double>(maxOffset) ));
        }
        ghostPlanes[i].setPosition(curve->getCurve()[index + indexOffset]);
        ghostPlanes[i].setOrientation(getNewOrientation(index + indexOffset));
    }

    update();
}

// Find the locations of the ghost planes from the distances from the planes in the mandible
void ViewerFibula::findGhostLocations(int nb, double distance[]){
    ghostLocation.clear();

    int index = curve->indexForLength(curveIndexL, distance[0]);
    ghostLocation.push_back(index);
    ghostLocation.push_back(index+25);      // NOTE 25 is just a temporary value
    for(unsigned int i=1; i<static_cast<unsigned int>(nb); i++){
        index = curve->indexForLength(ghostLocation[2*i-1], distance[i]);
        ghostLocation.push_back(index);
        int nbU = curve->getNbU();
        if((index+25)<nbU) ghostLocation.push_back(index+25);
        else ghostLocation.push_back(nbU-1);
    }
    curveIndexR = curve->indexForLength(ghostLocation[2*static_cast<unsigned int>(nb)-1], distance[nb]);
}

// Re-orientate the planes to correspond to the same angles as the jaw
void ViewerFibula::setPlaneOrientations(std::vector<Vec> angles){
    if(angles.size()==0) return;

    angleVectors.clear();
    angleVectors = angles;

    // Orientate the left plane
    Vec normal = leftPlane->getNormal();
    Quaternion s = Quaternion(-normal, angles[0]);  // -normal so it doesnt do a 180 flip
    leftPlane->setOrientation(s.normalized());

    // Orientate the ghost planes
    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        normal = ghostPlanes[i].getNormal();
        s = Quaternion(normal, angles[i]);
        ghostPlanes[i].setOrientation(s.normalized());
    }

    // Orientate the right plane
    normal = rightPlane->getNormal();
    s = Quaternion(normal, angles[angles.size()-1]);
    rightPlane->setOrientation(s.normalized());
}

// Don't wait for ghost planes, go ahead and cut
void ViewerFibula::noGhostPlanesToRecieve(){
    isPlanesRecieved = true;
    handleCut();
}

// Add ghost planes that correspond to the ghost planes in the jaw
void ViewerFibula::ghostPlanesRecieved(int nb, double distance[], std::vector<Vec> angles){
    // if no ghost planes were actually recieved
    if(nb==0){
        ghostPlanes.clear();        // TODO look at this (call noGhostPlanesToRecieve?)
        mesh.deleteGhostPlanes();
        return;
    }

    int oldNb = ghostPlanes.size() / 2;

    findGhostLocations(nb, distance);

    // doesn't work if its done before the curve is initialised (should never happen)
    // 2*nb ghost planes : there are 2 angles for each plane in the manible, so twice the number of ghost planes
    addGhostPlanes(2*nb);

    // Once everything is initialised, adjust the rotation
    setPlaneOrientations(angles);

    // If its cut and the number of planes has changed
    if(mesh.getIsCut() && nb!=oldNb){
        mesh.deleteGhostPlanes();
        cutMesh();
    }

    isPlanesRecieved = true;
    handleCut();
}

// When we want to move the right plane (the right plane is moved in the jaw)
void ViewerFibula::movePlaneDistance(double distance, std::vector<Vec> angles){
    int newIndex;

    if(ghostPlanes.size()==0) newIndex = curve->indexForLength(curveIndexL, distance);
    else newIndex = curve->indexForLength(ghostLocation[ghostPlanes.size()-1], distance);

    if(newIndex + indexOffset >= *nbU) return;      // This should never happen
    else curveIndexR = newIndex;

    rightPlane->setPosition(curve->getCurve()[curveIndexR + indexOffset]);
    rightPlane->setOrientation(getNewOrientation(curveIndexR + indexOffset));   // initial orientation

    angleVectors.clear();
    angleVectors = angles;

    setPlaneOrientations(angles);

    mesh.updatePlaneIntersections(rightPlane);
    update();
}

// When we want to move all the planes (except the static left) by a distance (the left plane is moved in the jaw)
void ViewerFibula::moveGhostPlaneDistance(double distance, std::vector<Vec> angles){
    int offset = 0;
    if(ghostPlanes.size()==0) movePlaneDistance(distance, angles);
    else offset = curve->indexForLength(curveIndexL, distance) - ghostLocation[0];

    // move all planes by the offset
    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        ghostLocation[i] += offset;
        ghostPlanes[i].setPosition(curve->getCurve()[ghostLocation[i] + indexOffset]);
        ghostPlanes[i].setOrientation(getNewOrientation(ghostLocation[i] + indexOffset));
    }

    // Gives a problem in extreme cases, so we check the index is still within bounds
    curveIndexR += offset;

    // If we're too far along the fibula, take it all back
    int overload = curveIndexR + indexOffset - curve->getNbU() + 1;   // The amount by which the actual index passes the end of the curve
    if(overload > 0){
        indexOffset -= overload;
        reinitialisePlanes(ghostPlanes.size()+1);
        Q_EMIT setPlaneSliderValue(static_cast<int>( (static_cast<double>(indexOffset)/static_cast<double>(*nbU)) * static_cast<double>(maxOffset) ));
    }

    rightPlane->setPosition(curve->getCurve()[curveIndexR + indexOffset]);
    rightPlane->setOrientation(getNewOrientation(curveIndexR + indexOffset));   // initial orientation

    angleVectors.clear();
    angleVectors = angles;

    setPlaneOrientations(angles);

    mesh.updatePlaneIntersections(rightPlane);
    update();
}

// One of the ghost planes is moved in the jaw
void ViewerFibula::middlePlaneMoved(int nb, double distances[], std::vector<Vec> angles){
    if(nb==0) return;

    findGhostLocations(nb, distances);

    // update the ghost planes
    for(unsigned int i=0; i<static_cast<unsigned int>(2*nb); i++){
        ghostPlanes[i].setPosition(curve->getCurve()[ghostLocation[i] + indexOffset]);
        ghostPlanes[i].setOrientation(getNewOrientation(ghostLocation[i] + indexOffset));
    }

    // If we're too far along the fibula, take it all back
    int overload = curveIndexR + indexOffset - curve->getNbU() + 1;   // The amount by which the actual index passes the end of the curve
    if(overload > 0){
        indexOffset -= overload;
        reinitialisePlanes(ghostPlanes.size()+1); // reinit everything but the right plane
        Q_EMIT setPlaneSliderValue(static_cast<int>( (static_cast<double>(indexOffset)/static_cast<double>(*nbU)) * static_cast<double>(maxOffset) ));
    }

    // update the right plane
    rightPlane->setPosition(curve->getCurve()[curveIndexR + indexOffset]);
    rightPlane->setOrientation(getNewOrientation(curveIndexR + indexOffset));

    angleVectors.clear();
    angleVectors = angles;

    setPlaneOrientations(angles);

    mesh.updatePlaneIntersections(rightPlane);

    update();
}

void ViewerFibula::reinitialisePlanes(unsigned int nbToInit){
    if(nbToInit==0) return;

    //Move the left plane and all the already initialised ghost planes back
    leftPlane->setPosition(curve->getCurve()[curveIndexL + indexOffset]);
    leftPlane->setOrientation(getNewOrientation(curveIndexL + indexOffset));

    if(nbToInit>=ghostPlanes.size()+2){
        rightPlane->setPosition(curve->getCurve()[curveIndexR + indexOffset]);
        rightPlane->setOrientation(getNewOrientation(curveIndexR + indexOffset));
        nbToInit = ghostPlanes.size() + 1;
    }

    for(unsigned int j=0; j<nbToInit-1; j++){
        ghostPlanes[j].setPosition(curve->getCurve()[ghostLocation[j] + indexOffset]);
        ghostPlanes[j].setOrientation(getNewOrientation(ghostLocation[j] + indexOffset));
    }
}

// Initialise the curve that the planes follow (to eventually be changed to automatically calculate the points)
void ViewerFibula::initCurve(){
    const long nbCP = 6;
    ControlPoint* control[nbCP];

    control[0] = new ControlPoint(108.241, 69.6891, -804.132);
    control[1] = new ControlPoint(97.122, 82.1788, -866.868);
    control[2] = new ControlPoint(93.5364, 90.1045, -956.126);
    control[3] = new ControlPoint(83.3966, 92.5807, -1069.7);
    control[4] = new ControlPoint(80.9, 90.1, -1155);
    control[5] = new ControlPoint(86.4811, 90.9929, -1199.7);

    curve = new Curve(nbCP, control);

    *nbU = 300;

    int nbSeg = nbCP-3;
    nbU -= *nbU%nbSeg;

    curve->generateCatmull(nbU);
    //curve->generateBSpline(nbU, 3);
    connect(curve, &Curve::curveReinitialised, this, &Viewer::updatePlanes);

    initPlanes(Movable::STATIC);
}

void ViewerFibula::cutMesh(){
    isCutSignal = true;
    handleCut();
    //update();
}

void ViewerFibula::handleCut(){
    if(isCutSignal && isPlanesRecieved){

        for(unsigned int i=0; i<ghostPlanes.size(); i++){
            mesh.addPlane(&ghostPlanes[i]);
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
    ghostPlanes.clear();        // NOTE To eventually be changed
    // Reset their orientations
    leftPlane->setOrientation(getNewOrientation(curveIndexL));
    rightPlane->setOrientation(getNewOrientation(curveIndexR));
    update();
}

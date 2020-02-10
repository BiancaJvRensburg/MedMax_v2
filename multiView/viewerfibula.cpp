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

    // Get the polyline vector in relation to the planes (in order of the planes)
    v = leftPlane->getPolylineVector(polyline[1]);
    polylineInPlanes.push_back(v);

    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        // +1 offset for the left plane
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
    for(unsigned int i=0; i<ghostPlanes.size(); i++) polyline.push_back(ghostPlanes[i]->getPosition());
    polyline.push_back(rightPlane->getPosition());
}

// Move all planes by the same offset (right plane INCLUDED) - when the slider is dragged
void ViewerFibula::movePlanes(int position){

    //std::cout << "Move planes called" << std::endl;

    int offset = static_cast<int>(static_cast<double>(position)/ static_cast<double>(maxOffset) * static_cast<double>(nbU));

    // Check that it this offset doesn't exceed the size of the fibula
    if(curveIndexL + offset < nbU && curveIndexL + offset > 0 && curveIndexR + offset < nbU && curveIndexR + offset > 0){
        indexOffset = offset;

        // Reset the position and orientations of ALL planes
        leftPlane->setPosition(curve->getPoint(curveIndexL + indexOffset));
        rightPlane->setPosition(curve->getPoint(curveIndexR + indexOffset));
        //leftPlane->setOrientation(getNewOrientation(curveIndexL + indexOffset));      // here
        //rightPlane->setOrientation(getNewOrientation(curveIndexR + indexOffset));

        for(unsigned int i=0; i<ghostPlanes.size(); i++){
            ghostPlanes[i]->setPosition(curve->getPoint(ghostLocation[i] + indexOffset));
            //ghostPlanes[i]->setOrientation(getNewOrientation(ghostLocation[i] + indexOffset));
        }
    }

    setPlaneOrientations(mandiblePolyline, mandibleAxes);
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

    //std::cout << "add ghost planes fibula called " << std::endl;

    for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];
    ghostPlanes.clear();

    for(unsigned int i=0; i<static_cast<unsigned int>(nb); i++){
        ghostPlanes.push_back(new Plane(25.0, Movable::STATIC));
        int index = ghostLocation[i];

        // If we're too far along the fibula, take it all back
        int overload = index + indexOffset - curve->getNbU() + 1;   // The amount by which the actual index passes the end of the curve
        if(overload > 0){
            indexOffset -= overload;
            reinitialisePlanes(i+1);
            // Q_EMIT setPlaneSliderValue(static_cast<int>( (static_cast<double>(indexOffset)/static_cast<double>(*nbU)) * static_cast<double>(maxOffset) ));
        }
        ghostPlanes[i]->setPosition(curve->getCurve()[index + indexOffset]);
        ghostPlanes[i]->setOrientation(getNewOrientation(index + indexOffset));       // here
        //std::cout << "Add ghost planes fibula called " << std::endl;
    }

    // std::cout << "Nb ghost planes fibula : " << ghostPlanes.size() << std::endl;
    //setPlaneOrientations(angleVectors);
    update();
}

// Find the locations of the ghost planes from the distances from the planes in the mandible
void ViewerFibula::findGhostLocations(int nb, double distance[]){
    ghostLocation.clear();
    //std::cout << "finding ghost locations " << std::endl;
    //std::cout << "*******************************";
    int index = curve->indexForLength(curveIndexL, distance[0]);
    ghostLocation.push_back(index);
    int nextIndex = curve->indexForLength(index, 30.0);        // 30 is a tempory security margin
    ghostLocation.push_back(nextIndex);
    for(unsigned int i=1; i<static_cast<unsigned int>(nb); i++){
        index = curve->indexForLength(ghostLocation[2*i-1], distance[i]);
        ghostLocation.push_back(index);
        int nbU = curve->getNbU();
        nextIndex = curve->indexForLength(index, 30.0);
        ghostLocation.push_back(nextIndex);
        if((nextIndex)<nbU) ghostLocation.push_back(nextIndex);
        else ghostLocation.push_back(nbU-1);
    }
    curveIndexR = curve->indexForLength(ghostLocation[2*static_cast<unsigned int>(nb)-1], distance[nb]);
}

// Re-orientate the planes to correspond to the same angles as the jaw
void ViewerFibula::setPlaneOrientations(std::vector<Vec> mandPolyline, std::vector<Vec> axes){
    if(mandPolyline.size()==0) return;

    // Orientate the left plane
    Vec normal = leftPlane->getNormal();
    Quaternion s = Quaternion(-normal, mandPolyline[0]);  // -normal so it doesnt do a 180 flip
    leftPlane->setOrientation(s.normalized());

    // Orientate the ghost planes
    for(unsigned int i=0; i<ghostPlanes.size(); i++){
        normal = ghostPlanes[i]->getNormal();
        s = Quaternion(normal, mandPolyline[i]);
        ghostPlanes[i]->setOrientation(s.normalized());
    }

    // Orientate the right plane
    normal = rightPlane->getNormal();
    s = Quaternion(normal, mandPolyline[mandPolyline.size()-1]);
    rightPlane->setOrientation(s.normalized());


    // Match the polylines

    if(mandPolyline.size()==0) return;

    std::vector<Vec> fibulaPolyline = getPolyline();

    if(ghostPlanes.size()==0){
        // Project the mand and the fib on the left plane
        Vec mandPoint = rightPlane->getLocalProjection(mandPolyline[1]);
        Vec fibPoint = rightPlane->getLocalProjection(fibulaPolyline[1]);
        // normalise them so they have the same length
        mandPoint.normalize();
        fibPoint.normalize();
        // Get the angle between the two
        double alpha = angle(mandPoint, fibPoint);
        //std::cout << "Angle : " << alpha << std::endl;
        // Rotate BOTH planes
        alpha += M_PI;
        //std::cout << " ROTATING" << std::endl;
        Vec axis = Vec(0,0,1);
        //double r = 0.25;
        //alpha = (M_PI*2.0)*r + M_PI;
        leftPlane->rotatePlane(axis, alpha);
        rightPlane->rotatePlane(axis, alpha);


    }

    else{
        Vec axis = Vec(0,0,1);

        // Left
        Vec mandPoint = leftPlane->getLocalProjection(mandPolyline[0]);
        Vec fibPoint = leftPlane->getLocalProjection(fibulaPolyline[0]);
        mandPoint.normalize();
        fibPoint.normalize();
        double alpha = angle(mandPoint, fibPoint) + M_PI;

        leftPlane->rotatePlane(axis, alpha);
        ghostPlanes[0]->rotatePlane(axis, alpha);

        // Ghost
        for(unsigned int i=1; i<ghostPlanes.size()-2; i+=2){
            mandPoint = ghostPlanes[i]->getLocalProjection(mandPolyline[i+1]);
            fibPoint = ghostPlanes[i]->getLocalProjection(fibulaPolyline[i+1]);
            mandPoint.normalize();
            fibPoint.normalize();
            alpha = angle(mandPoint, fibPoint) + M_PI;

            ghostPlanes[i]->rotatePlane(axis, alpha);
            ghostPlanes[i+1]->rotatePlane(axis, alpha);
        }

        // Right
        unsigned int lastIndex = mandPolyline.size()-1;
        mandPoint = rightPlane->getLocalProjection(mandPolyline[lastIndex]);
        fibPoint = rightPlane->getLocalProjection(fibulaPolyline[lastIndex]);
        mandPoint.normalize();
        fibPoint.normalize();
        alpha = angle(mandPoint, fibPoint) + M_PI;

        rightPlane->rotatePlane(axis, alpha);
        ghostPlanes[lastIndex-2]->rotatePlane(axis, alpha); // the last ghost plane
    }
}

void ViewerFibula::recieveTest(std::vector<Vec> axes){
    if(ghostPlanes.size()!=0) return;

    Vec rotationAxis = Vec(0,0,1);

    Vec xPlaneAxis = Vec(1,0,0);

    Vec xPlaneGoal = axes[0];       // our objective
    xPlaneGoal = rightPlane->getMeshVectorFromLocal(xPlaneGoal);        // convert to world coordinates
    xPlaneGoal = leftPlane->getLocalVector(xPlaneGoal);         // convert to left plane coordinates
    xPlaneGoal.normalize();

    double theta = angle(xPlaneAxis, xPlaneGoal) + M_PI;

    leftPlane->rotatePlane(rotationAxis, theta);

    /*Vec zPlaneAxis = Vec(0,0,1);
    zPlaneAxis = leftPlane->getMeshVectorFromLocal(zPlaneAxis);     // convert to world coordinates

    Vec rotationAxis = leftPlane->getPosition() - rightPlane->getPosition();   // get the position in the fibula space (leave it in world coordinates)
    rotationAxis.normalize();

    Vec zPlaneGoal = axes[0];       // our objective
    zPlaneGoal = rightPlane->getMeshVectorFromLocal(zPlaneGoal);        // convert to world coordinates
    zPlaneGoal.normalize();

    // Project onto a theoretical plane defined by the polyline
    Vec projZ = zPlaneAxis - rotationAxis * (zPlaneAxis * rotationAxis);
    Vec projGoal = zPlaneGoal - rotationAxis * (zPlaneGoal * rotationAxis);
    projZ.normalize();
    projGoal.normalize();
    double theta = angle(projZ, projGoal);      // get the angle between the two points projected on the plane

    //rotationAxis = leftPlane->getLocalVector(rotationAxis);         // get it in terms of the left plane (so left is at 0,0,0)
    //rotationAxis.normalize();

    std::cout << "Angle : " << (theta) * 180.0 / M_PI << std::endl;
    std::cout << "Polyline : " << rotationAxis.x << " " << rotationAxis.y  << " " << rotationAxis.z << std::endl;

    theta += M_PI;

    leftPlane->rotatePlane(rotationAxis, theta);*/



    /*Vec axis = Vec(0,0,1);

    if(ghostPlanes.size()!=0) return;
    // Rotate the plane on the polyline axis to match z
    Vec matchAxis = axes[0];
    matchAxis.normalize();
    Vec worldDirectionZ = rightPlane->getMeshVectorFromLocal(matchAxis);    // convert from right plane to the world
    worldDirectionZ.normalize();
    Vec fibDirectionZ = leftPlane->getLocalVector(worldDirectionZ);    // convert from world to left fibPlane
    fibDirectionZ.normalize();
    std::cout << "Where the z axis should be : " << fibDirectionZ.x << " " << fibDirectionZ.y  << " " << fibDirectionZ.z << std::endl;
    // get the angle to rotate the z axis
    Vec rotationAxe = rightPlane->getPosition();
    rotationAxe.normalize();

    // NEED TO CONVERT THEM TO WORLD COORDINATES
    axis = leftPlane->getMeshVectorFromLocal(axis);
    fibDirectionZ = leftPlane->getMeshVectorFromLocal(fibDirectionZ);
    axis.normalize();
    fibDirectionZ.normalize();

    Vec projAxis = axis - rotationAxe * (axis * rotationAxe);       // projection of the z axis of the plane projected onto a plane defined by the normal of the polyline
    Vec projMatch = fibDirectionZ - rotationAxe * (fibDirectionZ * rotationAxe);
    projAxis.normalize();
    projMatch.normalize();
    double theta = angle(projAxis, projMatch);
    std::cout << "Angle : " << (theta) * 180.0 / M_PI << std::endl;
    //std::cout << "Polyline : " << rotationAxe.x << " " << rotationAxe.y  << " " << rotationAxe.z << std::endl;
    // rotate around its OWN polyline
   // theta = M_PI / 2.0;
    theta += M_PI;
    leftPlane->rotatePlane(rotationAxe, theta);*/
}

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
void ViewerFibula::ghostPlanesRecieved(int nb, double distance[], std::vector<Vec> mandPolyline, std::vector<Vec> axes){
    // if no ghost planes were actually recieved
    if(nb==0){
        for(unsigned int i=0; i<ghostPlanes.size(); i++) delete ghostPlanes[i];
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
    setPlaneOrientations(mandPolyline, axes);
    mandiblePolyline.clear();
    mandiblePolyline = mandPolyline;
    mandibleAxes.clear();
    mandibleAxes = axes;

    // If its cut and the number of planes has changed
    if(mesh.getIsCut() && nb!=oldNb){
        mesh.deleteGhostPlanes();
        cutMesh();
    }

    isPlanesRecieved = true;
    handleCut();
}

// When we want to move the right plane (the right plane is moved in the jaw)
void ViewerFibula::movePlaneDistance(double distance, std::vector<Vec> mandPolyline, std::vector<Vec> axes){
    int newIndex;

    //if(isGhostPlanes) return;

    //std::cout << "move plane distance fibula called " << std::endl;
    //std::cout << "*******************************";
    if(ghostPlanes.size()==0) newIndex = curve->indexForLength(curveIndexL, distance);
    else newIndex = curve->indexForLength(ghostLocation[ghostPlanes.size()-1], distance);

    if(newIndex + indexOffset >= nbU) return;      // This should never happen
    else curveIndexR = newIndex;

    rightPlane->setPosition(curve->getCurve()[curveIndexR + indexOffset]);
    //rightPlane->setOrientation(getNewOrientation(curveIndexR + indexOffset));   // initial orientation

    setPlaneOrientations(mandPolyline, axes);
    mandiblePolyline.clear();
    mandiblePolyline = mandPolyline;
    mandibleAxes.clear();
    mandibleAxes = axes;

    mesh.updatePlaneIntersections(rightPlane);
    update();
}

// One of the ghost planes is moved in the jaw
void ViewerFibula::middlePlaneMoved(int nb, double distances[], std::vector<Vec> mandPolyline, std::vector<Vec> axes){
    if(nb==0) return;

    findGhostLocations(nb, distances);

    // update the ghost planes
    for(unsigned int i=0; i<static_cast<unsigned int>(2*nb); i++){
        ghostPlanes[i]->setPosition(curve->getCurve()[ghostLocation[i] + indexOffset]);
        // ghostPlanes[i]->setOrientation(getNewOrientation(ghostLocation[i] + indexOffset));
        // matchPlaneToFrenet(ghostPlanes[i], ghostLocation[i]+indexOffset);
    }

    // If we're too far along the fibula, take it all back
    int overload = curveIndexR + indexOffset - curve->getNbU() + 1;   // The amount by which the actual index passes the end of the curve
    if(overload > 0){
        indexOffset -= overload;
        reinitialisePlanes(ghostPlanes.size()+1); // reinit everything but the right plane
        Q_EMIT setPlaneSliderValue(static_cast<int>( (static_cast<double>(indexOffset)/static_cast<double>(nbU)) * static_cast<double>(maxOffset) ));
    }

    // update the right plane
    rightPlane->setPosition(curve->getCurve()[curveIndexR + indexOffset]);
    //rightPlane->setOrientation(getNewOrientation(curveIndexR + indexOffset));
    // matchPlaneToFrenet(rightPlane, curveIndexR+indexOffset);

    setPlaneOrientations(mandPolyline, axes);
    mandiblePolyline.clear();
    mandiblePolyline = mandPolyline;
    mandibleAxes.clear();
    mandibleAxes = axes;

    mesh.updatePlaneIntersections(rightPlane);

    update();
}

void ViewerFibula::reinitialisePlanes(unsigned int nbToInit){
    if(nbToInit==0) return;
    //std::cout << "reinit planes fibula called " << std::endl;
    //Move the left plane and all the already initialised ghost planes back
    leftPlane->setPosition(curve->getCurve()[curveIndexL + indexOffset]);
    leftPlane->setOrientation(getNewOrientation(curveIndexL + indexOffset));      // here

    if(nbToInit>=ghostPlanes.size()+2){
        rightPlane->setPosition(curve->getCurve()[curveIndexR + indexOffset]);
        rightPlane->setOrientation(getNewOrientation(curveIndexR + indexOffset));
        nbToInit = ghostPlanes.size() + 1;
    }

    for(unsigned int j=0; j<nbToInit-1; j++){
        ghostPlanes[j]->setPosition(curve->getCurve()[ghostLocation[j] + indexOffset]);
        ghostPlanes[j]->setOrientation(getNewOrientation(ghostLocation[j] + indexOffset));
        // matchPlaneToFrenet(ghostPlanes[j], ghostLocation[j]+indexOffset);
    }

    //setPlaneOrientations(angleVectors);
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

    nbU = 800;

    int nbSeg = nbCP-3;
    nbU -= nbU%nbSeg;

    curve->generateCatmull(nbU);
    //curve->generateBSpline(nbU, 3);
    connect(curve, &Curve::curveReinitialised, this, &Viewer::updatePlanes);

    initPlanes(Movable::STATIC);
}

void ViewerFibula::cutMesh(){
    //std::cout << "Cut mesh" << std::endl;
    /*Vec axis = Vec(1,0,0);
    leftPlane->rotatePlane(axis, M_PI);
    rightPlane->rotatePlane(axis, M_PI);*/
    isCutSignal = true;
    handleCut();
    //update();
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
    ghostPlanes.clear();        // NOTE To eventually be changed
    // Reset their orientations
    //leftPlane->setOrientation(getNewOrientation(curveIndexL));
    //rightPlane->setOrientation(getNewOrientation(curveIndexR));
    //matchPlaneToFrenet(leftPlane, curveIndexL);
    //matchPlaneToFrenet(rightPlane, curveIndexR);
    update();
}


void ViewerFibula::handleMovementStart(){

}

void ViewerFibula::handleMovementEnd(){

}

#include "polyline.h"
#include "Tools/point3.h"

Polyline::Polyline()
{
    frame = Frame();
}

// This is only called once. Used to set the reference frame.
void Polyline::init(const Frame *const refFrame, unsigned int nbPoints){
    //frame.setReferenceFrame(refFrame);
    reinit(nbPoints);
}

/* A simple reinitialisation where we set all the normals to the polyline normal
 * and create an evenly spaced polyline with the correct number of points.
 * Note: the points are not correctly positioned here, this is simply a initialisation.
*/
void Polyline::reinit(unsigned int nbPoints){
    points.clear();
    boxes.clear();
    segmentNormals.clear();
    segmentBinormals.clear();
    cuttingLines.clear();
    cuttingBinormals.clear();

    for(unsigned int i=0; i<nbPoints; i++) points.push_back(Vec(i, 0, 0));
    for(unsigned int i=0; i<points.size()-1; i++) segmentNormals.push_back(normal);
    for(unsigned int i=0; i<points.size()-1; i++) segmentBinormals.push_back(binormal);
    for(unsigned int i=1; i<points.size()-1; i++) cuttingLines.push_back(normal);
    for(unsigned int i=1; i<points.size()-1; i++) cuttingBinormals.push_back(binormal);
    for(unsigned int i=0; i<points.size()-1; i++){
        boxes.push_back(Box());
        boxes[i].init(frame.referenceFrame());
    }
    initManipulators();
    initCornerManipulators();
}

void Polyline::draw(){
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_DEPTH);

    glPushMatrix();
    glMultMatrixd(frame.matrix());

    glColor3f(1.,1.,1.);

    // QGLViewer::drawAxis(40.);

    if(isDrawLine){
        // The polyline
        glLineWidth(5.);
        glColor3f(0,0,1);
        glBegin(GL_LINES);
        for(unsigned int i=0; i<points.size()-1; i++){
            glVertex3d(points[i].x, points[i].y, points[i].z);
            glVertex3d(points[i+1].x, points[i+1].y, points[i+1].z);
        }
        glEnd();

        // The points
        glColor3f(0,1,0);
        glPointSize(10.);
        glBegin(GL_POINTS);
        for(unsigned int i=0; i<points.size(); i++) glVertex3d(points[i].x, points[i].y, points[i].z);
        glEnd();
    }

    if(isDrawBoxes){
        glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        for(unsigned int i=1; i<boxes.size()-1; i++){
            glColor4f(0,0,0, boxTransparency);
            boxes[i].draw(0);
        }
        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    }

    glPopMatrix();

    for(unsigned int i=0; i<boxManipulators.size(); i++) boxManipulators[i]->draw();
    for(unsigned int i=0; i<cornerManipulators.size(); i++) cornerManipulators[i]->draw();

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_DEPTH);
}

void Polyline::drawBox(unsigned int index){
    double size = 25.;
    const Vec& p0 = points[index];
    const Vec& p1 = points[index+1];
    Vec p2 = points[index] + segmentBinormals[index]*size;
    Vec p3 = points[index+1] + segmentBinormals[index]*size;
    Vec p4 = points[index] + segmentNormals[index]*size;
    Vec p5 = points[index+1] + segmentNormals[index]*size;
    Vec p6 = points[index] + (segmentNormals[index] +  segmentBinormals[index])*size;
    Vec p7 = points[index+1] + (segmentNormals[index] +  segmentBinormals[index])*size;

    glBegin(GL_LINES);
        glVertex3d(p0.x, p0.y, p0.z);
        glVertex3d(p1.x, p1.y, p1.z);

        glVertex3d(p0.x, p0.y, p0.z);
        glVertex3d(p4.x, p4.y, p4.z);

        glVertex3d(p0.x, p0.y, p0.z);
        glVertex3d(p2.x, p2.y, p2.z);

        glVertex3d(p5.x, p5.y, p5.z);
        glVertex3d(p1.x, p1.y, p1.z);

        glVertex3d(p3.x, p3.y, p3.z);
        glVertex3d(p1.x, p1.y, p1.z);

        glVertex3d(p7.x, p7.y, p7.z);
        glVertex3d(p5.x, p5.y, p5.z);

        glVertex3d(p7.x, p7.y, p7.z);
        glVertex3d(p3.x, p3.y, p3.z);

        glVertex3d(p7.x, p7.y, p7.z);
        glVertex3d(p6.x, p6.y, p6.z);

        glVertex3d(p4.x, p4.y, p4.z);
        glVertex3d(p6.x, p6.y, p6.z);

        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p6.x, p6.y, p6.z);

        glVertex3d(p4.x, p4.y, p4.z);
        glVertex3d(p5.x, p5.y, p5.z);

        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p3.x, p3.y, p3.z);

    glEnd();
}

void Polyline::deleteManipulators(){
    for(unsigned int i=0; i<boxManipulators.size(); i++) delete boxManipulators[i];
    boxManipulators.clear();
}

void Polyline::deleteCornerManipulators(){
    for(unsigned int i=0; i<cornerManipulators.size(); i++) delete cornerManipulators[i];
    cornerManipulators.clear();
}

void Polyline::initManipulators(){
    deleteManipulators();

    for(unsigned int i=0; i<boxes.size(); i++){
        boxManipulators.push_back(new SimpleManipulator);
        boxManipulators[i]->deactivate();
        boxManipulators[i]->setDisplayScale(manipulatorSize);
        boxManipulators[i]->setID(i);
    }

    setManipulatorsToBoxes();
}

void Polyline::initCornerManipulators(){
    deleteCornerManipulators();

    for(unsigned int i=0; i<boxes.size()*2; i++){
        cornerManipulators.push_back(new SimpleManipulator);
        cornerManipulators[i]->deactivate();
        cornerManipulators[i]->setDisplayScale(manipulatorSize);
        cornerManipulators[i]->setID(i);
        cornerManipulators[i]->setRotationActivated(false);
    }

    setCornerManipulatorsToBoxes();
}

void Polyline::setBoxToManipulator(unsigned int id, Vec manipulatorPosition, int s, Vec bX, Vec bY, Vec bZ){
    Vec x,y,z;
    boxManipulators[id]->getOrientation(x,y,z);

    // Set the orientation
    isAngleViolation(id, x, y, z, s, bX, bY, bZ);

   // Set box to the manipulator position - half the length
   manipulatorPosition = isDistanceViolation(id, manipulatorPosition);
   boxManipulators[id]->setOrigin(manipulatorPosition);
   Vec p = getLocalCoordinates(manipulatorPosition - (boxes[id].getLength()/2. * getWorldBoxTransform(id, boxes[id].getTangent()) + boxes[id].getHeight()/2. * getWorldBoxTransform(id, boxes[id].getBinormal()) + boxes[id].getWidth()/2. * getWorldBoxTransform(id, boxes[id].getNormal())));
   boxes[id].setPosition(p);

    setCornerManipulatorsToBoxes();
}

void Polyline::isAngleViolation(unsigned int id, Vec &x, Vec &y, Vec &z, int s, Vec &bX, Vec &bY, Vec &bZ){
     double maxRotation = M_PI / 8.;

    const Vec& n = segmentNormals[id];
    const Vec& b = segmentBinormals[id];
    const Vec t  =-cross(n,b);

    // x -> t, y -> b, z -> n
    bool isSet = false;

    switch (abs(s)) {
    case 6:     // xy
        if(abs(angle(t,x)) > maxRotation) {setToMaxRotation(id, maxRotation, 1, t, bX, x, bX, bY, bZ, s); isSet = true;}
        //if(abs(angle(b,y)) > maxRotation) {setToMaxRotation(id, maxRotation, 0, b, bY, y, bX, bY, bZ, s); isSet = true;}
    break;

    case 5:     // xz
        if(abs(angle(n,z)) > maxRotation) {setToMaxRotation(id, maxRotation, 0, n, bZ, z, bX, bY, bZ, s); isSet = true;}
        //if(abs(angle(t,x)) > maxRotation) {setToMaxRotation(id, maxRotation, 2, t, bX, x, bX, bY, bZ, s); isSet = true;}
        break;

    case 4:    // zy
        if(abs(angle(n,z)) > maxRotation) {setToMaxRotation(id, maxRotation, 1, n, bZ, z, bX, bY, bZ, s); isSet = true;}
        //if(abs(angle(b,y)) > maxRotation) {setToMaxRotation(id, maxRotation, 2, b, bY, y, bX, bY, bZ, s); isSet = true;}
        break;

    default:
        break;
    }

    if(!isSet) boxes[id].setFrameFromBasis(x,y,z);
}

void Polyline::setToMaxRotation(unsigned int id, const double &maxRotation, int comp, Vec a, Vec c, Vec &x, Vec &bX, Vec &bY, Vec &bZ, int s){
    double theta = angle(a, c);

    Quaternion q(a, x);     // The current rotation between the polyline and manipulator axis
    Frame f;
    Quaternion base;
    base.setFromRotatedBasis(bX, bY, bZ);       // Initialise a frame to the previous manipulator axis
    f.setOrientation(base);

    // compare c and x, get the direction
   /* Vec lx = f.localTransformOf(x);
    Vec la = f.localTransformOf(a);
    if(lx[comp] > 0) s *= -1;

    double alpha;
    int sSign = s/abs(s);
    int lSign = 1;
    if(la[comp] < 0) lSign = -1;
    std::cout << sSign << " , " << lSign << std::endl;
    if(sSign == lSign) alpha = maxRotation - theta;
    else alpha = maxRotation + theta;

    if(s < 0) alpha *= -1;*/

    double alpha = maxRotation - theta;

    Quaternion maxQ(q.axis(), alpha);     // The rotation between the polyline and the previous manipulator axis

    f.rotate(maxQ);         // Rotate the previous manipulator to its max
    Vec x0 = f.localInverseTransformOf(Vec(1,0,0));     // Get the x,y,z vectors of the max frame
    Vec y0 = f.localInverseTransformOf(Vec(0,1,0));
    Vec z0 = f.localInverseTransformOf(Vec(0,0,1));
    boxes[id].setFrameFromBasis(x0,y0,z0);
    boxManipulators[id]->setRepX(x0);
    boxManipulators[id]->setRepY(y0);
    boxManipulators[id]->setRepZ(z0);
}


Vec Polyline::isDistanceViolation(unsigned int id, const Vec &manipulatorPosition){
    double maxDistance = 25.;

    Vec midPoint = (points[id] + points[id+1]) / 2.0;
    double d = euclideanDistance(midPoint, manipulatorPosition);

    if(d <= maxDistance) return manipulatorPosition;

    Vec furtherestPosition = (manipulatorPosition - midPoint);
    furtherestPosition.normalize();

    return midPoint + furtherestPosition*maxDistance;
}

void Polyline::setBoxToManipulatorOrientation(unsigned int id){
    Vec x,y,z;
    boxManipulators[id]->getOrientation(x,y,z);
    boxes[id].setFrameFromBasis(x,y,z);
    setCornerManipulatorsToBoxes();
}

void Polyline::setBoxToCornerManipulator(unsigned int id, Vec manipulatorPosition){
    //if(id%2==0) boxes[id/2].setPosition(manipulatorPosition);       // TODO block the end and recalculate the box
    //else boxes[id/2].setPosition(p);
    unsigned int boxID = id/2;
    if(id%2==0){
        reorientateBox(boxID, manipulatorPosition, getMeshBoxEnd(boxID));
        boxes[boxID].setPosition(manipulatorPosition);
    }
    else{
        reorientateBox(boxID, getMeshBoxPoint(boxID), manipulatorPosition);
    }

    setManipulatorsToBoxes();
}

void Polyline::reorientateBox(unsigned int index, const Vec &start, const Vec &end){
    recalculateBinormal(index, start, end);

    const Vec& n = segmentNormals[index];
    const Vec& b = segmentBinormals[index];

    boxes[index].setFrameFromBasis(-cross(n,b),b,n);
}

void Polyline::setBoxToProjectionPoint(unsigned int id, Vec projPoint){
    Vec p = getLocalCoordinates(projPoint);
    boxes[id].setPosition(p);
}

void Polyline::setManipulatorsToBoxes(){
    Vec x,y,z;
    for(unsigned int i=0; i<boxes.size(); i++){
       boxes[i].getOrientation(x,y,z);
       boxManipulators[i]->setRepX(x);
       boxManipulators[i]->setRepY(y);
       boxManipulators[i]->setRepZ(z);

       Vec p = getMeshBoxPoint(i) + boxes[i].getLength()/2. * getWorldBoxTransform(i, boxes[i].getTangent()) + boxes[i].getHeight()/2. * getWorldBoxTransform(i, boxes[i].getBinormal()) + boxes[i].getWidth()/2. * getWorldBoxTransform(i, boxes[i].getNormal());
       boxManipulators[i]->setOrigin(p);
    }
}

void Polyline::setManipulatorOrientations(std::vector<Vec> &x, std::vector<Vec> &y, std::vector<Vec> &z){
    for(unsigned int i=0; i< x.size(); i++){
        boxManipulators[i]->setRepX(x[i]);
        boxManipulators[i]->setRepY(y[i]);
        boxManipulators[i]->setRepZ(z[i]);
    }
}

void Polyline::setCornerManipulatorsToBoxes(){
    Vec x,y,z;
    for(unsigned int i=0; i<boxes.size(); i++){
       boxes[i].getOrientation(x,y,z);
       cornerManipulators[i*2]->setRepX(x);
       cornerManipulators[i*2]->setRepY(y);
       cornerManipulators[i*2]->setRepZ(z);

       cornerManipulators[i*2+1]->setRepX(x);
       cornerManipulators[i*2+1]->setRepY(y);
       cornerManipulators[i*2+1]->setRepZ(z);

       cornerManipulators[i*2]->setOrigin(getMeshBoxPoint(i));

       //Vec p = getMeshBoxPoint(i) + boxes[i].getLength()/2. * getWorldBoxTransform(i, boxes[i].getTangent()) + boxes[i].getHeight()/2. * getWorldBoxTransform(i, boxes[i].getBinormal()) + boxes[i].getWidth()/2. * getWorldBoxTransform(i, boxes[i].getNormal());
       cornerManipulators[i*2+1]->setOrigin(getMeshBoxEnd(i));
    }

}

void Polyline::activateBoxManipulators(const bool &b){
    // Don't activate the first and the last - these boxes don't count
    for(unsigned int i=1; i<boxManipulators.size()-1; i++) boxManipulators[i]->setState(b);

}

void Polyline::activateFirstCornerManipulators(const bool &b){
    for(unsigned int i=2; i<cornerManipulators.size()-2; i+=2) cornerManipulators[i]->setState(b);
}

void Polyline::activateEndCornerManipulators(const bool &b){
    for(unsigned int i=3; i<cornerManipulators.size()-2; i+=2) cornerManipulators[i]->setState(b);
}

void Polyline::toggleBoxManipulators(unsigned int i, const bool &b){
    boxManipulators[i]->setState(b);
}

void Polyline::toggleFirstCornerManipulators(unsigned int i, const bool &b){
    cornerManipulators[i*2]->setState(b);
}

void Polyline::toggleEndCornerManipulators(unsigned int i, const bool &b){
    cornerManipulators[i*2+1]->setState(b);
}

// Update the points locations without updating their orientations
void Polyline::updatePoints(const std::vector<Vec> &newPoints){
    points.clear();
    for(unsigned int i=0; i<newPoints.size(); i++) points.push_back(newPoints[i]);
    for(unsigned int i=0; i<boxes.size(); i++){
        resetBox(i);
        boxes[i].restoreRotation();
        // setManipulatorsToBoxes();
    }
}

double Polyline::angle(const Vec &a, const Vec &b){
    double na = a.norm();
    double nb = b.norm();
    double ab = a*b;

    double val = ab / (na*nb);
    if(val >= static_cast<double>(1)) val = 1;          // protection from floating point errors (comparing it to an epsilon didn't work)
    else if(val < static_cast<double>(-1)) val = -1;
    return acos(val);
}

Vec Polyline::projection(Vec &a, const Vec &planeNormal){
    double alpha = (a * planeNormal);
    return a - planeNormal * alpha;
}

// Match the box to the current position and orientation of a point index
void Polyline::resetBox(unsigned int index){
    boxes[index].setPosition(points[index]);
    const Vec& n = segmentNormals[index];
    const Vec& b = segmentBinormals[index];

    boxes[index].setFrameFromBasis(-cross(n,b),b,n);
    double length = euclideanDistance(points[index], points[index+1]);
    boxes[index].setLength(length);
    setManipulatorsToBoxes();
    setCornerManipulatorsToBoxes();
}

void Polyline::bend(unsigned int index, Vec &newPosition, std::vector<Vec>& planeNormals, std::vector<Vec>& planeBinormals){
    if(index >= points.size()) return;

    bendFibula(index, newPosition);     // the fibula is bent in the same way but does not need to send info

    getCuttingAngles(planeNormals, planeBinormals);     // Get the normals and binormals of the planes which interset two boxes at a polyline point
}

void Polyline::bendFibula(unsigned int index, Vec &newPosition){
    if(index >= points.size()) return;

    points[index] = getLocalCoordinates(newPosition);

    // Recalculate the orientations of the two boxes attached to the point index (the box behind and the box infront)
    if(index!=0){
        recalculateBinormal(index-1, points[index-1], points[index]);
        resetBox(index-1);
        boxes[index-1].restoreRotation();
    }
    if(index!=points.size()-1){
        recalculateBinormal(index, points[index], points[index+1]);
        resetBox(index);
        boxes[index].restoreRotation();
    }
}

void Polyline::recalculateOrientations(){
    for(unsigned int i=0; i<points.size()-1; i++) recalculateBinormal(i, points[i], points[i+1]);
}

void Polyline::recalculateBinormal(unsigned int index, const Vec &origin, const Vec &newPosition){
    point3d  Nprev = normal;
    point3d  Tprev = tangent;

    point3d  T0 = newPosition - origin;        // the polyline tangent
    point3d  N0 = point3d::rotateVectorSimilarly( Nprev , Tprev , T0  );

    point3d  B0 = point3d::cross( T0 , N0 );

    segmentNormals[index] = Vec(N0);
    segmentBinormals[index] = -Vec(B0);
    segmentNormals[index].normalize();
    segmentBinormals[index].normalize();
}

Vec Polyline::vectorQuaternionRotation(double theta, const Vec &axis, const Vec &vectorToRotate){
    Quaternion r(cos(theta/2.0)*axis.x, cos(theta/2.0)*axis.y, cos(theta/2.0)*axis.z, sin(theta/2.0));      // rotation

    Frame f = Frame();
    initialiseFrame(f);
    f.rotate(r);

    return f.localInverseTransformOf(vectorToRotate);
}

void Polyline::initialiseFrame(Frame &f){
    Quaternion q = Quaternion();            // the base
    q.setFromRotatedBasis(getWorldTransform(Vec(1,0,0)), getWorldTransform(Vec(0,1,0)),getWorldTransform(Vec(0,0,1)));
    f.setOrientation(q);
}

// Get the normals and binormals of the planes which interset two boxes at a polyline point in order to send it to the fibula
void Polyline::getCuttingAngles(std::vector<Vec>& planeNormals, std::vector<Vec>& planeBinormals){
    cuttingLines.clear();
    cuttingBinormals.clear();
    planeNormals.clear();
    planeBinormals.clear();

    planeNormals.push_back(segmentNormals[0]);      // for the first plane (there's only one box attached to the first and last point)
    planeBinormals.push_back(segmentBinormals[0]);

    // Get the angle between the two boxes
    for(unsigned int i=0; i<segmentNormals.size()-1; i++){
        Vec v = (segmentNormals[i] + segmentNormals[i+1]);
        v.normalize();
        cuttingLines.push_back(v);

        Vec b = (segmentBinormals[i] + segmentBinormals[i+1]);
        b.z /= 2.;
        b.normalize();
        cuttingBinormals.push_back(b);
    }

    // Make sure the normal and binormal have a right angle between them
    for(unsigned int i=0; i<cuttingLines.size(); i++){
        double theta = angle(cuttingLines[i], cuttingBinormals[i]);     // Get the current angle between them
        double alpha = M_PI / 2.0 - theta + M_PI;           // Get the angle it needs to rotate in order for it to be a right angle
        Vec axis = cross(cuttingLines[i], cuttingBinormals[i]);     // Get the axis on which it has to rotate (which is orthogonal to the current normal and binormal)
        cuttingLines[i] = vectorQuaternionRotation(alpha, axis, cuttingLines[i]);       // Rotate the normal so its now at a right angle to the binormal

        planeNormals.push_back(cuttingLines[i]);            // save for the mandible
        planeBinormals.push_back(cuttingBinormals[i]);
    }

    planeNormals.push_back(segmentNormals.back());      // for the last plane
    planeBinormals.push_back(segmentBinormals.back());
}

// Get the distances between the polyline points
void Polyline::getDistances(std::vector<double> &distances){
    distances.clear();

    /*for(unsigned int i=0; i<points.size()-1; i++){
        distances.push_back(euclideanDistance(points[i], points[i+1]));
    }*/

    for(unsigned int i=0; i<boxes.size(); i++){
        distances.push_back(boxes[i].getLength());
    }

    distances[0] = 0.0001;       // We don't want an offet in the fibula, so set the first box to nearly zero (null vector if zero)
}

double Polyline::euclideanDistance(const Vec &a, const Vec &b){
    return sqrt(pow(a.x-b.x, 2.) + pow(a.y-b.y, 2.) + pow(a.z-b.z, 2.));
}

// Move the polyline point by the vector toLower
void Polyline::lowerPoint(unsigned int index, const Vec &toLower){
    points[index] += toLower;
}

// Reset the boxes to the current orientation, position and length of the polyline
void Polyline::resetBoxes(){
    recalculateOrientations();
    for(unsigned int i=0; i<boxes.size(); i++) resetBox(i);
}

// Rotate the box on its own axis
void Polyline::rotateBox(unsigned int i, double angle){
    boxes[i].rotateOnAxis(angle);
}

void Polyline::restoreBoxRotations(){
    for(unsigned int i=0; i<boxes.size(); i++) boxes[i].restoreRotation();
}

// Convert the vectors from the box frame to the world frame
void Polyline::getRelatvieNormals(std::vector<Vec> &relativeNorms){
    for(unsigned int i=0; i<relativeNorms.size(); i++){
        unsigned int boxID = i/4+1;
        relativeNorms[i] = getWorldTransform(boxes[boxID].worldTransform(relativeNorms[i]));
    }
}

// Get each plane in relation to each of its corresponding boxes
void Polyline::getRelativePlane(Plane &p, std::vector<Vec> &norms){
    // Get the plane norms in terms of the mesh
    Vec n(1,0,0);
    Vec b(0,1,0);
    n = p.getMeshVectorFromLocal(n);
    b = p.getMeshVectorFromLocal(b);

    // Now get it in terms of its boxes
    norms.clear();
    const unsigned int& id = p.getID();
    if(id!=0){
        norms.push_back(boxes[id-1].localTransform(n));
        norms.push_back(boxes[id-1].localTransform(b));
    }
    if(id<boxes.size()){
        norms.push_back(boxes[id].localTransform(n));
        norms.push_back(boxes[id].localTransform(b));
    }
}

// Get the vector direction of each segment of the polyline (in terms of the world)
void Polyline::getDirections(std::vector<Vec> &directions){
    directions.clear();
    for(unsigned int i=0; i<boxes.size(); i++) directions.push_back(getWorldTransform(boxes[i].worldTangent()));
}

// Move the entire polyline
void Polyline::lowerPolyline(Vec localDirection, double distance){
    Vec p = frame.position();
    Vec worldDirection = getWorldTransform(localDirection);     // convert to world coordinates
    frame.setPosition(p+worldDirection*distance);
}

void Polyline::adjustBoxLength(unsigned int i, double &distShift){
    boxes[i].setLength(distShift);
}

// Set the polyline points to correspond to the ghost planes. This doesn't update anything. Used exclusively for cntrl-z
void Polyline::setPointsToGhostPlanes(std::vector<Vec> &planePositions){
    for(unsigned int i=0; i<planePositions.size(); i++){
        points[i+2] = planePositions[i];
    }
}

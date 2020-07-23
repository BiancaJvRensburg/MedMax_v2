#include "plane.h"

Plane::Plane(double s, Movable status, float alpha, unsigned int id) : cp(id)
{
    size = s;
    rotationPercentage = 0;
    normal = Vec(0, 0, 1);
    constraint.setRotationConstraint(AxisPlaneConstraint::AXIS, Vec(0,0,1));
    constraintFree.setRotationConstraint(AxisPlaneConstraint::FREE, Vec(0,0,1));

    this->status = status;
    this->isVisible = true;
    this->alpha = alpha;
    this->isPoly = false;
    this->id = id;
    this->displayDimensions = Vec(s,s,0);
    manipulator.setDisplayScale(size/3.);
    manipulator.deactivate();
    manipulator.setID(id);
    manipulator.setRotationActivated(false);

    initBasePlane();
}

// Toggle how the plane is placed in terms of the polyline
void Plane::toggleIsPoly(){
    isPoly = !isPoly;
    initBasePlane();
}

// Initialise the plane points
void Plane::initBasePlane(){
    if(isPoly){     // If the polyline exits, put the curve point on the edge of the plane
        points[0] = Vec(cp.getPoint().x, cp.getPoint().y, cp.getPoint().z);
        points[1] = Vec(cp.getPoint().x, cp.getPoint().y + 2.*size, cp.getPoint().z);
        points[2] = Vec(cp.getPoint().x + 2.*size, cp.getPoint().y + 2.*size, cp.getPoint().z);
        points[3] = Vec(cp.getPoint().x + 2.*size, cp.getPoint().y, cp.getPoint().z);
    }
    else{       // Put it in the middle of the plane
        points[0] = Vec(cp.getPoint().x - size, cp.getPoint().y - size, cp.getPoint().z);
        points[1] = Vec(cp.getPoint().x - size, cp.getPoint().y + size, cp.getPoint().z);
        points[2] = Vec(cp.getPoint().x + size, cp.getPoint().y + size, cp.getPoint().z);
        points[3] = Vec(cp.getPoint().x + size, cp.getPoint().y - size, cp.getPoint().z);
    }

    reinitDisplayPoints();
}

void Plane::setDisplayDimensions(double height, double width){
    displayDimensions = Vec(width, height, 0);
    reinitDisplayPoints();
}

void Plane::reinitDisplayPoints(){
    const double& width = displayDimensions.x;
    const double& height = displayDimensions.y;

    if(isPoly){     // If the polyline exits, put the curve point on the edge of the plane
        displayPoints[0] = Vec(cp.getPoint().x, cp.getPoint().y, cp.getPoint().z);
        displayPoints[1] = Vec(cp.getPoint().x, cp.getPoint().y + 2.*height, cp.getPoint().z);
        displayPoints[2] = Vec(cp.getPoint().x + 2.*width, cp.getPoint().y + 2.*height, cp.getPoint().z);
        displayPoints[3] = Vec(cp.getPoint().x + 2.*width, cp.getPoint().y, cp.getPoint().z);
    }
    else{       // Put it in the middle of the plane
        displayPoints[0] = Vec(cp.getPoint().x - width, cp.getPoint().y - height, cp.getPoint().z);
        displayPoints[1] = Vec(cp.getPoint().x - width, cp.getPoint().y + height, cp.getPoint().z);
        displayPoints[2] = Vec(cp.getPoint().x + width, cp.getPoint().y + height, cp.getPoint().z);
        displayPoints[3] = Vec(cp.getPoint().x + width, cp.getPoint().y - height, cp.getPoint().z);
    }
}

// Gets the four corner points of the plane in the world coordinates
void Plane::getCorners(Vec &v0, Vec &v1, Vec &v2, Vec &v3){
    v0 = getMeshCoordinatesFromLocal(points[0]);
    v1 = getMeshCoordinatesFromLocal(points[1]);
    v2 = getMeshCoordinatesFromLocal(points[2]);
    v3 = getMeshCoordinatesFromLocal(points[3]);
}

void Plane::draw(){
    glPushMatrix();
    glMultMatrixd(cp.getFrame().matrix());

    if(isVisible){
        glEnable(GL_DEPTH);
        glEnable(GL_DEPTH_TEST);

        glBegin(GL_QUADS);
            glVertex3d(displayPoints[0].x, displayPoints[0].y, displayPoints[0].z);
            glVertex3d(displayPoints[1].x, displayPoints[1].y, displayPoints[1].z);
            glVertex3d(displayPoints[2].x, displayPoints[2].y, displayPoints[2].z);
            glVertex3d(displayPoints[3].x, displayPoints[3].y, displayPoints[3].z);
        glEnd();

        glDisable(GL_DEPTH);
        glDisable(GL_DEPTH_TEST);
    }

    glColor3f(1,1,1);
    //QGLViewer::drawAxis(size/2.);

    /*if(status==Movable::DYNAMIC){
        cp.toggleSwitchFrames();
        cp.draw();
        cp.toggleSwitchFrames();
    }*/

    glPopMatrix();

    manipulator.draw();
}

void Plane::rotatePlane(Vec axis, double theta){
    rotate(Quaternion(cos(theta/2.0)*axis.x, cos(theta/2.0)*axis.y, cos(theta/2.0)*axis.z, sin(theta/2.0)));
}

void Plane::setPlaneRotation(Vec axis, double theta){
    setRotation(Quaternion(cos(theta/2.0)*axis.x, cos(theta/2.0)*axis.y, cos(theta/2.0)*axis.z, sin(theta/2.0)));
}

void Plane::setPosition(Vec pos){
    cp.setPosition(pos);
    manipulator.setOrigin(pos);
}

void Plane::setManipulatorToOrientation(){
    Vec x,y,z;
    x = cp.getFrame().localInverseTransformOf(Vec(1,0,0));
    y = cp.getFrame().localInverseTransformOf(Vec(0,1,0));
    z = cp.getFrame().localInverseTransformOf(Vec(0,0,1));
    manipulator.setRepX(x);
    manipulator.setRepY(y);
    manipulator.setRepZ(z);
}

// Set the base from a basis x,y,z
Quaternion Plane::fromRotatedBasis(Vec x, Vec y, Vec z){
    Quaternion q = Quaternion();
    q.setFromRotatedBasis(x,y,z);
   return q;
}

// Actually set the orientation from the basis x,y,z
void Plane::setFrameFromBasis(Vec x, Vec y, Vec z){
    cp.getFrame().setOrientation(fromRotatedBasis(x,y,z));
    setManipulatorToOrientation();
}

void Plane::rotatePlaneXY(double percentage){
    double r = (percentage - rotationPercentage);       // Get the percentage to rotate it by
    rotationPercentage = percentage;

    double theta = (M_PI*2.0)*r + M_PI;     // Get the theta from the percentage
    Vec axis = Vec(0,0,1);

    rotatePlane(axis, theta);
}

// MESH INTERSECTION

/*
* Pass the vertices to check in the form of Vec
* The mesh will then use this information
* This is the only thing that the plane deals with
*/
bool Plane::isIntersection(Vec v0, Vec v1, Vec v2){

    // Put it all into local coordinates
    const Vec &tr0 = cp.getFrame().localCoordinatesOf(v0);
    const Vec &tr1 = cp.getFrame().localCoordinatesOf(v1);
    const Vec &tr2 = cp.getFrame().localCoordinatesOf(v2);

    Vec tr[3] = {tr0, tr1, tr2};

    if( (tr0.z < 0 && tr1.z < 0 && tr2.z < 0) || (tr0.z > 0 && tr1.z > 0 && tr2.z > 0) ) return false;  // if they all have the same sign
    else{
        for(int i=0; i<3; i++){
            Vec l = tr[(i+1)%3] - tr[i];

            if(l*normal == 0.0){
                if( tr[i]*normal == 0.0 ){
                    if( (abs(tr[i].x) < size && abs(tr[i].y) < size) || ( abs(tr[(i+1)%3].x) < size && abs(tr[(i+1)%3].y) < size)) return true;  // the plan contains the line
                }
                else continue;  // the line is parallel
            }

            double d = normal*(-tr[i]) / (normal*l);

            if(abs(d) > 1.0) continue;

            Vec intersection = d*l + tr[i];
            if(abs(intersection.x) < size*2. && abs(intersection.y) < size*2.) return true;
        }
    }

    return false;   // if we haven't found a line that meets the criteria
}

// On which sign of the plane is the point located
double Plane::getSign(Vec v){
    const Vec &tr0 = cp.getFrame().localCoordinatesOf(v);
    return tr0.z/(abs(tr0.z));
}

// Project coordinates onto the plane (coordinates in terms of the world)
Vec Plane::getProjection(Vec p){
    const Vec &localP = cp.getFrame().localCoordinatesOf(p);       // convert into local coordinates
    double alpha = (localP * normal);
    Vec newP = localP - normal *alpha;
    return cp.getFrame().localInverseCoordinatesOf(newP);   // convert back into original coordinate system
}

Vec Plane::getAxisProjection(Vec p, Vec axis){
    const Vec &localP = cp.getFrame().localCoordinatesOf(p);       // convert into local coordinates
    Vec localAxis = cp.getFrame().localTransformOf(axis);
    double alpha = localP.z / localAxis.z;
    Vec newP = localP - alpha*localAxis;
    return cp.getFrame().localInverseCoordinatesOf(newP);
}

// Project coordinates onto the plane (coordinates in terms of the plane)
Vec Plane::getLocalProjection(Vec localP){
    return localP - normal * (localP * normal);             // don't convert between coordinate systems
}

// Do two planes intersect? (may obsolete now)
bool Plane::isIntersectionPlane(Vec &v0, Vec &v1, Vec &v2, Vec &v3){

    // Put it all into local coordinates
    const Vec &tr0 = cp.getFrame().localCoordinatesOf(v0);
    const Vec &tr1 = cp.getFrame().localCoordinatesOf(v1);
    const Vec &tr2 = cp.getFrame().localCoordinatesOf(v2);
    const Vec &tr3 = cp.getFrame().localCoordinatesOf(v3);

    Vec tr[4] = {tr0, tr1, tr2, tr3};

    if( (tr0.z < 0 && tr1.z < 0 && tr2.z < 0 && tr3.z > 0) || (tr0.z > 0 && tr1.z > 0 && tr2.z > 0 && tr3.z > 0) ) return false;  // if they all have the same sign
    else{
        for(int i=0; i<4; i++){
            Vec l = tr[(i+1)%4] - tr[i];

            if(l*normal == 0.0){
                if( tr[i]*normal == 0.0 ){
                    if( (abs(tr[i].x) < size && abs(tr[i].y) < size) || ( abs(tr[(i+1)%3].x) < size && abs(tr[(i+1)%3].y) < size)) return true;  // the plan contains the line
                }
                else continue;  // the line is parallel
            }

            double d = normal*(-tr[i]) / (normal*l);

            if(abs(d) > 1.0) continue;

            Vec intersection = d*l + tr[i];
            //std::cout << "Intersection : " << intersection.x << " , " << intersection.y << " , " << intersection.z << std::endl;
            if(abs(intersection.x) < size && abs(intersection.y) < size) return true;
        }
    }

    return false;   // if we haven't found a line that meets the criteria
}

Quaternion Plane::getOrientation(){
    double a, x, y, z;
    cp.getFrame().getOrientation(a,x,y,z);
     return Quaternion(a,x,y,z);
}

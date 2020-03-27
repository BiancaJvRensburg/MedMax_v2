#include "plane.h"

Plane::Plane(double s, Movable status, Vec& pos, float alpha) : cp(pos)
{
    size = s;
    rotationPercentage = 0;
    normal = Vec(0, 0, 1);
    constraint.setRotationConstraint(AxisPlaneConstraint::AXIS, Vec(0,0,1));
    constraintFree.setRotationConstraint(AxisPlaneConstraint::FREE, Vec(0,0,1));

    this->status = status;
    this->isVisible = true;
    this->alpha = alpha;

    initBasePlane();
}

void Plane::initBasePlane(){
        points[0] = Vec(cp.getPoint().x - size, cp.getPoint().y - size, cp.getPoint().z);
        points[1] = Vec(cp.getPoint().x - size, cp.getPoint().y + size, cp.getPoint().z);
        points[2] = Vec(cp.getPoint().x + size, cp.getPoint().y + size, cp.getPoint().z);
        points[3] = Vec(cp.getPoint().x + size, cp.getPoint().y - size, cp.getPoint().z);
}

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
            glVertex3f(static_cast<float>(points[0].x), static_cast<float>(points[0].y), static_cast<float>(points[0].z));
            glVertex3f(static_cast<float>(points[1].x), static_cast<float>(points[1].y), static_cast<float>(points[1].z));
            glVertex3f(static_cast<float>(points[2].x), static_cast<float>(points[2].y), static_cast<float>(points[2].z));
            glVertex3f(static_cast<float>(points[3].x), static_cast<float>(points[3].y), static_cast<float>(points[3].z));
        glEnd();

        glDisable(GL_DEPTH);
        glDisable(GL_DEPTH_TEST);
    }

    glColor3f(1,1,1);
    //QGLViewer::drawAxis(size/2.0);

    if(status==Movable::DYNAMIC){
        cp.toggleSwitchFrames();
        cp.draw();
        cp.toggleSwitchFrames();
    }

    glPopMatrix();
}

void Plane::rotatePlane(Vec axis, double theta){
    rotate(Quaternion(cos(theta/2.0)*axis.x, cos(theta/2.0)*axis.y, cos(theta/2.0)*axis.z, sin(theta/2.0)));
}

void Plane::setPlaneRotation(Vec axis, double theta){
    setRotation(Quaternion(cos(theta/2.0)*axis.x, cos(theta/2.0)*axis.y, cos(theta/2.0)*axis.z, sin(theta/2.0)));
}

void Plane::rotatePlaneXY(double percentage){
    double r = (percentage - rotationPercentage);       // Get the percentage to rotate it by
    rotationPercentage = percentage;

    double theta = (M_PI*2.0)*r + M_PI;     // Get the theta from the percentage
    Vec axis = Vec(0,0,1);

    rotatePlane(axis, theta);
}

void Plane::setPosition(Vec pos){
    cp.setPosition(pos);
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
            if(abs(intersection.x) < size && abs(intersection.y) < size) return true;
        }
    }

    return false;   // if we haven't found a line that meets the criteria
}

double Plane::getSign(Vec v){
    const Vec &tr0 = cp.getFrame().localCoordinatesOf(v);
    return tr0.z/(abs(tr0.z));
}

Vec Plane::getProjection(Vec p){
    const Vec &localP = cp.getFrame().localCoordinatesOf(p);       // convert into local coordinates
    double alpha = (localP * normal);
    Vec newP = localP - normal *alpha;
    return cp.getFrame().localInverseCoordinatesOf(newP);   // convert back into original coordinate system
}

Vec Plane::getLocalProjection(Vec localP){
    return localP - normal * (localP * normal);             // don't convert between coordinate systems
}

Frame Plane::getFrameCopy(){
    return Frame(cp.getFrame());
}

void Plane::setOrientationFromOtherReference(std::vector<Vec> &frame, unsigned int startIndex, Plane *reference){
    Vec x = reference->getMeshVectorFromLocal(frame[startIndex]);
    x.normalize();
    Vec y = reference->getMeshVectorFromLocal(frame[startIndex+1]);
    y.normalize();
    Vec z = reference->getMeshVectorFromLocal(frame[startIndex+2]);
    z.normalize();

    setFrameFromBasis(x,y,z);
}

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
            std::cout << "Intersection : " << intersection.x << " , " << intersection.y << " , " << intersection.z << std::endl;
            if(abs(intersection.x) < size && abs(intersection.y) < size) return true;
        }
    }

    return false;   // if we haven't found a line that meets the criteria
}

void Plane::matchPlane(Plane *p){
    cp.matchCurvepoint(p->cp);
}

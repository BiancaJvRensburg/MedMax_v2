#include "box.h"

Box::Box()
{
    f = Frame();
    dimensions = Vec(100., 15., 15.);
    tangent = Vec(1.,0.,0.);
    binormal = Vec(0.,1.,0.);
    normal = Vec(0.,0.,1.);
    prevRotation = 0.;
    //manipulator = SimpleManipulator();
    //manipulator.setEtat(1);
}

// Set the reference frame
void Box::init(const Frame *ref){
    //f.setReferenceFrame(ref);
}

// Set the reference frame to x,y,z
void Box::setFrameFromBasis(Vec x, Vec y, Vec z){
    x.normalize();
    y.normalize();
    z.normalize();

    Quaternion q;
    q.setFromRotatedBasis(x,y,z);
    f.setOrientation(q);
}

void Box::draw(double offset){
    glPushMatrix();
    glMultMatrixd(f.matrix());

     QGLViewer::drawAxis(10.);

    //manipulator.draw();

    const double& length = getLength();
    const double& width = getWidth();
    const double& height = getHeight();

    const Vec& p0base = Vec(0,0,0);
    const Vec& p1base = p0base + length*tangent;
    Vec p0 = p0base - (normal+binormal)*offset;
    Vec p1 = p1base - (normal+binormal)*offset;
    Vec p2 = p0base + binormal*(width+offset) - normal*offset;
    Vec p3 = p1base + binormal*(width+offset) - normal*offset;
    Vec p4 = p0base + normal*(height+offset)  - binormal*offset;
    Vec p5 = p1base + normal*(height+offset)  - binormal*offset;
    Vec p6 = p0base + normal*height +  binormal*width + (normal+binormal)*offset;
    Vec p7 = p1base + normal*height +  binormal*width + (normal+binormal)*offset;

    glBegin(GL_QUADS);
        glVertex3d(p0.x, p0.y, p0.z);
        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p5.x, p5.y, p5.z);
        glVertex3d(p4.x, p4.y, p4.z);
    glEnd();

    glBegin(GL_QUADS);
        glVertex3d(p0.x, p0.y, p0.z);
        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p3.x, p3.y, p3.z);
        glVertex3d(p2.x, p2.y, p2.z);
    glEnd();

    glBegin(GL_QUADS);
        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p3.x, p3.y, p3.z);
        glVertex3d(p7.x, p7.y, p7.z);
        glVertex3d(p6.x, p6.y, p6.z);
    glEnd();

    glBegin(GL_QUADS);
        glVertex3d(p6.x, p6.y, p6.z);
        glVertex3d(p7.x, p7.y, p7.z);
        glVertex3d(p5.x, p5.y, p5.z);
        glVertex3d(p4.x, p4.y, p4.z);
    glEnd();

    glBegin(GL_QUADS);
        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p3.x, p3.y, p3.z);
        glVertex3d(p7.x, p7.y, p7.z);
        glVertex3d(p5.x, p5.y, p5.z);
    glEnd();

    glBegin(GL_QUADS);
        glVertex3d(p0.x, p0.y, p0.z);
        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p6.x, p6.y, p6.z);
        glVertex3d(p4.x, p4.y, p4.z);
    glEnd();

    glPopMatrix();
}

// Rotate the box on its own axis (so around its centre + tangent)
void Box::rotateOnAxis(double angle){
    double alpha = angle - prevRotation;
    prevRotation = angle;
    Vec centre = f.localInverseCoordinatesOf(dimensions/2.);

    Vec axis(1,0,0);
    Quaternion r(axis, alpha);

    f.rotateAroundPoint(r, centre);
}

// Restore the box's last rotation, in case the box was reset
void Box::restoreRotation(){
    Vec centre = f.localInverseCoordinatesOf(dimensions/2.);

    Vec axis(1,0,0);
    Quaternion r(axis, prevRotation);

    f.rotateAroundPoint(r, centre);
}

Vec Box::getLocation(){
    return f.position();
}

Vec Box::getEnd(){
    Vec t = tangent*getLength();

    return worldCoordinates(t);
}

Vec Box::getMidPoint(){
    Vec m = binormal*getWidth();

    return worldCoordinates(m);
}

Vec Box::getHighPoint(){
    Vec h = normal*getHeight();

    return worldCoordinates(h);
}

Vec Box::getHighEnd(){
    return getHighPoint() + getEnd();
}

void Box::getOrientation(Vec &x, Vec &y, Vec &z){
    Vec a(1,0,0);
    Vec b(0,1,0);
    Vec c(0,0,1);

    x = worldTransform(a);
    y = worldTransform(b);
    z = worldTransform(c);
}

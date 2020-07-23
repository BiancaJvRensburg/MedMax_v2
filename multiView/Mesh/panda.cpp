#include "panda.h"
#include "meshreader.h"

Panda::Panda()
{
    openOFF("effector.off", effector);
    openOFF("marker.off", marker);
    openOFF("Navex.off", navex);
    openOFF("WRJ.off", wrj);

    setFrames();
}

void Panda::draw(){
    glPushMatrix();
    glMultMatrixd(f.matrix());

    wrj.draw();
    marker.draw();
    navex.draw();
    effector.draw();

    QGLViewer::drawAxis(50.);

    //drawFrame(f);
    drawFrame(fsw);
    drawFrame(ftcp);

    glPopMatrix();
}

void Panda::drawFrame(Frame &frame){
    double size = 20;
    Vec centre = frame.localInverseCoordinatesOf(Vec(0,0,0));
    Vec x = frame.localInverseTransformOf(Vec(1,0,0))*size + centre;
    Vec y = frame.localInverseTransformOf(Vec(0,1,0))*size + centre;
    Vec z = frame.localInverseTransformOf(Vec(0,0,1))*size + centre;

    /*std::cout << "Coord x: " << x.x << "," << x.y << "," << x.z << std::endl;

    std::cout << "Coord y: " << y.x << "," << y.y << "," << y.z << std::endl;

    std::cout << "Coord : " << z.x << "," << z.y << "," << z.z << std::endl;*/

    glColor3f(1,0,0);
    glLineWidth(5.);
    glBegin(GL_LINES);
        glVertex3d(centre.x, centre.y, centre.z);
        glVertex3d(x.x, x.y, x.z);
    glEnd();

    glColor3f(0,1,0);
    glLineWidth(5.);
    glBegin(GL_LINES);
        glVertex3d(centre.x, centre.y, centre.z);
        glVertex3d(y.x, y.y, y.z);
    glEnd();

    glColor3f(0,0,1);
    glLineWidth(5.);
    glBegin(GL_LINES);
        glVertex3d(centre.x, centre.y, centre.z);
        glVertex3d(z.x, z.y, z.z);
    glEnd();
}

void Panda::openOFF(QString filename, SimpleMesh &mesh){
    std::vector<Vec3Df> &vertices = mesh.getVertices();
    std::vector<Triangle> &triangles = mesh.getTriangles();

    filename = QDir().relativeFilePath(filename);

    FileIO::openOFF(filename.toStdString(), vertices, triangles);

    mesh.update();
}

void Panda::setFrames(){
    fsw.setReferenceFrame(&f);
    ftcp.setReferenceFrame(&f);

    Vec pandaToSw[4] = {Vec(-1,0,0), Vec(0,-1,0), Vec(0,0,1), Vec(0,0,0)};

    Vec swToTCP[4] = {Vec(0.707106781186547, -0.707106781186548, 0), Vec(0.707106781186548, 0.707106781186547, -1.35585468084861E-31),
                     Vec(9.58734039131578E-32, 9.58734039131576E-32, 1), Vec(3.52273863634531, 3.52273863634531, 90)};

    // Create sw
    Quaternion qPandaSW;
    qPandaSW.setFromRotatedBasis(pandaToSw[0], pandaToSw[1], pandaToSw[2]);
    fsw.setPositionAndOrientation(pandaToSw[3], qPandaSW);

    // Get tcp in terms of the panda by passing through sw
    Vec tcpToPanda[4] = { fsw.localInverseTransformOf(swToTCP[0]), fsw.localInverseTransformOf(swToTCP[1]),
                          fsw.localInverseTransformOf(swToTCP[2]), fsw.localInverseCoordinatesOf(swToTCP[3]) };

    // Create tcp
    Quaternion qtcp;
    qtcp.setFromRotatedBasis(tcpToPanda[0], tcpToPanda[1], tcpToPanda[2]);
    ftcp.setPositionAndOrientation(tcpToPanda[3], qtcp);

    tcpPanda[0] = ftcp.localTransformOf(Vec(1,0,0));
    tcpPanda[1] = ftcp.localTransformOf(Vec(0,1,0));
    tcpPanda[2] = ftcp.localTransformOf(Vec(0,0,1));
    tcpPanda[3] = ftcp.localCoordinatesOf(Vec(0,0,0));
}

void Panda::checkPanda(){
    Vec x = fsw.localInverseTransformOf(Vec(1,0,0));
    Vec y = fsw.localInverseTransformOf(Vec(0,1,0));
    Vec z = fsw.localInverseTransformOf(Vec(0,0,1));

    std::cout << "Coord x: " << x.x << "," << x.y << "," << x.z << std::endl;
    std::cout << "Coord y: " << y.x << "," << y.y << "," << y.z << std::endl;
    std::cout << "Coord : " << z.x << "," << z.y << "," << z.z << std::endl;


    x = ftcp.localInverseTransformOf(Vec(1,0,0));
    y = ftcp.localInverseTransformOf(Vec(0,1,0));
    z = ftcp.localInverseTransformOf(Vec(0,0,1));

    std::cout << "Coord x: " << x.x << "," << x.y << "," << x.z << std::endl;
    std::cout << "Coord y: " << y.x << "," << y.y << "," << y.z << std::endl;
    std::cout << "Coord : " << z.x << "," << z.y << "," << z.z << std::endl;

}

void Panda::setLocation(const Vec &v){
    Frame ftemp;
    ftemp.setPositionAndOrientation(ftcp.position(), ftcp.rotation());
    ftemp.setPosition(v);

    Vec pos = ftemp.localInverseCoordinatesOf(tcpPanda[3]);
    f.setPosition(pos);
}

void Panda::setOrientation(const Quaternion &q){
    Frame ftemp;
    ftemp.setOrientation(q);

    Quaternion rotation;
    Vec x = ftemp.localInverseTransformOf(tcpPanda[0]);
    Vec y = ftemp.localInverseTransformOf(tcpPanda[1]);
    Vec z = ftemp.localInverseTransformOf(tcpPanda[2]);
    rotation.setFromRotatedBasis(x,y,z);

    f.setOrientation(rotation);
}

void Panda::setToPlane(const Vec &v, const Quaternion &q){
    Frame ftemp;
    ftemp.setPositionAndOrientation(v, q);      // set it to our goal
    ftemp.rotate(Quaternion(Vec(0,1,0), M_PI/2));       // rotate it to align to the correct axis

    Quaternion rotation;
    Vec x = ftemp.localInverseTransformOf(tcpPanda[0]);
    Vec y = ftemp.localInverseTransformOf(tcpPanda[1]);
    Vec z = ftemp.localInverseTransformOf(tcpPanda[2]);
    rotation.setFromRotatedBasis(x,y,z);

    Vec pos = ftemp.localInverseCoordinatesOf(tcpPanda[3]);

    f.setPositionAndOrientation(pos, rotation);
}

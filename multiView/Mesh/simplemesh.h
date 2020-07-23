#ifndef SIMPLEMESH_H
#define SIMPLEMESH_H

#include "Tools/vec3D.h"
#include "Tools/triangle.h"

class SimpleMesh
{

public:
    SimpleMesh(){}
    void draw();
    void recomputeNormals();
    void update();
    void clear();
    void invertNormal(){normalDirection *= -1;}

    std::vector<Vec3Df> &getVertices(){return vertices;}
    const std::vector<Vec3Df> &getVertices()const {return vertices;}

    std::vector<Triangle> &getTriangles(){return triangles;}
    const std::vector<Triangle> &getTriangles()const {return triangles;}

protected:
    Vec3Df computeTriangleNormal(unsigned int t);
    void computeVerticesNormals();
    void glTriangle(unsigned int i);

    std::vector <Vec3Df> vertices;      // starting verticies
    std::vector <Triangle> triangles;       // starting triangles
    std::vector<Vec3Df> verticesNormals;

    double normalDirection = 1.;
};

#endif // SIMPLEMESH_H

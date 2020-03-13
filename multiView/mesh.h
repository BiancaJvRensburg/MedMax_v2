#ifndef MESH_H
#define MESH_H

#include "Vec3D.h"
#include "Triangle.h"
#include "plane.h"
#include <queue>

enum Side {INTERIOR, EXTERIOR};

class Mesh : public QObject
{
    Q_OBJECT

public:

    Mesh():normalDirection(1.){}
    Mesh(std::vector<Vec3Df> &vertices, std::vector<Triangle> &triangles): vertices(vertices), triangles(triangles), normalDirection(1.){
        update();
    }
    ~Mesh(){}
    void init();
    Vec3Df& getBBCentre(){ return BBCentre; }

    std::vector<Vec3Df> &getVertices(){return vertices;}
    const std::vector<Vec3Df> &getVertices()const {return vertices;}

    std::vector<Triangle> &getTriangles(){return triangles;}
    const std::vector<Triangle> &getTriangles()const {return triangles;}

    std::vector< std::vector<unsigned int>> &getVertexNeighbours(){return vertexNeighbours;}
    const std::vector< std::vector<unsigned int>> &getVertexNeighbours()const {return vertexNeighbours;}

    std::vector< std::vector<unsigned int>> &getVertexTriangles(){return vertexTriangles;}
    const std::vector< std::vector<unsigned int>> &getVertexTriangles()const {return vertexTriangles;}

    void draw();

    void recomputeNormals();
    void update();
    void clear();

    float getBBRadius();

    void updatePlaneIntersections();    // need one for a single plane
    void updatePlaneIntersections(Plane *p);

    void addPlane(Plane *p);
    void deleteGhostPlanes();
    void setTransfer(bool isTransfer){ this->isTransfer = isTransfer; }
    void sendToManible();
    void setIsCut(Side s, bool isCut, bool isUpdate);
    void drawCut();
    bool getIsCut(){ return isCut; }

    typedef std::priority_queue< std::pair< float , int > , std::deque< std::pair< float , int > > , std::greater< std::pair< float , int > > > FacesQueue;

    void invertNormal(){normalDirection *= -1;}

    void readJSON(const QJsonObject &json, double &scale);

public Q_SLOTS:
    void recieveInfoFromFibula(std::vector<Vec>, std::vector<std::vector<int>>, std::vector<int>, std::vector<Vec>, int);

Q_SIGNALS:
    void sendInfoToManible(std::vector<int>, std::vector<Vec>, std::vector<std::vector<int>>, std::vector<int>, std::vector<Vec>, int);
    void updateViewer();

protected:
    void computeBB();

    void computeTriangleNormals();
    Vec3Df computeTriangleNormal(unsigned int t);
    void computeVerticesNormals();
    void glTriangle(unsigned int i);
    void glTriangleSmooth(unsigned int i);
    void glTriangleFibInMand(unsigned int i);
    void getColour(unsigned int vertex);

    void planeIntersection(unsigned int index);

    void floodNeighbour(unsigned int index, int id);     // flood the neighbours of the vertex index with the value id
    void mergeFlood();      // to be called after flooding; merges the regions between the planes

    void createSmoothedTriangles();
    void createSmoothedMandible();
    void createSmoothedFibula();

    void getSegmentsToKeep();   // Only for the fibula mesh (gets the segments between 2 planes that we want to keep)

    void cutMesh();
    void cutMandible(bool* truthTriangles);
    void cutFibula(bool* truthTriangles);
    void saveTrianglesToKeep(bool* truthTriangles, unsigned int i);
    void fillColours();

    void uniformScale(float s);

    std::vector <Vec3Df> vertices;      // starting verticies
    std::vector <Triangle> triangles;       // starting triangles
    std::vector <int> coloursIndicies;      // the value of each index allowing us to calculate its colour

    std::vector <Plane*> planes;
    std::vector <std::vector <unsigned int>> intersectionTriangles;    // Contains the index of the triangle instead of the actual triangle

    std::vector <int> flooding;
    std::vector< std::vector<unsigned int>> vertexNeighbours;       // each vertex's neighbours
    std::vector< std::vector<unsigned int>> vertexTriangles;        // the triangles each vertex belongs to
    std::vector<int> planeNeighbours;       // which planes are neighbours
    bool isCut = false;
    std::vector<unsigned int> trianglesCut;     // The list of triangles after the cutting (a list of triangle indicies)
    std::vector<unsigned int> trianglesExtracted;       // The list of triangles taken out (the complement of trianglesCut)
    std::vector<int> segmentsConserved; // filled with flooding values to keep

    std::vector<Vec3Df> smoothedVerticies;      // New verticies which line up with the cutting plane

    std::vector<Vec3Df> normals;
    std::vector<Vec3Df> verticesNormals;

    // The fibula in the manible
    std::vector<Vec3Df> fibInMandVerticies;
    std::vector<Triangle> fibInMandTriangles;
    std::vector<int> fibInMandColour;       // Only the fibula bones will be coloured
    std::vector<Vec3Df> fibInMandNormals;
    std::vector<Vec3Df> fibInMandVerticesNormals;
    int fibInMandNbColours;

    Side cuttingSide = Side::INTERIOR;

    Vec3Df BBMin;
    Vec3Df BBMax;
    Vec3Df BBCentre;
    float radius;

    bool isTransfer = true;

    int normalDirection;
};

#endif // MESH_H


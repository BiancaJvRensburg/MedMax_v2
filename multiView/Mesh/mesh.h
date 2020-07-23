#ifndef MESH_H
#define MESH_H

#include "Tools/vec3D.h"
#include "Tools/triangle.h"
#include "Planes/plane.h"
#include <queue>
#include <Eigenvalues>
#include <nanoflann.hpp>

typedef nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd> KDTree;

enum Side {INTERIOR, EXTERIOR};

class Mesh : public QObject
{
    Q_OBJECT

public:

    Mesh():normalDirection(1.),mat(0,3),tree(3, mat, 10){}
    Mesh(std::vector<Vec3Df> &vertices, std::vector<Triangle> &triangles):vertices(vertices), triangles(triangles), normalDirection(1.),mat(0,3),tree(3, mat, 10){
        update();
    }
    ~Mesh(){}
    void init();
    void computeBB(Vec3Df &centre, float& radius);

    std::vector<Vec3Df> &getVertices(){return vertices;}
    const std::vector<Vec3Df> &getVertices()const {return vertices;}

    std::vector<Triangle> &getTriangles(){return triangles;}
    const std::vector<Triangle> &getTriangles()const {return triangles;}

    std::vector<unsigned int> getVerticesOnPlane(unsigned int planeNb, Plane *p);
    Triangle& getTriangle(unsigned int i){ return triangles[i]; }
    Vec3Df& getSmoothVertex(unsigned int i){ return smoothedVerticies[i]; }

    void draw();
    void drawFragment(unsigned int i);

    void recomputeNormals();
    void update();
    void clear();

    void updatePlaneIntersections();    // need one for a single plane
    void updatePlaneIntersections(Plane *p);

    void addPlane(Plane *p);
    void deleteGhostPlanes();
    void setTransfer(bool isTransfer){ this->isTransfer = isTransfer; }
    void sendToMandible();
    void setIsCut(Side s, bool isCut, bool isUpdate);
    void drawCut();
    void drawCutMand();
    bool getIsCut(){ return isCut; }

    void setAlpha(float a){ alphaTransparency = a; }

    typedef std::priority_queue< std::pair< float , int > , std::deque< std::pair< float , int > > , std::greater< std::pair< float , int > > > FacesQueue;

    void invertNormal(){normalDirection *= -1;}

    void mlsProjection(const std::vector<Vec> &inputPoints, std::vector<Vec> &outputPoints);

    void getIntersectionForPlane(unsigned int index, std::vector <unsigned int> &intersectionTrianglesPlane);
    std::vector<Vec> getMinNormalForPlane(unsigned int index, unsigned int neighbourIndex);

public Q_SLOTS:
    void recieveInfoFromFibula(std::vector<int>&, std::vector<Vec>, std::vector<std::vector<int>>&, std::vector<int>&, std::vector<Vec>, int);

Q_SIGNALS:
    void sendInfoToManible(std::vector<int>&, std::vector<Vec>&, std::vector<std::vector<int>>&, std::vector<int>&, std::vector<Vec>&, int);
    void updateViewer();

protected:
    Vec3Df computeTriangleNormal(unsigned int t);
    void computeVerticesNormals();
    void glTriangle(unsigned int i);
    void glTriangleSmooth(unsigned int i, std::vector <int> &coloursIndicies);
    void glTriangleFragment(unsigned int index, unsigned int i, int colour);
    void getColour(unsigned int vertex, std::vector <int> &coloursIndicies);
    void collectOneRing(std::vector<std::vector<unsigned int>> &oneRing);
    void collectTriangleOneRing(std::vector<std::vector<unsigned int>> &oneTriangleRing);

    void planeIntersection(unsigned int index, std::vector <unsigned int> &intersectionTrianglesPlane, std::vector<bool>& isSharedVertex);
    void constructPlaneNeighbours(std::vector<int> &planeNeighbours);
    void setNeighbours(std::vector<int> &planeNeighbours, unsigned int a, unsigned int b);

    void floodNeighbour(unsigned int index, int id, std::vector<int> &planeNeighbours, std::vector<bool>& isSharedVertex);     // flood the neighbours of the vertex index with the value id
    void mergeFlood(const std::vector<int> &planeNeighbours);      // to be called after flooding; merges the regions between the planes

    void createSmoothedTriangles(std::vector <std::vector <unsigned int>> &intersectionTriangles, const std::vector<int> &planeNeighbours);
    void createSmoothedMandible(std::vector <std::vector <unsigned int>> &intersectionTriangles, const std::vector<int> &planeNeighbours);
    void createSmoothedFibula(std::vector <std::vector <unsigned int>> &intersectionTriangles, const std::vector<int> &planeNeighbours);
    void getSegmentsToKeep(const std::vector<int> &planeNeighbours);   // Only for the fibula mesh (gets the segments between 2 planes that we want to keep)

    void cutMesh(std::vector <std::vector <unsigned int>> &intersectionTriangles, const std::vector<int> &planeNeighbours);
    void cutMandible(bool* truthTriangles, const std::vector<int> &planeNeighbours);
    void cutFibula(bool* truthTriangles, std::vector <std::vector <unsigned int>> &intersectionTriangles, const std::vector<int> &planeNeighbours);
    void saveTrianglesToKeep(bool* truthTriangles, unsigned int i);
    void fillColours(std::vector <int> &coloursIndicies, const unsigned long long nbColours);

    Vec getPolylineProjectedVertex(unsigned int p1, unsigned int p2, unsigned int vertexIndex);

    Vec3Df& getVertex(unsigned int i){ return vertices[i]; }

    std::vector <Vec3Df> vertices;      // starting verticies
    std::vector <Triangle> triangles;       // starting triangles
    std::vector <Plane*> planes;

    std::vector<std::vector<unsigned int>> oneRing;
    std::vector<std::vector<unsigned int>> oneTriangleRing;
    std::vector<int> flooding;

    bool isCut = false;
    std::vector<unsigned int> trianglesCut;     // The list of triangles after the cutting (a list of triangle indicies)
    std::vector<unsigned int> trianglesExtracted;       // The list of triangles taken out (the complement of trianglesCut)
    std::vector<int> segmentsConserved; // filled with flooding values to keep

    std::vector<Vec3Df> smoothedVerticies;      // New verticies which line up with the cutting plane
    std::vector<Vec3Df> verticesNormals;

    // The fibula in the manible
    std::vector<Vec3Df> fibInMandVerticies;
    std::vector<int> fibInMandColour;       // Only the fibula bones will be coloured
    std::vector<Vec3Df> fibInMandNormals;
    int fibInMandNbColours;
    std::vector<std::vector<Triangle>> fibulaFragmentTriangles;

    Side cuttingSide = Side::INTERIOR;
    bool isTransfer = true;

    int normalDirection;
    float alphaTransparency = 1.f;
    Vec3Df BBMin;
    Vec3Df BBMax;

    Eigen::MatrixXd pointsToMatrix(const std::vector<Vec3Df> &basePoints, const int dimension);
    Eigen::MatrixXd mat;
    void HPSS(KDTree &tree, const int dimension, Vec inputPoint, Vec &outputPoint , Vec &outputNormal, double radius, unsigned int nbIterations=10, unsigned int knn=20, double s=1.0);
    void weightGauss(Vec x, std::vector<double> &weights, unsigned int knn, std::vector<size_t> const &id_nearest_neighbors, std::vector<double> const &square_distances_to_neighbors, double h);

    KDTree tree;
};

#endif // MESH_H


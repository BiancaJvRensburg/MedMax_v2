#include "mesh.h"
#include <algorithm>
#include <float.h>

void Mesh::init(){
    collectOneRing(oneRing);
    collectTriangleOneRing(oneTriangleRing);
    update();
}

void Mesh::computeBB(Vec3Df &centre, float &radius){

    Vec3Df BBMin( FLT_MAX, FLT_MAX, FLT_MAX );
    Vec3Df BBMax( -FLT_MAX, -FLT_MAX, -FLT_MAX );

    for( unsigned int i = 0 ; i < vertices.size() ; i ++ ){
        const Vec3Df &point = vertices[i];
        for( int v = 0 ; v < 3 ; v++ ){
            float value = point[v];
            if( BBMin[v] > value ) BBMin[v] = value;
            if( BBMax[v] < value ) BBMax[v] = value;
        }
    }

    radius = (BBMax - BBMin).norm() / 2.0f;
    centre = (BBMax + BBMin)/2.0f;
}

void Mesh::collectOneRing(std::vector<std::vector<unsigned int>> &oneRing) {
    oneRing.resize(vertices.size());

    for (unsigned int i = 0; i < triangles.size (); i++) {
        const Triangle &ti = triangles[i];
        for (unsigned int j = 0; j < 3; j++) {
            unsigned int vj = ti.getVertex(j);
            for (unsigned int k = 1; k < 3; k++) {
                unsigned int vk = ti.getVertex((j+k)%3);
                if (std::find (oneRing[vj].begin (), oneRing[vj].end (), vk) == oneRing[vj].end ())
                    oneRing[vj].push_back(vk);
            }
        }
    }
}

void Mesh::collectTriangleOneRing(std::vector<std::vector<unsigned int> > &oneTriangleRing){
    oneTriangleRing.resize(vertices.size());

    for (unsigned int i = 0; i < triangles.size (); i++) {
        const Triangle &t = triangles[i];
        for (unsigned int j = 0; j < 3; j++) {
            unsigned int v = t.getVertex(j);
            oneTriangleRing[v].push_back(i);
        }
    }
}

void Mesh::update(){
    recomputeNormals();
    updatePlaneIntersections();
}

void Mesh::clear(){
    vertices.clear();
    triangles.clear();
    verticesNormals.clear();
}

void Mesh::recomputeNormals () {
    computeVerticesNormals();
}

Vec3Df Mesh::computeTriangleNormal(unsigned int id ){
    const Triangle &t = triangles[id];
    Vec3Df normal = Vec3Df::crossProduct(vertices[t.getVertex (1)] - vertices[t.getVertex (0)], vertices[t.getVertex (2)]- vertices[t.getVertex (0)]);
    normal.normalize();
    return normal;
}

void Mesh::setIsCut(Side s, bool isCut, bool isUpdate){
    this->isCut = isCut;
    this->cuttingSide = s;
    if(!isCut) deleteGhostPlanes();
    if(isUpdate) updatePlaneIntersections();
}

void Mesh::computeVerticesNormals(){
    verticesNormals.clear();
    verticesNormals.resize( vertices.size() , Vec3Df(0.,0.,0.) );

    for( unsigned int t = 0 ; t<triangles.size(); ++t ){
        Vec3Df const &tri_normal = computeTriangleNormal(t);
        verticesNormals[ triangles[t].getVertex(0) ] += tri_normal;
        verticesNormals[ triangles[t].getVertex(1) ] += tri_normal;
        verticesNormals[ triangles[t].getVertex(2) ] += tri_normal;
    }

    for( unsigned int v = 0 ; v < verticesNormals.size() ; ++v ) verticesNormals[ v ].normalize();
}

// Access and colour each individual vertex here
void Mesh::glTriangle(unsigned int i){
    const Triangle &t = triangles[i];

    for(unsigned int j = 0 ; j < 3 ; j++ ){
        glNormal(verticesNormals[t.getVertex(j)]*normalDirection);
        glVertex(vertices[t.getVertex(j)]);
    }

    glColor4f(1.0, 1.0, 1.0, alphaTransparency);
}

void Mesh::glTriangleSmooth(unsigned int i, std::vector <int> &coloursIndicies){
    const Triangle &t = triangles[i];

    for(unsigned int j = 0 ; j < 3 ; j++ ){
        if(cuttingSide==Side::EXTERIOR) getColour(t.getVertex(j), coloursIndicies);
        glNormal(verticesNormals[t.getVertex(j)]*normalDirection);
        glVertex(smoothedVerticies[t.getVertex(j)]);
    }

    glColor4f(1.0, 1.0, 1.0, alphaTransparency);
}

void Mesh::glTriangleFibInMand(unsigned int i, std::vector <int> &coloursIndicies){
    const Triangle &t = fibInMandTriangles[i];

    for(unsigned int j=0 ; j<3 ; j++){
        getColour(t.getVertex(j), coloursIndicies);
        glNormal(fibInMandNormals[t.getVertex(j)]*normalDirection);
        glVertex(fibInMandVerticies[t.getVertex(j)]);
    }

    glColor4f(1.0, 1.0, 1.0, alphaTransparency);
}

void Mesh::getColour(unsigned int vertex, std::vector <int> &coloursIndicies){
    int colour;
    float nb;

    if(cuttingSide==Side::EXTERIOR && coloursIndicies.size()!=0){ // TODO protect this higher up
        colour = coloursIndicies[vertex];
        nb = static_cast<float>(planes.size())/2.f;
    }
    else if(cuttingSide==Side::INTERIOR && fibInMandColour.size()!=0){
        colour = fibInMandColour[vertex];
        nb = static_cast<float>(fibInMandNbColours);
    }
    else return;

    if(colour != -1){
        float r,g,b;
        float c = static_cast<float>(colour) + 1.f;
        r = c/nb;
        g = c/nb + (1.f/3.f);
        b = c/nb + (2.f/3.f);

        while(r>1.f) r -= 1.f;
        while(g>1.f) g -= 1.f;
        while(b>1.f) b -= 1.f;

        glColor4f(r,g,b, alphaTransparency);
    }
    else glColor4f(1.0, 1.0, 1.0, alphaTransparency);
}

void Mesh::addPlane(Plane *p){
    planes.push_back(p);
    updatePlaneIntersections(p);
}

void Mesh::deleteGhostPlanes(){
    if(planes.size()==2) return;        // Keep the left and right planes
    planes.erase(planes.begin()+2, planes.end());       // delete the ghost planes
}

void Mesh::updatePlaneIntersections(){
    if(isCut){
        std::vector <std::vector <unsigned int>> intersectionTriangles;
        intersectionTriangles.resize(planes.size());

        flooding.clear();
        for(unsigned int i=0; i<vertices.size(); i++) flooding.push_back(-1);       // reset the flooding values

        std::vector<int> planeNeighbours;
        for(unsigned int i=0; i<planes.size()*2; i++) planeNeighbours.push_back(-1);
        for(unsigned int i=0; i<planes.size(); i++) planeIntersection(i, intersectionTriangles[i]);

        for(unsigned int i=0; i<intersectionTriangles.size(); i++){
            std::vector<unsigned int> &triIndexes = intersectionTriangles[i];
            for(unsigned int k=0; k<triIndexes.size(); k++){
                for(unsigned int l=0; l<3; l++){
                    Triangle &t = triangles[triIndexes[k]];
                    unsigned int index = t.getVertex(l);
                    for(unsigned int j=0; j<oneRing[index].size(); j++){
                        floodNeighbour(oneRing[index][j], flooding[index], planeNeighbours);
                    }
                }
            }
        }

        mergeFlood(planeNeighbours);
        cutMesh(intersectionTriangles, planeNeighbours);
    }
}

void Mesh::cutMesh(std::vector<std::vector<unsigned int>> &intersectionTriangles, const std::vector<int> &planeNeighbours){
    trianglesCut.clear();

    bool truthTriangles[triangles.size()];  // keeps a record of the triangles who are already added
    for(unsigned int i=0; i<triangles.size(); i++) truthTriangles[i] = false;

    switch (cuttingSide) {
        case Side::INTERIOR:        // MANDIBLE
            cutMandible(truthTriangles, planeNeighbours);
        break;

        case Side::EXTERIOR:        // FIBULA
            cutFibula(truthTriangles, intersectionTriangles, planeNeighbours);
        break;
    }

    trianglesExtracted.clear();     // Store the rest of the triangles
    for(unsigned int i=0; i<triangles.size(); i++){
        if(!truthTriangles[i]) trianglesExtracted.push_back(i);
    }

    // ! Conserve this order
    createSmoothedTriangles(intersectionTriangles, planeNeighbours);

    if(cuttingSide == Side::EXTERIOR){      // send the segments to the mandible
        if(isTransfer){
            sendToMandible();
        }
    }
}

void Mesh::cutMandible(bool* truthTriangles, const std::vector<int> &planeNeighbours){
    for(unsigned int i=0; i<flooding.size(); i++){
        int flood = flooding[i];
        if(flood == -1) continue;
        int neighbour = planeNeighbours[static_cast<unsigned int>(flood)];
        if(neighbour ==-1){
            saveTrianglesToKeep(truthTriangles, i);
        }
    }
}

void Mesh::cutFibula(bool* truthTriangles, std::vector <std::vector <unsigned int>> &intersectionTriangles, const std::vector<int> &planeNeighbours){
    for(unsigned int j=0; j<intersectionTriangles.size(); j++){
            const std::vector<unsigned int> &v = intersectionTriangles[j];
            for(unsigned int k=0; k<v.size(); k++) {
                truthTriangles[v[k]] = true;
        }
    }

    getSegmentsToKeep(planeNeighbours);    // figure out what to keep (TODO can be done earlier)
    for(unsigned int i=0; i<flooding.size(); i++){
        bool isKeep = false;
        for(unsigned int k=0; k<segmentsConserved.size(); k++){      // Only keep it if it belongs to a kept segment
            if(segmentsConserved[k]==flooding[i]){
                isKeep = true;
                break;
            }
        }
        if(isKeep){
            for(unsigned int j=0; j<oneTriangleRing[i].size(); j++){        // Get the triangles they belong to
                saveTrianglesToKeep(truthTriangles, i);
            }
        }
    }
}

void Mesh::saveTrianglesToKeep(bool* truthTriangles, unsigned int i){
    for(unsigned int j=0; j<oneTriangleRing[i].size(); j++){        // Get the triangles the indicies belong to
        if(!truthTriangles[oneTriangleRing[i][j]]){     // If it's not already in the list
            trianglesCut.push_back(oneTriangleRing[i][j]);
            truthTriangles[oneTriangleRing[i][j]] = true;
        }
    }
}

void Mesh::fillColours(std::vector <int> &coloursIndicies, const unsigned long long nbColours){
    int tempColours[nbColours];
    coloursIndicies.clear();

    for(unsigned int i=0; i<nbColours; i++) tempColours[i] = -1;       // init all to -1

    for(unsigned int i=0; i<segmentsConserved.size(); i++) tempColours[segmentsConserved[i]] = static_cast<int>(i);       // change to the seg colours value

    for(unsigned int i=0; i<vertices.size(); i++){
        int index = flooding[i]; 
        coloursIndicies.push_back(tempColours[index]);      // Fill the colours
    }
}

// WARNING : this assumes that the left and right planes are the first planes added!
// Could search for the exterior planes beforehand using the fact that the other sides = -1
void Mesh::getSegmentsToKeep(const std::vector<int> &planeNeighbours){
    segmentsConserved.clear();

    // Find the non-discarded side of the left plane
    int planeToKeep;
    if(planeNeighbours[0]!=-1) planeToKeep = 0; // if it has a neighbour
    else planeToKeep = static_cast<int>(planes.size());   // keep the otherside if 0 is discared

    // if there are no ghost planes
    if(planes.size()==2){
        int rightPlaneKept;
        if(planeNeighbours[1]!=-1) rightPlaneKept = 1;
        else rightPlaneKept = 3;

        // Keep the smallest
        int max;
        if(planeToKeep<rightPlaneKept) max = planeToKeep;
        else max = rightPlaneKept;
        segmentsConserved.push_back(max);
        return;
    }

    // while we haven't found the right plane
    while(planeToKeep!=1 && planeToKeep!=static_cast<int>(planes.size())+1){
        int nextPlane = planeNeighbours[static_cast<unsigned int>(planeToKeep)];   // move on to the next plane

        // Keep the smaller of the two values to match the merge flood
        if(planeToKeep < nextPlane) segmentsConserved.push_back(planeToKeep);
        else segmentsConserved.push_back(nextPlane);

        // discard the other side
        int toDiscard;
        if( nextPlane < static_cast<int>(planes.size()) ) toDiscard = nextPlane + static_cast<int>(planes.size());
        else toDiscard = nextPlane - static_cast<int>(planes.size());

        if(toDiscard==1 || toDiscard==static_cast<int>(planes.size())+1) break;

        // move on to the next plane
        nextPlane = planeNeighbours[static_cast<unsigned int>(toDiscard)];

        // keep the other side
        if( nextPlane < static_cast<int>(planes.size()) ) planeToKeep = nextPlane + static_cast<int>(planes.size());
        else planeToKeep = nextPlane - static_cast<int>(planes.size());
    }
}

void Mesh::createSmoothedTriangles(std::vector <std::vector <unsigned int>> &intersectionTriangles, const std::vector<int> &planeNeighbours){
    smoothedVerticies.clear();

    for(unsigned int i=0; i<vertices.size(); i++) smoothedVerticies.push_back(vertices[i]);  // Copy the verticies table

    switch (cuttingSide) {
        case Side::INTERIOR:
            createSmoothedMandible(intersectionTriangles, planeNeighbours);
        break;

        case Side::EXTERIOR:
            createSmoothedFibula(intersectionTriangles, planeNeighbours);
        break;
    }
}

void Mesh::createSmoothedMandible(std::vector <std::vector <unsigned int>> &intersectionTriangles, const std::vector<int> &planeNeighbours){
    for(unsigned long long i=0; i<planes.size(); i++){
        for(unsigned long long j=0; j<intersectionTriangles[static_cast<unsigned long long>(i)].size(); j++){       // for each triangle cut
            for(unsigned int k=0; k<3; k++){    // find which verticies to keep
                const unsigned int &vertexIndex = triangles[intersectionTriangles[i][j]].getVertex(k);
                if(planeNeighbours[static_cast<unsigned int>(flooding[vertexIndex])] != -1){   // if we need to change it (here we only change it if it's outside of the cut (fine for mandible))
                    Vec newVertex = planes[i]->getProjection(Vec(static_cast<double>(vertices[vertexIndex][0]), static_cast<double>(vertices[vertexIndex][1]), static_cast<double>(vertices[vertexIndex][2])) );
                    smoothedVerticies[vertexIndex] = Vec3Df(static_cast<float>(newVertex.x), static_cast<float>(newVertex.y), static_cast<float>(newVertex.z)); // get the projection
                }
                // else don't change the original
            }

        }
    }
}

void Mesh::createSmoothedFibula(std::vector <std::vector <unsigned int>> &intersectionTriangles, const std::vector<int> &planeNeighbours){
    for(unsigned int i=0; i<planes.size(); i++){
        //verticesOnPlane[i].clear();
        for(unsigned long long j=0; j<intersectionTriangles[static_cast<unsigned long long>(i)].size(); j++){   // for each triangle cut
            int actualFlooding = -1;    //  Conserve the "real" flooding value (will never stay at -1)

            for(unsigned int k=0; k<3; k++){    // find which verticies are on the otherside of the cut
                const unsigned int &vertexIndex = triangles[intersectionTriangles[i][j]].getVertex(k);

                bool isOutlier = false;
                for(unsigned int l=0; l<segmentsConserved.size(); l++){
                    if(flooding[vertexIndex] == segmentsConserved[l]){
                        actualFlooding = flooding[vertexIndex];
                        isOutlier = true;
                    }
                }

                if(planeNeighbours[static_cast<unsigned int>(flooding[vertexIndex])]==-1 || isOutlier){        // if we need to change it
                    Vec newVertex;
                    const unsigned int &lastIndex = static_cast<unsigned int>(planes.size()-1);
                    if(i>2 && i<lastIndex){
                        if(i%2==0) newVertex = getPolylineProjectedVertex(i, i-1, vertexIndex);
                        else newVertex = getPolylineProjectedVertex(i, i+1, vertexIndex);
                    }
                    else if(i==2) newVertex = getPolylineProjectedVertex(i, 0, vertexIndex);
                    else if(i==lastIndex && lastIndex!=1) newVertex = getPolylineProjectedVertex(i, 1, vertexIndex);
                    else if(i==0){
                        if(lastIndex>1) newVertex = getPolylineProjectedVertex(i, 2, vertexIndex);
                        else newVertex = getPolylineProjectedVertex(i, 1, vertexIndex);
                    }
                    else if(i==1){
                        if(lastIndex>1) newVertex = getPolylineProjectedVertex(i, lastIndex, vertexIndex);
                        else newVertex = getPolylineProjectedVertex(i, 0, vertexIndex);
                    }
                    smoothedVerticies[vertexIndex] = Vec3Df(static_cast<float>(newVertex.x), static_cast<float>(newVertex.y), static_cast<float>(newVertex.z)); // get the projection
                }
                // else don't change the original
            }

            // Set the whole triangle to the correct flooding value
            for(unsigned int k=0; k<3; k++){
                const unsigned int &vertexIndex = triangles[intersectionTriangles[i][j]].getVertex(k);
                flooding[vertexIndex] = actualFlooding;
            }
        }
    }
}

Vec Mesh::getPolylineProjectedVertex(unsigned int p1, unsigned int p2, unsigned int vertexIndex){
    Vec n = planes[p2]->getPosition() - planes[p1]->getPosition();
    n.normalize();
    n = planes[p1]->getLocalVector(n);
    Vec p = Vec(static_cast<double>(vertices[vertexIndex][0]), static_cast<double>(vertices[vertexIndex][1]), static_cast<double>(vertices[vertexIndex][2]));
    p = planes[p1]->getLocalCoordinates(p);
    double alpha = p.z / n.z;
    Vec newVertex = p - alpha*n;
    return planes[p1]->getMeshCoordinatesFromLocal(newVertex);
}

void Mesh::updatePlaneIntersections(Plane *p){
    // Possible optimisation?

    updatePlaneIntersections();
}

/*void Mesh::floodNeighbour(unsigned int index, int id, std::vector<int> &planeNeighbours){

    if(flooding[index] == -1){      // Flood it
        flooding[index] = id;
        for(unsigned int i=0; i<oneRing[index].size(); i++){
            floodNeighbour(oneRing[index][i], id, planeNeighbours);
        }
    }

    else if(flooding[index] == id) return;      // stop if the vertex is already flooded with the same value

    else if(flooding[index]==id+static_cast<int>(planes.size()) || id==flooding[index]+static_cast<int>(planes.size())) return;     // stop if we've found our own neg/pos side

    else{       // else it already belongs to a different plane
        if(planeNeighbours[static_cast<unsigned int>(id)]== -1){       // They're not already neighbours
            planeNeighbours[static_cast<unsigned int>(id)] = flooding[index];     // equal to the old value
            planeNeighbours[static_cast<unsigned int>(flooding[index])] = id;
        }
        return;     // They're already neighbours
    }
}*/

void Mesh::floodNeighbour(unsigned int index, int id, std::vector<int> &planeNeighbours){
    std::queue<unsigned int> toFlood;
    toFlood.push(index);

    while(toFlood.size()!=0){
        index = toFlood.front();
        int flood = flooding[index];
        toFlood.pop();
        if(flood == -1){      // Flood it
            flooding[index] = id;
            for(unsigned int i=0; i<oneRing[index].size(); i++){
                toFlood.push(oneRing[index][i]);
            }
        }

        else if(flood == id) continue;      // stop if the vertex is already flooded with the same value

        else if(flood==id+static_cast<int>(planes.size()) || id==flood+static_cast<int>(planes.size())) continue;     // stop if we've found our own neg/pos side

        else{       // else it already belongs to a different plane
            if(planeNeighbours[static_cast<unsigned int>(id)]== -1){       // They're not already neighbours
                planeNeighbours[static_cast<unsigned int>(id)] = flooding[index];     // equal to the old value
                planeNeighbours[static_cast<unsigned int>(flooding[index])] = id;
            }
            //return;     // They're already neighbours
        }
    }

}

void Mesh::mergeFlood(const std::vector<int> &planeNeighbours){
    for(unsigned int i=0; i<flooding.size(); i++){
        int flood = flooding[i];
        if(flood != -1){
            int neighbour = planeNeighbours[static_cast<unsigned int>(flood)];
            if(neighbour != -1 && neighbour < flood){     // From the two neighbours, set them both to the lowest value
                flooding[i] = neighbour;
            }
        }
    }
}

// Finds all the intersecting triangles for plane nb index
void Mesh::planeIntersection(unsigned int index, std::vector <unsigned int> &intersectionTrianglesPlane){
    intersectionTrianglesPlane.clear();       // empty the list of intersections

    for(unsigned int i = 0 ; i < triangles.size(); i++){
        const unsigned int &t0 = triangles[i].getVertex(0);
        const unsigned int &t1 = triangles[i].getVertex(1);
        const unsigned int &t2 = triangles[i].getVertex(2);

        if(planes[index]->isIntersection(Vec(vertices[t0]), Vec(vertices[t1]), Vec(vertices[t2]) )){        // if the triangle intersects the plane
            intersectionTrianglesPlane.push_back(i);      // save the triangle index

            // For each vertex, get the apporiate sign
            for(unsigned int j=0; j<3; j++){
                double sign = planes[index]->getSign(Vec(vertices[triangles[i].getVertex(j)]));     // get which side of the plane the vertex is on
                if(sign >= 0 ){
                    flooding[triangles[i].getVertex(j)] = static_cast<int>(planes.size() + index);
                }
                else if(sign < 0){
                    flooding[triangles[i].getVertex(j)] =  static_cast<int>(index);
                }
            }
        }
    }
}

void Mesh::getIntersectionForPlane(unsigned int index, std::vector<unsigned int> &intersectionTrianglesPlane){
    for(unsigned int i = 0 ; i < triangles.size(); i++){
        const unsigned int &t0 = triangles[i].getVertex(0);
        const unsigned int &t1 = triangles[i].getVertex(1);
        const unsigned int &t2 = triangles[i].getVertex(2);

        if(planes[index]->isIntersection(Vec(vertices[t0]), Vec(vertices[t1]), Vec(vertices[t2]) ))        // if the triangle intersects the plane
            intersectionTrianglesPlane.push_back(i);
    }
}

void Mesh::sendToMandible(){
    std::vector<int> planeNb;       // the plane nb associated
    std::vector<Vec> convertedVerticies;    // the vertex coordinates in relation to the plane nb
    std::vector<std::vector<int>>convertedTriangles; // the new indicies of the triangles (3 indicies)
    std::vector<int> convertedColours;
    std::vector<Vec> convertedNormals;
    int tempVerticies[smoothedVerticies.size()];   // a temporary marker for already converted verticies

    std::vector <int> coloursIndicies;
    fillColours(coloursIndicies, planes.size()*2);

    for(unsigned int i=0; i<smoothedVerticies.size(); i++) tempVerticies[i] = -1;

    for(unsigned int i=0; i<trianglesCut.size(); i++){      // For every triangle we want to send (we've already filtered out the rest when cutting the mesh)
        std::vector<int> newTriangle;

        for(unsigned int j=0; j<3; j++){
            const unsigned int &triVert = triangles[trianglesCut[i]].getVertex(j);       // this must go into smoothedVerticies[triV....] (ie. index in the original)

            if(tempVerticies[triVert] != -1) newTriangle.push_back(tempVerticies[triVert]);     // If converted already

            else{       // convert to the corresponding plane
                int pNb = flooding[triVert];    // Get the plane nb
                if(pNb >= static_cast<int>(planes.size())) pNb -= planes.size();      // The plane nb is referenced by the smallest side
                planeNb.push_back(pNb);

                // Convert the vertex
                const Vec &unConverted = Vec(static_cast<double>(smoothedVerticies[triVert][0]), static_cast<double>(smoothedVerticies[triVert][1]), static_cast<double>(smoothedVerticies[triVert][2]));
                const Vec &convertedCoords = planes[static_cast<unsigned int>(pNb)]->getLocalCoordinates(unConverted);

                convertedVerticies.push_back(convertedCoords);      // add the converted vertex triVert
                convertedColours.push_back(coloursIndicies[triVert]); // add its corresponding colour and normal

                const Vec &unConvertedNorm = Vec(static_cast<double>(verticesNormals[triVert][0]), static_cast<double>(verticesNormals[triVert][1]), static_cast<double>(verticesNormals[triVert][2]));
                const Vec &convertedNorm = planes[static_cast<unsigned int>(pNb)]->getLocalVector(unConvertedNorm);
                convertedNormals.push_back(convertedNorm);

                const int &vertexIndex = static_cast<int>(convertedVerticies.size()) - 1;
                newTriangle.push_back(vertexIndex);  // set it to the last index of convertedVerticies

                tempVerticies[triVert] = vertexIndex;       // Store the corresponding index in tempVerticies
            }
        }      
        convertedTriangles.push_back(newTriangle);      // Add the triangle
    }

    Q_EMIT sendInfoToManible(planeNb, convertedVerticies, convertedTriangles, convertedColours, convertedNormals, (static_cast<int>(planes.size())/2));
}

void Mesh::recieveInfoFromFibula(const std::vector<Vec> &convertedVerticies, const std::vector<std::vector<int>> &convertedTriangles, const std::vector<int> &convertedColours, const std::vector<Vec> &convertedNormals, int nbColours){
    if(cuttingSide != Side::INTERIOR) return;

    fibInMandTriangles.clear();
    fibInMandVerticies.clear();
    fibInMandColour.clear();
    fibInMandNormals.clear();

    fibInMandNbColours = nbColours;

    for(unsigned int i=0; i<convertedVerticies.size(); i++){
        Vec3Df v = Vec3Df(static_cast<float>(convertedVerticies[i].x), static_cast<float>(convertedVerticies[i].y), static_cast<float>(convertedVerticies[i].z));
        fibInMandVerticies.push_back(v);
        fibInMandColour.push_back(convertedColours[i]);
        v = Vec3Df(static_cast<float>(convertedNormals[i].x), static_cast<float>(convertedNormals[i].y), static_cast<float>(convertedNormals[i].z));
        fibInMandNormals.push_back(v);
    }

    for(unsigned int i=0; i<convertedTriangles.size(); i++){
        const Triangle &t = Triangle(static_cast<unsigned int>(convertedTriangles[i][0]), static_cast<unsigned int>(convertedTriangles[i][1]), static_cast<unsigned int>(convertedTriangles[i][2]));
        fibInMandTriangles.push_back(t);
    }

    Q_EMIT updateViewer();
}

void Mesh::draw()
{

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_DEPTH);

    glBegin (GL_TRIANGLES);

    if(!isCut){
        for(unsigned int i=0; i<triangles.size(); i++){
            glTriangle(i);
        }
    }
    else{
        std::vector<int> coloursIndicies;
        fillColours(coloursIndicies, planes.size()*2);
        for(unsigned int i=0; i<trianglesCut.size(); i++){
            glTriangleSmooth(trianglesCut[i], coloursIndicies);
        }

        for(unsigned int i=0; i<fibInMandTriangles.size(); i++){
            glTriangleFibInMand(i, coloursIndicies);
        }
    }

    glEnd();

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_DEPTH);
}

void Mesh::drawCut(){
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_DEPTH);

    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

    glBegin (GL_TRIANGLES);
    std::vector<int> coloursIndicies;
    fillColours(coloursIndicies, planes.size()*2);

    for(unsigned int i = 0 ; i < trianglesExtracted.size(); i++){
        glTriangleSmooth(trianglesExtracted[i], coloursIndicies);
    }

    glEnd();

    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_DEPTH);
}

std::vector<unsigned int> Mesh::getVerticesOnPlane(unsigned int planeNb, Plane *p){
    std::vector<unsigned int> v;
    std::vector<unsigned int> intersectionTrianglesPlane;
    getIntersectionForPlane(planeNb, intersectionTrianglesPlane);

    for(unsigned int i=0; i<intersectionTrianglesPlane.size(); i++){
        for(unsigned int k=0; k<3; k++){
            const unsigned int &triangleNb = intersectionTrianglesPlane[i];
            const unsigned int &index = triangles[triangleNb].getVertex(k);
            if(abs(p->getLocalCoordinates(Vec(smoothedVerticies[index])).z) < 0.001){
                bool isFound = false;
                for(unsigned int j=0; j<v.size(); j++){
                    if(v[j]==index){
                        isFound = true;
                        break;
                    }
                }
                if(!isFound) v.push_back(index);
            }
        }
    }
    return v;
}

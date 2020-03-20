#include "mesh.h"
#include <algorithm>
#include <float.h>

void Mesh::computeBB(){

    BBMin = Vec3Df( FLT_MAX, FLT_MAX, FLT_MAX );
    BBMax = Vec3Df( -FLT_MAX, -FLT_MAX, -FLT_MAX );

    for( unsigned int i = 0 ; i < vertices.size() ; i ++ ){
        const Vec3Df & point = vertices[i];
        for( int v = 0 ; v < 3 ; v++ ){
            float value = point[v];
            if( BBMin[v] > value ) BBMin[v] = value;
            if( BBMax[v] < value ) BBMax[v] = value;
        }
    }

    radius = (BBMax - BBMin).norm() / 2.0f;

    BBCentre = (BBMax + BBMin)/2.0f;
}

void Mesh::update(){
    computeBB();
    recomputeNormals();
    updatePlaneIntersections();
}

void Mesh::clear(){
    vertices.clear();
    triangles.clear();
    normals.clear();
    verticesNormals.clear();
}

void Mesh::recomputeNormals () {
    computeTriangleNormals();
    computeVerticesNormals();
}

void Mesh::computeTriangleNormals(){
    normals.clear();

    for(unsigned int i = 0 ; i < triangles.size() ; i++){
        normals.push_back(computeTriangleNormal(i));
    }
}

Vec3Df Mesh::computeTriangleNormal(unsigned int id ){
    const Triangle & t = triangles[id];
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

    for( unsigned int t = 0 ; t < triangles.size(); ++t )
    {
        Vec3Df const & tri_normal = normals[t];
        verticesNormals[ triangles[t].getVertex(0) ] += tri_normal;
        verticesNormals[ triangles[t].getVertex(1) ] += tri_normal;
        verticesNormals[ triangles[t].getVertex(2) ] += tri_normal;
    }

    for( unsigned int v = 0 ; v < verticesNormals.size() ; ++v )
    {
        verticesNormals[ v ].normalize();
    }
}

// Access and colour each individual vertex here
void Mesh::glTriangle(unsigned int i){
    const Triangle & t = triangles[i];

    for(unsigned int j = 0 ; j < 3 ; j++ ){
        glNormal(verticesNormals[t.getVertex(j)]*normalDirection);
        glVertex(vertices[t.getVertex(j)]);
    }

    glColor4f(1.0, 1.0, 1.0, alphaTransparency);
}

void Mesh::glTriangleSmooth(unsigned int i){
    const Triangle & t = triangles[i];

    for(unsigned int j = 0 ; j < 3 ; j++ ){
        if(cuttingSide==Side::EXTERIOR) getColour(t.getVertex(j));
        glNormal(verticesNormals[t.getVertex(j)]*normalDirection);
        glVertex(smoothedVerticies[t.getVertex(j)]);
    }

    glColor4f(1.0, 1.0, 1.0, alphaTransparency);
}

void Mesh::glTriangleFibInMand(unsigned int i){
    const Triangle & t = fibInMandTriangles[i];

    for(unsigned int j = 0 ; j < 3 ; j++ ){
        getColour(t.getVertex(j));
        glNormal(fibInMandNormals[t.getVertex(j)]*normalDirection);
        glVertex(fibInMandVerticies[t.getVertex(j)]);
    }

    glColor4f(1.0, 1.0, 1.0, alphaTransparency);
}

void Mesh::getColour(unsigned int vertex){
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
    std::vector<unsigned int> init;
    intersectionTriangles.push_back(init);
    //verticesOnPlane.push_back(init);
    updatePlaneIntersections(p);
}

void Mesh::deleteGhostPlanes(){
    if(planes.size()==2) return;        // Keep the left and right planes
    planes.erase(planes.begin()+2, planes.end());       // delete the ghost planes
    intersectionTriangles.erase(intersectionTriangles.begin()+2, intersectionTriangles.end());
    //verticesOnPlane.erase(verticesOnPlane.begin()+2, verticesOnPlane.end());
}

void Mesh::updatePlaneIntersections(){
    if(isCut){
        flooding.clear();
        for(unsigned int i=0; i<vertices.size(); i++) flooding.push_back(-1);       // reset the flooding values

        planeNeighbours.clear();
        for(unsigned int i=0; i<planes.size()*2; i++) planeNeighbours.push_back(-1);        // reset the plane neighbours

        for(unsigned int i=0; i<planes.size(); i++) planeIntersection(i);

        for(unsigned int i=0; i<flooding.size(); i++){
            if(flooding[i] != -1){
                for(unsigned int j=0; j<vertexNeighbours[i].size(); j++){
                floodNeighbour(vertexNeighbours[i][j], flooding[i]);
                }
            }
        }
        mergeFlood();
        cutMesh();
    }
}

void Mesh::cutMesh(){
    trianglesCut.clear();

    bool truthTriangles[triangles.size()];  // keeps a record of the triangles who are already added
    for(unsigned int i=0; i<triangles.size(); i++) truthTriangles[i] = false;

    switch (cuttingSide) {
        case Side::INTERIOR:        // MANDIBLE
            cutMandible(truthTriangles);
        break;

        case Side::EXTERIOR:        // FIBULA
            cutFibula(truthTriangles);
        break;
    }

    trianglesExtracted.clear();     // Store the rest of the triangles
    for(unsigned int i=0; i<triangles.size(); i++){
        if(!truthTriangles[i]) trianglesExtracted.push_back(i);
    }

    // ! Conserve this order
    createSmoothedTriangles();

    if(cuttingSide == Side::EXTERIOR){      // fill the colours on the fibula and send the segments to the mandible
        fillColours();
        if(isTransfer){
            sendToManible();
        }
    }

}

void Mesh::cutMandible(bool* truthTriangles){
    for(unsigned int i=0; i<flooding.size(); i++){
        if(planeNeighbours[static_cast<unsigned int>(flooding[i])]==-1){
            saveTrianglesToKeep(truthTriangles, i);
        }
    }
}

void Mesh::cutFibula(bool* truthTriangles){
    for(unsigned int j=0; j<intersectionTriangles.size(); j++){
            const std::vector<unsigned int> &v = intersectionTriangles[j];
            for(unsigned int k=0; k<v.size(); k++) {
                truthTriangles[v[k]] = true;
        }
    }

    getSegmentsToKeep();    // figure out what to keep (TODO can be done earlier)
    for(unsigned int i=0; i<flooding.size(); i++){
        bool isKeep = false;
        for(unsigned int k=0; k<segmentsConserved.size(); k++){      // Only keep it if it belongs to a kept segment
            if(segmentsConserved[k]==flooding[i]){
                isKeep = true;
                break;
            }
        }
        if(isKeep){
            for(unsigned int j=0; j<vertexTriangles[i].size(); j++){        // Get the triangles they belong to
                saveTrianglesToKeep(truthTriangles, i);
            }
        }
    }
}

void Mesh::saveTrianglesToKeep(bool* truthTriangles, unsigned int i){
    for(unsigned int j=0; j<vertexTriangles[i].size(); j++){        // Get the triangles the indicies belong to
        if(!truthTriangles[vertexTriangles[i][j]]){     // If it's not already in the list
            trianglesCut.push_back(vertexTriangles[i][j]);
            truthTriangles[vertexTriangles[i][j]] = true;
        }
    }
}

void Mesh::fillColours(){
    int tempColours[planeNeighbours.size()];
    coloursIndicies.clear();

    for(unsigned int i=0; i<planeNeighbours.size(); i++) tempColours[i] = -1;       // init all to -1

    for(unsigned int i=0; i<segmentsConserved.size(); i++) tempColours[segmentsConserved[i]] = static_cast<int>(i);       // change to the seg colours value

    for(unsigned int i=0; i<vertices.size(); i++){
        int index = flooding[i]; 
        coloursIndicies.push_back(tempColours[index]);      // Fill the colours
    }
}

// WARNING : this assumes that the left and right planes are the first planes added!
// Could search for the exterior planes beforehand using the fact that the other sides = -1
void Mesh::getSegmentsToKeep(){
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

void Mesh::createSmoothedTriangles(){
    smoothedVerticies.clear();


    for(unsigned int i=0; i<vertices.size(); i++) smoothedVerticies.push_back(vertices[i]);  // Copy the verticies table

    switch (cuttingSide) {
        case Side::INTERIOR:
            createSmoothedMandible();
        break;

        case Side::EXTERIOR:
            createSmoothedFibula();
        break;
    }
}

void Mesh::createSmoothedMandible(){
    for(unsigned long long i=0; i<planes.size(); i++){
        for(unsigned long long j=0; j<intersectionTriangles[static_cast<unsigned long long>(i)].size(); j++){       // for each triangle cut
            for(unsigned int k=0; k<3; k++){    // find which verticies to keep
                unsigned int vertexIndex = triangles[intersectionTriangles[i][j]].getVertex(k);
                if(planeNeighbours[static_cast<unsigned int>(flooding[vertexIndex])] != -1){   // if we need to change it (here we only change it if it's outside of the cut (fine for mandible))
                    Vec newVertex = planes[i]->getProjection(Vec(static_cast<double>(vertices[vertexIndex][0]), static_cast<double>(vertices[vertexIndex][1]), static_cast<double>(vertices[vertexIndex][2])) );
                    smoothedVerticies[vertexIndex] = Vec3Df(static_cast<float>(newVertex.x), static_cast<float>(newVertex.y), static_cast<float>(newVertex.z)); // get the projection
                }
                // else don't change the original
            }

        }
    }
}

void Mesh::createSmoothedFibula(){
    for(unsigned int i=0; i<planes.size(); i++){
        //verticesOnPlane[i].clear();
        for(unsigned long long j=0; j<intersectionTriangles[static_cast<unsigned long long>(i)].size(); j++){   // for each triangle cut
            int actualFlooding = -1;    //  Conserve the "real" flooding value (will never stay at -1)

            for(unsigned int k=0; k<3; k++){    // find which verticies are on the otherside of the cut
                unsigned int vertexIndex = triangles[intersectionTriangles[i][j]].getVertex(k);

                bool isOutlier = false;
                for(unsigned int l=0; l<segmentsConserved.size(); l++){
                    if(flooding[vertexIndex] == segmentsConserved[l]){
                        actualFlooding = flooding[vertexIndex];
                        isOutlier = true;
                    }
                }

                if(planeNeighbours[static_cast<unsigned int>(flooding[vertexIndex])]==-1 || isOutlier){        // if we need to change it
                    Vec newVertex;
                    unsigned int lastIndex = planes.size()-1;
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
                    //else newVertex = planes[i]->getProjection(Vec(static_cast<double>(vertices[vertexIndex][0]), static_cast<double>(vertices[vertexIndex][1]), static_cast<double>(vertices[vertexIndex][2])) );
                    //verticesOnPlane[i].push_back(vertexIndex);
                    smoothedVerticies[vertexIndex] = Vec3Df(static_cast<float>(newVertex.x), static_cast<float>(newVertex.y), static_cast<float>(newVertex.z)); // get the projection
                   // Vec newVv = planes[i]->getLocalCoordinates(Vec(smoothedVerticies[vertexIndex]));
                    //if(newVv.z < -0.001) std::cout << " z : " << newVv.z << std::endl;
                }
                // else don't change the original
            }

            // Set the whole triangle to the correct flooding value
            for(unsigned int k=0; k<3; k++){
                unsigned int vertexIndex = triangles[intersectionTriangles[i][j]].getVertex(k);
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
    //if(newVertex.z < -0.001) std::cout << " z : " << newVertex.z << std::endl;
    return planes[p1]->getMeshCoordinatesFromLocal(newVertex);
}

void Mesh::updatePlaneIntersections(Plane *p){
    // Possible optimisation?

    updatePlaneIntersections();
}

void Mesh::floodNeighbour(unsigned int index, int id){

    if(flooding[index] == -1){      // Flood it
        flooding[index] = id;
        for(unsigned int i=0; i<vertexNeighbours[index].size(); i++){
            floodNeighbour(vertexNeighbours[index][i], id);
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
}

void Mesh::mergeFlood(){
    for(unsigned int i=0; i<flooding.size(); i++){
        if(planeNeighbours[static_cast<unsigned int>(flooding[i])] != -1 && planeNeighbours[static_cast<unsigned int>(flooding[i])] < flooding[i]){     // From the two neighbours, set them both to the lowest value
            flooding[i] = planeNeighbours[static_cast<unsigned int>(flooding[i])];
        }
    }
}

// Finds all the intersecting triangles for plane nb index
void Mesh::planeIntersection(unsigned int index){
    intersectionTriangles[index].clear();       // empty the list of intersections

    for(unsigned int i = 0 ; i < triangles.size(); i++){
        unsigned int t0 = triangles[i].getVertex(0);
        unsigned int t1 = triangles[i].getVertex(1);
        unsigned int t2 = triangles[i].getVertex(2);

        if(planes[index]->isIntersection(Vec(vertices[t0]), Vec(vertices[t1]), Vec(vertices[t2]) )){        // if the triangle intersects the plane
            intersectionTriangles[index].push_back(i);      // save the triangle index

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

void Mesh::sendToManible(){
    std::vector<int> planeNb;       // the plane nb associated
    std::vector<Vec> convertedVerticies;    // the vertex coordinates in relation to the plane nb
    std::vector<std::vector<int>>convertedTriangles; // the new indicies of the triangles (3 indicies)
    std::vector<int> convertedColours;
    std::vector<Vec> convertedNormals;
    int tempVerticies[smoothedVerticies.size()];   // a temporary marker for already converted verticies

    for(unsigned int i=0; i<smoothedVerticies.size(); i++) tempVerticies[i] = -1;


    for(unsigned int i=0; i<trianglesCut.size(); i++){      // For every triangle we want to send (we've already filtered out the rest when cutting the mesh)
        unsigned int triVert;       // Get the 3 verticies of the triangle
        std::vector<int> newTriangle;

        for(unsigned int j=0; j<3; j++){
            triVert = triangles[trianglesCut[i]].getVertex(j);       // this must go into smoothedVerticies[triV....] (ie. index in the original)

            if(tempVerticies[triVert] != -1) newTriangle.push_back(tempVerticies[triVert]);     // If converted already

            else{       // convert to the corresponding plane
                int pNb = flooding[triVert];    // Get the plane nb
                if(pNb >= static_cast<int>(planes.size())) pNb -= planes.size();      // The plane nb is referenced by the smallest side
                planeNb.push_back(pNb);

                // Convert the vertex
                Vec unConverted = Vec(static_cast<double>(smoothedVerticies[triVert][0]), static_cast<double>(smoothedVerticies[triVert][1]), static_cast<double>(smoothedVerticies[triVert][2]));
                Vec convertedCoords = planes[static_cast<unsigned int>(pNb)]->getLocalCoordinates(unConverted);

                convertedVerticies.push_back(convertedCoords);      // add the converted vertex triVert
                convertedColours.push_back(coloursIndicies[triVert]); // add its corresponding colour and normal

                Vec unConvertedNorm = Vec(static_cast<double>(verticesNormals[triVert][0]), static_cast<double>(verticesNormals[triVert][1]), static_cast<double>(verticesNormals[triVert][2]));
                Vec convertedNorm = planes[static_cast<unsigned int>(pNb)]->getLocalVector(unConvertedNorm);
                convertedNormals.push_back(convertedNorm);

                int vertexIndex = static_cast<int>(convertedVerticies.size()) - 1;
                newTriangle.push_back(vertexIndex);  // set it to the last index of convertedVerticies

                tempVerticies[triVert] = vertexIndex;       // Store the corresponding index in tempVerticies
            }
        }      
        convertedTriangles.push_back(newTriangle);      // Add the triangle
    }

    Q_EMIT sendInfoToManible(planeNb, convertedVerticies, convertedTriangles, convertedColours, convertedNormals, (static_cast<int>(planes.size())/2));
}

void Mesh::recieveInfoFromFibula(std::vector<Vec> convertedVerticies, std::vector<std::vector<int>> convertedTriangles, std::vector<int> convertedColours, std::vector<Vec> convertedNormals, int nbColours){
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
        Triangle t = Triangle(static_cast<unsigned int>(convertedTriangles[i][0]), static_cast<unsigned int>(convertedTriangles[i][1]), static_cast<unsigned int>(convertedTriangles[i][2]));
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
        for(unsigned int i = 0 ; i < triangles.size(); i++){
            glTriangle(i);
        }
    }
    else{
        for(unsigned int i = 0 ; i < trianglesCut.size(); i++){
            glTriangleSmooth(trianglesCut[i]);
        }

        for(unsigned int i=0; i<fibInMandTriangles.size(); i++){
            glTriangleFibInMand(i);
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
    for(unsigned int i = 0 ; i < trianglesExtracted.size(); i++){
        glTriangleSmooth(trianglesExtracted[i]);
    }

    glEnd();

    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_DEPTH);
}

float Mesh::getBBRadius(){
    computeBB();
    return radius;
}

void Mesh::readJSON(const QJsonObject &json, double &scale){
    if(json.contains("vertices") && json["vertices"].isArray()){
        vertices.clear();
        QJsonArray vArray = json["vertices"].toArray();
        for(int i=0; i<vArray.size(); i++){
            QJsonArray singleV = vArray[i].toArray();
            vertices.push_back(Vec3Df(singleV[0].toDouble(), singleV[1].toDouble(), singleV[2].toDouble()));
        }
    }

    if(json.contains("triangles") && json["triangles"].isArray()){
        triangles.clear();
        vertexNeighbours.clear();
        vertexTriangles.clear();

        for(unsigned int i=0; i<vertices.size(); i++){
            std::vector<unsigned int> init;
            vertexNeighbours.push_back(init);
            vertexTriangles.push_back(init);
        }

        QJsonArray tArray = json["triangles"].toArray();
        for(int i=0; i<tArray.size(); i++){
            QJsonArray singleT = tArray[i].toArray();
            unsigned int triIndex = static_cast<unsigned int>(triangles.size());
            unsigned int vert[3] = {static_cast<unsigned int>(singleT[0].toInt()), static_cast<unsigned int>(singleT[1].toInt()), static_cast<unsigned int>(singleT[2].toInt())};
            triangles.push_back(Triangle(vert[0], vert[1], vert[2]));

            // Add to neighbours
            for(int k=0; k<3; k++){
                bool found = false;
                for(unsigned int i=0; i<vertexNeighbours[vert[k]].size(); i++){
                    if(vertexNeighbours[vert[k]][i] == vert[(k+1)%3]){
                        found = true;
                        break;
                    }
                }
                if(found == false){
                    vertexNeighbours[vert[k]].push_back(vert[(k+1)%3]);
                    vertexNeighbours[vert[(k+1)%3]].push_back(vert[k]);
                }
            }

            // Add to vertexTriangles
            for(int k=0; k<3; k++){
                bool found = false;
                for(unsigned int i=0; i<vertexTriangles[vert[k]].size(); i++){
                    if(vertexTriangles[vert[k]][i] == triIndex){
                        found = true;
                        break;
                    }
                }
                if(found == false){
                    vertexTriangles[vert[k]].push_back(triIndex);
                }
            }
        }
    }

    if(json.contains("scale") && json["scale"].isDouble()){
        scale = 1.0/json["scale"].toDouble();
        uniformScale(scale);
    }

    update();
}

void Mesh::uniformScale(float s){
    for(unsigned int i=0; i<vertices.size(); i++) vertices[i] *= s;
    BBMax *= s;
    BBMin *= s;
    BBCentre *= s;
    radius *= s;
}

std::vector<unsigned int> Mesh::getVerticesOnPlane(unsigned int planeNb){
    std::vector<unsigned int> v;
    for(unsigned int i=0; i<intersectionTriangles[planeNb].size(); i++){
        unsigned int index = intersectionTriangles[planeNb][i];
        if(abs(planes[planeNb]->getLocalCoordinates(Vec(smoothedVerticies[index])).z) < 0.001){
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

    return v;
}

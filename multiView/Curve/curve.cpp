#include "curve.h"
#include "math.h"
#include <GL/gl.h>

Curve::Curve(){
    this->nbU = 0;
}

// Initialise the curve: set the control points / nb of control points
void Curve::init(unsigned int nbCP, std::vector<Vec>& cntrlPoints){
    nbControlPoint = nbCP;

    for(unsigned int i=0; i<nbCP; i++){
        TabControlPoint.push_back(new ControlPoint(cntrlPoints[i]));
    }

    initConnections();
}

// If a control point is moved, reinitialise the curve
void Curve::initConnections(){
    for(unsigned int i=0; i<nbControlPoint; i++){
        connect(TabControlPoint[i], &ControlPoint::cntrlPointTranslated, this, &Curve::reintialiseCurve);
    }
}

// Catmul-Rom generation process (this is the only public creation method)
void Curve::generateCatmull(unsigned int& n){
    unsigned int nbSeg = nbControlPoint-3;

    this->nbU = n - n%nbSeg;        // make sure the nbU can be divided up evenly between each segment
    n = nbU;    // update the nbU in viewer
    this->knotIndex = 0;
    this->degree = 3;

    generateCatmullKnotVector(0.3, this->knotVector);       // Generate the knot vector for the catmull rom
    catmullrom();       // Create the catmull rom
}

// Re-calculate the catmull rom and send a notification signal
void Curve::reintialiseCurve(){
    catmullrom();
    Q_EMIT curveReinitialised();
}

void Curve::generateCatmullKnotVector(double alpha, std::vector<double>& kv){
    kv.resize(static_cast<unsigned long long>(nbControlPoint));

    kv[0] = 0;

    for(unsigned int i=1; i<nbControlPoint; i++){
        Vec p = TabControlPoint[i]->getPoint() - TabControlPoint[i-1]->getPoint();
        kv[i] =  pow(p.norm(),alpha) + kv[i-1];
    }
}

/* Catmull rom */

// Get a point in the catmull rom
void Curve::calculateCatmullPoints(Vec& c, Vec& cp, Vec& cpp, double t){
    Vec p[4] = {TabControlPoint[knotIndex-1]->getPoint(), TabControlPoint[knotIndex]->getPoint(), TabControlPoint[knotIndex+1]->getPoint(), TabControlPoint[knotIndex+2]->getPoint()};

    const double &t0 = knotVector[knotIndex-1];
    const double &t1 = knotVector[knotIndex];
    const double &t2 = knotVector[knotIndex+1];
    const double &t3 = knotVector[knotIndex+2];

    Vec a1 = (t1-t)/(t1-t0)*p[0] + (t-t0)/(t1-t0)*p[1];
    Vec a2 = (t2-t)/(t2-t1)*p[1] + (t-t1)/(t2-t1)*p[2];
    Vec a3 = (t3-t)/(t3-t2)*p[2] + (t-t2)/(t3-t2)*p[3];

    Vec a1p = 1.0/(t1-t0)*(p[1]-p[0]);
    Vec a2p = 1.0/(t2-t1)*(p[2]-p[1]);
    Vec a3p = 1.0/(t3-t2)*(p[3]-p[2]);

    Vec b1 = (t2-t)/(t2-t0)*a1 + (t-t0)/(t2-t0)*a2;
    Vec b2 = (t3-t)/(t3-t1)*a2 + (t-t1)/(t3-t1)*a3;

    Vec b1p = 1.0/(t2-t0)*(a2-a1) + (t2-t)/(t2-t0)*a1p + (t-t0)/(t2-t0)*a2p;
    Vec b2p = 1.0/(t3-t1)*(a3-a2) + (t3-t)/(t3-t1)*a2p + (t-t1)/(t3-t1)*a3p;

    Vec b1pp = 1.0/(t2-t0)*(a2p-a1p);
    Vec b2pp = 1.0/(t3-t1)*(a3p-a2p);

    c = (t2-t)/(t2-t1)*b1 + (t-t1)/(t2-t1)*b2;      // Catmull rom
    cp = 1.0/(t2-t1)*(b2-b1) + (t2-t)/(t2-t1)*b1p + (t-t1)/(t2-t1)*b2p;     // First derivative
    cpp = 1.0/(t2-t1)*(b2p-b1p) + (t2-t)/(t2-t1)*b1pp + (t-t1)/(t2-t1)*b2pp;        // Second derivative
}

// Get the catmull rom
void Curve::catmullrom(){
    unsigned int nbSeg = nbControlPoint-3;
    unsigned int uPerSeg = nbU/nbSeg;       // the nbU per segment

    // Reset the curve vectors (curve, first and second derivative)
    curve.clear();
    curve.resize(nbU);
    dt.clear();
    dt.resize(nbU);
    d2t.clear();
    d2t.resize(nbU);

    // For each knot
    for(unsigned int j=1; j<=nbSeg; j++){
        unsigned int it=0;
        knotIndex = j;

        // Generate the segment
        for(double i=knotVector[j]; i<knotVector[j+1]; i+=((knotVector[j+1]-knotVector[j])/static_cast<double>(uPerSeg))){
            if((j-1)*uPerSeg+it >= nbU) return;
            curve[(j-1)*uPerSeg+it] = Vec();
            dt[(j-1)*uPerSeg+it] = Vec();
            d2t[(j-1)*uPerSeg+it] = Vec();

            calculateCatmullPoints(curve[(j-1)*uPerSeg+it], dt[(j-1)*uPerSeg+it], d2t[(j-1)*uPerSeg+it], i);
            it++;
        }
    }
}

// Euclidean distance between two points
double Curve::discreteLength(unsigned int indexS, unsigned int indexE){
    return sqrt( pow((curve[indexE].x - curve[indexS].x), 2.0) + pow((curve[indexE].y - curve[indexS].y), 2.0) + pow((curve[indexE].z - curve[indexS].z), 2.0));
}

// Length of the chord
double Curve::discreteChordLength(unsigned int indexS, unsigned int indexE){
    double sum = 0;

    for(unsigned int i=indexS; i<indexE; i++){
        sum += discreteLength(i, i+1);
    }

    return sum;
}

// Check which index is closest to the target distance from indexS
unsigned int Curve::getClosestDistance(double target, unsigned int indexS, unsigned int a, unsigned int b){
    double aDist = discreteLength(indexS, indexS+a);
    double bDist = discreteLength(indexS, indexS+b);
    double aTarget = abs(target - aDist);
    double bTarget = abs(target - bDist);

    if(aTarget < bTarget) return a;
    return b;
}

// Returns the index which is length away from indexS   (TODO optimise)
unsigned int Curve::indexForLength(unsigned int indexS, double length){
    unsigned int i=0;

    if(length > 0){
        while(indexS+i < nbU-1 && discreteLength(indexS, indexS+i) < length) i++;
        if(i!=0) i = getClosestDistance(length, indexS, i, i-1);
    }
    else{
        while(indexS+i > 0 && discreteLength(indexS, indexS+i) < abs(length)) i--;
        if(i!=nbU-1) i = getClosestDistance(length, indexS, i, i+1);
    }

    return indexS+i;
}

void Curve::draw(){
    /*glEnable(GL_DEPTH);
    glEnable(GL_DEPTH_TEST);*/

      glBegin(GL_LINE_STRIP);
      glColor3f(0.0, 1.0, 0.0);

      for(unsigned int i=0; i<nbU; i++){
        const Vec &p = curve[i];
        glVertex3f(static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z));
      }

      glEnd();

      /*glDisable(GL_DEPTH);
      glDisable(GL_DEPTH_TEST);*/
}

// Draw the control points
void Curve::drawControl(){
      glBegin(GL_LINE_STRIP);
      glColor3f(0.0, 0.0, 1.0);

      for(unsigned int i=0; i<nbControlPoint; i++){
        Vec *p = new Vec(TabControlPoint[i]->getX(), TabControlPoint[i]->getY(), TabControlPoint[i]->getZ());
        glVertex3f(static_cast<float>(p->x), static_cast<float>(p->y), static_cast<float>(p->z));
      }

      glEnd();

      for(unsigned int i=0; i<nbControlPoint; i++){
          TabControlPoint[i]->draw();
      }
}

// Frenet frame
Vec Curve::tangent(unsigned int index){
    Vec t(dt[index].x, dt[index].y, dt[index].z);
    t.normalize();

    return t;
}

Vec Curve::normal(unsigned int index){
    return cross(binormal(index), tangent(index));
}

Vec Curve::binormal(unsigned int index){
    Vec b = cross(Vec(dt[index].x, dt[index].y, dt[index].z), Vec(d2t[index].x, d2t[index].y, d2t[index].z));
    b.normalize();

    return b;
}

void Curve::getFrame(unsigned int index, Vec &t, Vec &n, Vec &b){
    t = tangent(index);
    b = binormal(index);
    n = cross(b, t);
}

void Curve::drawFrame(unsigned int index){
    Vec t,n,b;

    getFrame(index,t,n,b);

    glColor3f(0.0, 1.0, 1.0);
    glLineWidth(3);

    glBegin(GL_LINES);
      glVertex3f(static_cast<float>(curve[index].x), static_cast<float>(curve[index].y), static_cast<float>(curve[index].z));
      glVertex3f(static_cast<float>(curve[index].x + t.x*10), static_cast<float>(curve[index].y + t.y*10), static_cast<float>(curve[index].z + t.z*10));
    glEnd();

    glColor3f(1.0, 0.0, 1.0);
    glBegin(GL_LINES);
      glVertex3f(static_cast<float>(curve[index].x), static_cast<float>(curve[index].y), static_cast<float>(curve[index].z));
      glVertex3f(static_cast<float>(curve[index].x + n.x*10), static_cast<float>(curve[index].y + n.y*10), static_cast<float>(curve[index].z + n.z*10));
    glEnd();

     glColor3f(1.0, 1.0, 0.0);
     glBegin(GL_LINES);
      glVertex3f(static_cast<float>(curve[index].x), static_cast<float>(curve[index].y), static_cast<float>(curve[index].z));
      glVertex3f(static_cast<float>(curve[index].x + b.x*10), static_cast<float>(curve[index].y + b.y*10), static_cast<float>(curve[index].z + b.z*10));
    glEnd();

    glLineWidth(1);
}

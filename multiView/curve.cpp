#include "curve.h"
#include "math.h"
#include <GL/gl.h>

Curve::Curve(unsigned int nbCP, std::vector<Vec>& cntrlPoints){
    this->nbU = 0;
    nbControlPoint = nbCP;

    for(unsigned int i=0; i<nbCP; i++){
        TabControlPoint.push_back(new ControlPoint(cntrlPoints[i]));
    }

    initConnections();
}

void Curve::initConnections(){
    for(unsigned int i=0; i<nbControlPoint; i++){
        connect(TabControlPoint[i], &ControlPoint::cntrlPointTranslated, this, &Curve::reintialiseCurve);
    }
}

void Curve::generateBSpline(unsigned int& nbU, unsigned int degree){
    this->nbU = nbU;
    this->degree = degree;
    this->knotIndex = 0;

    generateUniformKnotVector(0, this->knotVector);
    splineDerivative(0, curve);
    splineDerivative(1, dt);
    splineDerivative(2, d2t);
}

void Curve::generateCatmull(unsigned int& n){
    unsigned int nbSeg = nbControlPoint-3;

    this->nbU = n - n%nbSeg;
    n = nbU;    // update the nbU in viewer
    this->knotIndex = 0;
    this->degree = 3;

    generateCatmullKnotVector(0.3, this->knotVector);
    catmullrom();
}


Vec Curve::deBoor(double u, unsigned int j, unsigned int r){

    if(r==0) return Vec(TabControlPoint[j]->getX(), TabControlPoint[j]->getY(),TabControlPoint[j]->getZ());

    double alpha = (u - knotVector[j]) / (knotVector[j+degree-(r-1)] - knotVector[j]);

    return (1.0 - alpha) * deBoor(u, j-1, r-1) + alpha * deBoor(u, j, r-1);
}


// Returns the kth derivative of the curve ( 0 <= k <= 2 )
void Curve::splineDerivative(unsigned int k, std::vector<Vec> &c){
    c.clear();
    c.resize(nbU);
    this->knotIndex = 0;

    for(unsigned int i=0; i<nbU; i++){
        double u = (1.0 / static_cast<double>(nbU-1)) * static_cast<double>(i);
        c[i] = Vec();

        while(u >= knotVector[knotIndex+1] && knotVector[knotIndex+1] != 1.0) knotIndex++;

        c[i] += Vec(deBoorDerivative(u, knotIndex, degree, k));
    }
}

Vec Curve::deBoorDerivative(double u, unsigned int j, unsigned int r, unsigned int k){

    if(r > degree - k){
        double denom = static_cast<double>(degree + nbControlPoint + 1) - 2.0* static_cast<double>(degree) - 1.0;
        double rnorm = static_cast<double>(r) / denom;
        if(k==2){
                double beta = rnorm / (knotVector[j+degree-(r-1)] - knotVector[j]);
                return - beta * deBoorDerivative(u, j-1, r-1, k) + beta * deBoorDerivative(u, j, r-1, k);
            }
        else if(k==1){
            double beta = rnorm / knotVector[j+1] - knotVector[j];
            return beta * (deBoor(u, j, r-1) - deBoor(u, j-1, r-1));
        }
    }
    return deBoor(u, j, r);
}

void Curve::generateUniformKnotVector(unsigned int a, std::vector<double>& kv){
    unsigned int k = degree - a;
    unsigned int n = nbControlPoint - a;
    unsigned int m = n + k + 1;

    kv.resize(static_cast<unsigned long long>(m));

    double denom = static_cast<double>(m) - 2.0* static_cast<double>(k) - 1.0;

    for(unsigned int i=0; i<=k; i++) kv[i] = 0;
    for(unsigned int i=k+1; i<m-k-1; i++) kv[i] = static_cast<double>(i-k) / denom;
    for(unsigned int i=m-k-1; i<m; i++) kv[i] = 1.0;
}

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

// Catmull rom
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

    c = (t2-t)/(t2-t1)*b1 + (t-t1)/(t2-t1)*b2;
    cp = 1.0/(t2-t1)*(b2-b1) + (t2-t)/(t2-t1)*b1p + (t-t1)/(t2-t1)*b2p;
    cpp = 1.0/(t2-t1)*(b2p-b1p) + (t2-t)/(t2-t1)*b1pp + (t-t1)/(t2-t1)*b2pp;
}

void Curve::catmullrom(){
    unsigned int nbSeg = nbControlPoint-3;
    unsigned int uPerSeg = nbU/nbSeg;

    curve.clear();
    curve.resize(nbU);
    dt.clear();
    dt.resize(nbU);
    d2t.clear();
    d2t.resize(nbU);

    for(unsigned int j=1; j<=nbSeg; j++){
        unsigned int it=0;
        knotIndex = j;

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

// Length as the crow flies
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

unsigned int Curve::getClosestDistance(double target, unsigned int indexS, unsigned int a, unsigned int b){
    double aDist = discreteLength(indexS, indexS+a);
    double bDist = discreteLength(indexS, indexS+b);
    double aTarget = abs(target - aDist);
    double bTarget = abs(target - bDist);

    if(aTarget < bTarget) return a;
    return b;
}

// Returns the index which is length away from indexS
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
    glEnable(GL_DEPTH);
    glEnable(GL_DEPTH_TEST);

      glBegin(GL_LINE_STRIP);
      glColor3f(0.0, 1.0, 0.0);

      for(unsigned int i=0; i<nbU; i++){
        const Vec &p = curve[i];
        glVertex3f(static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z));
      }

      glEnd();

      glDisable(GL_DEPTH);
      glDisable(GL_DEPTH_TEST);
}

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

void Curve::drawTangent(unsigned int index){
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

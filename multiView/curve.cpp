#include "curve.h"
#include "math.h"
#include <GL/gl.h>

/*Curve::Curve(long nbCP)
{
    this->nbU = new long();
    nbControlPoint = nbCP;
    TabControlPoint = new ControlPoint*[static_cast<unsigned long long>(nbControlPoint)];

    for(int i=0; i<nbCP; i++){
        TabControlPoint[i] = new ControlPoint();
    }

    initConnections();
}*/

Curve::Curve(long nbCP, std::vector<Vec>& cntrlPoints){
    this->nbU = 0;
    nbControlPoint = nbCP;

    //this->TabControlPoint = cntrlPoints;

    for(unsigned int i=0; i<nbCP; i++){
        TabControlPoint.push_back(new ControlPoint(static_cast<Vec>(cntrlPoints[i])));
    }

    initConnections();
}

void Curve::initConnections(){
    for(int i=0; i<nbControlPoint; i++){
        connect(TabControlPoint[i], &ControlPoint::cntrlPointTranslated, this, &Curve::reintialiseCurve);
    }
}

/*void Curve::updateConnections(ControlPoint* p){
    connect(p, &ControlPoint::cntrlPointTranslated, this, &Curve::reintialiseCurve);
}*/

void Curve::generateBSpline(long& nbU, int degree){
    this->nbU = nbU;
    this->degree = degree;
    this->knotIndex = 0;

    this->knotVector = generateUniformKnotVector(0);
    curve = splineDerivative(0);
    dt =  splineDerivative(1);
}

void Curve::generateCatmull(long& n){
    int nbSeg = nbControlPoint-3;

    this->nbU = n;
    this->nbU -= n%nbSeg;
    n = nbU;    // update the nbU in viewer
    this->knotIndex = 0;
    this->degree = 3;

    this->knotVector = generateCatmullKnotVector(0.3);
    catmullrom();
}


Vec Curve::deBoor(double u, int j, int r){

    if(r==0) return Vec(TabControlPoint[j]->getX(), TabControlPoint[j]->getY(),TabControlPoint[j]->getZ());

    double alpha = (u - knotVector[j]) / (knotVector[j+degree-(r-1)] - knotVector[j]);

    return (1.0 - alpha) * deBoor(u, j-1, r-1) + alpha * deBoor(u, j, r-1);
}


// Returns the kth derivative of the curve ( 0 <= k <= 2 )
Vec* Curve::splineDerivative(int k){

    this->knotIndex = 0;

    Vec *c = new Vec[static_cast<unsigned long long>(nbU)];

    for(int i=0; i<nbU; i++){
        double u = (1.0 / static_cast<double>(nbU-1)) * static_cast<double>(i);
        c[i] = Vec();

        while(u >= knotVector[knotIndex+1] && knotVector[knotIndex+1] != 1.0) knotIndex++;

        c[i] += Vec(deBoorDerivative(u, knotIndex, degree, k));
    }
    return c;
}

Vec Curve::deBoorDerivative(double u, int j, int r, int k){

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

std::vector<double> Curve::generateUniformKnotVector(int a){
    int k = degree - a;
    int n = nbControlPoint - a;
    int m = n + k + 1;
    std::vector<double> kv;
    kv.resize(static_cast<unsigned long long>(m));

    double denom = static_cast<double>(m) - 2.0* static_cast<double>(k) - 1.0;

    for(int i=0; i<=k; i++) kv[i] = 0;
    for(int i=k+1; i<m-k-1; i++) kv[i] = static_cast<double>(i-k) / denom;
    for(int i=m-k-1; i<m; i++) kv[i] = 1.0;

    return kv;
}

void Curve::reintialiseCurve(){
   catmullrom();

    /*curve = splineDerivative(0);
    dt =  splineDerivative(1);*/

    Q_EMIT curveReinitialised();
}

std::vector<double> Curve::generateCatmullKnotVector(double alpha){
    std::vector<double> kv;
    kv.resize(static_cast<unsigned long long>(nbControlPoint));

    kv[0] = 0;

    for(int i=1; i<nbControlPoint; i++){
        Vec p = TabControlPoint[i]->getPoint() - TabControlPoint[i-1]->getPoint();
        kv[i] =  pow(p.norm(),alpha) + kv[i-1];
    }

    return kv;
}

// Catmull rom
void Curve::calculateCatmullPoints(Vec& c, Vec& cp, double t){
    Vec p[4] = {TabControlPoint[knotIndex-1]->getPoint(), TabControlPoint[knotIndex]->getPoint(), TabControlPoint[knotIndex+1]->getPoint(), TabControlPoint[knotIndex+2]->getPoint()};

    double t0 = knotVector[knotIndex-1];
    double t1 = knotVector[knotIndex];
    double t2 = knotVector[knotIndex+1];
    double t3 = knotVector[knotIndex+2];

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

    c = (t2-t)/(t2-t1)*b1 + (t-t1)/(t2-t1)*b2;
    cp = 1.0/(t2-t1)*(b2-b1) + (t2-t)/(t2-t1)*b1p + (t-t1)/(t2-t1)*b2p;
}

void Curve::catmullrom(){
    int nbSeg = nbControlPoint-3;
    int uPerSeg = nbU/nbSeg;

    curve = new Vec[static_cast<unsigned long long>(nbU)];
    dt = new Vec[static_cast<unsigned long long>(nbU)];
   // controlPointIndicies = new bool[static_cast<unsigned long long>(*nbU)];

    //for(int i=0; i<*nbU; i++) controlPointIndicies[i] = false;

    for(int j=1; j<=nbSeg; j++){
        int it=0;
        knotIndex = j;
        //controlPointIndicies[(j-1)*uPerSeg] = true;
        for(double i=knotVector[j]; i<knotVector[j+1]; i+=((knotVector[j+1]-knotVector[j])/static_cast<double>(uPerSeg))){
            if((j-1)*uPerSeg+it >= nbU) return;
            curve[(j-1)*uPerSeg+it] = Vec();
            dt[(j-1)*uPerSeg+it] = Vec();

            calculateCatmullPoints(curve[(j-1)*uPerSeg+it], dt[(j-1)*uPerSeg+it], i);
            it++;
        }
    }
}

// Length as the crow flies
double Curve::discreteLength(int indexS, int indexE){
    return sqrt( pow((curve[indexE].x - curve[indexS].x), 2) + pow((curve[indexE].y - curve[indexS].y), 2) + pow((curve[indexE].z - curve[indexS].z), 2));
}

// Length of the chord
double Curve::discreteChordLength(int indexS, int indexE){
    double sum = 0;

    for(int i=indexS; i<indexE; i++){
        sum += discreteLength(i, i+1);
    }

    return sum;
}

// Returns the index which is length away from indexS
int Curve::indexForLength(int indexS, double length){
    int i=0;

    // First condition should never happen, just extra protection to prevent crashes
    if(length > 0){
        while(indexS+i < nbU-1 && discreteLength(indexS, indexS+i) < length) i++;
    }
    else{
        while(indexS+i > 0 && discreteLength(indexS, indexS+i) < abs(length)) i--;
    }

    return indexS+i;
}

void Curve::draw(){
    glEnable(GL_DEPTH);
    glEnable(GL_DEPTH_TEST);

      glBegin(GL_LINE_STRIP);
      glColor3f(0.0, 1.0, 0.0);

      for(int i=0; i<nbU; i++){
        Vec p = curve[i];
        glVertex3f(static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z));
      }

      glEnd();

      glDisable(GL_DEPTH);
      glDisable(GL_DEPTH_TEST);
}

void Curve::drawControl(){
      glBegin(GL_LINE_STRIP);
      glColor3f(0.0, 0.0, 1.0);

      for(int i=0; i<nbControlPoint; i++){
        Vec *p = new Vec(TabControlPoint[i]->getX(), TabControlPoint[i]->getY(), TabControlPoint[i]->getZ());
        glVertex3f(static_cast<float>(p->x), static_cast<float>(p->y), static_cast<float>(p->z));
      }

      glEnd();

      for(int i=0; i<nbControlPoint; i++){
          TabControlPoint[i]->draw();
      }
}

// Frenet frame
/*bool Curve::isControlPoint(int index){
    if(controlPointIndicies[index]) return true;
    return false;
}*/

Vec Curve::tangent(int index){
    Vec t = Vec(dt[index].x, dt[index].y, dt[index].z);
    t.normalize();

    return t;
}

void Curve::drawTangent(int index){
    Vec t = tangent(index);

    glColor3f(0.0, 1.0, 1.0);
    glLineWidth(3);

    glBegin(GL_LINES);
      glVertex3f(static_cast<float>(curve[index].x), static_cast<float>(curve[index].y), static_cast<float>(curve[index].z));
      glVertex3f(static_cast<float>(curve[index].x + t.x*10), static_cast<float>(curve[index].y + t.y*10), static_cast<float>(curve[index].z + t.z*10));
    glEnd();

    glLineWidth(1);
}

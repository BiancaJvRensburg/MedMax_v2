#ifndef SAVEDSTATE_H__
#define SAVEDSTATE_H__

#pragma once

#include "Planes/plane.h"
#include "Polyline/box.h"

enum Modification {PLANE, BOX};

class SavedState {
public:
  inline SavedState(Modification m){
      mod = m;
  }

  inline virtual ~SavedState () {}

  inline void addBoxes(std::vector<Box> &b){
      for(unsigned int i=0; i<b.size(); i++){
        boxPositions.push_back(b[i].getLocation());
        boxLengths.push_back(b[i].getLength());
        Vec x,y,z;
        b[i].getOrientation(x,y,z);
        boxOrientationsX.push_back(x);
        boxOrientationsY.push_back(y);
        boxOrientationsZ.push_back(z);
      }
  }

  inline void addPlanes(std::vector<Plane*> p){
      for(unsigned int i=0; i<p.size(); i++){
          planePositions.push_back(p[i]->getPosition());
          planeOrientations.push_back(p[i]->getOrientation());
      }
  }

  inline void getPlanePositions(std::vector<Vec> &v){ v = planePositions; }
  inline void getBoxPositions(std::vector<Vec> &v){ v = boxPositions; }
  inline void getBoxXOrientations(std::vector<Vec> &v){ v = boxOrientationsX; }
  inline void getBoxYOrientations(std::vector<Vec> &v){ v = boxOrientationsY; }
  inline void getBoxZOrientations(std::vector<Vec> &v){ v = boxOrientationsZ; }
  inline void getBoxLengths(std::vector<double> &v){ v = boxLengths; }
  inline void getPlaneOrientations(std::vector<Quaternion> &v){ v = planeOrientations; }
  inline Modification getModification() { return mod; }

  inline const Vec getBoxX(unsigned int id){ return boxOrientationsX[id]; }
  inline const Vec getBoxY(unsigned int id){ return boxOrientationsY[id]; }
  inline const Vec getBoxZ(unsigned int id){ return boxOrientationsZ[id]; }

private:
  std::vector<Vec> planePositions;
  std::vector<Vec> boxPositions;
  std::vector<Quaternion> planeOrientations;
  std::vector<Vec> boxOrientationsX;
  std::vector<Vec> boxOrientationsY;
  std::vector<Vec> boxOrientationsZ;
  std::vector<double> boxLengths;
  Modification mod;
};

#endif

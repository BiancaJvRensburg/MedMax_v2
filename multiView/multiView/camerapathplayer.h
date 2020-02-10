#ifndef CAMERAPATHPLAYER_H
#define CAMERAPATHPLAYER_H

#include <QGLViewer/qglviewer.h>
#include <QGLViewer/manipulatedFrame.h>

using namespace qglviewer;

class CameraPathPlayer : public qglviewer::MouseGrabber {
public:
  CameraPathPlayer(int nb) : pathNb(nb) {}
  void checkIfGrabsMouse(int x, int y, const qglviewer::Camera *const camera);
  int yPos() { return 25 * pathNb; }

protected:
  void mousePressEvent(QMouseEvent *const, qglviewer::Camera *const camera) {
    camera->playPath(static_cast<unsigned int>(pathNb));
  }

private:
  int pathNb;
};

#endif // CAMERAPATHPLAYER_H

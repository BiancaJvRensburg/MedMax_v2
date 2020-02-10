#include "camerapathplayer.h"

void CameraPathPlayer::checkIfGrabsMouse(int x, int y, const Camera *const) {
  setGrabsMouse((x < 80) && (y < yPos()) && ((yPos() - y) < 16));
}

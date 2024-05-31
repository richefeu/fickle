#pragma once

#include <cmath>
#include <fstream>
#include <vector>

struct DICParticle {

  int refcoord_xpix;
  int refcoord_ypix;
  double refrot;
  double radius_pix;
  double dx;
  double dy;
  double drot;
  double upix;
  double vpix;
  double rot_inc;
  double NCC;
  double NCC_rescue;
  double NCC_subpix;
};

class DICfile {
 public:
  std::vector<DICParticle> particles;
  std::vector<DICParticle> corners;
  std::vector<DICParticle> sticks;
  std::vector<DICParticle> fixed;
  std::vector<int> cornerID{0, 1, 2, 3};

  double meterPerPixel;

  void read(const char* fname, double actualDistance = 0.40) {

    std::ifstream file(fname);

    int nbParticles;
    file >> nbParticles;

    DICParticle P;
    for (int i = 0; i < nbParticles; ++i) {
      file >> P.refcoord_xpix >> P.refcoord_ypix >> P.refrot >> P.radius_pix >> P.dx >> P.dy >> P.drot >> P.upix >>
          P.vpix >> P.rot_inc >> P.NCC >> P.NCC_rescue >> P.NCC_subpix;
      if (P.radius_pix == 1.0) {
        corners.push_back(P);
      } else if (P.radius_pix == 2.1) {
        sticks.push_back(P);
      } else if (P.radius_pix == 2.2) {
        fixed.push_back(P);
      } else {
        particles.push_back(P);
      }
    }

    // Precomputations
    identifyCorners();
    indentifyPixelLength(actualDistance);
  }

  // CHECK DIFFERENCE BETWEEN pointers and references
  DICParticle& getCorner(int n) { return corners[cornerID[n]]; }

  void identifyCorners() {
    if (corners.size() != 4) {
      std::cout << "there should be 4 corners!" << std::endl;
    }
    double xmean = 0.0;
    double ymean = 0.0;
    for (size_t i = 0; i < corners.size(); ++i) {
      xmean += corners[i].refcoord_xpix;
      ymean += corners[i].refcoord_ypix;
    }
    xmean *= 0.25;
    ymean *= 0.25;

    for (size_t i = 0; i < corners.size(); ++i) {
      double dx = xmean - corners[i].refcoord_xpix;
      double dy = ymean - corners[i].refcoord_ypix;

      if (dx < 0.0 && dy < 0.0) {
        cornerID[i] = 3;
      } else if (dx > 0.0 && dy < 0.0) {
        cornerID[i] = 2;
      } else if (dx > 0.0 && dy > 0.0) {
        cornerID[i] = 1;
      } else {
        cornerID[i] = 0;
      }
    }
  }

  void indentifyPixelLength(double actualDistance) {
    double dx_pix = sticks[1].refcoord_xpix - sticks[0].refcoord_xpix;
    double dy_pix = sticks[1].refcoord_ypix - sticks[0].refcoord_ypix;
    double pixDistance = sqrt(dx_pix * dx_pix + dy_pix * dy_pix);
    meterPerPixel = actualDistance / pixDistance;
  }
};

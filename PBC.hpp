#pragma once

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>

#include "Interaction.hpp"
#include "Loading.hpp"
#include "Particle.hpp"
#include "PeriodicCell.hpp"

class PBC {
 public:
  std::vector<Particle> Particles;
  std::vector<Interaction> Interactions;
  Loading Load;
  PeriodicCell Cell;
  double Sigxx, Sigxy, Sigyx, Sigyy;  // Internal stress

  struct ftbak_t {
    size_t i, j;
    double ft;
  };

  // Parameters
  double t{0.0};
  double tmax{5.0};
  double dt{1e-6};
  double interCloseC{0.0}, interClose{0.01}, dVerlet;
  double interOutC{0.0}, interOut{0.1};
  double interHistC{0.0}, interHist{0.25};
  double density{2700.0};
  double kn{1.e4};
  double kt{1.e4};
  double dampRate{0.95};
  double vBarrier{0.0};
  double numericalDampingCoeff{0.0};
  double vBarrierExpo{2.0};
  double mu{0.8};
  int iconf{0};
  int constrainedInFrame{0};

  // Methods
  PBC();
  void setSample();
  void integrate();
  void accelerations();
  void ResetCloseList(double dmax);
  void saveConf(int i);
  void saveConf(const char* name);
  void loadConf(const char* name);
};

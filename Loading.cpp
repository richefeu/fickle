#include <cstdio>

#include "Loading.hpp"

Loading::Loading() {}

void Loading::BiaxialCompression(double pressure, double velocity) {
  snprintf(StoredCommand, 256, "BiaxialCompression %g %g", pressure, velocity);
  xxDrive = true;
  yyDrive = false;
  xyDrive = yxDrive = false;
  Sigxx = pressure;
  Sigxy = Sigyx = Sigyy = 0.0;
  vyy = -velocity;
  vxy = vyx = 0.0;
  vxx = 0.0;  // free in fact
  ServoFunction = nullptr;
}

void Loading::IsostaticCompression(double pressure) {
  snprintf(StoredCommand, 256, "IsostaticCompression %g", pressure);
  xxDrive = yyDrive = ForceDriven;
  xyDrive = yxDrive = VelocityDriven;
  Sigxx = Sigyy = pressure;
  Sigxy = Sigyx = 0.0;
  vxx = vyy = 0.0;  // free in fact
  vxy = vyx = 0.0;
  ServoFunction = nullptr;
}

void Loading::SimpleShear(double pressure, double gammaDot) {
  snprintf(StoredCommand, 256, "SimpleShear %g %g", pressure, gammaDot);
  xxDrive = xyDrive = yxDrive = VelocityDriven;
  yyDrive = ForceDriven;
  Sigxx = Sigxy = Sigyx = 0.0;
  Sigyy = pressure;
  vxx = vyx = vyy = 0.0;
  vxy = 0.0;  // will be driven by the servoFunction
  ServoFunction = [gammaDot](Loading& load, PeriodicCell& cell) -> void { load.vxy = gammaDot * cell.h.yy; };
}

void Loading::VelocityControl(double Vxx, double Vxy, double Vyx, double Vyy) {
  snprintf(StoredCommand, 256, "VelocityControl %g %g %g %g", Vxx, Vxy, Vyx, Vyy);
  xxDrive = xyDrive = yxDrive = yyDrive = VelocityDriven;
  Sigxx = Sigxy = Sigyx = Sigyy = 0.0;
  vxx = Vxx;
  vxy = Vxy;
  vyx = Vyx;
  vyy = Vyy;
  ServoFunction = nullptr;
}

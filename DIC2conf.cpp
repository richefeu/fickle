#include <iostream>

#include "DICfile.hpp"
#include "PBC.hpp"

int main() {

  // parameters that are hard coded =====
  double woodDensity = 800.0;
  double scalingLength = 0.4;
  double cornerRadius = 0.025;
  double kn = 1000.0;
  double pressure = 1.0;
  double mu = 0.8;

  DICfile myDICfile;
  myDICfile.read("dic_out_780.txt", scalingLength);

  std::cout << "Number of particles read = " << myDICfile.particles.size() << std::endl;
  std::cout << "Number of corners read = " << myDICfile.corners.size() << std::endl;
  std::cout << "Number of sticks read = " << myDICfile.sticks.size() << std::endl;
  std::cout << "Number of fixed read = " << myDICfile.fixed.size() << std::endl;

  std::cout << "meterPerPixel = " << myDICfile.meterPerPixel << std::endl;

  std::cout << "corner 0 : " << myDICfile.getCorner(0).refcoord_xpix << ", " << myDICfile.getCorner(0).refcoord_ypix
            << std::endl;
  std::cout << "corner 1 : " << myDICfile.getCorner(1).refcoord_xpix << ", " << myDICfile.getCorner(1).refcoord_ypix
            << std::endl;
  std::cout << "corner 2 : " << myDICfile.getCorner(2).refcoord_xpix << ", " << myDICfile.getCorner(2).refcoord_ypix
            << std::endl;
  std::cout << "corner 3 : " << myDICfile.getCorner(3).refcoord_xpix << ", " << myDICfile.getCorner(3).refcoord_ypix
            << std::endl;

  PBC conf;

  // compute h
  conf.Cell.h.xx = myDICfile.getCorner(1).refcoord_xpix - myDICfile.getCorner(0).refcoord_xpix;
  conf.Cell.h.yx = -(myDICfile.getCorner(1).refcoord_ypix - myDICfile.getCorner(0).refcoord_ypix);
  conf.Cell.h.xy = myDICfile.getCorner(3).refcoord_xpix - myDICfile.getCorner(0).refcoord_xpix;
  conf.Cell.h.yy = -(myDICfile.getCorner(3).refcoord_ypix - myDICfile.getCorner(0).refcoord_ypix);
  conf.Cell.h *= myDICfile.meterPerPixel;
  std::cout << "h = " << conf.Cell.h << std::endl;

  mat4r hinv = conf.Cell.h.get_inverse();

  Particle P;
  double totalMass = 0.0;
  for (size_t i = 0; i < myDICfile.particles.size(); ++i) {
    vec2r pos;
    pos.x = myDICfile.particles[i].refcoord_xpix - myDICfile.getCorner(0).refcoord_xpix;
    pos.y = -(myDICfile.particles[i].refcoord_ypix - myDICfile.getCorner(0).refcoord_ypix);
    pos *= myDICfile.meterPerPixel;
    pos = hinv * pos;
    P.pos = pos;
    P.radius = myDICfile.particles[i].radius_pix * myDICfile.meterPerPixel;
    P.mass = M_PI * P.radius * P.radius * woodDensity;
    totalMass += P.mass;
    P.inertia = 0.5 * P.mass * P.radius * P.radius;
    conf.Particles.push_back(P);
  }

  // corner fake particle
  P.pos.reset();
  P.radius = cornerRadius;
  P.mass = M_PI * P.radius * P.radius * woodDensity;
  totalMass += P.mass;
  P.inertia = 0.5 * P.mass * P.radius * P.radius;
  conf.Particles.push_back(P);

  conf.Cell.mass = totalMass / sqrt((double)myDICfile.particles.size());
  conf.constrainedInFrame = 1;
  conf.Load.IsostaticCompression(pressure);

  double meanMass = totalMass / (double)myDICfile.particles.size();
  conf.dt = M_PI * sqrt(meanMass / kn) / 100.0;
  std::cout << "time step = " << conf.dt << std::endl;
  conf.tmax = 5.0;
  conf.interClose = conf.dt * 10.0;
  conf.interOut = conf.tmax * 0.001;
  conf.interHist = conf.tmax * 0.01;
  conf.density = woodDensity;
  conf.kn = kn;
  conf.kt = kn;
  conf.dampRate = 0.95;
  conf.numericalDampingCoeff = 0.7;
  conf.mu = mu;
  conf.iconf = 0;
  conf.constrainedInFrame = 1;

  conf.saveConf("input.txt");

  return 0;
}
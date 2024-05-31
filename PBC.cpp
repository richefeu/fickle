// Periodic boundary conditions with disks
// This application is developped to investigate fluctuations

#include "PBC.hpp"

PBC::PBC() {
  // default values?
}

void PBC::saveConf(int i) {
  char fname[256];
  snprintf(fname, 256, "conf%d", i);
  saveConf(fname);
}

void PBC::saveConf(const char* fname) {
  std::ofstream conf(fname);

  conf << "PBC 27-11-2014" << std::endl;  // format: progName version-date
  conf << "t " << t << std::endl;
  conf << "tmax " << tmax << std::endl;
  conf << "dt " << dt << std::endl;
  conf << "interClose " << interClose << std::endl;
  conf << "interOut " << interOut << std::endl;
  conf << "interHist " << interHist << std::endl;
  conf << "dVerlet " << dVerlet << std::endl;
  conf << "constrainedInFrame " << constrainedInFrame << std::endl;
  conf << "density " << density << std::endl;
  conf << "kn " << kn << std::endl;
  conf << "kt " << kt << std::endl;
  conf << "dampRate " << dampRate << std::endl;
  conf << "mu " << mu << std::endl;
  conf << "iconf " << iconf << std::endl;
  conf << "h " << Cell.h << std::endl;
  conf << "vh " << Cell.vh << std::endl;
  conf << "ah " << Cell.ah << std::endl;
  conf << "hmass " << Cell.mass << std::endl;
  conf << "Load " << Load.StoredCommand << std::endl;
  conf << std::setprecision(15);
  conf << "Particles " << Particles.size() << std::endl;
  for (size_t i = 0; i < Particles.size(); i++) {
    conf << Particles[i].pos << " " << Particles[i].vel << " " << Particles[i].acc << " " << Particles[i].rot << " "
         << Particles[i].vrot << " " << Particles[i].arot << " " << Particles[i].radius << " " << Particles[i].inertia
         << " " << Particles[i].mass << std::endl;
  }
  conf << "Interactions " << Interactions.size() << std::endl;
  for (size_t i = 0; i < Interactions.size(); i++) {
    if (fabs(Interactions[i].fn) < 1.0e-20) continue;
    conf << Interactions[i].i << " " << Interactions[i].j << " " << Interactions[i].fn << " " << Interactions[i].ft
         << " " << Interactions[i].damp << std::endl;
  }
}

void PBC::loadConf(const char* name) {
  std::ifstream conf(name);
  if (!conf.is_open()) {
    std::cerr << "Cannot read " << name << std::endl;
  }

  // Check header
  std::string prog;
  conf >> prog;
  if (prog != "PBC") {
    std::cerr << "This is not file for PBC executable!" << std::endl;
  }
  std::string date;
  conf >> date;
  if (date != "27-11-2014") {
    std::cerr << "The version-date should be 27-11-2014!" << std::endl;
  }

  std::string token;
  conf >> token;
  while (conf.good()) {
    if (token == "t") {
      conf >> t;
    } else if (token == "tmax") {
      conf >> tmax;
    } else if (token == "dt") {
      conf >> dt;
    } else if (token == "interClose") {
      conf >> interClose;
    } else if (token == "interOut") {
      conf >> interOut;
    } else if (token == "interHist") {
      conf >> interHist;
    } else if (token == "dVerlet") {
      conf >> dVerlet;
    } else if (token == "density") {
      conf >> density;
    } else if (token == "constrainedInFrame") {
      conf >> constrainedInFrame;
    } else if (token == "kn") {
      conf >> kn;
    } else if (token == "kt") {
      conf >> kt;
    } else if (token == "dampRate") {
      conf >> dampRate;
    } else if (token == "vBarrier") {
      conf >> vBarrier;
      vBarrier = fabs(vBarrier);
    } else if (token == "vBarrierExpo") {
      conf >> vBarrierExpo;
      vBarrierExpo = fabs(vBarrierExpo);
    } else if (token == "numericalDampingCoeff") {
      conf >> numericalDampingCoeff;
    } else if (token == "mu") {
      conf >> mu;
    } else if (token == "iconf") {
      conf >> iconf;
    } else if (token == "h") {
      conf >> Cell.h;
    } else if (token == "vh") {
      conf >> Cell.vh;
    } else if (token == "ah") {
      conf >> Cell.ah;
    } else if (token == "hmass") {
      conf >> Cell.mass;
    } else if (token == "Load") {
      std::string command;
      conf >> command;
      if (command == "BiaxialCompression") {
        double pressure, velocity;
        conf >> pressure >> velocity;
        Load.BiaxialCompression(pressure, velocity);
      } else if (command == "IsostaticCompression") {
        double pressure;
        conf >> pressure;
        Load.IsostaticCompression(pressure);
      } else if (command == "VelocityControl") {
        double vhxx, vhxy, vhyx, vhyy;
        conf >> vhxx >> vhxy >> vhyx >> vhyy;
        Load.VelocityControl(vhxx, vhxy, vhyx, vhyy);
      } else if (command == "SimpleShear") {
        double pressure, gammaDot;
        conf >> pressure >> gammaDot;
        Load.SimpleShear(pressure, gammaDot);
      } else {
        std::cerr << "Unknown command for loading: " << command << std::endl;
      }
    } else if (token == "Particles") {
      size_t nb;
      conf >> nb;
      Particles.clear();
      Particle P;
      for (size_t i = 0; i < nb; i++) {
        conf >> P.pos >> P.vel >> P.acc >> P.rot >> P.vrot >> P.arot >> P.radius >> P.inertia >> P.mass;
        Particles.push_back(P);
      }
    } else if (token == "Interactions") {
      size_t nb;
      conf >> nb;
      Interactions.clear();
      Interaction I;
      for (size_t i = 0; i < nb; i++) {
        conf >> I.i >> I.j >> I.fn >> I.ft >> I.damp;
        Interactions.push_back(I);
      }
    } else {
      std::cerr << "Unknown token: " << token << std::endl;
    }

    conf >> token;
  }
}

void PBC::setSample() {
  Particle P;
  P.rot = 0.0;

  int ngw = 15;
  std::cout << std::endl;
  std::cout << "Nombre de spheres sur la largeur : ";
  std::cin >> ngw;
  double step = 1.0 / (2.0 * ngw);  // dans le domaine 01

  double radius = 1e-3;
  std::cout << "Rayon max : ";
  std::cin >> radius;
  double deltaR = 0.2e-3;
  std::cout << "Delta rayon : ";
  std::cin >> deltaR;
  vec2r from(0.0, 0.0);
  vec2r to(2.0 * (ngw + 1) * radius, 2.0 * (ngw + 1) * radius);
  dVerlet = 0.95 * (radius - deltaR);

  int i = 0;
  double massTot = 0.0;
  while (P.pos.y <= 1.0) {
    P.radius = radius - deltaR * (static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
    P.mass = M_PI * P.radius * P.radius * density;
    massTot += P.mass;
    P.inertia = 0.5 * P.mass * P.radius * P.radius;
    int column = i % ngw;
    int row = i / ngw;
    if (row % 2 == 0) {  // even row
      P.pos.x = step + 2 * column * step;
    } else {  // odd row
      P.pos.x = 2 * step + 2 * column * step;
    }
    P.pos.y = step + 2 * row * step;
    if (P.pos.y <= 1. - step) {
      Particles.push_back(P);
    }
    i++;
  }

  Cell.Define(to.x, from.y, from.x, to.y);

  Cell.mass = massTot;  // / sqrt((double)Particles.size());
  // Cell.mass corresponds to nearly the mass of a line or column of particles

  // Set ramdom velocities to particles
  double vmax = 0.1;
  std::cout << "Norme de la vitesse aleatoire : ";
  std::cin >> vmax;
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].vel.x = vmax * (static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
    Particles[i].vel.y = vmax * (static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
  }

  double press = 1000.0;
  std::cout << "Pression de confinement isostatique : ";
  std::cin >> press;
  Load.IsostaticCompression(press);

  // Parametres de simu

  density = 2700.0;
  std::cout << "densite des particules : ";
  std::cin >> density;

  kn = 1.e4;
  while (true) {
    char rep = 'o';
    std::cout << "kn : ";
    std::cin >> kn;
    std::cout << "ce qui donne kn/p = " << kn / press << std::endl;
    std::cout << "satisfait ? (o/n) : ";
    std::cin >> rep;
    if (rep == 'o' || rep == 'O') break;
  }

  double ktkn = 1.0;
  std::cout << "kt/kn : ";
  std::cin >> ktkn;
  kt = ktkn * kn;

  dampRate = 0.95;
  std::cout << "Taux d'amortissement : ";
  std::cin >> dampRate;

  mu = 0.8;
  std::cout << "Coefficient de frottement : ";
  std::cin >> mu;

  dt = 1e-6;
  while (true) {
    char rep = 'o';
    std::cout << "dt : ";
    std::cin >> dt;
    double m = M_PI * radius * radius * density;
    std::cout << "ce qui donne dt_crit/dt = " << (sqrt(m / kn)) / dt << std::endl;
    std::cout << "satisfait ? (o/n) : ";
    std::cin >> rep;
    if (rep == 'o' || rep == 'O') break;
  }

  t = 0.0;
  iconf = 0;
  interCloseC = 0.0;
  interOutC = 0.0;
  interHistC = 0.0;
  tmax = 5.0;
  std::cout << "durée : ";
  std::cin >> tmax;

  // TODO
  interClose = 0.01;
  interOut = 0.1;
  interHist = 0.25;

  return;
}

// The integration loop (Velocity Verlet)
void PBC::integrate() {
  double dt_2 = 0.5 * dt;
  double dt2_2 = 0.5 * dt * dt;

  ResetCloseList(dVerlet);
  saveConf(iconf);

  std::ofstream fileOut("output.txt");
  while (t <= tmax) {

    if (Load.ServoFunction != nullptr) {
      Load.ServoFunction(Load, Cell);
    }

    for (size_t i = 0; i < Particles.size(); i++) {
      Particles[i].pos += dt * Particles[i].vel + dt2_2 * Particles[i].acc;
      if (constrainedInFrame == 1) {
        while (Particles[i].pos.x < 0.0) {
          Particles[i].pos.x += 1.0;
        }
        while (Particles[i].pos.x > 1.0) {
          Particles[i].pos.x -= 1.0;
        }
        while (Particles[i].pos.y < 0.0) {
          Particles[i].pos.y += 1.0;
        }
        while (Particles[i].pos.y > 1.0) {
          Particles[i].pos.y -= 1.0;
        }
      }
      Particles[i].vel += dt_2 * Particles[i].acc;

      Particles[i].rot += dt * Particles[i].vrot + dt2_2 * Particles[i].arot;
      Particles[i].vrot += dt_2 * Particles[i].arot;
    }

    if (Load.xxDrive == ForceDriven) {
      Cell.h.xx += dt * Cell.vh.xx + dt2_2 * Cell.ah.xx;
      Cell.vh.xx += dt_2 * Cell.ah.xx;
    } else {
      Cell.h.xx += dt * Load.vxx;
      Cell.vh.xx = Load.vxx;
      Cell.ah.xx = 0.0;
    }

    if (Load.xyDrive == ForceDriven) {
      Cell.h.xy += dt * Cell.vh.xy + dt2_2 * Cell.ah.xy;
      Cell.vh.xy += dt_2 * Cell.ah.xy;
    } else {
      Cell.h.xy += dt * Load.vxy;
      Cell.vh.xy = Load.vxy;
      Cell.ah.xy = 0.0;
    }

    if (Load.yxDrive == ForceDriven) {
      Cell.h.yx += dt * Cell.vh.yx + dt2_2 * Cell.ah.yx;
      Cell.vh.yx += dt_2 * Cell.ah.yx;
    } else {
      Cell.h.yx += dt * Load.vyx;
      Cell.vh.yx = Load.vyx;
      Cell.ah.yx = 0.0;
    }

    if (Load.yyDrive == ForceDriven) {
      Cell.h.yy += dt * Cell.vh.yy + dt2_2 * Cell.ah.yy;
      Cell.vh.yy += dt_2 * Cell.ah.yy;
    } else {
      Cell.h.yy += dt * Load.vyy;
      Cell.vh.yy = Load.vyy;
      Cell.ah.yy = 0.0;
    }

    accelerations();

    // Limit the velocity of cell boundaries (vertical and horizontal only)
    // useful only for preparation
    if (vBarrier > 0.0) {
      double xxratio = pow(fabs(Cell.vh.xx / vBarrier), vBarrierExpo);
      Cell.ah.xx *= (1.0 - xxratio) / (1.0 + xxratio);
      double yyratio = pow(fabs(Cell.vh.yy / vBarrier), vBarrierExpo);
      Cell.ah.yy *= (1.0 - yyratio) / (1.0 + yyratio);
    }

    // Cundall damping (for preparation)
    if (numericalDampingCoeff > 0.0) {
      for (size_t i = 0; i < Particles.size(); i++) {
        double factor;
        double factorMinus = 1.0 - numericalDampingCoeff;
        double factorPlus = 1.0 + numericalDampingCoeff;

        vec2r vel = Cell.h * Particles[i].vel + Cell.vh * Particles[i].pos;
        factor = (Particles[i].force * vel > 0.0) ? factorMinus : factorPlus;
        Particles[i].force *= factor;
      }
    }

    vec2r vmean;
    for (size_t i = 0; i < Particles.size(); i++) {
      Particles[i].vel += dt_2 * Particles[i].acc;
      vmean += Particles[i].vel;
      Particles[i].vrot += dt_2 * Particles[i].arot;
    }
    vmean /= (double)(Particles.size());
    // substract the mean reduced velocity
    for (size_t i = 0; i < Particles.size(); i++) {
      Particles[i].vel -= vmean;
    }

    if (Load.xxDrive == ForceDriven) {
      Cell.vh.xx += dt_2 * Cell.ah.xx;
    }
    if (Load.xyDrive == ForceDriven) {
      Cell.vh.xy += dt_2 * Cell.ah.xy;
    }
    if (Load.yxDrive == ForceDriven) {
      Cell.vh.yx += dt_2 * Cell.ah.yx;
    }
    if (Load.yyDrive == ForceDriven) {
      Cell.vh.yy += dt_2 * Cell.ah.yy;
    }

    interHistC += dt;
    interOutC += dt;
    interCloseC += dt;
    t += dt;

    // ----

    if (interCloseC >= interClose - dt_2) {
      ResetCloseList(dVerlet);
      interCloseC = 0.0;
    }

    if (interOutC >= interOut - dt_2) {
      fileOut << t << " " << Cell.h << " " << Sigxx << " " << Sigxy << " " << Sigyx << " " << Sigyy << std::endl;
      interOutC = 0.0;
    }

    if (interHistC >= interHist - dt_2) {
      iconf++;
      std::cout << "iconf = " << iconf << ", Time = " << t << std::endl;
      saveConf(iconf);
      interHistC = 0.0;
    }
  }

  return;
}

void PBC::ResetCloseList(double dmax) {
  // store ft because the list will be cleared before being rebuilt
  std::vector<ftbak_t> ft_backup;
  ftbak_t I;
  for (size_t k = 0; k < Interactions.size(); k++) {
    I.i = Interactions[k].i;
    I.j = Interactions[k].j;
    I.ft = Interactions[k].ft;
    ft_backup.push_back(I);
  }

  // now rebuild the list
  Interactions.clear();
  for (size_t i = 0; i < Particles.size(); i++) {
    for (size_t j = i + 1; j < Particles.size(); j++) {

      vec2r sij = Particles[j].pos - Particles[i].pos;
      sij.x -= floor(sij.x + 0.5);
      sij.y -= floor(sij.y + 0.5);
      vec2r branch = Cell.h * sij;

      double sum = dmax + Particles[i].radius + Particles[j].radius;
      if (norm2(branch) <= sum * sum) {
        double m = (Particles[i].mass * Particles[j].mass) / (Particles[i].mass + Particles[j].mass);
        double Damp = dampRate * 2.0 * sqrt(kn * m);
        Interactions.push_back(Interaction(i, j, Damp));
      }
    }
  }

  // retrieve ft values
  size_t k, kold = 0;
  for (k = 0; k < Interactions.size(); ++k) {
    while (kold < ft_backup.size() && ft_backup[kold].i < Interactions[k].i) ++kold;
    if (kold == ft_backup.size()) break;

    while (kold < ft_backup.size() && ft_backup[kold].i == Interactions[k].i && ft_backup[kold].j < Interactions[k].j)
      ++kold;
    if (kold == ft_backup.size()) break;

    if (ft_backup[kold].i == Interactions[k].i && ft_backup[kold].j == Interactions[k].j) {
      Interactions[k].ft = ft_backup[kold].ft;
      ++kold;
    }
  }
}

void PBC::accelerations() {
  // Set forces and moments to zero
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].force.reset();
    Particles[i].moment = 0.0;
  }
  Sigxx = Sigxy = Sigyx = Sigyy = 0.0;

  size_t i, j;
  for (size_t k = 0; k < Interactions.size(); k++) {
    i = Interactions[k].i;
    j = Interactions[k].j;

    vec2r sij = Particles[j].pos - Particles[i].pos;
    sij.x -= floor(sij.x + 0.5);
    sij.y -= floor(sij.y + 0.5);
    vec2r branch = Cell.h * sij;

    double sum = Particles[i].radius + Particles[j].radius;
    if (norm2(branch) <= sum * sum) {  // it means that i and j are in contact

      vec2r n = branch;
      double len = n.normalize();
      vec2r t(-n.y, n.x);

      // real relative velocities
      vec2r vel = Particles[j].vel - Particles[i].vel;
      double dn = len - Particles[i].radius - Particles[j].radius;
      vec2r realVel = Cell.h * vel + Cell.vh * sij;
      double Bi = Particles[i].radius + 0.5 * dn;
      double Bj = Particles[j].radius + 0.5 * dn;

      // Normal force (elastic + viscuous)
      double vn = realVel * n;
      double fne = -kn * dn;
      double fnv = -Interactions[k].damp * vn;
      Interactions[k].fn = fne + fnv;

      // Tangential force (friction)
      double vijt = realVel * t - Particles[i].vrot * Bi - Particles[j].vrot * Bj;

      double ft = Interactions[k].ft - kt * dt * vijt;
      double ftest = mu * fne;
      if (fabs(ft) > ftest) {
        ft = (ft > 0.0) ? ftest : -ftest;
      }
      Interactions[k].ft = ft;

      // Resultant force and moment
      vec2r f = Interactions[k].fn * n + Interactions[k].ft * t;
      Particles[i].force -= f;
      Particles[j].force += f;
      Particles[i].moment -= ft * Bi;
      Particles[j].moment -= ft * Bj;

      // Internal stress
      Sigxx += f.x * branch.x;
      Sigxy += f.x * branch.y;
      Sigyx += f.y * branch.x;
      Sigyy += f.y * branch.y;
    }
  }  // Loop over interactions

  double invV = 1.0 / Cell.h.det();
  Sigxx *= invV;
  Sigxy *= invV;
  Sigyx *= invV;
  Sigyy *= invV;

  // Finally compute the accelerations (translation and rotation)
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].acc = Particles[i].force / Particles[i].mass;

    // The following 4 lines may (MUST!!!!) be removed in fact
    /*
    Particles[i].acc.x -= 2.0 * (Cell.vh.xx * Particles[i].vel.x + Cell.vh.xy * Particles[i].vel.y);
    Particles[i].acc.y -= 2.0 * (Cell.vh.yx * Particles[i].vel.x + Cell.vh.yy * Particles[i].vel.y);
    Particles[i].acc.x -= Cell.ah.xx * Particles[i].pos.x + Cell.ah.xy * Particles[i].pos.y;
    Particles[i].acc.y -= Cell.ah.yx * Particles[i].pos.x + Cell.ah.yy * Particles[i].pos.y;
    */

    // Compute inverse of the matrix h
    double invDet = 1.0 / (Cell.h.xx * Cell.h.yy - Cell.h.yx * Cell.h.xy);
    double hxxinv = invDet * Cell.h.yy;
    double hxyinv = -invDet * Cell.h.xy;
    double hyxinv = -invDet * Cell.h.yx;
    double hyyinv = invDet * Cell.h.xx;

    vec2r acc = Particles[i].acc;
    Particles[i].acc.x = hxxinv * acc.x + hxyinv * acc.y;
    Particles[i].acc.y = hyxinv * acc.x + hyyinv * acc.y;

    Particles[i].arot = Particles[i].moment / Particles[i].inertia;
  }

  if (Load.xxDrive == ForceDriven) {
    Cell.ah.xx = ((Sigxx - Load.Sigxx) * Cell.h.yy - (Sigyx - Load.Sigyx) * Cell.h.xy) / Cell.mass;
  }
  if (Load.xyDrive == ForceDriven) {
    Cell.ah.xy = ((Sigxy - Load.Sigxy) * Cell.h.yy - (Sigyy - Load.Sigyy) * Cell.h.xy) / Cell.mass;
  }
  if (Load.yxDrive == ForceDriven) {
    Cell.ah.yx = ((Sigyx - Load.Sigyx) * Cell.h.xx - (Sigxx - Load.Sigxx) * Cell.h.yx) / Cell.mass;
  }
  if (Load.yyDrive == ForceDriven) {
    Cell.ah.yy = ((Sigyy - Load.Sigyy) * Cell.h.xx - (Sigxy - Load.Sigxy) * Cell.h.yx) / Cell.mass;
  }
}

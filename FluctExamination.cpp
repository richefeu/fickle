#include "FluctExamination.hpp"

void concatFluctuations(int begin, int end, int step, std::vector<double>& xflucts, std::vector<double>& yflucts) {
  PBC confFrom;
  PBC confTo;
  char fnameFrom[128];
  char fnameTo[128];

  xflucts.clear();
  yflucts.clear();
  for (int ifile = begin; ifile <= end - step; ifile++) {
    snprintf(fnameFrom, 128, "conf%d", ifile);
    snprintf(fnameTo, 128, "conf%d", ifile + step);
    confFrom.loadConf(fnameFrom);
    confTo.loadConf(fnameTo);

    for (size_t i = 0; i < confFrom.Particles.size(); ++i) {
      vec2r delta = confTo.Particles[i].pos - confFrom.Particles[i].pos;
      delta.x -= floor(delta.x + 0.5);
      delta.y -= floor(delta.y + 0.5);
      vec2r fluct = confTo.Cell.h * delta;
      xflucts.push_back(fluct.x);
      yflucts.push_back(fluct.y);
    }
  }
}

void getDistribx(std::ostream& os, int begin, int end, int step) {
  std::vector<double> xflucts, yflucts;

  std::cout << "Concatenize the fluctuations from file " << begin << " to " << end << ", step " << step << std::endl;
  concatFluctuations(begin, end, step, xflucts, yflucts);

  CDF xcdf(xflucts, 1000);

  std::cout << "Optimisation..." << std::endl;
  size_t i25 = xcdf.cdens.size() / 4;
  size_t i50 = 2 * i25;
  size_t di = i25 / 2;
  std::cout << std::setprecision(15);

  double qini = 1.05;
  double Nini = (xcdf.cdens[i50 + di] - xcdf.cdens[i50 - di]) / (xcdf.value[i50 + di] - xcdf.value[i50 - di]);
  double Bini = Nini * Nini * 3.14159;
  std::cout << "q_ini = " << qini << ", beta_ini = " << Bini << ", normalizer_ini = " << Nini << std::endl;
  double fret = xcdf.findQGaussParameters(qini, Bini, Nini);
  std::cout << "q = " << xcdf.q << ", beta = " << xcdf.beta << ", normalizer = " << xcdf.normalizer << std::endl;

  leastSquaredDistance_QGauss func;
  func.xvalue = &xcdf.value[0];
  func.yvalue = &xcdf.cdens[0];
  func.minValue = std::min(2.0 * xcdf.value[0], -2.0 * xcdf.value.back());

  double integrale = func.qgauss_integrate(func.minValue, -func.minValue, xcdf.q, xcdf.beta, xcdf.normalizer, 1000);
  std::cout << "fret = " << fret << ", integrale = " << integrale << std::endl;

  os << step << ' ' << xcdf.q << ' ' << xcdf.beta << ' ' << xcdf.normalizer << ' ' << fret << ' ' << integrale
     << std::endl;

  char xfname[256];
  snprintf(xfname, 256, "cdfx_%d.txt", step);
  std::ofstream xfile(xfname);
  int deltai = 1;
  int ilast = (int)xcdf.value.size() - 1;
  for (size_t i = 0; i < xcdf.value.size(); i++) {
    int iprev = std::max(0, (int)i - deltai);
    int inext = std::min(ilast, (int)i + deltai);
    xfile << xcdf.value[i] << ' ' << xcdf.cdens[i] << ' '
          << func.qgauss_integrate(func.minValue, xcdf.value[i], qini, Bini, Nini, 1000) << ' '
          << func.qgauss_integrate(func.minValue, xcdf.value[i], xcdf.q, xcdf.beta, xcdf.normalizer, 1000) << ' '
          << (xcdf.cdens[inext] - xcdf.cdens[iprev]) / (xcdf.value[inext] - xcdf.value[iprev]) << ' '
          << func.qgauss(xcdf.value[i], xcdf.q, xcdf.beta, xcdf.normalizer) << std::endl;
  }
}

int main(int argc, char const* argv[]) {

  if (argc == 6) {
    int beg = atoi(argv[1]);
    int end = atoi(argv[2]);
    int step_beg = atoi(argv[3]);
    int step_end = atoi(argv[4]);
    int step_inc = atoi(argv[5]);

    std::ofstream file("qfitx.txt");
    file << std::setprecision(15);
    for (int step = step_beg; step <= step_end; step += step_inc) {
      getDistribx(file, beg, end, step);
    }

  } else {
    std::cout << "usage: " << argv[0] << " <beging number> <end number> <step_begin> <step_end> <step_inc>\n";
  }

  return 0;
}
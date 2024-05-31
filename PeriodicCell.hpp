#ifndef PERIODICCELL_HPP
#define PERIODICCELL_HPP

#include "mat4.hpp"

struct PeriodicCell {
  mat4r h;
  mat4r vh;
  mat4r ah;

  double mass;

  PeriodicCell();
  PeriodicCell(double a1x, double a1y, double a2x, double a2y);
  void Define(double a1x, double a1y, double a2x, double a2y);
};

#endif /* end of include guard: PERIODICCELL_HPP */

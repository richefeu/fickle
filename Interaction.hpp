#ifndef INTERACTION_HPP
#define INTERACTION_HPP

#include <cstddef>

struct Interaction {
  size_t i, j;
  double fn, ft;
  double damp;

  Interaction();
  Interaction(size_t I, size_t J, double Damp);
};

#endif /* end of include guard: INTERACTION_HPP */

#pragma once

#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>

#include "PBC.hpp"
#include "CDF.hpp"

void concatFluctuations(int begin, int end, int step, std::vector<vec2r>& flucts);
void getDistribx(std::ostream & os, int begin, int end, int step);
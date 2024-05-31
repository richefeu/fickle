#ifndef LOADING_HPP
#define LOADING_HPP

#include <functional>

#include "PeriodicCell.hpp"

const bool ForceDriven = true;
const bool VelocityDriven = false;

// Loading applied to collective degrees of freedom
struct Loading {
	// 1 for pressure and 0 for velocity
	bool xxDrive, xyDrive;
	bool yxDrive, yyDrive;

	// imposed external stress
	double Sigxx, Sigxy;
	double Sigyx, Sigyy;

	// imposed velocities (TODO rename vhxx...)
	double vxx, vxy;
	double vyx, vyy;

	char StoredCommand[256];

	// This function will be set to a lambda (c++11)
	std::function<void(Loading & load, PeriodicCell & cell)> ServoFunction;

	Loading();
	void BiaxialCompression(double pressure, double velocity);
	void IsostaticCompression(double pressure);
	void SimpleShear(double pressure, double gammaDot);
	void VelocityControl(double Vxx, double Vxy, double Vyx, double Vyy);
};


#endif /* end of include guard: LOADING_HPP */

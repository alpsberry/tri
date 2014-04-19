#ifndef TRI_DGPROBLEM_H
#define TRI_DGPROBLEM_H

#include <iostream>
#include "problem.h"

// the problem to be solved here is
// -\Delta u + u = 3\cos x \sin y, (x,y) \in \Omega = [\frac{\pi}{2}, \frac{3\pi}{2}] \times [0, \pi]
// u = 0, (x,y) \in \Gamma
class DGProblem: public Problem
{
public:
	double f(double x, double y){
		return 3 * cos(x) * sin(y);
		// return 1.0;
	}
	// read parameters from an input file
	int initProblem(int argc, char const *argv[]);
};

#endif

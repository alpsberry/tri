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
	double epsilon;
	double sigma0;
	double beta0;

	double f(double x, double y){
		return 3 * cos(x) * sin(y);
		// return 1.0;
	}
	double trueSol(double x, double y){
		return cos(x) * sin(y);
	}
	// read parameters from an input file
	int initDGProblem(int argc, char const *argv[]);
};

#endif

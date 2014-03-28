#ifndef TRI_FEMSOLVESYS_H
#define TRI_FEMSOLVESYS_H

#include "solveSys.h"
// the problem to be solved here is
//   -\Delta u + u = 3\cos x \sin y, (x,y) \in \Omega = [\frac{\pi}{2}, \frac{3\pi}{2}] \times [0, \pi]
//   u = 0, (x,y) \in \Gamma
// the variational form
//   \int\nolimits_{\Omega} \nabla u \nabla v + uv = \int\nolimits_{\Omega} fv
class FEMSolvingSystem: public BasicSolvingSystem
{
public:
	double integA(Element &ele, Mesh mesh, int vi, int vj, MyProblem prob);

	double integARH(Element ele, Mesh mesh, int vi, MyProblem prob);

	int getStiff(Element& ele, Mesh& mesh, MyProblem prob);

	int assembleStiff(Mesh &mesh, MyProblem prob);

	int consoleOutput(Mesh mesh);

	int fileOutput(Mesh mesh, MyProblem prob);

	int triOutput(MyProblem prob, Mesh mesh);

};

#endif /* TRI_FEMSOLVESYS_H */

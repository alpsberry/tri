#ifndef TRI_SOLVESYS_H
#define TRI_SOLVESYS_H

#include <ctime>
#include <vector>
#include <list>
#include "umfpack.h"
#include "mesh.h"
#include "../SuperLU_4.3/SRC/slu_ddefs.h"


struct maColEle
{
	int row;
	double value;
};

bool maCompareFunc(const maColEle t1,const maColEle t2);

class BasicSolvingSystem
{
public:
	std::vector< std::list<maColEle> > ma;
	double* rh; //right hand vector

	//for UMFPACK
	int* Ap; //Ap[0] = 0; Ap[k] num of nonzero entries in the first k columns
	int* Ai; //row of each nonzero entry, column-wise
	double* Ax; //value of each nonzero entry, column-wise
	
	double* x;

	virtual int assembleStiff(Mesh &mesh, MyProblem prob) = 0;
	virtual int getStiff(Element&, Mesh&, MyProblem prob) = 0;
	virtual double integElement(Element ele){
		return 0;
	}

	int solveSparse(Mesh mesh, MyProblem prob);

	int convertToCSC(Mesh mesh);

	int UMFSolve(Mesh mesh);

	int SuperLUSolve(Mesh mesh);

	int addToMA(double a, int row, int col);
	
	// ~BasicSolvingSystem() {} //need to deal with rh, ma, Ap... ?
};

// the problem to be solved here is
//   -\Delta u + u = 3\cos x \sin y, (x,y) \in \Omega = [\frac{\pi}{2}, \frac{3\pi}{2}] \times [0, \pi]
//   u = 0, (x,y) \in \Gamma
// the variational form
//   \int\nolimits_{\Omega} \nabla u \nabla v + uv = \int\nolimits_{\Omega} fv
class MySolvingSystem: public BasicSolvingSystem
{
public:
	double integA(Element &ele, Mesh mesh, int vi, int vj, MyProblem prob);

	double integARH(Element ele, Mesh mesh, int vi, MyProblem prob);

	int getStiff(Element& ele, Mesh& mesh, MyProblem prob);

	int assembleStiff(Mesh &mesh, MyProblem prob);
};

#endif

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

	int fileOutputTriplet(Mesh mesh, MyProblem prob);
	int fileOutputRH(Mesh mesh, MyProblem prob);
	int fileOutputMA(Mesh mesh, MyProblem prob);
	
	// ~BasicSolvingSystem() {} //need to deal with rh, ma, Ap... ?
};

#endif

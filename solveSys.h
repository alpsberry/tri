#ifndef TRI_SOLVESYS_H
#define TRI_SOLVESYS_H

#include <ctime>
#include <vector>
#include <list>
// #include "SuiteSparse_config/SuiteSparse_config.h"
#include "umfpack.h"
#include "mesh.h"

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

	int convertToUMF(Mesh mesh);

	int UMFSolve(Mesh mesh);

	int addToMA(double a, int row, int col);

};


class MySolvingSystem: public BasicSolvingSystem
{
public:
	double integA(Element &ele, Mesh mesh, int vi, int vj, MyProblem prob);

	double integARH(Element ele, Mesh mesh, int vi, MyProblem prob);

	int getStiff(Element& ele, Mesh& mesh, MyProblem prob);

	int assembleStiff(Mesh &mesh, MyProblem prob);

	~MySolvingSystem()
	{
		// delete [] rh;
		//need to deal with ma
	}
};

#endif

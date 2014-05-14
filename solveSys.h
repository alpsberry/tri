#ifndef TRI_SOLVESYS_H
#define TRI_SOLVESYS_H

#define __SOLVESYS_DEBUG

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

typedef std::vector< std::vector<double> > VECMATRIX;

bool maCompareFunc(const maColEle t1,const maColEle t2);

template <typename MyProblem>
class BasicSolvingSystem
{
public:
	int dof;

	std::vector< std::list<maColEle> > ma;
	double* rh; //right hand vector

	//for UMFPACK
	int* Ap; //Ap[0] = 0; Ap[k] num of nonzero entries in the first k columns
	int* Ai; //row of each nonzero entry, column-wise
	double* Ax; //value of each nonzero entry, column-wise
	
	double* x;

	// virtual int assembleStiff(Mesh<MyProblem> &mesh, MyProblem prob) = 0;
	// virtual int getStiff(Element&, Mesh&, MyProblem prob) = 0;
	virtual double integElement(Element ele){
		return 0;
	}

	int solveSparse(Mesh<MyProblem> mesh, MyProblem prob);

	int convertToCSC(Mesh<MyProblem> mesh);

	int UMFSolve(Mesh<MyProblem> mesh);

	int SuperLUSolve(Mesh<MyProblem> mesh);

	int addToMA(double a, int row, int col);

	int fileOutputTriplet(Mesh<MyProblem> mesh, MyProblem prob);
	int fileOutputRH(Mesh<MyProblem> mesh, MyProblem prob);
	int fileOutputMA(Mesh<MyProblem> mesh, MyProblem prob);
	
	// ~BasicSolvingSystem() {} //need to deal with rh, ma, Ap... ?
};

bool maCompareFunc(const maColEle t1, const maColEle t2){  
	return t1.row < t2.row;  
} 
template <typename MyProblem>
int BasicSolvingSystem<MyProblem>::solveSparse(Mesh<MyProblem> mesh, MyProblem prob)
{
	convertToCSC(mesh);
	if(prob.parameters.solPack == SolPack::UMFPACK)
		UMFSolve(mesh);
	else
		if (prob.parameters.solPack == SolPack::SuperLU)
			SuperLUSolve(mesh);

	return 0;
}

// convert the list-stored matrix to Compressed Sparse Column (CSC)
template <typename MyProblem>
int BasicSolvingSystem<MyProblem>::convertToCSC(Mesh<MyProblem> mesh)
{
#ifdef __SOLVESYS_DEBUG
	std::cout << "start converting to UMFPACK structure" << std::endl;
#endif

	//sort entries in each column by their row
	std::vector< std::list<maColEle> >::iterator it;
	for(it = ma.begin(); it != ma.end(); it++)
		it -> sort(maCompareFunc);
	
	Ap = new int [ma.size() + 1];
	Ap[0] = 0;
	int nnz(0), k(0);
	for( it = ma.begin(), k = 1;
		it != ma.end(); it++, k++)
	{
		nnz += it -> size();
		Ap[k] = nnz;
	}

	Ai = new int [nnz];
	Ax = new double [nnz];
	for(it = ma.begin(), k = 0; it != ma.end(); it++)
		for(std::list<maColEle>::iterator it1 = it -> begin();
			it1 != it -> end(); it1++, k++)
		{
			Ai[k] = it1 -> row;
			Ax[k] = it1 -> value;
		}

#ifdef __SOLVESYS_DEBUG
	std::cout << "finish converting to UMFPACK structure" << std::endl << std::endl;
#endif

	return 0;
}
template <typename MyProblem>
int BasicSolvingSystem<MyProblem>::UMFSolve(Mesh<MyProblem> mesh)
{

#ifdef __SOLVESYS_DEBUG
	std::cout << "start solving with UMFPACK" << std::endl;
#endif

	x = new double [dof];
	memset(x, 0, dof * sizeof(double));
	void *Symbolic, *Numeric ;
	(void) umfpack_di_symbolic (dof, dof, Ap, Ai, Ax, &Symbolic, NULL, NULL) ;
	(void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL) ;
	umfpack_di_free_symbolic (&Symbolic) ;
	(void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, rh, Numeric, NULL, NULL) ;
	umfpack_di_free_numeric (&Numeric) ;

#ifdef __SOLVESYS_DEBUG
	std::cout << "finish solving with UMFPACK" << std::endl << std::endl;
#endif		

#ifdef __SOLVESYS_DEBUG_LV2
	for(int i = 0; i < dof; i++)
	{
		std::cout << " x[" << i << "] = " << x[i] << std::endl << std::endl;
	}
#endif
	return 0;
}
template <typename MyProblem>
int BasicSolvingSystem<MyProblem>::SuperLUSolve(Mesh<MyProblem> mesh)
{

#ifdef __SOLVESYS_DEBUG
	std::cout << "start solving with SuperLU" << std::endl;
#endif

	superlu_options_t options;
	set_default_options(&options);

	int nnz = Ap[dof]; 

	SuperMatrix A, B;
	dCreate_CompCol_Matrix(&A, dof, dof, nnz, Ax, Ai, Ap, SLU_NC, SLU_D, SLU_GE);
	dCreate_Dense_Matrix(&B, dof, 1, rh, dof, SLU_DN, SLU_D, SLU_GE);

	set_default_options(&options);  
    options.ColPerm = NATURAL; 

    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
	if ( !(perm_c = intMalloc(dof)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(dof)) ) ABORT("Malloc fails for perm_r[].");

    SuperMatrix L, U;
    int info;
    SuperLUStat_t stat;
    StatInit(&stat);

    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    
    double *sol = (double*) ((DNformat*) B.Store)->nzval;
	x = new double [dof];
    memcpy(x, sol, dof * sizeof(double));

    StatFree(&stat);

    SUPERLU_FREE (rh);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U); 

#ifdef __SOLVESYS_DEBUG
	std::cout << "finish solving with SuperLU" << std::endl << std::endl;
#endif	

	return 0;
}
template <typename MyProblem>
int BasicSolvingSystem<MyProblem>::addToMA(double a, int row, int col)
{
	if(a == 0) return 0;

	maColEle* pmaColEle;
	std::list<maColEle>::iterator it;
	for(it = ma[col].begin(); it != ma[col].end(); it++)
		if(it -> row == row)
			break;

	if(it != ma[col].end()){
		it -> value += a;
	}
	else{
		pmaColEle = new maColEle;
		pmaColEle -> row = row;
		pmaColEle -> value = a;
		ma[col].push_back(*pmaColEle);
	}

	return 0;
}
template <typename MyProblem>
int BasicSolvingSystem<MyProblem>::fileOutputTriplet(Mesh<MyProblem> mesh, MyProblem prob)
{
	std::ofstream fout((prob.parameters.meshFilename + ".triplet").c_str());

	std::vector< std::list<maColEle> >::iterator it;
	std::list<maColEle>::iterator it1;
	int k;
	for(it = ma.begin(), k = 1; it != ma.end(); it++, k++)
		for(it1 = it -> begin(); it1 != it -> end(); it1++)
			fout << ((it1 -> row) + 1) << " " << k << " " << (it1 -> value) << std::endl;

	return 0;
}
template <typename MyProblem>
int BasicSolvingSystem<MyProblem>::fileOutputRH(Mesh<MyProblem> mesh, MyProblem prob)
{
	std::ofstream fout((prob.parameters.meshFilename + ".rh").c_str());

	for(int i = 0; i < dof; i++)
		fout << rh[i] << std::endl;

	return 0;		
}
template <typename MyProblem>
int BasicSolvingSystem<MyProblem>::fileOutputMA(Mesh<MyProblem> mesh, MyProblem prob)
{
	std::ofstream fout((prob.parameters.meshFilename + ".ma").c_str());

	for(int i = 0; i < ma.size() + 1; i++)
		fout << Ap[i] << " ";
	fout << std::endl;

	int nnz(0);
	for(std::vector< std::list<maColEle> >::iterator it = ma.begin();
		it != ma.end(); it++)
		nnz += it -> size();

	for(int i = 0; i < nnz; i++)
		fout << Ai[i] << " ";
	fout << std::endl;

	for(int i = 0; i < nnz; i++)
		fout << Ax[i] << " ";
	fout << std::endl;

	return 0;
}

#endif

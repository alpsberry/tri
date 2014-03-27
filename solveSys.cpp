#define __SOLVESYS_DEBUG
// #define __SOLVESYS_DEBUG_LV2m

#include "solveSys.h"

bool maCompareFunc(const maColEle t1, const maColEle t2){  
	return t1.row < t2.row;  
} 

int BasicSolvingSystem::solveSparse(Mesh mesh, MyProblem prob)
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
int BasicSolvingSystem::convertToCSC(Mesh mesh)
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

int BasicSolvingSystem::UMFSolve(Mesh mesh)
{

#ifdef __SOLVESYS_DEBUG
	std::cout << "start solving with UMFPACK" << std::endl;
#endif

	x = new double [mesh.kidof];
	memset(x, 0, mesh.kidof * sizeof(double));
	void *Symbolic, *Numeric ;
	(void) umfpack_di_symbolic (mesh.kidof, mesh.kidof, Ap, Ai, Ax, &Symbolic, NULL, NULL) ;
	(void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL) ;
	umfpack_di_free_symbolic (&Symbolic) ;
	(void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, rh, Numeric, NULL, NULL) ;
	umfpack_di_free_numeric (&Numeric) ;

#ifdef __SOLVESYS_DEBUG
	std::cout << "finish solving with UMFPACK" << std::endl << std::endl;
#endif		

#ifdef __SOLVESYS_DEBUG_LV2
	for(int i = 0; i < mesh.kidof; i++)
	{
		std::cout << " x[" << i << "] = " << x[i] << std::endl << std::endl;
	}
#endif
	return 0;
}

int BasicSolvingSystem::SuperLUSolve(Mesh mesh)
{

#ifdef __SOLVESYS_DEBUG
	std::cout << "start solving with SuperLU" << std::endl;
#endif

	superlu_options_t options;
	set_default_options(&options);

	int nnz = Ap[mesh.kidof]; 

	SuperMatrix A, B;
	dCreate_CompCol_Matrix(&A, mesh.kidof, mesh.kidof, nnz, Ax, Ai, Ap, SLU_NC, SLU_D, SLU_GE);
	dCreate_Dense_Matrix(&B, mesh.kidof, 1, rh, mesh.kidof, SLU_DN, SLU_D, SLU_GE);

	set_default_options(&options);  
    options.ColPerm = NATURAL; 

    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
	if ( !(perm_c = intMalloc(mesh.kidof)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(mesh.kidof)) ) ABORT("Malloc fails for perm_r[].");

    SuperMatrix L, U;
    int info;
    SuperLUStat_t stat;
    StatInit(&stat);

    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    
    double *sol = (double*) ((DNformat*) B.Store)->nzval;
	x = new double [mesh.kidof];
    memcpy(x, sol, mesh.kidof * sizeof(double));

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

int BasicSolvingSystem::addToMA(double a, int row, int col)
{
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

// gradInteg = \int_E \nabla \phi_i \nabla v
// timeInteg = \int_E uv
// nodal basis are used here
double MySolvingSystem::integA(Element &ele, Mesh mesh, int vi, int vj, MyProblem prob)
{
	if(vi > vj)
		std::swap(vi, vj);

	double x1(mesh.vertex[ele.vertex[0] - 1].x), y1(mesh.vertex[ele.vertex[0] - 1].y),
		   x2(mesh.vertex[ele.vertex[1] - 1].x), y2(mesh.vertex[ele.vertex[1] - 1].y),
		   x3(mesh.vertex[ele.vertex[2] - 1].x), y3(mesh.vertex[ele.vertex[2] - 1].y);

	if(ele.detBE == 0){
		ele.detBE = fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
	}

	double gradInteg = 0;
	if(vi == 0){
		switch(vj){
			case 0:{gradInteg = (y2 - y3) * (y2 - y3) + (x3 - x2) * (x3 - x2); break;}
			case 1:{gradInteg = (y2 - y3) * (y3 - y1) + (x3 - x2) * (x1 - x3); break;}
			case 2:{gradInteg = (y2 - y3) * (y1 - y2) + (x3 - x2) * (x2 - x1); break;}
		}
	}
	else if(vi == 1){
		switch(vj){
			case 1:{gradInteg = (y3 - y1) * (y3 - y1) + (x1 - x3) * (x1 - x3); break;}
			case 2:{gradInteg = (y3 - y1) * (y1 - y2) + (x1 - x3) * (x2 - x1); break;}
		}
	}
	else{
		gradInteg = (y1 - y2) * (y1 - y2) + (x2 - x1) * (x2 - x1);
	}

	gradInteg /= ele.detBE * 2;

	double timeInteg = 0;
	if(vi == vj)
		timeInteg = ele.detBE / 12.0;
	else
		timeInteg = ele.detBE / 24.0;

#ifdef __SOLVESYS_DEBUG_LV2
	std::cout << "  gradInteg = " << gradInteg << " timeInteg = " << timeInteg << std::endl;
#endif

	return gradInteg + timeInteg;
}

// calc \int_E f\phi
// \int_E f\phi = |detB_E| / 6 * (\sum_{i=1}^3 f(m_i)\phi(m_i))
// where m_i is the midpoint of each edge
double MySolvingSystem::integARH(Element ele, Mesh mesh, int vi, MyProblem prob)
{
	double x1(mesh.vertex[ele.vertex[vi] - 1].x), y1(mesh.vertex[ele.vertex[vi] - 1].y);

	double a = 0;
	for(int i = 0; i < ele.vertex.size(); i++){
		if (i == vi)
			continue;
		a += prob.f((x1 + mesh.vertex[ele.vertex[i] - 1].x) / 2, (y1 + mesh.vertex[ele.vertex[i] - 1].y) / 2);
	}

	return a * ele.detBE / 12.0;
}

int MySolvingSystem::getStiff(Element& ele, Mesh& mesh, MyProblem prob){

#ifdef __SOLVESYS_DEBUG_LV2
	std::cout << " assemble element " << ele.index << std::endl;
#endif	

	for(int i = 0; i < ele.vertex.size(); i++){
		if(mesh.vertex[ele.vertex[i] - 1].bctype > 0)
			continue;
		for(int j = 0; j < ele.vertex.size(); j++)
		{
			if(mesh.vertex[ele.vertex[j] - 1].bctype > 0)
				continue;
			double valIntegA = integA(ele, mesh, i, j, prob);
			addToMA(valIntegA, mesh.vertex[ele.vertex[i] - 1].index - 1,
				mesh.vertex[ele.vertex[j] - 1].index - 1);

#ifdef __SOLVESYS_DEBUG_LV2
				std::cout << "  global " << ele.vertex[i] << " local " << i
						  << " and global " << ele.vertex[j] << " local " << j
						  << " integ = " << valIntegA << std::endl;
#endif

		}
	}
	for(int i = 0; i < ele.vertex.size(); i++){
		if(mesh.vertex[ele.vertex[i] - 1].bctype > 0)
			continue;
		double valIntegARH = integARH(ele, mesh, i, prob);
		rh[mesh.vertex[ele.vertex[i] - 1].index - 1] += valIntegARH;

#ifdef __SOLVESYS_DEBUG_LV2
		std::cout << "  global " << ele.vertex[i] << " local " << i
				  << " rh" << mesh.vertex[ele.vertex[i] - 1].index - 1
				  << "= " << valIntegARH << std::endl;
#endif
	}
	return 0;
}

int MySolvingSystem::assembleStiff(Mesh &mesh, MyProblem prob)
{

#ifdef __SOLVESYS_DEBUG
	std::cout << "start forming system, kidof = " << mesh.kidof << std::endl;
#endif
	clock_t t = clock();

	rh = new double [mesh.kidof];
	memset(rh, 0, (mesh.kidof) * sizeof(double));
	ma.resize(mesh.kidof);

	for(std::vector<Element>::iterator it = mesh.element.begin();
		it != mesh.element.end(); it++)
	{
		getStiff(*it, mesh, prob);
	}

	t = clock() - t;

#ifdef __SOLVESYS_DEBUG
	std::cout << "finish forming system, t = "
			  << (double) t / CLOCKS_PER_SEC << "s"
			  << std::endl << std::endl;
#endif	

#ifdef __SOLVESYS_DEBUG_LV2
	for(int i = 0; i < mesh.kidof; i++)
	{
		std::cout<<" ma[" << i << "]" << std::endl;
		for(std::list<maColEle>::iterator it = ma[i].begin();
			it != ma[i].end(); it++)
			std::cout<<"    row " << it -> row << " value " << it -> value << std::endl;
	}
#endif

	return 0;
}
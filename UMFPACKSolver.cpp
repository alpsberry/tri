//
//  UMFPACKSolver.cpp
//  tri
//
//  Created by GBB on 21/2/15.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#include "UMFPACKSolver.h"
#include <umfpack.h>

std::vector<double> UMFPACKSolver::solveSparse()
{
    std::cout << "start solving with UMFPACK" << std::endl;
    double* x = new double [dof];
    memset(x, 0, dof * sizeof(double));
    void *Symbolic, *Numeric ;
    (void) umfpack_di_symbolic (dof, dof, Ap, Ai, Ax, &Symbolic, NULL, NULL) ;
    (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL) ;
    umfpack_di_free_symbolic (&Symbolic) ;
    (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, rh, Numeric, NULL, NULL) ;
    umfpack_di_free_numeric (&Numeric) ;

    std::vector<double> v(x, x + dof);

    delete [] x;

    std::cout << "finish solving with UMFPACK\n" << std::endl;

    return v;
}
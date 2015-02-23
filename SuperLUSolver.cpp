//
//  SuperLUSolver.cpp
//  tri
//
//  Created by GBB on 21/2/15.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#include "SuperLUSolver.h"
#include "../SuperLU_4.3/SRC/slu_ddefs.h"

std::vector<double> SuperLUSolver::solveSparse()
{
    std::cout << "start solving with SuperLU" << std::endl;
    
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

    double *sol = (double *) ((DNformat *) B.Store)->nzval;
    std::vector<double> v(sol, sol + dof);

    StatFree(&stat);

    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    // Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

    std::cout << "finish solving with SuperLU\n" << std::endl;

    return v;
}
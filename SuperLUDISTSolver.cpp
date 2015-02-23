//
//  SuperLUDISTSolver.cpp
//  tri
//
//  Created by GBB on 21/2/15.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#include "SuperLUDISTSolver.h"

SuperLUDISTSolver::SuperLUDISTSolver(std::vector< std::list<maColEle> > &ma, int femDof,
                                     double *femRH, gridinfo_t *superlu_grid, int femm_loc, int femfst_row)
    : LinearSolver(ma, femDof, femRH), m_loc(femm_loc), fst_row(femfst_row)
{
    grid = superlu_grid;

    nnz_loc = 0;
    for (auto &mal : ma) {
        mal.sort([](const maColEle& t1, const maColEle& t2) { return t1.row < t2.row; });
        nnz_loc += mal.size();
    }

    nzval_loc = new double [nnz_loc];
    rowptr = new int [m_loc + 1];
    colind = new int [nnz_loc];

    rowptr[0] = 0;

    for( int r = 0, current_nnz = 0; r < m_loc; ++r)
    {
        for (int c = 0; c < ma.size(); ++c)
        {
            std::list<maColEle>& ma_col = ma[c];
            std::list<maColEle>::iterator it = find_if(ma_col.begin(), ma_col.end(), [r] (const maColEle& m) { return m.row == r; } );
            if (it != ma_col.end())
            {
                nzval_loc[current_nnz] = it -> value;
                colind[current_nnz++] = c;
            }
        }
        rowptr[r + 1] = current_nnz;
    }

}

std::vector<double> SuperLUDISTSolver::solveSparse()
{
    if (grid -> iam == 0)
        std::cout << "start solving with SuperLU_DIST" << std::endl;

    superlu_options_t options;
    SuperLUStat_t stat;
    SuperMatrix A;
    ScalePermstruct_t ScalePermstruct;
    LUstruct_t LUstruct;
    SOLVEstruct_t SOLVEstruct;
    double   *berr;
    int      info, nrhs;

    nrhs = 1;   /* Number of right-hand side. */

    dCreate_CompRowLoc_Matrix_dist(&A, dof, dof, nnz_loc, m_loc, fst_row,
                                   nzval_loc, colind, rowptr,
                                   SLU_NR_loc, SLU_D, SLU_GE);

    if ( !(berr = doubleMalloc_dist(nrhs)) )
        ABORT("Malloc fails for berr[].");
    /* ------------------------------------------------------------
       NOW WE SOLVE THE LINEAR SYSTEM.
       ------------------------------------------------------------*/
    set_default_options_dist(&options);

    int m = A.nrow;
    int n = A.ncol;

    /* Initialize ScalePermstruct and LUstruct. */
    ScalePermstructInit(m, n, &ScalePermstruct);
    LUstructInit(m, n, &LUstruct);

    /* Initialize the statistics variables. */
    PStatInit(&stat);

    /* Call the linear equation solver. */
    pdgssvx(&options, &A, &ScalePermstruct, rh, m_loc, nrhs, grid,
            &LUstruct, &SOLVEstruct, berr, &stat, &info);
    
    std::vector<double> v(rh, rh + m_loc); // save the solution

    PStatPrint(&options, &stat, grid);        /* Print the statistics. */

    /* ------------------------------------------------------------
       DEALLOCATE STORAGE.
       ------------------------------------------------------------*/

    PStatFree(&stat);
    // Destroy_CompRowLoc_Matrix_dist(&A);
    ScalePermstructFree(&ScalePermstruct);
    Destroy_LU(n, grid, &LUstruct);
    LUstructFree(&LUstruct);
    if ( options.SolveInitialized )
    {
        dSolveFinalize(&options, &SOLVEstruct);
    }
    SUPERLU_FREE(berr);


    if (grid -> iam == 0)
        std::cout << "finish solving with SuperLU_DIST\n" << std::endl;

    return v;
}

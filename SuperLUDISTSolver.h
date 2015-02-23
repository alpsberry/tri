//
//  SuperLUDISTSolver.h
//  tri
//
//  Created by GBB on 23/2/15.
//  Copyright (c) 2015 Xiaolin Guo. All rights reserved.
//

#ifndef __tri__SuperLUDISTSolver__
#define __tri__SuperLUDISTSolver__

#include "LinearSolver.h"
#include "../SuperLU_DIST_3.3/SRC/superlu_ddefs.h"

class SuperLUDISTSolver: public LinearSolver
{
    gridinfo_t *grid;
    int m_loc, fst_row, nnz_loc;
    double *nzval_loc;
    int *colind, *rowptr;
public:
    SuperLUDISTSolver(std::vector< std::list<maColEle> > &ma,
                      int femDof, double *femRH, gridinfo_t *superlu_grid, int m_loc, int fst_row);

    std::vector<double> solveSparse();  // call SuperLU to solve the sparse linear system
    ~SuperLUDISTSolver()
    {
        delete [] nzval_loc;
        delete [] colind;
        delete [] rowptr;
    }
};


#endif /* defined(__tri__SuperLUDISTSolver__) */

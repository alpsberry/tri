//
//  SuperLUSolver.h
//  tri
//
//  Created by GBB on 21/2/15.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#ifndef __tri__SuperLUSolver__
#define __tri__SuperLUSolver__

#include "LinearSolver.h"
#include "../SuperLU_4.3/SRC/slu_ddefs.h"

class SuperLUSolver: public LinearSolver
{
public:
    SuperLUSolver(std::vector< std::list<maColEle> > &ma,
                  int femDof, double *femRH)
        : LinearSolver(ma, femDof, femRH) {};

    std::vector<double> solveSparse();  // call SuperLU to solve the sparse linear system
};


#endif /* defined(__tri__SuperLUSolver__) */

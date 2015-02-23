//
//  UMFPACKSolver.h
//  tri
//
//  Created by GBB on 21/2/15.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#ifndef __tri__UMFPACKSolver__
#define __tri__UMFPACKSolver__

#include "LinearSolver.h"

class UMFPACKSolver: public LinearSolver
{
public:
    UMFPACKSolver(std::vector< std::list<maColEle> > &ma,
                  int femDof, double *femRH)
        : LinearSolver(ma, femDof, femRH) {};

    std::vector<double> solveSparse();  // call UMFPACK to solve the sparse linear system
};


#endif /* defined(__tri__UMFPACKSolver__) */

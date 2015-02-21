//
//  LinearSolver.h
//  tri
//
//  Created by GBB on 21/2/15.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#ifndef __tri__LinearSolver__
#define __tri__LinearSolver__

#include <vector>
#include <list>
#include <iostream>

// #include "../SuperLU_4.3/SRC/slu_ddefs.h"


struct maColEle
{
    int row;
    double value;
    maColEle(int r, double v): row(r), value(v) {};
};


class LinearSolver
{
protected:
    //for Compressed Sparse Column (CSC) format
    int *Ap;    //Ap[0] = 0; Ap[k] num of nonzero entries in the first k columns
    int *Ai;    //row of each nonzero entry, column-wise
    double *Ax; //value of each nonzero entry, column-wise

    int dof;    // degrees of freedom
    double *rh; // right-hand side vector

    LinearSolver(std::vector< std::list<maColEle> > &ma,
                 int femDof, double *femRH);

public:
    virtual std::vector<double> solveSparse() = 0; // solve the sparse linear system

    virtual ~LinearSolver()
    {
        delete [] Ap;
        delete [] Ai;
        delete [] Ax;
    }
};


#endif /* defined(__tri__LinearSolver__) */

//
//  LinearSolver.cpp
//  tri
//
//  Created by GBB on 21/2/15.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#include "LinearSolver.h"

// convert the list-stored matrix to CSC
LinearSolver::LinearSolver(std::vector< std::list<maColEle> > &ma,
                           int femDof, double *femRH):
    dof(femDof), rh(femRH)
{
    std::cout << "start converting to CSC structure" << std::endl;
    //sort entries in each column by their row
    std::vector< std::list<maColEle> >::iterator it;
    for (auto& mal:ma)
        mal.sort([](const maColEle t1, const maColEle t2) {return t1.row < t2.row;});

    Ap = new int [ma.size() + 1];
    Ap[0] = 0;
    int nnz(0), k(0);
    for ( it = ma.begin(), k = 1;
            it != ma.end(); it++, k++)
    {
        nnz += it -> size();
        Ap[k] = nnz;
    }

    Ai = new int [nnz];
    Ax = new double [nnz];
    for (it = ma.begin(), k = 0; it != ma.end(); it++)
        for (std::list<maColEle>::iterator it1 = it -> begin();
                it1 != it -> end(); it1++, k++)
        {
            Ai[k] = it1 -> row;
            Ax[k] = it1 -> value;
        }
    std::cout << "finish converting to CSC structure" << std::endl << std::endl;
}
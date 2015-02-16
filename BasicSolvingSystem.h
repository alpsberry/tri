//
//  BasicSolvingSystem.h
//  tri
//
//  Created by GBB on 26/12/14.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#ifndef __tri__BasicSolvingSystem__
#define __tri__BasicSolvingSystem__

#define __SOLVESYS_DEBUG

#include <ctime>
#include <vector>
#include <list>
#include "mesh.h"
// #include "../SuperLU_4.3/SRC/slu_ddefs.h"


struct maColEle {
    int row;
    double value;
    maColEle(int r, double v):row(r), value(v) {};
};

typedef std::vector< std::vector<double> > VECMATRIX;

bool maCompareFunc(const maColEle t1, const maColEle t2);

class BasicSolvingSystem {
    //for Compressed Sparse Column (CSC) format
    int *Ap;    //Ap[0] = 0; Ap[k] num of nonzero entries in the first k columns
    int *Ai;    //row of each nonzero entry, column-wise
    double *Ax; //value of each nonzero entry, column-wise

    int convertToCSC(); // convert the list-stored matrix to CSC
    int UMFSolve();     // call UMFPACK to solve the sparse linear system
    // int SuperLUSolve(Mesh mesh); // call SuperLU to solve the sparse linear system
    
    int fileOutputTriplet(); // file-output stiffness matrix in triplet format
    int fileOutputRH();      // file-output right-hand side vector
    int fileOutputMA();      // file-output stiffness matrix in CSC format
    
protected:
    Mesh* mesh;
    Problem* prob;
    
    int dof;    // degrees of freedom
    std::vector< std::list<maColEle> > ma; // list-stored stiffness matrix
    double *rh; // right-hand side vector
    double *x;  // the numerical solution

    BasicSolvingSystem(Mesh* m, Problem* p):mesh(m), prob(p) {}

    int addToMA(double a, int row, int col); // add value to list-stored stiffness matrix ma
    
public:
    int solveSparse(); // solve the sparse linear system
    virtual void output();
    virtual void assembleStiff() = 0;
    
    virtual ~BasicSolvingSystem() {};
};


#endif /* defined(__tri__BasicSolvingSystem__) */

//
//  BasicSolvingSystem.h
//  tri
//
//  Created by GBB on 26/12/14.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#ifndef __tri__BasicSolvingSystem__
#define __tri__BasicSolvingSystem__

#include <ctime>
#include <vector>
#include <list>
#include "mesh.h"
#include "LinearSolver.h"
#include "UMFPACKSolver.h"
#include "SuperLUSolver.h"
#include "SuperLUDISTSolver.h"

typedef std::vector< std::vector<double> > VECMATRIX;

class BasicSolvingSystem {
  
    int fileOutputTriplet(); // file-output stiffness matrix in triplet format
    int fileOutputRH();      // file-output right-hand side vector
    int fileOutputMA();      // file-output stiffness matrix in CSC format
    
protected:
    Mesh* mesh;
    Problem* prob;
    
    int dof;    // degrees of freedom
    double *rh; // right-hand side vector
    std::vector<double> x;  // the numerical solution

    std::vector< std::list<maColEle> > ma; // list-stored stiffness matrix

    BasicSolvingSystem(Mesh* m, Problem* p):mesh(m), prob(p) {}

    int addToMA(double a, int row, int col); // add value to list-stored stiffness matrix ma
    
public:
    virtual void solveSparse(); // solve the sparse linear system
    virtual void output();
    virtual void assembleStiff() = 0;
    
    virtual ~BasicSolvingSystem() {
        delete[] rh;
    };
};


#endif /* defined(__tri__BasicSolvingSystem__) */

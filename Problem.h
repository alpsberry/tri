//
//  Problem.h
//  tri
//
//  Created by GBB on 26/12/14.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#ifndef __tri__Problem__
#define __tri__Problem__

#define DEFAULT_PARAM_FILE "tri.input"
#define DEFAULT_SOLVE_PACK 0
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>


const int constNonrefined = -1;

enum class SolPack {
    UMFPACK, SuperLU, SuperLUDist, Count
};

struct paramstruct {
    std::string meshFilename; // mesh filename
    int nRefine;              // number of refinement times
    SolPack solPack;          // solving package, UMFPACK or SuperLU
    int cprintMeshInfo;       // print mesh info in console
    int cprintError;          // compute and print error
    int printResults;         // output results in console
    int fprintResults;        // file output results, *.output
    int fprintMA;             // file output stiff matrix in compressed column form, *.ma
    int fprintRH;             // file output righ-hand side matrix, *.rh
    int fprintTriplet;        // file output stiff matrix in triplet form, *.triplet
};

class Problem {
public:
    paramstruct parameters;
    int dimension;
    double epsilon;
    double sigma0;
    double beta0;
    
    Problem(int argc, char const *argv[]);
    
    virtual double f(double x, double y) = 0;
    
    // Dirichlet boundary condition
    virtual double gd(double x, double y)
    {
        return 0;
    }
    
    virtual double trueSol(double x, double y)
    {
        return 0;
    }
    
};

#endif /* defined(__tri__Problem__) */

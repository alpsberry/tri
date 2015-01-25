//
//  DGProblem.cpp
//  tri
//
//  Created by GBB on 26/12/14.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#include "DGProblem.h"
using std::cout;
using std::endl;
using std::getline;
using std::string;

// read parameters from an input file
DGProblem::DGProblem(int argc, char const *argv[]):Problem(argc, argv)
{
    cout << " epsilon = " << epsilon << endl
    << " sigma0 = " << sigma0 << endl
    << " beta0 = " << beta0 << endl
    << " refinement level = " << parameters.nRefine << endl;
}

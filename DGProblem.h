//
//  DGProblem.h
//  tri
//
//  Created by GBB on 26/12/14.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#ifndef __tri__DGProblem__
#define __tri__DGProblem__

#include <iostream>
#include "problem.h"

// the problem to be solved here is
// -\Delta u + u = 3\cos x \sin y, (x,y) \in \Omega = [\frac{\pi}{2}, \frac{3\pi}{2}] \times [0, \pi]
// u = 0, (x,y) \in \Gamma
class DGProblem: public Problem {
public:
    // read parameters from an input file
    DGProblem(int argc, char const *argv[]);
    
    // double f(double x, double y)
    // {
    //     return 2 * cos(x) * sin(y);
    //     // return 1.0;
    // }
    // double trueSol(double x, double y)
    // {
    //     return cos(x) * sin(y);
    // }
    
    double f(double x, double y)
    {
        return -1;
    }
    double gd(double x, double y)
    {
        return 0.25 * (x * x + y * y) + 2;
    }
    double trueSol(double x, double y)
    {
        return 0.25 * (x * x + y * y) + 2;
    }
};

#endif /* defined(__tri__DGProblem__) */

//
//  main.cpp
//  tri
//
//  Created by GBB on 26/12/14.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#include <iostream>
#include "mesh.h"
#include "BasicSolvingSystem.h"
#include "DGSolvingSystem.h"
#include "DGProblem.h"
#include "problem.h"

using namespace std;

int main(int argc, const char * argv[]) {
    try {
        DGProblem prob(argc, argv);
        Mesh mesh(&prob);
        
        BasicSolvingSystem* solSys = new DGSolvingSystem(&mesh, &prob);
        solSys -> assembleStiff();
        solSys -> solveSparse();
        solSys -> output();
        
    } catch (std::runtime_error &e) {
        cout << e.what() << endl;
    }
    
    return 0;
}
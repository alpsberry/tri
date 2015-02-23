//
//  maintrimpi.cpp
//  tri
//
//  Created by GBB on 22/2/15.
//  Copyright (c) 2015 Xiaolin Guo. All rights reserved.
//

#include <iostream>
#include "mesh.h"
#include "DGSolvingSystemMPI.h"
#include "DGProblem.h"
#include "problem.h"

using namespace std;

int main(int argc, const char *argv[])
{
    gridinfo_t grid;

    try
    {
        MPI_Init(NULL, NULL);
        int nprow = 2, npcol = 2;

        superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);

        DGProblem prob(argc, argv);
        Mesh mesh(&prob);

        BasicSolvingSystem *solSys = new DGSolvingSystemMPI(&mesh, &prob, &grid);
        solSys -> assembleStiff();
        solSys -> solveSparse();
        solSys -> output();

        superlu_gridexit(&grid);
        MPI_Finalize();

    }
    catch (std::runtime_error &e)
    {
        cout << e.what() << endl;
    }

    return 0;
}
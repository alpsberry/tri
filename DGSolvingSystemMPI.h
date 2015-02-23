//
//  DGSolvingSystemMPI.h
//  tri
//
//  Created by GBB on 22/2/15.
//  Copyright (c) 2015 Xiaolin Guo. All rights reserved.
//
//  An MPI version of DGSolvingSystem
//  Divide the stiff matrix by row blocks, consisting with
// the data structure SuperLU_DIST uses, thus no explicit
// communication needed

#ifndef __tri__DGSolvingSystemMPI__
#define __tri__DGSolvingSystemMPI__

#include "DGSolvingSystem.h"
#include "../SuperLU_DIST_3.3/SRC/superlu_ddefs.h"
#include <mpi.h>

class DGSolvingSystemMPI: public DGSolvingSystem
{
    int iam;     // number of this processor
    int fst_row; // the first row this processor is in charge of in the stiff matrix
    int m_loc;   // number of rows in charge
    gridinfo_t *grid; // SuperLU_DIST grid

    int assembleElementMPI(Element ele);
    int assembleEdgeMPI(Edge edge);
	void addMiiToMAMPI(VECMATRIX M, Element E1, Element E2);
    int edgeIntegMPI(Edge edge, VECMATRIX &M11, VECMATRIX &M12, VECMATRIX &M21, VECMATRIX &M22);

    void gatherSolutions(); // gather solution vector x from all processors to root processor
public:
    DGSolvingSystemMPI(Mesh *m, Problem *p, gridinfo_t *superlu_grid): DGSolvingSystem(m, p)
    {
        grid = superlu_grid;
        iam = grid -> iam;
    }
    void solveSparse();
    void assembleStiff();
    void output(); 
};



#endif /* defined(__tri__DGSolvingSystemMPI__) */
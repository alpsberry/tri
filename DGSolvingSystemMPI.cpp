//
//  DGSolvingSystemMPI.cpp
//  tri
//
//  Created by GBB on 22/2/15.
//  Copyright (c) 2015 Xiaolin Guo. All rights reserved.
//

#include "DGSolvingSystemMPI.h"

using std::vector;
using std::cout;
using std::endl;

void DGSolvingSystemMPI::solveSparse() // now only solve with SuperLUDist
{
    LinearSolver *ls = nullptr;
    if (prob -> parameters.solPack == SolPack::SuperLUDist)
        ls = new SuperLUDISTSolver(ma, dof, rh, grid, m_loc, fst_row);
    // else havent implemented yet, need to gather ma and rh in order to solve with non-distributed solver

    if (ls != nullptr)
    {
        x = ls -> solveSparse();
        delete ls;
    }
}

int DGSolvingSystemMPI::assembleElementMPI(Element ele)
{
    
    VECMATRIX vecElementInteg = elementInteg(ele);
    
    for (int row = 0; row != vecElementInteg.size(); ++row)
        for (int col = 0; col != vecElementInteg[row].size(); ++col) {
            this -> addToMA(vecElementInteg[row][col], ele.dofIndex + row - fst_row, ele.dofIndex + col);
        }
    
    vector<double> vecElementIntegRhs;
    vecElementIntegRhs = elementIntegRhs(ele);
    for (int i = 0; i != vecElementIntegRhs.size(); ++i) {
        this -> rh[ele.dofIndex + i - fst_row] += vecElementIntegRhs[i];
    }
    
    return 0;
}

void DGSolvingSystemMPI::addMiiToMAMPI(VECMATRIX M, Element E1, Element E2)
{
    for (int row = 0; row != M.size(); ++row)
        for (int col = 0; col != M[row].size(); ++col)
            this -> addToMA(M[row][col], E1.dofIndex + row - fst_row, E2.dofIndex + col);
}

int DGSolvingSystemMPI::assembleEdgeMPI(Edge edge)
{
    VECMATRIX M11, M12, M21, M22;
    Element &E1 = mesh -> element[edge.neighborElement[0] - 1];
    
    if (edge.neighborElement.size() == 2)
    {
        Element &E2 = mesh -> element[edge.neighborElement[1] - 1];

        if ((fst_row <= E1.dofIndex  && E1.dofIndex < fst_row + m_loc) || (fst_row <= E2.dofIndex  && E2.dofIndex < fst_row + m_loc))
        {
            edgeIntegMPI(edge, M11, M12, M21, M22);

            if (fst_row <= E1.dofIndex  && E1.dofIndex < fst_row + m_loc)
            {
                addMiiToMAMPI(M11, E1, E1);
                addMiiToMAMPI(M12, E1, E2);
            }
            if (fst_row <= E2.dofIndex  && E2.dofIndex < fst_row + m_loc)
            {
                addMiiToMAMPI(M21, E2, E1);
                addMiiToMAMPI(M22, E2, E2);
            }
        }

    }
    else if ((fst_row <= E1.dofIndex  && E1.dofIndex < fst_row + m_loc))
    {

        vector<double> rhs;
        edgeInteg(edge, M11, rhs);

        addMiiToMAMPI(M11, E1, E1);

        for (int i = 0; i < 3; i++)
            rh[E1.dofIndex + i - fst_row] += rhs[i];
    }
    
    
    return 0;
}

int DGSolvingSystemMPI::edgeIntegMPI(Edge edge, VECMATRIX &M11, VECMATRIX &M12, VECMATRIX &M21, VECMATRIX &M22)
{
    if (edge.neighborElement.size() != 2) {
        cout << "invalid call in edge integ, element size not equal to 2" << endl;
        return 1;
    }
    
    Element &E1 = mesh -> element[edge.neighborElement[0] - 1];
    Element &E2 = mesh -> element[edge.neighborElement[1] - 1];
    
    double E1_x1(mesh -> vertex[ E1.vertex[0] - 1].x),
    E1_x2(mesh -> vertex[ E1.vertex[1] - 1].x),
    E1_x3(mesh -> vertex[ E1.vertex[2] - 1].x),
    E1_y1(mesh -> vertex[ E1.vertex[0] - 1].y),
    E1_y2(mesh -> vertex[ E1.vertex[1] - 1].y),
    E1_y3(mesh -> vertex[ E1.vertex[2] - 1].y);
    double E2_x1(mesh -> vertex[ E2.vertex[0] - 1].x),
    E2_x2(mesh -> vertex[ E2.vertex[1] - 1].x),
    E2_x3(mesh -> vertex[ E2.vertex[2] - 1].x),
    E2_y1(mesh -> vertex[ E2.vertex[0] - 1].y),
    E2_y2(mesh -> vertex[ E2.vertex[1] - 1].y),
    E2_y3(mesh -> vertex[ E2.vertex[2] - 1].y);
    
    vector<double> ne; // actually n_e^* as in the report
    VECMATRIX f_E1, f_E2;
    calc_ne_and_f_on_edge(edge, ne, f_E1, f_E2); // normal vector, not unit, actualy |e| / 4 * n_e
    
    vector<double> grad_ne_E1; // \nabla \phi_i^{E_1} \cdot n_e^*
    grad_ne_E1.push_back( ((E1_y2 - E1_y3) * ne[0] + (E1_x3 - E1_x2) * ne[1]) / E1.detBE );
    grad_ne_E1.push_back( ((E1_y3 - E1_y1) * ne[0] + (E1_x1 - E1_x3) * ne[1]) / E1.detBE );
    grad_ne_E1.push_back( ((E1_y1 - E1_y2) * ne[0] + (E1_x2 - E1_x1) * ne[1]) / E1.detBE );
    
    vector<double> grad_ne_E2; // \nabla \phi_i^{E_2} \cdot n_e^*
    grad_ne_E2.push_back( ((E2_y2 - E2_y3) * ne[0] + (E2_x3 - E2_x2) * ne[1]) / E2.detBE );
    grad_ne_E2.push_back( ((E2_y3 - E2_y1) * ne[0] + (E2_x1 - E2_x3) * ne[1]) / E2.detBE );
    grad_ne_E2.push_back( ((E2_y1 - E2_y2) * ne[0] + (E2_x2 - E2_x1) * ne[1]) / E2.detBE );
        
    initM(M11, M12, M21, M22, E1.localDof, E2.localDof);
    
    const double eps = prob->epsilon;
    
    if (fst_row <= E1.dofIndex  && E1.dofIndex < fst_row + m_loc)
    {
        getMii(edge, M11, E1, E1, f_E1, f_E1, eps, ne, grad_ne_E1, grad_ne_E1, -1,  1,  1);
        getMii(edge, M12, E1, E2, f_E1, f_E2, eps, ne, grad_ne_E1, grad_ne_E2, -1, -1, -1);
    }
    if (fst_row <= E2.dofIndex  && E2.dofIndex < fst_row + m_loc)
    {
        getMii(edge, M21, E2, E1, f_E2, f_E1, eps, ne, grad_ne_E2, grad_ne_E1,  1,  1, -1);
        getMii(edge, M22, E2, E2, f_E2, f_E2, eps, ne, grad_ne_E2, grad_ne_E2,  1, -1,  1);
    }
    return 0;
}

void DGSolvingSystemMPI::assembleStiff()
{
#ifdef __DGSOLVESYS_DEBUG
    if (iam == 0)
        cout << "start forming system" << endl;
#endif
    
    clock_t t = clock();
    
    // get dof, m_loc, fst_row
    this -> dof = retrieve_dof_count_element_dofIndex(*mesh); // get total dof
    m_loc = (dof / LocalDimension / (grid->nprow * grid->npcol)) * LocalDimension;
    fst_row = iam * m_loc;
    if ((m_loc * grid->nprow * grid->npcol) != dof)
    {
        if (iam == (grid->nprow * grid->npcol - 1)) /* last proc. gets all*/
            m_loc = dof - m_loc * (grid->nprow * grid->npcol - 1);
    }
#ifdef __DGSOLVESYS_DEBUG
    if (iam == 0)
        cout << " dof = " << this -> dof << endl;
    // cout << "I am processor " << iam << " m_loc = " << m_loc << " fst_row = " << fst_row << endl;
#endif
    
    // initialize rh, ma
    this -> rh = new double [m_loc];
    memset(this -> rh, 0, m_loc * sizeof(double));
    this -> ma.resize(dof);
    
    mesh -> calcDetBE(); //calculate det(B_E) for each element
    
    // assemble element integral related items
    int k = 1;
    for (auto it = mesh -> element.begin();
         it != mesh -> element.end(); ++it, ++k) {
        if (it -> reftype != constNonrefined)
            continue;
        
        if (fst_row <= it -> dofIndex  && it -> dofIndex < fst_row + m_loc) {
            assembleElementMPI(*it);
        }
    }
        
    // assemble edge integral related items
    k = 1;
    for (auto it = mesh -> edge.begin(); it != mesh -> edge.end(); ++it, ++k) {
        if (it -> reftype != constNonrefined)
            continue;
        
        assembleEdgeMPI(*it);
    }
    
    
    MPI_Barrier( grid -> comm);
    t = clock() - t;

#ifdef __DGSOLVESYS_DEBUG
    if (iam == 0)
        cout << "finish forming system, t = " << (double) t / CLOCKS_PER_SEC << endl << endl;
#endif

}


void DGSolvingSystemMPI::output()
{
    gatherSolutions();
    (prob -> parameters).fprintMA = 0;
    (prob -> parameters).fprintRH = 0;
    (prob -> parameters).fprintTriplet = 0;
    if (iam == 0)
        DGSolvingSystem::output();
    
}

void DGSolvingSystemMPI::gatherSolutions()
{
    int world_size;
    MPI_Comm_size(grid->comm, &world_size);

    int m_loc_tmp = (dof / LocalDimension / (grid->nprow * grid->npcol)) * LocalDimension;
    vector<int> recvcounts(world_size, m_loc_tmp);
    if ((m_loc_tmp * grid->nprow * grid->npcol) != dof)
        recvcounts[world_size - 1] = dof - m_loc_tmp * (grid->nprow * grid->npcol - 1);

    vector<int> displs(world_size);
    displs[0] = 0;
    for(int i = 1; i < world_size; ++i)
        displs[i] = recvcounts[i - 1] + displs[i - 1];

    double* tmp_x = new double [dof];
    MPI_Gatherv(x.data(), x.size(), MPI_DOUBLE, tmp_x, recvcounts.data(), displs.data(),
        MPI_DOUBLE, 0, grid->comm);

    if (iam == 0)
    {    
        std::vector<double> v(tmp_x, tmp_x + dof);
        x = v;
    }

}


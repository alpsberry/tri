//
//  BasicSolvingSystem.cpp
//  tri
//
//  Created by GBB on 26/12/14.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#include "BasicSolvingSystem.h"
#include <umfpack.h>

using namespace std;

void BasicSolvingSystem::solveSparse()
{
    LinearSolver *ls = nullptr;
    if (prob -> parameters.solPack == SolPack::UMFPACK)
        ls = new UMFPACKSolver(ma, dof, rh);
    else if (prob -> parameters.solPack == SolPack::SuperLU)
        ls = new SuperLUSolver(ma, dof, rh);

    if (ls != nullptr)
    {
        x = ls -> solveSparse();
        delete ls;
    }
}

int BasicSolvingSystem::addToMA(double a, int row, int col)
{
    if (a == 0) return 0;

    maColEle *pmaColEle;
    std::list<maColEle>::iterator it;
    for (it = ma[col].begin(); it != ma[col].end(); it++)
        if (it -> row == row)
            break;

    if (it != ma[col].end())
    {
        it -> value += a;
    }
    else
    {
        pmaColEle = new maColEle(row, a);
        ma[col].push_back(*pmaColEle);
    }

    return 0;
}

int BasicSolvingSystem::fileOutputTriplet()
{
    std::ofstream fout((prob->parameters.meshFilename + ".triplet").c_str());

    std::vector< std::list<maColEle> >::iterator it;
    std::list<maColEle>::iterator it1;
    int k;
    for (it = ma.begin(), k = 1; it != ma.end(); it++, k++)
        for (it1 = it -> begin(); it1 != it -> end(); it1++)
            fout << ((it1 -> row) + 1) << " " << k << " " << (it1 -> value) << std::endl;

    return 0;
}
int BasicSolvingSystem::fileOutputRH()
{
    std::ofstream fout((prob->parameters.meshFilename + ".rh").c_str());

    for (int i = 0; i < dof; i++)
        fout << rh[i] << std::endl;

    return 0;
}
int BasicSolvingSystem::fileOutputMA()
{
    std::ofstream fout((prob->parameters.meshFilename + ".ma").c_str());

    for (int i = 0; i < ma.size() + 1; i++)
        fout << Ap[i] << " ";
    fout << std::endl;

    int nnz(0);
    for (std::vector< std::list<maColEle> >::iterator it = ma.begin();
            it != ma.end(); it++)
        nnz += it -> size();

    for (int i = 0; i < nnz; i++)
        fout << Ai[i] << " ";
    fout << std::endl;

    for (int i = 0; i < nnz; i++)
        fout << Ax[i] << " ";
    fout << std::endl;

    return 0;
}

void BasicSolvingSystem::output()
{
    if (prob->parameters.fprintMA)
        fileOutputMA();
    if (prob->parameters.fprintRH)
        fileOutputRH();
    if (prob->parameters.fprintTriplet)
        fileOutputTriplet();
}


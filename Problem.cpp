//
//  Problem.cpp
//  tri
//
//  Created by GBB on 26/12/14.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#include "Problem.h"

using std::cout;
using std::endl;
using std::getline;
using std::string;


// read parameters from an input file
Problem::Problem(int argc, char const *argv[])
{
    string paramFile;
    if ( argc < 2 ) {
        paramFile = DEFAULT_PARAM_FILE;
        cout << "read parameters from default input file " << paramFile << endl;
    } else {
        paramFile = argv[1];
        cout << "read parameters from " << paramFile << endl;
    }
    
    std::ifstream fin(paramFile.c_str());
    if (!fin)
        throw std::runtime_error("input file open error!");
    
    string tempStr;
    fin >> parameters.meshFilename;
    getline(fin, tempStr);
    
    fin >> dimension;
    getline(fin, tempStr);
    
    fin >> epsilon;
    getline(fin, tempStr);
    
    fin >> sigma0;
    getline(fin, tempStr);
    
    fin >> beta0;
    getline(fin, tempStr);
    
    fin >> parameters.nRefine;
    getline(fin, tempStr);
    
    int intSolPack;
    fin >> intSolPack;
    if (intSolPack < static_cast<int>(SolPack::Count))
        parameters.solPack = static_cast<SolPack>(intSolPack);
    else
        parameters.solPack = static_cast<SolPack>(DEFAULT_SOLVE_PACK);
    
    getline(fin, tempStr);
    
    fin >> parameters.cprintMeshInfo;
    getline(fin, tempStr);
    
    fin >> parameters.cprintError;
    getline(fin, tempStr);
    
    fin >> parameters.printResults;
    getline(fin, tempStr);
    
    fin >> parameters.fprintResults;
    getline(fin, tempStr);
    
    fin >> parameters.fprintMA;
    getline(fin, tempStr);
    
    fin >> parameters.fprintRH;
    getline(fin, tempStr);
    
    fin >> parameters.fprintTriplet;
    getline(fin, tempStr);
    
}

#include "DGProblem.h"

using std::cout;
using std::endl;

// read parameters from an input file
int DGProblem::initProblem(int argc, char const *argv[])
{
    std::string paramFile;
    if ( argc < 2 ) {
        paramFile = DEFAULT_PARAM_FILE;
        cout << "read parameters from default input file " << paramFile << endl;
    } else {
        paramFile = argv[1];
        cout << "read parameters from " << paramFile << endl;
    }

    std::ifstream fin(paramFile.c_str());

    if (!fin) {
        cout << "input file open error!" << endl;
        return 1;
    }

    std::string tempStr;
    fin >> parameters.meshFilename;
    std::getline(fin, tempStr);

    fin >> dimension;
    std::getline(fin, tempStr);

    fin >> epsilon;
    std::getline(fin, tempStr);

    fin >> sigma0;
    std::getline(fin, tempStr);

    fin >> beta0;
    std::getline(fin, tempStr);

    fin >> parameters.nRefine;
    std::getline(fin, tempStr);

    int intSolPack;
    fin >> intSolPack;
    if (intSolPack < static_cast<int>(SolPack::Count))
        parameters.solPack = static_cast<SolPack>(intSolPack);
    else
        parameters.solPack = static_cast<SolPack>(DEFAULT_SOLVE_PACK);

    std::getline(fin, tempStr);

    fin >> parameters.cprintMeshInfo;
    std::getline(fin, tempStr);

    fin >> parameters.cprintError;
    std::getline(fin, tempStr);

    fin >> parameters.printResults;
    std::getline(fin, tempStr);

    fin >> parameters.fprintResults;
    std::getline(fin, tempStr);

    fin >> parameters.fprintMA;
    std::getline(fin, tempStr);

    fin >> parameters.fprintRH;
    std::getline(fin, tempStr);

    fin >> parameters.fprintTriplet;
    std::getline(fin, tempStr);

    cout << " epsilon = " << epsilon << endl
         << " sigma0 = " << sigma0 << endl
         << " beta0 = " << beta0 << endl
         << " refinement level = " << parameters.nRefine << endl;

    return 0;
}

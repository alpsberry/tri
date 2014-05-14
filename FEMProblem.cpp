#include "FEMProblem.h"

using std::cout;
using std::endl;

// read parameters from an input file
int FEMProblem::initFEMProblem(int argc, char const *argv[])
{
	std::string paramFile;
	if( argc < 2 ) {
		paramFile = DEFAULT_PARAM_FILE;
		cout << "read parameters from default input file " << paramFile << endl;
	}
	else{
		paramFile = argv[1];
		cout << "read parameters from " << paramFile << endl;
	}

	std::ifstream fin(paramFile.c_str());

	if(!fin){
		cout << "input file open error!" << endl;
		return 1;
	}

	std::string tempStr;
	fin >> parameters.meshFilename;
	std::getline(fin, tempStr);
	
	fin >> dimension;
	std::getline(fin, tempStr);

	parameters.nRefine = 0;

	initProblem(fin);

	return 0;
}

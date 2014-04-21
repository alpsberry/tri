#include "DGProblem.h"

using std::cout;
using std::endl;

// read parameters from an input file
int DGProblem::initDGProblem(int argc, char const *argv[])
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

	fin >> epsilon;
	std::getline(fin, tempStr);

	fin >> sigma0;
	std::getline(fin, tempStr);

	fin >> beta0;
	std::getline(fin, tempStr);

	initProblem(fin);

	cout << " epsilon = " << epsilon <<endl
		 << " sigma0 = " << sigma0 << endl
		 << " beta0 = " << beta0 << endl;

	return 0;
}

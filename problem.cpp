#include "problem.h"

// read parameters from an input file
int Problem::initProblem(int argc, char const *argv[])
{
	std::string paramFile;
	if( argc < 2 ) {
		paramFile = DEFAULT_PARAM_FILE;
		std::cout << "read parameters from default input file " << paramFile << std::endl;
	}
	else{
		paramFile = argv[1];
		std::cout << "read parameters from " << paramFile << std::endl;
	}

	std::ifstream fin(paramFile.c_str());

	if(!fin){
		std::cout << "input file open error!" << std::endl;
		return 1;
	}

	std::string tempStr;
	fin >> parameters.meshFilename;
	std::getline(fin, tempStr);
	
	fin >> dimension;
	std::getline(fin, tempStr);

	int intSolPack;
	fin >> intSolPack;
	if(intSolPack < static_cast<int>(SolPack::Count))
	// intSolPack = static_cast<std::underlying_type<SolPack>::type>(SolPack::SuperLU);
		parameters.solPack = static_cast<SolPack>(intSolPack);
	else
		parameters.solPack = static_cast<SolPack>(DEFAULT_SOLVE_PACK);
	// switch(intSolPack){
	// 	case 0:{parameters.solPack = UMFPACK; break;}
	// 	case 1:{parameters.solPack = SuperLU; break;}
	// 	default:{parameters.solPack = DEFAULT_SOLVE_PACK; break;}
	// }
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


	// std::cout << parameters.meshFilename << std::endl << dimension << std::endl;

	return 0;
}

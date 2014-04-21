#ifndef TRI_PROBLEM_H
#define TRI_PROBLEM_H
#define DEFAULT_PARAM_FILE "tri.input"
#define DEFAULT_SOLVE_PACK 0
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

struct Vertex{
	int index;
	double x, y;
	int bctype;
};

struct Element;

struct Edge{
	int index;
	std::vector<int> vertex;
	std::vector<int> neighborElement;
	int bctype;
};

struct Element{
	double detBE;
	int index;
	int localDof;
	std::vector<int> vertex;
	std::vector<int> edge;
};
// detBE is the determinant of B_E
// The mapping F_E from the reference element to any element E is
// F([x, y]^T) = B_E[x, y]^T + b_E
// where B_E = [x2 - x1, x3 - x1; y2 - y1. y3 - y1] and b_E = [x1, y1]^T

class Region{
public:
	std::vector<Element> element;
	std::vector<Edge> edge;
	std::vector<Vertex> vertex;
};

enum class SolPack {UMFPACK, SuperLU, Count};

struct paramstruct
{
	std::string meshFilename; // mesh filename
	SolPack solPack;          // solving package, UMFPACK or SuperLU 
	int cprintError;          // compute and print error
	int printResults;	      // output results in console
	int fprintResults;	      // file output results, *.output
	int fprintMA;	          // file output stiff matrix in compressed column form, *.ma
	int fprintRH;	          // file output righ-hand side matrix, *.rh
	int fprintTriplet;	      // file output stiff matrix in triplet form, *.triplet
};

class Problem
{
public:
	// Region region;
	paramstruct parameters;
	int dimension;

	virtual double f(double x, double y){
		return 0;
	}

	// Dirichlet boundary condition
	virtual double gd(double x, double y){
		return 0;
	}

	int initProblem(std::ifstream &fin)
	{
		std::string tempStr;
		
		int intSolPack;
		fin >> intSolPack;
		if(intSolPack < static_cast<int>(SolPack::Count))
			parameters.solPack = static_cast<SolPack>(intSolPack);
		else
			parameters.solPack = static_cast<SolPack>(DEFAULT_SOLVE_PACK);

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
	
		return 0;

	}


};

#endif

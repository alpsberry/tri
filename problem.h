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
	int dofIndex; // index in global dof
	double x, y;  // coordinates
	int bctype;   // boundary condition type, 0 for interior vertex, 1 for dirichlet boundary, 2 for neumann boundary
};

struct Element;

struct Edge{
	int index;
	int reftype;                       // refinement type, -1 for not refined, 0, 1, 2, ... for refinement level
	int bctype;                        // boundary condition type, 0 for interior vertex, 1 for dirichlet boundary, 2 for neumann boundary
	std::vector<int> vertex;           // indices of its vertices
	std::vector<int> neighborElement;  // indices of elements sharing the edge
};

struct Element{
	int index;
	double detBE;             // see below
	int localDof;             // local degrees of freedom
	int dofIndex;             // index in global dof for the first local dof
	int reftype;              // refinement type, -1 for not refined, 0, 1, 2, ... for refinement level
	std::vector<int> vertex;  // indices of its vertices
	std::vector<int> edge;    // indices of its edges
	int parent;  // indices of parent element for refinement
	std::vector<int> child;   // indices of child elements for refinement
};
// detBE is the determinant of B_E
// The mapping F_E from the reference element to any element E is
// F([x, y]^T) = B_E[x, y]^T + b_E
// where B_E = [x2 - x1, x3 - x1; y2 - y1. y3 - y1] and b_E = [x1, y1]^T

const int constNonrefined = -1;

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
	int nRefine;              // number of refinement times
	SolPack solPack;          // solving package, UMFPACK or SuperLU
	int cprintMeshInfo;       // print mesh info in console 
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
	double epsilon;
	double sigma0;
	double beta0;
	virtual double f(double x, double y){
		return 0;
	}

	// Dirichlet boundary condition
	virtual double gd(double x, double y){
		return 0;
	}

	virtual double trueSol(double x, double y){ return 0; }

	virtual int initProblem(int argc, char const *argv[]) = 0;

};

#endif

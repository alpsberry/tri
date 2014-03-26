//mesh.h
#ifndef TRI_MESH_H
#define TRI_MESH_H

#include <iostream>
#include <fstream>
#include <cstring>
#include "problem.h"

class Mesh : public Region
{
public:
	int kidof;
	int kbdof;
	Mesh(){
		kidof = 0;
		kbdof = 0;
	}
	int initElement(MyProblem prob, std::string filename);

	int initEdge(MyProblem prob, std::string filename);

	int initVertex(MyProblem prob, std::string filename);

	// void findElementEdge();

	int initMesh(MyProblem prob);
};

#endif /* TRI_MESH_H */
//mesh.h
#ifndef TRI_MESH_H
#define TRI_MESH_H

#define __MESH_DEBUG

#include <iostream>
#include <fstream>
#include <cstring>
#include "problem.h"
#include "DGProblem.h"

template < typename MyProblem > class Mesh : public Region
{
public:
	int kidof;
	int kbdof;
	Mesh(){
		kidof = 0;
		kbdof = 0;
	}
	int initElement(MyProblem prob);

	int initEdge(MyProblem prob);

	int initVertex(MyProblem prob);

	void findElementEdge();

	int initMesh(MyProblem prob);
};

template < typename MyProblem >
int Mesh<MyProblem>::initElement(MyProblem prob)
{
	std::ifstream fin((prob.parameters.meshFilename + ".ele").c_str());
		if(!fin){
		std::cout << "error opening " << prob.parameters.meshFilename + ".ele" << std::endl;
		return 1;
	}

	int numEle, numNode, numAttr;
	fin >> numEle >> numNode >> numAttr;

	if( numNode != prob.dimension + 1){
		std::cerr << "dimension does not match in "
			      << prob.parameters.meshFilename << std::endl;
		return 1;
	}

	element.resize(numEle);
	int tempVertex;
	for(std::vector<Element>::iterator it = element.begin();
		it != element.end(); it++)
	{
		it -> vertex.resize(prob.dimension + 1);
		fin >> it -> index;
		for(int j = 0; j < prob.dimension + 1; j++){
			fin >> tempVertex;
			it -> vertex[j] = tempVertex;
		}
		it -> detBE = 0;
	}
	return 0;
}

template < typename MyProblem >
int Mesh<MyProblem>::initEdge(MyProblem prob)
{
	std::ifstream fin((prob.parameters.meshFilename + ".edge").c_str());
	if(!fin){
		std::cout << "error opening " << prob.parameters.meshFilename + ".edge" << std::endl;
		return 1;
	}

	int numEdge, numBound;
	fin >> numEdge >> numBound;

	edge.resize(numEdge);
	for(std::vector<Edge>::iterator it = edge.begin(); it != edge.end(); it++){
		it -> vertex.resize(prob.dimension);
		fin >> it -> index;
		for(int j = 0; j < it -> vertex.size(); j++)
			fin >> it -> vertex[j];
		fin >> it -> bctype;
// it -> bctype = 0;
	}

	return 0;
}

template < typename MyProblem >
int Mesh<MyProblem>::initVertex(MyProblem prob)
{
	int numVer, dim, numAttr, numBound;
	std::ifstream fin((prob.parameters.meshFilename + ".node").c_str());
	if(!fin){
		std::cout << "error opening " << prob.parameters.meshFilename + ".node" << std::endl;
		return 1;
	}

	fin >> numVer >> dim >> numAttr >> numBound;

	if( dim != prob.dimension){
		std::cerr << "dimension does not match in "
			      << prob.parameters.meshFilename << std::endl;
		return 1;
	}

	vertex.resize(numVer);
	int tempBound;
	if(numBound > 0)
		for(std::vector<Vertex>::iterator it = vertex.begin();
			it != vertex.end(); it++){
			fin >> it -> index >> it -> x >> it -> y;
			fin >> tempBound;
			it -> bctype = tempBound;
// it -> bctype = 0;
			if(tempBound == 0){
				kidof++;
				it -> index = kidof;
			}
			else{
				kbdof++;
				it -> index = kbdof;
			}
		}

	return 0;
}

template < typename MyProblem >
void Mesh<MyProblem>::findElementEdge()
{
	for(std::vector<Edge>::iterator itEdge = edge.begin();
		itEdge != edge.end(); itEdge++)
		for(std::vector<Element>::iterator itEle = element.begin();
			itEle != element.end(); itEle++)
		{
			if(std::find(itEle -> vertex.begin(), itEle -> vertex.end(), itEdge -> vertex[0]) == itEle -> vertex.end())
				continue;
			if(std::find(itEle -> vertex.begin(), itEle -> vertex.end(), itEdge -> vertex[1]) == itEle -> vertex.end())
				continue;
			itEle -> edge.push_back(itEdge -> index);
			itEdge -> neighborElement.push_back(itEle -> index);
		}
}

template < typename MyProblem >
int Mesh<MyProblem>::initMesh(MyProblem prob)
{
#ifdef __MESH_DEBUG
	std::cout << "start initializing mesh" << std::endl;
#endif

	if(initEdge(prob) == 1)
		return 1;
#ifdef __MESH_DEBUG
	std::cout << " edge initialized" << std::endl;
#endif

	if(initElement(prob) == 1)
		return 1;
#ifdef __MESH_DEBUG
	std::cout << " element initialized" << std::endl;
#endif

	if(initVertex(prob) == 1)
		return 1;
#ifdef __MESH_DEBUG
	std::cout << " vertex initialized" << std::endl;
#endif
	
	findElementEdge();
#ifdef __MESH_DEBUG
	std::cout << " edge neighbor initialized" << std::endl;
#endif

#ifdef __MESH_DEBUG
	std::cout << "finish initializing mesh" << std::endl << std::endl;
#endif

	return 0;
}

#endif /* TRI_MESH_H */
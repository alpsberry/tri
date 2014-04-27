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

	int readRefinement(MyProblem prob);

	void findElementEdge(int previousLevelElementSize);

	int initMesh(MyProblem prob);

	void printVertex();
	void printEdge();
	void printElement();
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
		it -> reftype = -1;
		it -> localDof = 0;
		it -> detBE = 0;
		it -> parent = 0;
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
		it -> reftype = -1;
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
int Mesh<MyProblem>::readRefinement(MyProblem prob)
{
	for(int refineLevel = 0; refineLevel < prob.parameters.nRefine; refineLevel++)
	{
		std::string refFile = prob.parameters.meshFilename + ".ref" + std::to_string(refineLevel);
		std::ifstream fin(refFile.c_str());
		if(!fin){
			std::cout << "error opening " << refFile << std::endl;
			return 1;
		}

		// read refinement info on vertex
		int numNewVer;
		fin >> numNewVer;

		Vertex* pVer;
		for(int j = 0; j < numNewVer; j++){
			pVer = new Vertex;
			fin >> pVer -> index >> pVer -> x >> pVer -> y >> pVer -> bctype;
			vertex.push_back(*pVer);
		}

		// read refinement info on edge
		int numNewEdge, tempEdgeIndex, tempVertex, parentEdge, j(0);
		int sizeEdge = edge.size();
		fin >> numNewEdge;

		Edge* pEdge;
		while(j < numNewEdge){
			fin >> tempEdgeIndex;
			if(tempEdgeIndex <= sizeEdge){
				edge[tempEdgeIndex - 1].reftype = refineLevel;
				parentEdge = tempEdgeIndex;
			}
			else{
				pEdge = new Edge;
				pEdge -> index = tempEdgeIndex;
				for(int k = 0; k < prob.dimension; k++){
					fin >> tempVertex;
					pEdge -> vertex.push_back(tempVertex);
				}
				// fin >> pEdge -> bctype;
				pEdge -> reftype = -1;
				if(parentEdge > 0){
					pEdge -> neighborElement = edge[parentEdge - 1].neighborElement;
					pEdge -> bctype = edge[parentEdge - 1].bctype;
				}
				else
					pEdge -> bctype = 0;

				edge.push_back(*pEdge);

				++j;
			}
		}

		// read refinement info on element
		int numNewEle, tempEleIndex, parentIndex;
		int sizeEle = element.size();
		Element* pEle;
		fin >> numNewEle;

		j = 0;
		while(j < numNewEle)
		{
			fin >> tempEleIndex;
			if(tempEleIndex <= sizeEle){
				element[tempEleIndex - 1].reftype = refineLevel;
				parentIndex = tempEleIndex;
			}
			else{
				pEle = new Element;
				pEle -> index = tempEleIndex;
				for(int k = 0; k < prob.dimension + 1; k++){
					fin >> tempVertex;
					pEle -> vertex.push_back(tempVertex);
				}
				pEle -> reftype = -1;
				pEle -> localDof = 0;
				pEle -> detBE = 0;
				pEle -> parent = parentIndex;
				element.push_back(*pEle);
				element[parentIndex - 1].child.push_back(tempEleIndex);
				++j;
			}
		}

		findElementEdge(sizeEle);
	}

	return 0;
}

template < typename MyProblem >
void Mesh<MyProblem>::findElementEdge(int previousLevelElementSize)
{
	for(int i = previousLevelElementSize; i < element.size(); i++){
		Element& ele = element[i];
		for(Edge &ed : edge)
		{
			if(std::find(ele.vertex.begin(), ele.vertex.end(), ed.vertex[0]) == ele.vertex.end())
				continue;
			if(std::find(ele.vertex.begin(), ele.vertex.end(), ed.vertex[1]) == ele.vertex.end())
				continue;
			ele.edge.push_back(ed.index);

			if( (ed.bctype == 0 && ed.neighborElement.size() < 2) || (ed.bctype != 0 && ed.neighborElement.size() == 0))
				ed.neighborElement.push_back(ele.index);
			else{
				if(ed.bctype != 0){
					ed.neighborElement.clear();
					ed.neighborElement.push_back(ele.index);
				}
				else{
					std::vector<int>::iterator it = std::find(ed.neighborElement.begin(), ed.neighborElement.end(), ele.parent);
					if(it != ed.neighborElement.end())
						ed.neighborElement.erase(it);
					ed.neighborElement.push_back(ele.index);
				}
			}
		}
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
	findElementEdge(0);

	if(prob.parameters.nRefine > 0)
	{
		if(readRefinement(prob) == 1)
			return 1;
	#ifdef __MESH_DEBUG
		std::cout << " refinement initialized" << std::endl;
	#endif
	}

#ifdef __MESH_DEBUG
	std::cout << "finish initializing mesh" << std::endl << std::endl;
#endif

	if(prob.parameters.cprintMeshInfo){
		printVertex();
		printEdge();
		printElement();
	}

	return 0;
}

template < typename MyProblem >
void Mesh<MyProblem>::printVertex()
{
	std::cout << "vertex:" << std::endl;
	std::cout << "index\t" << "x\t" << "y\t" << "bctype" << std::endl;
	for(Vertex ver:vertex)
	{
		std::cout << ver.index << "\t" << ver.x << "\t" << ver.y << "\t" << ver.bctype << std::endl;
	}
	std::cout << std::endl;
	
}

template < typename MyProblem >
void Mesh<MyProblem>::printEdge()
{
	std::cout << "edge:" << std::endl;
	std::cout << "index\t" << "reftype\t" << "bctype\t" << "vertex\t" << "neighborElement" << std::endl;
	for(Edge ed:edge)
	{
		if(ed.neighborElement.size() == 1)
			std::cout << ed.index << "\t" << ed.reftype << "\t" << ed.bctype << "\t"
		        	  << ed.vertex[0] << "  " << ed.vertex[1] << "\t" // << std::endl;
		        	  << ed.neighborElement[0] << std::endl;
		else
			if(ed.neighborElement.size() == 2)
				std::cout << ed.index << "\t" << ed.reftype << "\t" << ed.bctype << "\t"
			        	  << ed.vertex[0] << "  " << ed.vertex[1] << "\t" // << std::endl;
			        	  << ed.neighborElement[0] << "  " << ed.neighborElement[1] << std::endl;
			else
				std::cout << ed.index << "\t" << ed.reftype << "\t" << ed.bctype << "\t"
			        	  << ed.vertex[0] << "  " << ed.vertex[1] << "\t" << std::endl;
	}
	std::cout << std::endl;
	
}

template < typename MyProblem >
void Mesh<MyProblem>::printElement()
{
	std::cout << "element:" << std::endl;
	std::cout << "index\t" << "reftype\t" << "vertex\t\t" << "edge\t\t" << "parent\t" << "child" << std::endl;
	for(Element ele:element)
	{
		std::cout << ele.index << "\t" << ele.reftype << "\t"
		          << ele.vertex[0] << "   " << ele.vertex[1] << "   " << ele.vertex[2] << "\t"
		          << ele.edge[0] << "   " << ele.edge[1] << "   " << ele.edge[2] << "\t"
		          << ele.parent << "\t";
		for(int jChild : ele.child)
			std::cout << jChild << " ";
		std::cout << std::endl;
	}

	std::cout << std::endl;

}


#endif /* TRI_MESH_H */
//
//  Mesh.cpp
//  tri
//
//  Created by GBB on 26/12/14.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#include "Mesh.h"

using std::vector;

Mesh::Mesh(Problem* prob)
{
    try {
#ifdef __MESH_DEBUG
        std::cout << "start initializing mesh" << std::endl;
#endif
        _meshFilename = prob -> parameters.meshFilename;
        _dimension = prob -> dimension;
        
        
        initEdge();
#ifdef __MESH_DEBUG
        std::cout << " edge initialized" << std::endl;
#endif
        
        initElement();
#ifdef __MESH_DEBUG
        std::cout << " element initialized" << std::endl;
#endif
        
        initVertex();
#ifdef __MESH_DEBUG
        std::cout << " vertex initialized" << std::endl;
#endif
        findElementEdge(0);
        
        if (prob -> parameters.nRefine > 0) {
            readRefinement(prob -> parameters.nRefine);
#ifdef __MESH_DEBUG
            std::cout << " refinement initialized" << std::endl;
#endif
        }
        
#ifdef __MESH_DEBUG
        std::cout << "finish initializing mesh" << std::endl << std::endl;
#endif
        
        if (prob -> parameters.cprintMeshInfo) {
            printVertex();
            printEdge();
            printElement();
        }
    } catch (std::runtime_error &e) {
        throw e;
    }
    
}

int Mesh::initElement()
{
    std::ifstream fin((_meshFilename + ".ele").c_str());
    if (!fin)
        throw std::runtime_error("error opening " + _meshFilename + ".ele");
    
    int numEle, numNode, numAttr;
    fin >> numEle >> numNode >> numAttr;
    
    if ( numNode != _dimension + 1)
        throw std::runtime_error("dimension does not match in "
                                 + _meshFilename + ".ele");
    
    element.resize(numEle);
    int tempVertex;
    for (std::vector<Element>::iterator it = element.begin();
         it != element.end(); it++) {
        it -> vertex.resize(_dimension + 1);
        fin >> it -> index;
        for (int j = 0; j < _dimension + 1; j++) {
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

int Mesh::initEdge()
{
    std::ifstream fin((_meshFilename + ".edge").c_str());
    if (!fin)
        throw std::runtime_error("error opening " + _meshFilename + ".edge");
    
    int numEdge, numBound;
    fin >> numEdge >> numBound;
    
    edge.resize(numEdge);
    for (std::vector<Edge>::iterator it = edge.begin(); it != edge.end(); it++) {
        it -> vertex.resize(_dimension);
        fin >> it -> index;
        for (int j = 0; j < it -> vertex.size(); j++)
            fin >> it -> vertex[j];
        fin >> it -> bctype;
        it -> reftype = -1;
    }
    
    return 0;
}

int Mesh::initVertex()
{
    int numVer, dim, numAttr, numBound;
    std::ifstream fin((_meshFilename + ".node").c_str());
    if (!fin)
        throw std::runtime_error("error opening " + _meshFilename + ".node");
    
    fin >> numVer >> dim >> numAttr >> numBound;
    
    if ( dim != _dimension)
        throw std::runtime_error("dimension does not match in "
                                 + _meshFilename + ".node");
    
    vertex.resize(numVer);
    int tempBound;
    if (numBound > 0)
        for (std::vector<Vertex>::iterator it = vertex.begin();
             it != vertex.end(); it++) {
            fin >> it -> index >> it -> x >> it -> y;
            fin >> tempBound;
            it -> bctype = tempBound;
        }
    
    return 0;
}

int Mesh::readRefinement(int nRefine)
{
    for (int refineLevel = 0; refineLevel < nRefine; refineLevel++) {
        std::string refFile = _meshFilename + ".ref" + std::to_string(refineLevel);
        std::ifstream fin(refFile.c_str());
        if (!fin)
            throw std::runtime_error("error opening " + refFile);
        
        // read refinement info on vertex
        int numNewVer;
        fin >> numNewVer;
        
        Vertex *pVer;
        for (int j = 0; j < numNewVer; j++) {
            pVer = new Vertex;
            fin >> pVer -> index >> pVer -> x >> pVer -> y >> pVer -> bctype;
            vertex.push_back(*pVer);
        }
        
        // read refinement info on edge
        int numNewEdge, tempEdgeIndex, tempVertex, parentEdge, j(0);
        vector<Edge>::size_type sizeEdge = edge.size();
        fin >> numNewEdge;
        
        Edge *pEdge;
        while (j < numNewEdge) {
            fin >> tempEdgeIndex;
            if (tempEdgeIndex <= sizeEdge) {
                edge[tempEdgeIndex - 1].reftype = refineLevel;
                parentEdge = tempEdgeIndex;
            } else {
                pEdge = new Edge;
                pEdge -> index = tempEdgeIndex;
                for (int k = 0; k < _dimension; k++) {
                    fin >> tempVertex;
                    pEdge -> vertex.push_back(tempVertex);
                }
                
                pEdge -> reftype = -1;
                if (parentEdge > 0) {
                    pEdge -> neighborElement = edge[parentEdge - 1].neighborElement;
                    pEdge -> bctype = edge[parentEdge - 1].bctype;
                } else
                    pEdge -> bctype = 0;
                
                edge.push_back(*pEdge);
                
                ++j;
            }
        }
        
        // read refinement info on element
        int numNewEle, tempEleIndex, parentIndex;
        vector<Edge>::size_type sizeEle = element.size();
        Element *pEle;
        fin >> numNewEle;
        
        j = 0;
        while (j < numNewEle) {
            fin >> tempEleIndex;
            if (tempEleIndex <= sizeEle) {
                element[tempEleIndex - 1].reftype = refineLevel;
                parentIndex = tempEleIndex;
            } else {
                pEle = new Element;
                pEle -> index = tempEleIndex;
                for (int k = 0; k < _dimension + 1; k++) {
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

void Mesh::findElementEdge(vector<Edge>::size_type previousLevelElementSize)
{
    for (auto i = previousLevelElementSize; i < element.size(); i++) {
        Element &ele = element[i];
        for (Edge &ed : edge) {
            if (std::find(ele.vertex.begin(), ele.vertex.end(), ed.vertex[0]) == ele.vertex.end())
                continue;
            if (std::find(ele.vertex.begin(), ele.vertex.end(), ed.vertex[1]) == ele.vertex.end())
                continue;
            ele.edge.push_back(ed.index);
            
            if ( (ed.bctype == 0 && ed.neighborElement.size() < 2) || (ed.bctype != 0 && ed.neighborElement.size() == 0))
                ed.neighborElement.push_back(ele.index);
            else {
                if (ed.bctype != 0) {
                    ed.neighborElement.clear();
                    ed.neighborElement.push_back(ele.index);
                } else {
                    std::vector<int>::iterator it = std::find(ed.neighborElement.begin(), ed.neighborElement.end(), ele.parent);
                    if (it != ed.neighborElement.end())
                        ed.neighborElement.erase(it);
                    ed.neighborElement.push_back(ele.index);
                }
            }
        }
    }
}

void Mesh::printVertex()
{
    std::cout << "vertex:" << std::endl;
    std::cout << "index\t" << "x\t" << "y\t" << "bctype" << std::endl;
    for (Vertex ver : vertex) {
        std::cout << ver.index << "\t" << ver.x << "\t" << ver.y << "\t" << ver.bctype << std::endl;
    }
    std::cout << std::endl;
    
}

void Mesh::printEdge()
{
    std::cout << "edge:" << std::endl;
    std::cout << "index\t" << "reftype\t" << "bctype\t" << "vertex\t" << "neighborElement" << std::endl;
    for (Edge ed : edge) {
        if (ed.neighborElement.size() == 1)
            std::cout << ed.index << "\t" << ed.reftype << "\t" << ed.bctype << "\t"
            << ed.vertex[0] << "  " << ed.vertex[1] << "\t" // << std::endl;
            << ed.neighborElement[0] << std::endl;
        else if (ed.neighborElement.size() == 2)
            std::cout << ed.index << "\t" << ed.reftype << "\t" << ed.bctype << "\t"
            << ed.vertex[0] << "  " << ed.vertex[1] << "\t" // << std::endl;
            << ed.neighborElement[0] << "  " << ed.neighborElement[1] << std::endl;
        else
            std::cout << ed.index << "\t" << ed.reftype << "\t" << ed.bctype << "\t"
            << ed.vertex[0] << "  " << ed.vertex[1] << "\t" << std::endl;
    }
    std::cout << std::endl;
    
}

void Mesh::printElement()
{
    std::cout << "element:" << std::endl;
    std::cout << "index\t" << "reftype\t" << "vertex\t\t" << "edge\t\t" << "parent\t" << "child" << std::endl;
    for (Element ele : element) {
        std::cout << ele.index << "\t" << ele.reftype << "\t"
        << ele.vertex[0] << "   " << ele.vertex[1] << "   " << ele.vertex[2] << "\t"
        << ele.edge[0] << "   " << ele.edge[1] << "   " << ele.edge[2] << "\t"
        << ele.parent << "\t";
        for (int jChild : ele.child)
            std::cout << jChild << " ";
        std::cout << std::endl;
    }
    
    std::cout << std::endl;
    
}

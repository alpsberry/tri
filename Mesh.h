//
//  Mesh.h
//  tri
//
//  Created by GBB on 26/12/14.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#ifndef __tri__Mesh__
#define __tri__Mesh__

//#define __MESH_DEBUG

#include <iostream>
#include <fstream>
#include <cstring>
#include "problem.h"

struct Vertex {
    int index;
    int dofIndex; // index in global dof
    double x, y;  // coordinates
    int bctype;   // boundary condition type, 0 for interior vertex, 1 for dirichlet boundary, 2 for neumann boundary
};

struct Element;

struct Edge {
    int index;
    int reftype;                       // refinement type, -1 for not refined, 0, 1, 2, ... for refinement level
    int bctype;                        // boundary condition type, 0 for interior vertex, 1 for dirichlet boundary, 2 for neumann boundary
    std::vector<int> vertex;           // indices of its vertices
    std::vector<int> neighborElement;  // indices of elements sharing the edge
};

struct Element {
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

class Mesh {
    std::string _meshFilename;
    int _dimension;

    
    int initElement();
    int initEdge();
    int initVertex();
    int readRefinement(int);
    void findElementEdge(std::vector<Edge>::size_type previousLevelElementSize);
public:
    void printVertex();
    void printEdge();
    void printElement();
public:
    void calcDetBE() // calculate detBE for each element on mesh
    {
        for (Element &ele : element) {
            if (ele.reftype != constNonrefined)
                continue;
            double x1(vertex[ele.vertex[0] - 1].x), y1(vertex[ele.vertex[0] - 1].y),
            x2(vertex[ele.vertex[1] - 1].x), y2(vertex[ele.vertex[1] - 1].y),
            x3(vertex[ele.vertex[2] - 1].x), y3(vertex[ele.vertex[2] - 1].y);
            
            ele.detBE = fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)); // absolute value needed here?
        }
    }

    std::vector<Element> element;
    std::vector<Edge> edge;
    std::vector<Vertex> vertex;
    
    Mesh(Problem* prob);
    
};

#endif /* defined(__tri__Mesh__) */

//
//  DGSolvingSystem.h
//  tri
//
//  Created by GBB on 26/12/14.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#ifndef __tri__DGSolvingSystem__
#define __tri__DGSolvingSystem__

#define __DGSOLVESYS_DEBUG
 #define __DGSOLVESYS_DEBUG_EDGE
// #define __DGSOLVESYS_DEBUG_LV2

#include "BasicSolvingSystem.h"
#include "DGProblem.h"

class DGSolvingSystem: public BasicSolvingSystem {
protected:
    const int LocalDimension = 3;
    double penaltyOver6;
    
    int retrieve_dof_count_element_dofIndex(Mesh &mesh); // assign dof to each element and return the total dof
    double innerProduct(std::vector<double> x, std::vector<double> y); // the inner product of two two-dimensional vector
    double dist(double x1, double y1, double x2, double y2);
    int vertexOnEdge(Vertex ver, Vertex v1, Vertex v2); // if ver is on the line of edge, edge is given by v1 and v2
    double penaltyTerm(Edge edge, int index_E1, int index_E2, int iver, int jver, VECMATRIX f_E1, VECMATRIX f_E2);
    
    VECMATRIX elementInteg(Element ele);
    std::vector<double> elementIntegRhs(Element ele);
    int assembleElement(Element ele);
    
    int calc_ne_and_f_on_edge(Edge edge, std::vector<double> &ne, VECMATRIX &integ_e);
    int calc_ne_and_f_on_edge(Edge edge, std::vector<double> &ne, VECMATRIX &f_E1, VECMATRIX &f_E2);
    
    void initM(VECMATRIX &M11, VECMATRIX &M12, VECMATRIX &M21, VECMATRIX &M22, int dofE1, int dofE2);
    void initM(VECMATRIX &M11, int dofE1);

    int getMii(Edge edge, VECMATRIX &M, Element E1, Element E2, VECMATRIX f_E1, VECMATRIX f_E2,
               double eps, std::vector<double> ne, std::vector<double> grad_ne_E1, std::vector<double> grad_ne_E2,
               int sign1, int sing2, int sing3);
    void addMiiToMA(VECMATRIX M, Element E1, Element E2);
    int edgeInteg(Edge edge, VECMATRIX &M11, VECMATRIX &M12, VECMATRIX &M21, VECMATRIX &M22);
    int edgeInteg(Edge edge, VECMATRIX &M11, std::vector<double> &rhs);
    int assembleEdge(Edge edge);
    
    void computeError(double &errL2, double &errH1); // compute error in L2 and H1 norm
    
    int consoleOutput();  // output the result in console
    int fileOutput();     // output the result in file *.output
public:
    DGSolvingSystem(Mesh* m, Problem* p);
    void assembleStiff(); // list-stored stiffness matrix saved in ma
    void output();      // output the result
};



#endif /* defined(__tri__DGSolvingSystem__) */

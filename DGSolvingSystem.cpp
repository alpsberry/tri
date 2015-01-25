//
//  DGSolvingSystem.cpp
//  tri
//
//  Created by GBB on 26/12/14.
//  Copyright (c) 2014 Xiaolin Guo. All rights reserved.
//

#include "DGSolvingSystem.h"

using std::vector;
using std::cout;
using std::endl;

VECMATRIX DGSolvingSystem::elementInteg(Element ele)
{
    VECMATRIX vecElementInteg;
    vecElementInteg.resize(ele.localDof);
    
    // #ifdef __DGSOLVESYS_DEBUG_LV2
    //     cout << "  local dof = " << ele.localDof << endl;
    // #endif
    
    double x1(mesh -> vertex[ele.vertex[0] - 1].x), y1(mesh -> vertex[ele.vertex[0] - 1].y),
    x2(mesh -> vertex[ele.vertex[1] - 1].x), y2(mesh -> vertex[ele.vertex[1] - 1].y),
    x3(mesh -> vertex[ele.vertex[2] - 1].x), y3(mesh -> vertex[ele.vertex[2] - 1].y);
    
    VECMATRIX vecGrad(3);
    vecGrad[0].push_back(y2 - y3); vecGrad[0].push_back(x3 - x2);
    vecGrad[1].push_back(y3 - y1); vecGrad[1].push_back(x1 - x3);
    vecGrad[2].push_back(y1 - y2); vecGrad[2].push_back(x2 - x1);
    
    for (int i = 0; i < 3; i++)
        for (int j = i; j < 3; j++) {
            vecElementInteg[i].push_back( innerProduct(vecGrad[i], vecGrad[j]) / ele.detBE / 2.0);
            if (i != j) {
                // vecElementInteg[i][j] += ele.detBE / 24.0; // comment this line to set alpha = 0, so this term is omitted
                vecElementInteg[j].push_back(vecElementInteg[i][j]);
            }
            // else {
            //     vecElementInteg[i][j] += ele.detBE / 12.0; // comment this line to set alpha = 0, so this term is omitted
            // }
        }
    
    return vecElementInteg;
}

vector<double> DGSolvingSystem::elementIntegRhs(Element ele)
{
    vector<double> vecElementIntegRhs;
    
    double x1(mesh -> vertex[ele.vertex[0] - 1].x), y1(mesh -> vertex[ele.vertex[0] - 1].y),
    x2(mesh -> vertex[ele.vertex[1] - 1].x), y2(mesh -> vertex[ele.vertex[1] - 1].y),
    x3(mesh -> vertex[ele.vertex[2] - 1].x), y3(mesh -> vertex[ele.vertex[2] - 1].y);
    
    vecElementIntegRhs.resize(LocalDimension);
    vecElementIntegRhs[0] = prob -> f(x1, y1) * ele.detBE / 6.0;
    vecElementIntegRhs[1] = prob -> f(x2, y2) * ele.detBE / 6.0;
    vecElementIntegRhs[2] = prob -> f(x3, y3) * ele.detBE / 6.0;
    
    return vecElementIntegRhs;
}

int DGSolvingSystem::assembleElement(Element ele)
{
    
    VECMATRIX vecElementInteg = elementInteg(ele);
    
    for (int row = 0; row != vecElementInteg.size(); ++row)
        for (int col = 0; col != vecElementInteg[row].size(); ++col) {
            this -> addToMA(vecElementInteg[row][col], ele.dofIndex + row, ele.dofIndex + col);
        }
    
    // why commenting off this part causes a segmentation fault?
    vector<double> vecElementIntegRhs;
    vecElementIntegRhs = elementIntegRhs(ele);
    for (int i = 0; i != vecElementIntegRhs.size(); ++i) {
        this -> rh[ele.dofIndex + i] += vecElementIntegRhs[i];
    }
    
    return 0;
}

double DGSolvingSystem::dist(double x1, double y1, double x2, double y2)
{
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

int DGSolvingSystem::vertexOnEdge(Vertex ver, Vertex v1, Vertex v2)
{
    double A = v2.y - v1.y;
    double B = v1.x - v2.x;
    double C = v1.y * (v2.x - v1.x) - v1.x * (v2.y - v1.y);
    double dis = fabs( (double)(A * ver.x + B * ver.y + C) / sqrt(A * A + B * B));
    if (dis < 0.0001)
        return 1;
    return 0;
}

int DGSolvingSystem::calc_ne_and_f_on_edge(Edge edge, vector<double> &ne,
                                           VECMATRIX &f_E1, VECMATRIX &f_E2)
{
    
    if (edge.neighborElement.size() != 2) {
        cout << "invalid call in computing normal vector, element size not equal to 2" << endl;
        return 1;
    }
    
    Element &E1 = mesh -> element[ edge.neighborElement[0] - 1];
    Element &E2 = mesh -> element[ edge.neighborElement[1] - 1];
    
    double x1, x2, x3, y1, y2, y3, a2, a3, b2, b3;
    
    int index_E1_v1 = -1;
    for (int ver : E1.vertex) {
        ++index_E1_v1;
        if (edge.vertex[0] == ver || edge.vertex[1] == ver)
            continue;
        if ( vertexOnEdge(mesh -> vertex[ver - 1], mesh -> vertex[edge.vertex[0] - 1], mesh -> vertex[edge.vertex[1] - 1]) )
            continue;
        x1 = mesh -> vertex[ver - 1].x;
        y1 = mesh -> vertex[ver - 1].y;
        break;
    }
    
    x2 = mesh -> vertex[E1.vertex[(index_E1_v1 + 1) % 3] - 1].x;
    y2 = mesh -> vertex[E1.vertex[(index_E1_v1 + 1) % 3] - 1].y;
    x3 = mesh -> vertex[E1.vertex[(index_E1_v1 + 2) % 3] - 1].x;
    y3 = mesh -> vertex[E1.vertex[(index_E1_v1 + 2) % 3] - 1].y;
    
    a2 = mesh -> vertex[edge.vertex[0] - 1].x;
    b2 = mesh -> vertex[edge.vertex[0] - 1].y;
    a3 = mesh -> vertex[edge.vertex[1] - 1].x;
    b3 = mesh -> vertex[edge.vertex[1] - 1].y;
    
    if ( ((a2 - x1) * (b3 - y1) - (a3 - x1) * (b2 - y1)) < 0 ) {
        std::swap(a2, a3);
        std::swap(b2, b3);
    }
    
    ne.clear();
    ne.push_back( (b3 - b2) / 4.0 );
    ne.push_back( (a2 - a3) / 4.0 );
    
    f_E1.resize(3);
    f_E1[index_E1_v1] = vector<double> {0, 0};
    
    double length_v2v3 = dist(x2, y2, x3, y3);
    f_E1[(index_E1_v1 + 1) % 3].push_back( dist(a2, b2, x3, y3) / length_v2v3 );
    f_E1[(index_E1_v1 + 1) % 3].push_back( dist(a3, b3, x3, y3) / length_v2v3 );
    f_E1[(index_E1_v1 + 2) % 3].push_back( dist(a2, b2, x2, y2) / length_v2v3 );
    f_E1[(index_E1_v1 + 2) % 3].push_back( dist(a3, b3, x2, y2) / length_v2v3 );
    
    int index_E2_v1 = -1;
    for (int ver : E2.vertex) {
        ++index_E2_v1;
        if (edge.vertex[0] == ver || edge.vertex[1] == ver)
            continue;
        if ( vertexOnEdge(mesh -> vertex[ver - 1], mesh -> vertex[edge.vertex[0] - 1], mesh -> vertex[edge.vertex[1] - 1]) )
            continue;
        x1 = mesh -> vertex[ver - 1].x;
        y1 = mesh -> vertex[ver - 1].y;
        break;
    }
    
    x2 = mesh -> vertex[E2.vertex[(index_E2_v1 + 1) % 3] - 1].x;
    y2 = mesh -> vertex[E2.vertex[(index_E2_v1 + 1) % 3] - 1].y;
    x3 = mesh -> vertex[E2.vertex[(index_E2_v1 + 2) % 3] - 1].x;
    y3 = mesh -> vertex[E2.vertex[(index_E2_v1 + 2) % 3] - 1].y;
    
    a2 = mesh -> vertex[edge.vertex[0] - 1].x;
    b2 = mesh -> vertex[edge.vertex[0] - 1].y;
    a3 = mesh -> vertex[edge.vertex[1] - 1].x;
    b3 = mesh -> vertex[edge.vertex[1] - 1].y;
    
    if ( ((a2 - x1) * (b3 - y1) - (a3 - x1) * (b2 - y1)) < 0 ) {
        std::swap(a2, a3);
        std::swap(b2, b3);
    }
    
    f_E2.resize(3);
    f_E2[index_E2_v1] = vector<double> {0, 0};
    
    length_v2v3 = dist(x2, y2, x3, y3);
    f_E2[(index_E2_v1 + 1) % 3].push_back( dist(a2, b2, x3, y3) / length_v2v3 );
    f_E2[(index_E2_v1 + 1) % 3].push_back( dist(a3, b3, x3, y3) / length_v2v3 );
    f_E2[(index_E2_v1 + 2) % 3].push_back( dist(a2, b2, x2, y2) / length_v2v3 );
    f_E2[(index_E2_v1 + 2) % 3].push_back( dist(a3, b3, x2, y2) / length_v2v3 );
    
    return 0;
}

int DGSolvingSystem::calc_ne_and_f_on_edge(Edge edge, vector<double> &ne,
                                           VECMATRIX &f_E1)
{
    if (edge.neighborElement.size() != 1) {
        cout << "invalid call in computing nomal vector, element size not equal to 1" << endl;
        return 1;
    }
    
    Element &E1 = mesh -> element[ edge.neighborElement[0] - 1];
    double x1, x2, x3, y1, y2, y3;
    
    int index_E1_v1 = -1;
    for (int ver : E1.vertex) {
        ++index_E1_v1;
        if (edge.vertex[0] == ver) {
            continue;
        }
        if (edge.vertex[1] == ver) {
            continue;
        }
        x1 = mesh -> vertex[ver - 1].x;
        y1 = mesh -> vertex[ver - 1].y;
        break;
    }
    
    x2 = mesh -> vertex[E1.vertex[(index_E1_v1 + 1) % 3] - 1].x;
    y2 = mesh -> vertex[E1.vertex[(index_E1_v1 + 1) % 3] - 1].y;
    x3 = mesh -> vertex[E1.vertex[(index_E1_v1 + 2) % 3] - 1].x;
    y3 = mesh -> vertex[E1.vertex[(index_E1_v1 + 2) % 3] - 1].y;
    
    ne.clear();
    ne.push_back( (y3 - y2) / 2.0 );
    ne.push_back( (x2 - x3) / 2.0 );
    
    f_E1.resize(3);
    f_E1[index_E1_v1] = vector<double> {0, 0};
    
    f_E1[(index_E1_v1 + 1) % 3].push_back( 1.0 );
    f_E1[(index_E1_v1 + 1) % 3].push_back( 0.0 );
    f_E1[(index_E1_v1 + 2) % 3].push_back( 0.0 );
    f_E1[(index_E1_v1 + 2) % 3].push_back( 1.0 );
    
    return 0;
}

double DGSolvingSystem::innerProduct(vector<double> x, vector<double> y)
{
    return x[0] * y[0] + x[1] * y[1];
}

void initM(VECMATRIX &M11, VECMATRIX &M12, VECMATRIX &M21, VECMATRIX &M22, int dofE1, int dofE2)
{
    M11.resize(dofE1);
    for (int i = 0; i != dofE1; ++i)
        M11[i].resize(dofE1);
    
    M12.resize(dofE1);
    for (int i = 0; i != dofE1; ++i)
        M12[i].resize(dofE2);
    
    M21.resize(dofE2);
    for (int i = 0; i != dofE2; ++i)
        M21[i].resize(dofE1);
    
    M22.resize(dofE2);
    for (int i = 0; i != dofE2; ++i)
        M22[i].resize(dofE2);
}

void initM(VECMATRIX &M11, int dofE1)
{
    M11.resize(dofE1);
    for (int i = 0; i != dofE1; ++i)
        M11[i].resize(dofE1);
}

double DGSolvingSystem::penaltyTerm(Edge edge, int index_E1, int index_E2, int iver, int jver, VECMATRIX f_E1, VECMATRIX f_E2)
{
    if (index_E1 == index_E2)
        return 2 * f_E2[jver][0] * f_E1[iver][0] + 2 * f_E2[jver][1] * f_E1[iver][1]
        + f_E2[jver][0] * f_E1[iver][1] + f_E2[jver][1] * f_E1[iver][0];
    else
        return 2 * f_E2[jver][0] * f_E1[iver][1] + 2 * f_E2[jver][1] * f_E1[iver][0]
        + f_E2[jver][0] * f_E1[iver][0] + f_E2[jver][1] * f_E1[iver][1];
    
}

int DGSolvingSystem::getMii(Edge edge, VECMATRIX &M, Element E1, Element E2, VECMATRIX f_E1, VECMATRIX f_E2, double eps,
                            vector<double> ne, vector<double> grad_ne_E1, vector<double> grad_ne_E2, int sign1, int sign2, int sign3)
{
    int row(0), col(0);
    for (int iver = 0; iver != E1.vertex.size(); ++iver) {
        // if(mesh -> vertex[ E1.vertex[iver] - 1 ].bctype > 0)
        //  continue;
        col = 0;
        for (int jver = 0; jver != E2.vertex.size(); ++jver) {
            // if(mesh -> vertex[ E2.vertex[jver] - 1 ].bctype > 0)
            //  continue;
            M[row][col] =  sign1 * (f_E1[iver][0] + f_E1[iver][1]) * grad_ne_E2[jver]
            + sign2 * eps * (f_E2[jver][0] + f_E2[jver][1]) * grad_ne_E1[iver];
            M[row][col] += sign3 * penaltyOver6 * penaltyTerm(edge, E1.index, E2.index, iver, jver, f_E1, f_E2);
            ++col;
        }
        ++row;
    }
    
    return 0;
}

int DGSolvingSystem::edgeInteg(Edge edge, VECMATRIX &M11, VECMATRIX &M12, VECMATRIX &M21, VECMATRIX &M22)
{
    if (edge.neighborElement.size() != 2) {
        cout << "invalid call in edge integ, element size not equal to 2" << endl;
        return 1;
    }
    
    Element &E1 = mesh -> element[edge.neighborElement[0] - 1];
    Element &E2 = mesh -> element[edge.neighborElement[1] - 1];
    
    double E1_x1(mesh -> vertex[ E1.vertex[0] - 1].x),
    E1_x2(mesh -> vertex[ E1.vertex[1] - 1].x),
    E1_x3(mesh -> vertex[ E1.vertex[2] - 1].x),
    E1_y1(mesh -> vertex[ E1.vertex[0] - 1].y),
    E1_y2(mesh -> vertex[ E1.vertex[1] - 1].y),
    E1_y3(mesh -> vertex[ E1.vertex[2] - 1].y);
    double E2_x1(mesh -> vertex[ E2.vertex[0] - 1].x),
    E2_x2(mesh -> vertex[ E2.vertex[1] - 1].x),
    E2_x3(mesh -> vertex[ E2.vertex[2] - 1].x),
    E2_y1(mesh -> vertex[ E2.vertex[0] - 1].y),
    E2_y2(mesh -> vertex[ E2.vertex[1] - 1].y),
    E2_y3(mesh -> vertex[ E2.vertex[2] - 1].y);
    
    vector<double> ne; // actually n_e^* as in the report
    VECMATRIX f_E1, f_E2;
    calc_ne_and_f_on_edge(edge, ne, f_E1, f_E2); // normal vector, not unit, actualy |e| / 4 * n_e
    
    vector<double> grad_ne_E1; // \nabla \phi_i^{E_1} \cdot n_e^*
    grad_ne_E1.push_back( ((E1_y2 - E1_y3) * ne[0] + (E1_x3 - E1_x2) * ne[1]) / E1.detBE );
    grad_ne_E1.push_back( ((E1_y3 - E1_y1) * ne[0] + (E1_x1 - E1_x3) * ne[1]) / E1.detBE );
    grad_ne_E1.push_back( ((E1_y1 - E1_y2) * ne[0] + (E1_x2 - E1_x1) * ne[1]) / E1.detBE );
    
    vector<double> grad_ne_E2; // \nabla \phi_i^{E_2} \cdot n_e^*
    grad_ne_E2.push_back( ((E2_y2 - E2_y3) * ne[0] + (E2_x3 - E2_x2) * ne[1]) / E2.detBE );
    grad_ne_E2.push_back( ((E2_y3 - E2_y1) * ne[0] + (E2_x1 - E2_x3) * ne[1]) / E2.detBE );
    grad_ne_E2.push_back( ((E2_y1 - E2_y2) * ne[0] + (E2_x2 - E2_x1) * ne[1]) / E2.detBE );
    
//     #ifdef __DGSOLVESYS_DEBUG_EDGE
//          cout << "  normal vector from element " << edge.neighborElement[0]
//               << " to element " << edge.neighborElement[1] << " = ("
//               << ne[0] << ", " << ne[1] << ")" << endl;
//     #endif
    
    initM(M11, M12, M21, M22, E1.localDof, E2.localDof);
    
    const double eps = prob->epsilon;
    
    getMii(edge, M11, E1, E1, f_E1, f_E1, eps, ne, grad_ne_E1, grad_ne_E1, -1,  1,  1);
    getMii(edge, M12, E1, E2, f_E1, f_E2, eps, ne, grad_ne_E1, grad_ne_E2, -1, -1, -1);
    getMii(edge, M21, E2, E1, f_E2, f_E1, eps, ne, grad_ne_E2, grad_ne_E1,  1,  1, -1);
    getMii(edge, M22, E2, E2, f_E2, f_E2, eps, ne, grad_ne_E2, grad_ne_E2,  1, -1,  1);
    return 0;
}

int DGSolvingSystem::edgeInteg(Edge edge, VECMATRIX &M11, vector<double> &rhs)
{
    if (edge.neighborElement.size() != 1) {
        cout << "invalid call in edge integ, element size not equal to 1" << endl;
        return 1;
    }
    
    Element &E1 = mesh -> element[edge.neighborElement[0] - 1];
    
    double E1_x1(mesh -> vertex[ E1.vertex[0] - 1].x),
    E1_x2(mesh -> vertex[ E1.vertex[1] - 1].x),
    E1_x3(mesh -> vertex[ E1.vertex[2] - 1].x),
    E1_y1(mesh -> vertex[ E1.vertex[0] - 1].y),
    E1_y2(mesh -> vertex[ E1.vertex[1] - 1].y),
    E1_y3(mesh -> vertex[ E1.vertex[2] - 1].y);
    
    vector<double> ne;
    VECMATRIX f_E1;
    calc_ne_and_f_on_edge(edge, ne, f_E1); // normal vector, not unit, actualy |e| / 2 * n_e
    
    vector<double> grad_ne_E1;
    grad_ne_E1.push_back( ((E1_y2 - E1_y3) * ne[0] + (E1_x3 - E1_x2) * ne[1]) / E1.detBE );
    grad_ne_E1.push_back( ((E1_y3 - E1_y1) * ne[0] + (E1_x1 - E1_x3) * ne[1]) / E1.detBE );
    grad_ne_E1.push_back( ((E1_y1 - E1_y2) * ne[0] + (E1_x2 - E1_x1) * ne[1]) / E1.detBE );
    
//     #ifdef __DGSOLVESYS_DEBUG_EDGE
//      cout << "  normal vector on boundary edge " << edge.index
//           << " = (" << ne[0] << ", " << ne[1] << ")" << endl;
//     #endif
    
    initM(M11, E1.localDof);
    const double eps = prob->epsilon;
    
    // M11
    getMii(edge, M11, E1, E1, f_E1, f_E1, eps, ne, grad_ne_E1, grad_ne_E1, -1, 1, 1);
    
    
    double eps_int_e_gd = 0; // actually 2 * \epsilon * int_e(g_D) / |e|, 2 / |e| is not divided here since the normal vector ne is not unified
    for (int iver : E1.vertex) {
        if (iver == edge.vertex[0] || iver == edge.vertex[1])
            eps_int_e_gd += prob->gd(mesh -> vertex[iver - 1].x, mesh -> vertex[iver - 1].y);
    }
    eps_int_e_gd *= eps;
    
    rhs.resize(3);
    for (int i = 0; i < 3; i++) {
        int iver = E1.vertex[i];
        if (iver != edge.vertex[0] && iver != edge.vertex[1])
            rhs[i] = grad_ne_E1[i] * eps_int_e_gd;
        else
            rhs[i] = grad_ne_E1[i] * eps_int_e_gd + prob->sigma0 / 2.0 * prob->gd(mesh -> vertex[iver - 1].x, mesh -> vertex[iver - 1].y);
    }
    return 0;
}

void DGSolvingSystem::addMiiToMA(VECMATRIX M, Element E1, Element E2)
{
    for (int row = 0; row != M.size(); ++row)
        for (int col = 0; col != M[row].size(); ++col)
            this -> addToMA(M[row][col], E1.dofIndex + row, E2.dofIndex + col);
}

int DGSolvingSystem::assembleEdge(Edge edge)
{
    VECMATRIX M11, M12, M21, M22;
    Element &E1 = mesh -> element[edge.neighborElement[0] - 1];
    
    if (edge.neighborElement.size() == 2) {
        Element &E2 = mesh -> element[edge.neighborElement[1] - 1];
        
        edgeInteg(edge, M11, M12, M21, M22);
        
        addMiiToMA(M11, E1, E1);
        addMiiToMA(M12, E1, E2);
        addMiiToMA(M21, E2, E1);
        addMiiToMA(M22, E2, E2);
        
    } else {
        vector<double> rhs;
        edgeInteg(edge, M11, rhs);
        
        addMiiToMA(M11, E1, E1);
        
        for (int i = 0; i < 3; i++)
            rh[E1.dofIndex + i] += rhs[i];
    }
    
    
    return 0;
}

int DGSolvingSystem::retrieve_dof_count_element_dofIndex(Mesh &mesh)
{
    int dof(0);
    for (auto itEle = mesh.element.begin(); itEle != mesh.element.end(); ++itEle)
        if (itEle -> reftype == constNonrefined) {
            itEle -> dofIndex = dof;
            itEle -> localDof = LocalDimension;
            dof += LocalDimension;
        }
    
    return dof;
}

void DGSolvingSystem::assembleStiff()
{
#ifdef __DGSOLVESYS_DEBUG
    cout << "start forming system" << endl;
#endif
    
    clock_t t = clock();
    
    this -> dof = retrieve_dof_count_element_dofIndex(*mesh); // get total dof
#ifdef __DGSOLVESYS_DEBUG
    cout << " dof = " << this -> dof << endl;
#endif
    
    // initialize rh, ma
    this -> rh = new double [this -> dof];
    memset(this -> rh, 0, (this -> dof) * sizeof(double));
    this -> ma.resize(this -> dof);
    
    mesh -> calcDetBE(); //calculate det(B_E) for each element
    
    // assemble element integral related items
    int k = 1;
    for (auto it = mesh -> element.begin();
         it != mesh -> element.end(); ++it, ++k) {
        if (it -> reftype != constNonrefined)
            continue;
        
        // #ifdef __DGSOLVESYS_DEBUG_LV2
        //      cout << " assemble element " << k << endl;
        // #endif
        
        assembleElement(*it);
    }
    t = clock() - t;
    
#ifdef __DGSOLVESYS_DEBUG
    cout << "finish assembling element, t = "
    << (double) t / CLOCKS_PER_SEC << "s"
    << endl;
#endif
    
    
    t = clock();
    
    //calc penalty
    penaltyOver6 = prob->sigma0 / 6.0;
    
    // assemble edge integral related items
    k = 1;
    for (auto it = mesh -> edge.begin(); it != mesh -> edge.end(); ++it, ++k) {
        if (it -> reftype != constNonrefined)
            continue;
        
//         #ifdef __DGSOLVESYS_DEBUG_EDGE
//              cout << " assemble edge " << k << endl;
//         #endif
        assembleEdge(*it);
    }
    
    t = clock() - t;
    
#ifdef __DGSOLVESYS_DEBUG
    cout << "finish assembling edge, t = "
    << (double) t / CLOCKS_PER_SEC << "s"
    << endl;
#endif
#ifdef __DGSOLVESYS_DEBUG
    cout << "finish forming system" << endl << endl;
#endif

}

int DGSolvingSystem::consoleOutput()
{
    std::vector<Element>::iterator it;
    int k(0);
    for (it = mesh -> element.begin(); it != mesh -> element.end(); it++) {
        if (it -> reftype != constNonrefined)
            continue;
        k = 0;
        for (int ver : it -> vertex)
            if (mesh -> vertex[ver - 1].bctype == 0)
                cout << mesh -> vertex[ver - 1].x << " " << mesh -> vertex[ver - 1].y << " " << this -> x[it -> dofIndex + (k++)] << std::endl;
            else
                cout << mesh -> vertex[ver - 1].x << " " << mesh -> vertex[ver - 1].y << " "
                << prob->gd(mesh -> vertex[ver - 1].x, mesh -> vertex[ver - 1].y) << std::endl;
    }
    
    return 0;
}

int DGSolvingSystem::fileOutput()
{
    std::ofstream fout((prob->parameters.meshFilename + ".output").c_str());
    
    std::vector<Element>::iterator it;
    int k(0);
    for (it = mesh -> element.begin(); it != mesh -> element.end(); it++) {
        if (it -> reftype != constNonrefined) {
            continue;
        }
        k = 0;
        for (int ver : it -> vertex) {
            // if (mesh -> vertex[ver - 1].bctype == 0)
            fout << mesh -> vertex[ver - 1].x << " " << mesh -> vertex[ver - 1].y << " " << this -> x[it -> dofIndex + (k++)] << std::endl;
            // else
            //     fout << mesh -> vertex[ver - 1].x << " " << mesh -> vertex[ver - 1].y << " "
            //          << prob.gd(mesh -> vertex[ver - 1].x, mesh -> vertex[ver - 1].y) << std::endl;
        }
    }
    
    
    return 0;
}

void DGSolvingSystem::output()
{
    if (prob->parameters.printResults)
        consoleOutput();
    if (prob->parameters.fprintResults)
        fileOutput();
    
    BasicSolvingSystem::output();
    
    if (prob->parameters.cprintError) {
        double errL2(0), errH1(0);
        computeError(errL2, errH1);
        std::cout << "error in L2 norm = " << errL2 << std::endl
                  << "error in H1 norm = " << errH1 << std::endl;
    }
    
}

void DGSolvingSystem::computeError(double &errL2, double &errH1)
{
    std::ofstream fout((prob->parameters.meshFilename + ".err").c_str());
    
    errL2 = 0;
    errH1 = 0;
    for (Element iEle : mesh -> element) {
        if (iEle.reftype != constNonrefined)
            continue;
        Vertex &v1 = mesh -> vertex[iEle.vertex[0] - 1];
        Vertex &v2 = mesh -> vertex[iEle.vertex[1] - 1];
        Vertex &v3 = mesh -> vertex[iEle.vertex[2] - 1];
        double x1(v1.x), y1(v1.y), x2(v2.x), y2(v2.y), x3(v3.x), y3(v3.y);
        double p1(0), p2(0), p3(0);
        double r1(0), r2(0), r3(0);
        
        if (v1.bctype > 0) {
            p1 = this -> x[iEle.dofIndex];
            r1 = prob->trueSol(x1, y1) - p1;
        }
        if (v2.bctype > 0) {
            p2 = this -> x[iEle.dofIndex + 1];
            r2 = prob->trueSol(x2, y2) - p2;
        }
        if (v3.bctype > 0) {
            p3 = this -> x[iEle.dofIndex + 2];
            r3 = prob->trueSol(x3, y3) - p3;
        }
        errL2 += (r1 * r1 + r2 * r2 + r3 * r3) * iEle.detBE / 6.0;
        
        fout << v1.x << " " << v1.y << " " << r1 << endl
        << v2.x << " " << v2.y << " " << r2 << endl
        << v3.x << " " << v3.y << " " << r3 << endl;
        
        errH1 += (  pow(r1 * (y2 - y3), 2) + pow(r2 * (y3 - y1), 2) + pow(r3 * (y1 - y2), 2)
                  + pow(r1 * (x3 - x2), 2) + pow(r2 * (x1 - x3), 2) + pow(r3 * (x2 - x1), 2)  ) / 2.0 / iEle.detBE;
    }
    errL2 = sqrt(errL2);
    errH1 = sqrt(errH1);
}
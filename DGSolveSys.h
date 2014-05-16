#ifndef TRI_DGSOLVESYS_H
#define TRI_DGSOLVESYS_H

#define __DGSOLVESYS_DEBUG
// #define __DGSOLVESYS_DEBUG_EDGE
// #define __DGSOLVESYS_DEBUG_LV2

#include "solveSys.h"
#include "DGProblem.h"

class DGSolvingSystem:public BasicSolvingSystem
{
	const int LocalDimension = 3;
	double penaltyOver3, penaltyOver6;
	
	int retrive_dof_count_element_dofIndex(Mesh &mesh); // assign dof to each element and return the total dof

	void calcDetBEOnMesh(Mesh &mesh); // calculate detBE for each element on mesh

	std::vector< std::vector<double> > elementInteg(Element ele, Mesh mesh);

	std::vector<double> elementIntegRhs(Element ele, Mesh mesh, Problem& prob);

	int edgeInteg(Edge edge, Mesh mesh, Problem& prob, VECMATRIX &M11, VECMATRIX &M12, VECMATRIX &M21, VECMATRIX &M22);

	int edgeInteg(Edge edge, Mesh mesh, Problem& prob, VECMATRIX &M11, std::vector<double> &rhs);

	int assembleElement(Element ele, Mesh mesh, Problem& prob);

	int assembleEdge(Edge edge, Mesh mesh, Problem& prob);

	double innerProduct(std::vector<double> x, std::vector<double> y); // the inner product of two two-dimensional vector
	
	double penaltyTerm(Edge edge, int v1, int v2);
	
	int calcNormalVectorAndIntegOnEdge(Edge edge, Mesh mesh, std::vector<double> &ne,
		std::vector<double> &integ_e);

	int calcNormalVectorAndIntegOnEdge(Edge edge, Mesh mesh, std::vector<double> &ne,
		std::vector<double> &integ_e1, std::vector<double> &integ_e2);

	int getMii(Mesh mesh, Edge edge, VECMATRIX &M, Element E1, Element E2, std::vector<double> integ_e1, std::vector<double> integ_e2, double eps,
		std::vector<double> ne, std::vector< std::vector<double> > grad_E1, std::vector< std::vector<double> > grad_E2, int sign1, int sing2, int sing3);

	void computeError(Mesh mesh, Problem& prob, double& errL2, double& errH1); // compute error in L2 and H1 norm

	int consoleOutput(Mesh mesh, Problem& prob);
	int fileOutput(Mesh mesh, Problem& prob);

public:

	int assembleStiff(Mesh &mesh, Problem& prob);

	int triOutput(Problem& prob, Mesh mesh);
};


using std::vector;
using std::cout;
using std::endl;

void DGSolvingSystem::calcDetBEOnMesh(Mesh &mesh)
{
	for(Element &ele : mesh.element){
		if(ele.reftype != constNonrefined)
			continue;
		double x1(mesh.vertex[ele.vertex[0] - 1].x), y1(mesh.vertex[ele.vertex[0] - 1].y),
		   x2(mesh.vertex[ele.vertex[1] - 1].x), y2(mesh.vertex[ele.vertex[1] - 1].y),
		   x3(mesh.vertex[ele.vertex[2] - 1].x), y3(mesh.vertex[ele.vertex[2] - 1].y);

		ele.detBE = fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)); // absolute value needed here?
	}
}

vector< vector<double> > DGSolvingSystem::elementInteg(Element ele, Mesh mesh)
{
	VECMATRIX vecElementInteg;
	vecElementInteg.resize(ele.localDof);

	#ifdef __DGSOLVESYS_DEBUG_LV2
				cout << "  local dof = " << ele.localDof << endl;
	#endif
	
	double x1(mesh.vertex[ele.vertex[0] - 1].x), y1(mesh.vertex[ele.vertex[0] - 1].y),
		   x2(mesh.vertex[ele.vertex[1] - 1].x), y2(mesh.vertex[ele.vertex[1] - 1].y),
		   x3(mesh.vertex[ele.vertex[2] - 1].x), y3(mesh.vertex[ele.vertex[2] - 1].y);

	double rec_2detBE = 0.5 / ele.detBE;
	// double detBE_over_12 = ele.detBE / 12.0;
	// double detBE_over_24 = ele.detBE / 24.0;

	double detBE_over_12 = 0;
	double detBE_over_24 = 0;

	vector< vector<double> > vecGrad(3);
	vecGrad[0].push_back(y2 - y3);
	vecGrad[0].push_back(x3 - x2);
	vecGrad[1].push_back(y3 - y1);	
	vecGrad[1].push_back(x1 - x3);
	vecGrad[2].push_back(y1 - y2);	
	vecGrad[2].push_back(x2 - x1);

	for(int i = 0; i < 3; i++)
		for(int j = i; j < 3; j++)
		{
			vecElementInteg[i].push_back( innerProduct(vecGrad[i], vecGrad[j]) * rec_2detBE);
			if(i == j)
				vecElementInteg[i][j] += detBE_over_12;
			else{
				vecElementInteg[i][j] += detBE_over_24;
				vecElementInteg[j].push_back(vecElementInteg[i][j]);
			}
		}

	return vecElementInteg;
}

vector<double> DGSolvingSystem::elementIntegRhs(Element ele, Mesh mesh, Problem& prob)
{
	vector<double> vecElementIntegRhs;
	
	double x1(mesh.vertex[ele.vertex[0] - 1].x), y1(mesh.vertex[ele.vertex[0] - 1].y),
		   x2(mesh.vertex[ele.vertex[1] - 1].x), y2(mesh.vertex[ele.vertex[1] - 1].y),
		   x3(mesh.vertex[ele.vertex[2] - 1].x), y3(mesh.vertex[ele.vertex[2] - 1].y);

	vector<double> F(LocalDimension, 0);
	
	F[0] = prob.f((x2 + x3) / 2.0, (y2 + y3) / 2.0);
	F[1] = prob.f((x1 + x3) / 2.0, (y1 + y3) / 2.0);
	F[2] = prob.f((x1 + x2) / 2.0, (y1 + y2) / 2.0);

	double detBE_over_12 = ele.detBE / 12.0;

	for(int vi = 0; vi < LocalDimension; ++vi){
		// if(mesh.vertex[ele.vertex[vi] - 1].bctype > 0)
		// 	continue;

		double a = 0;
		for(int i = 0; i < ele.vertex.size(); i++){
			if (i == vi)
				continue;
			a += F[i];
		}

		a *= detBE_over_12;
		vecElementIntegRhs.push_back(a);

	#ifdef __DGSOLVESYS_DEBUG_LV2
			cout << "  global " << ele.vertex[vi] << " local " << vi
					  << " rh = " << a << endl;
	#endif
	}

	return vecElementIntegRhs;
}

int DGSolvingSystem::assembleElement(Element ele, Mesh mesh, Problem& prob)
{

	VECMATRIX vecElementInteg = elementInteg(ele, mesh);
	
	for(int row = 0; row != vecElementInteg.size(); ++row)
		for(int col = 0; col != vecElementInteg[row].size(); ++col){

	#ifdef __DGSOLVESYS_DEBUG_LV2
				cout << "  add to " << ele.dofIndex + row << ", " << ele.dofIndex + col << endl;
	#endif
			this -> addToMA(vecElementInteg[row][col], ele.dofIndex + row, ele.dofIndex + col);
		}

	// why commenting off this part causes a segmentation fault?
	vector<double> vecElementIntegRhs;
	vecElementIntegRhs = elementIntegRhs(ele, mesh, prob);
	for(int i = 0; i != vecElementIntegRhs.size(); ++i){

	#ifdef __DGSOLVESYS_DEBUG_LV2
			cout << "  add rh to " << ele.dofIndex + i << endl;
	#endif
		this -> rh[ele.dofIndex + i] += vecElementIntegRhs[i];
	}

	return 0;
}

double dist(Vertex v1, Vertex v2)
{
	return sqrt((v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y));
}

int vertexOnEdge(Vertex ver, Vertex e1, Vertex e2)
{
	double A = e2.y - e1.y;
	double B = e1.x - e2.x;
	double C = e1.y * (e2.x - e1.x) - e1.x * (e2.y - e1.y);
	double dis = fabs( 1.0 * (A * ver.x + B * ver.y + C) / sqrt(A * A + B * B));
	if(dis < 0.0001)
		return 1;
	return 0;
}

int DGSolvingSystem::calcNormalVectorAndIntegOnEdge(Edge edge, Mesh mesh, vector<double> &ne,
	vector<double> &integ_e1, vector<double> &integ_e2)
{
	if(edge.neighborElement.size() != 2){
		cout << "invalid call in computing normal vector, element size not equal to 2" << endl;
		return 1;
	}

	Element &E1 = mesh.element[ edge.neighborElement[0] - 1];
	Element &E2 = mesh.element[ edge.neighborElement[1] - 1];

	double x1, x2, x3, y1, y2, y3;
	
	vector<double> integ_e1_tmp(LocalDimension, 1);
	vector<double> integ_e2_tmp(LocalDimension, 1);
	
	int k = -1;
	for(int ver:E1.vertex)
	{
		++k;
		if(edge.vertex[0] == ver){
			continue;
		}
		if(edge.vertex[1] == ver){
			continue;
		}
		if( vertexOnEdge(mesh.vertex[ver - 1], mesh.vertex[edge.vertex[0] - 1], mesh.vertex[edge.vertex[1] - 1]) )
			continue;
		x1 = mesh.vertex[ver - 1].x;
		y1 = mesh.vertex[ver - 1].y;

		integ_e1_tmp[k] = 0;
	}
	
	k = -1;
	for(int ver:E2.vertex)
	{
		++k;
		if(edge.vertex[0] == ver){
			continue;
		}
		if(edge.vertex[1] == ver){
			continue;
		}
		if( vertexOnEdge(mesh.vertex[ver - 1], mesh.vertex[edge.vertex[0] - 1], mesh.vertex[edge.vertex[1] - 1]) )
			continue;
		integ_e2_tmp[k] = 0;
	}

	x2 = mesh.vertex[edge.vertex[0] - 1].x;
	y2 = mesh.vertex[edge.vertex[0] - 1].y;
	x3 = mesh.vertex[edge.vertex[1] - 1].x;
	y3 = mesh.vertex[edge.vertex[1] - 1].y;

	if( ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) < 0 )
	{
		std::swap(x2, x3);
		std::swap(y2, y3);
	}

	ne.clear();
	ne.push_back( (y3 - y2) / 4.0 );
	ne.push_back( (x2 - x3) / 4.0 );

	integ_e1.clear();
	integ_e1 = integ_e1_tmp;
	integ_e2.clear();
	integ_e2 = integ_e2_tmp;

	// if(std::find(E1.vertex.begin(), E1.vertex.end(), edge.vertex[0]) == E1.vertex.end()
	// 	|| std::find(E1.vertex.begin(), E1.vertex.end(), edge.vertex[1]) == E1.vertex.end())
	// {
		for(int v2 = 0; v2 < E1.vertex.size(); v2++)
		{
			if(integ_e1[v2] == 0)
				continue;
			for(int v3 = 0; v3 < E1.vertex.size(); v3++){
				if(v2 == v3 || integ_e1[v3] == 0)
					continue;

				double length_v2v3 = dist(mesh.vertex[E1.vertex[v2] - 1], mesh.vertex[E1.vertex[v3] - 1]);
				// double length_edge = dist(mesh.vertex[edge.vertex[0] - 1], mesh.vertex[edge.vertex[1] - 1]);
				double f1 = dist(mesh.vertex[E1.vertex[v3] - 1], mesh.vertex[edge.vertex[0] - 1]);
				double f2 = dist(mesh.vertex[E1.vertex[v3] - 1], mesh.vertex[edge.vertex[1] - 1]);
				integ_e1[v2] = (f1 + f2) * 1.0 / length_v2v3;
			}
		}
	// }

	// if(std::find(E2.vertex.begin(), E2.vertex.end(), edge.vertex[0]) == E2.vertex.end()
	// 	|| std::find(E2.vertex.begin(), E2.vertex.end(), edge.vertex[1]) == E2.vertex.end())
	// {
		for(int v2 = 0; v2 < E2.vertex.size(); v2++)
		{
			if(integ_e2[v2] == 0)
				continue;
			for(int v3 = 0; v3 < E2.vertex.size(); v3++){
				if(v2 == v3 || integ_e2[v3] == 0)
					continue;

				double length_v2v3 = dist(mesh.vertex[E2.vertex[v2] - 1], mesh.vertex[E2.vertex[v3] - 1]);
				// double length_edge = dist(mesh.vertex[edge.vertex[0] - 1], mesh.vertex[edge.vertex[1] - 1]);
				double f1 = dist(mesh.vertex[E2.vertex[v3] - 1], mesh.vertex[edge.vertex[0] - 1]);
				double f2 = dist(mesh.vertex[E2.vertex[v3] - 1], mesh.vertex[edge.vertex[1] - 1]);
				integ_e2[v2] = (f1 + f2) * 1.0 / length_v2v3;
			}
		}
	// }
	return 0;
}

int DGSolvingSystem::calcNormalVectorAndIntegOnEdge(Edge edge, Mesh mesh, vector<double> &ne,
	vector<double> &integ_e)
{
	if(edge.neighborElement.size() != 1){
		cout << "invalid call in computing nomal vector, element size not equal to 1" << endl;
		return 1;
	}

	Element &ele = mesh.element[ edge.neighborElement[0] - 1];
	double x1, x2, x3, y1, y2, y3;
	
	vector<double> integ_e_tmp(LocalDimension, 1);
	
	int k = -1;
	for(int ver:ele.vertex)
	{
		++k;
		if(edge.vertex[0] == ver){
			continue;
		}
		if(edge.vertex[1] == ver){
			continue;
		}
		x1 = mesh.vertex[ver - 1].x;
		y1 = mesh.vertex[ver - 1].y;

		integ_e_tmp[k] = 0;
	}

	x2 = mesh.vertex[edge.vertex[0] - 1].x;
	y2 = mesh.vertex[edge.vertex[0] - 1].y;
	x3 = mesh.vertex[edge.vertex[1] - 1].x;
	y3 = mesh.vertex[edge.vertex[1] - 1].y;

	if( ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) < 0 )
	{
		std::swap(x2, x3);
		std::swap(y2, y3);
	}

	ne.clear();
	ne.push_back( (y3 - y2) / 2.0 );
	ne.push_back( (x2 - x3) / 2.0 );

	integ_e = integ_e_tmp;

	return 0;
}

double DGSolvingSystem::innerProduct(vector<double> x, vector<double> y)
{
	return x[0] * y[0] + x[1] * y[1];
}

void initM(VECMATRIX &M11, VECMATRIX &M12, VECMATRIX &M21, VECMATRIX &M22, int dofE1, int dofE2)
{
	M11.resize(dofE1);
	for(int i = 0; i != dofE1; ++i)
		M11[i].resize(dofE1);

	M12.resize(dofE1);
	for(int i = 0; i != dofE1; ++i)
		M12[i].resize(dofE2);

	M21.resize(dofE2);
	for(int i = 0; i != dofE2; ++i)
		M21[i].resize(dofE1);

	M22.resize(dofE2);
	for(int i = 0; i != dofE2; ++i)
		M22[i].resize(dofE2);
}

void initM(VECMATRIX &M11, int dofE1)
{
	M11.resize(dofE1);
	for(int i = 0; i != dofE1; ++i)
		M11[i].resize(dofE1);
}

double DGSolvingSystem::penaltyTerm(Edge edge, int v1, int v2)
{
	// if any of vi, v2 is not on the edge, then the penalty term is 0
	if((v1 != edge.vertex[0] && v1 != edge.vertex[1]) || (v2 != edge.vertex[0] && v2 != edge.vertex[1]))
		return 0;

	if(v1 == v2)
		return penaltyOver3; // i.e. sigma0 / 3.0
	else
		return penaltyOver6; // i.e. sigma0 / 6.0
}

int DGSolvingSystem::getMii(Mesh mesh, Edge edge, VECMATRIX &M, Element E1, Element E2, vector<double> integ_e1, vector<double> integ_e2, double eps,
	vector<double> ne, vector< vector<double> > grad_E1, vector< vector<double> > grad_E2, int sign1, int sign2, int sign3)
{
	int row(0), col(0);
	for(int iver = 0; iver != E1.vertex.size(); ++iver){
		// if(mesh.vertex[ E1.vertex[iver] - 1 ].bctype > 0)
		// 	continue;
		col = 0;
		for(int jver = 0; jver != E2.vertex.size(); ++jver)
		{
			// if(mesh.vertex[ E2.vertex[jver] - 1 ].bctype > 0)
			// 	continue;
			M[row][col] = sign1 * integ_e1[iver] * innerProduct(grad_E2[jver], ne) + sign2 * eps * integ_e2[jver] * innerProduct(grad_E1[iver], ne);
			M[row][col] += sign3 * penaltyTerm(edge, E2.vertex[jver], E1.vertex[iver]);
			++col;
		}
		++row;
	}

	return 0;
}

int DGSolvingSystem::edgeInteg(Edge edge, Mesh mesh, Problem& prob, VECMATRIX &M11, VECMATRIX &M12, VECMATRIX &M21, VECMATRIX &M22)
{
	if(edge.neighborElement.size() != 2){
		cout << "invalid call in edge integ, element size not equal to 2" << endl;
		return 1;
	}

	Element &E1 = mesh.element[edge.neighborElement[0] - 1];
	Element &E2 = mesh.element[edge.neighborElement[1] - 1];

	double E1_x1(mesh.vertex[ E1.vertex[0] - 1].x),
		   E1_x2(mesh.vertex[ E1.vertex[1] - 1].x),
		   E1_x3(mesh.vertex[ E1.vertex[2] - 1].x),
		   E1_y1(mesh.vertex[ E1.vertex[0] - 1].y),
		   E1_y2(mesh.vertex[ E1.vertex[1] - 1].y),
		   E1_y3(mesh.vertex[ E1.vertex[2] - 1].y);
	double E2_x1(mesh.vertex[ E2.vertex[0] - 1].x),
		   E2_x2(mesh.vertex[ E2.vertex[1] - 1].x),
		   E2_x3(mesh.vertex[ E2.vertex[2] - 1].x),
		   E2_y1(mesh.vertex[ E2.vertex[0] - 1].y),
		   E2_y2(mesh.vertex[ E2.vertex[1] - 1].y),
		   E2_y3(mesh.vertex[ E2.vertex[2] - 1].y);

	double rec_detBE1 = 1.0 / E1.detBE;	
	double rec_detBE2 = 1.0 / E2.detBE;
	
	vector< vector<double> > grad_E1(LocalDimension);
	grad_E1[0].push_back((E1_y2 - E1_y3) * rec_detBE1);
	grad_E1[0].push_back((E1_x3 - E1_x2) * rec_detBE1);
	grad_E1[1].push_back((E1_y3 - E1_y1) * rec_detBE1);	
	grad_E1[1].push_back((E1_x1 - E1_x3) * rec_detBE1);
	grad_E1[2].push_back((E1_y1 - E1_y2) * rec_detBE1);	
	grad_E1[2].push_back((E1_x2 - E1_x1) * rec_detBE1);

	vector< vector<double> > grad_E2(LocalDimension);
	grad_E2[0].push_back((E2_y2 - E2_y3) * rec_detBE2);
	grad_E2[0].push_back((E2_x3 - E2_x2) * rec_detBE2);
	grad_E2[1].push_back((E2_y3 - E2_y1) * rec_detBE2);	
	grad_E2[1].push_back((E2_x1 - E2_x3) * rec_detBE2);
	grad_E2[2].push_back((E2_y1 - E2_y2) * rec_detBE2);	
	grad_E2[2].push_back((E2_x2 - E2_x1) * rec_detBE2);

	vector<double> ne;
	vector<double> integ_e1, integ_e2;	
	calcNormalVectorAndIntegOnEdge(edge, mesh, ne, integ_e1, integ_e2); // normal vector, not unit, actualy |e| / 4 * n_e
	
	#ifdef __DGSOLVESYS_DEBUG_EDGE
			cout << "  normal vector from element " << edge.neighborElement[0]
				 << " to element " << edge.neighborElement[1] << " = ("
				 << ne[0] << ", " << ne[1] << ")" << endl;
	#endif
	#ifdef __DGSOLVESYS_DEBUG_EDGE
			cout << "  vertices of E1 = (" << E1.vertex[0]
				 << ", " << E1.vertex[1] << ", "
				 << E1.vertex[2] << "), integ on edge = ("
				 << integ_e1[0] << ", " << integ_e1[1]
				 << ", " << integ_e1[2] << ")" << endl;
	#endif
	
	#ifdef __DGSOLVESYS_DEBUG_EDGE
			cout << "  vertices of E2 = (" << E2.vertex[0]
				 << ", " << E2.vertex[1] << ", "
				 << E2.vertex[2] << "), integ on edge = ("
				 << integ_e2[0] << ", " << integ_e2[1]
				 << ", " << integ_e2[2] << ")" << endl;
	#endif

	initM(M11, M12, M21, M22, E1.localDof, E2.localDof);

	const double eps = prob.epsilon;
	
	getMii(mesh, edge, M11, E1, E1, integ_e1, integ_e1, eps, ne, grad_E1, grad_E1, -1,  1,  1);
	getMii(mesh, edge, M12, E1, E2, integ_e1, integ_e2, eps, ne, grad_E1, grad_E2, -1, -1, -1);
	getMii(mesh, edge, M21, E2, E1, integ_e2, integ_e1, eps, ne, grad_E2, grad_E1,  1,  1, -1);
	getMii(mesh, edge, M22, E2, E2, integ_e2, integ_e2, eps, ne, grad_E2, grad_E2,  1, -1,  1);


	return 0;
}

int DGSolvingSystem::edgeInteg(Edge edge, Mesh mesh, Problem& prob, VECMATRIX &M11, vector<double> &rhs)
{
	if(edge.neighborElement.size() != 1){
		cout << "invalid call in edge integ, element size not equal to 1" << endl;
		return 1;
	}

	Element &E1 = mesh.element[edge.neighborElement[0] - 1];

	double E1_x1(mesh.vertex[ E1.vertex[0] - 1].x),
		   E1_x2(mesh.vertex[ E1.vertex[1] - 1].x),
		   E1_x3(mesh.vertex[ E1.vertex[2] - 1].x),
		   E1_y1(mesh.vertex[ E1.vertex[0] - 1].y),
		   E1_y2(mesh.vertex[ E1.vertex[1] - 1].y),
		   E1_y3(mesh.vertex[ E1.vertex[2] - 1].y);

	double rec_detBE1 = 1.0 / E1.detBE;	
	
	vector< vector<double> > grad_E1(LocalDimension);
	grad_E1[0].push_back((E1_y2 - E1_y3) * rec_detBE1);
	grad_E1[0].push_back((E1_x3 - E1_x2) * rec_detBE1);
	grad_E1[1].push_back((E1_y3 - E1_y1) * rec_detBE1);	
	grad_E1[1].push_back((E1_x1 - E1_x3) * rec_detBE1);
	grad_E1[2].push_back((E1_y1 - E1_y2) * rec_detBE1);	
	grad_E1[2].push_back((E1_x2 - E1_x1) * rec_detBE1);

	vector<double> ne;
	vector<double> integ_e1;	
	calcNormalVectorAndIntegOnEdge(edge, mesh, ne, integ_e1); // normal vector, not unit, actualy |e| / 2 * n_e

	// #ifdef __DGSOLVESYS_DEBUG_EDGE
	// 	cout << "  normal vector on boundary edge " << edge.index
	// 		 << " = (" << ne[0] << ", " << ne[1] << ")" << endl;
	// #endif

	// #ifdef __DGSOLVESYS_DEBUG_EDGE
	// 		cout << "  vertices of E1 = (" << E1.vertex[0]
	// 			 << ", " << E1.vertex[1] << ", "
	// 			 << E1.vertex[2] << "), integ on edge = ("
	// 			 << integ_e1[0] << ", " << integ_e1[1]
	// 			 << ", " << integ_e1[2] << ")" << endl;
	// #endif

	initM(M11, E1.localDof);
	const double eps = prob.epsilon;
	// M11

	getMii(mesh, edge, M11, E1, E1, integ_e1, integ_e1, eps, ne, grad_E1, grad_E1, -1, 1, 1);

	double eps_int_e_gd = 0; // actually 2 * \epsilon * int_e(g_D) / |e|, 2 / |e| is not divided here since the normal vector ne is not unified
	for(int iver: E1.vertex){
		if(iver == edge.vertex[0] || iver == edge.vertex[1])
			eps_int_e_gd += prob.gd(mesh.vertex[iver - 1].x, mesh.vertex[iver - 1].y);
	}
	eps_int_e_gd *= eps;
	
	rhs.resize(3);
	for(int i = 0; i < 3; i++)
	{
		int iver = E1.vertex[i];
		if(iver != edge.vertex[0] && iver != edge.vertex[1])
			rhs[i] = innerProduct(grad_E1[i], ne) * eps_int_e_gd;
		else
			rhs[i] = innerProduct(grad_E1[i], ne) * eps_int_e_gd + prob.sigma0 / 2.0 * prob.gd(mesh.vertex[iver - 1].x, mesh.vertex[iver - 1].y);
	}
	return 0;
}

int DGSolvingSystem::assembleEdge(Edge edge, Mesh mesh, Problem& prob)
{
	VECMATRIX M11, M12, M21, M22;
	Element &E1 = mesh.element[edge.neighborElement[0] - 1];

	if(edge.neighborElement.size() == 2){
		Element &E2 = mesh.element[edge.neighborElement[1] - 1];

		edgeInteg(edge, mesh, prob, M11, M12, M21, M22);
		
		// M11
		// #ifdef __DGSOLVESYS_DEBUG_EDGE
		// 	cout << "  M11" << endl;
		// #endif
		for(int row = 0; row != M11.size(); ++row)
			for(int col = 0; col != M11[row].size(); ++col){
				this -> addToMA(M11[row][col], E1.dofIndex + row, E1.dofIndex + col);
			// #ifdef __DGSOLVESYS_DEBUG_EDGE
			// 		cout << "   ( " << E1.dofIndex + row
			// 			 << " , " << E1.dofIndex + col << ") = "
			// 			 << M11[row][col] << endl;
			// #endif

			}

		// M12
		// #ifdef __DGSOLVESYS_DEBUG_EDGE
		// 	cout << "  M12" << endl;
		// #endif
		for(int row = 0; row != M12.size(); ++row)
			for(int col = 0; col != M12[row].size(); ++col){
				this -> addToMA(M12[row][col], E1.dofIndex + row, E2.dofIndex + col);
			// #ifdef __DGSOLVESYS_DEBUG_EDGE
			// 		cout << "   ( " << E1.dofIndex + row
			// 			 << " , " << E2.dofIndex + col << ") = "
			// 			 << M12[row][col] << endl;
			// #endif
			}

		// M21
		// #ifdef __DGSOLVESYS_DEBUG_EDGE
		// 	cout << "  M21" << endl;
		// #endif
		for(int row = 0; row != M21.size(); ++row)
			for(int col = 0; col != M21[row].size(); ++col){
				this -> addToMA(M21[row][col], E2.dofIndex + row, E1.dofIndex + col);
			// #ifdef __DGSOLVESYS_DEBUG_EDGE
			// 		cout << "   ( " << E2.dofIndex + row
			// 			 << " , " << E1.dofIndex + col << ") = "
			// 			 << M21[row][col] << endl;
			// #endif
			}
		
		// M22
		// #ifdef __DGSOLVESYS_DEBUG_EDGE
		// 	cout << "  M22" << endl;
		// #endif
		for(int row = 0; row != M22.size(); ++row)
			for(int col = 0; col != M22[row].size(); ++col){
				this -> addToMA(M22[row][col], E2.dofIndex + row, E2.dofIndex + col);
			// #ifdef __DGSOLVESYS_DEBUG_EDGE
			// 		cout << "   ( " << E2.dofIndex + row
			// 			 << " , " << E2.dofIndex + col << ") = "
			// 			 << M22[row][col] << endl;
			// #endif
			}
	}
	else{
		vector<double> rhs;
		edgeInteg(edge, mesh, prob, M11, rhs);
		
		// // M11
		// #ifdef __DGSOLVESYS_DEBUG_EDGE
		// 	cout << "  M11" << endl;
		// #endif
		for(int row = 0; row != M11.size(); ++row)
			for(int col = 0; col != M11[row].size(); ++col){
				this -> addToMA(M11[row][col], E1.dofIndex + row, E1.dofIndex + col);
			// #ifdef __DGSOLVESYS_DEBUG_EDGE
			// 		cout << "   ( " << E1.dofIndex + row
			// 			 << " , " << E1.dofIndex + col << ") = "
			// 			 << M11[row][col] << endl;
			// #endif

			}
		for(int i = 0; i < 3; i++)
			rh[E1.dofIndex + i] += rhs[i];
	}


	return 0;
}

int DGSolvingSystem::retrive_dof_count_element_dofIndex(Mesh &mesh)
{
	int dof(0);
	for(auto itEle = mesh.element.begin(); itEle != mesh.element.end(); ++itEle)
	{
		if(itEle -> reftype == constNonrefined)
		{
			itEle -> dofIndex = dof;
			itEle -> localDof = LocalDimension;
			dof += LocalDimension;
		}
	}

	return dof;
}

int DGSolvingSystem::assembleStiff(Mesh &mesh, Problem& prob)
{
	#ifdef __DGSOLVESYS_DEBUG
		cout << "start forming system" << endl;
	#endif

	clock_t t = clock();

	this -> dof = retrive_dof_count_element_dofIndex(mesh); // get total dof
	#ifdef __DGSOLVESYS_DEBUG
		cout << " dof = " << this -> dof << endl;
	#endif

	// initialize rh, ma
	this -> rh = new double [this -> dof];
	memset(this -> rh, 0, (this -> dof) * sizeof(double));
	this -> ma.resize(this -> dof);

	calcDetBEOnMesh(mesh); //calculate det(B_E) for each element
	
	// assemble element integral related items
	int k = 1;
	for(auto it = mesh.element.begin();
		it != mesh.element.end(); ++it, ++k)
	{
		if(it -> reftype != constNonrefined)
			continue;

	// #ifdef __DGSOLVESYS_DEBUG_LV2
	// 		cout << " assemble element " << k << endl;
	// #endif

		assembleElement(*it, mesh, prob);
	}
	t = clock() - t;
	
	#ifdef __DGSOLVESYS_DEBUG
		cout << "finish assembling element, t = "
				  << (double) t / CLOCKS_PER_SEC << "s"
				  << endl;
	#endif


	t = clock();
	
	//calc penalty
	penaltyOver3 = prob.sigma0 / 3.0;
	penaltyOver6 = prob.sigma0 / 6.0;

	// assemble edge integral related items
	k = 1;
	for(auto it = mesh.edge.begin(); it != mesh.edge.end(); ++it, ++k)
	{
		if(it -> reftype != constNonrefined)
			continue;

	// #ifdef __DGSOLVESYS_DEBUG_EDGE
	// 		cout << " assemble edge " << k << endl;
	// #endif
		assembleEdge(*it, mesh, prob);
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
	return 0;
}

int DGSolvingSystem::consoleOutput(Mesh mesh, Problem& prob)
{
	std::vector<Element>::iterator it;
	int k(0);
	for(it = mesh.element.begin(); it != mesh.element.end(); it++){
		if(it -> reftype != constNonrefined)
			continue;
		k = 0;
		for(int ver : it -> vertex)
			if(mesh.vertex[ver - 1].bctype == 0)
				cout << mesh.vertex[ver - 1].x << " " << mesh.vertex[ver - 1].y << " " << this -> x[it -> dofIndex + (k++)] << std::endl;
			else
				cout << mesh.vertex[ver - 1].x << " " << mesh.vertex[ver - 1].y << " "
					 << prob.gd(mesh.vertex[ver - 1].x, mesh.vertex[ver - 1].y) << std::endl;
	}

	return 0;
}

int DGSolvingSystem::fileOutput(Mesh mesh, Problem& prob)
{
	std::ofstream fout((prob.parameters.meshFilename + ".output").c_str());

	std::vector<Element>::iterator it;
	int k(0);
	for(it = mesh.element.begin(); it != mesh.element.end(); it++){
		if(it -> reftype != constNonrefined){
			continue;
		}
		k = 0;
		for(int ver : it -> vertex){
			if(mesh.vertex[ver - 1].bctype == 0)
				fout << mesh.vertex[ver - 1].x << " " << mesh.vertex[ver - 1].y << " " << this -> x[it -> dofIndex + (k++)] << std::endl;
			else
				fout << mesh.vertex[ver - 1].x << " " << mesh.vertex[ver - 1].y << " "
					 << prob.gd(mesh.vertex[ver - 1].x, mesh.vertex[ver - 1].y) << std::endl;
		}
	}
		

	return 0;		
}

int DGSolvingSystem::triOutput(Problem& prob, Mesh mesh)
{
	if(prob.parameters.printResults)
		consoleOutput(mesh, prob);
	if(prob.parameters.fprintResults)
		fileOutput(mesh, prob);
	if(prob.parameters.fprintMA)
		this -> fileOutputMA(mesh, prob);
	if(prob.parameters.fprintRH)
		this -> fileOutputRH(mesh, prob);
	if(prob.parameters.fprintTriplet)
		this -> fileOutputTriplet(mesh, prob);

	if(prob.parameters.cprintError)
	{
		double errL2(0), errH1(0);
		computeError(mesh, prob, errL2, errH1);
		std::cout << "error in L2 norm = " << errL2 << std::endl
		          << "error in H1 norm = " << errH1 << std::endl;
	}

	return 0;
}

void DGSolvingSystem::computeError(Mesh mesh, Problem& prob, double& errL2, double& errH1)
{
	std::ofstream fout((prob.parameters.meshFilename + ".err").c_str());

	errL2 = 0;
	errH1 = 0;
	for(Element iEle:mesh.element){
		if(iEle.reftype != constNonrefined)
			continue;
		Vertex &v1 = mesh.vertex[iEle.vertex[0] - 1];
		Vertex &v2 = mesh.vertex[iEle.vertex[1] - 1];
		Vertex &v3 = mesh.vertex[iEle.vertex[2] - 1];
		double x1(v1.x), y1(v1.y), x2(v2.x), y2(v2.y), x3(v3.x), y3(v3.y);
		double p1(0), p2(0), p3(0); 
		double r1(0), r2(0), r3(0); 
		if(v1.bctype == 0){
			p1 = this -> x[iEle.dofIndex];
			r1 = prob.trueSol(x1, y1) - p1;
		}
		if(v2.bctype == 0){
			p2 = this -> x[iEle.dofIndex + 1];
			r2 = prob.trueSol(x2, y2) - p2;
		}
		if(v3.bctype == 0){
			p3 = this -> x[iEle.dofIndex + 2];
			r3 = prob.trueSol(x3, y3) - p3;
		}
	
		errL2 += (r1 * r1 + r2 * r2 + r3 * r3) * iEle.detBE / 6.0;

		fout << v1.x << " " << v1.y << " " << r1 << std::endl;
		fout << v2.x << " " << v2.y << " " << r2 << std::endl;
		fout << v3.x << " " << v3.y << " " << r3 << std::endl;
		
		errH1 += (  pow(r1 * (y2 - y3), 2) + pow(r2 * (y3 - y1), 2) + pow(r3 * (y1 - y2), 2)
	    	      + pow(r1 * (x3 - x2), 2) + pow(r2 * (x1 - x3), 2) + pow(r3 * (x2 - x1), 2)  ) / 2.0 / iEle.detBE;		
	}
	errL2 = sqrt(errL2);
	errH1 = sqrt(errH1);
}

#endif /* TRI_DGSOLVESYS_H */

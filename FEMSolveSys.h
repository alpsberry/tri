#ifndef TRI_FEMSOLVESYS_H
#define TRI_FEMSOLVESYS_H

#define __FEMSOLVESYS_DEBUG
// #define __FEMSOLVESYS_DEBUG_LV2

// the problem to be solved here is
//   -\Delta u + u = 3\cos x \sin y, (x,y) \in \Omega = [\frac{\pi}{2}, \frac{3\pi}{2}] \times [0, \pi]
//   u = 0, (x,y) \in \Gamma
// the variational form
//   \int\nolimits_{\Omega} \nabla u \nabla v + uv = \int\nolimits_{\Omega} fv

// here vertex.index is associated with its index in dof, not the index in all the vertices
// thus only vertex with bctype=0 is defined with this value

#include "solveSys.h"

using std::vector;

template <typename MyProblem>
class FEMSolvingSystem: public BasicSolvingSystem<MyProblem>
{
public:
	VECMATRIX integA(Element &ele, Mesh<MyProblem> mesh, MyProblem prob);

	double integARH(Element ele, Mesh<MyProblem> mesh, int vi, MyProblem prob);

	int getStiff(Element& ele, Mesh<MyProblem> &mesh, MyProblem prob);

	int retrive_dof_count_vertex_index(Mesh<MyProblem> &mesh);

	int assembleStiff(Mesh<MyProblem> &mesh, MyProblem prob);

	int consoleOutput(Mesh<MyProblem> mesh, MyProblem prob);

	int fileOutput(Mesh<MyProblem> mesh, MyProblem prob);

	int triOutput(MyProblem prob, Mesh<MyProblem> mesh);

	double computeError(Mesh<MyProblem> mesh, MyProblem prob);

};

double innerProduct(vector<double> x, vector<double> y)
{
	return x[0] * y[0] + x[1] * y[1];
}

template <typename MyProblem>
void calcDetBEOnMesh(Mesh<MyProblem> &mesh)
{
	for(Element &ele : mesh.element){
		double x1(mesh.vertex[ele.vertex[0] - 1].x), y1(mesh.vertex[ele.vertex[0] - 1].y),
		   x2(mesh.vertex[ele.vertex[1] - 1].x), y2(mesh.vertex[ele.vertex[1] - 1].y),
		   x3(mesh.vertex[ele.vertex[2] - 1].x), y3(mesh.vertex[ele.vertex[2] - 1].y);

		ele.detBE = fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)); // absolute value needed here?
	}
}


// gradInteg = \int_E \nabla \phi_i \nabla v
// timeInteg = \int_E uv
// nodal basis are used here
template <typename MyProblem>
VECMATRIX FEMSolvingSystem<MyProblem>::integA(Element &ele, Mesh<MyProblem> mesh, MyProblem prob)
{
	VECMATRIX vecIntegA;
	vecIntegA.resize(3);

	double x1(mesh.vertex[ele.vertex[0] - 1].x), y1(mesh.vertex[ele.vertex[0] - 1].y),
		   x2(mesh.vertex[ele.vertex[1] - 1].x), y2(mesh.vertex[ele.vertex[1] - 1].y),
		   x3(mesh.vertex[ele.vertex[2] - 1].x), y3(mesh.vertex[ele.vertex[2] - 1].y);

	// if(ele.detBE == 0){
	// 	ele.detBE = fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
	// }

	double rec_2detBE = 0.5 / ele.detBE;
	double detBE_over_12 = ele.detBE / 12.0;
	double detBE_over_24 = ele.detBE / 24.0;

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
			vecIntegA[i].push_back( innerProduct(vecGrad[i], vecGrad[j]) * rec_2detBE);
			if(i == j)
				vecIntegA[i][j] += detBE_over_12;
			else{
				vecIntegA[i][j] += detBE_over_24;
				vecIntegA[j].push_back(vecIntegA[i][j]);
			}
		}

	return vecIntegA;
}

// calc \int_E f\phi
// \int_E f\phi = |detB_E| / 6 * (\sum_{i=1}^3 f(m_i)\phi(m_i))
// where m_i is the midpoint of each edge
template <typename MyProblem>
double FEMSolvingSystem<MyProblem>::integARH(Element ele, Mesh<MyProblem> mesh, int vi, MyProblem prob)
{
	double x1(mesh.vertex[ele.vertex[vi] - 1].x), y1(mesh.vertex[ele.vertex[vi] - 1].y);

	// double a = 0;
	// for(int i = 0; i < ele.vertex.size(); i++){
	// 	if (i == vi)
	// 		continue;
	// 	a += prob.f((x1 + mesh.vertex[ele.vertex[i] - 1].x) / 2.0, (y1 + mesh.vertex[ele.vertex[i] - 1].y) / 2.0);
	// }

	// return a * ele.detBE / 12.0;

	return prob.f(x1, y1) * ele.detBE / 6.0;
}

template <typename MyProblem>
int FEMSolvingSystem<MyProblem>::getStiff(Element& ele, Mesh<MyProblem> &mesh, MyProblem prob){

#ifdef __FEMSOLVESYS_DEBUG_LV2
	std::cout << " assemble element " << ele.index << std::endl;
#endif	

	VECMATRIX vecIntegA = integA(ele, mesh, prob);
	for(int i = 0; i < ele.vertex.size(); i++){
		if(mesh.vertex[ele.vertex[i] - 1].bctype > 0)
			continue;
		for(int j = 0; j < ele.vertex.size(); j++)
		{
			if(mesh.vertex[ele.vertex[j] - 1].bctype > 0)
				continue;
			this -> addToMA(vecIntegA[i][j], mesh.vertex[ele.vertex[i] - 1].dofIndex,
				mesh.vertex[ele.vertex[j] - 1].dofIndex);
// #ifdef __FEMSOLVESYS_DEBUG_LV2
// 				std::cout << "  global " << ele.vertex[i] << " local " << i
// 						  << " and global " << ele.vertex[j] << " local " << j
// 						  << " integ = " << valIntegA << std::endl;
// #endif

		}
	}
	for(int i = 0; i < ele.vertex.size(); i++){
		if(mesh.vertex[ele.vertex[i] - 1].bctype > 0)
			continue;
		double valIntegARH = integARH(ele, mesh, i, prob);
		(this -> rh)[mesh.vertex[ele.vertex[i] - 1].dofIndex] += valIntegARH;

// #ifdef __FEMSOLVESYS_DEBUG_LV2
// 		std::cout << "  global " << ele.vertex[i] << " local " << i
// 				  << " rh" << mesh.vertex[ele.vertex[i] - 1].dofIndex
// 				  << "= " << valIntegARH << std::endl;
// 		std::cout << "  global " << ele.vertex[i] << " local " << i
// 				  << " rh" << mesh.vertex[ele.vertex[i] - 1].dofIndex
// 				  << "= " << (this -> rh)[mesh.vertex[ele.vertex[i] - 1].dofIndex] << std::endl;
// #endif
	}
	return 0;
}

template <typename MyProblem>
int FEMSolvingSystem<MyProblem>::retrive_dof_count_vertex_index(Mesh<MyProblem> &mesh)
{
	int dof(0);
	for(Vertex &iVer:mesh.vertex)
	{
		if(iVer.bctype == 0)
			iVer.dofIndex = dof++;
	}

	return dof;
}

template <typename MyProblem>
int FEMSolvingSystem<MyProblem>::assembleStiff(Mesh<MyProblem> &mesh, MyProblem prob)
{

	this -> dof = retrive_dof_count_vertex_index(mesh);

	calcDetBEOnMesh(mesh);

#ifdef __FEMSOLVESYS_DEBUG
	std::cout << "start forming system, dof = " << this -> dof << std::endl;
#endif
	clock_t t = clock();

	this -> rh = new double [this -> dof];
	memset(this -> rh, 0, (this -> dof) * sizeof(double));
	this -> ma.resize(this -> dof);

	for(std::vector<Element>::iterator it = mesh.element.begin();
		it != mesh.element.end(); it++)
	{
		getStiff(*it, mesh, prob);
	}

	t = clock() - t;

#ifdef __FEMSOLVESYS_DEBUG
	std::cout << "finish forming system, t = "
			  << (double) t / CLOCKS_PER_SEC << "s"
			  << std::endl << std::endl;
#endif	

// #ifdef __FEMSOLVESYS_DEBUG_LV2
// 	for(int i = 0; i < mesh.kidof; i++)
// 	{
// 		std::cout<<" ma[" << i << "]" << std::endl;
// 		for(std::list<maColEle>::iterator it = this -> ma[i].begin();
// 			it != this -> ma[i].end(); it++)
// 			std::cout<<"    row " << it -> row << " value " << it -> value << std::endl;
// 	}
// #endif

	return 0;
}

template <typename MyProblem>
int FEMSolvingSystem<MyProblem>::consoleOutput(Mesh<MyProblem> mesh, MyProblem prob)
{
	std::vector<Vertex>::iterator it;
	int k;
	for(k = 0, it = mesh.vertex.begin(); it != mesh.vertex.end(); it++){
		if(it -> bctype == 0)
		{
			std::cout << it -> x << " " << it -> y << " " << this -> x[it -> dofIndex] << std::endl;
		}
		else
			std::cout << it -> x << " " << it -> y << " " << prob.gd(it -> x, it -> y) << std::endl;
	}

	return 0;
}

template <typename MyProblem>
int FEMSolvingSystem<MyProblem>::fileOutput(Mesh<MyProblem> mesh, MyProblem prob)
{
	std::ofstream fout((prob.parameters.meshFilename + ".output").c_str());

	for(Element iEle:mesh.element){
		for(int iVer:iEle.vertex){
			Vertex &iver = mesh.vertex[iVer - 1];
			if(iver.bctype == 0)
			{
				fout << iver.x << " " << iver.y << " " << this -> x[iver.dofIndex] << std::endl;
			}
			else{
				fout << iver.x << " " << iver.y << " " << prob.gd(iver.x, iver.y) << std::endl;
			}
		}
	}

	return 0;		
}


template <typename MyProblem>
int FEMSolvingSystem<MyProblem>::triOutput(MyProblem prob, Mesh<MyProblem> mesh)
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
		double err = computeError(mesh, prob);
		std::cout << "error = " << err << std::endl;
	}

	return 0;
}

template <typename MyProblem>
double FEMSolvingSystem<MyProblem>::computeError(Mesh<MyProblem> mesh, MyProblem prob)
{
	std::ofstream fout((prob.parameters.meshFilename + ".err").c_str());

	double err = 0;
	for(Element iEle:mesh.element){
		Vertex &v1 = mesh.vertex[iEle.vertex[0] - 1];
		Vertex &v2 = mesh.vertex[iEle.vertex[1] - 1];
		Vertex &v3 = mesh.vertex[iEle.vertex[2] - 1];
		double p1(0), p2(0), p3(0); 
		double r1(0), r2(0), r3(0); 
		if(v1.bctype == 0){
			p1 = this -> x[v1.dofIndex];
			r1 = fabs(prob.trueSol(v1.x, v1.y) - p1);
		}
		if(v2.bctype == 0){
			p2 = this -> x[v2.dofIndex];
			r2 = fabs(prob.trueSol(v2.x, v2.y) - p2);
		}
		if(v3.bctype == 0){
			p3 = this -> x[v3.dofIndex];
			r3 = fabs(prob.trueSol(v3.x, v3.y) - p3);
		}
		
		fout << v1.x << " " << v1.y << " " << r1 << std::endl;
		fout << v2.x << " " << v2.y << " " << r2 << std::endl;
		fout << v3.x << " " << v3.y << " " << r3 << std::endl;

		err += (r1 * r1 + r2 * r2 + r3 * r3) * iEle.detBE / 6.0;
	}
	return sqrt(err);
}

#endif /* TRI_FEMSOLVESYS_H */

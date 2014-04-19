// the problem to be solved here is
//   -\Delta u + u = 3\cos x \sin y, (x,y) \in \Omega = [\frac{\pi}{2}, \frac{3\pi}{2}] \times [0, \pi]
//   u = 0, (x,y) \in \Gamma
// the variational form
//   \int\nolimits_{\Omega} \nabla u \nabla v + uv = \int\nolimits_{\Omega} fv

// here vertex.index is associated with its index in dof, not the index in all the vertices
// thus only vertex with bctype=0 is defined with this value

#ifndef TRI_FEMSOLVESYS_H
#define TRI_FEMSOLVESYS_H

#define __FEMSOLVESYS_DEBUG
// #define __FEMSOLVESYS_DEBUG_LV2

#include "solveSys.h"

template <typename MyProblem>
class FEMSolvingSystem: public BasicSolvingSystem<MyProblem>
{
public:
	double integA(Element &ele, Mesh<MyProblem> mesh, int vi, int vj, MyProblem prob);

	double integARH(Element ele, Mesh<MyProblem> mesh, int vi, MyProblem prob);

	int getStiff(Element& ele, Mesh<MyProblem> &mesh, MyProblem prob);

	int retrive_dof_count_vertex_index(Mesh<MyProblem> &mesh);

	int assembleStiff(Mesh<MyProblem> &mesh, MyProblem prob);

	int consoleOutput(Mesh<MyProblem> mesh);

	int fileOutput(Mesh<MyProblem> mesh, MyProblem prob);

	int triOutput(MyProblem prob, Mesh<MyProblem> mesh);

};


// gradInteg = \int_E \nabla \phi_i \nabla v
// timeInteg = \int_E uv
// nodal basis are used here
template <typename MyProblem>
double FEMSolvingSystem<MyProblem>::integA(Element &ele, Mesh<MyProblem> mesh, int vi, int vj, MyProblem prob)
{
	if(vi > vj)
		std::swap(vi, vj);

	double x1(mesh.vertex[ele.vertex[0] - 1].x), y1(mesh.vertex[ele.vertex[0] - 1].y),
		   x2(mesh.vertex[ele.vertex[1] - 1].x), y2(mesh.vertex[ele.vertex[1] - 1].y),
		   x3(mesh.vertex[ele.vertex[2] - 1].x), y3(mesh.vertex[ele.vertex[2] - 1].y);

	if(ele.detBE == 0){
		ele.detBE = fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
	}

	double gradInteg = 0;
	if(vi == 0){
		switch(vj){
			case 0:{gradInteg = (y2 - y3) * (y2 - y3) + (x3 - x2) * (x3 - x2); break;}
			case 1:{gradInteg = (y2 - y3) * (y3 - y1) + (x3 - x2) * (x1 - x3); break;}
			case 2:{gradInteg = (y2 - y3) * (y1 - y2) + (x3 - x2) * (x2 - x1); break;}
		}
	}
	else if(vi == 1){
		switch(vj){
			case 1:{gradInteg = (y3 - y1) * (y3 - y1) + (x1 - x3) * (x1 - x3); break;}
			case 2:{gradInteg = (y3 - y1) * (y1 - y2) + (x1 - x3) * (x2 - x1); break;}
		}
	}
	else{
		gradInteg = (y1 - y2) * (y1 - y2) + (x2 - x1) * (x2 - x1);
	}

	gradInteg /= ele.detBE * 2;

	double timeInteg = 0;
	if(vi == vj)
		timeInteg = ele.detBE / 12.0;
	else
		timeInteg = ele.detBE / 24.0;

#ifdef __FEMSOLVESYS_DEBUG_LV2
	std::cout << "  gradInteg = " << gradInteg << " timeInteg = " << timeInteg << std::endl;
#endif

	return gradInteg + timeInteg;
}

// calc \int_E f\phi
// \int_E f\phi = |detB_E| / 6 * (\sum_{i=1}^3 f(m_i)\phi(m_i))
// where m_i is the midpoint of each edge
template <typename MyProblem>
double FEMSolvingSystem<MyProblem>::integARH(Element ele, Mesh<MyProblem> mesh, int vi, MyProblem prob)
{
	double x1(mesh.vertex[ele.vertex[vi] - 1].x), y1(mesh.vertex[ele.vertex[vi] - 1].y);

	double a = 0;
	for(int i = 0; i < ele.vertex.size(); i++){
		if (i == vi)
			continue;
		a += prob.f((x1 + mesh.vertex[ele.vertex[i] - 1].x) / 2, (y1 + mesh.vertex[ele.vertex[i] - 1].y) / 2);
	}

	return a * ele.detBE / 12.0;
}

template <typename MyProblem>
int FEMSolvingSystem<MyProblem>::getStiff(Element& ele, Mesh<MyProblem> &mesh, MyProblem prob){

#ifdef __FEMSOLVESYS_DEBUG_LV2
	std::cout << " assemble element " << ele.index << std::endl;
#endif	

	for(int i = 0; i < ele.vertex.size(); i++){
		if(mesh.vertex[ele.vertex[i] - 1].bctype > 0)
			continue;
		for(int j = 0; j < ele.vertex.size(); j++)
		{
			if(mesh.vertex[ele.vertex[j] - 1].bctype > 0)
				continue;
			double valIntegA = integA(ele, mesh, i, j, prob);
			this -> addToMA(valIntegA, mesh.vertex[ele.vertex[i] - 1].index,
				mesh.vertex[ele.vertex[j] - 1].index);

#ifdef __FEMSOLVESYS_DEBUG_LV2
				std::cout << "  global " << ele.vertex[i] << " local " << i
						  << " and global " << ele.vertex[j] << " local " << j
						  << " integ = " << valIntegA << std::endl;
#endif

		}
	}
	for(int i = 0; i < ele.vertex.size(); i++){
		if(mesh.vertex[ele.vertex[i] - 1].bctype > 0)
			continue;
		double valIntegARH = integARH(ele, mesh, i, prob);
		this -> rh[mesh.vertex[ele.vertex[i] - 1].index] += valIntegARH;

#ifdef __FEMSOLVESYS_DEBUG_LV2
		std::cout << "  global " << ele.vertex[i] << " local " << i
				  << " rh" << mesh.vertex[ele.vertex[i] - 1].index
				  << "= " << valIntegARH << std::endl;
#endif
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
			iVer.index = dof++;
	}

	return dof;
}

template <typename MyProblem>
int FEMSolvingSystem<MyProblem>::assembleStiff(Mesh<MyProblem> &mesh, MyProblem prob)
{

	this -> dof = retrive_dof_count_vertex_index(mesh);

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

#ifdef __FEMSOLVESYS_DEBUG_LV2
	for(int i = 0; i < mesh.kidof; i++)
	{
		std::cout<<" ma[" << i << "]" << std::endl;
		for(std::list<maColEle>::iterator it = this -> ma[i].begin();
			it != this -> ma[i].end(); it++)
			std::cout<<"    row " << it -> row << " value " << it -> value << std::endl;
	}
#endif

	return 0;
}

template <typename MyProblem>
int FEMSolvingSystem<MyProblem>::consoleOutput(Mesh<MyProblem> mesh)
{
	std::vector<Vertex>::iterator it;
	int k;
	for(k = 0, it = mesh.vertex.begin(); it != mesh.vertex.end(); it++)
		if(it -> bctype == 0)
		{
			std::cout << it -> x << " " << it -> y << " " << this -> x[k++] << std::endl;
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
				fout << iver.x << " " << iver.y << " " << this -> x[iver.index] << std::endl;
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
		consoleOutput(mesh);
	if(prob.parameters.fprintResults)
		fileOutput(mesh, prob);
	if(prob.parameters.fprintMA)
		this -> fileOutputMA(mesh, prob);
	if(prob.parameters.fprintRH)
		this -> fileOutputRH(mesh, prob);
	if(prob.parameters.fprintTriplet)
		this -> fileOutputTriplet(mesh, prob);

	return 0;
}

#endif /* TRI_FEMSOLVESYS_H */

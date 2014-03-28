#define __FEMSOLVESYS_DEBUG
// #define __FEMSOLVESYS_DEBUG_LV2

#include "FEMSolveSys.h"

// gradInteg = \int_E \nabla \phi_i \nabla v
// timeInteg = \int_E uv
// nodal basis are used here
double FEMSolvingSystem::integA(Element &ele, Mesh mesh, int vi, int vj, MyProblem prob)
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
double FEMSolvingSystem::integARH(Element ele, Mesh mesh, int vi, MyProblem prob)
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

int FEMSolvingSystem::getStiff(Element& ele, Mesh& mesh, MyProblem prob){

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
			addToMA(valIntegA, mesh.vertex[ele.vertex[i] - 1].index - 1,
				mesh.vertex[ele.vertex[j] - 1].index - 1);

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
		rh[mesh.vertex[ele.vertex[i] - 1].index - 1] += valIntegARH;

#ifdef __FEMSOLVESYS_DEBUG_LV2
		std::cout << "  global " << ele.vertex[i] << " local " << i
				  << " rh" << mesh.vertex[ele.vertex[i] - 1].index - 1
				  << "= " << valIntegARH << std::endl;
#endif
	}
	return 0;
}

int FEMSolvingSystem::assembleStiff(Mesh &mesh, MyProblem prob)
{

#ifdef __FEMSOLVESYS_DEBUG
	std::cout << "start forming system, kidof = " << mesh.kidof << std::endl;
#endif
	clock_t t = clock();

	rh = new double [mesh.kidof];
	memset(rh, 0, (mesh.kidof) * sizeof(double));
	ma.resize(mesh.kidof);

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
		for(std::list<maColEle>::iterator it = ma[i].begin();
			it != ma[i].end(); it++)
			std::cout<<"    row " << it -> row << " value " << it -> value << std::endl;
	}
#endif

	return 0;
}

int FEMSolvingSystem::consoleOutput(Mesh mesh)
{
	std::vector<Vertex>::iterator it;
	int k;
	for(k = 0, it = mesh.vertex.begin(); it != mesh.vertex.end(); it++)
		if(it -> bctype == 0)
		{
			std::cout << it -> x << " " << it -> y << " " << x[k++] << std::endl;
		}

	return 0;
}

int FEMSolvingSystem::fileOutput(Mesh mesh, MyProblem prob)
{
	std::ofstream fout((prob.parameters.meshFilename + ".output").c_str());

	std::vector<Vertex>::iterator it;
	int k;
	for(k = 0, it = mesh.vertex.begin(); it != mesh.vertex.end(); it++)
		if(it -> bctype == 0)
		{
			fout << it -> x << " " << it -> y << " " << x[k++] << std::endl;
		}

	return 0;		
}


int FEMSolvingSystem::triOutput(MyProblem prob, Mesh mesh)
{
	if(prob.parameters.printResults)
		consoleOutput(mesh);
	if(prob.parameters.fprintResults)
		fileOutput(mesh, prob);
	if(prob.parameters.fprintMA)
		fileOutputMA(mesh, prob);
	if(prob.parameters.fprintRH)
		fileOutputRH(mesh, prob);
	if(prob.parameters.fprintTriplet)
		fileOutputTriplet(mesh, prob);

	return 0;
}
#include "util.h"

int consoleOutput(Mesh mesh, MySolvingSystem solSys)
{
	std::vector<Vertex>::iterator it;
	int k;
	for(k = 0, it = mesh.vertex.begin(); it != mesh.vertex.end(); it++)
		if(it -> bctype == 0)
		{
			std::cout << it -> x << " " << it -> y << " " << solSys.x[k++] << std::endl;
		}

	return 0;
}

int fileOutput(Mesh mesh, MySolvingSystem solSys)
{
	std::ofstream fout("dat/square.1.output");

	std::vector<Vertex>::iterator it;
	int k;
	for(k = 0, it = mesh.vertex.begin(); it != mesh.vertex.end(); it++)
		if(it -> bctype == 0)
		{
			fout << it -> x << " " << it -> y << " " << solSys.x[k++] << std::endl;
		}

	return 0;		
}

int fileOutputRH(Mesh mesh, MySolvingSystem solSys)
{
	std::ofstream fout("dat/square.1.rh");

	for(int i = 0; i < mesh.kidof; i++)
		fout << solSys.rh[i] << std::endl;

	return 0;		
}

int fileOutputMA(Mesh mesh, MySolvingSystem solSys)
{
	std::ofstream fout("dat/square.1.ma");

	for(int i = 0; i < solSys.ma.size() + 1; i++)
		fout << solSys.Ap[i] << " ";
	fout << std::endl;

	int nnz(0);
	for(std::vector< std::list<maColEle> >::iterator it = solSys.ma.begin();
		it != solSys.ma.end(); it++)
		nnz += it -> size();

	for(int i = 0; i < nnz; i++)
		fout << solSys.Ai[i] << " ";
	fout << std::endl;

	for(int i = 0; i < nnz; i++)
		fout << solSys.Ax[i] << " ";
	fout << std::endl;

	return 0;
}

int fileOutputTriplet(Mesh mesh, MySolvingSystem solSys)
{
	std::ofstream fout("dat/square.1.triplet");

	std::vector< std::list<maColEle> >::iterator it;
	std::list<maColEle>::iterator it1;
	int k;
	for(it = solSys.ma.begin(), k = 1; it != solSys.ma.end(); it++, k++)
		for(it1 = it -> begin(); it1 != it -> end(); it1++)
			fout << ((it1 -> row) + 1) << " " << k << " " << (it1 -> value) << std::endl;

	return 0;
}
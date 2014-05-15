#include "tri_DG.h"

int main(int argc, char const *argv[])
{
	DGProblem prob;
	if(prob.initProblem(argc, argv) == 1)
		return 1;

	Mesh mesh;
	if(mesh.initMesh(prob) == 1)
		return 1;

	DGSolvingSystem solSys;
	solSys.assembleStiff(mesh, prob);
	solSys.solveSparse(mesh, prob);
	solSys.triOutput(prob, mesh);
	return 0;
}
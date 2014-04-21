#include "tri.h"

int main(int argc, char const *argv[])
{
	FEMProblem prob;
	if(prob.initFEMProblem(argc, argv) == 1)
		return 1;

	Mesh<FEMProblem> mesh;
	if(mesh.initMesh(prob) == 1)
		return 1;

	FEMSolvingSystem<FEMProblem> solSys;
	solSys.assembleStiff(mesh, prob);
	solSys.solveSparse(mesh, prob);
	solSys.triOutput(prob, mesh);
	return 0;
}

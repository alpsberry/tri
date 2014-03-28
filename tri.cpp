#include "tri.h"

int main(int argc, char const *argv[])
{
	MyProblem prob;
	if(prob.initProblem(argc, argv) == 1)
		return 1;

	Mesh mesh;
	if(mesh.initMesh(prob) == 1)
		return 1;

	FEMSolvingSystem solSys;
	solSys.assembleStiff(mesh, prob);
	solSys.solveSparse(mesh, prob);
	solSys.triOutput(prob, mesh);
	return 0;
}

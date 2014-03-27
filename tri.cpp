// tri 0.02 Mar. 26
#include "tri.h"

int main(int argc, char const *argv[])
{
	MyProblem prob;
	if(prob.initProblem(argc, argv) == 1)
		return 1;

	Mesh mesh;
	if(mesh.initMesh(prob) == 1)
		return 1;

	MySolvingSystem solSys;
	solSys.assembleStiff(mesh, prob);
	solSys.solveSparse(mesh, prob);

	triOutput(prob, mesh, solSys);
	return 0;
}

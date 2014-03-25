// tri 0.01 Mar. 17
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
	solSys.convertToUMF(mesh);
	solSys.UMFSolve(mesh);

	// consoleOutput(mesh, solSys);
	fileOutput(mesh, solSys);
	fileOutputMA(mesh, solSys);
	fileOutputRH(mesh, solSys);
	fileOutputTriplet(mesh, solSys);
	return 0;
}

#include "tri_DG.h"

int main(int argc, char const *argv[])
{
    DGProblem prob;
    Mesh mesh;
    DGSolvingSystem solSys;
    try {
        prob.initProblem(argc, argv);
        mesh.initMesh(prob);

        solSys.assembleStiff(mesh, prob);
        solSys.solveSparse(mesh, prob);
        solSys.triOutput(prob, mesh);
        
    } catch (std::runtime_error &e) {
        cout << e.what() << endl;
    }

    return 0;
}
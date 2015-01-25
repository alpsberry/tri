#include "tri_DG.h"

int main(int argc, char const *argv[])
{
    try {
        DGProblem prob(argc, argv);
        Mesh mesh(prob);

        DGSolvingSystem solSys;
        solSys.assembleStiff(mesh, prob);
        solSys.solveSparse(mesh, prob);
        solSys.triOutput(mesh, prob);

    } catch (std::runtime_error &e) {
        cout << e.what() << endl;
    }

    return 0;
}
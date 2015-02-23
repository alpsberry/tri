CC = clang++
MPICC = mpic++
CFLAGS = -std=c++11 -Wall
# compile with UMFPACK, SuperLU, and SuperLUDIST
LDFLAGS =  -lumfpack -lamd -lsuitesparseconfig -lcholmod -lcolamd  -framework Accelerate ../SuperLU_4.3/lib/libblas.a ../SuperLU_4.3/lib/libsuperlu_4.3.a ../SuperLU_DIST_3.3/lib/libsuperlu_dist_3.3.a -lmetis -lparmetis
# compile with UMFPACK
# LDFLAGS =  -lumfpack -lamd -lsuitesparseconfig -lcholmod -lcolamd  -framework Accelerate
# compile with SuperLU4.3
# LDFLAGS =  ../SuperLU_4.3/lib/libblas.a ../SuperLU_4.3/lib/libsuperlu_4.3.a

all: tri

release:
	(make CFLAGS="-Wall -O2 -std=c++11" all;)

tri: main.o Mesh.o DGSolvingSystem.o DGProblem.o BasicSolvingSystem.o Problem.o LinearSolver.o UMFPACKSolver.o SuperLUSolver.o
	$(CC) $(CFLAGS) $(LDFLAGS) main.o Mesh.o DGSolvingSystem.o DGProblem.o BasicSolvingSystem.o LinearSolver.o UMFPACKSolver.o SuperLUSolver.o Problem.o -o tri

trimpi: maintrimpi.o Mesh.o DGSolvingSystemMPI.o DGSolvingSystem.o DGProblem.o BasicSolvingSystem.o Problem.o LinearSolver.o UMFPACKSolver.o SuperLUSolver.o SuperLUDISTSolver.o
	$(MPICC) $(CFLAGS) $(LDFLAGS) maintrimpi.o DGSolvingSystemMPI.o Mesh.o DGSolvingSystem.o DGProblem.o BasicSolvingSystem.o LinearSolver.o UMFPACKSolver.o SuperLUSolver.o SuperLUDISTSolver.o Problem.o -o trimpi

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

maintrimpi.o: maintrimpi.cpp
	$(MPICC) $(CFLAGS) -c maintrimpi.cpp

Mesh.o: Mesh.cpp
	$(CC) $(CFLAGS) -c Mesh.cpp

DGSolvingSystem.o: DGSolvingSystem.cpp
	$(CC) $(CFLAGS) -c DGSolvingSystem.cpp

Problem.o: Problem.cpp
	$(CC) $(CFLAGS) -c Problem.cpp

DGProblem.o: DGProblem.cpp
	$(CC) $(CFLAGS) -c DGProblem.cpp

BasicSolvingSystem.o: BasicSolvingSystem.cpp
	$(CC) $(CFLAGS) -c BasicSolvingSystem.cpp

DGSolvingSystemMPI.o: DGSolvingSystemMPI.cpp
	$(MPICC) $(CFLAGS) -c DGSolvingSystemMPI.cpp

LinearSolver.o: LinearSolver.cpp
	$(CC) $(CFLAGS) -c LinearSolver.cpp

UMFPACKSolver.o: UMFPACKSolver.cpp
	$(CC) $(CFLAGS) -c UMFPACKSolver.cpp	

SuperLUSolver.o: SuperLUSolver.cpp
	$(CC) $(CFLAGS) -c SuperLUSolver.cpp

SuperLUDISTSolver.o: SuperLUDISTSolver.cpp
	$(CC) $(CFLAGS) -c SuperLUDISTSolver.cpp


clean:
	rm -rf *o tri


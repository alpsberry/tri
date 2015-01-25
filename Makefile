CC = clang++
CFLAGS = -std=c++11 -Wall
# compile with both UMFPACK and SuperLU 
LDFLAGS =  -lumfpack -lamd -lsuitesparseconfig -lcholmod -lcolamd  -framework Accelerate ../SuperLU_4.3/lib/libblas.a ../SuperLU_4.3/lib/libsuperlu_4.3.a
# compile with UMFPACK
# LDFLAGS =  -lumfpack -lamd -lsuitesparseconfig -lcholmod -lcolamd  -framework Accelerate
# compile with SuperLU4.3
# LDFLAGS =  ../SuperLU_4.3/lib/libblas.a ../SuperLU_4.3/lib/libsuperlu_4.3.a

all: tri

release:
	(make CFLAGS="-Wall -O2 -std=c++11" all;)

tri: main.o Mesh.o DGSolvingSystem.o DGProblem.o BasicSolvingSystem.o Problem.o
	$(CC) $(CFLAGS) $(LDFLAGS) main.o Mesh.o DGSolvingSystem.o DGProblem.o BasicSolvingSystem.o Problem.o -o tri

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

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

clean:
	rm -rf *o tri


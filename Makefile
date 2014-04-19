CC = clang++
#CFLAGS = -Wall -O2
CFLAGS = -std=c++11 -Wall
# compile with both UMFPACK and SuperLU 
LDFLAGS =  -lumfpack -lamd -lsuitesparseconfig -lcholmod -lcolamd  -framework Accelerate ../SuperLU_4.3/lib/libblas.a ../SuperLU_4.3/lib/libsuperlu_4.3.a
# compile with UMFPACK
# LDFLAGS =  -lumfpack -lamd -lsuitesparseconfig -lcholmod -lcolamd  -framework Accelerate
# compile with SuperLU4.3
# LDFLAGS =  ../SuperLU_4.3/lib/libblas.a ../SuperLU_4.3/lib/libsuperlu_4.3.a

TRI_TARGET = tri
TRIDG_TARGET = tri_DG

TRI_SRCS = tri.cpp FEMProblem.cpp
TRIDG_SRCS = tri_DG.cpp DGProblem.cpp

TRI_OBJS = ${TRI_SRCS:.cpp=.o}
TRIDG_OBJS = ${TRIDG_SRCS:.cpp=.o}

CLEANFILES = 	$(TRI_OBJS) $(TRI_TARGET) $(TRIDG_OBJS) $(TRIDG_TARGET) 

all: $(TRI_TARGET) $(TRIDG_TARGET)

release:
	(make CFLAGS="-Wall -O2 -std=c++11" all;)

$(TRI_TARGET): $(TRI_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(INCLDIRS)

$(TRIDG_TARGET): $(TRIDG_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(INCLDIRS)


%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLDIRS) -c $<

clean:
	rm -f $(CLEANFILES)

tri.o: tri.h FEMSolveSys.h solveSys.h mesh.h problem.h FEMProblem.h
tri_DG.o: tri_DG.h DGSolveSys.h solveSys.h mesh.h problem.h DGProblem.h

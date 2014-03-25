CC = g++
#CFLAGS = -Wall -O2
# CFLAGS = -g -Wall
CFLAGS = 
LDFLAGS =  -lumfpack -lamd -lsuitesparseconfig -lcholmod -lcolamd  -framework Accelerate
INCLDIRS = 

TRI_TARGET = tri

TRI_SRCS = tri.cpp solveSys.cpp problem.cpp mesh.cpp util.cpp

TRI_OBJS = ${TRI_SRCS:.cpp=.o}

CLEANFILES = 	$(TRI_OBJS) $(TRI_TARGET) 

all: $(TRI_TARGET)

release:
	(make CFLAGS="-Wall -O2" all;)

$(TRI_TARGET): $(TRI_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(INCLDIRS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLDIRS) -c $<

clean:
	rm -f $(CLEANFILES)

# mesh.o: mesh.h
# problem.o: mesh.h problem.h solveSys.h
# solveSys.o: mesh.h problem.h solveSys.h
tri.o: tri.h mesh.h problem.h solveSys.h util.h

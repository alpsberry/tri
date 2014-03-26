#ifndef TRI_UTIL_H
#define TRI_UTIL_H

#include "solveSys.h"

int consoleOutput(Mesh mesh, MySolvingSystem solSys);

int fileOutput(Mesh mesh, MySolvingSystem solSys);

int fileOutputRH(Mesh mesh, MySolvingSystem solSys);

int fileOutputMA(Mesh mesh, MySolvingSystem solSys);

int fileOutputTriplet(Mesh mesh, MySolvingSystem solSys);

#endif /* TRI_UTIL_H */
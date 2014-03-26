#ifndef TRI_PROBLEM_H
#define TRI_PROBLEM_H
#include <cmath>
#include <vector>

struct Vertex{
	int index;
	double x, y;
	int bctype;
};

struct Element;

struct Edge{
	int index;
	std::vector<int> vertex;
	std::vector<int> neighborElement;
	int bctype;
};

struct Element{
	double detBE;
	int index;
	std::vector<int> vertex;
	std::vector<int> edge;
};
// detBE is the absolute value of the determinant of B_E
// The mapping F_E from the reference element to any element E is
// F([x, y]^T) = B_E[x, y]^T + b_E
// where B_E = [x2 - x1, x3 - x1; y2 - y1. y3 - y1] and b_E = [x1, y1]^T

class Region{
public:
	std::vector<Element> element;
	std::vector<Edge> edge;
	std::vector<Vertex> vertex;
};

class Problem
{
public:
	// Region region;
	int dimension;
	virtual double f(double x, double y){
		return 0;
	}
	virtual double gd(double x, double y){
		return 0;
	}
};

// the problem to be solved here is
// -\Delta u + u = 3\cos x \sin y, (x,y) \in \Omega = [\frac{\pi}{2}, \frac{3\pi}{2}] \times [0, \pi]
// u = 0, (x,y) \in \Gamma
class MyProblem: public Problem
{
public:
	double f(double x, double y){
		return 3 * cos(x) * sin(y);
	}
	int initProblem(int argc, char const *argv[])
	{
		dimension = 2;
		return 0;
	}
};

#endif

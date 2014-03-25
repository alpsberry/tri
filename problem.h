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
	int initProblem(int argc, char const *argv[])
	{
		dimension = 2;
		return 0;
	}
};

class MyProblem: public Problem
{
public:
	double f(double x, double y){
		return 3 * cos(x) * sin(y);
	}
};

#endif

#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

const string filename = "square";

// const int N = 10;
int N = 10;
const int M = 2;

const double pi = 3.14159;
const double Width = pi;
const double Height = pi;
const double Left = 0.5 * pi;
const double Bottom = 0;

// const double Width = 1;
// const double Height = 1;
// const double Left = 0;
// const double Bottom = 0;

double hx;
double hy;

void outputPoly()
{
	ofstream fout((filename + ".poly").c_str());
	fout << "4  2  0  1" << endl
	     << "   1\t" << Left << "\t" << Bottom << "\t" << "1" << endl
	     << "   2\t" << Left + Width << "\t" << Bottom << "\t" << "1" << endl
	     << "   3\t" << Left + Width << "\t" << Bottom + Height << "\t" << "1" << endl
	     << "   4\t" << Left << "\t" << Bottom + Height << "\t" << "1" << endl;
	fout << "4  1" << endl
	     << "   1\t1\t2\t1" << endl
	     << "   2\t2\t3\t1" << endl
	     << "   3\t3\t4\t1" << endl
	     << "   4\t4\t1\t1" << endl;
	fout << "0" << endl;

	ofstream fout1((filename + ".1.poly").c_str());
	fout1 << "4  2  0  1" << endl
	     << "   1\t" << Left << "\t" << Bottom << "\t" << "1" << endl
	     << "   2\t" << Left + Width << "\t" << Bottom << "\t" << "1" << endl
	     << "   3\t" << Left + Width << "\t" << Bottom + Height << "\t" << "1" << endl
	     << "   4\t" << Left << "\t" << Bottom + Height << "\t" << "1" << endl;
	fout1 << "4  1" << endl
	     << "   1\t1\t2\t1" << endl
	     << "   2\t2\t3\t1" << endl
	     << "   3\t3\t4\t1" << endl
	     << "   4\t4\t1\t1" << endl;
	fout1 << "0" << endl;
}

void outputNode()
{
	ofstream fout((filename + ".1.node").c_str());
	fout << N * N << "  " << "2  0  1" << endl;

	for(int i = 1; i <= N; ++i)
		fout << "   " << i << "\t" << Left + (i - 1) * hx << "\t" << Bottom << "\t1" << endl;
	for(int i = 2; i <= N - 1; ++i)
	{
		fout << "   " << (i - 1) * N + 1 << "\t" << Left << "\t" << Bottom + (i - 1) * hy << "\t1" << endl;
		for(int j = 2; j <= N-1; ++j)
			fout << "   " << (i - 1) * N + j << "\t" << Left + (j - 1) * hx << "\t" << Bottom + (i - 1) * hy << "\t0" << endl;
		fout << "   " << i * N << "\t" << Left + Width << "\t" << Bottom + (i - 1) * hy << "\t1" << endl;

	}
	for(int i = 1; i <= N; ++i)
		fout << "   " << (N - 1) * N + i << "\t" << Left + (i - 1) * hx << "\t" << Bottom + Height << "\t1" << endl;

}

void outputElement()
{
	ofstream fout((filename + ".1.ele").c_str());
	fout << 2 * (N - 1) * (N - 1) << "  " << "3  0" << endl;
	for(int i = 1; i <= N - 1; ++i)
		for(int j = 1; j <= N-1; ++j)
		{
			fout << "   " << (2 * N - 2) * (i - 1) + 2 * j - 1 << "\t" << (i - 1) * N + j << "\t" << i * N + j + 1 << "\t" << i * N + j << endl;
			fout << "   " << (2 * N - 2) * (i - 1) + 2 * j << "\t" << (i - 1) * N + j << "\t" << (i - 1) * N + j + 1 << "\t" << i * N + j + 1 << endl;
		}

}

void outputEdge()
{
	ofstream fout((filename + ".1.edge").c_str());
	fout << (N - 1) * (3 * N - 1) << "  1" << endl;

	// row
	for(int j = 1; j <= N-1; j++)
		fout << "   " << j << "\t" << j << "\t" << j + 1 << "\t" << "1" << endl;
	for(int i = 2; i <= N - 1; i++)
		for(int j = 1; j <= N - 1; j++)
			fout << "   " << (i - 1) * (N - 1) + j << "\t" << (i - 1) * N + j << "\t" << (i - 1) * N + j + 1 << "\t" << "0" << endl;
	for(int j = 1; j <= N - 1; j++)
		fout << "   " << (N - 1) * (N - 1) + j << "\t" << (N - 1) * N + j << "\t" << (N - 1) * N + j + 1 << "\t" << "1" << endl;

	// col
	for(int i = 1; i <= N - 1; i++)
		fout << "   " << N * (N - 1) + i << "\t" << (i - 1) * N + 1 << "\t" << i * N + 1 << "\t" << "1" << endl;
	for(int j = 2; j <= N - 1; j++)
		for(int i = 1; i <= N - 1; i++)
			fout << "   " << N * (N - 1) + i + (j - 1) * (N - 1) << "\t" << (i - 1) * N + j << "\t" << i * N + j << "\t" << "0" << endl;
	for(int i = 1; i <= N - 1; i++)
		fout << "   " << N * (N - 1) + i + (N - 1) * (N - 1) << "\t" << (i - 1) * N + N << "\t" << i * N + N << "\t" << "1" << endl;


	// diag
	for(int i = 1; i <= N-1; i++)
		for(int j = 1; j <= N-1; j++)
			fout << "   " << 2 * N * (N - 1) + (i - 1) * (N - 1) + j << "\t" << (i - 1) * N + j << "\t" << i * N + j + 1 << "\t0" << endl;

}

void outputRefineLeft()
{
	ofstream fout((filename + ".1.ref0").c_str());

	// node
	const int originNodeSize = N * N;
	fout << 4 * N - 3 << endl;

	fout << "   " << originNodeSize + 1 << "\t" << Left + 0.5 * hx << "\t" << Bottom << "\t1" << endl;
	fout << "   " << originNodeSize + 2 << "\t" << Left << "\t" << Bottom + 0.5 * hy << "\t1" << endl;
	fout << "   " << originNodeSize + 3 << "\t" << Left + 0.5 * hx << "\t" << Bottom + 0.5 * hy << "\t0" << endl;
	fout << "   " << originNodeSize + 4 << "\t" << Left + hx << "\t" << Bottom + 0.5 * hy << "\t0" << endl;
	for(int i = 2; i <= N - 1; i++)
	{
		fout << "   " << originNodeSize + 4 * (i - 1) + 1 << "\t" << Left + 0.5 * hx << "\t" << Bottom + (i - 1) * hy << "\t0" << endl;
		fout << "   " << originNodeSize + 4 * (i - 1) + 2 << "\t" << Left << "\t" << Bottom + (i - 1) * hy + 0.5 * hy << "\t1" << endl;
		fout << "   " << originNodeSize + 4 * (i - 1) + 3 << "\t" << Left + 0.5 * hx << "\t" << Bottom + (i - 1) * hy + 0.5 * hy << "\t0" << endl;
		fout << "   " << originNodeSize + 4 * (i - 1) + 4 << "\t" << Left + hx << "\t" << Bottom + (i - 1) * hy + 0.5 * hy << "\t0" << endl;
	}
	fout << "   " << originNodeSize + 4 * (N - 1) + 1 << "\t" << Left + 0.5 * hx << "\t" << Bottom + Height << "\t1" << endl;

	// edge
	const int originEdgeSize = (N - 1) * (3 * N - 1);
	fout << 14 * N - 12 << endl;
	// edge row
	for(int j = 1; j <= N-1; j++)
	{
		fout << "\t" << (j - 1) * (N - 1) + 1 << endl
		     << "\t\t" << originEdgeSize + (j - 1) * 4 + 1 << "\t" << (j - 1) * N + 1 << "\t" << originNodeSize + (j - 1) * 4 + 1 << endl
		     << "\t\t" << originEdgeSize + (j - 1) * 4 + 2 << "\t" << originNodeSize + (j - 1) * 4 + 1 << "\t" << (j - 1) * N + 2 << endl;
		fout << "\t" << "0" << endl
		     << "\t\t" << originEdgeSize + (j - 1) * 4 + 3 << "\t" << originNodeSize + (j - 1) * 4 + 2 <<"\t" << originNodeSize + (j - 1) * 4 + 3 << endl
		     << "\t\t" << originEdgeSize + (j - 1) * 4 + 4 << "\t" << originNodeSize + (j - 1) * 4 + 3 <<"\t" << originNodeSize + (j - 1) * 4 + 4 << endl;
	}
	fout << "\t" << (N - 1) * (N - 1) + 1 << endl
	     << "\t\t" << originEdgeSize + (N - 1) * 4 + 1 << "\t" << (N - 1) * N + 1 << "\t" << originNodeSize + (N - 1) * 4 + 1 << endl
	     << "\t\t" << originEdgeSize + (N - 1) * 4 + 2 << "\t" << originNodeSize + (N - 1) * 4 + 1 << "\t" << (N - 1) * N + 2 << endl;

	// edge col
	const int newRowSize = 4 * N - 2;
	for(int i = 1; i <= N - 1; i++)
		fout << "\t" << N * (N - 1) + i << endl
	         << "\t\t" << originEdgeSize + newRowSize + 2 * (i - 1) + 1 << "\t" << (i - 1) * N + 1 << "\t" << originNodeSize + 4 * (i - 1) + 2 << endl
	         << "\t\t" << originEdgeSize + newRowSize + 2 * (i - 1) + 2 << "\t" << originNodeSize + 4 * (i - 1) + 2 << "\t" << i * N + 1 << endl;
	fout << "\t0" << endl;
	for(int i = 1; i <= N - 1; i++)
		fout << "\t\t" << originEdgeSize + newRowSize + (N - 1) * 2 + 2 * (i - 1) + 1 << "\t" << originNodeSize + 4 * (i - 1) + 1 << "\t" << originNodeSize + 4 * (i - 1) + 3 << endl
	         << "\t\t" << originEdgeSize + newRowSize + (N - 1) * 2 + 2 * (i - 1) + 2 << "\t" << originNodeSize + 4 * (i - 1) + 3 << "\t" << originNodeSize + 4 * i + 1 << endl;
	for(int i = 1; i <= N - 1; i++)
		fout << "\t" << N * (N - 1) + (N - 1) + i << endl
	         << "\t\t" << originEdgeSize + newRowSize + (N - 1) * 4 + 2 * (i - 1) + 1 << "\t" << (i - 1) * N + 2 << "\t" << originNodeSize + 4 * (i - 1) + 4 << endl
	         << "\t\t" << originEdgeSize + newRowSize + (N - 1) * 4 + 2 * (i - 1) + 2 << "\t" << originNodeSize + 4 * (i - 1) + 4 << "\t" << i * N + 2 << endl;

	// edge diag
	const int newColSize = 6 * (N - 1);
	for(int i = 1; i <= N - 1; i++)
		fout << "\t" << 2 * N * (N - 1) + (i - 1) * (N - 1) + 1 << endl
	         << "\t\t" << originEdgeSize + newRowSize + newColSize + 4 * (i - 1) + 1 << "\t" << (i - 1) * N + 1 << "\t" << originNodeSize + 4 * (i - 1) + 3 << endl
	         << "\t" << "0" << endl
	         << "\t\t" << originEdgeSize + newRowSize + newColSize + 4 * (i - 1) + 2 << "\t" << originNodeSize + 4 * (i - 1) + 1 << "\t" << originNodeSize + 4 * (i - 1) + 4 << endl
	         << "\t" << "0" << endl
	         << "\t\t" << originEdgeSize + newRowSize + newColSize + 4 * (i - 1) + 3 << "\t" << originNodeSize + 4 * (i - 1) + 2 << "\t" << originNodeSize + 4 * i + 1 << endl
	         << "\t" << 2 * N * (N - 1) + (i - 1) * (N - 1) + 1 << endl
	         << "\t\t" << originEdgeSize + newRowSize + newColSize + 4 * (i - 1) + 4 << "\t" << originNodeSize + 4 * (i - 1) + 3 << "\t" << i * N + 2 << endl;



	// element
	const int originElementSize = 2 * (N - 1) * (N - 1);
	fout << 8 * (N - 1) << endl;
	for(int j = 1; j <= N - 1; j++)
	{
		fout << "\t" << 2 * (N - 1) * (j - 1) + 1<< endl
		     << "\t\t" << originElementSize + 8 * (j - 1) + 1 << "\t" << (j - 1) * N + 1 << "\t" << originNodeSize + 4 * (j - 1) + 3 << "\t" << originNodeSize + 4 * (j - 1) + 2 << endl
		     << "\t" << 2 * (N - 1) * (j - 1) + 2<< endl
		     << "\t\t" << originElementSize + 8 * (j - 1) + 2 << "\t" << (j - 1) * N + 1 << "\t" << originNodeSize + 4 * (j - 1) + 1 << "\t" << originNodeSize + 4 * (j - 1) + 3 << endl
		     << "\t\t" << originElementSize + 8 * (j - 1) + 3 << "\t" << originNodeSize + 4 * (j - 1) + 1 << "\t" << originNodeSize + 4 * (j - 1) + 4 << "\t" << originNodeSize + 4 * (j - 1) + 3 << endl
		     << "\t\t" << originElementSize + 8 * (j - 1) + 4 << "\t" << originNodeSize + 4 * (j - 1) + 1 << "\t" << (j - 1) * N + 2 << "\t" << originNodeSize + 4 * (j - 1) + 4 << endl
		     << "\t" << 2 * (N - 1) * (j - 1) + 1<< endl
		     << "\t\t" << originElementSize + 8 * (j - 1) + 5 << "\t" << originNodeSize + 4 * (j - 1) + 2 << "\t" << originNodeSize + 4 * j + 1 << "\t" << j * N + 1 << endl
		     << "\t\t" << originElementSize + 8 * (j - 1) + 6 << "\t" << originNodeSize + 4 * (j - 1) + 2 << "\t" << originNodeSize + 4 * (j - 1) + 3 << "\t" << originNodeSize + 4 * j + 1 << endl
		     << "\t\t" << originElementSize + 8 * (j - 1) + 7 << "\t" << originNodeSize + 4 * (j - 1) + 3 << "\t" << j * N + 2 << "\t" << originNodeSize + 4 * j + 1 << endl
		     << "\t" << 2 * (N - 1) * (j - 1) + 2<< endl
		     << "\t\t" << originElementSize + 8 * (j - 1) + 8 << "\t" << originNodeSize + 4 * (j - 1) + 3 << "\t" << originNodeSize + 4 * (j - 1) + 4 << "\t" << j * N + 2 << endl;
	}

}

void outputRefineCentral0()
{
	ofstream fout((filename + ".1.ref0").c_str());

	const double startx = N / 2;
	const double starty = N / 2;

	// node
	const int originNodeSize = N * N;
	fout << 5 << endl;
	fout << "   " << originNodeSize + 1 << "\t" << Left + (startx - 1) * hx + 0.5 * hx << "\t" << Bottom + (starty - 1) * hy << "\t0"<< endl
	     << "   " << originNodeSize + 2 << "\t" << Left + (startx - 1) * hx << "\t" << Bottom + (starty - 1) * hy + 0.5 * hy << "\t0"<< endl
	     << "   " << originNodeSize + 3 << "\t" << Left + (startx - 1) * hx + 0.5 * hx << "\t" << Bottom + (starty - 1) * hy + 0.5 * hy << "\t0"<< endl
	     << "   " << originNodeSize + 4 << "\t" << Left + startx * hx << "\t" << Bottom + (starty - 1) * hy + 0.5 * hy << "\t0"<< endl
	     << "   " << originNodeSize + 5 << "\t" << Left + (startx - 1) * hx + 0.5 * hx << "\t" << Bottom + starty * hy << "\t0"<< endl;

	// edge
	const int originEdgeSize = (N - 1) * (3 * N - 1);
	fout << 16 << endl;

	// edge row
	fout << "\t" << (starty - 1) * (N - 1) + startx << endl
	     << "\t\t" << originEdgeSize + 1 << "\t" << (starty - 1) * N + startx << "\t" << originNodeSize + 1 << endl
	     << "\t\t" << originEdgeSize + 2 << "\t" << originNodeSize + 1 << "\t" << (starty - 1) * N + startx + 1 << endl
	     << "\t0" << endl
	     << "\t\t" << originEdgeSize + 3 << "\t" << originNodeSize + 2 << "\t" << originNodeSize + 3 << endl
	     << "\t\t" << originEdgeSize + 4 << "\t" << originNodeSize + 3 << "\t" << originNodeSize + 4 << endl
	     << "\t" << starty * (N - 1) + startx << endl
	     << "\t\t" << originEdgeSize + 5 << "\t" << starty * N + startx << "\t" << originNodeSize + 5 << endl
	     << "\t\t" << originEdgeSize + 6 << "\t" << originNodeSize + 5 << "\t" << starty * N + startx + 1 << endl;

	// edge col
	const int newRowSize = 6;
	fout << "\t" << N * (N - 1) + (startx - 1) * (N - 1) + starty << endl
	     << "\t\t" << originEdgeSize + newRowSize + 1 << "\t" << (starty - 1) * N + startx << "\t" << originNodeSize + 2 << endl
	     << "\t\t" << originEdgeSize + newRowSize + 2 << "\t" << originNodeSize + 2 << "\t" << starty * N + startx << endl
	     << "\t0" << endl
	     << "\t\t" << originEdgeSize + newRowSize + 3 << "\t" << originNodeSize + 1 << "\t" << originNodeSize + 3 << endl
	     << "\t\t" << originEdgeSize + newRowSize + 4 << "\t" << originNodeSize + 3 << "\t" << originNodeSize + 5 << endl
	     << "\t" << N * (N - 1) + startx * (N - 1) + starty << endl
	     << "\t\t" << originEdgeSize + newRowSize + 5 << "\t" << (starty - 1) * N + startx + 1 << "\t" << originNodeSize + 4 << endl
	     << "\t\t" << originEdgeSize + newRowSize + 6 << "\t" << originNodeSize + 4 << "\t" << starty * N + startx + 1 << endl;	

	// edge diag
	const int newColSize = 6;
	fout << "\t" << 2 * N * (N - 1) + (starty - 1) * (N - 1) + startx << endl
	     << "\t\t" << originEdgeSize + newRowSize + newColSize + 1 << "\t" << (starty - 1) * N + startx << "\t" << originNodeSize + 3 << endl
	     << "\t0" << endl
	     << "\t\t" << originEdgeSize + newRowSize + newColSize + 2 << "\t" << originNodeSize + 1 << "\t" << originNodeSize + 4 << endl
	     << "\t\t" << originEdgeSize + newRowSize + newColSize + 3 << "\t" << originNodeSize + 2 << "\t" << originNodeSize + 5 << endl
	     << "\t" << 2 * N * (N - 1) + (starty - 1) * (N - 1) + startx << endl
	     << "\t\t" << originEdgeSize + newRowSize + newColSize + 4 << "\t" << originNodeSize + 3 << "\t" << starty * N + startx + 1 << endl;


	// element
	const int originElementSize = 2 * (N - 1) * (N - 1);
	fout << 8 << endl
	     << "\t" << (starty - 1) * 2 * (N - 1) + 2 * (startx - 1) + 1 << endl
	     << "\t\t" << originElementSize + 1 << "\t" << (starty - 1) * N + startx << "\t" << originNodeSize + 3 << "\t" << originNodeSize + 2 << endl
	     << "\t" << (starty - 1) * 2 * (N - 1) + 2 * (startx - 1) + 2<< endl
	     << "\t\t" << originElementSize + 2 << "\t" << (starty - 1) * N + startx << "\t" << originNodeSize + 1 << "\t" << originNodeSize + 3 << endl
	     << "\t\t" << originElementSize + 3 << "\t" << originNodeSize + 1 << "\t" << originNodeSize + 4 << "\t" << originNodeSize + 3 << endl
	     << "\t\t" << originElementSize + 4 << "\t" << originNodeSize + 1 << "\t" << (starty - 1) * N + startx + 1 << "\t" << originNodeSize + 4 << endl
	     << "\t" << (starty - 1) * 2 * (N - 1) + 2 * (startx - 1) + 1<< endl
	     << "\t\t" << originElementSize + 5 << "\t" << originNodeSize + 2 << "\t" << originNodeSize + 5 << "\t" << starty * N + startx << endl
	     << "\t\t" << originElementSize + 6 << "\t" << originNodeSize + 2 << "\t" << originNodeSize + 3 << "\t" << originNodeSize + 5 << endl
	     << "\t\t" << originElementSize + 7 << "\t" << originNodeSize + 3 << "\t" << starty * N + startx + 1 << "\t" << originNodeSize + 5 << endl
	     << "\t" << (starty - 1) * 2 * (N - 1) + 2 * (startx - 1) + 2<< endl
	     << "\t\t" << originElementSize + 8 << "\t" << originNodeSize + 3 << "\t" << originNodeSize + 4 << "\t" << starty * N + startx + 1 << endl;
}

void outputRefineCentral1()
{
	ofstream fout((filename + ".1.ref1").c_str());

	const double startx = N / 2;
	const double starty = N / 2;	
	
	// node
	const int originNodeSize = N * N + 5;
	fout << 1 << endl;
	fout << "   " << originNodeSize + 1 << "\t" << Left + (startx - 1) * hx + 0.25 * hx << "\t" << Bottom + (starty - 1) * hy + 0.25 * hy << "\t0"<< endl;

	// edge
	const int originEdgeSize = (N - 1) * (3 * N - 1) + 16;
	fout << 2 << endl
	     << "\t" << originEdgeSize - 3 << endl
	     << "\t\t" << originEdgeSize + 1 << "\t" << (starty - 1) * N + startx << "\t" << originNodeSize + 1 << endl
	     << "\t\t" << originEdgeSize + 2 << "\t" << originNodeSize + 1 << "\t" << originNodeSize - 2 << endl;

	// element
	const int originElementSize = 2 * (N - 1) * (N - 1) + 8;
	fout << 2 << endl
	     << "\t" << originElementSize - 6 << endl
	     << "\t\t" << originElementSize + 1 << "\t" << (starty - 1) * N + startx << "\t" << originNodeSize - 4 << "\t" << originNodeSize + 1 << endl
	     << "\t\t" << originElementSize + 2 << "\t" << originNodeSize - 2 << "\t" << originNodeSize - 4 << "\t" << originNodeSize + 1 << endl;

}

int main(int argc, char **argv)
{	
	if( argc > 2 )
		return 1;

	N = atoi(argv[1]);      /* Number of processors requested */
	cout << "N = " << N << endl;

	if(N < 3){
		cout << "N cannot be less than 3" << endl;
		return 1;
	}

	hx = Width / (N - 1);
	hy = Height / (N - 1);

	outputPoly();
	outputNode();
	outputElement();
	outputEdge();
	outputRefineLeft();
	// outputRefineCentral0();
	// outputRefineCentral1();
}
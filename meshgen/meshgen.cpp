#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>

using namespace std;
using std::vector;

const string filename = "square";

// const int N = 10;
int N = 10;
int M = 5;
int Z = 3;

const double pi = 3.14159;
const double Width = pi;
const double Height = pi;
const double Left = 0.5 * pi;
const double Bottom = 0;

// const double Width = 1;
// const double Height = 1;
// const double Left = 0;
// const double Bottom = 0;

double hx, hy, mx, my, zx, zy;

class Node {
public:
    int index;
    double x;
    double y;
    int bctype;
    Node(int i, double xx, double yy, int bct): index(i), x(xx), y(yy), bctype(bct) {}
};

class Edge {
public:
    int index;
    int v1, v2;
    int bctype;
    int parent;
    Edge(int i, int v11, int v22, int bct, int pa): index(i), v1(v11), v2(v22), bctype(bct), parent(pa) {}
};

class Element {
public:
    int index;
    int v1, v2, v3;
    int parent;
    Element(int i, int v11, int v22, int v33, int pa): index(i), v1(v11), v2(v22), v3(v33), parent(pa) {}
};

vector<Node> vecNode;
vector<Edge> vecEdge;
vector<Element> vecElement;

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

void genNode()
{
    Node *pNode;

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            if (j == 1 || j == N || i == 1 || i == N)
                pNode = new Node((i - 1) * N + j, Left + (j - 1) * hx, Bottom + (i - 1) * hy, 1);
            else
                pNode = new Node((i - 1) * N + j, Left + (j - 1) * hx, Bottom + (i - 1) * hy, 0);
            vecNode.push_back(*pNode);
        }
    }
}

void outputNode()
{
    genNode();
    ofstream fout((filename + ".1.node").c_str());
    fout << N *N << "  " << "2  0  1" << endl;
    for (int i = 0; i < N * N; i++)
        fout << "   " << vecNode[i].index << "\t" << vecNode[i].x << "\t" << vecNode[i].y << "\t" << vecNode[i].bctype << endl;

}

void outputElement()
{
    ofstream fout((filename + ".1.ele").c_str());
    fout << 2 * (N - 1) * (N - 1) << "  " << "3  0" << endl;
    for (int i = 1; i <= N - 1; ++i)
        for (int j = 1; j <= N - 1; ++j) {
            fout << "   " << (2 * N - 2) * (i - 1) + 2 * j - 1 << "\t" << (i - 1) * N + j << "\t" << i *N + j + 1 << "\t" << i *N + j << endl;
            fout << "   " << (2 * N - 2) * (i - 1) + 2 * j << "\t" << (i - 1) * N + j << "\t" << (i - 1) * N + j + 1 << "\t" << i *N + j + 1 << endl;
        }

}

void genEdge()
{
    // row
    for (int i = 1; i <= N; i++)
        for (int j = 1; j <= N - 1; j++) {
            if (i == 1 || i == N)
                vecEdge.push_back(Edge((i - 1) * (N - 1) + j, (i - 1) * N + j, (i - 1) * N + j + 1, 1, 0));
            else
                vecEdge.push_back(Edge((i - 1) * (N - 1) + j, (i - 1) * N + j, (i - 1) * N + j + 1, 0, 0));
        }

    // col
    for (int j = 1; j <= N; j++)
        for (int i = 1; i <= N - 1; i++)
            if (j == 1 || j == N)
                vecEdge.push_back(Edge(N * (N - 1) + i + (j - 1) * (N - 1), (i - 1) * N + j, i * N + j, 1, 0));
            else
                vecEdge.push_back(Edge(N * (N - 1) + i + (j - 1) * (N - 1), (i - 1) * N + j, i * N + j, 0, 0));

    // diag
    for (int i = 1; i <= N - 1; i++)
        for (int j = 1; j <= N - 1; j++)
            vecEdge.push_back(Edge(2 * N * (N - 1) + (i - 1) * (N - 1) + j, (i - 1) * N + j, i * N + j + 1, 0, 0));

}
void outputEdge()
{
    genEdge();

    ofstream fout((filename + ".1.edge").c_str());
    fout << (N - 1) * (3 * N - 1) << "  1" << endl;
    for (int i = 0; i < (N - 1) * (3 * N - 1); i++)
        fout << "   " << vecEdge[i].index << "\t" << vecEdge[i].v1 << "\t" << vecEdge[i].v2 << "\t" << vecEdge[i].bctype << endl;

}

void outputRefineRound();
void outputRefineRound1();

int main(int argc, char **argv)
{
    if ( argc < 3 || argc > 4)
        return 1;

    N = atoi(argv[1]);
    M = atoi(argv[2]);
    if(argc == 4)
        Z = atoi(argv[3]);
    cout << "N = " << N << "   M = " << M << "   Z = " << Z << endl;

    if (N < 3 || M < 3 || Z < 3) {
        cout << "N, M, Z cannot be less than 3" << endl;
        return 1;
    }
    if (M % 2 == 0 || Z % 2 == 0) {
        cout << "M and Z must be odd" << endl;
        return 1;
    }

    hx = Width / (N - 1);
    hy = Height / (N - 1);

    outputPoly();
    outputNode();
    outputElement();
    outputEdge();
    // genNode();
    // genEdge();
    outputRefineRound();
    outputRefineRound1();
}

void genRefineRoundNode()
{
    int nNode = N * N;

    // new node on row edge
    for (int i = 1; i <= N; i++) {
        if (i > 2 && i < N - 1) {
            for (int j = 2; j <= M - 1; j ++)
                vecNode.push_back(Node(++nNode, Left + (j - 1) * mx, Bottom + (i - 1) * hy, 0));
            for (int j = 2; j <= M - 1; j ++)
                vecNode.push_back(Node(++nNode, Left + Width - hx + (j - 1) * mx, Bottom + (i - 1) * hy, 0));
        } else {
            for (int j = 2; j <= (M - 1) * (N - 1); j++) {
                if ((j - 1) % (M - 1) == 0) continue;
                if (i == 1 || i == N)
                    vecNode.push_back(Node(++nNode, Left + (j - 1) * mx, Bottom + (i - 1) * hy, 1));
                else
                    vecNode.push_back(Node(++nNode, Left + (j - 1) * mx, Bottom + (i - 1) * hy, 0));
            }
        }
    }

    // new node on col edge
    for (int j = 1; j <= N; j++) {
        if (j > 2 && j < N - 1) {
            for (int i = 2; i <= M - 1; i++)
                vecNode.push_back(Node(++nNode, Left + (j - 1) * hx, Bottom + (i - 1) * my, 0));
            for (int i = 2; i <= M - 1; i++)
                vecNode.push_back(Node(++nNode, Left + (j - 1) * hx, Bottom + Height - hy + (i - 1) * my, 0));
        } else {
            for (int i = 2; i <= (M - 1) * (N - 1); i++) {
                if ((i - 1) % (M - 1) == 0) continue;
                if (j == 1 || j == N)
                    vecNode.push_back(Node(++nNode, Left + (j - 1) * hx, Bottom + (i - 1) * my, 1));
                else
                    vecNode.push_back(Node(++nNode, Left + (j - 1) * hx, Bottom + (i - 1) * my, 0));
            }
        }
    }

    // inside each square
    for (int i = 1; i <= N - 1; i++)
        for (int j = 1; j <= N - 1; j++) {
            if (i > 1 && i < N - 1 && j > 1 && j < N - 1 ) continue;
            double startx = Left + (j - 1) * hx;
            double starty = Bottom + (i - 1) * hy;
            for (int mi = 2; mi <= M - 1; mi++)
                for (int mj = 2; mj <= M - 1; mj++)
                    vecNode.push_back(Node(++nNode, startx + (mj - 1) * mx, starty + (mi - 1) * my, 0));
        }
}

int findNode(double x, double y)
{
    for (Node node : vecNode)
        if (fabs(node.x - x) < 0.0001 && fabs(node.y - y) < 0.0001)
            return node.index;
    return 0;
}

int findElement(int v1, int v2, int v3)
{
    for (Element ele : vecElement)
        if (ele.v1 == v1 && ele.v2 == v2 && ele.v3 == v3)
            return ele.index;
    return 0;
}

int findEdge(int v1, int v2)
{
    for (Edge edge : vecEdge)
        if (edge.v1 == v1 && edge.v2 == v2)
            return edge.index;
    return 0;
}

void genRefindEdgeOnSquareEdge(int indexEdge, int indexNode, int m)
{
    int nEdge = vecEdge.size();

    vecEdge.push_back(Edge(++nEdge, vecEdge[indexEdge - 1].v1, indexNode, vecEdge[indexEdge - 1].bctype, indexEdge));
    for (int i = 0; i < m - 3; i++)
        vecEdge.push_back(Edge(++nEdge, indexNode + i, indexNode + i + 1, vecEdge[indexEdge - 1].bctype, indexEdge));
    vecEdge.push_back(Edge(++nEdge, indexNode + m - 3, vecEdge[indexEdge - 1].v2, vecEdge[indexEdge - 1].bctype, indexEdge));
}

void genRefindEdgeOnSquare(int j, int i, int flagBottom, int flagRight, int flagTop, int flagLeft)
{
    int indexEdgeBottom = (i - 1) * (N - 1) + j;
    int indexEdgeTop = i * (N - 1) + j;
    int indexEdgeLeft = N * (N - 1) + i + (j - 1) * (N - 1);
    int indexEdgeRight = N * (N - 1) + i + j * (N - 1);
    int indexEdgeDiag = 2 * N * (N - 1) + (i - 1) * (N - 1) + j;

    int indexNodeBottom = findNode(Left + (j - 1) * hx + mx, Bottom + (i - 1) * hy);
    int indexNodeTop = findNode(Left + (j - 1) * hx + mx, Bottom + i * hy);
    int indexNodeLeft = findNode(Left + (j - 1) * hx, Bottom + (i - 1) * hy + my);
    int indexNodeRight = findNode(Left + j * hx, Bottom + (i - 1) * hy + my);
    // int indexNodeInterior = findNode(Left + (j - 1) * hx + mx, Bottom + (i - 1) * hy + my);

    if (flagBottom)
        genRefindEdgeOnSquareEdge(indexEdgeBottom, indexNodeBottom, M);

    if (flagTop)
        genRefindEdgeOnSquareEdge(indexEdgeTop, indexNodeTop, M);

    if (flagLeft)
        genRefindEdgeOnSquareEdge(indexEdgeLeft, indexNodeLeft, M);

    if (flagRight)
        genRefindEdgeOnSquareEdge(indexEdgeRight, indexNodeRight, M);

    int nEdge = vecEdge.size();
    double startx = Left + (j - 1) * hx;
    double starty = Bottom + (i - 1) * hy;
    // interior row
    for (int ii = 2; ii <= M - 1; ii++) {
        for (int jj = 1; jj <= M - 1; jj++) {
            int left = findNode(startx + (jj - 1) * mx, starty + (ii - 1) * my);
            int right = findNode(startx + jj * mx, starty + (ii - 1) * my);
            vecEdge.push_back(Edge(++nEdge, left, right, 0, 0));
        }
    }

    // interior col
    for (int jj = 2; jj <= M - 1; jj++) {
        for (int ii = 1; ii <= M - 1; ii++) {
            int bottom = findNode(startx + (jj - 1) * mx, starty + (ii - 1) * my);
            int top = findNode(startx + (jj - 1) * mx, starty + ii * my);
            vecEdge.push_back(Edge(++nEdge, bottom, top, 0, 0));
        }
    }

    // diag
    for (int ii = 1; ii <= M - 1; ii++) {
        for (int jj = 1; jj <= M - 1; jj++) {
            int v1 = findNode(startx + (jj - 1) * mx, starty + (ii - 1) * my);
            int v2 = findNode(startx + jj * mx, starty + ii * my);
            if (ii == jj)
                vecEdge.push_back(Edge(++nEdge, v1, v2, 0, indexEdgeDiag));
            else
                vecEdge.push_back(Edge(++nEdge, v1, v2, 0, 0));
        }
    }
}

void genRefineRoundElement()
{
    int nElement = vecElement.size();
    for (int i = 1; i <= N - 1; i++)
        for (int j = 1; j <= N - 1; j++) {
            if (i > 1 && i < N - 1 && j > 1 && j < N - 1 ) continue;
            double startx = Left + (j - 1) * hx;
            double starty = Bottom + (i - 1) * hy;
            int element1 = (2 * N - 2) * (i - 1) + 2 * j - 1;
            int element2 = element1 + 1;
            for (int ii = 1; ii <= M - 1; ii++)
                for (int jj = 1; jj <= M - 1; jj ++) {
                    int v1 = findNode(startx + (jj - 1) * mx, starty + (ii - 1) * my);
                    int v2 = findNode(startx + jj * mx, starty + (ii - 1) * my);
                    int v3 = findNode(startx + (jj - 1) * mx, starty + ii * my);
                    int v4 = findNode(startx + jj * mx, starty + ii * my);

                    if (ii >= jj)
                        vecElement.push_back(Element(++nElement, v1, v4, v3, element1));
                    else
                        vecElement.push_back(Element(++nElement, v1, v4, v3, element2));

                    if (ii <= jj)
                        vecElement.push_back(Element(++nElement, v1, v2, v4, element2));
                    else
                        vecElement.push_back(Element(++nElement, v1, v2, v4, element1));
                }
        }
}

void outputRefineRound()
{
    ofstream fout((filename + ".1.ref0").c_str());

    mx = hx / (M - 1);
    my = hy / (M - 1);

    genRefineRoundNode();

    fout << vecNode.size() - N *N << endl;
    for (int i = N * N; i < vecNode.size(); i++)
        fout << "   " << vecNode[i].index << "\t" << vecNode[i].x << "\t" << vecNode[i].y << "\t" << vecNode[i].bctype << endl;

    genRefindEdgeOnSquare(1, 1, 1, 1, 1, 1);
    for (int j = 2; j <= N - 1; j++)
        genRefindEdgeOnSquare(j, 1, 1, 1, 1, 0);
    for (int i = 2; i <= N - 2; i++) {
        genRefindEdgeOnSquare(1, i, 0, 1, 1, 1);
        genRefindEdgeOnSquare(N - 1, i, 0, 1, 1, 1);
    }
    genRefindEdgeOnSquare(1, N - 1, 0, 1, 1, 1);
    for (int j = 2; j <= N - 2; j++)
        genRefindEdgeOnSquare(j, N - 1, 1, 1, 1, 0);
    genRefindEdgeOnSquare(N - 1, N - 1, 0, 1, 1, 0);

    fout << vecEdge.size() - (N - 1) * (3 * N - 1) << endl;
    int currentParent = 0;
    for (int i = (N - 1) * (3 * N - 1); i < vecEdge.size(); i++) {
        if (vecEdge[i].parent != currentParent) {
            currentParent = vecEdge[i].parent;
            fout << "\t" << vecEdge[i].parent << endl;
        }
        fout << "\t\t" << vecEdge[i].index << "\t" << vecEdge[i].v1 << "\t" << vecEdge[i].v2 << endl;
        // fout << "   " << vecEdge[i].index << "\t" << vecEdge[i].v1 << "\t" << vecEdge[i].v2 << "\t" << vecEdge[i].bctype << endl;
    }

    genRefineRoundElement();
    fout << vecElement.size() << endl;
    currentParent = 0;
    for (int i = 0; i < vecElement.size(); i++) {
        if (vecElement[i].parent != currentParent) {
            currentParent = vecElement[i].parent;
            fout << "\t" << vecElement[i].parent << endl;
        }
        fout << "\t\t" << 2 * (N - 1) * (N - 1) + vecElement[i].index << "\t" << vecElement[i].v1 << "\t" << vecElement[i].v2 << "\t" << vecElement[i].v3 << endl;
    }

}

void genRefineRoundNode1()
{
    int nNode = vecNode.size();

    // new node on row edge
    for (int i = 1; i <= (N - 1) * (M - 1) + 1; i++) {
        if (i > 2 && i < (N - 1) * (M - 1) ) {
            for (int j = 2; j <= Z - 1; j ++)
                vecNode.push_back(Node(++nNode, Left + (j - 1) * zx, Bottom + (i - 1) * my, 0));
            for (int j = 2; j <= Z - 1; j ++)
                vecNode.push_back(Node(++nNode, Left + Width - mx + (j - 1) * zx, Bottom + (i - 1) * my, 0));
        } else {
            for (int j = 2; j <= (M - 1) * (N - 1) * (Z - 1); j++) {
                if ((j - 1) % (Z - 1) == 0) continue;
                if (i == 1 || i == (N - 1) * (M - 1) + 1)
                    vecNode.push_back(Node(++nNode, Left + (j - 1) * zx, Bottom + (i - 1) * my, 1));
                else
                    vecNode.push_back(Node(++nNode, Left + (j - 1) * zx, Bottom + (i - 1) * my, 0));
            }
        }
    }

    // new node on col edge
    for (int j = 1; j <= (N - 1) * (M - 1) + 1; j++) {
        if (j > 2 && j < (N - 1) * (M - 1) ) {
            for (int i = 2; i <= Z - 1; i++)
                vecNode.push_back(Node(++nNode, Left + (j - 1) * mx, Bottom + (i - 1) * zy, 0));
            for (int i = 2; i <= Z - 1; i++)
                vecNode.push_back(Node(++nNode, Left + (j - 1) * mx, Bottom + Height - my + (i - 1) * zy, 0));
        } else {
            for (int i = 2; i <= (M - 1) * (N - 1) * (Z - 1); i++) {
                if ((i - 1) % (Z - 1) == 0) continue;
                if (j == 1 || j == (N - 1) * (M - 1) + 1)
                    vecNode.push_back(Node(++nNode, Left + (j - 1) * mx, Bottom + (i - 1) * zy, 1));
                else
                    vecNode.push_back(Node(++nNode, Left + (j - 1) * mx, Bottom + (i - 1) * zy, 0));
            }
        }
    }

    // inside each square
    for (int i = 1; i <= (N - 1) * (M - 1); i++)
        for (int j = 1; j <= (N - 1) * (M - 1); j++) {
            if (i > 1 && i < (N - 1) * (M - 1) && j > 1 && j < (N - 1) * (M - 1) ) continue;
            double startx = Left + (j - 1) * mx;
            double starty = Bottom + (i - 1) * my;
            for (int zi = 2; zi <= Z - 1; zi++)
                for (int zj = 2; zj <= Z - 1; zj++)
                    vecNode.push_back(Node(++nNode, startx + (zj - 1) * zx, starty + (zi - 1) * zy, 0));
        }
}

void genRefineRoundElement1()
{
    int nElement = vecElement.size();
    for (int i = 1; i <= (N - 1) * (M - 1); i++)
        for (int j = 1; j <= (N - 1) * (M - 1); j++) {
            if (i > 1 && i < (N - 1) * (M - 1) && j > 1 && j < (N - 1) * (M - 1) ) continue;
            double startx = Left + (j - 1) * mx;
            double starty = Bottom + (i - 1) * my;
            int startV1 = findNode(startx, starty);
            int startV2 = findNode(startx + mx, starty);
            int startV3 = findNode(startx, starty + my);
            int startV4 = findNode(startx + mx, starty + my);
            int element1 = findElement(startV1, startV4, startV3) + 2 * (N - 1) * (N - 1);
            int element2 = findElement(startV1, startV2, startV4) + 2 * (N - 1) * (N - 1);

            for (int ii = 1; ii <= Z - 1; ii++)
                for (int jj = 1; jj <= Z - 1; jj ++) {
                    int v1 = findNode(startx + (jj - 1) * zx, starty + (ii - 1) * zy);
                    int v2 = findNode(startx + jj * zx, starty + (ii - 1) * zy);
                    int v3 = findNode(startx + (jj - 1) * zx, starty + ii * zy);
                    int v4 = findNode(startx + jj * zx, starty + ii * zy);

                    if (ii >= jj)
                        vecElement.push_back(Element(++nElement, v1, v4, v3, element1));
                    else
                        vecElement.push_back(Element(++nElement, v1, v4, v3, element2));

                    if (ii <= jj)
                        vecElement.push_back(Element(++nElement, v1, v2, v4, element2));
                    else
                        vecElement.push_back(Element(++nElement, v1, v2, v4, element1));
                }
        }
}

void genRefindEdgeOnSquare1(int j, int i, int flagBottom, int flagRight, int flagTop, int flagLeft)
{
    double startx = Left + (j - 1) * mx;
    double starty = Bottom + (i - 1) * my;
    int startV1 = findNode(startx, starty);
    int startV2 = findNode(startx + mx, starty);
    int startV3 = findNode(startx, starty + my);
    int startV4 = findNode(startx + mx, starty + my);

    int indexEdgeBottom = findEdge(startV1, startV2);
    int indexEdgeTop = findEdge(startV3, startV4);
    int indexEdgeLeft = findEdge(startV1, startV3);
    int indexEdgeRight = findEdge(startV2, startV4);
    int indexEdgeDiag = findEdge(startV1, startV4);

    int indexNodeBottom = findNode(Left + (j - 1) * mx + zx, Bottom + (i - 1) * my);
    int indexNodeTop = findNode(Left + (j - 1) * mx + zx, Bottom + i * my);
    int indexNodeLeft = findNode(Left + (j - 1) * mx, Bottom + (i - 1) * my + zy);
    int indexNodeRight = findNode(Left + j * mx, Bottom + (i - 1) * my + zy);

    if (flagBottom)
        genRefindEdgeOnSquareEdge(indexEdgeBottom, indexNodeBottom, Z);

    if (flagTop)
        genRefindEdgeOnSquareEdge(indexEdgeTop, indexNodeTop, Z);

    if (flagLeft)
        genRefindEdgeOnSquareEdge(indexEdgeLeft, indexNodeLeft, Z);

    if (flagRight)
        genRefindEdgeOnSquareEdge(indexEdgeRight, indexNodeRight, Z);

    int nEdge = vecEdge.size();
    // interior row
    for (int ii = 2; ii <= Z - 1; ii++) {
        for (int jj = 1; jj <= Z - 1; jj++) {
            int left = findNode(startx + (jj - 1) * zx, starty + (ii - 1) * zy);
            int right = findNode(startx + jj * zx, starty + (ii - 1) * zy);
            if(left == 0 || right == 0)
                cout << i << " " << j << " " << ii << " " << jj << " " << left << " " << right << endl;
            vecEdge.push_back(Edge(++nEdge, left, right, 0, 0));
        }
    }

    // interior col
    for (int jj = 2; jj <= Z - 1; jj++) {
        for (int ii = 1; ii <= Z - 1; ii++) {
            int bottom = findNode(startx + (jj - 1) * zx, starty + (ii - 1) * zy);
            int top = findNode(startx + (jj - 1) * zx, starty + ii * zy);
            vecEdge.push_back(Edge(++nEdge, bottom, top, 0, 0));
        }
    }

    // diag
    for (int ii = 1; ii <= Z - 1; ii++) {
        for (int jj = 1; jj <= Z - 1; jj++) {
            int v1 = findNode(startx + (jj - 1) * zx, starty + (ii - 1) * zy);
            int v2 = findNode(startx + jj * zx, starty + ii * zy);
            if (ii == jj)
                vecEdge.push_back(Edge(++nEdge, v1, v2, 0, indexEdgeDiag));
            else
                vecEdge.push_back(Edge(++nEdge, v1, v2, 0, 0));
        }
    }
}

void outputRefineRound1()
{
    ofstream fout((filename + ".1.ref1").c_str());

    zx = mx / (Z - 1);
    zy = my / (Z - 1);

    int nNodeRef0 = vecNode.size();
    genRefineRoundNode1();

    fout << vecNode.size() - nNodeRef0 << endl;
    for (int i = nNodeRef0; i < vecNode.size(); i++)
        fout << "   " << vecNode[i].index << "\t" << vecNode[i].x << "\t" << vecNode[i].y << "\t" << vecNode[i].bctype << endl;

    int nEdgeRef0 = vecEdge.size();
    genRefindEdgeOnSquare1(1, 1, 1, 1, 1, 1);
    for (int j = 2; j <= (N - 1) * (M - 1); j++)
        genRefindEdgeOnSquare1(j, 1, 1, 1, 1, 0);
    for (int i = 2; i <= (N - 1) * (M - 1) - 1; i++) {
        genRefindEdgeOnSquare1(1, i, 0, 1, 1, 1);
        genRefindEdgeOnSquare1((N - 1) * (M - 1), i, 0, 1, 1, 1);
    }
    genRefindEdgeOnSquare1(1, (N - 1) * (M - 1), 0, 1, 1, 1);
    for (int j = 2; j <= (N - 1) * (M - 1) - 1; j++)
        genRefindEdgeOnSquare1(j, (N - 1) * (M - 1), 1, 1, 1, 0);
    genRefindEdgeOnSquare1((N - 1) * (M - 1), (N - 1) * (M - 1), 0, 1, 1, 0);

    fout << vecEdge.size() - nEdgeRef0 << endl;
    int currentParent = -1;
    for (int i = nEdgeRef0; i < vecEdge.size(); i++) {
        if (vecEdge[i].parent != currentParent) {
            currentParent = vecEdge[i].parent;
            fout << "\t" << vecEdge[i].parent << endl;
        }
        fout << "\t\t" << vecEdge[i].index << "\t" << vecEdge[i].v1 << "\t" << vecEdge[i].v2 << endl;
        // fout << "   " << vecEdge[i].index << "\t" << vecEdge[i].v1 << "\t" << vecEdge[i].v2 << "\t" << vecEdge[i].bctype << endl;
    }

    int nElementRef0 = vecElement.size();
    genRefineRoundElement1();
    fout << vecElement.size() - nElementRef0 << endl;
    currentParent = -1;
    for (int i = nElementRef0; i < vecElement.size(); i++) {
        if (vecElement[i].parent != currentParent) {
            currentParent = vecElement[i].parent;
            fout << "\t" << vecElement[i].parent << endl;
        }
        fout << "\t\t" << 2 * (N - 1) * (N - 1) + vecElement[i].index << "\t" << vecElement[i].v1 << "\t" << vecElement[i].v2 << "\t" << vecElement[i].v3 << endl;
    }

}
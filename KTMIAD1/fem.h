#pragma once
#ifndef FEM_H
#define FEM_H

#include "grid.h"

typedef std::function<double(double, double)> ScalFunc2D;
typedef std::function<double(double)> ScalFunc1D;

struct FE
{
   std::vector<int> nodes;  // � ������� �� 8 �����
   std::vector<int> edges;  // � ������� �� 12 ����
   std::vector<std::pair<int, int>> nodes_edges; //� ������� ����� 2 ����
};

struct edge
{
   std::vector<int> nodes;  //� ������ ����� 2 ����
};

struct bc
{
   int node;
   double value;
};

struct bc_2
{
   int ielem;
   int node;
   double value;
};

class FEM
{
public:
   void init();
   void Compute();
   void generate_FE();
   void edge_from_nodes(int node1, int node2);
   int FE_from_edge(int edge);
   void buildPortraitOfMatrix();
   bool isOrdered(const std::vector<int>& v);
   void assemblyGlobalMatrix();
   void addElementToGlobal(int i, int j, double elem);
   void getLocalG(double hx, double hy, double hz);
   void getLocalM(double hx, double hy, double hz);
   void getLocalRightPart(point p, double hx, double hy, double hz, int ielem, int cur_tLayer);
   void bc1(int cur_tLayer);
   double bc1_value(int iedge, int cur_tLayer);
   double RP_value(int iedge, int cur_tLayer);
   int globalToLocalEdge(int ielem, int gl_num);
   int globalToLocalNode(int ielem, int gl_num);
   bool isInBoundary(int num);
   void CheckSol(int cur_tLayer);


   std::vector<double> q;

private:
   grid g;
   std::vector<FE> _elems;                // ������ � ��������� ����������
   std::vector<edge> _edges;              // ������ � �������
   int nFE, nUz;                          // ����� �������� ��������� � �����
   std::vector<int> ig, jg;               // ������� �������� ��� ��������
   std::vector<double> points, weights;   // ������� ��� ����� � ���� ������-3
   std::vector<double> alpha;             // ������� ������������� ��� �������� �� � ���������� ��������
   std::vector<double> beta;
   std::vector<double> x_l, y_l;          // ������� � ������������ ��������� ����� ��

   std::vector<double> di, gg, F;      // ������� ��� ����
   std::vector<bc> _bc1;
   std::vector<bc_2> _bc2;
   double Pi = 3.1415926535;
   std::vector<int> result;
   int n_edges;
   std::vector<std::vector<double>> M_loc, G_loc;
   std::vector<double> F_loc;
   double mu = 1.0;
   double gamma = 1.0;
   int num_test;
   int nx, ny, nz;

   //��������� �����
   std::vector<double> time_layer;
   int time_size;
   std::vector<double> q0; //t-1
   std::vector<double> q1; //t-2
   double dt, dt0, dt1;
   double d1, d2, d3;
   double t, t_1, t_2;
   int cur_tLayer;

   std::vector<std::vector<double>> G1 =
   {
	   {7,-8,1},
	   {-8,16,-8},
	   {1,-8,7}
   };

   std::vector<std::vector<double>> M1 =
   {
	   {4,2,-1},
	   {2,16,2},
	   {-1,2,4}
   };

   std::vector<std::vector<double>> M =
   {
	   {4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0},
	   {2, 4, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0},
	   {2, 1, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0},
	   {1, 2, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0},
	   {0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0},
	   {0, 0, 0, 0, 2, 4, 1, 2, 0, 0, 0, 0},
	   {0, 0, 0, 0, 2, 1, 4, 2, 0, 0, 0, 0},
	   {0, 0, 0, 0, 1, 2, 2, 4, 0, 0, 0, 0},
	   {0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 2, 1},
	   {0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 1, 2},
	   {0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 4, 2},
	   {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 4}
   };

  /* std::vector<std::vector<double>> G1 =
   {
	   {2, 1, -2, -1},
	   {1, 2, -1, -2},
	   {-2, -1, 2, 1},
	   {-1, -2, 1, 2}
   };*/

   std::vector<std::vector<double>> G2 =
   {
	   {2, -2, 1, -1},
	   {-2, 2, -1, 1},
	   {1, -1, 2, -2},
	   {-1, 1, -2, 2}
   };

   std::vector<std::vector<double>> G3 =
   {
	   {-2, 2, -1, 1},
	   {-1, 1, -2, 2},
	   {2, -2, 1, -1},
	   {1, -1, 2, -2}
   };

   std::vector<std::vector<double>> G3_T =
   {
	   {-2, -1, 2, 1},
	   {2, 1, -2, -1},
	   {-1, -2, 1, 2},
	   {1, 2, -1, -2}
   };
};
#endif
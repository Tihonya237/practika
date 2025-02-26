#pragma once
#ifndef GRID_H
#define GRID_H

#include "resource.h"

class grid
{
public:
   double angle_step(int n, double k);
   double coord_step(double x0, double x1 ,int n, double k);
   void get_grid();

   std::vector<point> MeshXYZ;
   int n_x, n_y, n_z;

private:
   std::vector<point> st;                          // Ключевые точки
   std::vector<point> points;                      // Набор точек сетки
   std::vector<double> _x, _y, _z;
   double Pi = 3.1415926535;
   double k_x, k_y, k_z;
};

#endif
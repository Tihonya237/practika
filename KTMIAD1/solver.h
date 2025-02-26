#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include "resource.h"

class solver
{
public:
   std::vector<double> BCGSTAB(std::vector<double> &gg, std::vector<double>& di, std::vector<double>& F, std::vector<int> &_ig, std::vector<int>& _jg);
   void mult_reverse(std::vector<double>& ggu, std::vector<double>& di, std::vector<double>& res, std::vector<double>& b);
   void mult_direct(std::vector<double>& ggl, std::vector<double>& di, std::vector<double>& res, std::vector<double>& b);
   void LU();
   std::vector<double> mult(const std::vector<double>& v);
   std::vector<double> LOS_LU(std::vector<double>& _gg, std::vector<double>& _di, std::vector<double>& F, std::vector<int>& _ig, std::vector<int>& _jg);

   std::vector<double> q;

private:
   std::vector<double> di_n, ggu_n, ggl_n, r, p, di, gg;
   std::vector<int> ig, jg;
   int maxIter = 1000;
   double eps = 1e-15;
};

inline std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b)
{
   std::vector<double> res = a;
   for (int i = 0; i < res.size(); i++)
      res[i] += b[i];
   return res;
}
inline std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b)
{
   std::vector<double> res = a;
   for (int i = 0; i < res.size(); i++)
      res[i] -= b[i];
   return res;
}
inline double operator*(const std::vector<double>& a, const std::vector<double>& b)
{
   double scalar = 0.0;
   for (int i = 0; i < a.size(); i++)
      scalar += a[i] * b[i];
   return scalar;
}
inline std::vector<double> operator*(double c, const std::vector<double>& a)
{
   std::vector<double> res = a;
   for (int i = 0; i < res.size(); i++)
      res[i] *= c;
   return res;
}
#endif
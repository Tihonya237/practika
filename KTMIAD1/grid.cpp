#include "grid.h"

double grid::angle_step(int n, double k)
{
   if (k != 1)
      return 90.0 * (1 - k) / (1 - pow(k, n));
   else
      return 90.0 / n;
}

double grid::coord_step(double x0, double x1, int n, double k)
{
   if (k != 1)
      return (x1 - x0) * (1 - k) / (1 - pow(k, n));
   else
      return (x1 - x0) / n;
}

void grid::get_grid()
{
   st.resize(2);
   std::ifstream input("files/input.json ");
   nlohmann::json inGrid{};
   if (input.is_open())
   {
      input >> inGrid;

      input.close();

      st[0].x = inGrid["area"]["start"][0];
      st[0].y = inGrid["area"]["start"][1];
      st[0].z = inGrid["area"]["start"][2];

      st[1].x = inGrid["area"]["end"][0];
      st[1].y = inGrid["area"]["end"][1];
      st[1].z = inGrid["area"]["end"][2];

      n_x = inGrid["area"]["n"][0];
      n_y = inGrid["area"]["n"][1];
      n_z = inGrid["area"]["n"][2];

      k_x = inGrid["area"]["k"][0];
      k_y = inGrid["area"]["k"][1];
      k_z = inGrid["area"]["k"][2];
   }
   else throw "Can't open file input.json\n";

   double tmp, h;

   //Обработка по X
   if (k_x != 1)
      h = (st[1].x - st[0].x) * (1. - k_x) / (1. - pow(k_x, n_x));
   else
      h = (st[1].x - st[0].x) / n_x;

   _x.push_back(st[0].x);
   for (int i = 0; i < n_x; i++)
   {
      _x.push_back(_x.back() + h);
      h *= k_x;
   }

   //Обработка по y
   if (k_y != 1)
      h = (st[1].y - st[0].y) * (1. - k_y) / (1. - pow(k_y, n_y));
   else
      h = (st[1].y - st[0].y) / n_y;

   _y.push_back(st[0].y);
   for (int i = 0; i < n_y; i++)
   {
      _y.push_back(_y.back() + h);
      h *= k_y;
   }

   //Обработка по z
   if (k_z != 1)
      h = (st[1].z - st[0].z) * (1. - k_z) / (1. - pow(k_z, n_z));
   else
      h = (st[1].z - st[0].z) / n_z;

   _z.push_back(st[0].z);
   for (int i = 0; i < n_z; i++)
   {
      _z.push_back(_z.back() + h);
      h *= k_z;
   }

   //Заполнение структуры сетки
   for (int i = 0; i < _z.size(); i++)
      for (int j = 0; j < _y.size(); j++)
         for (int k = 0; k < _x.size(); k++)
            MeshXYZ.push_back({_x[k], _y[j], _z[i]});
}
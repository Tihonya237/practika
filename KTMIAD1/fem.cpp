#include "fem.h"
#include "solver.h"

void FEM::Compute()
{
   double start_time, end_time, prog_time;

   g.get_grid();
   init();
   generate_FE();
   std::cout << "Сборка СЛАУ началась." << std::endl;
   buildPortraitOfMatrix();
   assemblyGlobalMatrix(0);
   bc1();
   std::cout << "Сборка СЛАУ завершена." << std::endl;
   std::cout << "Решение СЛАУ началось." << std::endl;
   solver s;
   start_time = clock();
   q = s.BCGSTAB(gg, di, F, ig, jg);
   end_time = clock();
   std::cout << "Решение СЛАУ завершено." << std::endl;
   CheckSol(0.0);
   prog_time = end_time - start_time;

   std::cout << prog_time / 1000 << std::endl;

   for (int i = 1; i < time_size; i++)
   {
       start_time = clock();
       std::cout << i;
       assemblyGlobalMatrix(i);
       std::cout << "Решение СЛАУ началось." << std::endl;
       solver s;
       q = s.BCGSTAB(gg, di, F, ig, jg);
       std::cout << "Решение СЛАУ завершено." << std::endl;
       end_time = clock();
       CheckSol(i);

       prog_time = end_time - start_time;

       std::cout << prog_time / 1000 << std::endl;
   }
   system("pause");
}

void FEM::init()
{
   std::ifstream input("files/input.json");
   nlohmann::json inParam{};

   if (input.is_open())
   {
      input >> inParam;

      input.close();

      mu = inParam["parameters"]["mu"];
      gamma = inParam["parameters"]["gamma"];
      num_test = inParam["parameters"]["num_test"];

      time_size = inParam["time_grid"]["time_size"];

      time_layer.resize(time_size);

      for (int i = 0; i < time_size; i++)
          time_layer[i] = inParam["time_grid"]["time_st_en"][i];
   }
   else throw "Can't open file input.json\n";
   nUz = g.MeshXYZ.size();
}

// Получение значения первого краевого условия
double FEM::bc1_value(int iedge, int cur_tLayer)
{
    double t = time_layer[cur_tLayer];
   int ielem = FE_from_edge(iedge);
   int loc_edge = globalToLocalEdge(ielem, iedge);
   int loc_node1 = globalToLocalNode(ielem, _elems[ielem].nodes_edges[loc_edge].first);
   int loc_node2 = globalToLocalNode(ielem, _elems[ielem].nodes_edges[loc_edge].second);

   point p1, p2;
   p1.x = g.MeshXYZ[_elems[ielem].nodes[loc_node1]].x;
   p1.y = g.MeshXYZ[_elems[ielem].nodes[loc_node1]].y;
   p1.z = g.MeshXYZ[_elems[ielem].nodes[loc_node1]].z;

   p2.x = g.MeshXYZ[_elems[ielem].nodes[loc_node2]].x;
   p2.y = g.MeshXYZ[_elems[ielem].nodes[loc_node2]].y;
   p2.z = g.MeshXYZ[_elems[ielem].nodes[loc_node2]].z;

   double hx = abs(p1.x - p2.x);
   double hy = abs(p1.y - p2.y);
   double hz = abs(p1.z - p2.z);

   int var;

   if (hx > 0) var = 0;
   else 
      if (hy > 0) var = 1;
      else
         if (hz > 0) var = 2;

   switch(var)
   {
   case 0:
      p1.x = (g.MeshXYZ[_elems[ielem].nodes[loc_node1]].x + g.MeshXYZ[_elems[ielem].nodes[loc_node2]].x) / 2.0;
   case 1:
      p1.y = (g.MeshXYZ[_elems[ielem].nodes[loc_node1]].y + g.MeshXYZ[_elems[ielem].nodes[loc_node2]].y) / 2.0;
   case 2:
      p1.z = (g.MeshXYZ[_elems[ielem].nodes[loc_node1]].z + g.MeshXYZ[_elems[ielem].nodes[loc_node2]].z) / 2.0;
   }

   switch (num_test)
   {
   case 1:
      switch (var)
      {
      case 0:
         return p1.x + p1.y + p1.z + 2 *  t;
      case 1:
         return p1.x + 2 * p1.y + p1.z + t;
      case 2:
         return p1.x + p1.y + 3 * p1.z + 3 * t;
      }
   case 2:
      switch (var)
      {
      case 0:
         return pow(p1.x, 2) + pow(p1.y, 2) + pow(p1.z, 2) + pow(t, 2);
      case 1:
         return pow(p1.x, 2) + 2 * pow(p1.y, 2) + pow(p1.z, 2) + pow(t, 2);
      case 2:
         return pow(p1.x, 2) + pow(p1.y, 2) + 3 * pow(p1.z, 2) + pow(t, 2);
      }
   case 3:
      switch (var)
      {
      case 0:
         return pow(p1.x, 3) + pow(p1.y, 3) + pow(p1.z, 3) + pow(t, 3);
      case 1:
         return pow(p1.x, 3) + 2 * pow(p1.y, 3) + pow(p1.z, 3) + pow(t, 3);
      case 2:
         return pow(p1.x, 3) + pow(p1.y, 3) + 3 * pow(p1.z, 3) + pow(t, 3);
      }
   case 4:
      switch (var)
      {
      case 0:
         return pow(p1.x, 4) + pow(p1.y, 4) + pow(p1.z, 4) + pow(t, 4);
      case 1:
         return pow(p1.x, 4) + 2 * pow(p1.y, 4) + pow(p1.z, 4) + pow(t, 4);
      case 2:
         return pow(p1.x, 4) + pow(p1.y, 4) + 3 * pow(p1.z, 4) + pow(t, 4);
      }
   case 5:
      switch (var)
      {
      case 0:
         return sin(p1.y + p1.z);
      case 1:
         return sin(p1.x + p1.z);
      case 2:
         return sin(p1.x + p1.y);
      }
   }
   throw "Wrong number of test!";
}

double FEM::RP_value(int iedge)
{
   int ielem = FE_from_edge(iedge);
   int loc_edge = globalToLocalEdge(ielem, iedge);
   int loc_node1 = globalToLocalNode(ielem, _elems[ielem].nodes_edges[loc_edge].first);
   int loc_node2 = globalToLocalNode(ielem, _elems[ielem].nodes_edges[loc_edge].second);

   point p1, p2;
   p1.x = g.MeshXYZ[_elems[ielem].nodes[loc_node1]].x;
   p1.y = g.MeshXYZ[_elems[ielem].nodes[loc_node1]].y;
   p1.z = g.MeshXYZ[_elems[ielem].nodes[loc_node1]].z;

   p2.x = g.MeshXYZ[_elems[ielem].nodes[loc_node2]].x;
   p2.y = g.MeshXYZ[_elems[ielem].nodes[loc_node2]].y;
   p2.z = g.MeshXYZ[_elems[ielem].nodes[loc_node2]].z;

   double hx = abs(p1.x - p2.x);
   double hy = abs(p1.y - p2.y);
   double hz = abs(p1.z - p2.z);

   int var;

   if (hx > 0) var = 0;
   else
      if (hy > 0) var = 1;
      else
         if (hz > 0) var = 2;

   switch (var)
   {
   case 0:
      p1.x = (g.MeshXYZ[_elems[ielem].nodes[loc_node1]].x + g.MeshXYZ[_elems[ielem].nodes[loc_node2]].x) / 2.0;
   case 1:
      p1.y = (g.MeshXYZ[_elems[ielem].nodes[loc_node1]].y + g.MeshXYZ[_elems[ielem].nodes[loc_node2]].y) / 2.0;
   case 2:
      p1.z = (g.MeshXYZ[_elems[ielem].nodes[loc_node1]].z + g.MeshXYZ[_elems[ielem].nodes[loc_node2]].z) / 2.0;
   }

   switch (num_test)
   {
   case 1:
      switch (var)
      {
      case 0:
          return 1.0 / mu * (0);
      case 1:
          return 1.0 / mu * (0);
      case 2:
          return 1.0 / mu * (0);
      }
   case 2:
      switch (var)
      {
      case 0:
         return 1.0 / mu * (-4);
      case 1:
         return 1.0 / mu * (-4);
      case 2:
         return 1.0 / mu * (-4);
      }
   case 3:
      switch (var)
      {
      case 0:
         return 1.0 / mu * (-6 * (p1.y + p1.z));
      case 1:
         return 1.0 / mu * (-6 * (p1.x + p1.z));
      case 2:
         return 1.0 / mu * (-6 * (p1.x + p1.y));
      }
   case 4:
      switch (var)
      {
      case 0:
         return 1.0 / mu * (-12 * (pow(p1.y, 2) + pow(p1.z, 2)));
      case 1:
         return 1.0 / mu * (-12 * (pow(p1.x, 2) + pow(p1.z, 2)));
      case 2:
         return 1.0 / mu * (-12 * (pow(p1.x, 2) + pow(p1.y, 2)));
      }
   case 5:
      switch (var)
      {
      case 0:
         return 1.0 / mu * (2 * sin(p1.y + p1.z));
      case 1:
         return 1.0 / mu * (2 * sin(p1.x + p1.z));
      case 2:
         return 1.0 / mu * (2 * sin(p1.x + p1.y));
      }
   }
}
double FEM::dAdt(int iedge, double t)
{
    int ielem = FE_from_edge(iedge);
    int loc_edge = globalToLocalEdge(ielem, iedge);
    int loc_node1 = globalToLocalNode(ielem, _elems[ielem].nodes_edges[loc_edge].first);
    int loc_node2 = globalToLocalNode(ielem, _elems[ielem].nodes_edges[loc_edge].second);

    point p1, p2;
    p1.x = g.MeshXYZ[_elems[ielem].nodes[loc_node1]].x;
    p1.y = g.MeshXYZ[_elems[ielem].nodes[loc_node1]].y;
    p1.z = g.MeshXYZ[_elems[ielem].nodes[loc_node1]].z;

    p2.x = g.MeshXYZ[_elems[ielem].nodes[loc_node2]].x;
    p2.y = g.MeshXYZ[_elems[ielem].nodes[loc_node2]].y;
    p2.z = g.MeshXYZ[_elems[ielem].nodes[loc_node2]].z;

    double hx = abs(p1.x - p2.x);
    double hy = abs(p1.y - p2.y);
    double hz = abs(p1.z - p2.z);

    int var;

    if (hx > 0) var = 0;
    else
        if (hy > 0) var = 1;
        else
            if (hz > 0) var = 2;

    switch (var)
    {
    case 0:
        p1.x = (g.MeshXYZ[_elems[ielem].nodes[loc_node1]].x + g.MeshXYZ[_elems[ielem].nodes[loc_node2]].x) / 2.0;
    case 1:
        p1.y = (g.MeshXYZ[_elems[ielem].nodes[loc_node1]].y + g.MeshXYZ[_elems[ielem].nodes[loc_node2]].y) / 2.0;
    case 2:
        p1.z = (g.MeshXYZ[_elems[ielem].nodes[loc_node1]].z + g.MeshXYZ[_elems[ielem].nodes[loc_node2]].z) / 2.0;
    }

    switch (num_test)
    {
    case 1:
        switch (var)
        {
        case 0:
            return 2;
        case 1:
            return 1;
        case 2:
            return 3;
        }
    case 2:
        switch (var)
        {
        case 0:
            return 2 * t;
        case 1:
            return 2 * t;
        case 2:
            return 2 * t;
        }
    case 3:
        switch (var)
        {
        case 0:
            return 3 * pow(t, 2);
        case 1:
            return 3 * pow(t, 2);
        case 2:
            return 3 * pow(t, 2);
        }
    case 4:
        switch (var)
        {
        case 0:
            return 4 * pow(t, 3);
        case 1:
            return 4 * pow(t, 3);
        case 2:
            return 4 * pow(t, 3);
        }
    case 5:
        switch (var)
        {
        case 0:
            return 1.0 / mu * (2 * sin(p1.y + p1.z));
        case 1:
            return 1.0 / mu * (2 * sin(p1.x + p1.z));
        case 2:
            return 1.0 / mu * (2 * sin(p1.x + p1.y));
        }
    }
}

void FEM::generate_FE()
{
   nFE = g.n_x * g.n_y * g.n_z;
   std::vector<int> tmp_nodes(4);
   _elems.resize(nFE);
   for (int i = 0; i < nFE; i++)
   {
      _elems[i].nodes.resize(8);
      _elems[i].edges.resize(12);
      _elems[i].nodes_edges.resize(12);
   }

   nx = g.n_x + 1;
   ny = g.n_y + 1;
   nz = g.n_z + 1;
   // Глобальная нумерация узлов для каждого конечного элемента
   for (int i = 0, k = 0, m = 1; i < nz - 1; i++, m += nx)
   {
      for (int j = 0; j < ny - 1; j++, m++)
      {
         for (int p = 0; p < nx - 1; p++, m++, k++)
         {
            _elems[k].nodes[0] = m - 1;
            _elems[k].nodes[1] = m;
            _elems[k].nodes[2] = m + nx - 1;
            _elems[k].nodes[3] = m + nx;
            _elems[k].nodes[4] = m + ny * nx - 1;
            _elems[k].nodes[5] = m + ny * nx;
            _elems[k].nodes[6] = m + nx + ny * nx - 1;
            _elems[k].nodes[7] = m + nx + ny * nx;
         }
      }
   }

   // Глобальная нумерация ребер для каждого конечного элемента
   int m2 = (nx - 1) * ny + (ny - 1) * nx;
   n_edges = m2 * nz + nx * ny * (nz - 1);
   int r = m2 - nx * ny - 1;
   int m3 = m2 - r;
   int m4 = m2 + nx * ny;
   //for (int i = 1, k = 0, m = 0; i < nz; i++)
   //{
   //   m = (i - 1) * m4;
   //   for (int j = 0; j < ny - 1; j++, m += nx)
   //   {
   //      for (int p = 0; p < nx - 1; p++, m++, k++)
   //      {
   //         _elems[k].edges[0] = m;
   //         _elems[k].edges[4] = m + nx - 1;
   //         _elems[k].edges[5] = m + nx;
   //         _elems[k].edges[1] = m + 2 * nx - 1;
   //         _elems[k].edges[8] = (i - 1) * m4 + m2 + j * nx + p;
   //         _elems[k].edges[9] = (i - 1) * m4 + m2 + j * nx + p + 1;
   //         _elems[k].edges[10] = (i - 1) * m4 + m2 + j * nx + p + nx;
   //         _elems[k].edges[11] = (i - 1) * m4 + m2 + j * nx + p + nx + 1;
   //         _elems[k].edges[2] = m + m4;
   //         _elems[k].edges[6] = m + m4 + nx - 1;
   //         _elems[k].edges[7] = m + m4 + nx;
   //         _elems[k].edges[3] = m + m4 + 2 * nx - 1;
   //      }
   //   }
   //}

   for(int k = 0, tmp = 0; k < nz - 1; k++)
       for(int i = 0; i < ny - 1; i++)
          for (int j = 0; j < nx - 1; j++, tmp++)
          {
             int gr = nx * (ny - 1) + ny * (nx - 1);
             int pop = nx * ny;
             _elems[tmp].edges[0] = k * (gr + pop) + i * (2 * nx - 1) + j;
             _elems[tmp].nodes_edges[0].first = _elems[tmp].nodes[0];
             _elems[tmp].nodes_edges[0].second = _elems[tmp].nodes[1];

             _elems[tmp].edges[1] = k * (gr + pop) + (i + 1) * (2 * nx - 1) + j;
             _elems[tmp].nodes_edges[1].first = _elems[tmp].nodes[2];
             _elems[tmp].nodes_edges[1].second = _elems[tmp].nodes[3];

             _elems[tmp].edges[2] = (k + 1) * (gr + pop) + i * (2 * nx - 1) + j;
             _elems[tmp].nodes_edges[2].first = _elems[tmp].nodes[4];
             _elems[tmp].nodes_edges[2].second = _elems[tmp].nodes[5];

             _elems[tmp].edges[3] = (k + 1) * (gr + pop) + (i + 1) * (2 * nx - 1) + j;
             _elems[tmp].nodes_edges[3].first = _elems[tmp].nodes[6];
             _elems[tmp].nodes_edges[3].second = _elems[tmp].nodes[7];

             _elems[tmp].edges[4] = k * (gr + pop) + i * (2 * nx - 1) + j + (nx - 1);
             _elems[tmp].nodes_edges[4].first = _elems[tmp].nodes[0];
             _elems[tmp].nodes_edges[4].second = _elems[tmp].nodes[2];

             _elems[tmp].edges[5] = k * (gr + pop) + i * (2 * nx - 1) + j + 1 + (nx - 1);
             _elems[tmp].nodes_edges[5].first = _elems[tmp].nodes[1];
             _elems[tmp].nodes_edges[5].second = _elems[tmp].nodes[3];

             _elems[tmp].edges[6] = (k + 1) * (gr + pop) + i * (2 * nx - 1) + (nx - 1) + j;
             _elems[tmp].nodes_edges[6].first = _elems[tmp].nodes[4];
             _elems[tmp].nodes_edges[6].second = _elems[tmp].nodes[6];

             _elems[tmp].edges[7] = (k + 1) * (gr + pop) + i * (2 * nx - 1) + (nx - 1) + j + 1;
             _elems[tmp].nodes_edges[7].first = _elems[tmp].nodes[5];
             _elems[tmp].nodes_edges[7].second = _elems[tmp].nodes[7];

             _elems[tmp].edges[8] = k * (gr + pop) + gr + j + i * nx;
             _elems[tmp].nodes_edges[8].first = _elems[tmp].nodes[0];
             _elems[tmp].nodes_edges[8].second = _elems[tmp].nodes[4];

             _elems[tmp].edges[9] = k * (gr + pop) + gr + j + i * nx + 1;
             _elems[tmp].nodes_edges[9].first = _elems[tmp].nodes[1];
             _elems[tmp].nodes_edges[9].second = _elems[tmp].nodes[5];

             _elems[tmp].edges[10] = k * (gr + pop) + gr + j + i * nx + nx;
             _elems[tmp].nodes_edges[10].first = _elems[tmp].nodes[2];
             _elems[tmp].nodes_edges[10].second = _elems[tmp].nodes[6];

             _elems[tmp].edges[11] = k * (gr + pop) + gr + j + i * nx + nx + 1;
             _elems[tmp].nodes_edges[11].first = _elems[tmp].nodes[3];
             _elems[tmp].nodes_edges[11].second = _elems[tmp].nodes[7];
          }

   di.resize(n_edges);
   F.resize(n_edges);
   q.resize(n_edges);
}

void FEM::edge_from_nodes(int node1, int node2)
{
   int result = -1;
   for (int i = 0; i < _edges.size(); i++)
   {
      if (node1 == _edges[i].nodes[0] && node2 == _edges[i].nodes[1])
      {
         result = i;
         break;
      }
   }

   if (result != -1)
   {
      std::cout << "Номер ребра для данных узлов: " << result << std::endl << std::endl;
   }
   else
   {
      std::cout << "Введены некорректные номера узлов" << std::endl << std::endl;
   }
}

int FEM::FE_from_edge(int edge)
{
   result.clear();
   for (int i = 0; i < _elems.size(); i++)
   {
      for (int j = 0; j < 12; j++)
      {
         if (_elems[i].edges[j] == edge)
         {
            return i;
         }
      }
   }
}

void FEM::buildPortraitOfMatrix()
{
   std::vector<std::vector<int>> list;
   list.resize(n_edges);
   std::vector<int> tmp_edge;
   // Идем по всем КЭ
   for (int ielem = 0; ielem < nFE; ielem++)
   {
      tmp_edge = _elems[ielem].edges;
      sort(tmp_edge.begin(), tmp_edge.end());
       for (int i = 0; i < _elems[ielem].edges.size() - 1; i++)
       {
           for (int j = i + 1; j < _elems[ielem].edges.size(); j++)
           {
               int insertPos = tmp_edge[j];
               int element = tmp_edge[i];
               bool isIn = false;
               for (int k = 0; k < list[insertPos].size(); k++)
                  if (element == list[insertPos][k])
                  {
                     isIn = true;
                     break;
                  }
               if (!isIn)
                   list[insertPos].push_back(element);
           }
       }
       tmp_edge.clear();
   }
   // Сортируем все получившиеся списки (по возрастанию номеров)
   for (int i = 0; i < list.size(); i++)
      if (!isOrdered(list[i]))
         sort(list[i].begin(), list[i].end());

   // Формируем массив ig
   ig.resize(n_edges + 1);
   // 1-ый и 2-ой элементы всегда равны 1, но мы будем нумеровать с 0
   ig[0] = 0;
   ig[1] = 0;
   for (int i = 1; i < list.size(); i++)
         ig[i + 1] = ig[i] + list[i].size();

   // Формируем массив jg
   jg.resize(ig.back());
   for (int i = 1, j = 0; i < n_edges; i++)
      for (int k = 0; k < list[i].size(); k++)
         jg[j++] = list[i][k];
}
// Проверка списка на упорядоченность по возрастанию
bool FEM::isOrdered(const std::vector<int>& v)
{
    if (v.size() > 0)
    {
        for (int i = 0; i < v.size() - 1; i++)
            if (v[i + 1] < v[i])
                return false;
    }
   return true;
}

// Получить значение локальной матрицы жесткости
void FEM::getLocalG(double hx, double hy, double hz)
{
    double hy_hz_6hx = hy * hz / (hx * 6.0);
    double hx_hz_6hy = hx * hz / (hy * 6.0);
    double hx_hy_6hz = hx * hy / (hz * 6.0);
    double hx_6 = hx / 6.0;
    double hy_6 = hy / 6.0;
    double hz_6 = hz / 6.0;
    double coef = 1.0 / mu;

    for(int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            G_loc[i][j] = coef * (hx_hy_6hz * G1[i][j] + hx_hz_6hy * G2[i][j]);

    for (int i = 0; i < 4; i++)
        for (int j = 4; j < 8; j++)
            G_loc[i][j] = coef * (-hz_6) * G2[i][j - 4];

    for (int i = 0; i < 4; i++)
        for (int j = 8; j < 12; j++)
            G_loc[i][j] = coef * (hy_6)*G3[i][j - 8];

    for (int i = 4; i < 8; i++)
        for (int j = 0; j < 4; j++)
            G_loc[i][j] = coef * (-hz_6) * G2[i - 4][j];

    for (int i = 4; i < 8; i++)
        for (int j = 4; j < 8; j++)
            G_loc[i][j] = coef * (hx_hy_6hz * G1[i - 4][j - 4] + hy_hz_6hx * G2[i - 4][j - 4]);

    for (int i = 4; i < 8; i++)
        for (int j = 8; j < 12; j++)
            G_loc[i][j] = coef * (-hx_6) * G1[i - 4][j - 8];

    for (int i = 8; i < 12; i++)
        for (int j = 0; j < 4; j++)
            G_loc[i][j] = coef * hy_6 * G3_T[i - 8][j];

    for (int i = 8; i < 12; i++)
        for (int j = 4; j < 8; j++)
            G_loc[i][j] = coef * (-hx_6) * G1[i - 8][j - 4];

    for (int i = 8; i < 12; i++)
        for (int j = 8; j < 12; j++)
            G_loc[i][j] = coef * (hx_hz_6hy * G1[i - 8][j - 8] + hy_hz_6hx * G2[i - 8][j - 8]);

}

void FEM::getLocalM(double hx, double hy, double hz)
{
    double coef = gamma * (hx * hy * hz) / 36.0;
    for (int i = 0; i < M_loc.size(); i++)
        for (int j = 0; j < M_loc[i].size(); j++)
            M_loc[i][j] = coef * M[i][j];
}

void FEM::getLocalRightPart(point p, double hx, double hy, double hz, int ielem, int cur_tLayer)
{
    point tmp;
    int k = 0;
    F_loc.clear();
    F_loc.resize(12);
    q0.resize(12);
    q1.resize(12);

    std::vector<double> F_loc_tmp = F_loc;

    // b^j
    for (int i = 0; i < F_loc_tmp.size(); i++)
        F_loc_tmp[i] = RP_value(_elems[ielem].edges[i]) + gamma * dAdt(_elems[ielem].edges[i], time_layer[cur_tLayer]);

    // q^(j-1)
    for (int i = 0; i < q0.size(); i++)
        q0[i] = q[_elems[ielem].edges[i]];

    if (cur_tLayer > 0)
        dt = 1.0 / (time_layer[cur_tLayer] - time_layer[cur_tLayer - 1]);

    double coef = (hx * hy * hz) / 36.0;

    //b^j + 1/delT * M * q^(j-1)
    for (int i = 0; i < M.size(); i++)
        for (int j = 0; j < M.size(); j++)
            F_loc[i] += coef * M[i][j] * F_loc_tmp[j] + dt * M_loc[i][j] * q0[j];
}

// Добавить элемент в глобальную матрицу
void FEM::addElementToGlobal(int i, int j, double elem)
{
   if (i == j)
   {
      di[i] += elem;
      return;
   }
   else
   {
      for (int ind = ig[i]; ind < ig[i + 1]; ind++)
         if (jg[ind] == j)
         {
            gg[ind] += elem;
            return;
         }
   }
}

// Сборка глобальной матрицы
void FEM::assemblyGlobalMatrix(int cur_tLayer)
{
   gg.clear();
   gg.resize(ig.back(), 0);
   M_loc.resize(12);
   G_loc.resize(12);
   for (int i = 0; i < 12; i++)
   {
       M_loc[i].resize(12);
       G_loc[i].resize(12);
   }
   for (int ielem = 0; ielem < nFE; ielem++)
   {
      point p1 = g.MeshXYZ[_elems[ielem].nodes[0]];
      point p2 = g.MeshXYZ[_elems[ielem].nodes[7]];
      double hx = p2.x - p1.x;
      double hy = p2.y - p1.y;
      double hz = p2.z - p1.z;
      getLocalG(hx, hy, hz);
      getLocalM(hx, hy, hz);
      getLocalRightPart(p1, hx, hy, hz, ielem, cur_tLayer);
      for (int i = 0; i < _elems[ielem].edges.size(); i++)
      {
         for (int j = 0; j < _elems[ielem].edges.size(); j++)
         {
            addElementToGlobal(_elems[ielem].edges[i], _elems[ielem].edges[j], G_loc[i][j]);
            addElementToGlobal(_elems[ielem].edges[i], _elems[ielem].edges[j], M_loc[i][j]);
         }

         F[_elems[ielem].edges[i]] += F_loc[i];
      }
   }
}

// Установка первых краевых условий
void FEM::bc1()
{
    //Левая грань
   int gr = nx * (ny - 1) + ny * (nx - 1);
   int pop = nx * ny;

    for(int i = 0; i < nz - 1; i++)
       for (int j = 0; j < ny - 1; j++)
       {

          std::vector<int> n(4, 0);

          n[0] = i * (gr + pop) + j * (2 * nx - 1) + (nx - 1);
          n[1] = i * (gr + pop) + gr + j * nx;
          n[2] = i * (gr + pop) + gr + j * nx + nx;
          n[3] = (i + 1) * (gr + pop) + j * (2 * nx - 1) + (nx - 1);

          for (int i = 0; i < 4; i+=3)
          {
             if (!isInBoundary(n[i]))
             {
                _bc1.push_back({ {n[i]}, {bc1_value(n[i], 0)}});
             }
          }

          for (int i = 1; i < 3; i ++)
          {
             if (!isInBoundary(n[i]))
             {
                _bc1.push_back({ {n[i]}, {bc1_value(n[i], 0)} });
             }
          }
       }

    //Правая грань

    for (int i = 0; i < nz - 1; i++)
       for (int j = 0; j < ny - 1; j++)
       {
             std::vector<int> n(4, 0);

             n[0] = i * (gr + pop) + j * (2 * nx - 1) + (nx - 2) + 1 + (nx - 1);
             n[1] = i * (gr + pop) + gr + (nx - 2) + j * nx + 1;
             n[2] = i * (gr + pop) + gr + (nx - 2) + j * nx + nx + 1;
             n[3] = (i + 1) * (gr + pop) + j * (2 * nx - 1) + (nx - 1) + (nx - 2) + 1;

          for (int i = 0; i < 4; i += 3)
          {
             if (!isInBoundary(n[i]))
             {
                _bc1.push_back({ {n[i]}, {bc1_value(n[i], 0)} });
             }
          }

          for (int i = 1; i < 3; i++)
          {
             if (!isInBoundary(n[i]))
             {
                _bc1.push_back({ {n[i]}, {bc1_value(n[i], 0)} });
             }
          }
       }


    //Нижняя грань

    for (int i = 0; i < ny - 1; i++)
       for (int j = 0; j < nx - 1; j++)
       {
          std::vector<int> n(4, 0);

          n[0] = i * (2 * nx - 1) + j;
          n[1] = i * (2 * nx - 1) + j + (nx - 1);
          n[2] = i * (2 * nx - 1) + j + (nx - 1) + 1;
          n[3] = (i + 1) * (2 * nx - 1) + j;

          for (int i = 0; i < 4; i += 3)
          {
             if (!isInBoundary(n[i]))
             {
                _bc1.push_back({ {n[i]}, {bc1_value(n[i], 0)} });
             }
          }

          for (int i = 1; i < 3; i++)
          {
             if (!isInBoundary(n[i]))
             {
                _bc1.push_back({ {n[i]}, {bc1_value(n[i], 0)} });
             }
          }
       }


    //Верхняя грань

    for (int i = 0; i < ny - 1; i++)
       for (int j = 0; j < nx - 1; j++)
       {
          std::vector<int> n(4, 0);

          n[0] = (nz - 1) * (gr + pop) + i * (2*nx-1) + j;
          n[1] = (nz - 1) * (gr + pop) + i * (2*nx-1) + (nx - 1) + j;
          n[2] = (nz - 1) * (gr + pop) + i * (2*nx-1) + (nx - 1) + j + 1;
          n[3] = (nz - 1) * (gr + pop) + (i+1) * (2*nx-1) + j;

          for (int i = 0; i < 4; i += 3)
          {
             if (!isInBoundary(n[i]))
             {
                _bc1.push_back({ {n[i]}, {bc1_value(n[i], 0)} });
             }
          }

          for (int i = 1; i < 3; i++)
          {
             if (!isInBoundary(n[i]))
             {
                _bc1.push_back({ {n[i]}, {bc1_value(n[i], 0)} });
             }
          }
       }


    //Передняя грань

    for (int i = 0; i < nz - 1; i++)
       for (int j = 0; j < nx - 1; j++)
       {
          std::vector<int> n(4, 0);

          n[0] = i * (gr + pop) + j;
          n[1] = i * (gr + pop) + gr + j;
          n[2] = i * (gr + pop) + gr + j + 1;
          n[3] = (i + 1) * (gr + pop) + j;

          for (int i = 0; i < 4; i += 3)
          {
             if (!isInBoundary(n[i]))
             {
                _bc1.push_back({ {n[i]}, {bc1_value(n[i], 0)} });
             }
          }

          for (int i = 1; i < 3; i++)
          {
             if (!isInBoundary(n[i]))
             {
                _bc1.push_back({ {n[i]}, {bc1_value(n[i], 0)} });
             }
          }
       }

    //Задняя грань

    for (int i = 0; i < nz - 1; i++)
       for(int j = 0; j < nx - 1; j++)
       {
          std::vector<int> n(4, 0);

          n[0] = i * (gr + pop) + ((ny - 2) + 1) * (2 * nx - 1) + j;
          n[1] = i * (gr + pop) + gr + j + (ny - 2) * nx + nx;
          n[2] = i * (gr + pop) + gr + j + (ny - 2) * nx + nx + 1;
          n[3] = (i + 1) * (gr + pop) + ((ny - 2) + 1) * (2 * nx - 1) + j; 

          for (int i = 0; i < 4; i += 3)
          {
             if (!isInBoundary(n[i]))
             {
                _bc1.push_back({ {n[i]}, {bc1_value(n[i], 0)} });
             }
          }

          for (int i = 1; i < 3; i++)
          {
             if (!isInBoundary(n[i]))
             {
                _bc1.push_back({ {n[i]}, {bc1_value(n[i], 0)} });
             }
          }
       }

    std::sort(_bc1.begin(), _bc1.end(), [](const bc& b1, const bc& b2) -> bool
       {
          return b1.node < b2.node;
       }
    );

    std::vector<int> bc1_nodes(n_edges, -1);
    for (int i = 0; i < _bc1.size(); i++)
    {
       bc1_nodes[_bc1[i].node] = i; // В узле задано краевое
    }
    int k;
    for (int i = 0; i < n_edges; i++)
    {
       if (bc1_nodes[i] != -1) //задано
       {
          di[i] = 1.0;
          F[i] = _bc1[bc1_nodes[i]].value;
          for (int j = ig[i]; j < ig[i + 1]; j++)
          {
             k = jg[j];
             if (bc1_nodes[k] == -1)
             {
                F[k] -= gg[j] * F[i];
             }
             gg[j] = 0.0;
          }
       }
       else //не задано
       {
          for (int j = ig[i]; j < ig[i + 1]; j++)
          {
             k = jg[j];
             if (bc1_nodes[k] != -1)
             {
                F[i] -= gg[j] * F[k];
                gg[j] = 0.0;
             }
          }
       }
    }
}

int FEM::globalToLocalEdge(int ielem, int gl_num)
{
   for (int i = 0; i < 12; i++)
      if (gl_num == _elems[ielem].edges[i])
         return i;
}

int FEM::globalToLocalNode(int ielem, int gl_num)
{
   for (int i = 0; i < 8; i++)
      if (gl_num == _elems[ielem].nodes[i])
         return i;
}

bool FEM::isInBoundary(int num)
{
   for (int i = 0; i < _bc1.size(); i++)
      if (_bc1[i].node == num)
         return true;
   return false;
}

void FEM::CheckSol(int cur_tLayer)
{
   std::ofstream out("files/Solution"+ std::to_string(cur_tLayer) +".txt");
   std::vector<double> q_true;

   out << "iedge" << std::setw(5) << "q" << std::setw(21) << "q*" << std::setw(27) << "|q - q*|" << std::setw(29) << "||q - q*||/||q*||" << std::endl;
   
   for(int iedge = 0; iedge < n_edges; iedge++)
      q_true.push_back(bc1_value(iedge, cur_tLayer));

   std::vector<double> pogr = q;

   for (int i = 0; i < pogr.size(); i++)
      pogr[i] = abs(q[i] - q_true[i]);

   double norm = 0;
   for (int i = 0; i < pogr.size(); i++)
      norm += pogr[i] * pogr[i];

   double norm_true = 0;
   for (int i = 0; i < pogr.size(); i++)
      norm_true += q_true[i] * q_true[i];

   //out << std::scientific << sqrt(norm) / sqrt(norm_true);

   out << std::scientific << 0 << std::setw(20) << q[0] << std::setw(20) << q_true[0] << std::setw(20) << pogr[0] << std::setw(20) << sqrt(norm) / sqrt(norm_true) << std::endl;

   for(int i = 1; i < 10; i++)
      out << std::scientific << i << std::setw(20) << q[i] << std::setw(20) << q_true[i] << std::setw(20) << pogr[i] << std::endl;

   for (int i = 11; i < q.size(); i++)
      out << std::scientific << i << std::setw(19) << q[i] << std::setw(20) << q_true[i] << std::setw(20) << pogr[i] << std::endl;
}
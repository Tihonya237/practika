#pragma once
#ifndef RESOURCE_H
#define RESOURCE_H

#include <vector>
#include <functional>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <locale.h>
#include <math.h>
#include "../KTMIAD1/nlohmann/json.hpp"

const std::string directory = "../";

struct point
{
   double x, y, z;
};

#endif

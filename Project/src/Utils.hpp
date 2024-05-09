#pragma once
#include <iostream>
#include "Eigen/Eigen"
#include <Data.hpp>

using namespace std;
using namespace Eigen;

namespace FracturesTraces
{
bool ImportData(const string &filepath, DFN& data);

bool ImportAll(const string &filename, DFN& data);

}

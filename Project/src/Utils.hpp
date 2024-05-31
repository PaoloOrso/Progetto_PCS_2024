#pragma once
#include <iostream>
#include "Eigen/Eigen"
#include <Data.hpp>


using namespace std;
using namespace Eigen;

namespace FracturesTraces
{
bool FinalTest(DFN& data);

bool ImportAll(const string &filename, DFN& data);

bool CreateSpheres(DFN &data);

bool CreateNormals(DFN &data);

bool FindIntersections(DFN &data);

bool PrintResults(DFN &data);

bool SubPolygons(DFN &data);

}

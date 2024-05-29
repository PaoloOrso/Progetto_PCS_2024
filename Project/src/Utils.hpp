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

bool Testsfera(DFN &data);

bool Testpianiparalleli(DFN &data);

bool Testintersezione(DFN &data);

bool Stampa(DFN &data);


}

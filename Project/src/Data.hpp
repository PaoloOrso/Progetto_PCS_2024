#pragma once
#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace FracturesTraces{

struct DFN
{
    unsigned int NumberFractures = {};
    vector<unsigned int> FracturesId = {};
    vector<unsigned int> NumberVertices = {};
    map<unsigned int, Matrix<double, 3, Dynamic>> Vertices ={};


    unsigned int NumberTraces = 0;
    vector<unsigned int> TracesId = {};
    map<unsigned int, vector<unsigned int>> GeneratingFractures = {};
    map<unsigned int, Matrix<double, 3,2>> GeneratingVertices = {};

    vector<bool> Tips = {};
    vector<double> LenghtTraces = {};

};

}

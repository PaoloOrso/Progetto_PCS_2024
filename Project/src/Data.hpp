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
    map<unsigned int, vector<vector<double>>> Vertices ={};


    unsigned int NumberTraces = 0;
    vector<unsigned int> TracesId = {};
    map<unsigned int, vector<unsigned int>> GeneratingFractures = {};
    map<unsigned int, vector<vector<double>>> GeneratingVertices = {};

    vector<bool> Tips = {};
    vector<double> LenghtTraces = {};

};

}

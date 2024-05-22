#pragma once
#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace FracturesTraces{

struct DFN
{
    unsigned int NumberFractures = 0;
    vector<unsigned int> FracturesId = {};
    vector<unsigned int> NumberVertices = {};
    map<unsigned int, vector<vector<double>>> Vertices ={};

    map<unsigned int, double> raggi = {};

    map<unsigned int, vector<double>> Normals = {};

    map<unsigned int, double> Directors;



    unsigned int NumberTraces = 0;
    vector<unsigned int> TracesId = {};
    map<unsigned int, vector<unsigned int>> GeneratingFractures = {};
    map<unsigned int, vector<vector<double>>> GeneratingVertices = {};


    map<unsigned int, vector<double>> Baricentri = {};


    vector<bool> Tips = {};
    vector<double> LenghtTraces = {};

};

}

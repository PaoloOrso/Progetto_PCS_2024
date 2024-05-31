#pragma once
#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace FracturesTraces{

struct DFN
{
    double Tol = 1e-10;

    unsigned int N_Fractures;
    unsigned int MaxId;
    vector<unsigned int> N_Vertices;
    vector<vector<Vector3d>> Vertices;

    vector<Vector3d> Barycentres;
    vector<double> Rays;

    vector<Vector3d> Normals;  
    vector<double> Ds;

    unsigned int NumberTraces;
    vector<Vector2i> GeneratingFractures;
    vector<vector<Vector3d>> GeneratingPoints;
    vector<double> LenghtTraces;
    vector<vector<bool>> Tips;
    vector<unsigned int> TracesinFigures;
    vector<vector<unsigned int>> IdTraces;
    vector<vector<bool>> BoolTraces;
    vector<vector<pair<unsigned int,double>>> Id_Lenght_Fractures;

};    

}

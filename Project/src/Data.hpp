#pragma once
#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace FracturesTraces{

struct DFN
{
    unsigned int NumberFractures;
    unsigned int MaxId;
    vector<unsigned int> NumberVertices;
    map<unsigned int, vector<Vector3d>> Vertices;

    vector<Vector3d> Baricentri;
    vector<double> raggi;

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

    map<unsigned int, double> Idpassanti;
    map<unsigned int, double> Idnonpassanti;

   // vector<pair<unsigned int,double>> Id_Lenght;
    vector<vector<pair<unsigned int,double>>> Id_Lenght_Fractures;




};    

}

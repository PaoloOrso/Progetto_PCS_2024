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

    vector<vector<unsigned int>> Id_Traces_Sorted;
    vector<vector<bool>> BoolTraces_Sorted;




    vector<unsigned int> Number_Points; //vettore numero dei punti per ogni poligono
    vector<vector<unsigned int>> Id_Points; //vettore di vettori degli id dei punti
    vector<vector<Vector3d>> Id_Points_Coordinates; //vettore di vettori di coordinate dei punti

    vector<unsigned int> Number_Segments; // vettore di numero di segmenti
    vector<vector<unsigned int>> Id_Segments; // vettore di vettori degli dei segmenti
    vector<vector<Vector2i>> Id_Points_Extremes; // vettore di vettori degli id degli estremi

    vector<unsigned int> Number_Polygons; // vettore di numero di celle
    vector<vector<unsigned int>> Id_Polygons; //vettore di vettori di id delle celle
    vector<vector<vector<unsigned int>>> Id_Polygons_Vertices; // vettore di vettori degli id dei punti delle celle
    vector<vector<vector<unsigned int>>> Id_Polygons_Edges;  // vettore di vettori degli id dei segmenti delle celle
    vector<vector<unsigned int>> NumberVerticesEges;



};

}



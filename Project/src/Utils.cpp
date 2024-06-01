// 3 10 50 82 200 362
// 2 25 481 1 8985 1
#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <Eigen/Eigen>
#include <algorithm>

using namespace Eigen;

namespace FracturesTraces
{

//-------CONTROLLO-FINALE---------------------------------------------------------------------------

bool FinalTest(DFN& data)
{
    // 3 10 50 82 200 362
    if(!ImportAll("./FR10_data.txt", data))
    {
        return false;
    }

    if(!CreateSpheres(data))
    {
        return false;
    }

    if(!CreateNormals(data))
    {
        return false;
    }

    if(!FindIntersections(data))
    {
        return false;
    }

    if(!PrintResults(data))
    {
        return false;
    }

    if(!SubPolygons(data))
    {
        return false;
    }

    return true;
}

//-------INIZIALIZZO-N.FRATTURE--MAX.ID--VETTORE-DI-N.VERTICI--VETTORE-DEI-VERTICI------------------

bool ImportAll(const string &filename, DFN &data)
{
    ifstream file;
    file.open(filename);

    if(file.fail())
    {
        cerr << "errore" << endl;
        return false;
    }

    string line;
    unsigned int NFractures = 0;
    istringstream convertN(line);

    getline(file, line);
    getline(file, line);
    convertN.str(line);
    convertN >> NFractures;
    data.N_Fractures = NFractures;
    getline(file,line);

    for(unsigned int i = 0; i != NFractures; i++)
    {
        istringstream convertN(line);

        unsigned int id, vertices;

        char tmp;

        double coordx, coordy, coordz;

        vector<double> coordinatesx, coordinatesy, coordinatesz;

        Vector3d vert;

        vector<Vector3d> Vertices;

        getline(file, line);
        convertN.str(line);
        convertN >> id >> tmp >> vertices;

        data.N_Vertices.push_back(vertices);

        getline(file, line);

        getline(file, line);
        replace(line.begin(),line.end(), ';' ,' ');
        istringstream convert1(line);

        for(unsigned int j = 0; j != vertices; j++)
        {
            convert1 >> coordx;
            coordinatesx.push_back(coordx);
        }

        getline(file, line);
        replace(line.begin(),line.end(), ';' ,' ');
        istringstream convert2(line);

        for(unsigned int j = 0; j != vertices; j++)
        {
            convert2 >> coordy;
            coordinatesy.push_back(coordy);
        }

        getline(file, line);
        replace(line.begin(),line.end(), ';' ,' ');
        istringstream convert3(line);

        for(unsigned int j = 0; j != vertices; j++)
        {
            convert3 >> coordz;
            coordinatesz.push_back(coordz);
        }

        for(unsigned int i = 0; i != vertices;i++)
        {
            vert[0] = coordinatesx[i];
            vert[1] = coordinatesy[i];
            vert[2] = coordinatesz[i];

            Vertices.push_back(vert);
            vert = {};
        }

        data.MaxId = id;
        data.Vertices.push_back(Vertices);
        getline(file, line);
    }

    return true;
}

//-------INIZIALIZZO-BARICENTRI-E-RAGGI-------------------------------------------------------------

bool CreateSpheres(DFN &data)
{
    for(unsigned int id = 0; id != (data.MaxId + 1); id++)
    {
        Vector3d bari;

        double tmpx = 0.0, tmpy = 0.0, tmpz = 0.0;

        double distance = 0.0, maxdistance = 0.0;

        for(unsigned int i = 0; i != data.N_Vertices[id]; i++)
        {
            Vector3d point = data.Vertices[id][i];

            tmpx += point[0];
            tmpy += point[1];
            tmpz += point[2];
        }

        bari[0] = tmpx / data.N_Vertices[id];
        bari[1] = tmpy / data.N_Vertices[id];
        bari[2] = tmpz / data.N_Vertices[id];

        data.Barycentres.push_back(bari);

        for(unsigned int i = 0; i != data.N_Vertices[id]; i++)
        {
            Vector3d point = data.Vertices[id][i];

            distance = sqrt(pow(bari[0] - point[0],2) + pow(bari[1] - point[1],2) + pow(bari[2] - point[2],2));
            if(distance > maxdistance)
            {
                maxdistance = distance;
            }
        }

        data.Rays.push_back(maxdistance);

    }

    return true;

}

//-------INIZIALIZZO-COEFFICIENTI-DEI-PIANI-E-VERSORI-NORMALI---------------------------------------

bool CreateNormals(DFN &data)
{   
    Vector3d point0, point1, point2;

    Vector3d vector1, vector2;
    
    for(unsigned int id = 0; id != (data.MaxId + 1); id++)
    {
        Vector3d normal;

        double d;

        point0 = data.Vertices[id][0];
        point1 = data.Vertices[id][1];
        point2 = data.Vertices[id][2];

        vector1 = point2-point0;
        vector2 = point1-point0;

        normal = (vector1.cross(vector2)).normalized();

        d = -normal[0]*point0[0] - normal[1]*point0[1] - normal[2]*point0[2];

        data.Normals.push_back(normal);
        data.Ds.push_back(d);
    }

    return true;

}

//-------TROVO-INTERSEZIONI-------------------------------------------------------------------------

bool FindIntersections(DFN &data)
{

    for(unsigned int k = 0; k!= data.N_Fractures; k++)
    {
        data.IdTraces.push_back({});
        data.BoolTraces.push_back({});
        data.Id_Lenght_Fractures_Tips.push_back({});
    }

    unsigned int N_Intersections = 0;

    Vector3d bari1, bari2;

    Vector3d normal1, normal2;

    double dist_bari, sum_rays;

    Vector3d director;

    vector<unsigned int> N_traces_fractures(data.N_Fractures);

    Vector3d point0_trace, point1_trace;

    Vector3d vert0, vert1;

    Vector3d diff1, diff2, diff3;

    Vector3d det;

    double alfa, beta;

    Vector2d alfa_beta;

    double trace_lenght, delta_alfa;

    Vector2i gen_frac;

    vector<Vector3d> gen_points;

    Vector3d trace_end1, trace_end2;

    vector<double> tmp;

    bool cond1, cond2, cond3, cond4;

    for (unsigned int i = 0; i != data.N_Fractures; i++)
        for ( unsigned int j = i+1; j != data.N_Fractures; j++)
        {

            bari1 = data.Barycentres[i];
            bari2 = data.Barycentres[j];

            normal1 = data.Normals[i];
            normal2 = data.Normals[j];

            dist_bari = sqrt(pow(bari1[0]-bari2[0],2) + pow(bari1[1]-bari2[1],2) + pow(bari1[2]-bari2[2],2));
            sum_rays = data.Rays[i] + data.Rays[j];

            if(( dist_bari - sum_rays < data.Tol ) && ( (normal1.cross(normal2)).norm() > data.Tol ))
            {

                vector<double> alfa_vector;

                director = normal1.cross(normal2);

                Matrix3d A;

                A << normal1[0], normal1[1], normal1[2],
                     normal2[0], normal2[1], normal2[2],
                     director[0], director[1], director[2];

                Vector3d b;

                b << -data.Ds[i], -data.Ds[j], 0;

                point0_trace = A.colPivHouseholderQr().solve(b);

                point1_trace = point0_trace + director;

                for(unsigned int w = 0; w < data.N_Vertices[i]; w++)
                {
                    vert0 = data.Vertices[i][w];
                    vert1 = data.Vertices[i][ (w + 1) % data.N_Vertices[i] ];

                    diff1 = point1_trace-point0_trace;
                    diff2 = vert1-vert0;
                    diff3 = vert0-point0_trace;

                    det = diff1.cross(diff2);

                    if(det.norm() > data.Tol)
                    {
                        Matrix<double,3,2> AA;

                        AA << diff1[0], diff2[0],
                              diff1[1], diff2[1],
                              diff1[2], diff2[2];


                        Vector3d bb;

                        bb(0) = diff3[0]; bb(1) = diff3[1]; bb(2) = diff3[2];

                        alfa_beta = AA.colPivHouseholderQr().solve(bb);

                        alfa = alfa_beta(0);
                        beta = alfa_beta(1);

                        if(-beta > 0.0 && -beta < 1.0)
                        {
                            alfa_vector.push_back(alfa);
                        }
                    }
                }

                for(unsigned int w = 0; w < data.N_Vertices[j]; w++)
                {
                    vert0 = data.Vertices[j][w];
                    vert1 = data.Vertices[j][(w + 1) % data.N_Vertices[j] ];

                    diff1 = point1_trace-point0_trace;
                    diff2 = vert1-vert0;
                    diff3 = vert0-point0_trace;

                    det = diff1.cross(diff2);

                    if(det.norm() > data.Tol)
                    {
                        Matrix<double,3,2> AA;

                        AA << diff1[0], diff2[0],
                              diff1[1], diff2[1],
                              diff1[2], diff2[2];


                        Vector3d bb;

                        bb(0) = diff3[0]; bb(1) = diff3[1]; bb(2) = diff3[2];

                        alfa_beta = AA.colPivHouseholderQr().solve(bb);

                        alfa = alfa_beta(0);
                        beta = alfa_beta(1);

                        if(-beta > 0.0 && -beta < 1.0)
                        {
                            alfa_vector.push_back(alfa);
                        }
                    }
                }

                if( size(alfa_vector) == 4)
                {
                    tmp = alfa_vector;

                    sort(tmp.begin(),tmp.begin()+4);

                    delta_alfa = abs(tmp[1]-tmp[2]);

                    trace_lenght = delta_alfa * director.norm();

                    gen_frac = {i,j};

                    trace_end1 = point0_trace + alfa_vector[1]*(point1_trace-point0_trace);
                    trace_end2 = point0_trace + alfa_vector[2]*(point1_trace-point0_trace);

                    gen_points = {trace_end1,trace_end2};

                    tuple<unsigned int, double, bool> Id_trace_lenght_tips_false (N_Intersections,trace_lenght,false);
                    tuple<unsigned int, double, bool> Id_trace_lenght_tips_true (N_Intersections,trace_lenght,true);

                    cond1 = abs(alfa_vector[0] - alfa_vector[2]) < data.Tol && abs(alfa_vector[1] - alfa_vector[3]) < data.Tol;

                    cond2 = ((max(alfa_vector[0],alfa_vector[1]) >= max(alfa_vector[2],alfa_vector[3])) && (min(alfa_vector[0],alfa_vector[1]) <= min(alfa_vector[2],alfa_vector[3]))) ||
                            ((max(alfa_vector[2],alfa_vector[3]) >= max(alfa_vector[0],alfa_vector[1])) && (min(alfa_vector[2],alfa_vector[3]) <= min(alfa_vector[0],alfa_vector[1])));

                    cond3 = ((max(alfa_vector[0],alfa_vector[1]) > min(alfa_vector[2],alfa_vector[3])) && (min(alfa_vector[0],alfa_vector[1]) < min(alfa_vector[2],alfa_vector[3]))) ||
                            ((max(alfa_vector[2],alfa_vector[3]) > min(alfa_vector[0],alfa_vector[1])) && (min(alfa_vector[2],alfa_vector[3]) < min(alfa_vector[0],alfa_vector[1])));

                    cond4 = (alfa_vector[0] >= min(alfa_vector[2],alfa_vector[3]) && alfa_vector[0] <= max(alfa_vector[2],alfa_vector[3]))  &&
                            (alfa_vector[1] >= min(alfa_vector[2],alfa_vector[3]) && alfa_vector[1] <= max(alfa_vector[2],alfa_vector[3]));

                    if (cond1 || cond2 || cond3)
                    {                       
                        data.GeneratingFractures.push_back(gen_frac);
                        data.GeneratingPoints.push_back(gen_points);
                        data.LenghtTraces.push_back(trace_lenght);

                        N_traces_fractures[i]++;
                        N_traces_fractures[j]++;

                        data.IdTraces[i].push_back(N_Intersections);
                        data.IdTraces[j].push_back(N_Intersections);

                        N_Intersections++;

                        if(cond1)
                        {                         
                            data.Id_Lenght_Fractures_Tips[i].push_back(Id_trace_lenght_tips_true);
                            data.Id_Lenght_Fractures_Tips[j].push_back(Id_trace_lenght_tips_true);

                            data.Tips.push_back({true,true});
                            data.BoolTraces[i].push_back(true);
                            data.BoolTraces[j].push_back(true);

                        }
                        else if(cond2)
                        {
                            if(cond4)
                            {
                                data.Id_Lenght_Fractures_Tips[i].push_back(Id_trace_lenght_tips_true);
                                data.Id_Lenght_Fractures_Tips[j].push_back(Id_trace_lenght_tips_false);

                                data.Tips.push_back({true,false});
                                data.BoolTraces[i].push_back(true);
                                data.BoolTraces[j].push_back(false);
                            }
                            else
                            {
                                data.Id_Lenght_Fractures_Tips[i].push_back(Id_trace_lenght_tips_false);
                                data.Id_Lenght_Fractures_Tips[j].push_back(Id_trace_lenght_tips_true);

                                data.Tips.push_back({false,true});
                                data.BoolTraces[i].push_back(false);
                                data.BoolTraces[j].push_back(true);
                            }
                        }
                        else if(cond3)
                        {
                            data.Id_Lenght_Fractures_Tips[i].push_back(Id_trace_lenght_tips_false);
                            data.Id_Lenght_Fractures_Tips[j].push_back(Id_trace_lenght_tips_false);

                            data.Tips.push_back({false,false});
                            data.BoolTraces[i].push_back(false);
                            data.BoolTraces[j].push_back(false);
                        }
                    }
                }
            }
        }

    data.NumberTraces = N_Intersections;
    data.TracesinFigures = N_traces_fractures;

    for(unsigned int k = 0; k!= data.N_Fractures;k++)
    {
        sort(data.Id_Lenght_Fractures_Tips[k].begin(), data.Id_Lenght_Fractures_Tips[k].end(), [](const tuple<unsigned int, double,bool>& a, const tuple<unsigned int,double,bool>& b)
        {
            return get<1>(a) > get<1>(b);
        });
    }

    return true;
}

//-------STAMPA-RISULTATI---------------------------------------------------------------------------

bool PrintResults(DFN &data)
{
    ofstream out("result.txt");

    out << "# Number of Traces" << endl << data.NumberTraces << endl << endl;
    for(unsigned int t = 0; t != data.NumberTraces; t++)
    {
        out << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
        out << t << " ; " << data.GeneratingFractures[t][0] << " , " << data.GeneratingFractures[t][1] << " ; ";
        out << data.GeneratingPoints[t][0][0] << " , " << data.GeneratingPoints[t][0][1] << " , " << data.GeneratingPoints[t][0][2] << " ; ";
        out << data.GeneratingPoints[t][1][0] << " , " << data.GeneratingPoints[t][1][1] << " , " << data.GeneratingPoints[t][1][2] << endl;
    }
    out << endl;
    for(unsigned int id = 0; id != data.N_Fractures; id++)
    {
        vector<unsigned int> Sorted_ids;

        if(size(data.IdTraces[id]) != 0)
        {
        out << "# FractureId; NumTraces" << endl << id << " ; " << data.TracesinFigures[id] << endl;
        out << "# TraceId; Tips; Lenght" << endl;
        for(unsigned int l = 0; l != size(data.Id_Lenght_Fractures_Tips[id]); l++)
        {
            if(data.BoolTraces[id][l] == true)
            {
                out << get<0>(data.Id_Lenght_Fractures_Tips[id][l]) << " ; " << boolalpha << data.BoolTraces[id][l] << noboolalpha << " ; " << get<1>(data.Id_Lenght_Fractures_Tips[id][l]) << endl;
                Sorted_ids.push_back(get<1>(data.Id_Lenght_Fractures_Tips[id][l]));
            }
        }

        for(unsigned int l = 0; l != size(data.Id_Lenght_Fractures_Tips[id]); l++)
        {
            if(data.BoolTraces[id][l] == false)
            {
                out << get<0>(data.Id_Lenght_Fractures_Tips[id][l]) << " ; " << boolalpha << data.BoolTraces[id][l] << noboolalpha << " ; " << get<1>(data.Id_Lenght_Fractures_Tips[id][l]) << endl;
                Sorted_ids.push_back(get<1>(data.Id_Lenght_Fractures_Tips[id][l]));
            }
        }
        out << endl;
        }

        data.Id_Traces_Sorted.push_back(Sorted_ids);
    }
    out.close();

    return true;
}

//-------CREO-SOTTOPOLIGONI-------------------------------------------------------------------------

bool SubPolygons(DFN &data)
{
    vector<Vector3d> Polygon;
    vector<vector<Vector3d> Polygroup;
    unsigned int N_traces;
    vector<unsigned int> Id_traces_sorted;
    vector<Vector3d> Trace_points;
    Vector3d ver0, ver1;
    Vector3d dif1, dif2, dif3;
    Vector2d coeffs;
    Vector3d dett;
    Vector3d new_vert0, new_vert1;
    double coef0, coef1;


    for(unsigned int i = 0; i != data.N_Fractures; i++)
    {
        Polygon = data.Vertices[i];

        N_traces = data.TracesinFigures[i];

        Id_traces_sorted = data.Id_Traces_Sorted[i];

        for(unsigned int j = 0; j != N_traces ; j++)  // qui inizio con i tagli
        {
            Trace_points = data.GeneratingPoints[Id_traces_sorted[j]];  // punti generatori della traccia

            for(unsigned int w = 0; w < data.N_Vertices[i]; w++)
            {
                vector<double> prova;

                //se la frattura è passante, i punti generatori sono i nuovi vertici

                if(get<2>data.Id_Lenght_Fractures_Tips[i][j] == true)
                {
                    new_vert0 = data.GeneratingPoints[Id_traces_sorted][0];
                    new_vert1 = data.GeneratingPoints[Id_traces_sorted][1];

                    // devo capire che lati le tracce dividono





                }

                ver0 = data.Vertices[i][w];
                ver1 = data.Vertices[i][ (w + 1) % data.N_Vertices[i] ];

                dif1 = Trace_points[1]-Trace_points[0];
                dif2 = ver1-ver0;
                dif3 = ver0-Trace_points[0];

                dett = dif1.cross(dif2);

                if(dett.norm() > data.Tol)
                {
                    Matrix<double,3,2> E;

                    E << dif1[0], dif2[0],
                         dif1[1], dif2[1],
                         dif1[2], dif2[2];


                    Vector3d f;

                    f(0) = dif3[0]; f(1) = dif3[1]; f(2) = dif3[2];

                    coeffs = E.colPivHouseholderQr().solve(f);  // trovo coeffs punto int traccia-lati

                    coef0 = coeffs(0);
                    coef1 = coeffs(1);

                    if(-coef1 > 0.0 && -coef1 < 1.0)
                    {
                        prova.push_back(coef0);
                    }

                }


            }

        }

    }

    return true;

}

}

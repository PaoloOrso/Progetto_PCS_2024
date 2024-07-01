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
    if(!ImportAll("./FR3_data.txt", data))
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
    unsigned int id, vertices;
    char tmp;
    double coordx, coordy, coordz;

    getline(file, line);
    getline(file, line);
    convertN.str(line);
    convertN >> NFractures;
    data.N_Fractures = NFractures;
    getline(file,line);

    for(unsigned int i = 0; i != NFractures; i++)
    {
        istringstream convertN(line);

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

        data.Vertices.push_back(Vertices);
        getline(file, line);
    }

    data.MaxId = id;


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
        data.Id_Lenght_Fractures.push_back({});
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

                    pair<unsigned int, double> Id_trace_lenght(N_Intersections,trace_lenght);

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

                        data.Id_Lenght_Fractures[i].push_back(Id_trace_lenght);
                        data.Id_Lenght_Fractures[j].push_back(Id_trace_lenght);

                        N_traces_fractures[i]++;
                        N_traces_fractures[j]++;

                        data.IdTraces[i].push_back(N_Intersections);
                        data.IdTraces[j].push_back(N_Intersections);
                        N_Intersections++;

                        if(cond1)
                        {
                            data.Tips.push_back({true,true});
                            data.BoolTraces[i].push_back(true);
                            data.BoolTraces[j].push_back(true);
                        }
                        else if(cond2)
                        {
                            if(cond4)
                            {
                                data.Tips.push_back({true,false});
                                data.BoolTraces[i].push_back(true);
                                data.BoolTraces[j].push_back(false);
                            }
                            else
                            {
                                data.Tips.push_back({false,true});
                                data.BoolTraces[i].push_back(false);
                                data.BoolTraces[j].push_back(true);
                            }
                        }
                        else if(cond3)
                        {
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
        sort(data.Id_Lenght_Fractures[k].begin(), data.Id_Lenght_Fractures[k].end(), [](const pair<unsigned int, double>& a, const pair<unsigned int,double>& b)
             {
                 return a.second > b.second;
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
            for(unsigned int l = 0; l != size(data.Id_Lenght_Fractures[id]); l++)
            {
                if(data.BoolTraces[id][l] == true)
                {
                    out << (data.Id_Lenght_Fractures[id][l]).first << " ; " << boolalpha << data.BoolTraces[id][l] << noboolalpha << " ; " << (data.Id_Lenght_Fractures[id][l]).second << endl;
                    Sorted_ids.push_back(data.Id_Lenght_Fractures[id][l].first);
                }
            }

            for(unsigned int l = 0; l != size(data.Id_Lenght_Fractures[id]); l++)
            {
                if(data.BoolTraces[id][l] == false)
                {
                    out << (data.Id_Lenght_Fractures[id][l]).first << " ; " << boolalpha << data.BoolTraces[id][l] << noboolalpha << " ; " << (data.Id_Lenght_Fractures[id][l]).second << endl;
                    Sorted_ids.push_back(data.Id_Lenght_Fractures[id][l].first);
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
    data.BoolTraces_Sorted = data.BoolTraces;

    for(unsigned int g = 0; g != size(data.BoolTraces_Sorted); g++)
    {
        partition(data.BoolTraces_Sorted[g].begin(), data.BoolTraces_Sorted[g].end(), [](bool b) { return b; });
    }

    unsigned int id_trace;
    vector<Vector3d> id_points_coordinates;
    vector<Vector2i> id_points_extremes;
    Vector3d vert0, vert1, trace0, trace1, diff1, diff2, diff3, det;
    Vector2d coeff;
    double coef0, coef1;
    Vector3d new_point;
    bool condition, condition0, condition1, condition2, condition3, condition4, condition5;
    bool control = false;
    bool control0 = false;
    bool control1 = false;
    vector<int> id_new_points;
    vector<int> ids_new_segment;
    int id_new_seg;
    vector<int> ids_new_segs;
    vector<vector<Vector3d>> Subpolygons, Subpolygons_checked;
    vector<vector<unsigned int>> SubpolygonsById;
    vector<int> ids_old_verts;
    Vector3d N,D,P;
    vector<Vector3d> Polygon1, Polygon2;
    vector<double> alfas;
    int poly1_start, poly1_end, poly2_start, poly2_end;
    vector<Vector3d> new_points;
    bool new_poly;
    vector<vector<vector<unsigned int>>> Export_points;
    vector<vector<Vector2i>> Export_extremes;
    vector<vector<Vector3d>> Export_coordinates_id;


    for(unsigned int i=0; i!= data.N_Fractures; i++) // Fratture
    {
        id_points_coordinates = data.Vertices[i];
        Subpolygons.push_back(data.Vertices[i]);

        if(data.TracesinFigures[i] != 0)
        {
            for(unsigned int t=0; t != size(id_points_coordinates); t++)
            {
                id_points_extremes.push_back({t,(t+1) % size(id_points_coordinates)});
            }

            for(unsigned int j=0; j!= data.TracesinFigures[i]; j++)  // Tracce
            {
                id_trace = data.Id_Traces_Sorted[i][j];
                trace0 = data.GeneratingPoints[id_trace][0];
                trace1 = data.GeneratingPoints[id_trace][1];


                for(unsigned int p = 0; p!=size(Subpolygons) ; p++ )  // Sottopoligoni
                {
                    unsigned int limit = size(Subpolygons[p]);

//-------------------------------------------------------------------------TRACCIA-PASSANTE---------

                    if(data.BoolTraces_Sorted[i][j] == true)
                    {
                        for(unsigned int v = 0; v!= limit; v++)  // vertici dei sottopoligoni
                        {
                            vert0 = Subpolygons[p][v];
                            vert1 = Subpolygons[p][ (v+1) % limit ];

                            diff1 = trace1-trace0;
                            diff2 = vert1-vert0;
                            diff3 = vert0-trace0;

                            det = diff1.cross(diff2);

                            if(det.norm() > data.Tol)
                            {
                                Matrix<double,3,2> AA;

                                AA << diff1[0], diff2[0],
                                      diff1[1], diff2[1],
                                      diff1[2], diff2[2];


                                Vector3d bb;

                                bb(0) = diff3[0]; bb(1) = diff3[1]; bb(2) = diff3[2];

                                coeff = AA.colPivHouseholderQr().solve(bb);

                                coef0 = coeff(0);
                                coef1 = coeff(1);

                                if(-coef1 > 0.0 && -coef1 < 1.0)
                                {
                                    new_point = (trace0 + coef0*(trace1-trace0));

                                    for(unsigned int u = 0; u != size(id_points_coordinates); u++)
                                    {
                                        condition = abs(id_points_coordinates[u][0] - new_point[0] ) < data.Tol &&
                                                    abs(id_points_coordinates[u][1] - new_point[1] ) < data.Tol &&
                                                    abs(id_points_coordinates[u][2] - new_point[2] ) < data.Tol;

                                        if(condition == true) // controllo se il punto cè già
                                        {
                                            ids_new_segment.push_back(u);
                                            id_new_seg = u;
                                            control = true;
                                        }
                                    }

                                    if(control == false) // devo aggiungere il punto
                                    {
                                        Subpolygons[p].push_back(new_point);
                                        id_points_coordinates.push_back(new_point);
                                        ids_new_segment.push_back(size(id_points_coordinates) - 1);
                                        id_new_seg = size(id_points_coordinates) - 1;
                                    }

                                    for(unsigned int u = 0; u != size(id_points_coordinates); u++)
                                    {
                                        condition0 = abs(id_points_coordinates[u][0] - vert0[0] ) < data.Tol &&
                                                    abs(id_points_coordinates[u][1] - vert0[1] ) < data.Tol &&
                                                    abs(id_points_coordinates[u][2] - vert0[2] ) < data.Tol;

                                        condition1 = abs(id_points_coordinates[u][0] - vert1[0] ) < data.Tol &&
                                                     abs(id_points_coordinates[u][1] - vert1[1] ) < data.Tol &&
                                                     abs(id_points_coordinates[u][2] - vert1[2] ) < data.Tol;

                                        if(condition0 == true)
                                        {
                                            ids_old_verts.push_back(u);
                                        }
                                        if(condition1 == true)
                                        {
                                            ids_old_verts.push_back(u);
                                        }

                                    }

                                    id_points_extremes.push_back({ids_old_verts[0], id_new_seg});
                                    id_points_extremes.push_back({id_new_seg, ids_old_verts[1]});

                                    for(unsigned int u = 0; u!= size(id_points_extremes); u++)
                                    {
                                        if( (id_points_extremes[u][0] == ids_old_verts[0] && id_points_extremes[u][1] == ids_old_verts[1]) ||
                                            (id_points_extremes[u][0] == ids_old_verts[1] && id_points_extremes[u][1] == ids_old_verts[0]) )
                                        {
                                            id_points_extremes.erase(id_points_extremes.begin() + u);
                                        }

                                    }

                                    ids_old_verts = {};

                                }
                            }
                        }

                        id_points_extremes.push_back({ids_new_segment[0],ids_new_segment[1]});
                        ids_new_segment = {};
                        new_poly = true;

                    }

//-----------------------------------------------------------------TRACCIA-NON-PASSANTE-------------

                    else if (data.BoolTraces_Sorted[i][j] == false)
                    {
                        for(unsigned int v = 0; v!= size(Subpolygons[p]); v++)  // vertici dei sottopoligoni
                        {
                            vert0 = Subpolygons[p][v];
                            vert1 = Subpolygons[p][ (v+1) % limit ];

                            diff1 = trace1-trace0;
                            diff2 = vert1-vert0;
                            diff3 = vert0-trace0;

                            det = diff1.cross(diff2);

                            if(det.norm() > data.Tol)
                            {
                                Matrix<double,3,2> AA;

                                AA << diff1[0], diff2[0],
                                    diff1[1], diff2[1],
                                    diff1[2], diff2[2];


                                Vector3d bb;

                                bb(0) = diff3[0]; bb(1) = diff3[1]; bb(2) = diff3[2];

                                coeff = AA.colPivHouseholderQr().solve(bb);

                                coef0 = coeff(0);
                                coef1 = coeff(1);

                                if((-coef1 > 0.0 && -coef1 < 1.0))
                                {
                                    new_points.push_back(trace0 + coef0*(trace1-trace0));

                                    for(unsigned int u = 0; u != size(id_points_coordinates); u++)
                                    {
                                        condition4 = abs(id_points_coordinates[u][0] - vert0[0] ) < data.Tol &&
                                                     abs(id_points_coordinates[u][1] - vert0[1] ) < data.Tol &&
                                                     abs(id_points_coordinates[u][2] - vert0[2] ) < data.Tol;

                                        condition5 = abs(id_points_coordinates[u][0] - vert1[0] ) < data.Tol &&
                                                     abs(id_points_coordinates[u][1] - vert1[1] ) < data.Tol &&
                                                     abs(id_points_coordinates[u][2] - vert1[2] ) < data.Tol;

                                        if(condition4 == true)
                                        {
                                            ids_old_verts.push_back(u);
                                        }
                                        if(condition5 == true)
                                        {
                                            ids_old_verts.push_back(u);
                                        }

                                    }
                                }
                            }
                        }

                        if(size(new_points)==2)
                        {
                            alfas.push_back(new_points[0][0] / diff1[0]);
                            alfas.push_back(new_points[1][0] / diff1[0]);
                        }

                        if(size(alfas) == 2)
                        {
                            if(!((alfas[0] > 1.0 && alfas[1] > 1.0 ) || (alfas[0] < 0.0 && alfas[1] < 0.0 )) )
                            {
                                for(unsigned int u = 0; u != size(id_points_coordinates); u++)
                                {
                                    condition2 = abs(id_points_coordinates[u][0] - new_points[0][0] ) < data.Tol &&
                                                 abs(id_points_coordinates[u][1] - new_points[0][1] ) < data.Tol &&
                                                 abs(id_points_coordinates[u][2] - new_points[0][2] ) < data.Tol;

                                    condition3 = abs(id_points_coordinates[u][0] - new_points[1][0] ) < data.Tol &&
                                                 abs(id_points_coordinates[u][1] - new_points[1][1] ) < data.Tol &&
                                                 abs(id_points_coordinates[u][2] - new_points[1][2] ) < data.Tol;

                                    if(condition2 == true) // controllo se il primo punto cè già
                                    {
                                        ids_new_segment.push_back(u);
                                        ids_new_segs.push_back(u);
                                        control0 = true;
                                    }
                                    if(condition3 == true) // controllo se il secondo punto cè già
                                    {
                                        ids_new_segment.push_back(u);
                                        ids_new_segs.push_back(u);
                                        control1 = true;
                                    }
                                }

                                if(control0 == false) // devo aggiungere il primo punto
                                {
                                    Subpolygons[p].push_back(new_points[0]);
                                    id_points_coordinates.push_back(new_points[0]);
                                    ids_new_segment.push_back(size(id_points_coordinates) - 1);
                                    ids_new_segs.push_back(size(id_points_coordinates) - 1);

                                }

                                if(control1 == false) // devo aggiungere il secondo punto
                                {
                                    Subpolygons[p].push_back(new_points[1]);
                                    id_points_coordinates.push_back(new_points[1]);
                                    ids_new_segment.push_back(size(id_points_coordinates) - 1);
                                    ids_new_segs.push_back(size(id_points_coordinates) - 1);
                                }

                                id_points_extremes.push_back({ids_old_verts[0], ids_new_segs[0]});
                                id_points_extremes.push_back({ids_new_segs[0], ids_old_verts[1]});
                                id_points_extremes.push_back({ids_old_verts[2], ids_new_segs[1]});
                                id_points_extremes.push_back({ids_new_segs[1], ids_old_verts[3]});

                                for(unsigned int u = 0; u!= size(id_points_extremes); u++)
                                {
                                    if( (id_points_extremes[u][0] == ids_old_verts[0] && id_points_extremes[u][1] == ids_old_verts[1]) ||
                                        (id_points_extremes[u][0] == ids_old_verts[1] && id_points_extremes[u][1] == ids_old_verts[0]) )
                                    {
                                        id_points_extremes.erase(id_points_extremes.begin() + u);
                                    }
                                    if( (id_points_extremes[u][0] == ids_old_verts[2] && id_points_extremes[u][1] == ids_old_verts[3]) ||
                                        (id_points_extremes[u][0] == ids_old_verts[3] && id_points_extremes[u][1] == ids_old_verts[2]) )
                                    {
                                        id_points_extremes.erase(id_points_extremes.begin() + u);
                                    }
                                }


                                id_points_extremes.push_back({ids_new_segment[0],ids_new_segment[1]});

                                new_poly = true;
                            }
                        }

                        new_points = {};
                        ids_old_verts = {};
                        alfas = {};
                        ids_new_segs = {};
                        ids_new_segment = {};

                    }

                    N = data.Normals[i];
                    D = diff1;

                    if(new_poly == true)
                    {
                        for(unsigned int x = 0; x!= size(Subpolygons[p]); x++)
                        {
                            if(abs(trace1.norm() - Subpolygons[p][x].norm() ) > data.Tol ||
                                abs(trace0.norm() - Subpolygons[p][x].norm() ) > data.Tol)
                            {
                                P = Subpolygons[p][x] - trace0;

                                if( ((D).cross(P)).dot(N) > data.Tol)
                                {
                                    Polygon1.push_back(Subpolygons[p][x]);
                                }
                                else if( ((D).cross(P)).dot(N) < -data.Tol )
                                {
                                    Polygon2.push_back(Subpolygons[p][x]);
                                }
                            }

                        }

                        for(unsigned int y = 0; y!= size(id_points_coordinates); y++)
                        {
                            if(id_points_coordinates[y] == Polygon1[0])
                            {
                                poly1_start = y;
                            }
                            else if(id_points_coordinates[y] == Polygon1[size(Polygon1) - 1])
                            {
                                poly1_end = y;
                            }
                            else if(id_points_coordinates[y] == Polygon2[0])
                            {
                                poly2_start = y;
                            }
                            else if(id_points_coordinates[y] == Polygon2[size(Polygon2) - 1])
                            {
                                poly2_end = y;
                            }
                        }

                        for(unsigned int z = 0; z!= size(id_points_extremes); z++)
                        {
                            if(id_points_extremes[z][0] == poly1_start && ( id_points_extremes[z][1] == ids_new_segment[0] || id_points_extremes[z][1] == ids_new_segment[1]) )
                            {
                                Polygon1.insert(Polygon1.begin(),id_points_coordinates[id_points_extremes[z][1]]);
                            }
                            else if(id_points_extremes[z][1] == poly1_start && ( id_points_extremes[z][0] == ids_new_segment[0] || id_points_extremes[z][0] == ids_new_segment[1]))
                            {
                                Polygon1.insert(Polygon1.begin(),id_points_coordinates[id_points_extremes[z][0]]);
                            }
                            else if(id_points_extremes[z][0] == poly2_start && ( id_points_extremes[z][1] == ids_new_segment[0] || id_points_extremes[z][1] == ids_new_segment[1]))
                            {
                                Polygon2.insert(Polygon2.begin(),id_points_coordinates[id_points_extremes[z][1]]);
                            }
                            else if(id_points_extremes[z][1] == poly2_start && ( id_points_extremes[z][0] == ids_new_segment[0] || id_points_extremes[z][0] == ids_new_segment[1]))
                            {
                                Polygon2.insert(Polygon2.begin(),id_points_coordinates[id_points_extremes[z][0]]);
                            }
                            else if(id_points_extremes[z][0] == poly1_end && ( id_points_extremes[z][1] == ids_new_segment[0] || id_points_extremes[z][1] == ids_new_segment[1]))
                            {
                                Polygon1.push_back(id_points_coordinates[id_points_extremes[z][1]]);
                            }
                            else if(id_points_extremes[z][1] == poly1_end && ( id_points_extremes[z][0] == ids_new_segment[0] || id_points_extremes[z][0] == ids_new_segment[1]))
                            {
                                Polygon1.push_back(id_points_coordinates[id_points_extremes[z][0]]);
                            }
                            else if(id_points_extremes[z][0] == poly2_end && ( id_points_extremes[z][1] == ids_new_segment[0] || id_points_extremes[z][1] == ids_new_segment[1]))
                            {
                                Polygon2.push_back(id_points_coordinates[id_points_extremes[z][1]]);
                            }
                            else if(id_points_extremes[z][1] == poly2_end && ( id_points_extremes[z][0] == ids_new_segment[0] || id_points_extremes[z][0] == ids_new_segment[1]))
                            {
                                Polygon2.push_back(id_points_coordinates[id_points_extremes[z][0]]);
                            }
                        }

                        Subpolygons_checked.push_back(Polygon1);
                        Subpolygons_checked.push_back(Polygon2);

                    }
                    else if(new_poly == false)
                    {
                        Subpolygons_checked.push_back(Subpolygons[p]);
                    }

                    new_poly = false;

                    Polygon1 = {};
                    Polygon2 = {};

                }

//-----------------------------FINE-CICLO-SUI-SOTTOPOLIGONI------------------------------

                Subpolygons = Subpolygons_checked;
                Subpolygons_checked = {};

            }

            for (const auto& subpolygon : Subpolygons)
            {
                vector<unsigned int> subpolygonById;
                for (const auto& point : subpolygon) {
                    for (unsigned int i = 0; i < id_points_coordinates.size(); ++i) {
                        if (point == id_points_coordinates[i]) {
                            subpolygonById.push_back(i);
                            break;
                        }
                    }
                }
                SubpolygonsById.push_back(subpolygonById);
            }

            Export_points.push_back(SubpolygonsById);
            Export_extremes.push_back(id_points_extremes);
            Export_coordinates_id.push_back(id_points_coordinates);
            Subpolygons = {};
            Subpolygons_checked = {};
            id_points_extremes = {};
            id_new_points = {};
            SubpolygonsById = {};

        }
        else
        {
            for (const auto& subpolygon : Subpolygons)
            {
                vector<unsigned int> subpolygonById;
                for (const auto& point : subpolygon) {
                    for (unsigned int i = 0; i < id_points_coordinates.size(); ++i) {
                        if (point == id_points_coordinates[i]) {
                            subpolygonById.push_back(i);
                            break;
                        }
                    }
                }
                SubpolygonsById.push_back(subpolygonById);
            }

            Export_points.push_back(SubpolygonsById);
            Export_extremes.push_back(id_points_extremes);
            Export_coordinates_id.push_back(id_points_coordinates);
        }

        Subpolygons = {};
        SubpolygonsById = {};


    }

    data.Export_Points = Export_points;
    data.Export_Extremes = Export_extremes;
    data.Export_Coordinates_Id = Export_coordinates_id;

    return true;

}

}

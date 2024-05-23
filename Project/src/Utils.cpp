#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <Eigen/Eigen>

using namespace Eigen;

namespace FracturesTraces
{

//--------------------------------CONTROLLO-GENERICO----------------------------------------------------------------------------

bool ImportData(DFN& data)
{
    if(!ImportAll("./FR3_data.txt", data))
    {
        return false;
    }

    if(!Testsfera(data))
    {
        return false;
    }

    if(!Testpianiparalleli(data))
    {
        return false;
    }

    if(!Testintersezione(data))
    {
        return false;
    }



    return true;


}

//--------------------------------INIZIALIZZO-LA-MAPPA-CHE-ASSOCIA-ID-A-UN-VETTORE-DI-VETTORI-CHE-DEFINISCE-LE-COORDINATE-------

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
    data.NumberFractures = NFractures;
    getline(file,line);

    for(unsigned int i = 0; i != NFractures; i++)
    {
        istringstream convertN(line);
        unsigned int id;
        unsigned int vertices;
        char tmp;
        double coordx;
        double coordy;
        double coordz;

        vector<double> coordinatesx;
        vector<double> coordinatesy;
        vector<double> coordinatesz;

        Vector3d vert;

        vector<Vector3d> Vertices;

        getline(file, line);
        convertN.str(line);
        convertN >> id >> tmp >> vertices;
        data.FracturesId.push_back(id);
        data.NumberVertices.push_back(vertices);

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

        data.Vertices.insert({id,Vertices});

        getline(file, line);

    }

    return true;

}

//--------------------------------TEST-DELLA-SFERA-E-DEL-BARICENTRO-------------------------------------------------------------

bool Testsfera(DFN &data)
{

    for(unsigned int id = 0; id != 3; id++)
    {
        double x_bari = 0.0;
        double y_bari = 0.0;
        double z_bari = 0.0;
        double tmpx = 0.0;
        double tmpy = 0.0;
        double tmpz = 0.0;
        vector<double> tmp_bari = {};

        for(unsigned int i = 0; i != 4; i++)
        {
            tmpx += data.Vertices[id][i][0];

        }

        x_bari = tmpx / 4;

        for(unsigned int i = 0; i != 4; i++ )
        {
            tmpy += data.Vertices[id][i][1];
        }

        y_bari = tmpy / 4;

        for(unsigned int i = 0; i != 4; i++ )
        {
            tmpz += data.Vertices[id][i][2];
        }

        z_bari = tmpz / 4;

        tmp_bari.push_back(x_bari);
        tmp_bari.push_back(y_bari);
        tmp_bari.push_back(z_bari);

        data.Baricentri.insert({id,tmp_bari});

    }

    double distance = 0.0;

    for(unsigned int id = 0; id != 3; id++)
    {
        double maxdistance = 0.0;

        for(unsigned int i = 0; i !=4; i++)
        {
            distance = pow(data.Baricentri[id][0] - data.Vertices[id][i][0],2) + pow(data.Baricentri[id][1] - data.Vertices[id][i][1],2) + pow(data.Baricentri[id][2] - data.Vertices[id][i][2],2);
            if(distance > maxdistance)
            {
                maxdistance = distance;
            }

        }

        data.raggi.insert({id,maxdistance});

    }

    return true;

}

//--------------------------------TEST-PIANI-PARALLELI--------------------------------------------------------------------------

bool Testpianiparalleli(DFN &data)
{

    Vector3d point0;
    Vector3d point1;
    Vector3d point2;
    Vector3d u;
    Vector3d v;

    for(unsigned int id = 0; id != 3; id++)
    {

        Vector3d normal;
        double d = 0.0;

        point0 = data.Vertices[id][0];
        point1 = data.Vertices[id][1];
        point2 = data.Vertices[id][2];

        u = point2-point0;
        v = point1-point0;

        normal = (u.cross(v));
        normal = normal.normalized();

        d = -normal[0]*point0[0] - normal[1]*point0[1] - normal[2]*point0[2];

        data.Normals.push_back(normal);
        data.Directors.insert({id,d});
    }

    return true;

}

//--------------------------------TEST-INTERSEZIONE-----------------------------------------------------------------------------

bool Testintersezione(DFN &data)
{

    for (unsigned int i = 0; i != data.NumberFractures;i++ )
        for ( unsigned int j = i +1; j != data.NumberFractures; j++)
        {

            // test sfera e piani paralleli



            Vector3d normale1;
            Vector3d normale2;
            Vector3d director;

            normale1 = data.Normals[i];
            normale2 = data.Normals[j];

            director = normale1.cross(normale2);

            Matrix3d A;
            A << normale1[0],normale1[1],normale1[2],
                normale2[0],normale2[1],normale2[2],
                director[0],director[1],director[2];

            Vector3d b;

            b << -data.Directors[i] , -data.Directors[j],0;

            Vector3d P0 = A.colPivHouseholderQr().solve(b);  //primo punto sulla retta di intersezione

            Vector3d P1 = P0+director;   //secondo punto sulla retta di intersezione


            const unsigned int numVertices = data.Vertices[i].size();
            for(unsigned int w = 0; w < numVertices; w++)
            {
                const Vector3d punto0 = data.Vertices[i][w];
                const Vector3d punto1 = data.Vertices[i][(w + 1) % numVertices];

                Vector3d diff1 = P1-P0;
                Vector3d diff2 = punto1-punto0;
                Vector3d diff3 = punto0-P0;

                Vector3d sos = diff1.cross(diff2);


                if(sos.norm() > 1e-10)
                {
                    Matrix<double,3,2> AA;

                    AA << diff1[0],diff2[0],
                        diff1[1],diff2[1],
                        diff1[2],diff2[2];


                    Vector3d bb;
                    bb(0) = diff3[0];
                    bb(1) = diff3[1];
                    bb(2) = diff3[2];


                    Vector2d intersection = AA.colPivHouseholderQr().solve(bb);

                    double alfa = intersection(0);
                    double beta = intersection(1);

                    cout << alfa << " " << beta << endl;

                    cout << "punto:" << endl << P0+alfa*(P1-P0) << endl << "punto:" << endl << punto0-beta*(punto1-punto0) << endl;


                }

            }

        }

    return true;

}

}

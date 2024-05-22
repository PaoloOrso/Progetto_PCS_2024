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

        vector<double> vert;

        vector<vector<double>> Vertices;

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
            vert.push_back(coordinatesx[i]);
            vert.push_back(coordinatesy[i]);
            vert.push_back(coordinatesz[i]);

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

    vector<double> normal;
    vector<double> point0;
    vector<double> point1;
    vector<double> point2;

    for(unsigned int id = 0; id != 3; id++)
    {
        double normalx = 0.0;
        double normaly = 0.0;
        double normalz = 0.0;
        double d = 0.0;

        point0 = data.Vertices[id][0];
        point1 = data.Vertices[id][1];
        point2 = data.Vertices[id][2];

        normalx = (point1[1]-point0[1])*(point2[2]-point0[2])-(point1[2]-point0[2])*(point2[1]-point0[1]);
        normaly = (point1[2]-point0[2])*(point2[0]-point0[0])-(point1[0]-point0[0])*(point2[2]-point0[2]);
        normalz = (point1[0]-point0[0])*(point2[1]-point0[1])-(point1[1]-point0[1])*(point2[0]-point0[0]);

        d = -normalx*point0[0] - normaly*point0[1] - normalz*point0[2];

        normal.push_back(normalx);
        normal.push_back(normaly);
        normal.push_back(normalz);

        data.Normals.insert({id,normal});
        data.Directors.insert({id,d});

        normal = {};

    }

    return true;


}

//--------------------------------TEST-INTERSEZIONE-----------------------------------------------------------------------------

bool Testintersezione(DFN &data)
{

     // for (unsigned int i = 0; i != data.NumberFractures;i++ )
     //     for ( unsigned int j = i +1; j != data.NumberFractures; j++)
     //     {

                unsigned int i = 0;
                unsigned int j = 1;

                vector<double> vettore1;
                vector<double> vettore2;
                vector<double> director;

                vettore1 = data.Normals[i];
                vettore2 = data.Normals[j];

                director.push_back((vettore1[1]*vettore2[2])-(vettore1[2]*vettore2[1]));
                director.push_back((vettore1[2]*vettore2[0])-(vettore1[0]*vettore2[2]));
                director.push_back((vettore1[0]*vettore2[1])-(vettore1[1]*vettore2[0]));

                Matrix3d A;
                A << vettore1[0],vettore1[1],vettore1[2],
                    vettore2[0],vettore2[1],vettore2[2],
                    director[0],director[1],director[2];

                Vector3d b;

                b << -data.Directors[i] , -data.Directors[j],0;

                Vector3d solution = A.colPivHouseholderQr().solve(b);

                double x0 = solution(0);
                double y0 = solution(1);
                double z0 = solution(2);

                // cout << "Punto sulla retta: " << "X:" << x0 << " Y:" << y0 << " Z:" << z0 << endl;

                vector<vector<double>> vertici1;
                vector<vector<double>> vertici2;

                for(unsigned int k = 0; k != data.NumberVertices[i]; k++)
                {
                    vertici1.push_back(data.Vertices[i][k]);
                }
                vertici1.push_back(data.Vertices[i][0]);

                for(unsigned int l = 0; l != data.NumberVertices[j]; l++)
                {
                    vertici2.push_back(data.Vertices[j][l]);
                }
                vertici2.push_back(data.Vertices[j][0]);




                for(unsigned int w = 0; w != data.NumberVertices[i]; w++)
                {




                Matrix2d AA;
                AA(0, 0) = director[i];            AA(0, 1) = -(vertici1[w+1][0] - vertici1[w][0]);
                AA(1, 0) = director[j];            AA(1, 1) = -(vertici1[w+1][1] - vertici1[w][1]);

                Vector2d bb;
                bb(0) = vertici1[w][0] - x0;
                bb(1) = vertici1[w][1] - y0;


                Vector2d x = AA.colPivHouseholderQr().solve(bb);

                double t = x(0);
                double u = x(1);


                if(u>=0 && u <=1)
                {
                    cout << "Intersezione con lato" << endl;
                }

                }



       //  }











    return true;

}

}

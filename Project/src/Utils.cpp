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

using namespace Eigen;

namespace FracturesTraces
{

//--------------------------------CONTROLLO-GENERICO----------------------------------------------------------------------------

bool ImportData(DFN& data)
{
    // 3 10 50 82 200 362
    if(!ImportAll("./FR10_data.txt", data))
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
        data.MaxId = id;
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

    for(unsigned int id = 0; id != (data.MaxId + 1); id++)
    {
        Vector3d bari;
        unsigned int num = data.NumberVertices[id];

        double tmpx = 0.0;
        double tmpy = 0.0;
        double tmpz = 0.0;

        for(unsigned int i = 0; i != num; i++)
        {
            Vector3d point = data.Vertices[id][i];
            tmpx += point[0];
            tmpy += point[1];
            tmpz += point[2];

        }

        bari[0] = tmpx / num;
        bari[1] = tmpy / num;
        bari[2] = tmpz / num;

        data.Baricentri.push_back(bari);

        double distance = 0.0;
        double maxdistance = 0.0;

        for(unsigned int i = 0; i != num; i++)
        {
            Vector3d point = data.Vertices[id][i];
            distance = sqrt(pow(bari[0] - point[0],2) + pow(bari[1] - point[1],2) + pow(bari[2] - point[2],2));
            if(distance > maxdistance)
            {
                maxdistance = distance;
            }

        }

        data.raggi.push_back(maxdistance);

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
    
    for(unsigned int id = 0; id != (data.MaxId + 1); id++)
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
        data.Ds.push_back(d);
    }

    return true;

}

//--------------------------------TEST-INTERSEZIONE-----------------------------------------------------------------------------

bool Testintersezione(DFN &data)
{
    unsigned int NUMEROINTERSEZIONI = 0;
    vector<unsigned int> Frig_frac(data.NumberFractures);

    for (unsigned int i = 0; i != data.NumberFractures;i++ )
        for ( unsigned int j = i +1; j != data.NumberFractures; j++)
        {

            Vector3d bari1 = data.Baricentri[i];
            Vector3d bari2 = data.Baricentri[j];

            Vector3d normale1;
            Vector3d normale2;

            normale1 = data.Normals[i];
            normale2 = data.Normals[j];

            double distanza_bari = sqrt(pow(bari1[0]-bari2[0],2) + pow(bari1[1]-bari2[1],2) + pow(bari1[2]-bari2[2],2));
            double somma_raggi = data.raggi[i] + data.raggi[j];

            if(( distanza_bari - somma_raggi < 1e-10 ) && ( (normale1.cross(normale2)).norm() > 1e-10 ))  // test baricentro e piani paralleli
            {

            vector<double> test;
            Vector3d director;

            director = normale1.cross(normale2);     // vettore direttore della retta di intersezione tra i due piani

            Matrix3d A;
            A << normale1[0],normale1[1],normale1[2],
                normale2[0],normale2[1],normale2[2],
                director[0],director[1],director[2];

            Vector3d b;

            b << -data.Ds[i] , -data.Ds[j],0;

            Vector3d P0 = A.colPivHouseholderQr().solve(b);  //primo punto sulla retta di intersezione

            Vector3d P1 = P0+director;   //secondo punto sulla retta di intersezione

            unsigned int numVertices = data.NumberVertices[i];
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

                    // Vector3d punto = punto0-beta*(punto1-punto0);

                    if(-beta > 0.0 && -beta < 1.0)
                    {
                        test.push_back(alfa);
                        //cout << "punto poligono " << i << ": " << punto[0] << " , " << punto[1] << " , " << punto[2] << endl;
                    }
                }
            }

            numVertices = data.NumberVertices[j];
            for(unsigned int w = 0; w < numVertices; w++)
            {
                const Vector3d punto0 = data.Vertices[j][w];
                const Vector3d punto1 = data.Vertices[j][(w + 1) % numVertices];

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

                    // Vector3d Punto = punto0-beta*(punto1-punto0);

                    if(-beta > 0.0 && -beta < 1.0)
                    {
                        test.push_back(alfa);
                        //cout << "punto poligono " << j << ": " << Punto[0] << " , " << Punto[1] << " , " << Punto[2] << endl;
                    }
                }
            }

            cout << endl;

            if( size(test) == 4)
            {
                double lunghezza = 0.0;
                double delta_alfa = 0.0;

                Vector2i gen_frac;
                vector<Vector3d> gen_points;
                Vector3d PT0;
                Vector3d PT1;

                delta_alfa = abs(test[1]-test[2]);

                lunghezza = delta_alfa * director.norm();

                gen_frac = {i,j};

                PT0 = P0 + test[1]*(P1-P0);
                PT1 = P0 + test[2]*(P1-P0);

                gen_points = {PT0,PT1};

                if(abs(test[0] - test[2]) < 1e-10 && abs(test[1] - test[3]) < 1e-10)
                {
                    NUMEROINTERSEZIONI++;
                    data.GeneratingFractures.push_back(gen_frac);
                    data.GeneratingPoints.push_back(gen_points);
                    data.LenghtTraces.push_back(lunghezza);
                    data.Tips.push_back({1,1});
                    Frig_frac[i]++;
                    Frig_frac[j]++;
                    cout << "Due fratture passanti tra poligoni " << i << " e " << j << " di lunghezza " << lunghezza << endl << endl;
                }
                else if(((max(test[0],test[1]) >= max(test[2],test[3])) && (min(test[0],test[1]) <= min(test[2],test[3]))) || ((max(test[2],test[3]) >= max(test[0],test[1])) && (min(test[2],test[3]) <= min(test[0],test[1]))))
                {
                    NUMEROINTERSEZIONI++;
                    data.GeneratingFractures.push_back(gen_frac);
                    data.GeneratingPoints.push_back(gen_points);
                    data.LenghtTraces.push_back(lunghezza);
                    data.Tips.push_back({0,1});  // sistema l'ordine di 1 e 0 nei due casi
                    Frig_frac[i]++;
                    Frig_frac[j]++;
                    cout << "Una frattura passante e una non passante tra poligoni " << i << " e " << j << " di lunghezza " << lunghezza << endl << endl;
                }
                else if(((max(test[0],test[1]) > min(test[2],test[3])) && (min(test[0],test[1]) < min(test[2],test[3]))) || ((max(test[2],test[3]) > min(test[0],test[1])) && (min(test[2],test[3]) < min(test[0],test[1]))))
                {
                    NUMEROINTERSEZIONI++;
                    data.GeneratingFractures.push_back(gen_frac);
                    data.GeneratingPoints.push_back(gen_points);
                    data.LenghtTraces.push_back(lunghezza);
                    data.Tips.push_back({0,0});
                    Frig_frac[i]++;
                    Frig_frac[j]++;
                    cout << "Due fratture non passanti tra poligoni " << i << " e " << j << " di lunghezza " << lunghezza << endl << endl;
                }
                // else if((max(test[0],test[1]) < min(test[2],test[3])) || ((max(test[2],test[3]) < min(test[0],test[1]))))
                // {
                //     cout << "No intersezioni" << endl << endl;
                // }
            }
            // else
            // {
            //     cout << "No intersezioni" << endl << endl;
            // }
        }
        }

    data.NumberTraces = NUMEROINTERSEZIONI;
    data.TracesinFigures = Frig_frac;
    cout << "Numero intersezioni: " << NUMEROINTERSEZIONI << endl;

    return true;

}

//---------------------------------STAMPA-RISULTATI-----------------------------------------------------------------------------


}

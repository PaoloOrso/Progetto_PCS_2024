#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>


namespace FracturesTraces
{

//---------------------------------------------------------------CONTROLLO-GENERICO--------------------------------------------------------------------------------

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

    return true;


}

//-------------------------------------------------INIZIALIZZO-LA-MAPPA-CHE-ASSOCIA-ID-A-UN-VETTORE-DI-VETTORI-CHE-DEFINISCE-LE-COORDINATE-------------------------

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

//-------------------------------------------------------TEST-DELLA-SFERA-E-DEL-BARICENTRO-------------------------------------------------------------------------

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
            cout << distance;
            if(distance > maxdistance)
            {
                maxdistance = distance;
            }

        }

        data.raggi.insert({id,maxdistance});

    }

    return true;

}

//-------------------------------------------------------------TEST-PIANI-PARALLELI-------------------------------------------------------------------------

bool Testpianiparalleli(DFN &data);


}

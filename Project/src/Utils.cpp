#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>


namespace FracturesTraces
{

//---------------------------------------------------------------CONTROLLO-GENERICO--------------------------------------------------------------------------------

bool ImportData(const string& filepath, DFN& data)
{
    if(!ImportAll(filepath + "/FR3_data.txt", data))
    {
        return false;
    }



    if(!Baricentri(data))
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
        double coord;

        vector<double> coordinates;
        vector<vector<double>> Vertices;

        getline(file, line);
        convertN.str(line);
        convertN >> id >> tmp >> vertices;
        data.FracturesId.push_back(id);
        data.NumberVertices.push_back(vertices);

        getline(file, line);

        for(unsigned int i = 0; i!=3; i++)
        {
            getline(file, line);
            replace(line.begin(),line.end(), ';' ,' ');
            istringstream convert(line);

            for(unsigned int j = 0; j != vertices; j++)
            {
                convert >> coord;
                //cout << coord << ' ';                             // STAMPO LE COORDINATE COME CONTROLLO
                coordinates.push_back(coord);

            }
            Vertices.push_back(coordinates);
            coordinates = {};
            //cout << endl;                                         // STAMPO LE COORDINATE COME CONTROLLO
        }

        data.Vertices.insert({id,Vertices});

        getline(file, line);

    }

    return true;

}

//-------------------------------------------------------TEST-DELLA-SFERA-E-DEL-BARICENTRO------------------------------------------------------------------------

bool Baricentri(DFN &data)

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
            tmpx += data.Vertices[id][0][i];
            //cout << "x " << tmpx << endl;

        }

        x_bari = tmpx / 4;
        //cout << fixed << setprecision(6) << "barix " << x_bari << endl;


        for(unsigned int i = 0; i != 4; i++ )
        {
            tmpy += data.Vertices[id][1][i];
            //cout << "y " << tmpy << endl;
        }

        y_bari = tmpy / 4;
        //cout << fixed << setprecision(6) << "bariy " <<  y_bari << endl;


        for(unsigned int i = 0; i != 4; i++ )
        {
            tmpz += data.Vertices[id][2][i];
            //cout << "z " << tmpz << endl;
        }

        z_bari = tmpz / 4;
        //cout << fixed << setprecision(6) << "bariz " << z_bari<< endl;

        tmp_bari.push_back(x_bari);
        tmp_bari.push_back(y_bari);
        tmp_bari.push_back(z_bari);

        data.Baricentri.insert({id,tmp_bari});

    }

    return true;

}

}

// vector<double> calcolaBaricentro(const vector<vector<double>>& vertici)
// {
//     vector<double> baricentro(3, 0.0);
//     int numVertici = 4;


//     for (const auto& vertice : vertici)
//     {
//         for (int i = 0; i < 3; ++i)
//         {
//             baricentro[i] += vertice[i];
//         }
//     }


//     for (int i = 0; i < 3; ++i)
//     {
//         baricentro[i] /= numVertici;
//     }

//     return baricentro;
// }

// int main() {
//     Data data;


//     // Calcolo e stampa dei baricentri
//     for (const auto& pair : data.Vertices) {
//         std::vector<double> baricentro = calcolaBaricentro(pair.second);
//         std::cout << "Baricentro del gruppo " << pair.first << ": ";
//         for (const auto& coord : baricentro) {
//             std::cout << coord << " ";
//         }
//         std::cout << std::endl;
//     }

//     return 0;
// }


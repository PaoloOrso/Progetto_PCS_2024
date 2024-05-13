#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>


namespace FracturesTraces
{


bool ImportData(const string& filepath, DFN& data)
{
    if(!ImportAll(filepath + "/FR3_data.txt", data))
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool ImportAll(const string &filename, DFN &data)
{
    ifstream file;
    file.open(filename);

    if(file.fail())
    {
        cerr << "errore" << endl;
        return false;
    }

    // string line;

    // unsigned int NFractures = 0;

    // istringstream convertN(line);

    // getline(file, line);
    // getline(file, line);
    // convertN.str(line);
    // convertN >> NFractures;
    // getline(file,line);

    // for(unsigned int i = 0; i != NFractures; i++)
    // {
    //     istringstream convertN(line);
    //     unsigned int id;
    //     unsigned int vertices;
    //     char tmp;
    //     double coord;

    //     vector<double> coordinates;
    //     vector<vector<double>> Vertices;

    //     getline(file, line);
    //     convertN.str(line);
    //     convertN >> id >> tmp >> vertices;
    //     data.FracturesId.push_back(id);
    //     data.NumberVertices.push_back(vertices);


    //     getline(file, line);

    //     for(unsigned int i = 0; i!=3; i++)
    //     {
    //         getline(file, line);
    //         replace(line.begin(),line.end(), ';' ,' ');


    //         for(unsigned int j = 0; j != vertices; j++)  // il problema è qui
    //         {
    //             istringstream convert(line);
    //             convert >> coord;                         //
    //             cout << coord << endl;
    //             coordinates.push_back(coord);            //
    //         }                                            //

    //         Vertices.push_back(coordinates);
    //     }



    //     getline(file, line);

    vector<vector<double>> vettore1 = {{0.0,0.0,0.0},{1.0,0.0,0.0},{1.0,1.0,0.0},{0.0,1.0,0.0}};
    data.Vertices.insert({0,vettore1});
    vector<vector<double>> vettore2 = {{0.8,0.0,-0.1},{0.8,0.0,0.3},{0.8,1.0,0.3},{0.8,1.0,-0.1}};
    data.Vertices.insert({1,vettore2});
    vector<vector<double>> vettore3 = {{-0.238,0.5,-0.344},{0.301,0.5,-0.344},{0.316,0.5,0.453},{-0.238,0.5,0.453}};
    data.Vertices.insert({2,vettore3});

    return true;

}

bool TestSfera(DFN &data)

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
        tmpx += data.Vertices[0][i][0];
        cout << tmpx << endl;

    }

    x_bari = tmpx / 4;
    cout << x_bari;


    for(unsigned int i = 0; i != 4; i++ )
    {
        tmpy += data.Vertices[0][i][1];
        cout << tmpy << endl;
    }

    y_bari = tmpy / 4;
    cout << fixed << setprecision(6) << y_bari;


    for(unsigned int i = 0; i != 4; i++ )
    {
        tmpz += data.Vertices[0][i][2];
        cout << tmpz << endl;
    }

    z_bari = tmpz / 4;
    cout << fixed << setprecision(6) << z_bari;


    tmp_bari.push_back(x_bari);
    tmp_bari.push_back(y_bari);
    tmp_bari.push_back(z_bari);

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


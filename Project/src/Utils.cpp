#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


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
        istringstream convert(line);

        getline(file, line);

        for(unsigned int i = 0; i!=3; i++)
        {
            getline(file, line);
            replace(line.begin(),line.end(), ';' ,' ');

            for(unsigned int j = 0; j != vertices; j++)  // il problema è qui
            {
                convert >> coord;                        //
                coordinates.push_back(coord);            //
            }                                            //

            Vertices.push_back(coordinates);
        }

        getline(file, line);

    }

    return true;

}



}

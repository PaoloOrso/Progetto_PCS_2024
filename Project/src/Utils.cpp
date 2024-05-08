#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

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
    bool end = false;

    istringstream convertN(line);


    getline(file, line);
    getline(file, line);
    cout << line << endl;
    convertN.str(line);
    convertN >> NFractures;
    getline(file,line);

    for(unsigned int i = 0; i != NFractures; i++)
    {
        unsigned int id;
        unsigned int vertices;
        char tmp;
        istringstream convertN(line);

        getline(file, line);
        cout << line << endl;

        convertN >> id >> tmp >> vertices;
        data.FracturesId.push_back(id);
        data.NumberVertices.push_back(vertices);
        getline(file, line);
        getline(file, line);

        convertN >> data.Vertices(0,0) >> data.Vertices(0,1) >> data.Vertices(0,2) >> data.Vertices(0,3);
        getline(file, line);

        convertN >> data.Vertices(1,0) >> data.Vertices(1,1) >> data.Vertices(1,2) >> data.Vertices(1,3);
        getline(file, line);

        convertN >> data.Vertices(2,0) >> data.Vertices(2,1) >> data.Vertices(2,2) >> data.Vertices(2,3);
        getline(file, line);

    }

    end = true;


    return true;

}



}













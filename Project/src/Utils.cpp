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
    char coma;
    unsigned int id = 0;
    unsigned int vertices = 0;
    unsigned int NFractures = 0;
    bool end = false;

    istringstream convertN;


    getline(file, line);
    getline(file, line);
    cout << line << endl;
    convertN.str(line);
    convertN >> NFractures;
    getline(file,line);

    for(unsigned int i = 0; i != NFractures; i++)
    {
        getline(file, line);
        cout << line << endl;
        convertN.str(line);
        cout << line << endl;
        convertN >> id >> coma >> vertices;
        data.FracturesId.push_back(id);
        data.NumberVertices.push_back(vertices);
        getline(file, line);
        getline(file, line);
        getline(file, line);
        getline(file, line);
        getline(file, line);

    }

    end = true;


    return true;

}



}













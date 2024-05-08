#include <iostream>
#include "Data.hpp"
#include "Utils.hpp"
#include <sstream>

using namespace std;
using namespace Eigen;
using namespace FracturesTraces;

int main()
{
    DFN data;

    string filepath = "DFN";
    if(!ImportData(filepath, data))
    {
        return 1;
    }
}

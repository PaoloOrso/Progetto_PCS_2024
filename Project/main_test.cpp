#include <iostream>
#include <gtest/gtest.h>
#include "TestLibrary.hpp"

using namespace std;
using namespace Eigen;
using namespace FracturesTraces;

int main(int argc, char ** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

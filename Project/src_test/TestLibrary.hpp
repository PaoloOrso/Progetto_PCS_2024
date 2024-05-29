#ifndef __TESTPOLYGONS_H
#define __TESTPOLYGONS_H

#include <gtest/gtest.h>
#include "Data.hpp"
#include "Utils.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include "UCDUtilities.hpp"

using namespace Eigen;
using namespace std;

namespace FracturesTraces {
//------------------------------------TEST-FINALE--------------
// TEST(TEST_FINALE, TestFinalTest){

//     Matrix3d vertices = Matrix3d::Zero();

//     vertices << 0.0, 1.0, 0.0,
//         0.0, 0.0, 1.0,
//         0.0, 0.0, 0.0;

//     Triangle t(vertices);

//     double area = t.computeArea();
//     EXPECT_EQ(area, 0.5);
// }

//--------------------------------TEST-IMPORT-------------------

TEST(IMPORT_TEST, TestImportAll)
{

    DFN test1;
    ImportAll("./Test_data.txt", test1);


    vector<unsigned int> num_vertices_true = {3,4,5,4};
    vector<vector<Vector3d>> vertices_true = {  {{0,0,0},{1,0,0},{0,1,0}},
                                              {{2,0,0},{2,2,0},{0,2,2},{0,0,2}},
                                              {{1,-1,1},{3,0,1},{2,3,1},{0,3,1},{-1,0,1}},
                                              {{0.2,0.2,0.1},{0.2,0.4,0.1},{0.2,0.2,-0.1},{0.2,0.4,-0.1}}
                                             };

    EXPECT_EQ(test1.NumberFractures, 4);
    EXPECT_EQ(test1.MaxId, 3);
    EXPECT_EQ(test1.NumberVertices, num_vertices_true);
    EXPECT_EQ(test1.Vertices, vertices_true);

}

//------------------------------TEST-SFERA-------------------------

TEST(SFERA_TEST, TestTestsfera)
{

    DFN test2;
    test2.NumberFractures = 4;
    test2.MaxId = 3;
    test2.NumberVertices = {4,4,3,5};
    test2.Vertices = {   {{0,0,0},{1,0,0},{1,1,0},{0,1,0}},
                         {{0,1,2},{0,1,3},{0,2,4},{0,2,1}},
                         {{3.5,2.5,0},{2.5,3.5,0},{0,0,9}},
                         {{0.2,0,0},{0.6,0,0},{0.8,0,3},{0.4,0,5},{0,0,2}}
                      };

    Testsfera(test2);

    vector<Vector3d> bari_true = {{0.5,0.5,0},{0,1.5,2.5},{2,2,3},{0.4,0,2}};
    vector<double> raggi_true = {0.707106,0.632455,7.011419,1.0770329615};


    EXPECT_EQ(test2.Baricentri,bari_true);
    EXPECT_EQ(test2.raggi, raggi_true );
}


// //------------------------------TEST-PIANI-PARALLELI----------------

// TEST(PIANI_TEST, TestTestpianiparalleli){

//             DFN test = Testsfera();

//     EXPECTED_EQ()
//     EXPECTED_EQ()

// }

// //-------------------------------TESTO-VARIE-INTERSEZIONI-------------------------

// TEST(INTERSEZIONI_TEST, TestTestintersezione{

//     DFN test = Testsfera();

//     EXPECTED_EQ()
//     EXPECTED_EQ()

// }

// //-------------------------------TEST-STAMPA---------------------------------

// TEST(SPAMPA_TEST,TestStampa){



// }

// }
}


#endif

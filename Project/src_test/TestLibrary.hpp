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
    vector<double> raggi_true = {0.707107,1.581138,6.633249,3};


    for(unsigned int i = 0; i != 4; i++)
    {
    EXPECT_TRUE(test2.Baricentri[i].norm() - bari_true[i].norm() < 1e-6);
    EXPECT_TRUE(abs(test2.raggi[i]- raggi_true[i]) < 1e-6);
    }

}


//------------------------------TEST-PIANI-PARALLELI----------------

TEST(PIANI_TEST, TestTestpianiparalleli)
{

    DFN test3;
    test3.NumberFractures = 4;
    test3.MaxId = 3;
    test3.NumberVertices = {4,4,3,3};
    test3.Vertices = {   {{0,0,0},{1,0,0},{1,1,0},{0,1,0}},
                        {{0,0,2},{1,0,2},{1,1,2},{0,1,2}},
                        {{5,0,0},{5,2,0},{5,1,2}},
                        {{-1,1,0},{-1,3,0},{-3,2,4}},
                     };

    Testpianiparalleli(test3);

    vector<double> d_true = {0,2,5,-0.894427};


    for(unsigned int i = 0; i!= 4; i++)
    {
    EXPECT_TRUE(abs(test3.Normals[i].norm()) -1 < 1e-6);
    EXPECT_TRUE(abs(test3.Ds[i] - d_true[i]) < 1e-6 );
    }
 }

//-------------------------------TESTO-VARIE-INTERSEZIONI-------------------------

 TEST(INTERSEZIONI_TEST, TestTestintersezione)
{

    DFN test4;
    ImportAll("./Test_data2.txt", test4);
    Testsfera(test4);
    Testpianiparalleli(test4);
    Testintersezione(test4);

    vector<Vector2i> gen_frac_true = {{0,1},{0,2},{3,4}};
    vector<double> Len_traces_true = {0.5,0.4,0.4};
    vector<vector<bool>> tips_true = {{false,false},{false,true},{true,true}};
    vector<unsigned int> traces_figures_true = {2,1,1,1,1};
    vector<vector<unsigned int>> id_traces_true = {{0,1},{0},{1},{2},{2}};


    EXPECT_EQ(test4.GeneratingFractures, gen_frac_true);
    EXPECT_EQ(test4.Tips, tips_true);
    EXPECT_EQ(test4.TracesinFigures, traces_figures_true);
    EXPECT_EQ(test4.IdTraces, id_traces_true);

    for(unsigned int i = 0; i != test4.NumberTraces;i++)
    {
        EXPECT_TRUE(abs(test4.LenghtTraces[i] - Len_traces_true[i]) < 1e-6);
    }
}

//--TEST_STAMPA?----


//------------------------------------TEST-FINALE--------------
TEST(TEST_FINALE, TestFinalTest){

    DFN test6;
    ImportAll("./Test_data2.txt", test6);
    Testsfera(test6);
    Testpianiparalleli(test6);
    Testintersezione(test6);
    Stampa(test6);

    bool final_test_true = false;

    EXPECT_EQ(ImportAll("./Test_data2.txt", test6) && !(Testsfera(test6)) && Testpianiparalleli(test6)
                  && Testintersezione(test6) && Stampa(test6), final_test_true);


}

}


#endif

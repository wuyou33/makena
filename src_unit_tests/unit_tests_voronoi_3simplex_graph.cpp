#include "gtest/gtest.h"
#include <vector>
#include "primitives.hpp"
#include "manifold.hpp"
#include "quaternion.hpp"
#include "voronoi_3simplex_graph.hpp"

#include <iostream>
#include <fstream>
#include <sstream>


namespace Makena {


class Voronoi3SimplexGraphTest : public ::testing::Test {

  protected:  


    Voronoi3SimplexGraphTest(){;};
    virtual ~Voronoi3SimplexGraphTest(){;}
    virtual void SetUp()    {;};
    virtual void TearDown() {;};

};


static Mat3x3 randomRotMat()
{
    Vec3 v(0.0, 0.0, 0.0);
    while (v == Vec3(0.0, 0.0, 0.0)) {
        v = Vec3( (rand()%100)/100.0 - 0.5, 
                  (rand()%100)/100.0 - 0.5, 
                  (rand()%100)/100.0 - 0.5 );
    }
    
    Quaternion q(v, (rand()%100)/50.0 * M_PI);
    q.normalize();
    return q.rotationMatrix();
}

/*
static void printPredicate(const enum predicate pred) {
    switch (pred) {

      case VORONOI_INSIDE_TETRAHEDRON:
        cerr << "VORONOI_INSIDE_TETRAHEDRON\n";
        break;

      case VORONOI_OVER_TRIANGLE_1_3_2:
        cerr << "VORONOI_OVER_TRIANGLE_1_3_2\n";
        break;

      case VORONOI_OVER_TRIANGLE_1_2_4:
        cerr << "VORONOI_OVER_TRIANGLE_1_2_4\n";
        break;

      case VORONOI_OVER_TRIANGLE_1_4_3:
        cerr << "VORONOI_OVER_TRIANGLE_1_4_3\n";
        break;

      case VORONOI_OVER_TRIANGLE_2_3_4:
        cerr << "VORONOI_OVER_TRIANGLE_2_3_4\n";
        break;

      case VORONOI_OVER_EDGE_1_2:
        cerr << "VORONOI_OVER_EDGE_1_2\n";
        break;

      case VORONOI_OVER_EDGE_2_3:
        cerr << "VORONOI_OVER_EDGE_2_3\n";
        break;

      case VORONOI_OVER_EDGE_3_1:
        cerr << "VORONOI_OVER_EDGE_3_1\n";
        break;

      case VORONOI_OVER_EDGE_1_4:
        cerr << "VORONOI_OVER_EDGE_1_4\n";
        break;

      case VORONOI_OVER_EDGE_2_4:
        cerr << "VORONOI_OVER_EDGE_2_4\n";
        break;

      case VORONOI_OVER_EDGE_3_4:
        cerr << "VORONOI_OVER_EDGE_3_4\n";
        break;

      case VORONOI_OVER_VERTEX_1:
        cerr << "VORONOI_OVER_VERTEX_1\n";
        break;

      case VORONOI_OVER_VERTEX_2:
        cerr << "VORONOI_OVER_VERTEX_2\n";
        break;

      case VORONOI_OVER_VERTEX_3:
        cerr << "VORONOI_OVER_VERTEX_3\n";
        break;

      case VORONOI_OVER_VERTEX_4:
        cerr << "VORONOI_OVER_VERTEX_4\n";
        break;

      case VORONOI_INSIDE_TRIANGLE:
        cerr << "VORONOI_INSIDE_TRIANGLE\n";
        break;

      case BETWEEN_1_AND_2:
        cerr << "BETWEEN_1_AND_2\n";
        break;

      case ON_POINT1:
        cerr << "ON_POINT1\n";
        break;

      case ON_POINT2:
        cerr << "ON_POINT2\n";
        break;

      case NONE:
        cerr << "NONE\n";
        break;

      default:
        cerr << "<Unknown>: " << pred << "\n";

    }
}
*/

 
/**  @brief check the graph structure after construction
 */
TEST_F(Voronoi3SimplexGraphTest, Test01) {
    Voronoi3SimplexGraph g01;


    EXPECT_EQ(g01.V1.mType, Voronoi3SimplexNode::VERTEX);
    EXPECT_EQ(g01.V1.mVertexSet, 1);
    EXPECT_EQ(g01.V1.mNeighbors[0], &g01.E12_V1);
    EXPECT_EQ(g01.V1.mNeighbors[1], &g01.E31_V1);
    EXPECT_EQ(g01.V1.mNeighbors[2], &g01.E14_V1);
    EXPECT_EQ(g01.V1.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.V1.mInward, 0);
    EXPECT_EQ(g01.V1.mInvisible, false);
    EXPECT_EQ(g01.V1.mQueued, false);


    EXPECT_EQ(g01.V2.mType, Voronoi3SimplexNode::VERTEX);
    EXPECT_EQ(g01.V2.mVertexSet, 2);
    EXPECT_EQ(g01.V2.mNeighbors[0], &g01.E12_V2);
    EXPECT_EQ(g01.V2.mNeighbors[1], &g01.E23_V2);
    EXPECT_EQ(g01.V2.mNeighbors[2], &g01.E24_V2);
    EXPECT_EQ(g01.V2.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.V2.mInward, 0);
    EXPECT_EQ(g01.V2.mInvisible, false);
    EXPECT_EQ(g01.V2.mQueued, false);


    EXPECT_EQ(g01.V3.mType, Voronoi3SimplexNode::VERTEX);
    EXPECT_EQ(g01.V3.mVertexSet, 3);
    EXPECT_EQ(g01.V3.mNeighbors[0], &g01.E23_V3);
    EXPECT_EQ(g01.V3.mNeighbors[1], &g01.E31_V3);
    EXPECT_EQ(g01.V3.mNeighbors[2], &g01.E34_V3);
    EXPECT_EQ(g01.V3.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.V3.mInward, 0);
    EXPECT_EQ(g01.V3.mInvisible, false);
    EXPECT_EQ(g01.V3.mQueued, false);


    EXPECT_EQ(g01.V4.mType, Voronoi3SimplexNode::VERTEX);
    EXPECT_EQ(g01.V4.mVertexSet, 4);
    EXPECT_EQ(g01.V4.mNeighbors[0], &g01.E14_V4);
    EXPECT_EQ(g01.V4.mNeighbors[1], &g01.E24_V4);
    EXPECT_EQ(g01.V4.mNeighbors[2], &g01.E34_V4);
    EXPECT_EQ(g01.V4.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.V4.mInward, 0);
    EXPECT_EQ(g01.V4.mInvisible, false);
    EXPECT_EQ(g01.V4.mQueued, false);


    EXPECT_EQ(g01.E12.mType, Voronoi3SimplexNode::EDGE);
    EXPECT_EQ(g01.E12.mVertexSet, 12);
    EXPECT_EQ(g01.E12.mNeighbors[0], &g01.F132_E12);
    EXPECT_EQ(g01.E12.mNeighbors[1], &g01.F124_E12);
    EXPECT_EQ(g01.E12.mNeighbors[2], &g01.E12_V1);
    EXPECT_EQ(g01.E12.mNeighbors[3], &g01.E12_V2);
    EXPECT_EQ(g01.E12.mInward, 0);
    EXPECT_EQ(g01.E12.mInvisible, false);
    EXPECT_EQ(g01.E12.mQueued, false);


    EXPECT_EQ(g01.E23.mType, Voronoi3SimplexNode::EDGE);
    EXPECT_EQ(g01.E23.mVertexSet, 23);
    EXPECT_EQ(g01.E23.mNeighbors[0], &g01.F132_E23);
    EXPECT_EQ(g01.E23.mNeighbors[1], &g01.F234_E23);
    EXPECT_EQ(g01.E23.mNeighbors[2], &g01.E23_V2);
    EXPECT_EQ(g01.E23.mNeighbors[3], &g01.E23_V3);
    EXPECT_EQ(g01.E23.mInward, 0);
    EXPECT_EQ(g01.E23.mInvisible, false);
    EXPECT_EQ(g01.E23.mQueued, false);


    EXPECT_EQ(g01.E31.mType, Voronoi3SimplexNode::EDGE);
    EXPECT_EQ(g01.E31.mVertexSet, 31);
    EXPECT_EQ(g01.E31.mNeighbors[0], &g01.F132_E31);
    EXPECT_EQ(g01.E31.mNeighbors[1], &g01.F143_E31);
    EXPECT_EQ(g01.E31.mNeighbors[2], &g01.E31_V3);
    EXPECT_EQ(g01.E31.mNeighbors[3], &g01.E31_V1);
    EXPECT_EQ(g01.E31.mInward, 0);
    EXPECT_EQ(g01.E31.mInvisible, false);
    EXPECT_EQ(g01.E31.mQueued, false);


    EXPECT_EQ(g01.E14.mType, Voronoi3SimplexNode::EDGE);
    EXPECT_EQ(g01.E14.mVertexSet, 14);
    EXPECT_EQ(g01.E14.mNeighbors[0], &g01.F124_E14);
    EXPECT_EQ(g01.E14.mNeighbors[1], &g01.F143_E14);
    EXPECT_EQ(g01.E14.mNeighbors[2], &g01.E14_V1);
    EXPECT_EQ(g01.E14.mNeighbors[3], &g01.E14_V4);
    EXPECT_EQ(g01.E14.mInward, 0);
    EXPECT_EQ(g01.E14.mInvisible, false);
    EXPECT_EQ(g01.E14.mQueued, false);


    EXPECT_EQ(g01.E24.mType, Voronoi3SimplexNode::EDGE);
    EXPECT_EQ(g01.E24.mVertexSet, 24);
    EXPECT_EQ(g01.E24.mNeighbors[0], &g01.F124_E24);
    EXPECT_EQ(g01.E24.mNeighbors[1], &g01.F234_E24);
    EXPECT_EQ(g01.E24.mNeighbors[2], &g01.E24_V2);
    EXPECT_EQ(g01.E24.mNeighbors[3], &g01.E24_V4);
    EXPECT_EQ(g01.E24.mInward, 0);
    EXPECT_EQ(g01.E24.mInvisible, false);
    EXPECT_EQ(g01.E24.mQueued, false);


    EXPECT_EQ(g01.E34.mType, Voronoi3SimplexNode::EDGE);
    EXPECT_EQ(g01.E34.mVertexSet, 34);
    EXPECT_EQ(g01.E34.mNeighbors[0], &g01.F143_E34);
    EXPECT_EQ(g01.E34.mNeighbors[1], &g01.F234_E34);
    EXPECT_EQ(g01.E34.mNeighbors[2], &g01.E34_V3);
    EXPECT_EQ(g01.E34.mNeighbors[3], &g01.E34_V4);
    EXPECT_EQ(g01.E34.mInward, 0);
    EXPECT_EQ(g01.E34.mInvisible, false);
    EXPECT_EQ(g01.E34.mQueued, false);


    EXPECT_EQ(g01.F132.mType, Voronoi3SimplexNode::FACE);
    EXPECT_EQ(g01.F132.mVertexSet, 132);
    EXPECT_EQ(g01.F132.mNeighbors[0], &g01.F132_E12);
    EXPECT_EQ(g01.F132.mNeighbors[1], &g01.F132_E23);
    EXPECT_EQ(g01.F132.mNeighbors[2], &g01.F132_E31);
    EXPECT_EQ(g01.F132.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.F132.mInward, 0);
    EXPECT_EQ(g01.F132.mInvisible, false);
    EXPECT_EQ(g01.F132.mQueued, false);


    EXPECT_EQ(g01.F124.mType, Voronoi3SimplexNode::FACE);
    EXPECT_EQ(g01.F124.mVertexSet, 124);
    EXPECT_EQ(g01.F124.mNeighbors[0], &g01.F124_E12);
    EXPECT_EQ(g01.F124.mNeighbors[1], &g01.F124_E24);
    EXPECT_EQ(g01.F124.mNeighbors[2], &g01.F124_E14);
    EXPECT_EQ(g01.F124.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.F124.mInward, 0);
    EXPECT_EQ(g01.F124.mInvisible, false);
    EXPECT_EQ(g01.F124.mQueued, false);


    EXPECT_EQ(g01.F143.mType, Voronoi3SimplexNode::FACE);
    EXPECT_EQ(g01.F143.mVertexSet, 143);
    EXPECT_EQ(g01.F143.mNeighbors[0], &g01.F143_E31);
    EXPECT_EQ(g01.F143.mNeighbors[1], &g01.F143_E14);
    EXPECT_EQ(g01.F143.mNeighbors[2], &g01.F143_E34);
    EXPECT_EQ(g01.F143.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.F143.mInward, 0);
    EXPECT_EQ(g01.F143.mInvisible, false);
    EXPECT_EQ(g01.F143.mQueued, false);


    EXPECT_EQ(g01.F234.mType, Voronoi3SimplexNode::FACE);
    EXPECT_EQ(g01.F234.mVertexSet, 234);
    EXPECT_EQ(g01.F234.mNeighbors[0], &g01.F234_E23);
    EXPECT_EQ(g01.F234.mNeighbors[1], &g01.F234_E34);
    EXPECT_EQ(g01.F234.mNeighbors[2], &g01.F234_E24);
    EXPECT_EQ(g01.F234.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.F234.mInward, 0);
    EXPECT_EQ(g01.F234.mInvisible, false);
    EXPECT_EQ(g01.F234.mQueued, false);


}


/**  @brief test reset()
 */
TEST_F(Voronoi3SimplexGraphTest, Test02) {

    Voronoi3SimplexGraph g01;

    g01.V1.mInward    = 99;
    g01.V1.mInvisible = true;
    g01.V1.mQueued    = true;

    g01.V2.mInward    = 99;
    g01.V2.mInvisible = true;
    g01.V2.mQueued    = true;

    g01.V3.mInward    = 99;
    g01.V3.mInvisible = true;
    g01.V3.mQueued    = true;

    g01.V4.mInward    = 99;
    g01.V4.mInvisible = true;
    g01.V4.mQueued    = true;

    g01.E12.mInward    = 99;
    g01.E12.mInvisible = true;
    g01.E12.mQueued    = true;

    g01.E23.mInward    = 99;
    g01.E23.mInvisible = true;
    g01.E23.mQueued    = true;

    g01.E31.mInward    = 99;
    g01.E31.mInvisible = true;
    g01.E31.mQueued    = true;

    g01.E14.mInward    = 99;
    g01.E14.mInvisible = true;
    g01.E14.mQueued    = true;

    g01.E24.mInward    = 99;
    g01.E24.mInvisible = true;
    g01.E24.mQueued    = true;

    g01.E34.mInward    = 99;
    g01.E34.mInvisible = true;
    g01.E34.mQueued    = true;

    g01.F132.mInward    = 99;
    g01.F132.mInvisible = true;
    g01.F132.mQueued    = true;

    g01.F124.mInward    = 99;
    g01.F124.mInvisible = true;
    g01.F124.mQueued    = true;

    g01.F143.mInward    = 99;
    g01.F143.mInvisible = true;
    g01.F143.mQueued    = true;

    g01.F234.mInward    = 99;
    g01.F234.mInvisible = true;
    g01.F234.mQueued    = true;

    Vec3 p01(1.0, 0.0, 0.0);
    Vec3 p02(0.0, 0.0, 0.0);
    Vec3 p03(0.0, 0.0, 1.0);
    Vec3 p04(0.0, 1.0, 0.0);

    Vec3 pt01(2.0, 2.0, 2.0);

    g01.reset(p01, p02, p03, p04, pt01);

    EXPECT_EQ(g01.mP1, p01);
    EXPECT_EQ(g01.mP2, p02);
    EXPECT_EQ(g01.mP3, p03);
    EXPECT_EQ(g01.mP4, p04);
    EXPECT_EQ(g01.mPT, pt01);

    EXPECT_EQ(g01.V1.mType, Voronoi3SimplexNode::VERTEX);
    EXPECT_EQ(g01.V1.mVertexSet, 1);
    EXPECT_EQ(g01.V1.mNeighbors[0], &g01.E12_V1);
    EXPECT_EQ(g01.V1.mNeighbors[1], &g01.E31_V1);
    EXPECT_EQ(g01.V1.mNeighbors[2], &g01.E14_V1);
    EXPECT_EQ(g01.V1.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.V1.mInward, 0);
    EXPECT_EQ(g01.V1.mInvisible, false);
    EXPECT_EQ(g01.V1.mQueued, false);


    EXPECT_EQ(g01.V2.mType, Voronoi3SimplexNode::VERTEX);
    EXPECT_EQ(g01.V2.mVertexSet, 2);
    EXPECT_EQ(g01.V2.mNeighbors[0], &g01.E12_V2);
    EXPECT_EQ(g01.V2.mNeighbors[1], &g01.E23_V2);
    EXPECT_EQ(g01.V2.mNeighbors[2], &g01.E24_V2);
    EXPECT_EQ(g01.V2.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.V2.mInward, 0);
    EXPECT_EQ(g01.V2.mInvisible, false);
    EXPECT_EQ(g01.V2.mQueued, false);


    EXPECT_EQ(g01.V3.mType, Voronoi3SimplexNode::VERTEX);
    EXPECT_EQ(g01.V3.mVertexSet, 3);
    EXPECT_EQ(g01.V3.mNeighbors[0], &g01.E23_V3);
    EXPECT_EQ(g01.V3.mNeighbors[1], &g01.E31_V3);
    EXPECT_EQ(g01.V3.mNeighbors[2], &g01.E34_V3);
    EXPECT_EQ(g01.V3.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.V3.mInward, 0);
    EXPECT_EQ(g01.V3.mInvisible, false);
    EXPECT_EQ(g01.V3.mQueued, false);


    EXPECT_EQ(g01.V4.mType, Voronoi3SimplexNode::VERTEX);
    EXPECT_EQ(g01.V4.mVertexSet, 4);
    EXPECT_EQ(g01.V4.mNeighbors[0], &g01.E14_V4);
    EXPECT_EQ(g01.V4.mNeighbors[1], &g01.E24_V4);
    EXPECT_EQ(g01.V4.mNeighbors[2], &g01.E34_V4);
    EXPECT_EQ(g01.V4.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.V4.mInward, 0);
    EXPECT_EQ(g01.V4.mInvisible, false);
    EXPECT_EQ(g01.V4.mQueued, false);


    EXPECT_EQ(g01.E12.mType, Voronoi3SimplexNode::EDGE);
    EXPECT_EQ(g01.E12.mVertexSet, 12);
    EXPECT_EQ(g01.E12.mNeighbors[0], &g01.F132_E12);
    EXPECT_EQ(g01.E12.mNeighbors[1], &g01.F124_E12);
    EXPECT_EQ(g01.E12.mNeighbors[2], &g01.E12_V1);
    EXPECT_EQ(g01.E12.mNeighbors[3], &g01.E12_V2);
    EXPECT_EQ(g01.E12.mInward, 0);
    EXPECT_EQ(g01.E12.mInvisible, false);
    EXPECT_EQ(g01.E12.mQueued, false);


    EXPECT_EQ(g01.E23.mType, Voronoi3SimplexNode::EDGE);
    EXPECT_EQ(g01.E23.mVertexSet, 23);
    EXPECT_EQ(g01.E23.mNeighbors[0], &g01.F132_E23);
    EXPECT_EQ(g01.E23.mNeighbors[1], &g01.F234_E23);
    EXPECT_EQ(g01.E23.mNeighbors[2], &g01.E23_V2);
    EXPECT_EQ(g01.E23.mNeighbors[3], &g01.E23_V3);
    EXPECT_EQ(g01.E23.mInward, 0);
    EXPECT_EQ(g01.E23.mInvisible, false);
    EXPECT_EQ(g01.E23.mQueued, false);


    EXPECT_EQ(g01.E31.mType, Voronoi3SimplexNode::EDGE);
    EXPECT_EQ(g01.E31.mVertexSet, 31);
    EXPECT_EQ(g01.E31.mNeighbors[0], &g01.F132_E31);
    EXPECT_EQ(g01.E31.mNeighbors[1], &g01.F143_E31);
    EXPECT_EQ(g01.E31.mNeighbors[2], &g01.E31_V3);
    EXPECT_EQ(g01.E31.mNeighbors[3], &g01.E31_V1);
    EXPECT_EQ(g01.E31.mInward, 0);
    EXPECT_EQ(g01.E31.mInvisible, false);
    EXPECT_EQ(g01.E31.mQueued, false);


    EXPECT_EQ(g01.E14.mType, Voronoi3SimplexNode::EDGE);
    EXPECT_EQ(g01.E14.mVertexSet, 14);
    EXPECT_EQ(g01.E14.mNeighbors[0], &g01.F124_E14);
    EXPECT_EQ(g01.E14.mNeighbors[1], &g01.F143_E14);
    EXPECT_EQ(g01.E14.mNeighbors[2], &g01.E14_V1);
    EXPECT_EQ(g01.E14.mNeighbors[3], &g01.E14_V4);
    EXPECT_EQ(g01.E14.mInward, 0);
    EXPECT_EQ(g01.E14.mInvisible, false);
    EXPECT_EQ(g01.E14.mQueued, false);


    EXPECT_EQ(g01.E24.mType, Voronoi3SimplexNode::EDGE);
    EXPECT_EQ(g01.E24.mVertexSet, 24);
    EXPECT_EQ(g01.E24.mNeighbors[0], &g01.F124_E24);
    EXPECT_EQ(g01.E24.mNeighbors[1], &g01.F234_E24);
    EXPECT_EQ(g01.E24.mNeighbors[2], &g01.E24_V2);
    EXPECT_EQ(g01.E24.mNeighbors[3], &g01.E24_V4);
    EXPECT_EQ(g01.E24.mInward, 0);
    EXPECT_EQ(g01.E24.mInvisible, false);
    EXPECT_EQ(g01.E24.mQueued, false);


    EXPECT_EQ(g01.E34.mType, Voronoi3SimplexNode::EDGE);
    EXPECT_EQ(g01.E34.mVertexSet, 34);
    EXPECT_EQ(g01.E34.mNeighbors[0], &g01.F143_E34);
    EXPECT_EQ(g01.E34.mNeighbors[1], &g01.F234_E34);
    EXPECT_EQ(g01.E34.mNeighbors[2], &g01.E34_V3);
    EXPECT_EQ(g01.E34.mNeighbors[3], &g01.E34_V4);
    EXPECT_EQ(g01.E34.mInward, 0);
    EXPECT_EQ(g01.E34.mInvisible, false);
    EXPECT_EQ(g01.E34.mQueued, false);


    EXPECT_EQ(g01.F132.mType, Voronoi3SimplexNode::FACE);
    EXPECT_EQ(g01.F132.mVertexSet, 132);
    EXPECT_EQ(g01.F132.mNeighbors[0], &g01.F132_E12);
    EXPECT_EQ(g01.F132.mNeighbors[1], &g01.F132_E23);
    EXPECT_EQ(g01.F132.mNeighbors[2], &g01.F132_E31);
    EXPECT_EQ(g01.F132.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.F132.mInward, 0);
    EXPECT_EQ(g01.F132.mInvisible, false);
    EXPECT_EQ(g01.F132.mQueued, false);


    EXPECT_EQ(g01.F124.mType, Voronoi3SimplexNode::FACE);
    EXPECT_EQ(g01.F124.mVertexSet, 124);
    EXPECT_EQ(g01.F124.mNeighbors[0], &g01.F124_E12);
    EXPECT_EQ(g01.F124.mNeighbors[1], &g01.F124_E24);
    EXPECT_EQ(g01.F124.mNeighbors[2], &g01.F124_E14);
    EXPECT_EQ(g01.F124.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.F124.mInward, 0);
    EXPECT_EQ(g01.F124.mInvisible, false);
    EXPECT_EQ(g01.F124.mQueued, false);


    EXPECT_EQ(g01.F143.mType, Voronoi3SimplexNode::FACE);
    EXPECT_EQ(g01.F143.mVertexSet, 143);
    EXPECT_EQ(g01.F143.mNeighbors[0], &g01.F143_E31);
    EXPECT_EQ(g01.F143.mNeighbors[1], &g01.F143_E14);
    EXPECT_EQ(g01.F143.mNeighbors[2], &g01.F143_E34);
    EXPECT_EQ(g01.F143.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.F143.mInward, 0);
    EXPECT_EQ(g01.F143.mInvisible, false);
    EXPECT_EQ(g01.F143.mQueued, false);


    EXPECT_EQ(g01.F234.mType, Voronoi3SimplexNode::FACE);
    EXPECT_EQ(g01.F234.mVertexSet, 234);
    EXPECT_EQ(g01.F234.mNeighbors[0], &g01.F234_E23);
    EXPECT_EQ(g01.F234.mNeighbors[1], &g01.F234_E34);
    EXPECT_EQ(g01.F234.mNeighbors[2], &g01.F234_E24);
    EXPECT_EQ(g01.F234.mNeighbors[3], nullptr);
    EXPECT_EQ(g01.F234.mInward, 0);
    EXPECT_EQ(g01.F234.mInvisible, false);
    EXPECT_EQ(g01.F234.mQueued, false);


}


/** @brief FindNormalsAndTestVisibility()
 */
TEST_F(Voronoi3SimplexGraphTest, Test03) {

    Voronoi3SimplexGraph g01;

    for (long i = 0; i < 1000; i++) {

        Vec3 trans((rand()%100)/10.0, (rand()%100)/10.0, (rand()%100)/10.0);
        Mat3x3 M = randomRotMat();

        Vec3 p01( 0.0, 0.0, 0.0);
        Vec3 p02( ((rand()%100)+1)/10.0, 0.0, 0.0);
        Vec3 p03( 0.0, ((rand()%100)+1)/10.0, 0.0);
        Vec3 p04( 0.0, 0.0, ((rand()%100)+1)/10.0);

        Vec3 pt01(p02.x()/4.0, p03.y()/4.0, p04.z()/4.0);// Inside.
        Vec3 pt02(p02.x()/1.5, p03.y()/1.5, p04.z()/1.5);// Over 234 near
        Vec3 pt03(p02.x()*1000.0,p03.y()*1000.0,p04.z()*1000.0);//Over 234far
        Vec3 pt04(p02.x()/4.0,p03.y()/4.0,-0.0001); //Over 132 near
        Vec3 pt05(p02.x()/4.0,p03.y()/4.0,-10000.0);//Over 132 far
        Vec3 pt06(p02.x()/4.0,-0.0001, p04.z()/4.0);//Over 124 near
        Vec3 pt07(p02.x()/4.0,-10000.0,p04.z()/4.0);//Over 124 far
        Vec3 pt08(-0.0001, p03.y()/4.0,p04.z()/4.0);//Over 143 near
        Vec3 pt09(-10000.0,p03.y()/4.0,p04.z()/4.0);//Over 143 far
        Vec3 pt10(-0.0001, -0.0001, -0.0001);       //Over 132,124,143 near
        Vec3 pt11(-10000.0, -10000.0, -10000.0);    //Over 132,124,143 near
        Vec3 pt12(-0.0001,  p03.y()/2.0, -0.0001);  //Over 132,143 near
        Vec3 pt13(-10000.0, p03.y()/2.0, -10000.0); //Over 132,143 far
        Vec3 pt14(p02.x()/2.0, -0.0001, -0.0001);   //Over 132,124 near
        Vec3 pt15(p02.x()/2.0, -10000.0,-10000.0);  //Over 132,124 far
        Vec3 pt16(-0.0001, -0.0001,   p04.z()/2.0); //Over 143,124 near
        Vec3 pt17(-10000.0,-10000.0,  p04.z()/2.0); //Over 143,124 far
        Vec3 pt18(p02.x(), p03.y(), -0.1*p04.z());  //Over 132,234
        Vec3 pt19(p02.x(), -0.1*p03.y(), p04.z());  //Over 124,234
        Vec3 pt20(-0.1*p02.x(), p03.y(), p04.z());  //Over 143,234
        Vec3 pt21(2.0*p02.x(),-0.1*p03.y(),-0.1*p04.z());//Over 132,124,234
        Vec3 pt22(-0.1*p02.x(),2.0*p03.y(),-0.1*p04.z());//Over 132,143,234
        Vec3 pt23(-0.1*p02.x(),-0.1*p03.y(),2.0*p04.z());//Over 124,143,234

        p01  = M * p01  + trans;
        p02  = M * p02  + trans;
        p03  = M * p03  + trans;
        p04  = M * p04  + trans;
        pt01 = M * pt01 + trans;
        pt02 = M * pt02 + trans;
        pt03 = M * pt03 + trans;
        pt04 = M * pt04 + trans;
        pt05 = M * pt05 + trans;
        pt06 = M * pt06 + trans;
        pt07 = M * pt07 + trans;
        pt08 = M * pt08 + trans;
        pt09 = M * pt09 + trans;
        pt10 = M * pt10 + trans;
        pt11 = M * pt11 + trans;
        pt12 = M * pt12 + trans;
        pt13 = M * pt13 + trans;
        pt14 = M * pt14 + trans;
        pt15 = M * pt15 + trans;
        pt16 = M * pt16 + trans;
        pt17 = M * pt17 + trans;
        pt18 = M * pt18 + trans;
        pt19 = M * pt19 + trans;
        pt20 = M * pt20 + trans;
        pt21 = M * pt21 + trans;
        pt22 = M * pt22 + trans;
        pt23 = M * pt23 + trans;

        g01.reset(p01, p02, p03, p04, pt01);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, true);
        EXPECT_EQ(g01.V2.mInvisible, true);
        EXPECT_EQ(g01.V3.mInvisible, true);
        EXPECT_EQ(g01.V4.mInvisible, true);
        EXPECT_EQ(g01.E12.mInvisible, true);
        EXPECT_EQ(g01.E23.mInvisible, true);
        EXPECT_EQ(g01.E31.mInvisible, true);
        EXPECT_EQ(g01.E14.mInvisible, true);
        EXPECT_EQ(g01.E24.mInvisible, true);
        EXPECT_EQ(g01.E34.mInvisible, true);
        EXPECT_EQ(g01.F132.mInvisible, true);
        EXPECT_EQ(g01.F124.mInvisible, true);
        EXPECT_EQ(g01.F143.mInvisible, true);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt02);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, true);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, true);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, true);
        EXPECT_EQ(g01.E14.mInvisible, true);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, true);
        EXPECT_EQ(g01.F124.mInvisible, true);
        EXPECT_EQ(g01.F143.mInvisible, true);
        EXPECT_EQ(g01.F234.mInvisible, false);

        g01.reset(p01, p02, p03, p04, pt03);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, true);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, true);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, true);
        EXPECT_EQ(g01.E14.mInvisible, true);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, true);
        EXPECT_EQ(g01.F124.mInvisible, true);
        EXPECT_EQ(g01.F143.mInvisible, true);
        EXPECT_EQ(g01.F234.mInvisible, false);

        g01.reset(p01, p02, p03, p04, pt04);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, true);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, true);
        EXPECT_EQ(g01.E24.mInvisible, true);
        EXPECT_EQ(g01.E34.mInvisible, true);
        EXPECT_EQ(g01.F132.mInvisible, false);
        EXPECT_EQ(g01.F124.mInvisible, true);
        EXPECT_EQ(g01.F143.mInvisible, true);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt05);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, true);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, true);
        EXPECT_EQ(g01.E24.mInvisible, true);
        EXPECT_EQ(g01.E34.mInvisible, true);
        EXPECT_EQ(g01.F132.mInvisible, false);
        EXPECT_EQ(g01.F124.mInvisible, true);
        EXPECT_EQ(g01.F143.mInvisible, true);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt06);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, true);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, true);
        EXPECT_EQ(g01.E31.mInvisible, true);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, true);
        EXPECT_EQ(g01.F132.mInvisible, true);
        EXPECT_EQ(g01.F124.mInvisible, false);
        EXPECT_EQ(g01.F143.mInvisible, true);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt07);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, true);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, true);
        EXPECT_EQ(g01.E31.mInvisible, true);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, true);
        EXPECT_EQ(g01.F132.mInvisible, true);
        EXPECT_EQ(g01.F124.mInvisible, false);
        EXPECT_EQ(g01.F143.mInvisible, true);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt08);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, true);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, true);
        EXPECT_EQ(g01.E23.mInvisible, true);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, true);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, true);
        EXPECT_EQ(g01.F124.mInvisible, true);
        EXPECT_EQ(g01.F143.mInvisible, false);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt09);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, true);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, true);
        EXPECT_EQ(g01.E23.mInvisible, true);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, true);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, true);
        EXPECT_EQ(g01.F124.mInvisible, true);
        EXPECT_EQ(g01.F143.mInvisible, false);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt10);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, false);
        EXPECT_EQ(g01.F124.mInvisible, false);
        EXPECT_EQ(g01.F143.mInvisible, false);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt11);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, false);
        EXPECT_EQ(g01.F124.mInvisible, false);
        EXPECT_EQ(g01.F143.mInvisible, false);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt12);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, true);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, false);
        EXPECT_EQ(g01.F124.mInvisible, true);
        EXPECT_EQ(g01.F143.mInvisible, false);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt13);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, true);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, false);
        EXPECT_EQ(g01.F124.mInvisible, true);
        EXPECT_EQ(g01.F143.mInvisible, false);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt14);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, true);
        EXPECT_EQ(g01.F132.mInvisible, false);
        EXPECT_EQ(g01.F124.mInvisible, false);
        EXPECT_EQ(g01.F143.mInvisible, true);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt15);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, true);
        EXPECT_EQ(g01.F132.mInvisible, false);
        EXPECT_EQ(g01.F124.mInvisible, false);
        EXPECT_EQ(g01.F143.mInvisible, true);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt16);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, true);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, true);
        EXPECT_EQ(g01.F124.mInvisible, false);
        EXPECT_EQ(g01.F143.mInvisible, false);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt17);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, true);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, true);
        EXPECT_EQ(g01.F124.mInvisible, false);
        EXPECT_EQ(g01.F143.mInvisible, false);
        EXPECT_EQ(g01.F234.mInvisible, true);

        g01.reset(p01, p02, p03, p04, pt18);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, true);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, false);
        EXPECT_EQ(g01.F124.mInvisible, true);
        EXPECT_EQ(g01.F143.mInvisible, true);
        EXPECT_EQ(g01.F234.mInvisible, false);

        g01.reset(p01, p02, p03, p04, pt19);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, true);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, true);
        EXPECT_EQ(g01.F124.mInvisible, false);
        EXPECT_EQ(g01.F143.mInvisible, true);
        EXPECT_EQ(g01.F234.mInvisible, false);

        g01.reset(p01, p02, p03, p04, pt20);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, true);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, true);
        EXPECT_EQ(g01.F124.mInvisible, true);
        EXPECT_EQ(g01.F143.mInvisible, false);
        EXPECT_EQ(g01.F234.mInvisible, false);

        g01.reset(p01, p02, p03, p04, pt21);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, false);
        EXPECT_EQ(g01.F124.mInvisible, false);
        EXPECT_EQ(g01.F143.mInvisible, true);
        EXPECT_EQ(g01.F234.mInvisible, false);

        g01.reset(p01, p02, p03, p04, pt22);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, false);
        EXPECT_EQ(g01.F124.mInvisible, true);
        EXPECT_EQ(g01.F143.mInvisible, false);
        EXPECT_EQ(g01.F234.mInvisible, false);

        g01.reset(p01, p02, p03, p04, pt23);
        g01.FindNormalsAndTestVisibility();
        EXPECT_EQ(g01.V1.mInvisible, false);
        EXPECT_EQ(g01.V2.mInvisible, false);
        EXPECT_EQ(g01.V3.mInvisible, false);
        EXPECT_EQ(g01.V4.mInvisible, false);
        EXPECT_EQ(g01.E12.mInvisible, false);
        EXPECT_EQ(g01.E23.mInvisible, false);
        EXPECT_EQ(g01.E31.mInvisible, false);
        EXPECT_EQ(g01.E14.mInvisible, false);
        EXPECT_EQ(g01.E24.mInvisible, false);
        EXPECT_EQ(g01.E34.mInvisible, false);
        EXPECT_EQ(g01.F132.mInvisible, true);
        EXPECT_EQ(g01.F124.mInvisible, false);
        EXPECT_EQ(g01.F143.mInvisible, false);
        EXPECT_EQ(g01.F234.mInvisible, false);

    }

}


static void generateRandomAlphas(
    double& a01,
    double& a02,
    double& a03,
    double& a04
) {
    double v1 = ((rand()%10000)+1)/10000.0;
    double v2 = ((rand()%10000)+1)/10000.0;
    double v3 = ((rand()%10000)+1)/10000.0;
    double v4 = ((rand()%10000)+1)/10000.0;
    double sum =v1 + v2 + v3 + v4;    
    a01 = v1/sum;
    a02 = v2/sum;
    a03 = v3/sum;
    a04 = v4/sum;
}


static void generateRandomAlphas(
    double& a01,
    double& a02,
    double& a03
) {
    double v1 = ((rand()%10000)+1)/10000.0;
    double v2 = ((rand()%10000)+1)/10000.0;
    double v3 = ((rand()%10000)+1)/10000.0;
    double sum =v1 + v2 + v3;
    a01 = v1/sum;
    a02 = v2/sum;
    a03 = v3/sum;
}


static void generateRandomAlphas(
    double& a01,
    double& a02
) {
    double v1 = ((rand()%10000)+1)/10000.0;
    double v2 = ((rand()%10000)+1)/10000.0;
    double sum =v1 + v2;
    a01 = v1/sum;
    a02 = v2/sum;
}


static double rand100()
{
    return ((rand()%10000)+1)/100.0;
}


static bool isAnyOf(
    enum predicate ptest,
    enum predicate p1,
    enum predicate p2,
    enum predicate p3,
    enum predicate p4    )
{ return ptest == p1 || ptest == p2 || ptest == p3 || ptest == p4; }

/*
static bool isAnyOf(
    enum predicate ptest,
    enum predicate p1,
    enum predicate p2,
    enum predicate p3 )
{ return ptest == p1 || ptest == p2 || ptest == p3; }
*/

static bool isAnyOf(
    enum predicate ptest,
    enum predicate p1,
    enum predicate p2 )
{ return ptest == p1 || ptest == p2; }



/** @brief tests find()
 */
TEST_F(Voronoi3SimplexGraphTest, Test04) {

    Voronoi3SimplexGraph g01;

    for (long i = 0; i < 100000; i++) {

        Vec3 trans((rand()%100)/10.0, (rand()%100)/10.0, (rand()%100)/10.0);
        Mat3x3 M = randomRotMat();

        Vec3 p01( 0.0, 0.0, 0.0);
        Vec3 p02( ((rand()%100)+1)/10.0, 0.0, 0.0);
        Vec3 p03( 0.0, ((rand()%100)+1)/10.0, 0.0);
        Vec3 p04( 0.0, 0.0, ((rand()%100)+1)/10.0);

        Vec3 n132 = (p03 - p01).cross(p02 - p01);
        n132.normalize();

        Vec3 n124 = (p02 - p01).cross(p04 - p01);
        n124.normalize();

        Vec3 n143 = (p04 - p01).cross(p03 - p01);
        n143.normalize();

        Vec3 n234 = (p03 - p02).cross(p04 - p02);
        n234.normalize();

        double a01, a02, a03, a04;
        generateRandomAlphas(a01, a02, a03, a04);
        Vec3 pt01 = p01 * a01 + p02 * a02 + p03 * a03 + p04 * a04;

        generateRandomAlphas(a01, a02, a03);
        Vec3 pt02 = p02*a01 + p03*a02 + p04*a03 + n234*rand100();

        generateRandomAlphas(a01, a02, a03);
        Vec3 pt03 = p01*a01 + p03*a02 + p02*a03 + n132*rand100();

        generateRandomAlphas(a01, a02, a03);
        Vec3 pt04 = p01*a01 + p02*a02 + p04*a03 + n124*rand100();

        generateRandomAlphas(a01, a02, a03);
        Vec3 pt05 = p01*a01 + p04*a02 + p03*a03 + n143*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt06 = p01*a01 + p02*a02 + n132*rand100() + n124*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt07 = p02*a01 + p03*a02 + n132*rand100() + n234*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt08 = p03*a01 + p01*a02 + n132*rand100() + n143*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt09 = p01*a01 + p04*a02 + n124*rand100() + n143*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt10 = p02*a01 + p04*a02 + n124*rand100() + n234*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt11 = p03*a01 + p04*a02 + n143*rand100() + n234*rand100();

        Vec3 pt12 = p01 + n132*rand100() + n124*rand100() + n143*rand100();
        Vec3 pt13 = p02 + n132*rand100() + n124*rand100() + n234*rand100();
        Vec3 pt14 = p03 + n132*rand100() + n143*rand100() + n234*rand100();
        Vec3 pt15 = p04 + n124*rand100() + n143*rand100() + n234*rand100();

        Vec3 pt16 = p01 + n132*rand100();
        Vec3 pt17 = p01 + n124*rand100();
        Vec3 pt18 = p01 + n143*rand100();

        Vec3 pt19 = p02 + n132*rand100();
        Vec3 pt20 = p02 + n124*rand100();
        Vec3 pt21 = p02 + n234*rand100();

        Vec3 pt22 = p03 + n132*rand100();
        Vec3 pt23 = p03 + n143*rand100();
        Vec3 pt24 = p03 + n234*rand100();

        Vec3 pt25 = p04 + n124*rand100();
        Vec3 pt26 = p04 + n143*rand100();
        Vec3 pt27 = p04 + n234*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt28 = p01*a01 + p02*a02 + n132*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt29 = p01*a01 + p02*a02 + n124*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt30 = p02*a01 + p03*a02 + n132*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt31 = p02*a01 + p03*a02 + n234*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt32 = p03*a01 + p01*a02 + n132*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt33 = p03*a01 + p01*a02 + n143*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt34 = p01*a01 + p04*a02 + n124*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt35 = p01*a01 + p04*a02 + n143*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt36 = p02*a01 + p04*a02 + n124*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt37 = p02*a01 + p04*a02 + n234*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt38 = p03*a01 + p04*a02 + n143*rand100();

        generateRandomAlphas(a01, a02);
        Vec3 pt39 = p03*a01 + p04*a02 + n234*rand100();

        p01  = M * p01  + trans;
        p02  = M * p02  + trans;
        p03  = M * p03  + trans;
        p04  = M * p04  + trans;
        pt01 = M * pt01 + trans;
        pt02 = M * pt02 + trans;
        pt03 = M * pt03 + trans;
        pt04 = M * pt04 + trans;
        pt05 = M * pt05 + trans;
        pt06 = M * pt06 + trans;
        pt07 = M * pt07 + trans;
        pt08 = M * pt08 + trans;
        pt09 = M * pt09 + trans;
        pt10 = M * pt10 + trans;
        pt11 = M * pt11 + trans;
        pt12 = M * pt12 + trans;
        pt13 = M * pt13 + trans;
        pt14 = M * pt14 + trans;
        pt15 = M * pt15 + trans;
        pt16 = M * pt16 + trans;
        pt17 = M * pt17 + trans;
        pt18 = M * pt18 + trans;
        pt19 = M * pt19 + trans;
        pt20 = M * pt20 + trans;
        pt21 = M * pt21 + trans;
        pt22 = M * pt22 + trans;
        pt23 = M * pt23 + trans;
        pt24 = M * pt24 + trans;
        pt25 = M * pt25 + trans;
        pt26 = M * pt26 + trans;
        pt27 = M * pt27 + trans;
        pt28 = M * pt28 + trans;
        pt29 = M * pt29 + trans;
        pt30 = M * pt30 + trans;
        pt31 = M * pt31 + trans;
        pt32 = M * pt32 + trans;
        pt33 = M * pt33 + trans;
        pt34 = M * pt34 + trans;
        pt35 = M * pt35 + trans;
        pt36 = M * pt36 + trans;
        pt37 = M * pt37 + trans;
        pt38 = M * pt38 + trans;
        pt39 = M * pt39 + trans;

        auto res01 = g01.find(p01, p02, p03, p04, pt01);
        EXPECT_EQ(res01, VORONOI_INSIDE_TETRAHEDRON);

        auto res02 = g01.find(p01, p02, p03, p04, pt02);
        EXPECT_EQ(res02, VORONOI_OVER_TRIANGLE_2_3_4);

        auto res03 = g01.find(p01, p02, p03, p04, pt03);
        EXPECT_EQ(res03, VORONOI_OVER_TRIANGLE_1_3_2);

        auto res04 = g01.find(p01, p02, p03, p04, pt04);
        EXPECT_EQ(res04, VORONOI_OVER_TRIANGLE_1_2_4);

        auto res05 = g01.find(p01, p02, p03, p04, pt05);
        EXPECT_EQ(res05, VORONOI_OVER_TRIANGLE_1_4_3);

        auto res06 = g01.find(p01, p02, p03, p04, pt06);
        EXPECT_EQ(res06, VORONOI_OVER_EDGE_1_2);

        auto res07 = g01.find(p01, p02, p03, p04, pt07);
        EXPECT_EQ(res07, VORONOI_OVER_EDGE_2_3);

        auto res08 = g01.find(p01, p02, p03, p04, pt08);
        EXPECT_EQ(res08, VORONOI_OVER_EDGE_3_1);

        auto res09 = g01.find(p01, p02, p03, p04, pt09);
        EXPECT_EQ(res09, VORONOI_OVER_EDGE_1_4);

        auto res10 = g01.find(p01, p02, p03, p04, pt10);
        EXPECT_EQ(res10, VORONOI_OVER_EDGE_2_4);

        auto res11 = g01.find(p01, p02, p03, p04, pt11);
        EXPECT_EQ(res11, VORONOI_OVER_EDGE_3_4);

        auto res12 = g01.find(p01, p02, p03, p04, pt12);
        EXPECT_EQ(res12, VORONOI_OVER_VERTEX_1);

        auto res13 = g01.find(p01, p02, p03, p04, pt13);
        EXPECT_EQ(res13, VORONOI_OVER_VERTEX_2);

        auto res14 = g01.find(p01, p02, p03, p04, pt14);
        EXPECT_EQ(res14, VORONOI_OVER_VERTEX_3);

        auto res15 = g01.find(p01, p02, p03, p04, pt15);
        EXPECT_EQ(res15, VORONOI_OVER_VERTEX_4);

        auto res16 = g01.find(p01, p02, p03, p04, pt16);

        EXPECT_EQ(isAnyOf(res16, 
               VORONOI_OVER_VERTEX_1,
               VORONOI_OVER_EDGE_1_2,
               VORONOI_OVER_EDGE_3_1,
               VORONOI_OVER_TRIANGLE_1_3_2), true);


        auto res17 = g01.find(p01, p02, p03, p04, pt17);

        EXPECT_EQ(isAnyOf(res17, 
               VORONOI_OVER_VERTEX_1,
               VORONOI_OVER_EDGE_1_2,
               VORONOI_OVER_EDGE_1_4,
               VORONOI_OVER_TRIANGLE_1_2_4), true);


        auto res18 = g01.find(p01, p02, p03, p04, pt18);

        EXPECT_EQ(isAnyOf(res18, 
               VORONOI_OVER_VERTEX_1,
               VORONOI_OVER_EDGE_3_1,
               VORONOI_OVER_EDGE_1_4,
               VORONOI_OVER_TRIANGLE_1_4_3), true);

        auto res19 = g01.find(p01, p02, p03, p04, pt19);

        EXPECT_EQ(isAnyOf(res19, 
               VORONOI_OVER_VERTEX_2,
               VORONOI_OVER_EDGE_1_2,
               VORONOI_OVER_EDGE_2_3,
               VORONOI_OVER_TRIANGLE_1_3_2), true);

        auto res20 = g01.find(p01, p02, p03, p04, pt20);

        EXPECT_EQ(isAnyOf(res20, 
               VORONOI_OVER_VERTEX_2,
               VORONOI_OVER_EDGE_1_2,
               VORONOI_OVER_EDGE_2_4,
               VORONOI_OVER_TRIANGLE_1_2_4), true);

        auto res21 = g01.find(p01, p02, p03, p04, pt21);

        EXPECT_EQ(isAnyOf(res21, 
               VORONOI_OVER_VERTEX_2,
               VORONOI_OVER_EDGE_2_3,
               VORONOI_OVER_EDGE_2_4,
               VORONOI_OVER_TRIANGLE_2_3_4), true);

        auto res22 = g01.find(p01, p02, p03, p04, pt22);

        EXPECT_EQ(isAnyOf(res22, 
               VORONOI_OVER_VERTEX_3,
               VORONOI_OVER_EDGE_2_3,
               VORONOI_OVER_EDGE_3_1,
               VORONOI_OVER_TRIANGLE_1_3_2), true);

        auto res23 = g01.find(p01, p02, p03, p04, pt23);

        EXPECT_EQ(isAnyOf(res23, 
               VORONOI_OVER_VERTEX_3,
               VORONOI_OVER_EDGE_3_1,
               VORONOI_OVER_EDGE_3_4,
               VORONOI_OVER_TRIANGLE_1_4_3), true);

        auto res24 = g01.find(p01, p02, p03, p04, pt24);

        EXPECT_EQ(isAnyOf(res24, 
               VORONOI_OVER_VERTEX_3,
               VORONOI_OVER_EDGE_2_3,
               VORONOI_OVER_EDGE_3_4,
               VORONOI_OVER_TRIANGLE_2_3_4), true);

        auto res25 = g01.find(p01, p02, p03, p04, pt25);

        EXPECT_EQ(isAnyOf(res25, 
               VORONOI_OVER_VERTEX_4,
               VORONOI_OVER_EDGE_1_4,
               VORONOI_OVER_EDGE_2_4,
               VORONOI_OVER_TRIANGLE_1_2_4), true);

        auto res26 = g01.find(p01, p02, p03, p04, pt26);

        EXPECT_EQ(isAnyOf(res26, 
               VORONOI_OVER_VERTEX_4,
               VORONOI_OVER_EDGE_1_4,
               VORONOI_OVER_EDGE_3_4,
               VORONOI_OVER_TRIANGLE_1_4_3), true);

        auto res27 = g01.find(p01, p02, p03, p04, pt27);

        EXPECT_EQ(isAnyOf(res27, 
               VORONOI_OVER_VERTEX_4,
               VORONOI_OVER_EDGE_2_4,
               VORONOI_OVER_EDGE_3_4,
               VORONOI_OVER_TRIANGLE_2_3_4), true);

        auto res28 = g01.find(p01, p02, p03, p04, pt28);

        EXPECT_EQ(isAnyOf(res28, 
               VORONOI_OVER_EDGE_1_2,
               VORONOI_OVER_TRIANGLE_1_3_2), true);


        auto res29 = g01.find(p01, p02, p03, p04, pt29);

        EXPECT_EQ(isAnyOf(res29, 
               VORONOI_OVER_EDGE_1_2,
               VORONOI_OVER_TRIANGLE_1_2_4), true);

        auto res30 = g01.find(p01, p02, p03, p04, pt30);

        EXPECT_EQ(isAnyOf(res30, 
               VORONOI_OVER_EDGE_2_3,
               VORONOI_OVER_TRIANGLE_1_3_2), true);

        auto res31 = g01.find(p01, p02, p03, p04, pt31);

        EXPECT_EQ(isAnyOf(res31, 
               VORONOI_OVER_EDGE_2_3,
               VORONOI_OVER_TRIANGLE_2_3_4), true);

        auto res32 = g01.find(p01, p02, p03, p04, pt32);

        EXPECT_EQ(isAnyOf(res32, 
               VORONOI_OVER_EDGE_3_1,
               VORONOI_OVER_TRIANGLE_1_3_2), true);

        auto res33 = g01.find(p01, p02, p03, p04, pt33);

        EXPECT_EQ(isAnyOf(res33, 
               VORONOI_OVER_EDGE_3_1,
               VORONOI_OVER_TRIANGLE_1_4_3), true);

        auto res34 = g01.find(p01, p02, p03, p04, pt34);

        EXPECT_EQ(isAnyOf(res34, 
               VORONOI_OVER_EDGE_1_4,
               VORONOI_OVER_TRIANGLE_1_2_4), true);

        auto res35 = g01.find(p01, p02, p03, p04, pt35);

        EXPECT_EQ(isAnyOf(res35, 
               VORONOI_OVER_EDGE_1_4,
               VORONOI_OVER_TRIANGLE_1_4_3), true);

        auto res36 = g01.find(p01, p02, p03, p04, pt36);

        EXPECT_EQ(isAnyOf(res36, 
               VORONOI_OVER_EDGE_2_4,
               VORONOI_OVER_TRIANGLE_1_2_4), true);

        auto res37 = g01.find(p01, p02, p03, p04, pt37);

        EXPECT_EQ(isAnyOf(res37, 
               VORONOI_OVER_EDGE_2_4,
               VORONOI_OVER_TRIANGLE_2_3_4), true);

        auto res38 = g01.find(p01, p02, p03, p04, pt38);

        EXPECT_EQ(isAnyOf(res38, 
               VORONOI_OVER_EDGE_3_4,
               VORONOI_OVER_TRIANGLE_1_4_3), true);

        auto res39 = g01.find(p01, p02, p03, p04, pt39);

        EXPECT_EQ(isAnyOf(res39, 
               VORONOI_OVER_EDGE_3_4,
               VORONOI_OVER_TRIANGLE_2_3_4), true);
    }
}

/** @brief tests find() flat tetrahedron
 */
TEST_F(Voronoi3SimplexGraphTest, Test05) {

    Voronoi3SimplexGraph g01;

    for (long i = 0; i < 100000; i++) {

        Vec3 trans((rand()%100)/10.0, (rand()%100)/10.0, (rand()%100)/10.0);
        Mat3x3 M = randomRotMat();

        Vec3 p01( 1000.0,   0.0,   500.0);
        Vec3 p02(    0.0,   0.0,  -500.0);
        Vec3 p03(-1000.0,   0.0,   500.0);
        Vec3 p04(    0.0,   0.1,     0.0);

        Vec3 n132 = (p03 - p01).cross(p02 - p01);
        n132.normalize();

        Vec3 n124 = (p02 - p01).cross(p04 - p01);
        n124.normalize();

        Vec3 n143 = (p04 - p01).cross(p03 - p01);
        n143.normalize();

        Vec3 n234 = (p03 - p02).cross(p04 - p02);
        n234.normalize();

        double a01, a02, a03;

        Vec3 pt01 = p04 + Vec3(0.0,    0.01, 0.0);
        Vec3 pt02 = p04 + Vec3(0.0,10000.0, 0.0);
        Vec3 pt03 = p04 + n143 * rand100();
        generateRandomAlphas(a01, a02);
        Vec3 pt04_1 = p04 + n143 * rand100();
        Vec3 pt04_2 = p04 + Vec3(0.0, 1.0, 0.0) * rand100();
        Vec3 pt04   = pt04_1 * a01 + pt04_2 * a02;
        Vec3 pt05 = p04 + n143 * rand100() + Vec3(0.0, 0.0, +0.000001);

        Vec3 pt06_1 = p04 + n143 * rand100();
        Vec3 pt06_2 = p04 + n124 * rand100();
        Vec3 pt06_3 = p04 + n234 * rand100();


        generateRandomAlphas(a01, a02, a03);
        Vec3 pt06 = pt06_1 * a01 +  pt06_2 * a02 +  pt06_3 * a03;

        generateRandomAlphas(a01, a02);
        Vec3 pt07 = pt06_1 * a01 +  pt06_2 * a02;
        Vec3 pt08 = pt06_2 * a01 +  pt06_3 * a02;
        Vec3 pt09 = pt06_3 * a01 +  pt06_1 * a02;

        generateRandomAlphas(a01, a02);
        Vec3 pt10_1 = p01 * a01 +  p04 * a02;
        Vec3 pt10_2 = pt10_1 + n143 * rand100();
        Vec3 pt10_3 = pt10_1 + n124 * rand100();
        generateRandomAlphas(a01, a02);
        Vec3 pt10 = pt10_2 * a01 +  pt10_3 * a02;
        Vec3 pt11 = pt10_2;
        Vec3 pt12 = pt10_3;

        p01  = M * p01  + trans;
        p02  = M * p02  + trans;
        p03  = M * p03  + trans;
        p04  = M * p04  + trans;
        pt01 = M * pt01 + trans;
        pt02 = M * pt02 + trans;
        pt03 = M * pt03 + trans;
        pt04 = M * pt04 + trans;
        pt05 = M * pt05 + trans;
        pt06 = M * pt06 + trans;
        pt07 = M * pt07 + trans;
        pt08 = M * pt08 + trans;
        pt09 = M * pt09 + trans;
        pt10 = M * pt10 + trans;
        pt11 = M * pt11 + trans;
        pt12 = M * pt12 + trans;


        auto res01 = g01.find(p01, p02, p03, p04, pt01);
        EXPECT_EQ(res01, VORONOI_OVER_VERTEX_4);

        auto res02 = g01.find(p01, p02, p03, p04, pt02);
        EXPECT_EQ(res02, VORONOI_OVER_VERTEX_4);

        auto res03 = g01.find(p01, p02, p03, p04, pt03);
        EXPECT_EQ(isAnyOf(res03, VORONOI_OVER_VERTEX_4,
                                 VORONOI_OVER_EDGE_1_4,
                                 VORONOI_OVER_EDGE_3_4,
                                 VORONOI_OVER_TRIANGLE_1_4_3), true);
                
        auto res04 = g01.find(p01, p02, p03, p04, pt04);
        EXPECT_EQ(res04, VORONOI_OVER_VERTEX_4);

        auto res05 = g01.find(p01, p02, p03, p04, pt05);
        EXPECT_EQ(res05, VORONOI_OVER_TRIANGLE_1_4_3);

        auto res06 = g01.find(p01, p02, p03, p04, pt06);
        EXPECT_EQ(res06, VORONOI_OVER_VERTEX_4);

        auto res07 = g01.find(p01, p02, p03, p04, pt07);
        EXPECT_EQ(isAnyOf(res07, VORONOI_OVER_VERTEX_4,
                                 VORONOI_OVER_EDGE_1_4), true);
                
        auto res08 = g01.find(p01, p02, p03, p04, pt08);
        EXPECT_EQ(isAnyOf(res08, VORONOI_OVER_VERTEX_4,
                                 VORONOI_OVER_EDGE_2_4), true);

        auto res09 = g01.find(p01, p02, p03, p04, pt09);
        EXPECT_EQ(isAnyOf(res09, VORONOI_OVER_VERTEX_4,
                                 VORONOI_OVER_EDGE_3_4), true);

//cerr << "p01: " << p01 << "\n";
//cerr << "p02: " << p02 << "\n";
//cerr << "p03: " << p03 << "\n";
//cerr << "p04: " << p04 << "\n";
//cerr << "pt10_1: " << pt10_1 << "\n";

        auto res10 = g01.find(p01, p02, p03, p04, pt10);
        EXPECT_EQ(res10, VORONOI_OVER_EDGE_1_4);

        auto res11 = g01.find(p01, p02, p03, p04, pt11);
        EXPECT_EQ(isAnyOf(res11, VORONOI_OVER_EDGE_1_4,
                                 VORONOI_OVER_TRIANGLE_1_4_3),true);
        auto res12 = g01.find(p01, p02, p03, p04, pt12);

        EXPECT_EQ(isAnyOf(res12, VORONOI_OVER_EDGE_1_4,
                                 VORONOI_OVER_TRIANGLE_1_2_4),true);

    }
}

/** @brief tests find() needle tetrahedron
 */
TEST_F(Voronoi3SimplexGraphTest, Test06) {

    Voronoi3SimplexGraph g01;

    for (long i = 0; i < 10000; i++) {


        Vec3 trans((rand()%100)/10.0, (rand()%100)/10.0, (rand()%100)/10.0);
        Mat3x3 M = randomRotMat();

        Vec3 p01(    0.1,   0.0,    0.05);
        Vec3 p02(    0.0,   0.0,   -0.05);
        Vec3 p03(   -0.1,   0.0,    0.05);
        Vec3 p04(    0.0, 1000.0,     0.0);

        Vec3 n132 = (p03 - p01).cross(p02 - p01);
        n132.normalize();

        Vec3 n124 = (p02 - p01).cross(p04 - p01);
        n124.normalize();

        Vec3 n143 = (p04 - p01).cross(p03 - p01);
        n143.normalize();

        Vec3 n234 = (p03 - p02).cross(p04 - p02);
        n234.normalize();

        double a01, a02, a03;

        Vec3 pt01 = p04 + Vec3(0.0,    0.01, 0.0);
        Vec3 pt02 = p04 + Vec3(0.0,10000.0, 0.0);
        Vec3 pt03 = p04 + n143 * rand100();


        // Not used as this exeeds numerical precision.
        generateRandomAlphas(a01, a02);
        Vec3 pt04_1 = p04 + n143 * rand100();
        Vec3 pt04_2 = p04 + Vec3(0.0, 1.0, 0.0) * rand100();
        Vec3 pt04   = pt04_1 * a01 + pt04_2 * a02;
        Vec3 pt05 = p04 + n143 * rand100() + Vec3(0.0, 0.0, 1.0);

        Vec3 pt06_1 = p04 + n143 * rand100();
        Vec3 pt06_2 = p04 + n124 * rand100();
        Vec3 pt06_3 = p04 + n234 * rand100();


        generateRandomAlphas(a01, a02, a03);
        Vec3 pt06 = pt06_1 * a01 +  pt06_2 * a02 +  pt06_3 * a03;

        generateRandomAlphas(a01, a02);
        Vec3 pt07 = pt06_1 * a01 +  pt06_2 * a02;
        Vec3 pt08 = pt06_2 * a01 +  pt06_3 * a02;
        Vec3 pt09 = pt06_3 * a01 +  pt06_1 * a02;

        generateRandomAlphas(a01, a02);
        Vec3 pt10_1 = p01 * a01 +  p04 * a02;
        Vec3 pt10_2 = pt10_1 + n143 * rand100();
        Vec3 pt10_3 = pt10_1 + n124 * rand100();
        generateRandomAlphas(a01, a02);
        Vec3 pt10 = pt10_2 * a01 +  pt10_3 * a02;
        Vec3 pt11 = pt10_2;
        Vec3 pt12 = pt10_3;

        p01  = M * p01  + trans;
        p02  = M * p02  + trans;
        p03  = M * p03  + trans;
        p04  = M * p04  + trans;
        pt01 = M * pt01 + trans;
        pt02 = M * pt02 + trans;
        pt03 = M * pt03 + trans;
        pt04 = M * pt04 + trans;
        pt05 = M * pt05 + trans;
        pt06 = M * pt06 + trans;
        pt07 = M * pt07 + trans;
        pt08 = M * pt08 + trans;
        pt09 = M * pt09 + trans;
        pt10 = M * pt10 + trans;
        pt11 = M * pt11 + trans;
        pt12 = M * pt12 + trans;


        auto res01 = g01.find(p01, p02, p03, p04, pt01);
        EXPECT_EQ(res01, VORONOI_OVER_VERTEX_4);

        auto res02 = g01.find(p01, p02, p03, p04, pt02);
        EXPECT_EQ(res02, VORONOI_OVER_VERTEX_4);

        auto res03 = g01.find(p01, p02, p03, p04, pt03);
        EXPECT_EQ(isAnyOf(res03, VORONOI_OVER_VERTEX_4,
                                 VORONOI_OVER_EDGE_1_4,
                                 VORONOI_OVER_EDGE_3_4,
                                 VORONOI_OVER_TRIANGLE_1_4_3), true);
                
        auto res04 = g01.find(p01, p02, p03, p04, pt04);
        EXPECT_EQ(res04, VORONOI_OVER_VERTEX_4);

        auto res05 = g01.find(p01, p02, p03, p04, pt05);

        EXPECT_EQ(res05, VORONOI_OVER_TRIANGLE_1_4_3);

        auto res06 = g01.find(p01, p02, p03, p04, pt06);
        EXPECT_EQ(res06, VORONOI_OVER_VERTEX_4);

        auto res07 = g01.find(p01, p02, p03, p04, pt07);
        EXPECT_EQ(isAnyOf(res07, VORONOI_OVER_VERTEX_4,
                                 VORONOI_OVER_EDGE_1_4), true);
                
        auto res08 = g01.find(p01, p02, p03, p04, pt08);
        EXPECT_EQ(isAnyOf(res08, VORONOI_OVER_VERTEX_4,
                                 VORONOI_OVER_EDGE_2_4), true);

        auto res09 = g01.find(p01, p02, p03, p04, pt09);
        EXPECT_EQ(isAnyOf(res09, VORONOI_OVER_VERTEX_4,
                                 VORONOI_OVER_EDGE_3_4), true);

//cerr << "p01: " << p01 << "\n";
//cerr << "p02: " << p02 << "\n";
//cerr << "p03: " << p03 << "\n";
//cerr << "p04: " << p04 << "\n";
//cerr << "pt10_1: " << pt10_1 << "\n";

        auto res10 = g01.find(p01, p02, p03, p04, pt10);
        EXPECT_EQ(res10, VORONOI_OVER_EDGE_1_4);

        auto res11 = g01.find(p01, p02, p03, p04, pt11);
        EXPECT_EQ(isAnyOf(res11, VORONOI_OVER_EDGE_1_4,
                                 VORONOI_OVER_TRIANGLE_1_4_3),true);
        auto res12 = g01.find(p01, p02, p03, p04, pt12);

        EXPECT_EQ(isAnyOf(res12, VORONOI_OVER_EDGE_1_4,
                                 VORONOI_OVER_TRIANGLE_1_2_4),true);

    }
}












} // namespace Makena

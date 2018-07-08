#include "gtest/gtest.h"
#include <vector>
#include "primitives.hpp"
#include "manifold.hpp"
#include "quaternion.hpp"
#include "bd_boundary_simplex_finder.hpp"

#include <iostream>
#include <fstream>
#include <sstream>


namespace Makena {


class BDBoundarySimplexFinderTest : public ::testing::Test {

  protected:  

    BDBoundarySimplexFinderTest(){;};
    virtual ~BDBoundarySimplexFinderTest(){;};
    virtual void SetUp() {;};
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


static Mat3x3 randomRotMatXY()
{
    Vec3 v(0.0, 0.0, 0.0);
    while (v == Vec3(0.0, 0.0, 0.0)) {
        v = Vec3( 0.0,
                  0.0,
                  (rand()%100)/100.0 - 0.5 );
    }
    
    Quaternion q(v, (rand()%100)/50.0 * M_PI);
    q.normalize();
    return q.rotationMatrix();
}

/*
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
*/
/*
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
*/
/*
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
*/

static double rand100()
{
    return ((rand()%10000)+1)/100.0;
}

/*
static bool isAnyOf(
    enum predicate pTest,
    enum predicate p1,
    enum predicate p2,
    enum predicate p3,
    enum predicate p4
) {
    return pTest==p1 || pTest==p2 || pTest==p3 ||pTest==p4;
}
*/
/*
static bool isAnyOf(
    enum predicate pTest,
    enum predicate p1,
    enum predicate p2,
    enum predicate p3
) {
    return pTest==p1 || pTest==p2 || pTest==p3;
}
*/
/*
static bool isAnyOf(
    enum predicate pTest,
    enum predicate p1,
    enum predicate p2
) {
    return pTest==p1 || pTest==p2;
}
*/
/*
static bool isAnyOf(
    enum predicate pTest,
    enum predicate p1
) {
    return pTest==p1;
}
*/

/**  @brief isOriginInsideTriangleXY()
 */
TEST_F(BDBoundarySimplexFinderTest, Test01) {


    for (long i = 0; i < 1000; i++) {

        auto M = randomRotMatXY();

        Vec3 p01( -1.0, -1.0, rand100());
        Vec3 p02(  5.0, -1.0, rand100());
        Vec3 p03( -1.0,  5.0, rand100());
        Vec3 p04(  1.0,  1.0, rand100());
        Vec3 p05( -1.0,  4.0, rand100());
        Vec3 p06(  1.0, -4.0, rand100());
        Vec3 p07( 16.0,-10.0, rand100());
        Vec3 p08(  0.0,  0.0, rand100());

        p01 = M * p01;
        p02 = M * p02;
        p03 = M * p03;
        p04 = M * p04;
        p05 = M * p05;
        p06 = M * p06;
        p07 = M * p07;
        p08 = M * p08;

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);
        vector<BDVertex> Q;
        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                  std::cerr);

        EXPECT_EQ(f.isOriginInsideTriangleXY(p01, p02, p03), true);
        EXPECT_EQ(f.isOriginInsideTriangleXY(p02, p03, p01), true);
        EXPECT_EQ(f.isOriginInsideTriangleXY(p03, p01, p02), true);

        EXPECT_EQ(f.isOriginInsideTriangleXY(p04, p02, p03), false);
        EXPECT_EQ(f.isOriginInsideTriangleXY(p02, p03, p04), false);
        EXPECT_EQ(f.isOriginInsideTriangleXY(p03, p04, p02), false);

        EXPECT_EQ(f.isOriginInsideTriangleXY(p05, p06, p07), true);
        EXPECT_EQ(f.isOriginInsideTriangleXY(p06, p07, p05), true);
        EXPECT_EQ(f.isOriginInsideTriangleXY(p07, p05, p06), true);

        EXPECT_EQ(f.isOriginInsideTriangleXY(p08, p06, p07), true);
        EXPECT_EQ(f.isOriginInsideTriangleXY(p06, p07, p08), true);
        EXPECT_EQ(f.isOriginInsideTriangleXY(p07, p08, p06), true);

    }
}


/**  @brief isOriginInsideTriangleXYrelaxed()
 */
TEST_F(BDBoundarySimplexFinderTest, Test02) {


    for (long i = 0; i < 1000; i++) {

        auto M = randomRotMatXY();

        Vec3 p01( -1.0, -1.0, rand100());
        Vec3 p02(  5.0, -1.0, rand100());
        Vec3 p03( -1.0,  5.0, rand100());
        Vec3 p04(  1.0,  1.0, rand100());
        Vec3 p05( -1.0,  4.0, rand100());
        Vec3 p06(  1.0, -4.0, rand100());
        Vec3 p07( 16.0,-10.0, rand100());
        Vec3 p08(  0.0,  0.0, rand100());
        Vec3 p09( 1.0e-11, -1.0e-11, rand100());
        Vec3 p10( 1.0e-9, -1.0e-9, rand100());
        Vec3 p11( -1.0 + 1.0e-11,  4.0, rand100());
        Vec3 p12(  1.0 + 1.0e-11, -4.0, rand100());
        Vec3 p13( -1.0 + 1.0e-9,  4.0, rand100());
        Vec3 p14(  1.0 + 1.0e-9, -4.0, rand100());
        Vec3 p15( -1.0 - 1.0e-9,  4.0, rand100());
        Vec3 p16(  1.0 - 1.0e-9, -4.0, rand100());

        p01 = M * p01;
        p02 = M * p02;
        p03 = M * p03;
        p04 = M * p04;
        p05 = M * p05;
        p06 = M * p06;
        p07 = M * p07;
        p08 = M * p08;
        p09 = M * p09;
        p10 = M * p10;
        p11 = M * p11;
        p12 = M * p12;
        p13 = M * p13;
        p14 = M * p14;
        p15 = M * p15;
        p16 = M * p16;

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);
        vector<BDVertex> Q;

        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);

        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p01, p02, p03), true);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p02, p03, p01), true);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p03, p01, p02), true);

        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p04, p02, p03), false);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p02, p03, p04), false);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p03, p04, p02), false);

        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p05, p06, p07), true);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p06, p07, p05), true);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p07, p05, p06), true);

        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p08, p06, p07), true);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p06, p07, p08), true);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p07, p08, p06), true);

        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p06, p07, p09), true);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p07, p09, p06), true);

        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p06, p07, p10), false);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p07, p10, p06), false);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p10, p06, p07), false);

        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p11, p12, p07), true);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p12, p07, p11), true);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p07, p11, p12), true);

        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p13, p14, p07), false);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p14, p07, p13), false);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p07, p13, p14), false);

        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p15, p16, p07), true);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p16, p07, p15), true);
        EXPECT_EQ(f.isOriginInsideTriangleXYrelaxed(p07, p15, p16), true);

    }
}


/**  @brief findInitial2SimplexFrom3Simplex() #1
 */
TEST_F(BDBoundarySimplexFinderTest, Test03) {


    for (long i = 0; i < 1000; i++) {

        auto M = randomRotMatXY();

        Vec3 p01( -5.0, -5.0, -1.0);
        Vec3 p02(  5.0, -5.0, -1.0);
        Vec3 p03(  0.0,  5.0, -1.0);
        Vec3 p04(  1.0,  2.0,  5.0);

        p01 = M * p01;
        p02 = M * p02;
        p03 = M * p03;
        p04 = M * p04;

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);
        vector<BDVertex> Q;
        VertexIt vit01;

        BDVertex mv01(vit01, vit01);
        mv01.mP = p01;
        BDVertex mv02(vit01, vit01);
        mv02.mP = p02;
        BDVertex mv03(vit01, vit01);
        mv03.mP = p03;
        BDVertex mv04(vit01, vit01);
        mv04.mP = p04;

        Q.push_back(mv01);
        Q.push_back(mv02);
        Q.push_back(mv03);
        Q.push_back(mv04);

        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);

        EXPECT_EQ(f.findInitial2SimplexFrom3Simplex(), true);
        EXPECT_EQ(Q.size(),3);
        EXPECT_EQ(Q[0].mP, p01);
        EXPECT_EQ(Q[1].mP, p02);
        EXPECT_EQ(Q[2].mP, p04);
    }
}


/**  @brief findInitial2SimplexFrom3Simplex() #2
 */
TEST_F(BDBoundarySimplexFinderTest, Test04) {


    for (long i = 0; i < 1000; i++) {

        auto M = randomRotMatXY();

        Vec3 p01( -5.0, -5.0, -1.0);
        Vec3 p02(  5.0, -5.0, -1.0);
        Vec3 p03(  0.0,  5.0, -1.0);
        Vec3 p04(  2.0,  1.0,  5.0);

        p01 = M * p01;
        p02 = M * p02;
        p03 = M * p03;
        p04 = M * p04;

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);
        vector<BDVertex> Q;
        VertexIt vit01;

        BDVertex mv01(vit01, vit01);
        mv01.mP = p01;
        BDVertex mv02(vit01, vit01);
        mv02.mP = p02;
        BDVertex mv03(vit01, vit01);
        mv03.mP = p03;
        BDVertex mv04(vit01, vit01);
        mv04.mP = p04;

        Q.push_back(mv01);
        Q.push_back(mv02);
        Q.push_back(mv03);
        Q.push_back(mv04);

        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);

        EXPECT_EQ(f.findInitial2SimplexFrom3Simplex(), true);
        EXPECT_EQ(Q.size(),3);
        EXPECT_EQ(Q[0].mP, p01);
        EXPECT_EQ(Q[1].mP, p04);
        EXPECT_EQ(Q[2].mP, p03);
    }
}


/**  @brief findInitial2SimplexFrom3Simplex() #3
 */
TEST_F(BDBoundarySimplexFinderTest, Test05) {


    for (long i = 0; i < 1000; i++) {

        auto M = randomRotMatXY();

        Vec3 p01( -5.0, -5.0, -1.0);
        Vec3 p02(  5.0, -5.0, -1.0);
        Vec3 p03(  0.0,  5.0, -1.0);
        Vec3 p04( -2.0,  1.0,  5.0);

        p01 = M * p01;
        p02 = M * p02;
        p03 = M * p03;
        p04 = M * p04;

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);
        vector<BDVertex> Q;
        VertexIt vit01;

        BDVertex mv01(vit01, vit01);
        mv01.mP = p01;
        BDVertex mv02(vit01, vit01);
        mv02.mP = p02;
        BDVertex mv03(vit01, vit01);
        mv03.mP = p03;
        BDVertex mv04(vit01, vit01);
        mv04.mP = p04;

        Q.push_back(mv01);
        Q.push_back(mv02);
        Q.push_back(mv03);
        Q.push_back(mv04);

        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);

        EXPECT_EQ(f.findInitial2SimplexFrom3Simplex(), true);
        EXPECT_EQ(Q.size(),3);
        EXPECT_EQ(Q[0].mP, p02);
        EXPECT_EQ(Q[1].mP, p03);
        EXPECT_EQ(Q[2].mP, p04);
    }
}


/**  @brief findInitial2SimplexFrom3Simplex() #4
 */
TEST_F(BDBoundarySimplexFinderTest, Test06) {


    for (long i = 0; i < 1000; i++) {

        auto M = randomRotMatXY();

        Vec3 p01( -5.0, -5.0, -1.0);
        Vec3 p02(  5.0, -5.0, -1.0);
        Vec3 p03(  0.0,  5.0, -1.0);
        Vec3 p04(  1.0,  1.0,  5.0);

        p01 = M * p01;
        p02 = M * p02;
        p03 = M * p03;
        p04 = M * p04;

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);
        vector<BDVertex> Q;
        VertexIt vit01;

        BDVertex mv01(vit01, vit01);
        mv01.mP = p01;
        BDVertex mv02(vit01, vit01);
        mv02.mP = p02;
        BDVertex mv03(vit01, vit01);
        mv03.mP = p03;
        BDVertex mv04(vit01, vit01);
        mv04.mP = p04;

        Q.push_back(mv01);
        Q.push_back(mv02);
        Q.push_back(mv03);
        Q.push_back(mv04);

        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                   std::cerr);

        EXPECT_EQ(f.findInitial2SimplexFrom3Simplex(), true);
        EXPECT_EQ(Q.size(),3);
        EXPECT_EQ( ( Q[0].mP==p01 && Q[1].mP==p02 && Q[2].mP==p04 )||
                   ( Q[0].mP==p01 && Q[1].mP==p04 && Q[2].mP==p03 )  , true);
    }
}


/**  @brief findInitial2SimplexFrom3Simplex() #5
 */
TEST_F(BDBoundarySimplexFinderTest, Test07) {


    for (long i = 0; i < 1000; i++) {

        auto M = randomRotMatXY();

        Vec3 p01( -5.0, -5.0, -1.0);
        Vec3 p02(  5.0, -5.0, -1.0);
        Vec3 p03(  0.0,  5.0, -1.0);
        Vec3 p04( -1.0,  1.0,  5.0);

        p01 = M * p01;
        p02 = M * p02;
        p03 = M * p03;
        p04 = M * p04;

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);
        vector<BDVertex> Q;
        VertexIt vit01;

        BDVertex mv01(vit01, vit01);
        mv01.mP = p01;
        BDVertex mv02(vit01, vit01);
        mv02.mP = p02;
        BDVertex mv03(vit01, vit01);
        mv03.mP = p03;
        BDVertex mv04(vit01, vit01);
        mv04.mP = p04;

        Q.push_back(mv01);
        Q.push_back(mv02);
        Q.push_back(mv03);
        Q.push_back(mv04);

        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);

        EXPECT_EQ(f.findInitial2SimplexFrom3Simplex(), true);
        EXPECT_EQ(Q.size(),3);
        EXPECT_EQ( ( Q[0].mP==p01 && Q[1].mP==p02 && Q[2].mP==p04 )||
                   ( Q[0].mP==p02 && Q[1].mP==p03 && Q[2].mP==p04 )  , true);
    }
}


/**  @brief findInitial2SimplexFrom3Simplex() #6
 */
TEST_F(BDBoundarySimplexFinderTest, Test08) {


    for (long i = 0; i < 1000; i++) {

        auto M = randomRotMatXY();

        Vec3 p01( -5.0, -5.0, -1.0);
        Vec3 p02(  5.0, -5.0, -1.0);
        Vec3 p03(  0.0,  5.0, -1.0);
        Vec3 p04(  0.0, -1.0,  5.0);

        p01 = M * p01;
        p02 = M * p02;
        p03 = M * p03;
        p04 = M * p04;

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);
        vector<BDVertex> Q;
        VertexIt vit01;

        BDVertex mv01(vit01, vit01);
        mv01.mP = p01;
        BDVertex mv02(vit01, vit01);
        mv02.mP = p02;
        BDVertex mv03(vit01, vit01);
        mv03.mP = p03;
        BDVertex mv04(vit01, vit01);
        mv04.mP = p04;

        Q.push_back(mv01);
        Q.push_back(mv02);
        Q.push_back(mv03);
        Q.push_back(mv04);

        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);

        EXPECT_EQ(f.findInitial2SimplexFrom3Simplex(), true);
        EXPECT_EQ(Q.size(),3);
        EXPECT_EQ( ( Q[0].mP==p01 && Q[1].mP==p04 && Q[2].mP==p03 )||
                   ( Q[0].mP==p02 && Q[1].mP==p03 && Q[2].mP==p04 )  , true);
    }
}


/**  @brief findInitial2SimplexFrom3Simplex() #7
 */
TEST_F(BDBoundarySimplexFinderTest, Test09) {

    for (long i = 0; i < 1000; i++) {

        auto M = randomRotMatXY();

        Vec3 p01( -5.0, -5.0, -1.0);
        Vec3 p02(  5.0, -5.0, -1.0);
        Vec3 p03(  0.0,  5.0, -1.0);
        Vec3 p04(  0.0,  0.0,  5.0);

        p01 = M * p01;
        p02 = M * p02;
        p03 = M * p03;
        p04 = M * p04;

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);
        vector<BDVertex> Q;
        VertexIt vit01;

        BDVertex mv01(vit01, vit01);
        mv01.mP = p01;
        BDVertex mv02(vit01, vit01);
        mv02.mP = p02;
        BDVertex mv03(vit01, vit01);
        mv03.mP = p03;
        BDVertex mv04(vit01, vit01);
        mv04.mP = p04;

        Q.push_back(mv01);
        Q.push_back(mv02);
        Q.push_back(mv03);
        Q.push_back(mv04);

        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);

        EXPECT_EQ(f.findInitial2SimplexFrom3Simplex(), true);
        EXPECT_EQ(Q.size(),3);
        EXPECT_EQ( ( Q[0].mP==p01 && Q[1].mP==p02 && Q[2].mP==p04 )||
                   ( Q[0].mP==p01 && Q[1].mP==p04 && Q[2].mP==p03 )||
                   ( Q[0].mP==p02 && Q[1].mP==p03 && Q[2].mP==p04 )  , true);
    }
}


/**  @brief findInitial2SimplexFrom3Simplex() #8
 */
TEST_F(BDBoundarySimplexFinderTest, Test10) {

    for (long i = 0; i < 1000; i++) {

        auto M = randomRotMatXY();

        Vec3 p01(  5.0, -5.0,  1.0);
        Vec3 p02( -5.0, -5.0,  1.0);
        Vec3 p03(  0.0,  5.0,  1.0);
        Vec3 p04( rand100()-50.0, rand100()-50.0, -5.0);

        p01 = M * p01;
        p02 = M * p02;
        p03 = M * p03;
        p04 = M * p04;

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);
        vector<BDVertex> Q;
        VertexIt vit01;

        BDVertex mv01(vit01, vit01);
        mv01.mP = p01;
        BDVertex mv02(vit01, vit01);
        mv02.mP = p02;
        BDVertex mv03(vit01, vit01);
        mv03.mP = p03;
        BDVertex mv04(vit01, vit01);
        mv04.mP = p04;

        Q.push_back(mv01);
        Q.push_back(mv02);
        Q.push_back(mv03);
        Q.push_back(mv04);

        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);

        EXPECT_EQ(f.findInitial2SimplexFrom3Simplex(), true);
        EXPECT_EQ(Q.size(),3);
        EXPECT_EQ( (Q[0].mP==p01 && Q[1].mP==p03 && Q[2].mP==p02), true);
    }
}


/**  @brief findInitial2SimplexFrom1Simplex()
 */
TEST_F(BDBoundarySimplexFinderTest, Test11) {

    for (long i = 0; i < 1000; i++) {

        Manifold m01, m02;
        m01.constructCuboid(
            Vec3(  1.0,  -1.0,  1.0),
            Vec3(  1.0,   1.0,  1.0),
            Vec3( -1.0,   1.0,  1.0),
            Vec3( -1.0,  -1.0,  1.0),
            Vec3(  1.0,  -1.0, -1.0),
            Vec3(  1.0,   1.0, -1.0),
            Vec3( -1.0,   1.0, -1.0),
            Vec3( -1.0,  -1.0, -1.0) );

        m02.constructCuboid(
            Vec3(  1.0,  -1.0,  1.0),
            Vec3(  1.0,   1.0,  1.0),
            Vec3( -1.0,   1.0,  1.0),
            Vec3( -1.0,  -1.0,  1.0),
            Vec3(  1.0,  -1.0, -1.0),
            Vec3(  1.0,   1.0, -1.0),
            Vec3( -1.0,   1.0, -1.0),
            Vec3( -1.0,  -1.0, -1.0) );

        Manifold::Martialled mm01 = m01.exportData();
        Manifold::Martialled mm02 = m02.exportData();
       
        VertexIt vit01;
        for (vit01 = m01.vertices().first;
             vit01 != m01.vertices().second && 
             (*vit01)->pLCS()!=Vec3(1.0, 1.0, 1.0);
             vit01++) {;}

        VertexIt vit02;
        for (vit02 = m01.vertices().first;
             vit02 != m01.vertices().second && 
             (*vit02)->pLCS()!=Vec3(-1.0, -1.0, -1.0);
             vit02++) {;}

        VertexIt vit03;
        for (vit03 = m02.vertices().first;
             vit03 != m02.vertices().second && 
             (*vit03)->pLCS()!=Vec3(-1.0, -1.0, -1.0);
             vit03++) {;}

        VertexIt vit04;
        for (vit04 = m02.vertices().first;
             vit04 != m02.vertices().second && 
             (*vit04)->pLCS()!=Vec3( 1.0,  1.0,  1.0);
             vit04++) {;}

        Vec3 v_ref01( 2.0, -2.0, -2.0);
        Vec3 v_ref02(-2.0,  2.0, -2.0);
        Vec3 v_ref03(-2.0, -2.0,  2.0);
        Vec3 v_ref04(-2.0,  2.0,  2.0);
        Vec3 v_ref05( 2.0, -2.0,  2.0);
        Vec3 v_ref06( 2.0,  2.0, -2.0);

        Vec3 v_ref07( 2.0,  2.0,  2.0);
        Vec3 v_ref08(-2.0, -2.0, -2.0);

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);

        Mat3x3     U = randomRotMat();
        Vec3       zv(rand100()-50.0,rand100()-50.0,rand100()-50.0);
        Quaternion zq(U);

        v_ref01 = U * v_ref01;
        v_ref02 = U * v_ref02;
        v_ref03 = U * v_ref03;
        v_ref04 = U * v_ref04;
        v_ref05 = U * v_ref05;
        v_ref06 = U * v_ref06;
        v_ref07 = U * v_ref07;
        v_ref08 = U * v_ref08;

        b01.setObjectAttributes(mm01, 1.0, U, 0.0, 0.0, false);
        b02.setObjectAttributes(mm02, 1.0, U, 0.0, 0.0, false);
        b01.setGeomConfig(zv, zq, zv, zv);
        b02.setGeomConfig(zv, zq, zv, zv);
        b01.updateGeomConfig(0.0);
        b02.updateGeomConfig(0.0);
        b01.commitGeomConfig();
        b02.commitGeomConfig();

        vector<BDVertex> Q;

        BDVertex mv01(vit01, vit03);
        mv01.updateP(m01, U, zv, m02, U, zv);
        BDVertex mv02(vit02, vit04);
        mv02.updateP(m01, U, zv, m02, U, zv);
        Q.push_back(mv01);
        Q.push_back(mv02);

        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                   std::cerr);

        f.mRotMat1 = U;
        f.mRotMat2 = U;
        f.mCoM1    = zv;
        f.mCoM2    = zv;

        f.findInitial2SimplexFrom1Simplex();
        EXPECT_EQ(Q.size(), 3);
        EXPECT_EQ(Q[0].v1(), vit01);
        EXPECT_EQ(Q[0].v2(), vit03);
        EXPECT_EQ(Q[0].p(),  v_ref07);
        EXPECT_EQ(Q[1].v1(), vit02);
        EXPECT_EQ(Q[1].v2(), vit04);
        EXPECT_EQ(Q[1].p(),  v_ref08);

        EXPECT_EQ( (Q[0].v1()!=Q[2].v1() || Q[0].v2()!=Q[2].v2()) &&
                   (Q[1].v1()!=Q[2].v1() || Q[1].v2()!=Q[2].v2())   , true);

        EXPECT_EQ( Q[2].p() == v_ref01 ||
                   Q[2].p() == v_ref02 ||
                   Q[2].p() == v_ref03 ||
                   Q[2].p() == v_ref04 ||
                   Q[2].p() == v_ref05 ||
                   Q[2].p() == v_ref06    , true);
    }
}



/**  @brief findInitial2SimplexFrom0Simplex()
 */
TEST_F(BDBoundarySimplexFinderTest, Test12) {

    for (long i = 0; i < 1000; i++) {

        Manifold m01, m02;
        m01.constructCuboid(
            Vec3(  1.0,  -1.0,  1.0),
            Vec3(  1.0,   1.0,  1.0),
            Vec3( -1.0,   1.0,  1.0),
            Vec3( -1.0,  -1.0,  1.0),
            Vec3(  1.0,  -1.0, -1.0),
            Vec3(  1.0,   1.0, -1.0),
            Vec3( -1.0,   1.0, -1.0),
            Vec3( -1.0,  -1.0, -1.0) );

        m02.constructCuboid(
            Vec3(  1.0,  -1.0,  1.0),
            Vec3(  1.0,   1.0,  1.0),
            Vec3( -1.0,   1.0,  1.0),
            Vec3( -1.0,  -1.0,  1.0),
            Vec3(  1.0,  -1.0, -1.0),
            Vec3(  1.0,   1.0, -1.0),
            Vec3( -1.0,   1.0, -1.0),
            Vec3( -1.0,  -1.0, -1.0) );

        Manifold::Martialled mm01 = m01.exportData();
        Manifold::Martialled mm02 = m02.exportData();
       
        VertexIt vit01;
        for (vit01 = m01.vertices().first;
             vit01 != m01.vertices().second && 
             (*vit01)->pLCS()!=Vec3(1.0, 1.0, 1.0);
             vit01++) {;}

        VertexIt vit02;
        for (vit02 = m02.vertices().first;
             vit02 != m02.vertices().second && 
             (*vit02)->pLCS()!=Vec3(1.0, 1.0, 1.0);
             vit02++) {;}

        Vec3 v_ref01( 0.0,  0.0,  0.0);
        Vec3 v_ref02( 2.0,  2.0,  2.0);
        Vec3 v_ref03( 2.0,  2.0, -2.0);
        Vec3 v_ref04( 2.0, -2.0,  2.0);
        Vec3 v_ref05( 2.0, -2.0, -2.0);
        Vec3 v_ref06(-2.0,  2.0,  2.0);
        Vec3 v_ref07(-2.0,  2.0, -2.0);
        Vec3 v_ref08(-2.0, -2.0,  2.0);
        Vec3 v_ref09(-2.0, -2.0, -2.0);

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);

        Mat3x3     U = randomRotMat();
        Vec3       zv(rand100()-50.0,rand100()-50.0,rand100()-50.0);
//        Mat3x3     U(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
//        Vec3       zv(0.0, 0.0, 0.0);

        Quaternion zq(U);

        v_ref01 = U * v_ref01;
        v_ref02 = U * v_ref02;
        v_ref03 = U * v_ref03;
        v_ref04 = U * v_ref04;
        v_ref05 = U * v_ref05;
        v_ref06 = U * v_ref06;
        v_ref07 = U * v_ref07;
        v_ref08 = U * v_ref08;
        v_ref09 = U * v_ref09;

        b01.setObjectAttributes(mm01, 1.0, U, 0.0, 0.0, false);
        b02.setObjectAttributes(mm02, 1.0, U, 0.0, 0.0, false);
        b01.setGeomConfig(zv, zq, zv, zv);
        b02.setGeomConfig(zv, zq, zv, zv);
        b01.updateGeomConfig(0.0);
        b02.updateGeomConfig(0.0);
        b01.commitGeomConfig();
        b02.commitGeomConfig();

        vector<BDVertex> Q;

        BDVertex mv01(vit01, vit02);
        mv01.updateP(m01, U, zv, m02, U, zv);
        Q.push_back(mv01);

        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                   std::cerr);

        f.mRotMat1 = U;
        f.mRotMat2 = U;
        f.mCoM1    = zv;
        f.mCoM2    = zv;
        f.findInitial2SimplexFrom0Simplex();
        EXPECT_EQ(Q.size(), 3);
        EXPECT_EQ(Q[0].v1(), vit01);
        EXPECT_EQ(Q[0].v2(), vit02);
        EXPECT_EQ(Q[0].p(),  v_ref01);
        EXPECT_EQ( (Q[0].v1()!=Q[1].v1() || Q[0].v2()!=Q[1].v2()) &&
                   (Q[0].v1()!=Q[2].v1() || Q[0].v2()!=Q[2].v2()) &&
                   (Q[1].v1()!=Q[2].v1() || Q[1].v2()!=Q[2].v2())    , true);

        EXPECT_EQ(  Q[1].p()==v_ref02 ||
                    Q[1].p()==v_ref03 ||
                    Q[1].p()==v_ref04 ||
                    Q[1].p()==v_ref05 ||
                    Q[1].p()==v_ref06 ||
                    Q[1].p()==v_ref07 ||
                    Q[1].p()==v_ref08 ||
                    Q[1].p()==v_ref09    , true );

        EXPECT_EQ(  Q[2].p()==v_ref02 ||
                    Q[2].p()==v_ref03 ||
                    Q[2].p()==v_ref04 ||
                    Q[2].p()==v_ref05 ||
                    Q[2].p()==v_ref06 ||
                    Q[2].p()==v_ref07 ||
                    Q[2].p()==v_ref08 ||
                    Q[2].p()==v_ref09    , true );

        EXPECT_EQ(  Q[1].p()!=Q[2].p(), true);

    }
}


/**  @brief findSupportDirectionOfTriangle()
 */
TEST_F(BDBoundarySimplexFinderTest, Test13) {

    for (long i = 0; i < 1000; i++) {

        auto M = randomRotMat();
        Vec3 trans(rand100()-50.0, rand100()-50.0, rand100()-50.0);

        Vec3 p01(  0.0,  0.0,  0.0);
        Vec3 p02(  rand100(),  0.0,  0.0);
        Vec3 p03(  0.0,  rand100(),  0.0);
        Vec3 v04(  0.0,  0.0,  1.0);

        p01 = M * p01 + trans;
        p02 = M * p02 + trans;
        p03 = M * p03 + trans;
        v04 = M * v04;

        ConvexRigidBody  b01(1);
        ConvexRigidBody  b02(2);
        vector<BDVertex> Q;
        VertexIt vit01;

        BDVertex mv01(vit01, vit01);
        mv01.mP = p01;
        BDVertex mv02(vit01, vit01);
        mv02.mP = p02;
        BDVertex mv03(vit01, vit01);
        mv03.mP = p03;

        Q.push_back(mv01);
        Q.push_back(mv02);
        Q.push_back(mv03);

        long             maxNumPivots = 100;
        long             maxNumCycles = 3;
        double           epsilonZero  = 1.0e-10;

        BDBoundarySimplexFinder 
           f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                   std::cerr);

        Vec3 v_res01 = f.findSupportDirectionOfTriangle();
        v_res01.normalize();
        if (v04.z() >= 0.0) {
            EXPECT_EQ(v_res01, v04);
        }
        else {
            EXPECT_EQ(v_res01, v04*-1.0);    
        }
    }
}


} // namespace Makena

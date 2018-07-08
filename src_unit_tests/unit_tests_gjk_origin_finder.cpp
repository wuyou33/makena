#include "gtest/gtest.h"
#include <vector>
#include "primitives.hpp"
#include "manifold.hpp"
#include "quaternion.hpp"
#include "gjk_origin_finder.hpp"

#include <iostream>
#include <fstream>
#include <sstream>


namespace Makena {


class GJKOriginFinderTest : public ::testing::Test {

  protected:  

    GJKOriginFinderTest(){;};
    virtual ~GJKOriginFinderTest(){;};
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
*/

static double rand100()
{
    return ((rand()%10000)+1)/100.0;
}


static bool isAnyOf(
    enum predicate pTest,
    enum predicate p1,
    enum predicate p2,
    enum predicate p3,
    enum predicate p4
) {
    return pTest==p1 || pTest==p2 || pTest==p3 ||pTest==p4;
}

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

static bool isAnyOf(
    enum predicate pTest,
    enum predicate p1,
    enum predicate p2
) {
    return pTest==p1 || pTest==p2;
}

/*
static bool isAnyOf(
    enum predicate pTest,
    enum predicate p1
) {
    return pTest==p1;
}
*/

/**  @brief isInWrongOrderForTetrahedron()
 */
TEST_F(GJKOriginFinderTest, Test01) {


    for (long i = 0; i < 1000; i++) {

        Vec3 trans((rand()%100)/10.0, (rand()%100)/10.0, (rand()%100)/10.0);
        Mat3x3 M = randomRotMat();

        Vec3 p01( 0.0, 0.0, 0.0);
        Vec3 p02( ((rand()%100)+1)/10.0, 0.0, 0.0);
        Vec3 p03( 0.0, ((rand()%100)+1)/10.0, 0.0);
        Vec3 p04( 0.0, 0.0, ((rand()%100)+1)/10.0);

        p01 = M * p01 + trans;
        p02 = M * p02 + trans;
        p03 = M * p03 + trans;
        p04 = M * p04 + trans;

        ConvexRigidBody         b01(1);
        ConvexRigidBody         b02(2);
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

        long   maxNumPivots = 100;
        long   maxNumCycles = 3;
        double epsilonZero   = 1.0e-10;

        GJKOriginFinder 
            f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);

        EXPECT_EQ(f.isInWrongOrderForTetrahedron(), false);
        std::swap(Q[1], Q[2]);
        EXPECT_EQ(f.isInWrongOrderForTetrahedron(), true);
    }
}


/**  @brief testPointAgainstTriangle(() random triangle some fixed test points.
 */
TEST_F(GJKOriginFinderTest, Test02) {


    for (long i = 0; i < 1000; i++) {

        Vec3 trans((rand()%100)/10.0, (rand()%100)/10.0, (rand()%100)/10.0);
        Mat3x3 M = randomRotMat();

        Vec3 p01( 0.0, 0.0, 0.0); // Origin
        Vec3 p02( ((rand()%100)+1)/10.0, 0.0, 0.0); // X
        Vec3 p03( 0.0, ((rand()%100)+1)/10.0, 0.0); // Y
        Vec3 p04( 0.0, 0.0, ((rand()%100)+1)/10.0); // Z
                                                    // (test against (Ori,X,Y))
        Vec3 p05( p02.x(), 0.0, ((rand()%100)+1)/10.0); // Over X
        Vec3 p06( 0.0, p03.y(), ((rand()%100)+1)/10.0); // Over Y
        Vec3 p07( p02.x()*0.25, p03.y()*0.25, ((rand()%100)+1)/10.0); 
                                                          // Over (Ori, X, Y)
        double scale01 = (rand()%98 + 1)/100.0;
        Vec3 p08 = p01 * scale01 + p02 * (1.0 - scale01); // Over (Orig, X)
        Vec3 p09 = p05 * scale01 + p06 * (1.0 - scale01); // Over (X,Y)
        Vec3 p10 = p06 * scale01 + p01 * (1.0 - scale01); // Over (Y, Orig)


        p01 = M * p01 + trans;
        p02 = M * p02 + trans;
        p03 = M * p03 + trans;
        p04 = M * p04 + trans;
        p05 = M * p05 + trans;
        p06 = M * p06 + trans;
        p07 = M * p07 + trans;
        p08 = M * p08 + trans;
        p09 = M * p09 + trans;
        p10 = M * p10 + trans;

        ConvexRigidBody         b01(1);
        ConvexRigidBody         b02(2);
        vector<BDVertex>        Q;
        long                    maxNumPivots = 100;
        long                    maxNumCycles = 3;
        double                  epsilonZero  = 1.0e-10;

        GJKOriginFinder
            f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                   std::cerr);
        f.setLogLevel(GJKOriginFinder::ALL);

        enum predicate pred01 = f.testPointAgainstTriangle(p01, p02, p03, p04);
        EXPECT_EQ( pred01 == VORONOI_INSIDE_TRIANGLE ||
                   pred01 == VORONOI_OVER_EDGE_3_1   ||
                   pred01 == VORONOI_OVER_EDGE_1_2   ||
                   pred01 == VORONOI_OVER_VERTEX_1      , true );

        enum predicate pred02 = f.testPointAgainstTriangle(p01, p02, p03, p05);
        EXPECT_EQ( pred02 == VORONOI_INSIDE_TRIANGLE ||
                   pred02 == VORONOI_OVER_EDGE_2_3   ||
                   pred02 == VORONOI_OVER_EDGE_1_2   ||
                   pred02 == VORONOI_OVER_VERTEX_2      , true );

        enum predicate pred03 = f.testPointAgainstTriangle(p01, p02, p03, p06);
        EXPECT_EQ( pred03 == VORONOI_INSIDE_TRIANGLE ||
                   pred03 == VORONOI_OVER_EDGE_2_3   ||
                   pred03 == VORONOI_OVER_EDGE_3_1   ||
                   pred03 == VORONOI_OVER_VERTEX_3      , true );

        enum predicate pred04 = f.testPointAgainstTriangle(p01, p02, p03, p07);
        EXPECT_EQ( pred04 == VORONOI_INSIDE_TRIANGLE , true );

        enum predicate pred05 = f.testPointAgainstTriangle(p01, p02, p03, p08);

        EXPECT_EQ( pred05 == VORONOI_INSIDE_TRIANGLE ||
                   pred05 == VORONOI_OVER_EDGE_1_2   , true );

        enum predicate pred06 = f.testPointAgainstTriangle(p01, p02, p03, p09);

        EXPECT_EQ( pred06 == VORONOI_INSIDE_TRIANGLE ||
                   pred06 == VORONOI_OVER_EDGE_2_3   , true );

        enum predicate pred07 = f.testPointAgainstTriangle(p01, p02, p03, p10);

        EXPECT_EQ( pred07 == VORONOI_INSIDE_TRIANGLE ||
                   pred07 == VORONOI_OVER_EDGE_3_1   , true );

        //f.logPredicate(GJKOriginFinder::INFO, pred07);

    }
}


static bool Test03Test3Points(
    Vec3             p1, 
    Vec3             p2, 
    Vec3             p3, 
    Vec3             pt, 
    const Mat3x3&    M,
    const Vec3&      trans,
    GJKOriginFinder& f,
    enum predicate   pred01,
    enum predicate   pred02,
    enum predicate   pred03,
    enum predicate   pred04
) {

    p1 = M * p1 + trans;
    p2 = M * p2 + trans;
    p3 = M * p3 + trans;

    pt.setZ(-100.0);
    Vec3 pt1 = M * pt + trans;
//    Vec3 pt1 = pt;
    pt.setZ(0.0);
    Vec3 pt2 = M * pt + trans;
//    Vec3 pt2 = pt;
    pt.setZ(100.0);
    Vec3 pt3 = M * pt + trans;
//    Vec3 pt3 = pt;

//    cerr << "p1: " << p1 << "\n";
//    cerr << "p2: " << p2 << "\n";
//    cerr << "p3: " << p3 << "\n";
//    cerr << "pt1: " << pt1 << "\n";
//    cerr << "pt2: " << pt2 << "\n";
//    cerr << "pt3: " << pt3 << "\n";

    enum predicate test1 = f.testPointAgainstTriangle(p1, p2, p3, pt1);
    enum predicate test2 = f.testPointAgainstTriangle(p1, p2, p3, pt2);
    enum predicate test3 = f.testPointAgainstTriangle(p1, p2, p3, pt3);

//    f.logPredicate(GJKOriginFinder::INFO, test1);
//    f.logPredicate(GJKOriginFinder::INFO, test2);
//    f.logPredicate(GJKOriginFinder::INFO, test3);

    bool res1 = false;
    bool res2 = false;
    bool res3 = false;

    if (pred01 != NONE&& test1==pred01) {
        res1 = true;
    }
    if (pred02 != NONE&& test1==pred02) {
        res1 = true;
    }
    if (pred03 != NONE&& test1==pred03) {
        res1 = true;
    }
    if (pred04 != NONE&& test1==pred04) {
        res1 = true;
    }

    if (pred01 != NONE&& test2==pred01) {
        res2 = true;
    }
    if (pred02 != NONE&& test2==pred02) {
        res2 = true;
    }
    if (pred03 != NONE&& test2==pred03) {
        res2 = true;
    }
    if (pred04 != NONE&& test1==pred04) {
        res2 = true;
    }

    if (pred01 != NONE&& test3==pred01) {
        res3 = true;
    }
    if (pred02 != NONE&& test3==pred02) {
        res3 = true;
    }
    if (pred03 != NONE&& test3==pred03) {
        res3 = true;
    }
    if (pred04 != NONE&& test3==pred04) {
        res3 = true;
    }

    return res1 && res2 && res3;
}


/**  @brief testPointAgainstTriangle(() one flat triangle with 
 *          various test points.
 */
TEST_F(GJKOriginFinderTest, Test03) {


    ConvexRigidBody         b01(1);
    ConvexRigidBody         b02(2);
    vector<BDVertex> Q;
    long                    maxNumPivots = 100;
    long                    maxNumCycles = 3;
    double                  epsilonZero  = 1.0e-10;
    GJKOriginFinder 
            f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0, 
                                                                   std::cerr);
    f.setLogLevel(GJKOriginFinder::ALL);

    for (long i = 0; i < 1000; i++) {

        Vec3 trans((rand()%100)/10.0, (rand()%100)/10.0, (rand()%100)/10.0);
        Mat3x3 M = randomRotMat();

        Vec3 p01(-100.0, 0.0, 0.0);
        Vec3 p02( 100.0, 0.0, 0.0);
        Vec3 p03(   0.0, 0.1, 0.0);

        Vec3 t01(-100.0, 0.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t01, M, trans, f,
                      VORONOI_INSIDE_TRIANGLE,
                      VORONOI_OVER_EDGE_1_2,
                      VORONOI_OVER_EDGE_3_1,
                      VORONOI_OVER_VERTEX_1), true );

        Vec3 t02( 100.0, 0.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t02, M, trans, f,
                      VORONOI_INSIDE_TRIANGLE,
                      VORONOI_OVER_EDGE_1_2,
                      VORONOI_OVER_EDGE_2_3,
                      VORONOI_OVER_VERTEX_2), true );

        Vec3 t03( 0.0, 0.1,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t03, M, trans, f,
                      VORONOI_INSIDE_TRIANGLE,
                      VORONOI_OVER_EDGE_2_3,
                      VORONOI_OVER_EDGE_3_1,
                      VORONOI_OVER_VERTEX_3), true );

        Vec3 t04 = p01 * 0.3 + p02 * 0.7;
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t04, M, trans, f,
                      VORONOI_INSIDE_TRIANGLE,
                      VORONOI_OVER_EDGE_1_2,
                      NONE,
                      NONE                 ), true );

        Vec3 t05 = p02 * 0.3 + p03 * 0.7;
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t05, M, trans, f,
                      VORONOI_INSIDE_TRIANGLE,
                      VORONOI_OVER_EDGE_2_3,
                      NONE,
                      NONE                 ), true );

        Vec3 t06 = p03 * 0.3 + p01 * 0.7;
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t06, M, trans, f,
                      VORONOI_INSIDE_TRIANGLE,
                      VORONOI_OVER_EDGE_3_1,
                      NONE,
                      NONE                 ), true );

        Vec3 t07 = p01 * 0.5 + p02 * 0.4 + p03 * 0.1;
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t07, M, trans, f,
                      VORONOI_INSIDE_TRIANGLE,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t08(-100.0, -1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t08, M, trans, f,
                      VORONOI_OVER_EDGE_1_2,
                      VORONOI_OVER_VERTEX_1,
                      NONE,
                      NONE                 ), true );

        Vec3 t09(-99.0, -1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t09, M, trans, f,
                      VORONOI_OVER_EDGE_1_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t10(-100.1, -1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t10, M, trans, f,
                      VORONOI_OVER_VERTEX_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t11( -101.0, -1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t11, M, trans, f,
                      VORONOI_OVER_VERTEX_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t12( -101.0, 0.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t12, M, trans, f,
                      VORONOI_OVER_VERTEX_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t13( -101.0, 1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t13, M, trans, f,
                      VORONOI_OVER_VERTEX_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t14( -100.001, 1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t14, M, trans, f,
                      VORONOI_OVER_EDGE_3_1,
                      VORONOI_OVER_VERTEX_1,
                      NONE,
                      NONE                 ), true );

        Vec3 t15( -100.0, 1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t15, M, trans, f,
                      VORONOI_OVER_EDGE_3_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t16(  99.9,-1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t16, M, trans, f,
                      VORONOI_OVER_EDGE_1_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t17(  100.0,-1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t17, M, trans, f,
                      VORONOI_OVER_EDGE_1_2,
                      VORONOI_OVER_VERTEX_2,
                      NONE,
                      NONE                 ), true );

        Vec3 t18(  100.001,-1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t18, M, trans, f,
                      VORONOI_OVER_VERTEX_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t19(  101.0,-1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t19, M, trans, f,
                      VORONOI_OVER_VERTEX_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t20(  101.0, 0.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t20, M, trans, f,
                      VORONOI_OVER_VERTEX_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t21(  101.0, 1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t21, M, trans, f,
                      VORONOI_OVER_VERTEX_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t22(  100.001, 1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t22, M, trans, f,
                      VORONOI_OVER_EDGE_2_3,
                      VORONOI_OVER_VERTEX_2,
                      NONE,
                      NONE                 ), true );

        Vec3 t23(  100.0, 1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t23, M, trans, f,
                      VORONOI_OVER_EDGE_2_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t24(  98.0, 1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t24, M, trans, f,
                      VORONOI_OVER_EDGE_2_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t25(   0.0, 1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t25, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t26(   0.0, 0.2,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t26, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t27(   0.0,1000.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t27, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t28( 0.001, 1.1,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t28, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      VORONOI_OVER_EDGE_2_3,
                      NONE,
                      NONE                 ), true );

        Vec3 t29( 0.002, 1.1,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t29, M, trans, f,
                      VORONOI_OVER_EDGE_2_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t30(-0.001, 1.1,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t30, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      VORONOI_OVER_EDGE_3_1,
                      NONE,
                      NONE                 ), true );

        Vec3 t31(-0.002, 1.1,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t31, M, trans, f,
                      VORONOI_OVER_EDGE_3_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t32( 1.0, 2.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t32, M, trans, f,
                      VORONOI_OVER_EDGE_2_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t33(-1.0, 2.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t33, M, trans, f,
                      VORONOI_OVER_EDGE_3_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

    }
}


/**  @brief testPointAgainstTriangle(() one sharp triangle with 
 *          various test points.
 */
TEST_F(GJKOriginFinderTest, Test04) {


    ConvexRigidBody         b01(1);
    ConvexRigidBody         b02(2);
    vector<BDVertex> Q;
    long                    maxNumPivots = 100;
    long                    maxNumCycles = 3;
    double                  epsilonZero  = 1.0e-10;
    GJKOriginFinder 
             f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);
    f.setLogLevel(GJKOriginFinder::OFF);

    for (long i = 0; i < 1000; i++) {

        Vec3 trans((rand()%100)/10.0, (rand()%100)/10.0, (rand()%100)/10.0);
        Mat3x3 M = randomRotMat();

        Vec3 p01(-0.1,   0.0, 0.0);
        Vec3 p02( 0.1,   0.0, 0.0);
        Vec3 p03( 0.0,1000.0, 0.0);

        Vec3 t01(-0.1,  -1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t01, M, trans, f,
                      VORONOI_OVER_EDGE_1_2,
                      VORONOI_OVER_VERTEX_1,
                      NONE,
                      NONE                 ), true );

        Vec3 t02(-0.1001,  -1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t02, M, trans, f,
                      VORONOI_OVER_VERTEX_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t03(-0.0999,  -1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t03, M, trans, f,
                      VORONOI_OVER_EDGE_1_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t04(-2.0,  -1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t04, M, trans, f,
                      VORONOI_OVER_VERTEX_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t05(-2.0,   0.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t05, M, trans, f,
                      VORONOI_OVER_VERTEX_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t06(-1.01,   0.0001,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t06, M, trans, f,
                      VORONOI_OVER_VERTEX_1,
                      VORONOI_OVER_EDGE_3_1,
                      NONE,
                      NONE                 ), true );

        Vec3 t07(-1.01,   0.00011,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t07, M, trans, f,
                      VORONOI_OVER_EDGE_3_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t08(-1.01,   0.00009,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t08, M, trans, f,
                      VORONOI_OVER_VERTEX_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t09(   0.0,  -0.5,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t09, M, trans, f,
                      VORONOI_OVER_EDGE_1_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t10(   0.0,   0.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t10, M, trans, f,
                      VORONOI_OVER_EDGE_1_2,
                      VORONOI_INSIDE_TRIANGLE,
                      NONE,
                      NONE                 ), true );

        Vec3 t11(   0.0, -0.00001,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t11, M, trans, f,
                      VORONOI_OVER_EDGE_1_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t12(   0.0, 0.00001,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t12, M, trans, f,
                      VORONOI_INSIDE_TRIANGLE,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t13(   0.1, -1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t13, M, trans, f,
                      VORONOI_OVER_EDGE_1_2,
                      VORONOI_OVER_VERTEX_2,
                      NONE,
                      NONE                 ), true );

        Vec3 t14(   0.09, -1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t14, M, trans, f,
                      VORONOI_OVER_EDGE_1_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t15(   0.101, -1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t15, M, trans, f,
                      VORONOI_OVER_VERTEX_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t16(   2.0, -1.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t16, M, trans, f,
                      VORONOI_OVER_VERTEX_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t17(   2.0, 0.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t17, M, trans, f,
                      VORONOI_OVER_VERTEX_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t18(   1.1, 0.0001,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t18, M, trans, f,
                      VORONOI_OVER_EDGE_2_3,
                      VORONOI_OVER_VERTEX_2,
                      NONE,
                      NONE                 ), true );

        Vec3 t19(   1.1, 0.00009,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t19, M, trans, f,
                      VORONOI_OVER_VERTEX_2,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t20(   1.1, 0.00011,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t20, M, trans, f,
                      VORONOI_OVER_EDGE_2_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t21(   0.05, 500.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t21, M, trans, f,
                      VORONOI_OVER_EDGE_2_3,
                      VORONOI_INSIDE_TRIANGLE,
                      NONE,
                      NONE                 ), true );

        Vec3 t22( 0.049,  500.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t22, M, trans, f,
                      VORONOI_INSIDE_TRIANGLE,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t23( 0.051,  500.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t23, M, trans, f,
                      VORONOI_OVER_EDGE_2_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t24( 1.1, 1000.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t24, M, trans, f,
                      VORONOI_OVER_EDGE_2_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t25( 1.0, 1000.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t25, M, trans, f,
                      VORONOI_OVER_EDGE_2_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t26( 1.0, 1000.0001, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t26, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      VORONOI_OVER_EDGE_2_3,
                      NONE,
                      NONE                 ), true );

        Vec3 t27( 1.0, 1000.00011, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t27, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t28( 0.5, 1002.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t28, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t29( 0.0, 1100.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t29, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t30( 0.1, 2000.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t30, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t31( 0.0, 2000.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t31, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t32(-0.1, 2000.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t32, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t33(-0.5, 1002.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t33, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t34(-1.0, 1000.00011, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t34, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t35(-1.0, 1000.0001, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t35, M, trans, f,
                      VORONOI_OVER_VERTEX_3,
                      VORONOI_OVER_EDGE_3_1,
                      NONE,
                      NONE                 ), true );

        Vec3 t36(-1.0, 1000.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t36, M, trans, f,
                      VORONOI_OVER_EDGE_3_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t37(-1.1, 1000.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t37, M, trans, f,
                      VORONOI_OVER_EDGE_3_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t38(-0.051,  500.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t38, M, trans, f,
                      VORONOI_OVER_EDGE_3_1,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t39( -0.05, 500.0,  0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t39, M, trans, f,
                      VORONOI_OVER_EDGE_3_1,
                      VORONOI_INSIDE_TRIANGLE,
                      NONE,
                      NONE                 ), true );

        Vec3 t40(-0.049,  500.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t40, M, trans, f,
                      VORONOI_INSIDE_TRIANGLE,
                      NONE,
                      NONE,
                      NONE                 ), true );

        Vec3 t41( 0.0,  500.0, 0.0);
        EXPECT_EQ(Test03Test3Points(p01, p02, p03, t41, M, trans, f,
                      VORONOI_INSIDE_TRIANGLE,
                      NONE,
                      NONE,
                      NONE                 ), true );


    }
}

/**  @brief projectPointOnLine()
 */
TEST_F(GJKOriginFinderTest, Test05) {

    ConvexRigidBody         b01(1);
    ConvexRigidBody         b02(2);
    vector<BDVertex> Q;
    long   maxNumPivots = 100;
    long   maxNumCycles = 3;
    double epsilonZero   = 1.0e-10;
    GJKOriginFinder 
            f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                   std::cerr);

    for (long i = 0; i < 1000; i++) {

        Vec3 trans((rand()%100)/10.0, (rand()%100)/10.0, (rand()%100)/10.0);
        Mat3x3 M = randomRotMat();

        Vec3 p01(       ((rand()%100)+1)/10.0, 0.0, 0.0);
        Vec3 p02(-1.0 * ((rand()%100)+1)/10.0, 0.0, 0.0);

        Vec3 pt01  ( 0.0, ((rand()%100)+1)/10.0, 0.0);
        Vec3 exp01 ( 0.0, 0.0, 0.0);

        Vec3 pt02  ( 0.0, 0.0, ((rand()%100)+1)/10.0);
        Vec3 exp02 ( 0.0, 0.0, 0.0);

        Vec3 pt03  ( 0.0, 0.0, 0.0);
        Vec3 exp03 ( 0.0, 0.0, 0.0);

        Vec3 pt04  ( p01.x(), ((rand()%100)+1)/10.0, 0.0);
        Vec3 exp04 ( p01.x(), 0.0, 0.0);

        Vec3 pt05  ( p01.x(), 0.0, ((rand()%100)+1)/10.0);
        Vec3 exp05 ( p01.x(), 0.0, 0.0);

        Vec3 pt06  ( p01.x(), 0.0, 0.0);
        Vec3 exp06 ( p01.x(), 0.0, 0.0);

        Vec3 pt07  ( p02.x(), ((rand()%100)+1)/10.0, 0.0);
        Vec3 exp07 ( p02.x(), 0.0, 0.0);

        Vec3 pt08  ( p02.x(), 0.0, ((rand()%100)+1)/10.0);
        Vec3 exp08 ( p02.x(), 0.0, 0.0);

        Vec3 pt09  ( p02.x(), 0.0, 0.0);
        Vec3 exp09 ( p02.x(), 0.0, 0.0);

        Vec3 pt10  ( ((rand()%100)+1)/10.0, 
                     ((rand()%100)+1)/10.0, ((rand()%100)+1)/10.0  );

        Vec3 exp10 ( pt10.x(), 0.0, 0.0);
                    

        p01 = M * p01 + trans;
        p02 = M * p02 + trans;

        pt01 = M * pt01 + trans;
        exp01 = M * exp01 + trans;

        pt02 = M * pt02 + trans;
        exp02 = M * exp02 + trans;

        pt03 = M * pt03 + trans;
        exp03 = M * exp03 + trans;

        pt04 = M * pt04 + trans;
        exp04 = M * exp04 + trans;

        pt05 = M * pt05 + trans;
        exp05 = M * exp05 + trans;

        pt06 = M * pt06 + trans;
        exp06 = M * exp06 + trans;

        pt07 = M * pt07 + trans;
        exp07 = M * exp07 + trans;

        pt08 = M * pt08 + trans;
        exp08 = M * exp08 + trans;

        pt09 = M * pt09 + trans;
        exp09 = M * exp09 + trans;

        pt10 = M * pt10 + trans;
        exp10 = M * exp10 + trans;

        Vec3 proj01 = f.projectPointOnLine(p01, p02, pt01);
        EXPECT_EQ(proj01, exp01);

        Vec3 proj02 = f.projectPointOnLine(p01, p02, pt02);
        EXPECT_EQ(proj02, exp02);

        Vec3 proj03 = f.projectPointOnLine(p01, p02, pt03);
        EXPECT_EQ(proj03, exp03);

        Vec3 proj04 = f.projectPointOnLine(p01, p02, pt04);
        EXPECT_EQ(proj04, exp04);

        Vec3 proj05 = f.projectPointOnLine(p01, p02, pt05);
        EXPECT_EQ(proj05, exp05);

        Vec3 proj06 = f.projectPointOnLine(p01, p02, pt06);
        EXPECT_EQ(proj06, exp06);

        Vec3 proj07 = f.projectPointOnLine(p01, p02, pt07);
        EXPECT_EQ(proj07, exp07);

        Vec3 proj08 = f.projectPointOnLine(p01, p02, pt08);
        EXPECT_EQ(proj08, exp08);

        Vec3 proj09 = f.projectPointOnLine(p01, p02, pt09);
        EXPECT_EQ(proj09, exp09);

        Vec3 proj10 = f.projectPointOnLine(p01, p02, pt10);
        EXPECT_EQ(proj10, exp10);

    }
}


/**  @brief  testPointAgainstEdge()
 */

TEST_F(GJKOriginFinderTest, Test06) {

    ConvexRigidBody         b01(1);
    ConvexRigidBody         b02(2);
    vector<BDVertex> Q;
    long   maxNumPivots = 100;
    long   maxNumCycles = 3;
    double epsilonZero   = 1.0e-10;
    GJKOriginFinder 
            f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                   std::cerr);

    for (long i = 0; i < 1000; i++) {

        Vec3 trans((rand()%100)/10.0, (rand()%100)/10.0, (rand()%100)/10.0);
        Mat3x3 M = randomRotMat();

        Vec3 p01(       ((rand()%100)+1)/10.0, 0.0, 0.0);
        Vec3 p02(-1.0 * ((rand()%100)+1)/10.0, 0.0, 0.0);

        Vec3 pt01  ( 0.0, ((rand()%100)+1)/10.0, 0.0);
        Vec3 exp01 ( 0.0, 0.0, 0.0);

        Vec3 pt02  ( 0.0, 0.0, ((rand()%100)+1)/10.0);
        Vec3 exp02 ( 0.0, 0.0, 0.0);

        Vec3 pt03  ( 0.0, 0.0, 0.0);
        Vec3 exp03 ( 0.0, 0.0, 0.0);

        Vec3 pt04  ( p01.x(), ((rand()%100)+1)/10.0, 0.0);
        Vec3 exp04 ( p01.x(), 0.0, 0.0);

        Vec3 pt05  ( p01.x(), 0.0, ((rand()%100)+1)/10.0);
        Vec3 exp05 ( p01.x(), 0.0, 0.0);

        Vec3 pt06  ( p01.x(), 0.0, 0.0);
        Vec3 exp06 ( p01.x(), 0.0, 0.0);

        Vec3 pt07  ( p02.x(), ((rand()%100)+1)/10.0, 0.0);
        Vec3 exp07 ( p02.x(), 0.0, 0.0);

        Vec3 pt08  ( p02.x(), 0.0, ((rand()%100)+1)/10.0);
        Vec3 exp08 ( p02.x(), 0.0, 0.0);

        Vec3 pt09  ( p02.x(), 0.0, 0.0);
        Vec3 exp09 ( p02.x(), 0.0, 0.0);

        Vec3 pt10  ( ((rand()%100)+1)/10.0, 
                     ((rand()%100)+1)/10.0, ((rand()%100)+1)/10.0  );

        Vec3 exp10 ( std::max(std::min(pt10.x(), p01.x()), p02.x()),
                     0.0, 0.0);
        if (fabs(pt10.x() - p01.x()) < 0.0000001){
             pt10.setX(pt10.x() + 0.0000001);
             exp10.setX(p01.x());
        }

        if (fabs(pt10.x() - p02.x()) < 0.0000001){
             pt10.setX(pt10.x() - 0.0000001);
             exp10.setX(p02.x());
        }

        enum predicate expPred10;
        if (pt10.x() > p01.x()) {
            expPred10 = ON_POINT1;
        }
        else if (pt10.x() == p01.x()) {
            expPred10 = BETWEEN_1_AND_2;
        }
        else if (pt10.x() == p02.x()) {
            expPred10 = BETWEEN_1_AND_2;
        }
        else if (pt10.x() < p02.x()) {
            expPred10 = ON_POINT2;
        }
        else {
            expPred10 = BETWEEN_1_AND_2;
        }

        p01 = M * p01 + trans;
        p02 = M * p02 + trans;

        pt01 = M * pt01 + trans;
        exp01 = M * exp01 + trans;

        pt02 = M * pt02 + trans;
        exp02 = M * exp02 + trans;

        pt03 = M * pt03 + trans;
        exp03 = M * exp03 + trans;

        pt04 = M * pt04 + trans;
        exp04 = M * exp04 + trans;

        pt05 = M * pt05 + trans;
        exp05 = M * exp05 + trans;

        pt06 = M * pt06 + trans;
        exp06 = M * exp06 + trans;

        pt07 = M * pt07 + trans;
        exp07 = M * exp07 + trans;

        pt08 = M * pt08 + trans;
        exp08 = M * exp08 + trans;

        pt09 = M * pt09 + trans;
        exp09 = M * exp09 + trans;

        pt10 = M * pt10 + trans;
        exp10 = M * exp10 + trans;

        Vec3 proj01; 
        enum predicate pred01=f.testPointAgainstEdge(p01, p02, pt01, proj01);
        EXPECT_EQ(pred01, BETWEEN_1_AND_2);
        EXPECT_EQ(proj01, exp01);

        Vec3 proj02; 
        enum predicate pred02=f.testPointAgainstEdge(p01, p02, pt02, proj02);
        EXPECT_EQ(pred02, BETWEEN_1_AND_2);
        EXPECT_EQ(proj02, exp03);

        Vec3 proj03; 
        enum predicate pred03=f.testPointAgainstEdge(p01, p02, pt03, proj03);
        EXPECT_EQ(pred03, BETWEEN_1_AND_2);
        EXPECT_EQ(proj03, exp03);

        Vec3 proj04; 
        enum predicate pred04=f.testPointAgainstEdge(p01, p02, pt04, proj04);
        EXPECT_EQ(pred04==BETWEEN_1_AND_2 ||pred04==ON_POINT1, true);
        EXPECT_EQ(proj04, exp04);

        Vec3 proj05; 
        enum predicate pred05=f.testPointAgainstEdge(p01, p02, pt05, proj05);
        EXPECT_EQ(pred05==BETWEEN_1_AND_2 ||pred05==ON_POINT1, true);
        EXPECT_EQ(proj05, exp05);

        Vec3 proj06; 
        enum predicate pred06=f.testPointAgainstEdge(p01, p02, pt06, proj06);
        EXPECT_EQ(pred06==BETWEEN_1_AND_2 ||pred06==ON_POINT1, true);
        EXPECT_EQ(proj06, exp06);

        Vec3 proj07; 
        enum predicate pred07=f.testPointAgainstEdge(p01, p02, pt07, proj07);
        EXPECT_EQ(pred07==BETWEEN_1_AND_2 ||pred07==ON_POINT2, true);
        EXPECT_EQ(proj07, exp07);

        Vec3 proj08; 
        enum predicate pred08=f.testPointAgainstEdge(p01, p02, pt08, proj08);
        EXPECT_EQ(pred08==BETWEEN_1_AND_2 ||pred08==ON_POINT2, true);
        EXPECT_EQ(proj08, exp08);

        Vec3 proj09; 
        enum predicate pred09=f.testPointAgainstEdge(p01, p02, pt09, proj09);
        EXPECT_EQ(pred09==BETWEEN_1_AND_2 ||pred09==ON_POINT2, true);
        EXPECT_EQ(proj09, exp09);

        Vec3 proj10; 
        enum predicate pred10=f.testPointAgainstEdge(p01, p02, pt10, proj10);
        EXPECT_EQ(pred10==expPred10, true);
        EXPECT_EQ(proj10, exp10);

    }
}


/**  @brief  projectPointOnLine()
 */
TEST_F(GJKOriginFinderTest, Test07) {

    ConvexRigidBody         b01(1);
    ConvexRigidBody         b02(2);
    vector<BDVertex> Q;
    long   maxNumPivots = 100;
    long   maxNumCycles = 3;
    double epsilonZero   = 1.0e-10;
    GJKOriginFinder
            f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0, 
                                                                    std::cerr);
    f.setLogLevel(GJKOriginFinder::ALL);
    for (long i = 0; i < 1000; i++) {

        Vec3 trans((rand()%100)/10.0, (rand()%100)/10.0, (rand()%100)/10.0);
        Mat3x3 M = randomRotMat();

        Vec3 p01(       (rand()%100)/10.0, 0.0, 0.0);
        Vec3 p02( -1.0* (rand()%100)/10.0, 0.0, 0.0);

        auto v12 = p02 - p01;
        Vec3 pt01((rand()%100)/10.0, 
                  (rand()%100)/10.0,
                  (rand()%100)/10.0  );

        Vec3 exp01;
        if (v12.squaredNorm2() < EPSILON_SQUARED) {
            exp01 = p01;
        }
        else {
            exp01 = Vec3( pt01.x(), 0.0, 0.0);
        }

        p01 = M * p01 + trans;
        p02 = M * p02 + trans;
        pt01 = M * pt01 + trans;
        exp01 = M * exp01 + trans;

        Vec3 proj01 = f.projectPointOnLine(p01, p02, pt01);

        EXPECT_EQ(proj01, exp01);

    }
}


/**  @brief  findClosestPointToOrigin0Simplex(
 */
TEST_F(GJKOriginFinderTest, Test08) {

    ConvexRigidBody         b01(1);
    ConvexRigidBody         b02(2);
    vector<BDVertex> Q;
    Vec3 zero(0.0, 0.0, 0.0);
    VertexIt vit01;
    BDVertex mv01(vit01, vit01);
    mv01.mP = zero;
    Q.push_back(mv01);

    long   maxNumPivots = 100;
    long   maxNumCycles = 3;
    double epsilonZero  = 1.0e-10;

    GJKOriginFinder
            f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);

    f.setLogLevel(GJKOriginFinder::ALL);
    for (long i = 0; i < 1000; i++) {

        Vec3 trans((rand()%100)/10.0, (rand()%100)/10.0, (rand()%100)/10.0);
        Mat3x3 M = randomRotMat();

        Vec3 p01(rand100(),  rand100(),  rand100());
        p01 = M * p01 + trans;

        f.mQ[0].mP = p01;


        enum predicate pred01;
        Vec3 res01;

        f.findClosestPointToOrigin0Simplex(res01, pred01);

        EXPECT_EQ(pred01, ON_POINT1);
        EXPECT_EQ(res01,  p01);

    }
}


/**  @brief  findClosestPointToOrigin1Simplex(
 */
TEST_F(GJKOriginFinderTest, Test09) {

    ConvexRigidBody         b01(1);
    ConvexRigidBody         b02(2);
    vector<BDVertex> Q;
    Vec3 zero(0.0, 0.0, 0.0);
    VertexIt vit01;
    BDVertex mv01(vit01, vit01);
    mv01.mP = zero;
    Q.push_back(mv01);
    BDVertex mv02(vit01, vit01);
    mv02.mP = zero;
    Q.push_back(mv02);

    long   maxNumPivots = 100;
    long   maxNumCycles = 3;
    double epsilonZero  = 1.0e-10;

    GJKOriginFinder
            f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);

    f.setLogLevel(GJKOriginFinder::ALL);
    for (long i = 0; i < 100000; i++) {

        Mat3x3 M = randomRotMat();
        Vec3 p01( rand100(),  rand100(), rand100() );
        Vec3 p02( p01.x() + rand100(),  p01.y(), p01.z());
        p01 = M * p01;
        p02 = M * p02;
        f.mQ[0].mP = p01;
        f.mQ[1].mP = p02;
        enum predicate pred01;
        Vec3 res01;

        f.findClosestPointToOrigin1Simplex(res01, pred01);

        EXPECT_EQ(pred01, ON_POINT1);
        EXPECT_EQ(res01,  p01);

        Vec3 p03(    -1.0 * rand100(),  rand100(), rand100() );
        Vec3 p04( p03.x() - rand100(),    p03.y(),   p03.z() );
        p03 = M * p03;
        p04 = M * p04;
        f.mQ[0].mP = p04;
        f.mQ[1].mP = p03;

        enum predicate pred02;
        Vec3 res02;

        f.findClosestPointToOrigin1Simplex(res02, pred02);

        EXPECT_EQ(pred02, ON_POINT2);
        EXPECT_EQ(res02,  p03);


        Vec3 p05(     rand100(),  rand100(), rand100() );
        Vec3 p06(-1.0*rand100(),    p05.y(),   p05.z() );
        Vec3 exp01(     0.0,    p05.y(),   p05.z() );
        p05   = M * p05;
        p06   = M * p06;
        exp01 = M * exp01;
        f.mQ[0].mP = p05;
        f.mQ[1].mP = p06;

        enum predicate pred03;
        Vec3 res03;

        f.findClosestPointToOrigin1Simplex(res03, pred03);

        EXPECT_EQ(pred03, BETWEEN_1_AND_2);
        EXPECT_EQ(res03,  exp01);


    }
}


/**  @brief  findClosestPointToOrigin2Simplex(
 */
TEST_F(GJKOriginFinderTest, Test10) {

    ConvexRigidBody         b01(1);
    ConvexRigidBody         b02(2);
    vector<BDVertex> Q;
    Vec3 zero(0.0, 0.0, 0.0);
    VertexIt vit01;
    BDVertex mv01(vit01, vit01);
    mv01.mP = zero;
    Q.push_back(mv01);
    BDVertex mv02(vit01, vit01);
    mv02.mP = zero;
    Q.push_back(mv02);
    BDVertex mv03(vit01, vit01);
    mv03.mP = zero;
    Q.push_back(mv03);

    long   maxNumPivots = 100;
    long   maxNumCycles = 3;
    double epsilonZero  = 1.0e-10;

    GJKOriginFinder
            f(b01, b02, Q, maxNumPivots, maxNumCycles, epsilonZero, true, 1.0,
                                                                    std::cerr);

    f.setLogLevel(GJKOriginFinder::ALL);
    for (long i = 0; i < 100000; i++) {

        Mat3x3 M = randomRotMat();

        Vec3 b01(       0.0,       0.0, rand100());
        Vec3 b02( rand100(),       0.0,   b01.z());
        Vec3 b03(       0.0, rand100(),   b01.z());

        ////

        Vec3 shift01(-1.0*b02.x()/4.0, -1.0*b03.y()/4.0, 0.0);
        Vec3 exp01(       0.0,       0.0,   b01.z());

        f.mQ[0].mP = M * (b01+shift01);
        f.mQ[1].mP = M * (b02+shift01);
        f.mQ[2].mP = M * (b03+shift01);
        exp01 = M * exp01;

        enum predicate pred01;
        Vec3 res01;
        f.findClosestPointToOrigin2Simplex(res01, pred01);
        EXPECT_EQ(pred01, VORONOI_INSIDE_TRIANGLE);
        EXPECT_EQ(res01,  exp01);

        ////

        Vec3 shift02(-1.0*b02.x()/2.0, 0.0, 0.0);
        Vec3 exp02(       0.0,       0.0,   b01.z());

        f.mQ[0].mP = M * (b01+shift02);
        f.mQ[1].mP = M * (b02+shift02);
        f.mQ[2].mP = M * (b03+shift02);
        exp02 = M * exp02;

        enum predicate pred02;
        Vec3 res02;
        f.findClosestPointToOrigin2Simplex(res02, pred02);
        //        f.logPredicate(GJKOriginFinder::INFO, pred02);
        EXPECT_EQ(isAnyOf(pred02, VORONOI_INSIDE_TRIANGLE,
                                  VORONOI_OVER_EDGE_1_2), true);
        EXPECT_EQ(res02,  exp02);

        ////

        Vec3 shift03(-1.0*b02.x()/2.0, -1.0*b03.y()/2.0, 0.0);
        Vec3 exp03(       0.0,       0.0,   b01.z());

        f.mQ[0].mP = M * (b01+shift03);
        f.mQ[1].mP = M * (b02+shift03);
        f.mQ[2].mP = M * (b03+shift03);
        exp03 = M * exp03;

        enum predicate pred03;
        Vec3 res03;
        f.findClosestPointToOrigin2Simplex(res03, pred03);
        EXPECT_EQ(isAnyOf(pred03, VORONOI_INSIDE_TRIANGLE,
                                  VORONOI_OVER_EDGE_2_3), true);
        EXPECT_EQ(res03,  exp03);

        ////

        Vec3 shift04(  0.0, -1.0*b03.y()/2.0, 0.0);
        Vec3 exp04(       0.0,       0.0,   b01.z());

        f.mQ[0].mP = M * (b01+shift04);
        f.mQ[1].mP = M * (b02+shift04);
        f.mQ[2].mP = M * (b03+shift04);
        exp04 = M * exp04;

        enum predicate pred04;
        Vec3 res04;
        f.findClosestPointToOrigin2Simplex(res04, pred04);
        EXPECT_EQ(isAnyOf(pred04, VORONOI_INSIDE_TRIANGLE,
                                  VORONOI_OVER_EDGE_3_1), true);
        EXPECT_EQ(res04,  exp04);

        ////

        Vec3 shift05(0.0, 0.0, 0.0);
        Vec3 exp05(       0.0,       0.0,   b01.z());

        f.mQ[0].mP = M * (b01+shift05);
        f.mQ[1].mP = M * (b02+shift05);
        f.mQ[2].mP = M * (b03+shift05);
        exp05 = M * exp05;

        enum predicate pred05;
        Vec3 res05;
        f.findClosestPointToOrigin2Simplex(res05, pred05);
        EXPECT_EQ(isAnyOf(pred05, VORONOI_INSIDE_TRIANGLE,
                                  VORONOI_OVER_EDGE_3_1, 
                                  VORONOI_OVER_EDGE_1_2, 
                                  VORONOI_OVER_VERTEX_1  ), true);
        //f.logPredicate(GJKOriginFinder::INFO, pred05);
        EXPECT_EQ(res05,  exp05);

        ////

        Vec3 shift06(-1.0*b02.x(), 0.0, 0.0);
        Vec3 exp06(       0.0,       0.0,   b01.z());

        f.mQ[0].mP = M * (b01+shift06);
        f.mQ[1].mP = M * (b02+shift06);
        f.mQ[2].mP = M * (b03+shift06);
        exp06 = M * exp06;

        enum predicate pred06;
        Vec3 res06;
        f.findClosestPointToOrigin2Simplex(res06, pred06);
        EXPECT_EQ(isAnyOf(pred06, VORONOI_INSIDE_TRIANGLE,
                                  VORONOI_OVER_EDGE_1_2, 
                                  VORONOI_OVER_EDGE_2_3, 
                                  VORONOI_OVER_VERTEX_2  ), true);
        //f.logPredicate(GJKOriginFinder::INFO, pred05);
        EXPECT_EQ(res06,  exp06);

        ////

        Vec3 shift07( 0.0, -1.0*b03.y(), 0.0);
        Vec3 exp07(       0.0,       0.0,   b01.z());

        f.mQ[0].mP = M * (b01+shift07);
        f.mQ[1].mP = M * (b02+shift07);
        f.mQ[2].mP = M * (b03+shift07);
        exp07 = M * exp07;

        enum predicate pred07;
        Vec3 res07;
        f.findClosestPointToOrigin2Simplex(res07, pred07);
        EXPECT_EQ(isAnyOf(pred07, VORONOI_INSIDE_TRIANGLE,
                                  VORONOI_OVER_EDGE_2_3, 
                                  VORONOI_OVER_EDGE_3_1, 
                                  VORONOI_OVER_VERTEX_3  ), true);
        //f.logPredicate(GJKOriginFinder::INFO, pred05);
        EXPECT_EQ(res07,  exp07);

        ////

        Vec3 shift08( -1.0*b02.x()/2.0, rand100(),   0.0);
        Vec3 exp08  (              0.0, shift08.y(), b01.z());

        f.mQ[0].mP = M * (b01+shift08);
        f.mQ[1].mP = M * (b02+shift08);
        f.mQ[2].mP = M * (b03+shift08);
        exp08 = M * exp08;

        enum predicate pred08;
        Vec3 res08;
        f.findClosestPointToOrigin2Simplex(res08, pred08);
        EXPECT_EQ(pred08, VORONOI_OVER_EDGE_1_2);
        //f.logPredicate(GJKOriginFinder::INFO, pred08);
        EXPECT_EQ(res08,  exp08);

        EXPECT_EQ(res07,  exp07);

        ////
        double alpha09 = rand100();
        Vec3 shift09(
            -1.0*b02.x()/2.0 - alpha09 * b03.y(),
            -1.0*b03.y()/2.0 - alpha09 * b02.x(), 
            0.0
        );

        Vec3 exp09  (-1.0*alpha09*b03.y(), -1.0*alpha09*b02.x(), b01.z());

        f.mQ[0].mP = M * (b01+shift09);
        f.mQ[1].mP = M * (b02+shift09);
        f.mQ[2].mP = M * (b03+shift09);
        exp09 = M * exp09;

        enum predicate pred09;
        Vec3 res09;
        f.findClosestPointToOrigin2Simplex(res09, pred09);
        EXPECT_EQ(pred09, VORONOI_OVER_EDGE_2_3);
        //f.logPredicate(GJKOriginFinder::INFO, pred09);
        EXPECT_EQ(res09,  exp09);

        ////

        Vec3 shift10(rand100(),   -1.0*b03.y()/2.0, 0.0);
        Vec3 exp10  (shift10.x(),              0.0, b01.z());

        f.mQ[0].mP = M * (b01+shift10);
        f.mQ[1].mP = M * (b02+shift10);
        f.mQ[2].mP = M * (b03+shift10);
        exp10 = M * exp10;

        enum predicate pred10;
        Vec3 res10;
        f.findClosestPointToOrigin2Simplex(res10, pred10);
        EXPECT_EQ(pred10, VORONOI_OVER_EDGE_3_1);
        //f.logPredicate(GJKOriginFinder::INFO, pred10);
        EXPECT_EQ(res10,  exp10);

        ////

        Vec3 shift11(rand100(),   0.0, 0.0);
        Vec3 exp11  (shift11.x(), 0.0, b01.z());

        f.mQ[0].mP = M * (b01+shift11);
        f.mQ[1].mP = M * (b02+shift11);
        f.mQ[2].mP = M * (b03+shift11);
        exp11 = M * exp11;

        enum predicate pred11;
        Vec3 res11;
        f.findClosestPointToOrigin2Simplex(res11, pred11);
        EXPECT_EQ(isAnyOf(pred11, VORONOI_OVER_EDGE_3_1,
                                  VORONOI_OVER_VERTEX_1  ), true);
        //f.logPredicate(GJKOriginFinder::INFO, pred11);
        EXPECT_EQ(res11,  exp11);

        ////

        Vec3 shift12(0.0, rand100(),   0.0);
        Vec3 exp12  (0.0, shift12.y(), b01.z());

        f.mQ[0].mP = M * (b01+shift12);
        f.mQ[1].mP = M * (b02+shift12);
        f.mQ[2].mP = M * (b03+shift12);
        exp12 = M * exp12;

        enum predicate pred12;
        Vec3 res12;
        f.findClosestPointToOrigin2Simplex(res12, pred12);
        EXPECT_EQ(isAnyOf(pred12, VORONOI_OVER_EDGE_1_2,
                                  VORONOI_OVER_VERTEX_1  ), true);
        //f.logPredicate(GJKOriginFinder::INFO, pred12);
        EXPECT_EQ(res12,  exp12);

        ////

        Vec3 shift13(rand100(),     rand100(),     0.0);
        Vec3 exp13  (shift13.x(), shift13.y(), b01.z());

        f.mQ[0].mP = M * (b01+shift13);
        f.mQ[1].mP = M * (b02+shift13);
        f.mQ[2].mP = M * (b03+shift13);
        exp13 = M * exp13;

        enum predicate pred13;
        Vec3 res13;
        f.findClosestPointToOrigin2Simplex(res13, pred13);
        EXPECT_EQ(pred13, VORONOI_OVER_VERTEX_1);
        //f.logPredicate(GJKOriginFinder::INFO, pred13);
        EXPECT_EQ(res13,  exp13);

        ////

        Vec3 shift14(-1.0* b02.x(),     rand100(),     0.0);
        Vec3 exp14  (          0.0,   shift14.y(), b01.z());

        f.mQ[0].mP = M * (b01+shift14);
        f.mQ[1].mP = M * (b02+shift14);
        f.mQ[2].mP = M * (b03+shift14);
        exp14 = M * exp14;

        enum predicate pred14;
        Vec3 res14;
        f.findClosestPointToOrigin2Simplex(res14, pred14);
        EXPECT_EQ(isAnyOf(pred14, VORONOI_OVER_EDGE_1_2,
                                  VORONOI_OVER_VERTEX_2  ), true);
        //f.logPredicate(GJKOriginFinder::INFO, pred14);
        EXPECT_EQ(res14,  exp14);

        ////
        double alpha15 = rand100();
        Vec3 shift15(
            -1.0*b02.x() - alpha15 * b03.y(),
                         - alpha15 * b02.x(), 
            0.0
        );
        Vec3 exp15  (
           - alpha15 * b03.y(),
           - alpha15 * b02.x(), 
           b01.z()
        );

        f.mQ[0].mP = M * (b01+shift15);
        f.mQ[1].mP = M * (b02+shift15);
        f.mQ[2].mP = M * (b03+shift15);
        exp15 = M * exp15;

        enum predicate pred15;
        Vec3 res15;
        f.findClosestPointToOrigin2Simplex(res15, pred15);
        EXPECT_EQ(isAnyOf(pred15, VORONOI_OVER_EDGE_2_3,
                                  VORONOI_OVER_VERTEX_2  ), true);
        //f.logPredicate(GJKOriginFinder::INFO, pred15);
        EXPECT_EQ(res15,  exp15);

        ////

        double alpha16_1 = rand100();
        double alpha16_2 = rand100();
        Vec3 shift16(
            -1.0*b02.x() - alpha16_1 * b03.y() - alpha16_2*b02.x(),
                         - alpha16_1 * b02.x() + alpha16_2*b03.y(), 
            0.0
        );
        Vec3 exp16  (
           - alpha16_1 * b03.y() - alpha16_2*b02.x(),
           - alpha16_1 * b02.x() + alpha16_2*b03.y(),
           b01.z()
        );

        f.mQ[0].mP = M * (b01+shift16);
        f.mQ[1].mP = M * (b02+shift16);
        f.mQ[2].mP = M * (b03+shift16);
        exp16 = M * exp16;

        enum predicate pred16;
        Vec3 res16;
        f.findClosestPointToOrigin2Simplex(res16, pred16);
        EXPECT_EQ(pred16, VORONOI_OVER_VERTEX_2);
        //f.logPredicate(GJKOriginFinder::INFO, pred16);
        EXPECT_EQ(res16,  exp16);

        ////

        Vec3 shift17(rand100(), -1.0*b03.y(), 0.0);
        Vec3 exp17  (shift17.x(), 0.0, b01.z());

        f.mQ[0].mP = M * (b01+shift17);
        f.mQ[1].mP = M * (b02+shift17);
        f.mQ[2].mP = M * (b03+shift17);
        exp17 = M * exp17;

        enum predicate pred17;
        Vec3 res17;
        f.findClosestPointToOrigin2Simplex(res17, pred17);
        EXPECT_EQ(isAnyOf(pred17, VORONOI_OVER_EDGE_3_1,
                                  VORONOI_OVER_VERTEX_3  ), true);
        //f.logPredicate(GJKOriginFinder::INFO, pred17);
        EXPECT_EQ(res17,  exp17);

        ////
        double alpha18 = rand100();
        Vec3 shift18(
           -1.0*alpha18*b03.y(),
           -1.0*b03.y() - alpha18*b02.x(),
           0.0
        );
        Vec3 exp18  (
           -1.0*alpha18*b03.y(),
           -1.0*alpha18*b02.x(),
           b01.z()
        );


        f.mQ[0].mP = M * (b01+shift18);
        f.mQ[1].mP = M * (b02+shift18);
        f.mQ[2].mP = M * (b03+shift18);
        exp18 = M * exp18;

        enum predicate pred18;
        Vec3 res18;
        f.findClosestPointToOrigin2Simplex(res18, pred18);
        EXPECT_EQ(isAnyOf(pred18, VORONOI_OVER_EDGE_2_3,
                                  VORONOI_OVER_VERTEX_3  ), true);
        //f.logPredicate(GJKOriginFinder::INFO, pred18);
        EXPECT_EQ(res18,  exp18);

        ////
        double alpha19_1 = rand100();
        double alpha19_2 = rand100();

        Vec3 shift19(
                        -1.0*alpha19_1*b03.y() + alpha19_2,
           -1.0*b03.y() -1.0*alpha19_1*b02.x(),
           0.0
        );
        Vec3 exp19  (
           -1.0*alpha19_1*b03.y() + alpha19_2,
           -1.0*alpha19_1*b02.x(),
           b01.z()
        );


        f.mQ[0].mP = M * (b01+shift19);
        f.mQ[1].mP = M * (b02+shift19);
        f.mQ[2].mP = M * (b03+shift19);
        exp19 = M * exp19;

        enum predicate pred19;
        Vec3 res19;
        f.findClosestPointToOrigin2Simplex(res19, pred19);
        EXPECT_EQ(pred19, VORONOI_OVER_VERTEX_3);
        //f.logPredicate(GJKOriginFinder::INFO, pred19);
        EXPECT_EQ(res19,  exp19);


    }
}




} // namespace Makena

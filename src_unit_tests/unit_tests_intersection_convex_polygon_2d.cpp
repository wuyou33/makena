#include "gtest/gtest.h"
#include <vector>

#include <iostream>
#include <fstream>
#include <sstream>
#include "intersection_convex_polygon_2d.hpp"

namespace Makena {



class IntersectionFinderConvexPolygon2DTests : public ::testing::Test {

  protected:  

    IntersectionFinderConvexPolygon2DTests(){;};
    virtual ~IntersectionFinderConvexPolygon2DTests(){;};
    virtual void SetUp() {;};
    virtual void TearDown() {;};

    IntersectionFinderConvexPolygon2D mFinder;

};


static void createPartialChainsCV(
    Vec2 a1, Vec2 a2, Vec2 a3, 
    Vec2 b1, Vec2 b2, Vec2 b3, 
    CHNode** na, CHNode** nb
) {
    CHNode* n1 = new CHNode(a1, true,  1, -1, CHNode::NT_ORIGINAL_VERTEX);
    CHNode* n2 = new CHNode(a2, true,  2, -1, CHNode::NT_ORIGINAL_VERTEX);
    CHNode* n3 = new CHNode(a3, true,  3, -1, CHNode::NT_ORIGINAL_VERTEX);
    CHNode* n4 = new CHNode(b1, false, 4, -1, CHNode::NT_ORIGINAL_VERTEX);
    CHNode* n5 = new CHNode(b2, false, 5, -1, CHNode::NT_ORIGINAL_VERTEX);
    CHNode* n6 = new CHNode(b3, false, 6, -1, CHNode::NT_ORIGINAL_VERTEX);

    CHEdge* e12 = new CHEdge();
    CHEdge* e23 = new CHEdge();
    CHEdge* e45 = new CHEdge();
    CHEdge* e56 = new CHEdge();

    n1->mNext = e12;
    e12->mSrc = n1;
    e12->mDst = n2;
    n2->mPrev = e12;
    n2->mNext = e23;
    e23->mSrc = n2;
    e23->mDst = n3;
    n3->mPrev = e23;

    n4->mNext = e45;
    e45->mSrc = n4;
    e45->mDst = n5;
    n5->mPrev = e45;
    n5->mNext = e56;
    e56->mSrc = n5;
    e56->mDst = n6;
    n6->mPrev = e56;

    (*na) = n2;
    (*nb) = n5;
}


static void createPartialChainsOE(
      Vec2 a1, Vec2 a2, Vec2 b1, Vec2 b2, Vec2 b3, CHEdge** e, CHNode** nb ) {

    CHNode* n1 = new CHNode(a1, true,  1, -1, CHNode::NT_ORIGINAL_VERTEX);
    CHNode* n2 = new CHNode(a2, true,  2, -1, CHNode::NT_ORIGINAL_VERTEX);
    CHNode* n4 = new CHNode(b1, false,  3,-1, CHNode::NT_ORIGINAL_VERTEX);
    CHNode* n5 = new CHNode(b2, false,  4,-1, CHNode::NT_ORIGINAL_VERTEX);
    CHNode* n6 = new CHNode(b3, false,  5,-1, CHNode::NT_ORIGINAL_VERTEX);

    CHEdge* e12 = new CHEdge();
    CHEdge* e45 = new CHEdge();
    CHEdge* e56 = new CHEdge();

    n1->mNext = e12;
    e12->mSrc = n1;
    e12->mDst = n2;
    n2->mPrev = e12;

    n4->mNext = e45;
    e45->mSrc = n4;
    e45->mDst = n5;
    n5->mPrev = e45;
    n5->mNext = e56;
    e56->mSrc = n5;
    e56->mDst = n6;
    n6->mPrev = e56;

    (*e) =  e12;
    (*nb) = n5;
}



static void createPartialChainsIS(
               Vec2 a1, Vec2 a2, Vec2 b1, Vec2 b2, CHEdge** ea, CHEdge** eb)
{
    CHNode* n1 = new CHNode(a1, true,  1, -1, CHNode::NT_ORIGINAL_VERTEX);
    CHNode* n2 = new CHNode(a2, true,  2, -1, CHNode::NT_ORIGINAL_VERTEX);
    CHNode* n3 = new CHNode(b1, false,  3,-1,  CHNode::NT_ORIGINAL_VERTEX);
    CHNode* n4 = new CHNode(b2, false,  4,-1,  CHNode::NT_ORIGINAL_VERTEX);

    CHEdge* e12 = new CHEdge();
    CHEdge* e34 = new CHEdge();

    n1->mNext = e12;
    e12->mSrc = n1;
    e12->mDst = n2;
    n2->mPrev = e12;

    n3->mNext = e34;
    e34->mSrc = n3;
    e34->mDst = n4;
    n4->mPrev = e34;

    (*ea) =  e12;
    (*eb) =  e34;
}


double TORR = EPSILON_SQUARED * 0.5;
static double pi = 3.1415926;
static Vec2 rotatedVec(double x, double y, double theta)
{
    Vec2 v(cos(theta) * x - sin(theta) * y, sin(theta) * x + cos(theta) * y);
    return v;
}


static bool compareResultsIntsec(
    vector<IntersectionFinderConvexPolygon2D::OutputElem>& v1,
    vector<IntersectionFinderConvexPolygon2D::OutputElem>& v2
) {
    if (v1.size() != v2.size()) {
        return false;
    }
    for (size_t i = 0; i < v1.size() ; i++) {
        auto offset = i;
        bool res = true;
        for (size_t j = 0; j < v1.size() ; j++) {
            auto index1 = j;
            auto index2 = (j+offset)%v1.size();
            Vec2 v = v1[index1].mP - v2[index2].mP;
            if (v.norm2() > 0.01) {
                res = false;
                break;
            }
            if (v1[index1].mType!=v2[index2].mType) {
                res = false;
                break;
            }
            if (((v1[index1].mIndexA    == v2[index2].mIndexA    &&
                  v1[index1].mIndexAaux == v2[index2].mIndexAaux    ) ||
                 (v1[index1].mIndexA    == v2[index2].mIndexAaux &&    
                  v1[index1].mIndexAaux == v2[index2].mIndexA       )   ) &&
                ((v1[index1].mIndexB    == v2[index2].mIndexB    &&
                  v1[index1].mIndexBaux == v2[index2].mIndexBaux    ) ||
                (v1[index1].mIndexB    == v2[index2].mIndexBaux &&  
                  v1[index1].mIndexBaux == v2[index2].mIndexB       )   )   ) {
               ;
            }
            else {
                res = false;
                break;
            }
        }
        if (res) {
            return true;
        }
    }

    return false;
}



/**  @brief CHNode
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test01) {
    Vec2 v01(0.01, 0.01);
    CHNode* n01 = new CHNode(v01, true,  1, -1, CHNode::NT_ORIGINAL_VERTEX);
    Vec2 v02(0.01000001, 0.01000001);
    CHNode* n02 = new CHNode(v02, true,  2, -1, CHNode::NT_ORIGINAL_VERTEX);
    EXPECT_EQ(n01->isCoincident(n02), true);

    Vec2 v03(0.01, 0.02);
    CHNode* n03 = new CHNode(v03, true,  3, -1, CHNode::NT_ORIGINAL_VERTEX);
    EXPECT_EQ(n01->isCoincident(n03), false);

    Vec2 v04(1.2, 3.4);
    CHNode* n04 = new CHNode(v04, true,  4, -1, CHNode::NT_ORIGINAL_VERTEX);
    Vec2 v05(1.2, 3.3);
    CHNode* n05 = new CHNode(v05, true,  5, -1, CHNode::NT_ORIGINAL_VERTEX);

    EXPECT_EQ(n04->isHigherThan(n05), true);
    EXPECT_EQ(n05->isHigherThan(n04), false);

    Vec2 v06(-1.0, -1.0);
    CHNode* n06 = new CHNode(v06, true,  6, -1, CHNode::NT_ORIGINAL_VERTEX);
    Vec2 v07(1.0, 1.0);
    CHNode* n07 = new CHNode(v07, true,  7, -1, CHNode::NT_ORIGINAL_VERTEX);

    CHEdge* e01 = new CHEdge();
    e01->mSrc = n06;
    e01->mDst = n07;
    EXPECT_EQ(n06->isOnEdge(e01), true);
    EXPECT_EQ(n07->isOnEdge(e01), true);
    Vec2 v08(0.0, 0.0);
    CHNode* n08 = new CHNode(v08, true,  8, -1, CHNode::NT_ORIGINAL_VERTEX);
    EXPECT_EQ(n08->isOnEdge(e01), true);    
    Vec2 v09(0.0000001, 0.0);
    CHNode* n09 = new CHNode(v09, true,  9, -1, CHNode::NT_ORIGINAL_VERTEX);
    EXPECT_EQ(n09->isOnEdge(e01), false);    

    Vec2 v10(0.5, -1.0);
    CHNode* n10 = new CHNode(v10, true,  10, -1, CHNode::NT_ORIGINAL_VERTEX);
    Vec2 v11(-0.5, 1.0);
    CHNode* n11 = new CHNode(v11, true,  11, -1, CHNode::NT_ORIGINAL_VERTEX);
    CHEdge* e02 = new CHEdge();
    e02->mSrc = n10;
    e02->mDst = n11;

    Vec2 v12(0.0, 0.0);
    CHNode* n12 = new CHNode(v12, true,  12, -1, CHNode::NT_ORIGINAL_VERTEX);
    EXPECT_EQ(n12->isOnEdge(e02), true);
    Vec2 v13(0.1, -0.200000001);
    CHNode* n13 = new CHNode(v13, true,  13, -1, CHNode::NT_ORIGINAL_VERTEX);
    EXPECT_EQ(n13->isOnEdge(e02), true);
}



/**  @brief CHEdge
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test02) {
    Vec2 v01(0.0, 0.0);
    CHNode* n01 = new CHNode(v01, true,  1, -1, CHNode::NT_ORIGINAL_VERTEX);
    Vec2 v02(0.0, 1.0);
    CHNode* n02 = new CHNode(v02, true,  2, -1, CHNode::NT_ORIGINAL_VERTEX);
    Vec2 v03(1.0, 1.0);
    CHNode* n03 = new CHNode(v03, true,  3, -1, CHNode::NT_ORIGINAL_VERTEX);
    CHEdge* e01 = new CHEdge();
    e01->mSrc = n01;
    e01->mDst = n02;

    CHEdge* e02 = new CHEdge();
    e02->mSrc = n01;
    e02->mDst = n03;

    EXPECT_EQ(e01->isInTheSameDirection(e02), true);

    CHEdge* e03 = new CHEdge();
    e03->mSrc = n02;
    e03->mDst = n01;

    EXPECT_EQ(e01->isInOppositeDirection(e03), true);

    Vec2 v04(1.5, 1.7);
    CHNode* n04 = new CHNode(v04, true,  4, -1, CHNode::NT_ORIGINAL_VERTEX);

    Vec2 v05(3.3, 2.1);
    CHNode* n05 = new CHNode(v05, true,  5, -1, CHNode::NT_ORIGINAL_VERTEX);

    Vec2 v06(-7.5, -0.3);
    CHNode* n06 = new CHNode(v06, true,  6, -1, CHNode::NT_ORIGINAL_VERTEX);

    CHEdge* e04 = new CHEdge();
    e04->mSrc = n04;
    e04->mDst = n05;

    CHEdge* e05 = new CHEdge();
    e05->mSrc = n04;
    e05->mDst = n06;

    EXPECT_EQ(e04->isColinear(e05), true);
    EXPECT_EQ(e04->isColinear(e01), false);

    Vec2 v07(0.0, 0.0);
    CHNode* n07 = new CHNode(v07, true,  7, -1, CHNode::NT_ORIGINAL_VERTEX);
    Vec2 v08(1.0, -1.0);
    CHNode* n08 = new CHNode(v08, true,  8, -1, CHNode::NT_ORIGINAL_VERTEX);

    CHEdge* e06 = new CHEdge();
    e06->mSrc = n07;
    e06->mDst = n08;
    
    Vec2 v09(0.5, 0.0);
    EXPECT_EQ(e06->isCCWAroundSrcTo(v09), true);
    Vec2 v10(0.0, 1.0);
    EXPECT_EQ(e06->isCCWAroundSrcTo(v10), true);
    Vec2 v11(-0.5, 0.0);
    EXPECT_EQ(e06->isCCWAroundSrcTo(v11), false);
    Vec2 v12(0.0, -0.5);
    EXPECT_EQ(e06->isCCWAroundSrcTo(v12), false);
    Vec2 v13(1.0, -1.00001);
    EXPECT_EQ(e06->isCCWAroundSrcTo(v13), false);
    Vec2 v14(1.000001, -1.0);
    EXPECT_EQ(e06->isCCWAroundSrcTo(v14), true);

    Vec2 v15(1.0, 1.0);
    CHNode* n15 = new CHNode(v15, true,  15, -1, CHNode::NT_ORIGINAL_VERTEX);
    Vec2 v16(2.0, 2.0);
    CHNode* n16 = new CHNode(v16, true,  16, -1, CHNode::NT_ORIGINAL_VERTEX);
    Vec2 v17(2.0, 1.0);
    CHNode* n17 = new CHNode(v17, true,  17, -1, CHNode::NT_ORIGINAL_VERTEX);
    Vec2 v18(1.0, 2.0);
    CHNode* n18 = new CHNode(v18, true,  18, -1, CHNode::NT_ORIGINAL_VERTEX);

    CHEdge* e07 = new CHEdge();
    e07->mSrc = n15;
    e07->mDst = n16;
    CHEdge* e08 = new CHEdge();
    e08->mSrc = n17;
    e08->mDst = n18;
    EXPECT_EQ(e07->doesIntersect(e08), true);

    CHEdge* e09 = new CHEdge();
    e09->mSrc = n15;
    e09->mDst = n17;
    CHEdge* e10 = new CHEdge();
    e10->mSrc = n16;
    e10->mDst = n18;
    EXPECT_EQ(e09->doesIntersect(e10), false);

}
        

/**  @brief IntersectionFinderConvexPolygon::areInCCW()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test03) {


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 v01 = rotatedVec(0.0, 0.0, i);
        Vec2 v02 = rotatedVec(0.0, 1.0, i);
        Vec2 v03 = rotatedVec(1.0, 1.0, i);
        Vec2 v04 = rotatedVec(1.0, 1.5, i);

        CHNode* n01 = new CHNode(v01, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n02 = new CHNode(v02, true, 2, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n03 = new CHNode(v03, true, 3, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n04 = new CHNode(v04, true, 4, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e01 = new CHEdge();
        e01->mSrc = n01;
        e01->mDst = n02;

        CHEdge* e02 = new CHEdge();
        e02->mSrc = n01;
        e02->mDst = n03;

        CHEdge* e03 = new CHEdge();
        e03->mSrc = n01;
        e03->mDst = n04;

        CHEdge* e04 = new CHEdge();
        e04->mSrc = n02;
        e04->mDst = n01;


        CHEdge* e05 = new CHEdge();
        e05->mSrc = n03;
        e05->mDst = n01;


        CHEdge* e06 = new CHEdge();
        e06->mSrc = n04;
        e06->mDst = n01;

        EXPECT_EQ(mFinder.areInCCW(
                                    e01, true,  e02, true, e03, true), true);
        EXPECT_EQ(mFinder.areInCCW(
                                    e04, false, e02, true, e03, true), true);
        EXPECT_EQ(mFinder.areInCCW(
                                    e01, true, e05, false, e03, true), true);
        EXPECT_EQ(mFinder.areInCCW(
                                    e01, true, e02, true, e06, false), true);
        EXPECT_EQ(mFinder.areInCCW(
                                    e04, false, e02, true, e06, false), true);
        EXPECT_EQ(mFinder.areInCCW(
                                    e01, true, e03, true, e02, true), false);

    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 v01 = rotatedVec(0.0, 0.0, i);
        Vec2 v02 = rotatedVec(-1.0,0.0, i);
        Vec2 v03 = rotatedVec(1.0, 0.0, i);
        Vec2 v04 = rotatedVec(0.0, 1.0, i);

        CHNode* n01 = new CHNode(v01, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n02 = new CHNode(v02, true, 2, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n03 = new CHNode(v03, true, 3, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n04 = new CHNode(v04, true, 4, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e01 = new CHEdge();
        e01->mSrc = n01;
        e01->mDst = n02;

        CHEdge* e02 = new CHEdge();
        e02->mSrc = n01;
        e02->mDst = n03;

        CHEdge* e03 = new CHEdge();
        e03->mSrc = n01;
        e03->mDst = n04;

        CHEdge* e04 = new CHEdge();
        e04->mSrc = n02;
        e04->mDst = n01;


        CHEdge* e05 = new CHEdge();
        e05->mSrc = n03;
        e05->mDst = n01;


        CHEdge* e06 = new CHEdge();
        e06->mSrc = n04;
        e06->mDst = n01;

        EXPECT_EQ(mFinder.areInCCW(
                                    e01, true, e02, true, e03, true), true);
        EXPECT_EQ(mFinder.areInCCW(
                                    e04, false, e02, true, e03, true), true);
        EXPECT_EQ(mFinder.areInCCW(
                                    e01, true, e05, false, e03, true), true);
        EXPECT_EQ(mFinder.areInCCW(
                                    e01, true, e02, true, e06, false), true);
        EXPECT_EQ(mFinder.areInCCW(
                                    e04, false, e02, true, e06, false), true);
        EXPECT_EQ(mFinder.areInCCW(
                                    e01, true, e03, true, e02, true), false);

    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 v01 = rotatedVec( 0.0,  0.0, i);
        Vec2 v02 = rotatedVec( 0.0,  1.0, i);
        Vec2 v03 = rotatedVec(-1.0, -1.0, i);
        Vec2 v04 = rotatedVec( 1.0, -1.0, i);

        CHNode* n01 = new CHNode(v01, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n02 = new CHNode(v02, true, 2, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n03 = new CHNode(v03, true, 3, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n04 = new CHNode(v04, true, 4, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e01 = new CHEdge();
        e01->mSrc = n01;
        e01->mDst = n02;

        CHEdge* e02 = new CHEdge();
        e02->mSrc = n01;
        e02->mDst = n03;

        CHEdge* e03 = new CHEdge();
        e03->mSrc = n01;
        e03->mDst = n04;

        CHEdge* e04 = new CHEdge();
        e04->mSrc = n02;
        e04->mDst = n01;


        CHEdge* e05 = new CHEdge();
        e05->mSrc = n03;
        e05->mDst = n01;


        CHEdge* e06 = new CHEdge();
        e06->mSrc = n04;
        e06->mDst = n01;

        EXPECT_EQ(mFinder.areInCCW(
                                     e01, true, e02, true, e03, true), true);
        EXPECT_EQ(mFinder.areInCCW(
                                     e04, false, e02, true, e03, true), true);
        EXPECT_EQ(mFinder.areInCCW(
                                     e01, true, e05, false, e03, true), true);
        EXPECT_EQ(mFinder.areInCCW(
                                     e01, true, e02, true, e06, false), true);
        EXPECT_EQ(mFinder.areInCCW(
                                     e04, false, e02, true, e06, false), true);
        EXPECT_EQ(mFinder.areInCCW(
                                     e01, true, e03, true, e02, true), false);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 v01 = rotatedVec( 0.0,  0.0, i);
        Vec2 v02 = rotatedVec( 0.0,  1.0, i);
        Vec2 v03 = rotatedVec(-1.0,  1.0, i);
        Vec2 v04 = rotatedVec(-1.0, -1.0, i);

        CHNode* n01 = new CHNode(v01, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n02 = new CHNode(v02, true, 2, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n03 = new CHNode(v03, true, 3, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n04 = new CHNode(v04, true, 4, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e01 = new CHEdge();
        e01->mSrc = n01;
        e01->mDst = n02;

        CHEdge* e02 = new CHEdge();
        e02->mSrc = n01;
        e02->mDst = n03;

        CHEdge* e03 = new CHEdge();
        e03->mSrc = n01;
        e03->mDst = n04;

        CHEdge* e04 = new CHEdge();
        e04->mSrc = n02;
        e04->mDst = n01;


        CHEdge* e05 = new CHEdge();
        e05->mSrc = n03;
        e05->mDst = n01;


        CHEdge* e06 = new CHEdge();
        e06->mSrc = n04;
        e06->mDst = n01;

        EXPECT_EQ(mFinder.areInCCW(
                                     e01, true, e02, true, e03, true), true);
        EXPECT_EQ(mFinder.areInCCW(
                                     e04, false, e02, true, e03, true), true);
        EXPECT_EQ(mFinder.areInCCW(
                                     e01, true, e05, false, e03, true), true);
        EXPECT_EQ(mFinder.areInCCW(
                                     e01, true, e02, true, e06, false), true);
        EXPECT_EQ(mFinder.areInCCW(
                                     e04, false, e02, true, e06, false), true);
        EXPECT_EQ(mFinder.areInCCW(
                                     e01, true, e03, true, e02, true), false);

    }

}


/**  @brief IntersectionFinderConvexPolygon2D::removeDegeneracy()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test04) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 v01  = rotatedVec(0.0, 0.0, i);
        Vec2 v02  = rotatedVec(1.0, 0.0, i);
        Vec2 v021 = rotatedVec(1.0+TORR, 0.0-TORR, i);
        Vec2 v022 = rotatedVec(1.0-TORR, 0.0+TORR, i);
        Vec2 v03  = rotatedVec(2.0, 1.0, i);
        Vec2 v04  = rotatedVec(2.0, 2.0, i);
        Vec2 v041 = rotatedVec(2.0+TORR, 2.0, i);
        Vec2 v045 = rotatedVec(1.5, 2.5, i);
        Vec2 v046 = rotatedVec(1.5+TORR, 2.5-TORR, i);
        Vec2 v05  = rotatedVec(1.0, 3.0, i);
        Vec2 v06  = rotatedVec(-1.0, 3.0, i);
        Vec2 v061 = rotatedVec(-1.3, 2.7, i);
        Vec2 v062 = rotatedVec(-1.6, 2.4, i);
        Vec2 v063 = rotatedVec(-1.7, 2.3, i);
        Vec2 v07  = rotatedVec(-2.0, 2.0, i);
        Vec2 v08  = rotatedVec(-2.0, 1.0, i);
        Vec2 v09  = rotatedVec(-1.0, 0.0, i);
        Vec2 v091 = rotatedVec(-1.0+TORR, 0.0, i);
        Vec2 v095 = rotatedVec(-0.8, 0.0, i);
        Vec2 v096 = rotatedVec(-0.6, 0.0, i);
        Vec2 v097 = rotatedVec(-0.4, 0.0, i);
        Vec2 v098 = rotatedVec(-0.2, 0.0, i);


        IntersectionFinderConvexPolygon2D::InputElem n01 (v01,   1);
        IntersectionFinderConvexPolygon2D::InputElem n02 (v02,   2);
        IntersectionFinderConvexPolygon2D::InputElem n021(v021, 21);
        IntersectionFinderConvexPolygon2D::InputElem n022(v022, 22);
        IntersectionFinderConvexPolygon2D::InputElem n03 (v03,   3);
        IntersectionFinderConvexPolygon2D::InputElem n04 (v04,   4);
        IntersectionFinderConvexPolygon2D::InputElem n041(v041, 41);
        IntersectionFinderConvexPolygon2D::InputElem n045(v045, 45);
        IntersectionFinderConvexPolygon2D::InputElem n046(v046, 46);
        IntersectionFinderConvexPolygon2D::InputElem n05 (v05,   5);
        IntersectionFinderConvexPolygon2D::InputElem n06 (v06,   6);
        IntersectionFinderConvexPolygon2D::InputElem n061(v061, 61);
        IntersectionFinderConvexPolygon2D::InputElem n062(v062, 62);
        IntersectionFinderConvexPolygon2D::InputElem n063(v063, 63);
        IntersectionFinderConvexPolygon2D::InputElem n07 (v07,   7);
        IntersectionFinderConvexPolygon2D::InputElem n08 (v08,   8);
        IntersectionFinderConvexPolygon2D::InputElem n09 (v09,   9);
        IntersectionFinderConvexPolygon2D::InputElem n091(v091, 91);
        IntersectionFinderConvexPolygon2D::InputElem n095(v095, 95);
        IntersectionFinderConvexPolygon2D::InputElem n096(v096, 96);
        IntersectionFinderConvexPolygon2D::InputElem n097(v097, 97);
        IntersectionFinderConvexPolygon2D::InputElem n098(v098, 98);

        vector<IntersectionFinderConvexPolygon2D::InputElem> vec01; 

        vec01.push_back(n01);
        vec01.push_back(n02);
        vec01.push_back(n021);
        vec01.push_back(n022);
        vec01.push_back(n03);
        vec01.push_back(n04);
        vec01.push_back(n045);
        vec01.push_back(n046);
        vec01.push_back(n05);
        vec01.push_back(n06);
        vec01.push_back(n061);
        vec01.push_back(n062);
        vec01.push_back(n063);
        vec01.push_back(n07);
        vec01.push_back(n08);
        vec01.push_back(n09);
        vec01.push_back(n091);
        vec01.push_back(n095);
        vec01.push_back(n096);
        vec01.push_back(n098);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> vec02; 
        mFinder.removeDegeneracy(vec01, vec02);

        EXPECT_EQ(vec02.size(), 8);

    }
}



/**  @brief IntersectionFinderConvexPolygon2D::removeDegeneracy()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test04out) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 v01  = rotatedVec(0.0, 0.0, i);
        Vec2 v02  = rotatedVec(1.0, 0.0, i);
        Vec2 v021 = rotatedVec(1.0+TORR, 0.0-TORR, i);
        Vec2 v022 = rotatedVec(1.0-TORR, 0.0+TORR, i);
        Vec2 v03  = rotatedVec(2.0, 1.0, i);
        Vec2 v04  = rotatedVec(2.0, 2.0, i);
        Vec2 v041 = rotatedVec(2.0+TORR, 2.0, i);
        Vec2 v045 = rotatedVec(1.5, 2.5, i);
        Vec2 v046 = rotatedVec(1.5+TORR, 2.5-TORR, i);
        Vec2 v05  = rotatedVec(1.0, 3.0, i);
        Vec2 v06  = rotatedVec(-1.0, 3.0, i);
        Vec2 v061 = rotatedVec(-1.3, 2.7, i);
        Vec2 v062 = rotatedVec(-1.6, 2.4, i);
        Vec2 v063 = rotatedVec(-1.7, 2.3, i);
        Vec2 v07  = rotatedVec(-2.0, 2.0, i);
        Vec2 v08  = rotatedVec(-2.0, 1.0, i);
        Vec2 v09  = rotatedVec(-1.0, 0.0, i);
        Vec2 v091 = rotatedVec(-1.0+TORR, 0.0, i);
        Vec2 v095 = rotatedVec(-0.8, 0.0, i);
        Vec2 v096 = rotatedVec(-0.6, 0.0, i);
        Vec2 v097 = rotatedVec(-0.4, 0.0, i);
        Vec2 v098 = rotatedVec(-0.2, 0.0, i);


        IntersectionFinderConvexPolygon2D::OutputElem n01 (
        v01, 1, -1, 1, -1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n02 (
        v02,   2,-1, 2,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n021(
        v021, 21,-1, 21,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n022(
        v022, 22, -1,22,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n03 (
        v03,   3,-1, 3,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n04 (
        v04,   4,-1, 4,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n041(
        v041, 41,-1, 41,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n045(
        v045, 45,-1, 45,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n046(
        v046, 46,-1, 46,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n05 (
        v05,   5,-1, 5,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n06 (
        v06,   6,-1, 6,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n061(
        v061, 61,-1, 61,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n062(
        v062, 62,-1, 62,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n063(
        v063, 63,-1, 63,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n07 (
        v07,   7,-1, 7,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n08 (
        v08,   8,-1, 8,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n09 (
        v09,   9,-1, 9,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n091(
        v091, 91, -1,91,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n095(
        v095, 95,-1, 95,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n096(
        v096, 96,-1, 96,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n097(
        v097, 97,-1, 97,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);
        IntersectionFinderConvexPolygon2D::OutputElem n098(
        v098, 98,-1, 98,-1, IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> vec01; 

        vec01.push_back(n01);
        vec01.push_back(n02);
        vec01.push_back(n021);
        vec01.push_back(n022);
        vec01.push_back(n03);
        vec01.push_back(n04);
        vec01.push_back(n045);
        vec01.push_back(n046);
        vec01.push_back(n05);
        vec01.push_back(n06);
        vec01.push_back(n061);
        vec01.push_back(n062);
        vec01.push_back(n063);
        vec01.push_back(n07);
        vec01.push_back(n08);
        vec01.push_back(n09);
        vec01.push_back(n091);
        vec01.push_back(n095);
        vec01.push_back(n096);
        vec01.push_back(n098);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> vec02; 
        mFinder.removeDegeneracy(vec01, vec02);

        EXPECT_EQ(vec02.size(), 8);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::
 *          findPerpendicularProjectionOfPointOntoEdge()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test05) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 v05((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 v01(0.0, 0.0);
        Vec2 v02(2.0, 1.0);
        Vec2 v03(1.0, 1.0);
        Vec2 v04(1.2, 0.6);

        v01 += v05;
        v02 += v05;
        v03 += v05;
        v04 += v05;

        v01  = rotatedVec(v01.x(), v01.y(), i);
        v02  = rotatedVec(v02.x(), v02.y(), i);
        v03  = rotatedVec(v03.x(), v03.y(), i);
        v04  = rotatedVec(v04.x(), v04.y(), i);

        vector<Vec2> pts;
        pts.push_back(v01);
        pts.push_back(v02);
        double s;
        Vec2 v06 = mFinder.
                findPerpendicularProjectionOfPointOntoEdge(v03, v01, v02, s);
        //cerr << setprecision(20) << "v04: " << v04.x() << "," 
        //                                    << v04.y() << "\n";
        //cerr << setprecision(20) << "v06: " << v06.x() << "," 
        //                                    << v06.y() << "\n";
        EXPECT_EQ(v04==v06, true);

        Vec2 v07 = mFinder.
                findPerpendicularProjectionOfPointOntoEdge(v04, v01, v02, s);
        EXPECT_EQ(v04==v07, true);
    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection_edge_edge()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test06) {

    //////////////////////
    // COINCIDENT CASES //
    //////////////////////

    // a1 == b1
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {
        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 1.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(1.3, 0.6);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a1, true);
        EXPECT_EQ(output_01[0].mIndexA==1, true);
        EXPECT_EQ(output_01[0].mIndexB==3, true);
        EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
    }


    // a1 == b1 colinear opposite
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 1.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(-2.0, -1.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a1, true);
        EXPECT_EQ(output_01[0].mIndexA==1, true);
        EXPECT_EQ(output_01[0].mIndexB==3, true);
        EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

    }


    // a1 == b1 colinear b2
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 1.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(1.2, 0.6);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a1&&output_01[1].mP==b2)||
                  (output_01[1].mP==a1&&output_01[0].mP==b2),   true);

        if (output_01[0].mP==a1) {
            EXPECT_EQ(output_01[0].mIndexA==1, true);
            EXPECT_EQ(output_01[0].mIndexB==3, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(
                   (output_01[1].mIndexA==1&&output_01[1].mIndexAaux==2)||
                   (output_01[1].mIndexA==2&&output_01[1].mIndexAaux==1),true);
            EXPECT_EQ(output_01[1].mIndexB==4, true);
            EXPECT_EQ(output_01[1].mType==
                        IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==1, true);
            EXPECT_EQ(output_01[1].mIndexB==3, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
            EXPECT_EQ(output_01[0].mIndexB==4, true);
            EXPECT_EQ(output_01[0].mType==
                        IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
    }

    // a1 == b1 colinear a2
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 1.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(3.0, 1.5);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);

        EXPECT_EQ((output_01[0].mP==a1&&output_01[1].mP==a2)||
                  (output_01[1].mP==a1&&output_01[0].mP==a2),   true);
        if (output_01[0].mP==a1) {
            EXPECT_EQ(output_01[0].mIndexA==1, true);
            EXPECT_EQ(output_01[0].mIndexB==3, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[1].mIndexA==2, true);
            EXPECT_EQ(
                   (output_01[1].mIndexB==3&&output_01[1].mIndexBaux==4)||
                   (output_01[1].mIndexB==4&&output_01[1].mIndexBaux==3),true);
            EXPECT_EQ(output_01[1].mType==
                        IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==1, true);
            EXPECT_EQ(output_01[1].mIndexB==3, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[0].mIndexA==2, true);
            EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);
            EXPECT_EQ(output_01[0].mType==
                        IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);
        }

    }

    // a1 == b1 colinear a2 == b2
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 1.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(2.0, 1.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a1&&output_01[1].mP==a2)||
                  (output_01[1].mP==a1&&output_01[0].mP==a2),   true);

        if (output_01[0].mP==a1) {
            EXPECT_EQ(output_01[0].mIndexA==1, true);
            EXPECT_EQ(output_01[0].mIndexB==3, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[1].mIndexA==2, true);
            EXPECT_EQ(output_01[1].mIndexB==4, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==1, true);
            EXPECT_EQ(output_01[1].mIndexB==3, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[0].mIndexA==2, true);
            EXPECT_EQ(output_01[0].mIndexB==4, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
        }
    }

    // a1 == b2
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 1.0);
        Vec2 b1(1.3, 0.6);
        Vec2 b2(0.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a1, true);
        EXPECT_EQ(output_01[0].mIndexA==1, true);
        EXPECT_EQ(output_01[0].mIndexB==4, true);
        EXPECT_EQ(output_01[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
    }


    // a1 == b2 colinear opposite
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 1.0);
        Vec2 b1(-2.0, -1.0);
        Vec2 b2(0.0, 0.0);


        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a1, true);
        EXPECT_EQ(output_01[0].mIndexA==1, true);
        EXPECT_EQ(output_01[0].mIndexB==4, true);
        EXPECT_EQ(output_01[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

    }


    // a1 == b2 colinear b1
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 1.0);
        Vec2 b1(1.2, 0.6);
        Vec2 b2(0.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a1&&output_01[1].mP==b1)||
                  (output_01[1].mP==a1&&output_01[0].mP==b1),   true);

        if (output_01[0].mP==a1) {
            EXPECT_EQ(output_01[0].mIndexA==1, true);
            EXPECT_EQ(output_01[0].mIndexB==4, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(
                   (output_01[1].mIndexA==1&&output_01[1].mIndexAaux==2)||
                   (output_01[1].mIndexA==2&&output_01[1].mIndexAaux==1),true);
            EXPECT_EQ(output_01[1].mIndexB==3  , true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==1, true);
            EXPECT_EQ(output_01[1].mIndexB==4, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
            EXPECT_EQ(output_01[0].mIndexB==3   , true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
    }

    // a1 == b2 colinear a2
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 1.0);
        Vec2 b1(3.0, 1.5);
        Vec2 b2(0.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a1&&output_01[1].mP==a2)||
                  (output_01[1].mP==a1&&output_01[0].mP==a2),   true);

        if (output_01[0].mP==a1) {
            EXPECT_EQ(output_01[0].mIndexA==1, true);
            EXPECT_EQ(output_01[0].mIndexB==4, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[1].mIndexA==2, true);
            EXPECT_EQ(
                   (output_01[1].mIndexB==3&&output_01[1].mIndexBaux==4)||
                   (output_01[1].mIndexB==4&&output_01[1].mIndexBaux==3),true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==1, true);
            EXPECT_EQ(output_01[1].mIndexB==4, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[0].mIndexA==2 , true);
            EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);
        }
    }


    // a1 == b2 colinear a2 == b1
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 1.0);
        Vec2 b1(2.0, 1.0);
        Vec2 b2(0.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a1&&output_01[1].mP==a2)||
                  (output_01[1].mP==a1&&output_01[0].mP==a2),   true);

        if (output_01[0].mP==a1) {
            EXPECT_EQ(output_01[0].mIndexA==1, true);
            EXPECT_EQ(output_01[0].mIndexB==4, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[1].mIndexA==2, true);
            EXPECT_EQ(output_01[1].mIndexB==3, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==1, true);
            EXPECT_EQ(output_01[1].mIndexB==4, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[0].mIndexA==2, true);
            EXPECT_EQ(output_01[0].mIndexB==3, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
        }
    }


    // a2 == b1
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 1.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(1.3, 0.6);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a2, true);
        EXPECT_EQ(output_01[0].mIndexA==2, true);
        EXPECT_EQ(output_01[0].mIndexB==3, true);
        EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

    }


    // a2 == b1 colinear opposite
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 1.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(-2.0, -1.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a2, true);
        EXPECT_EQ(output_01[0].mIndexA==2, true);
        EXPECT_EQ(output_01[0].mIndexB==3, true);
        EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

    }


    // a2 == b1 colinear b2
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 1.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(1.2, 0.6);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a2&&output_01[1].mP==b2)||
                  (output_01[1].mP==a2&&output_01[0].mP==b2),   true);

        if (output_01[0].mP==a2) {
            EXPECT_EQ(output_01[0].mIndexA==2, true);
            EXPECT_EQ(output_01[0].mIndexB==3, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(
                   (output_01[1].mIndexA==1&&output_01[1].mIndexAaux==2)||
                   (output_01[1].mIndexA==2&&output_01[1].mIndexAaux==1),true);
            EXPECT_EQ(output_01[1].mIndexB==4, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==2, true);
            EXPECT_EQ(output_01[1].mIndexB==3, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
            EXPECT_EQ(output_01[0].mIndexB==4, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
    }

    // a2 == b1 colinear a1
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 1.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(3.0, 1.5);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a2&&output_01[1].mP==a1)||
                  (output_01[1].mP==a2&&output_01[0].mP==a1),   true);

        if (output_01[0].mP==a2) {
            EXPECT_EQ(output_01[0].mIndexA==2, true);
            EXPECT_EQ(output_01[0].mIndexB==3, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[1].mIndexA==1  , true);
            EXPECT_EQ(
                   (output_01[1].mIndexB==3&&output_01[1].mIndexBaux==4)||
                   (output_01[1].mIndexB==4&&output_01[1].mIndexBaux==3),true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==2, true);
            EXPECT_EQ(output_01[1].mIndexB==3, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[0].mIndexA==1  , true);
            EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);
        }
    }


    // a2 == b1 colinear a1 == b2
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 1.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(2.0, 1.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a2&&output_01[1].mP==a1)||
                  (output_01[1].mP==a2&&output_01[0].mP==a1),   true);

        if (output_01[0].mP==a2) {
            EXPECT_EQ(output_01[0].mIndexA==2, true);
            EXPECT_EQ(output_01[0].mIndexB==3, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[1].mIndexA==1  , true);
            EXPECT_EQ(output_01[1].mIndexB==4  , true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==2, true);
            EXPECT_EQ(output_01[1].mIndexB==3, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[0].mIndexA==1  , true);
            EXPECT_EQ(output_01[0].mIndexB==4  , true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
        }
    }

    // a2 == b2
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 1.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(1.3, 0.6);
        Vec2 b2(0.0, 0.0);


        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a2, true);
        EXPECT_EQ(output_01[0].mIndexA==2, true);
        EXPECT_EQ(output_01[0].mIndexB==4, true);
        EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

    }

    // a2 == b2 colinear opposite
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 1.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(-2.0, -1.0);
        Vec2 b2(0.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a2, true);
        EXPECT_EQ(output_01[0].mIndexA==2, true);
        EXPECT_EQ(output_01[0].mIndexB==4, true);
        EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
    }


    // a2 == b2 colinear b1
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 1.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(1.2, 0.6);
        Vec2 b2(0.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a2&&output_01[1].mP==b1)||
                  (output_01[1].mP==a2&&output_01[0].mP==b1),   true);

        if (output_01[0].mP==a2) {
            EXPECT_EQ(output_01[0].mIndexA==2, true);
            EXPECT_EQ(output_01[0].mIndexB==4, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(
                   (output_01[1].mIndexA==1&&output_01[1].mIndexAaux==2)||
                   (output_01[1].mIndexA==2&&output_01[1].mIndexAaux==1),true);
            EXPECT_EQ(output_01[1].mIndexB==3  , true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==2, true);
            EXPECT_EQ(output_01[1].mIndexB==4, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
            EXPECT_EQ(output_01[0].mIndexB==3  , true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
    }


    // a2 == b2 colinear a1
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 1.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(3.0, 1.5);
        Vec2 b2(0.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a2&&output_01[1].mP==a1)||
                  (output_01[1].mP==a2&&output_01[0].mP==a1),   true);

        if (output_01[0].mP==a2) {
            EXPECT_EQ(output_01[0].mIndexA==2, true);
            EXPECT_EQ(output_01[0].mIndexB==4, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[1].mIndexA==1, true);
            EXPECT_EQ(
                   (output_01[1].mIndexB==3&&output_01[1].mIndexBaux==4)||
                   (output_01[1].mIndexB==4&&output_01[1].mIndexBaux==3),true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==2, true);
            EXPECT_EQ(output_01[1].mIndexB==4, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[0].mIndexA==1, true);
            EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);
        }
    }

    // a2 == b2 colinear a1 == b1
    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 1.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(2.0, 1.0);
        Vec2 b2(0.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a2&&output_01[1].mP==a1)||
                  (output_01[1].mP==a2&&output_01[0].mP==a1),   true);

        if (output_01[0].mP==a2) {
            EXPECT_EQ(output_01[0].mIndexA==2, true);
            EXPECT_EQ(output_01[0].mIndexB==4, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[1].mIndexA==1, true);
            EXPECT_EQ(output_01[1].mIndexB==3, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==2, true);
            EXPECT_EQ(output_01[1].mIndexB==4, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
            EXPECT_EQ(output_01[0].mIndexA==1, true);
            EXPECT_EQ(output_01[0].mIndexB==3, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
        }
    }

    //////////////////////
    //  COLINEAR CASES  //
    //////////////////////

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-2.0, 0.0);
        Vec2 a2(-1.0, 0.0);
        Vec2 b1(1.0, 0.0);
        Vec2 b2(2.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(output_01.size(), 0);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-2.0, 0.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(2.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a2, true);
        EXPECT_EQ(output_01[0].mIndexA==2, true);
        EXPECT_EQ(output_01[0].mIndexB==3, true);
        EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);
    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-2.0, 0.0);
        Vec2 a2(1.0, 0.0);
        Vec2 b1(-1.0, 0.0);
        Vec2 b2(2.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a2&&output_01[1].mP==b1)||
                  (output_01[1].mP==a2&&output_01[0].mP==b1),   true);

        if (output_01[0].mP==a2) {
            EXPECT_EQ(output_01[0].mIndexA==2, true);
            EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);
            EXPECT_EQ(
                   (output_01[1].mIndexA==1&&output_01[1].mIndexAaux==2)||
                   (output_01[1].mIndexA==2&&output_01[1].mIndexAaux==1),true);
            EXPECT_EQ(output_01[1].mIndexB==3, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==2, true);
            EXPECT_EQ(
                   (output_01[1].mIndexB==3&&output_01[1].mIndexBaux==4)||
                   (output_01[1].mIndexB==4&&output_01[1].mIndexBaux==3),true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);
            EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
            EXPECT_EQ(output_01[0].mIndexB==3, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-2.0, 0.0);
        Vec2 a2(2.0, 0.0);
        Vec2 b1(-1.0, 0.0);
        Vec2 b2(1.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==b1&&output_01[1].mP==b2)||
                  (output_01[1].mP==b1&&output_01[0].mP==b2),   true);

        if (output_01[0].mP==b1) {
            EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
            EXPECT_EQ(output_01[0].mIndexB==3, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
            EXPECT_EQ(
                   (output_01[1].mIndexA==1&&output_01[1].mIndexAaux==2)||
                   (output_01[1].mIndexA==2&&output_01[1].mIndexAaux==1),true);
            EXPECT_EQ(output_01[1].mIndexB==4, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
        else {
            EXPECT_EQ(
                   (output_01[1].mIndexA==1&&output_01[1].mIndexAaux==2)||
                   (output_01[1].mIndexA==2&&output_01[1].mIndexAaux==1),true);
            EXPECT_EQ(output_01[1].mIndexB==3, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
            EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
            EXPECT_EQ(output_01[0].mIndexB==4, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-1.0, 0.0);
        Vec2 a2(2.0, 0.0);
        Vec2 b1(-2.0, 0.0);
        Vec2 b2(1.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 2);
        EXPECT_EQ((output_01[0].mP==a1&&output_01[1].mP==b2)||
                  (output_01[1].mP==a1&&output_01[0].mP==b2),   true);

        if (output_01[0].mP==a1) {
            EXPECT_EQ(output_01[0].mIndexA==1, true);
            EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);
            EXPECT_EQ(
                   (output_01[1].mIndexA==1&&output_01[1].mIndexAaux==2)||
                   (output_01[1].mIndexA==2&&output_01[1].mIndexAaux==1),true);
            EXPECT_EQ(output_01[1].mIndexB==4, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
        else {
            EXPECT_EQ(output_01[1].mIndexA==1, true);
            EXPECT_EQ(
                   (output_01[1].mIndexB==3&&output_01[1].mIndexBaux==4)||
                   (output_01[1].mIndexB==4&&output_01[1].mIndexBaux==3),true);
            EXPECT_EQ(output_01[1].mIndexB==3, true);
            EXPECT_EQ(output_01[1].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);
            EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
            EXPECT_EQ(output_01[0].mIndexB==4, true);
            EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);
        }
    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 0.0);
        Vec2 b1(-2.0, 0.0);
        Vec2 b2(0.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a1, true);

        EXPECT_EQ(output_01[0].mIndexA==1, true);
        EXPECT_EQ(output_01[0].mIndexB==4, true);
        EXPECT_EQ(output_01[0].mType==
                      IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(1.0, 0.0);
        Vec2 a2(2.0, 0.0);
        Vec2 b1(-2.0, 0.0);
        Vec2 b2(-1.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(output_01.size(), 0);

    }



    ///////////////////////////
    // PARALLEL NOT COLINEAR // 
    ///////////////////////////

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-2.0, 0.0);
        Vec2 a2(-1.0, 0.0);
        Vec2 b1(1.0, 0.5);
        Vec2 b2(2.0, 0.5);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(output_01.size(), 0);

    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-2.0, 0.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(0.0, 0.5);
        Vec2 b2(2.0, 0.5);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(output_01.size(), 0);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-2.0, 0.0);
        Vec2 a2(1.0, 0.0);
        Vec2 b1(-1.0, 0.5);
        Vec2 b2(2.0, 0.5);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(output_01.size(), 0);

    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-2.0, 0.0);
        Vec2 a2(2.0, 0.0);
        Vec2 b1(-1.0, 0.5);
        Vec2 b2(1.0, 0.5);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(output_01.size(), 0);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-1.0, 0.0);
        Vec2 a2(1.0, 0.0);
        Vec2 b1(-2.0, 0.5);
        Vec2 b2(2.0, 0.5);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(output_01.size(), 0);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-1.0, 0.0);
        Vec2 a2(2.0, 0.0);
        Vec2 b1(-2.0, 0.5);
        Vec2 b2(1.0, 0.5);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(output_01.size(), 0);

    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 0.0);
        Vec2 b1(-2.0, 0.5);
        Vec2 b2(0.0, 0.5);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(output_01.size(), 0);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(1.0, 0.0);
        Vec2 a2(2.0, 0.0);
        Vec2 b1(-2.0, 0.5);
        Vec2 b2(-1.0, 0.5);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(output_01.size(), 0);

    }


    ////////////////////
    // VERTEX ON EDGE //
    ////////////////////

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-1.0, 1.0);
        Vec2 a2(0.5, 0.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(2.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a2, true);

        EXPECT_EQ(output_01[0].mIndexA==2, true);
        EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);
        EXPECT_EQ(output_01[0].mType==
                       IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.5, 0.0);
        Vec2 a2(-1.0, 1.0);
        Vec2 b1(0.0, 0.0);
        Vec2 b2(2.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a1, true);

        EXPECT_EQ(output_01[0].mIndexA==1, true);
        EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);

        EXPECT_EQ(output_01[0].mType==
                       IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);

    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-1.0, 1.0);
        Vec2 a2(0.5, 0.0);
        Vec2 b1(2.0, 0.0);
        Vec2 b2(0.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==a2, true);

        EXPECT_EQ(output_01[0].mIndexA==2, true);
        EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);
        EXPECT_EQ(output_01[0].mType==
                       IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 0.0);
        Vec2 b1(-1.0, 1.0);
        Vec2 b2(0.5, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==b2, true);

        EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
        EXPECT_EQ(output_01[0].mIndexB==4, true);
        EXPECT_EQ(output_01[0].mType==
                       IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 0.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(-1.0, 1.0);
        Vec2 b2(0.5, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==b2, true);

        EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
        EXPECT_EQ(output_01[0].mIndexB==4, true);
        EXPECT_EQ(output_01[0].mType==
                       IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 0.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(0.5, 0.0);
        Vec2 b2(-1.0, 1.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==b1, true);

        EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);

        EXPECT_EQ(output_01[0].mIndexB==3, true);
        EXPECT_EQ(output_01[0].mType==
                       IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);

    }


    /////////////////////
    // PROPER CROSSING //
    /////////////////////

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(0.0, 0.0);
        Vec2 a2(2.0, 1.0);
        Vec2 b1(1.0, 1.0);
        Vec2 b2(2.0, 0.5);
        Vec2 v1(1.5, 0.75);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;
        v1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        v1  = rotatedVec(v1.x(), v1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==v1, true);

        EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
        EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);
        EXPECT_EQ(output_01[0].mType==
                       IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE, true);

    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 1.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(1.0, 1.0);
        Vec2 b2(2.0, 0.5);
        Vec2 v1(1.5, 0.75);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;
        v1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        v1  = rotatedVec(v1.x(), v1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==v1, true);

        EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
        EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);
        EXPECT_EQ(output_01[0].mType==
                       IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE, true);

    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 1.0);
        Vec2 a2(0.0, 0.0);
        Vec2 b1(2.0, 0.5);
        Vec2 b2(1.0, 1.0);
        Vec2 v1(1.5, 0.75);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;
        v1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        v1  = rotatedVec(v1.x(), v1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==v1, true);

        EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
        EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);

        EXPECT_EQ(output_01[0].mType==
                       IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE, true);

    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(2.0, 0.5);
        Vec2 a2(1.0, 1.0);
        Vec2 b1(2.0, 1.0);
        Vec2 b2(0.0, 0.0);
        Vec2 v1(1.5, 0.75);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;
        v1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        v1  = rotatedVec(v1.x(), v1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_a1);
        input_a01.push_back(&ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_b1);
        input_b01.push_back(&ie_b2);
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        bool res_01 = mFinder.
                  findIntersection_edge_edge(input_a01, input_b01, output_01);

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size(), 1);
        EXPECT_EQ(output_01[0].mP==v1, true);

        EXPECT_EQ(
                   (output_01[0].mIndexA==1&&output_01[0].mIndexAaux==2)||
                   (output_01[0].mIndexA==2&&output_01[0].mIndexAaux==1),true);
        EXPECT_EQ(
                   (output_01[0].mIndexB==3&&output_01[0].mIndexBaux==4)||
                   (output_01[0].mIndexB==4&&output_01[0].mIndexBaux==3),true);

        EXPECT_EQ(output_01[0].mType==
                       IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE, true);
    }
}


/**  @brief IntersectionFinderConvexPolygon2D::
 *          findIntersection_point_polygon()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test07) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(0.0, -1.0);
        Vec2 p2(3.0, -1.0);
        Vec2 p3(4.0, -0.6);
        Vec2 p4(4.0,  1.0);
        Vec2 p5(1.5,  3.0);
        Vec2 p6(0.0,  3.0);
        Vec2 p7(-1.0, 2.0);
        Vec2 p8(-1.0, 0.0);

        Vec2 t01(0.0, 0.0);
        Vec2 t02(-0.5, -0.5);
        Vec2 t03(2.0, -1.0);
        Vec2 t04(4.0, 0.0);
        Vec2 t05(-1.0, 0.0001);
        Vec2 t06(1.0, 1.0);
        Vec2 t07(5.0, 1.0);
        Vec2 t08(4.0, 1.0001);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;
        p7 += trans;
        p8 += trans;

        t01 += trans;
        t02 += trans;
        t03 += trans;
        t04 += trans;
        t05 += trans;
        t06 += trans;
        t07 += trans;
        t08 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);
        p7  = rotatedVec(p7.x(), p7.y(), i);
        p8  = rotatedVec(p8.x(), p8.y(), i);

        t01  = rotatedVec(t01.x(), t01.y(), i);
        t02  = rotatedVec(t02.x(), t02.y(), i);
        t03  = rotatedVec(t03.x(), t03.y(), i);
        t04  = rotatedVec(t04.x(), t04.y(), i);
        t05  = rotatedVec(t05.x(), t05.y(), i);
        t06  = rotatedVec(t06.x(), t06.y(), i);
        t07  = rotatedVec(t07.x(), t07.y(), i);
        t08  = rotatedVec(t08.x(), t08.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_p1(p1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_p2(p2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_p3(p3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_p4(p4, 4);
        IntersectionFinderConvexPolygon2D::InputElem ie_p5(p5, 5);
        IntersectionFinderConvexPolygon2D::InputElem ie_p6(p6, 6);
        IntersectionFinderConvexPolygon2D::InputElem ie_p7(p7, 7);
        IntersectionFinderConvexPolygon2D::InputElem ie_p8(p8, 8);

        IntersectionFinderConvexPolygon2D::InputElem ie_t01(t01, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_t02(t02, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_t03(t03, 13);
        IntersectionFinderConvexPolygon2D::InputElem ie_t04(t04, 14);
        IntersectionFinderConvexPolygon2D::InputElem ie_t05(t05, 15);
        IntersectionFinderConvexPolygon2D::InputElem ie_t06(t06, 16);
        IntersectionFinderConvexPolygon2D::InputElem ie_t07(t07, 17);
        IntersectionFinderConvexPolygon2D::InputElem ie_t08(t08, 18);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
        input_a01.push_back(&ie_t01);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a02;
        input_a02.push_back(&ie_t02);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a03;
        input_a03.push_back(&ie_t03);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a04;
        input_a04.push_back(&ie_t04);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a05;
        input_a05.push_back(&ie_t05);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a06;
        input_a06.push_back(&ie_t06);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a07;
        input_a07.push_back(&ie_t07);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a08;
        input_a08.push_back(&ie_t08);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a09;
        input_a09.push_back(&ie_p1);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a10;
        input_a10.push_back(&ie_p2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a11;
        input_a11.push_back(&ie_p3);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a12;
        input_a12.push_back(&ie_p4);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a13;
        input_a13.push_back(&ie_p5);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a14;
        input_a14.push_back(&ie_p6);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a15;
        input_a15.push_back(&ie_p7);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a16;
        input_a16.push_back(&ie_p8);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
        input_b01.push_back(&ie_p1);
        input_b01.push_back(&ie_p2);
        input_b01.push_back(&ie_p3);
        input_b01.push_back(&ie_p4);
        input_b01.push_back(&ie_p5);
        input_b01.push_back(&ie_p6);
        input_b01.push_back(&ie_p7);
        input_b01.push_back(&ie_p8);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_03;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_04;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_05;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_06;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_07;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_08;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_09;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_10;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_11;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_12;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_13;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_14;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_15;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_16;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_17;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_18;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_19;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_20;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_21;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_22;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_23;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_24;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_25;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_26;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_27;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_28;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_29;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_30;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_31;
        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_32;

        bool res_01 = mFinder.
               findIntersection_point_polygon(input_a01, input_b01, output_01);
                                               
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(output_01.size()==1, true);
        EXPECT_EQ(output_01[0].mP==t01, true);
        EXPECT_EQ(output_01[0].mIndexA==11, true);
        EXPECT_EQ(output_01[0].mIndexB==-1, true);
        EXPECT_EQ(output_01[0].mType==
                   IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR, true);


        bool res_02 = mFinder.
               findIntersection_polygon_point(input_b01, input_a01, output_02);
                                               
        EXPECT_EQ(res_02, true);
        EXPECT_EQ(output_02.size()==1, true);
        EXPECT_EQ(output_02[0].mP==t01, true);
        EXPECT_EQ(output_02[0].mIndexA==-1, true);
        EXPECT_EQ(output_02[0].mIndexB==11, true);
        EXPECT_EQ(output_02[0].mType==
                    IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX, true);

        bool res_03 = mFinder.
               findIntersection_point_polygon(input_a02, input_b01, output_03);
                                               
        EXPECT_EQ(res_03, true);
        EXPECT_EQ(output_03.size()==1, true);
        EXPECT_EQ(output_03[0].mP==t02, true);
        EXPECT_EQ(output_03[0].mIndexA==12, true);
        EXPECT_EQ((output_03[0].mIndexB==1 && output_03[0].mIndexBaux==8)||
                  (output_03[0].mIndexB==8 && output_03[0].mIndexBaux==1)
                  , true);
        EXPECT_EQ(output_03[0].mType==
                        IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);

        bool res_04 = mFinder.
               findIntersection_polygon_point(input_b01, input_a02, output_04);
                                               
        EXPECT_EQ(res_04, true);
        EXPECT_EQ(output_04.size()==1, true);
        EXPECT_EQ(output_04[0].mP==t02, true);
        EXPECT_EQ((output_04[0].mIndexA==1 && output_04[0].mIndexAaux==8)||
                  (output_04[0].mIndexA==8 && output_04[0].mIndexAaux==1)
                  , true);
        EXPECT_EQ(output_04[0].mIndexB==12, true);
        EXPECT_EQ(output_04[0].mType==
                        IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);

        bool res_05 = mFinder.
               findIntersection_point_polygon(input_a03, input_b01, output_05);
                                               
        EXPECT_EQ(res_05, true);
        EXPECT_EQ(output_05.size()==1, true);
        EXPECT_EQ(output_05[0].mP==t03, true);
        EXPECT_EQ(output_05[0].mIndexA==13, true);
        EXPECT_EQ((output_05[0].mIndexB==1 && output_05[0].mIndexBaux==2)||
                  (output_05[0].mIndexB==2 && output_05[0].mIndexBaux==1)
                  , true);
        EXPECT_EQ(output_05[0].mType==
                        IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);

        bool res_06 = mFinder.
               findIntersection_polygon_point(input_b01, input_a03, output_06);
                                               
        EXPECT_EQ(res_06, true);
        EXPECT_EQ(output_06.size()==1, true);
        EXPECT_EQ(output_06[0].mP==t03, true);
        EXPECT_EQ((output_06[0].mIndexA==1 && output_06[0].mIndexAaux==2)||
                  (output_06[0].mIndexA==2 && output_06[0].mIndexAaux==1)
                  , true);
        EXPECT_EQ(output_06[0].mIndexB==13, true);
        EXPECT_EQ(output_06[0].mType==
                        IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);

        bool res_07 = mFinder.
               findIntersection_point_polygon(input_a04, input_b01, output_07);
                                               
        EXPECT_EQ(res_07, true);
        EXPECT_EQ(output_07.size()==1, true);
        EXPECT_EQ(output_07[0].mP==t04, true);
        EXPECT_EQ(output_07[0].mIndexA==14, true);
        EXPECT_EQ((output_07[0].mIndexB==3 && output_07[0].mIndexBaux==4)||
                  (output_07[0].mIndexB==4 && output_07[0].mIndexBaux==3)
                  , true);
        EXPECT_EQ(output_07[0].mType==
                        IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);

        bool res_08 = mFinder.
               findIntersection_polygon_point(input_b01, input_a04, output_08);
                                               
        EXPECT_EQ(res_08, true);
        EXPECT_EQ(output_08.size()==1, true);
        EXPECT_EQ(output_08[0].mP==t04, true);
        EXPECT_EQ((output_08[0].mIndexA==3 && output_08[0].mIndexAaux==4)||
                  (output_08[0].mIndexA==4 && output_08[0].mIndexAaux==3)
                  , true);
        EXPECT_EQ(output_08[0].mIndexB==14, true);
        EXPECT_EQ(output_08[0].mType==
                        IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);

        bool res_09 = mFinder.
               findIntersection_point_polygon(input_a05, input_b01, output_09);
                                               
        EXPECT_EQ(res_09, true);
        EXPECT_EQ(output_09.size()==1, true);
        EXPECT_EQ(output_09[0].mP==t05, true);
        EXPECT_EQ(output_09[0].mIndexA==15, true);
        EXPECT_EQ((output_09[0].mIndexB==7 && output_09[0].mIndexBaux==8)||
                  (output_09[0].mIndexB==8 && output_09[0].mIndexBaux==7)
                  , true);
        EXPECT_EQ(output_09[0].mType==
                        IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE, true);

        bool res_10 = mFinder.
               findIntersection_polygon_point(input_b01, input_a05, output_10);
                                               
        EXPECT_EQ(res_10, true);
        EXPECT_EQ(output_10.size()==1, true);
        EXPECT_EQ(output_10[0].mP==t05, true);
        EXPECT_EQ((output_10[0].mIndexA==7 && output_10[0].mIndexAaux==8)||
                  (output_10[0].mIndexA==8 && output_10[0].mIndexAaux==7)
                  , true);
        EXPECT_EQ(output_10[0].mIndexB==15, true);
        EXPECT_EQ(output_10[0].mType==
                        IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX, true);

        bool res_11 = mFinder.
               findIntersection_point_polygon(input_a06, input_b01, output_11);
                                               
        EXPECT_EQ(res_11, true);
        EXPECT_EQ(output_11.size()==1, true);
        EXPECT_EQ(output_11[0].mP==t06, true);
        EXPECT_EQ(output_11[0].mIndexA==16, true);
        EXPECT_EQ(output_11[0].mIndexB==-1, true);
        EXPECT_EQ(output_11[0].mType==
                    IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR, true);

        bool res_12 = mFinder.
               findIntersection_polygon_point(input_b01, input_a06, output_12);
                                               
        EXPECT_EQ(res_12, true);
        EXPECT_EQ(output_12.size()==1, true);
        EXPECT_EQ(output_12[0].mP==t06, true);
        EXPECT_EQ(output_12[0].mIndexA==-1, true);
        EXPECT_EQ(output_12[0].mIndexB==16, true);
        EXPECT_EQ(output_12[0].mType==
                   IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX, true);

        bool res_13 = mFinder.
               findIntersection_point_polygon(input_a07, input_b01, output_13);
                                               
        EXPECT_EQ(res_13, false);
        EXPECT_EQ(output_13.size()==0, true);

        bool res_14 = mFinder.
               findIntersection_polygon_point(input_b01, input_a07, output_14);
                                               
        EXPECT_EQ(res_14, false);
        EXPECT_EQ(output_14.size()==0, true);

        bool res_15 = mFinder.
               findIntersection_point_polygon(input_a08, input_b01, output_15);
                                               
        EXPECT_EQ(res_15, false);
        EXPECT_EQ(output_15.size()==0, true);

        bool res_16 = mFinder.
               findIntersection_polygon_point(input_b01, input_a08, output_16);
                                               
        EXPECT_EQ(res_16, false);
        EXPECT_EQ(output_16.size()==0, true);

        bool res_17 = mFinder.
               findIntersection_point_polygon(input_a09, input_b01, output_17);
                                               
        EXPECT_EQ(res_17, true);
        EXPECT_EQ(output_17.size()==1, true);
        EXPECT_EQ(output_17[0].mP==p1, true);
        EXPECT_EQ(output_17[0].mIndexA==1, true);
        EXPECT_EQ(output_17[0].mIndexB==1, true);
        EXPECT_EQ(output_17[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_18 = mFinder.
               findIntersection_polygon_point(input_b01, input_a09, output_18);
                                               
        EXPECT_EQ(res_18, true);
        EXPECT_EQ(output_18.size()==1, true);
        EXPECT_EQ(output_18[0].mP==p1, true);
        EXPECT_EQ(output_18[0].mIndexA==1, true);
        EXPECT_EQ(output_18[0].mIndexB==1, true);
        EXPECT_EQ(output_18[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_19 = mFinder.
               findIntersection_point_polygon(input_a10, input_b01, output_19);
                                               
        EXPECT_EQ(res_19, true);
        EXPECT_EQ(output_19.size()==1, true);
        EXPECT_EQ(output_19[0].mP==p2, true);
        EXPECT_EQ(output_19[0].mIndexA==2, true);
        EXPECT_EQ(output_19[0].mIndexB==2, true);
        EXPECT_EQ(output_19[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_20 = mFinder.
               findIntersection_polygon_point(input_b01, input_a10, output_20);
                                               
        EXPECT_EQ(res_20, true);
        EXPECT_EQ(output_20.size()==1, true);
        EXPECT_EQ(output_20[0].mP==p2, true);
        EXPECT_EQ(output_20[0].mIndexA==2, true);
        EXPECT_EQ(output_20[0].mIndexB==2, true);
        EXPECT_EQ(output_20[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_21 = mFinder.
               findIntersection_point_polygon(input_a11, input_b01, output_21);
                                               
        EXPECT_EQ(res_21, true);
        EXPECT_EQ(output_21.size()==1, true);
        EXPECT_EQ(output_21[0].mP==p3, true);
        EXPECT_EQ(output_21[0].mIndexA==3, true);
        EXPECT_EQ(output_21[0].mIndexB==3, true);
        EXPECT_EQ(output_21[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_22 = mFinder.
               findIntersection_polygon_point(input_b01, input_a11, output_22);
                                               
        EXPECT_EQ(res_22, true);
        EXPECT_EQ(output_22.size()==1, true);
        EXPECT_EQ(output_22[0].mP==p3, true);
        EXPECT_EQ(output_22[0].mIndexA==3, true);
        EXPECT_EQ(output_22[0].mIndexB==3, true);
        EXPECT_EQ(output_22[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_23 = mFinder.
               findIntersection_point_polygon(input_a12, input_b01, output_23);
                                               
        EXPECT_EQ(res_23, true);
        EXPECT_EQ(output_23.size()==1, true);
        EXPECT_EQ(output_23[0].mP==p4, true);
        EXPECT_EQ(output_23[0].mIndexA==4, true);
        EXPECT_EQ(output_23[0].mIndexB==4, true);
        EXPECT_EQ(output_23[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_24 = mFinder.
               findIntersection_polygon_point(input_b01, input_a12, output_24);
                                               
        EXPECT_EQ(res_24, true);
        EXPECT_EQ(output_24.size()==1, true);
        EXPECT_EQ(output_24[0].mP==p4, true);
        EXPECT_EQ(output_24[0].mIndexA==4, true);
        EXPECT_EQ(output_24[0].mIndexB==4, true);
        EXPECT_EQ(output_24[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_25 = mFinder.
               findIntersection_point_polygon(input_a13, input_b01, output_25);
                                               
        EXPECT_EQ(res_25, true);
        EXPECT_EQ(output_25.size()==1, true);
        EXPECT_EQ(output_25[0].mP==p5, true);
        EXPECT_EQ(output_25[0].mIndexA==5, true);
        EXPECT_EQ(output_25[0].mIndexB==5, true);
        EXPECT_EQ(output_25[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_26 = mFinder.
               findIntersection_polygon_point(input_b01, input_a13, output_26);
                                               
        EXPECT_EQ(res_26, true);
        EXPECT_EQ(output_26.size()==1, true);
        EXPECT_EQ(output_26[0].mP==p5, true);
        EXPECT_EQ(output_26[0].mIndexA==5, true);
        EXPECT_EQ(output_26[0].mIndexB==5, true);
        EXPECT_EQ(output_26[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_27 = mFinder.
               findIntersection_point_polygon(input_a14, input_b01, output_27);
                                               
        EXPECT_EQ(res_27, true);
        EXPECT_EQ(output_27.size()==1, true);
        EXPECT_EQ(output_27[0].mP==p6, true);
        EXPECT_EQ(output_27[0].mIndexA==6, true);
        EXPECT_EQ(output_27[0].mIndexB==6, true);
        EXPECT_EQ(output_27[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_28 = mFinder.
               findIntersection_polygon_point(input_b01, input_a14, output_28);
                                               
        EXPECT_EQ(res_28, true);
        EXPECT_EQ(output_28.size()==1, true);
        EXPECT_EQ(output_28[0].mP==p6, true);
        EXPECT_EQ(output_28[0].mIndexA==6, true);
        EXPECT_EQ(output_28[0].mIndexB==6, true);
        EXPECT_EQ(output_28[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_29 = mFinder.
               findIntersection_point_polygon(input_a15, input_b01, output_29);
                                               
        EXPECT_EQ(res_29, true);
        EXPECT_EQ(output_29.size()==1, true);
        EXPECT_EQ(output_29[0].mP==p7, true);
        EXPECT_EQ(output_29[0].mIndexA==7, true);
        EXPECT_EQ(output_29[0].mIndexB==7, true);
        EXPECT_EQ(output_29[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_30 = mFinder.
               findIntersection_polygon_point(input_b01, input_a15, output_30);
                                               
        EXPECT_EQ(res_30, true);
        EXPECT_EQ(output_30.size()==1, true);
        EXPECT_EQ(output_30[0].mP==p7, true);
        EXPECT_EQ(output_30[0].mIndexA==7, true);
        EXPECT_EQ(output_30[0].mIndexB==7, true);
        EXPECT_EQ(output_30[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_31 = mFinder.
               findIntersection_point_polygon(input_a16, input_b01, output_31);
                                               
        EXPECT_EQ(res_31, true);
        EXPECT_EQ(output_31.size()==1, true);
        EXPECT_EQ(output_31[0].mP==p8, true);
        EXPECT_EQ(output_31[0].mIndexA==8, true);
        EXPECT_EQ(output_31[0].mIndexB==8, true);
        EXPECT_EQ(output_31[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

        bool res_32 = mFinder.
               findIntersection_polygon_point(input_b01, input_a16, output_32);
                                               
        EXPECT_EQ(res_32, true);
        EXPECT_EQ(output_32.size()==1, true);
        EXPECT_EQ(output_32[0].mP==p8, true);
        EXPECT_EQ(output_32[0].mIndexA==8, true);
        EXPECT_EQ(output_32[0].mIndexB==8, true);
        EXPECT_EQ(output_32[0].mType==
                     IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX, true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::constructCycle()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test08) {

    {

        Vec2 p1(0.0, -1.0);
        Vec2 p2(3.0, -1.0);
        Vec2 p3(4.0, -0.6);
        Vec2 p4(4.0,  1.0);
        Vec2 p5(1.5,  3.0);
        Vec2 p6(0.0,  3.0);
        Vec2 p7(-1.0, 2.0);
        Vec2 p8(-1.0, 0.0);

        IntersectionFinderConvexPolygon2D::InputElem ie_p1(p1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_p2(p2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_p3(p3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_p4(p4, 4);
        IntersectionFinderConvexPolygon2D::InputElem ie_p5(p5, 5);
        IntersectionFinderConvexPolygon2D::InputElem ie_p6(p6, 6);
        IntersectionFinderConvexPolygon2D::InputElem ie_p7(p7, 7);
        IntersectionFinderConvexPolygon2D::InputElem ie_p8(p8, 8);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> pts_01;
        pts_01.push_back(&ie_p1);
        pts_01.push_back(&ie_p2);
        pts_01.push_back(&ie_p3);
        pts_01.push_back(&ie_p4);
        pts_01.push_back(&ie_p5);
        pts_01.push_back(&ie_p6);
        pts_01.push_back(&ie_p7);
        pts_01.push_back(&ie_p8);

        CHNode* n1;
        CHNode* n2;
        long index01, index02;
        mFinder.constructCycle(pts_01,&n1,&n2,true, index01, index02);

        EXPECT_EQ(index01, 5);
        EXPECT_EQ(index02, 1);

        EXPECT_EQ(n1->mP == p6, true);
        EXPECT_EQ(n2->mP == p2, true);

        CHNode* n_06 = n1;
        EXPECT_EQ(n_06->mP == p6, true);
        CHEdge* e_06_07 = n_06->mNext;
        EXPECT_EQ(e_06_07->mSrc, n_06);
        CHNode* n_07 = e_06_07->mDst;
        EXPECT_EQ(n_07->mP == p7, true);
        CHEdge* e_07_08 = n_07->mNext;
        EXPECT_EQ(e_07_08->mSrc, n_07);
        CHNode* n_08 = e_07_08->mDst;
        EXPECT_EQ(n_08->mP == p8, true);
        CHEdge* e_08_01 = n_08->mNext;
        EXPECT_EQ(e_08_01->mSrc, n_08);
        CHNode* n_01 = e_08_01->mDst;
        EXPECT_EQ(n_01->mP == p1, true);
        CHEdge* e_01_02 = n_01->mNext;
        EXPECT_EQ(e_01_02->mSrc, n_01);
        CHNode* n_02 = e_01_02->mDst;
        EXPECT_EQ(n_02->mP == p2, true);
        CHEdge* e_02_03 = n_02->mNext;
        EXPECT_EQ(e_02_03->mSrc, n_02);
        CHNode* n_03 = e_02_03->mDst;
        EXPECT_EQ(n_03->mP == p3, true);
        CHEdge* e_03_04 = n_03->mNext;
        EXPECT_EQ(e_03_04->mSrc, n_03);
        CHNode* n_04 = e_03_04->mDst;
        EXPECT_EQ(n_04->mP == p4, true);
        CHEdge* e_04_05 = n_04->mNext;
        EXPECT_EQ(e_04_05->mSrc, n_04);
        CHNode* n_05 = e_04_05->mDst;
        EXPECT_EQ(n_05->mP == p5, true);
        CHEdge* e_05_06 = n_05->mNext;
        EXPECT_EQ(e_05_06->mSrc, n_05);
        CHNode* n_06_2 = e_05_06->mDst;
        EXPECT_EQ(n_06_2->mP == p6, true);
        EXPECT_EQ(n_06, n_06_2);

        EXPECT_EQ(e_06_07->mDst, n_07);
        EXPECT_EQ(e_07_08->mDst, n_08);
        EXPECT_EQ(e_08_01->mDst, n_01);
        EXPECT_EQ(e_01_02->mDst, n_02);
        EXPECT_EQ(e_02_03->mDst, n_03);
        EXPECT_EQ(e_03_04->mDst, n_04);
        EXPECT_EQ(e_04_05->mDst, n_05);
        EXPECT_EQ(e_05_06->mDst, n_06);

        EXPECT_EQ(n_06->mPrev, e_05_06);
        EXPECT_EQ(n_07->mPrev, e_06_07);
        EXPECT_EQ(n_08->mPrev, e_07_08);
        EXPECT_EQ(n_01->mPrev, e_08_01);
        EXPECT_EQ(n_02->mPrev, e_01_02);
        EXPECT_EQ(n_03->mPrev, e_02_03);
        EXPECT_EQ(n_04->mPrev, e_03_04);
        EXPECT_EQ(n_05->mPrev, e_04_05);

        EXPECT_EQ(n_06->mNext, e_06_07);
        EXPECT_EQ(n_07->mNext, e_07_08);
        EXPECT_EQ(n_08->mNext, e_08_01);
        EXPECT_EQ(n_01->mNext, e_01_02);
        EXPECT_EQ(n_02->mNext, e_02_03);
        EXPECT_EQ(n_03->mNext, e_03_04);
        EXPECT_EQ(n_04->mNext, e_04_05);
        EXPECT_EQ(n_05->mNext, e_05_06);

        EXPECT_EQ(n_06->mPeer, nullptr);
        EXPECT_EQ(n_07->mPeer, nullptr);
        EXPECT_EQ(n_08->mPeer, nullptr);
        EXPECT_EQ(n_01->mPeer, nullptr);
        EXPECT_EQ(n_02->mPeer, nullptr);
        EXPECT_EQ(n_03->mPeer, nullptr);
        EXPECT_EQ(n_04->mPeer, nullptr);
        EXPECT_EQ(n_05->mPeer, nullptr);

        EXPECT_EQ(e_06_07->mSrcInterior, false);
        EXPECT_EQ(e_07_08->mSrcInterior, false);
        EXPECT_EQ(e_08_01->mSrcInterior, false);
        EXPECT_EQ(e_01_02->mSrcInterior, false);
        EXPECT_EQ(e_02_03->mSrcInterior, false);
        EXPECT_EQ(e_03_04->mSrcInterior, false);
        EXPECT_EQ(e_04_05->mSrcInterior, false);
        EXPECT_EQ(e_05_06->mSrcInterior, false);

        EXPECT_EQ(e_06_07->mDstInterior, false);
        EXPECT_EQ(e_07_08->mDstInterior, false);
        EXPECT_EQ(e_08_01->mDstInterior, false);
        EXPECT_EQ(e_01_02->mDstInterior, false);
        EXPECT_EQ(e_02_03->mDstInterior, false);
        EXPECT_EQ(e_03_04->mDstInterior, false);
        EXPECT_EQ(e_04_05->mDstInterior, false);
        EXPECT_EQ(e_05_06->mDstInterior, false);

        EXPECT_EQ(n_01->mType, CHNode::NT_ORIGINAL_VERTEX);
        EXPECT_EQ(n_02->mType, CHNode::NT_ORIGINAL_VERTEX);
        EXPECT_EQ(n_03->mType, CHNode::NT_ORIGINAL_VERTEX);
        EXPECT_EQ(n_04->mType, CHNode::NT_ORIGINAL_VERTEX);
        EXPECT_EQ(n_05->mType, CHNode::NT_ORIGINAL_VERTEX);
        EXPECT_EQ(n_06->mType, CHNode::NT_ORIGINAL_VERTEX);
        EXPECT_EQ(n_07->mType, CHNode::NT_ORIGINAL_VERTEX);
        EXPECT_EQ(n_08->mType, CHNode::NT_ORIGINAL_VERTEX);

        EXPECT_EQ(n_01->mIsTouchingPoint, false);
        EXPECT_EQ(n_02->mIsTouchingPoint, false);
        EXPECT_EQ(n_03->mIsTouchingPoint, false);
        EXPECT_EQ(n_04->mIsTouchingPoint, false);
        EXPECT_EQ(n_05->mIsTouchingPoint, false);
        EXPECT_EQ(n_06->mIsTouchingPoint, false);
        EXPECT_EQ(n_07->mIsTouchingPoint, false);
        EXPECT_EQ(n_08->mIsTouchingPoint, false);

        EXPECT_EQ(n_01->mThisIsA, true);
        EXPECT_EQ(n_02->mThisIsA, true);
        EXPECT_EQ(n_03->mThisIsA, true);
        EXPECT_EQ(n_04->mThisIsA, true);
        EXPECT_EQ(n_05->mThisIsA, true);
        EXPECT_EQ(n_06->mThisIsA, true);
        EXPECT_EQ(n_07->mThisIsA, true);
        EXPECT_EQ(n_08->mThisIsA, true);

        EXPECT_EQ(n_01->mIndex, 1);
        EXPECT_EQ(n_02->mIndex, 2);
        EXPECT_EQ(n_03->mIndex, 3);
        EXPECT_EQ(n_04->mIndex, 4);
        EXPECT_EQ(n_05->mIndex, 5);
        EXPECT_EQ(n_06->mIndex, 6);
        EXPECT_EQ(n_07->mIndex, 7);
        EXPECT_EQ(n_08->mIndex, 8);

        EXPECT_EQ(n_01->mIndexAux, -1);
        EXPECT_EQ(n_02->mIndexAux, -1);
        EXPECT_EQ(n_03->mIndexAux, -1);
        EXPECT_EQ(n_04->mIndexAux, -1);
        EXPECT_EQ(n_05->mIndexAux, -1);
        EXPECT_EQ(n_06->mIndexAux, -1);
        EXPECT_EQ(n_07->mIndexAux, -1);
        EXPECT_EQ(n_08->mIndexAux, -1);

    }
}

/**  @brief IntersectionFinderConvexPolygon2D::constructCycle()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test09) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(0.0, -1.0);
        Vec2 p2(3.0, -1.0);


        IntersectionFinderConvexPolygon2D::InputElem ie_p1(p1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_p2(p2, 2);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> pts_01;
        pts_01.push_back(&ie_p1);
        pts_01.push_back(&ie_p2);

        CHNode* n1;
        CHNode* n2;
        long index_01, index_02;
        mFinder.constructCycle(pts_01,&n1,&n2,false, index_01, index_02);

        EXPECT_EQ(index_01, 0);
        EXPECT_EQ(index_02, 1);
        EXPECT_EQ(n1->mP == p1, true);
        EXPECT_EQ(n2->mP == p2, true);

        CHNode* n_01 = n1;
        EXPECT_EQ(n_01->mP == p1, true);
        CHEdge* e_01_02 = n_01->mNext;
        EXPECT_EQ(e_01_02->mSrc, n_01);
        CHNode* n_02 = e_01_02->mDst;
        EXPECT_EQ(n_02->mP == p2, true);
        CHEdge* e_02_01 = n_02->mNext;
        EXPECT_EQ(e_02_01->mSrc, n_02);
        CHNode* n_01_2 = e_02_01->mDst;
        EXPECT_EQ(n_01, n_01_2);

        EXPECT_EQ(e_01_02->mDst, n_02);
        EXPECT_EQ(n_01->mPrev, e_02_01);
        EXPECT_EQ(n_01->mNext, e_01_02);
        EXPECT_EQ(n_02->mPrev, e_01_02);
        EXPECT_EQ(n_02->mNext, e_02_01);
        EXPECT_EQ(n_01->mPeer, nullptr);
        EXPECT_EQ(n_02->mPeer, nullptr);
        EXPECT_EQ(e_01_02->mSrcInterior, false);
        EXPECT_EQ(e_02_01->mSrcInterior, false);
        EXPECT_EQ(e_01_02->mDstInterior, false);
        EXPECT_EQ(e_02_01->mDstInterior, false);


        EXPECT_EQ(n_01->mType, CHNode::NT_ORIGINAL_VERTEX);
        EXPECT_EQ(n_02->mType, CHNode::NT_ORIGINAL_VERTEX);

        EXPECT_EQ(n_01->mIsTouchingPoint, false);
        EXPECT_EQ(n_02->mIsTouchingPoint, false);

        EXPECT_EQ(n_01->mThisIsA, false);
        EXPECT_EQ(n_02->mThisIsA, false);

        EXPECT_EQ(n_01->mIndex, 1);
        EXPECT_EQ(n_02->mIndex, 2);


        EXPECT_EQ(n_01->mIndexAux, -1);
        EXPECT_EQ(n_02->mIndexAux, -1);
    }
}


/**  @brief IntersectionFinderConvexPolygon2D::destructCycle()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test10) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(0.0, -1.0);
        Vec2 p2(3.0, -1.0);
        Vec2 p3(4.0, -0.6);
        Vec2 p4(4.0,  1.0);
        Vec2 p5(1.5,  3.0);
        Vec2 p6(0.0,  3.0);
        Vec2 p7(-1.0, 2.0);
        Vec2 p8(-1.0, 0.0);

        IntersectionFinderConvexPolygon2D::InputElem ie_p1(p1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_p2(p2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_p3(p3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_p4(p4, 4);
        IntersectionFinderConvexPolygon2D::InputElem ie_p5(p5, 5);
        IntersectionFinderConvexPolygon2D::InputElem ie_p6(p6, 6);
        IntersectionFinderConvexPolygon2D::InputElem ie_p7(p7, 7);
        IntersectionFinderConvexPolygon2D::InputElem ie_p8(p8, 8);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> pts_01;
        pts_01.push_back(&ie_p1);
        pts_01.push_back(&ie_p2);
        pts_01.push_back(&ie_p3);
        pts_01.push_back(&ie_p4);
        pts_01.push_back(&ie_p5);
        pts_01.push_back(&ie_p6);
        pts_01.push_back(&ie_p7);
        pts_01.push_back(&ie_p8);

        CHNode* n1;
        CHNode* n2;
        long index_01, index_02;
        mFinder.constructCycle(pts_01,&n1,&n2,true, index_01, index_02);

        mFinder.destructCycle(n1);

    }
}


/**  @brief static Vec2 IntersectionFinderConvexPolygon2D::
 *                      findIntersectionOfTwoLines()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test11) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(-1.0, 0.0);
        Vec2 p2(1.0,  0.0);
        Vec2 p3(0.0, -1.0);
        Vec2 p4(0.0,  1.0);

        Vec2 p5(0.0,  0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
       
        Vec2 n1 = mFinder.
                              findIntersectionOfTwoLines(p1, p2-p1, p3, p4-p3);
        EXPECT_EQ(n1,p5);

    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(0.0,  0.0);
        Vec2 p2(4.0,  1.0);
        Vec2 p3(3.0,  0.0);
        Vec2 p4(4.0,  2.0);
        Vec2 p5(24.0/7.0, 6.0/7.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        Vec2 n1 = mFinder.
                             findIntersectionOfTwoLines(p1, p2-p1, p3, p4-p3);
        EXPECT_EQ(n1, p5);
    }
}


/**  @brief IntersectionFinderConvexPolygon2D::splitEdge() 1
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test12) {


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(0.0, -1.0);
        Vec2 p2(3.0, -1.0);
        Vec2 p3(4.0, -0.6);
        Vec2 p4(4.0,  1.0);
        Vec2 p5(1.5,  3.0);
        Vec2 p6(0.0,  3.0);
        Vec2 p7(-1.0, 2.0);
        Vec2 p8(-1.0, 0.0);

        Vec2 p9(4.5, 0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;
        p7 += trans;
        p8 += trans;
        p9 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);
        p7  = rotatedVec(p7.x(), p7.y(), i);
        p8  = rotatedVec(p8.x(), p8.y(), i);
        p9  = rotatedVec(p9.x(), p9.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_p1(p1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_p2(p2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_p3(p3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_p4(p4, 4);
        IntersectionFinderConvexPolygon2D::InputElem ie_p5(p5, 5);
        IntersectionFinderConvexPolygon2D::InputElem ie_p6(p6, 6);
        IntersectionFinderConvexPolygon2D::InputElem ie_p7(p7, 7);
        IntersectionFinderConvexPolygon2D::InputElem ie_p8(p8, 8);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> pts_01;
        pts_01.push_back(&ie_p1);
        pts_01.push_back(&ie_p2);
        pts_01.push_back(&ie_p3);
        pts_01.push_back(&ie_p4);
        pts_01.push_back(&ie_p5);
        pts_01.push_back(&ie_p6);
        pts_01.push_back(&ie_p7);
        pts_01.push_back(&ie_p8);

        CHNode* n1;
        CHNode* n2;
        long index_01, index_02;
        mFinder.constructCycle(pts_01, &n1, &n2,true, index_01, index_02);
        CHEdge* e;
        CHNode* cur;
        for (cur = n1; cur->mP != p3; cur=cur->mNext->mDst) {;}

        e = cur->mNext;

        CHNode* n3 = mFinder.splitEdge(
                                           e, p9, true,CHNode::NT_SPLIT_POINT);
        CHEdge* en = n3->mNext;
        CHEdge* ep = n3->mPrev;

        EXPECT_EQ(n3->mP,               p9);
        EXPECT_EQ(n3->mPeer,            nullptr);
        EXPECT_EQ(n3->mType,            CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(n3->mIsTouchingPoint, false);
        EXPECT_EQ(n3->mTop,             false);
        EXPECT_EQ(n3->mBottom,          false);
        EXPECT_EQ(n3->mThisIsA,         true);
        EXPECT_EQ(n3->mIndex,           3);
        EXPECT_EQ(n3->mIndexAux,        4);

        EXPECT_EQ(en->mDst->mP,    p4);
        EXPECT_EQ(en->mDst->mPrev, en);
        EXPECT_EQ(en->mSrc,        n3);

        EXPECT_EQ(ep->mSrc->mP,    p3);
        EXPECT_EQ(ep->mSrc->mNext, ep);
        EXPECT_EQ(ep->mDst,        n3);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::splitEdge() 2
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test13) {


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(0.0, -1.0);
        Vec2 p2(3.0, -1.0);
        Vec2 p3(4.0, -0.6);
        Vec2 p4(4.0,  1.0);
        Vec2 p5(1.5,  3.0);
        Vec2 p6(0.0,  3.0);
        Vec2 p7(-1.0, 2.0);
        Vec2 p8(-1.0, 0.0);

        Vec2 p9(4.5, 0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;
        p7 += trans;
        p8 += trans;
        p9 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);
        p7  = rotatedVec(p7.x(), p7.y(), i);
        p8  = rotatedVec(p8.x(), p8.y(), i);
        p9  = rotatedVec(p9.x(), p9.y(), i);


        IntersectionFinderConvexPolygon2D::InputElem ie_p1(p1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_p2(p2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_p3(p3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_p4(p4, 4);
        IntersectionFinderConvexPolygon2D::InputElem ie_p5(p5, 5);
        IntersectionFinderConvexPolygon2D::InputElem ie_p6(p6, 6);
        IntersectionFinderConvexPolygon2D::InputElem ie_p7(p7, 7);
        IntersectionFinderConvexPolygon2D::InputElem ie_p8(p8, 8);

        IntersectionFinderConvexPolygon2D::InputElem ie_p9(p9, 9);

        vector<IntersectionFinderConvexPolygon2D::InputElem*> pts_01;
        pts_01.push_back(&ie_p1);
        pts_01.push_back(&ie_p2);
        pts_01.push_back(&ie_p3);
        pts_01.push_back(&ie_p4);
        pts_01.push_back(&ie_p5);
        pts_01.push_back(&ie_p6);
        pts_01.push_back(&ie_p7);
        pts_01.push_back(&ie_p8);

        CHNode* n1;
        CHNode* n2;
        long index_01, index_02;
        mFinder.constructCycle(pts_01, &n1, &n2,true, index_01, index_02);
        CHEdge* e;
        CHNode* cur;
        for (cur = n1; cur->mP != p3; cur=cur->mNext->mDst) {;}
        e = cur->mNext;

        CHNode* n9 = new CHNode(p9, false, 10, -1, CHNode::NT_INTERSECTION);
        CHNode* n3 = mFinder.splitEdge(
                                            e, n9, CHNode::NT_INTERSECTION);
        CHEdge* en = n3->mNext;
        CHEdge* ep = n3->mPrev;

        EXPECT_EQ(n3->mP,               p9);
        EXPECT_EQ(n3->mPeer,            n9);
        EXPECT_EQ(n3->mType,            CHNode::NT_INTERSECTION);
        EXPECT_EQ(n3->mIsTouchingPoint, false);
        EXPECT_EQ(n3->mTop,             false);
        EXPECT_EQ(n3->mBottom,          false);
        EXPECT_EQ(n3->mThisIsA,         true);
        EXPECT_EQ(n3->mIndex,           3);
        EXPECT_EQ(n3->mIndexAux,        4);

        EXPECT_EQ(n9->mPeer,            n3);


        EXPECT_EQ(en->mDst->mP,    p4);
        EXPECT_EQ(en->mDst->mPrev, en);
        EXPECT_EQ(en->mSrc,        n3);

        EXPECT_EQ(ep->mSrc->mP,    p3);
        EXPECT_EQ(ep->mSrc->mNext, ep);
        EXPECT_EQ(ep->mDst,        n3);

        EXPECT_EQ(n3->mPeer, n9);
        EXPECT_EQ(n9->mPeer, n3);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::classifyCoincidentNodes()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test14) {


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(0.0, -1.0);
        Vec2 p2(0.0,  0.0);
        Vec2 p3(-1.0, 1.0);
        Vec2 p4(0.0, -2.0);
        Vec2 p5(0.0,  0.0);
        Vec2 p6(-2.0, 2.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);


        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                            classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_1_1);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(0.0, -1.0);
        Vec2 p2(0.0,  0.0);
        Vec2 p3(-1.0, 0.0);
        Vec2 p4(0.0, -2.0);
        Vec2 p5(0.0,  0.0);
        Vec2 p6(-2.0, 2.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                              classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_1_2);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);


        Vec2 p1(0.0, -1.0);
        Vec2 p2(0.0,  0.0);
        Vec2 p3(-1.0, 1.0);
        Vec2 p4(0.0, -2.0);
        Vec2 p5(0.0,  0.0);
        Vec2 p6(-2.0, 0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                              classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_1_3);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(0.0, -1.0);
        Vec2 p2(0.0,  0.0);
        Vec2 p3(0.0, 1.0);
        Vec2 p4(0.0, 2.0);
        Vec2 p5(0.0,  0.0);
        Vec2 p6(0.0, -2.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                              classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_2_1);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(0.0, -1.0);
        Vec2 p2(0.0,  0.0);
        Vec2 p3(-2.0, 2.0);
        Vec2 p4(1.0,  1.0);
        Vec2 p5(0.0,  0.0);
        Vec2 p6(0.0, -2.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                              classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_2_2);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(1.0,  1.0);
        Vec2 p2(0.0,  0.0);
        Vec2 p3(0.0, -2.0);
        Vec2 p4(0.0, -1.0);
        Vec2 p5(0.0,  0.0);
        Vec2 p6(-2.0, 2.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                              classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_3_1);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(1.0,  -1.0);
        Vec2 p2(0.0,  0.0);
        Vec2 p3(-2.0, 0.0);
        Vec2 p4(0.5, -1.0);
        Vec2 p5(0.0,  0.0);
        Vec2 p6(-1.0, 0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                               classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_4_1);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(0.5, -1.0);
        Vec2 p2(0.0,  0.0);
        Vec2 p3(-1.0, 0.0);
        Vec2 p4(1.0,  -1.0);
        Vec2 p5(0.0,  0.0);
        Vec2 p6(-2.0, 0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                              classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_4_2);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(1.0, -1.0);
        Vec2 p2(0.0,  0.0);
        Vec2 p3(-1.0,-1.0);
        Vec2 p4(-1.0, 1.0);
        Vec2 p5(0.0,  0.0);
        Vec2 p6(1.0,  1.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                              classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_5);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(-0.5,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.0);
        Vec2 p4( 1.0,  0.5);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-0.5, -1.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                               classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_7);

    }



    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1( 1.0,  0.5);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-0.5, -1.0);
        Vec2 p4(-0.5,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                               classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_8);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1( 1.0,  0.5);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -0.5);
        Vec2 p4( 1.0,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                              classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_9);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1( 1.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.0);
        Vec2 p4( 1.0,  0.5);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -0.5);
        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;
        p6 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);
        p6  = rotatedVec(p6.x(), p6.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n6 = new CHNode(p6, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();
        CHEdge* e56 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;
        n5->mNext = e56;
        e56->mSrc = n5;
        e56->mDst = n6;
        n6->mPrev = e56;

        enum _classifier c1 = mFinder.
                                             classifyCoincidentNodes(n2, n5);
        EXPECT_EQ(c1, TYPE_CV_10);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::classifyVertexOnEdge();
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test15) {


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1( 1.0,  0.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  0.0);
        Vec2 p4(-2.0,  0.0);
        Vec2 p5( 2.0,  0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;

        enum _classifier c1 = mFinder.
                                                 classifyVertexOnEdge(n2, e45);
        EXPECT_EQ(c1, TYPE_OE_1_1);
    }

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1( 1.0,  0.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0, -1.0);
        Vec2 p4(-2.0,  0.0);
        Vec2 p5( 2.0,  0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;

        enum _classifier c1 = mFinder.
                                                 classifyVertexOnEdge(n2, e45);
        EXPECT_EQ(c1, TYPE_OE_1_2);
    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1( 0.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  0.0);
        Vec2 p4(-2.0,  0.0);
        Vec2 p5( 2.0,  0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;

        enum _classifier c1 = mFinder.
                                                 classifyVertexOnEdge(n2, e45);
        EXPECT_EQ(c1, TYPE_OE_2_1);
    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(-1.0,  0.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0,  0.0);
        Vec2 p4(-2.0,  0.0);
        Vec2 p5( 2.0,  0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;

        enum _classifier c1 = mFinder.
                                                 classifyVertexOnEdge(n2, e45);
        EXPECT_EQ(c1, TYPE_OE_3_1);
    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(-1.0,  0.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0,  1.0);
        Vec2 p4(-2.0,  0.0);
        Vec2 p5( 2.0,  0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;

        enum _classifier c1 = mFinder.
                                                 classifyVertexOnEdge(n2, e45);
        EXPECT_EQ(c1, TYPE_OE_3_2);
    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(-1.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0,  0.0);
        Vec2 p4(-2.0,  0.0);
        Vec2 p5( 2.0,  0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;

        enum _classifier c1 = mFinder.
                                                 classifyVertexOnEdge(n2, e45);
        EXPECT_EQ(c1, TYPE_OE_4_1);
    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1( 1.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0, -1.0);
        Vec2 p4(-2.0,  0.0);
        Vec2 p5( 2.0,  0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;

        enum _classifier c1 = mFinder.
                                                 classifyVertexOnEdge(n2, e45);
        EXPECT_EQ(c1, TYPE_OE_5_1);
    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(-1.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0,  1.0);
        Vec2 p4(-2.0,  0.0);
        Vec2 p5( 2.0,  0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;

        enum _classifier c1 = mFinder.
                                                 classifyVertexOnEdge(n2, e45);
        EXPECT_EQ(c1, TYPE_OE_5_2);
    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(-1.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0,  1.0);
        Vec2 p4(-2.0,  0.0);
        Vec2 p5( 2.0,  0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;

        enum _classifier c1 = mFinder.
                                                classifyVertexOnEdge(n2, e45);
        EXPECT_EQ(c1, TYPE_OE_6_1);
    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1( 1.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0, -0.5);
        Vec2 p4(-2.0,  0.0);
        Vec2 p5( 2.0,  0.0);

        p1 += trans;
        p2 += trans;
        p3 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p3  = rotatedVec(p3.x(), p3.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n3 = new CHNode(p3, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e23 = new CHEdge();
        CHEdge* e45 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;
        n2->mNext = e23;
        e23->mSrc = n2;
        e23->mDst = n3;
        n3->mPrev = e23;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;

        enum _classifier c1 = mFinder.
                                                classifyVertexOnEdge(n2, e45);
        EXPECT_EQ(c1, TYPE_OE_6_2);
    }

}


/**  @brief IntersectionFinderConvexPolygon2D::classifyCrossing();
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test16) {


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

       Vec2 p1(-0.5, -1.0);
        Vec2 p2( 0.5,  1.0);
        Vec2 p4(-2.0,  0.0);
        Vec2 p5( 2.0,  0.0);

        p1 += trans;
        p2 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e45 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;

        enum _classifier c1 = mFinder.
                                                    classifyCrossing(e12, e45);
        EXPECT_EQ(c1, TYPE_IS_4);

    }


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 p1(-0.5, -1.0);
        Vec2 p2( 0.5,  1.0);
        Vec2 p4( 2.0,  0.0);
        Vec2 p5(-2.0,  0.0);

        p1 += trans;
        p2 += trans;
        p4 += trans;
        p5 += trans;

        p1  = rotatedVec(p1.x(), p1.y(), i);
        p2  = rotatedVec(p2.x(), p2.y(), i);
        p4  = rotatedVec(p4.x(), p4.y(), i);
        p5  = rotatedVec(p5.x(), p5.y(), i);

        CHNode* n1 = new CHNode(p1, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n2 = new CHNode(p2, true, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n4 = new CHNode(p4, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);
        CHNode* n5 = new CHNode(p5, false, 1, -1, CHNode::NT_ORIGINAL_VERTEX);

        CHEdge* e12 = new CHEdge();
        CHEdge* e45 = new CHEdge();

        n1->mNext = e12;
        e12->mSrc = n1;
        e12->mDst = n2;
        n2->mPrev = e12;

        n4->mNext = e45;
        e45->mSrc = n4;
        e45->mDst = n5;
        n5->mPrev = e45;

        enum _classifier c1 = mFinder.
                                                    classifyCrossing(e12, e45);
        EXPECT_EQ(c1, TYPE_IS_3);

    }
}


/**  @brief CV_1_1
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test17) {

    // AL-BL sweeping down

    for (long i = 0; i < 2; i++) {
        //       *
        //      / \
        //     /   \
        //    /     \
        //  na/nb pa/pb

        Vec2 p1( 1.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0, -1.0);

        Vec2 p4( 1.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0, -1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_1(
                                                  &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
    }


    for (long i = 0; i < 2; i++) {

        //           *
        //          / \
        //         /   \
        //        /     \
        //      na/nb    *pa
        //                \
        //                 *pb

        Vec2 p1( 1.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-2.0, -2.0);
        Vec2 p4( 1.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0, -1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_1(
                                                   &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                   &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
    }


    for (long i = 0; i < 2; i++) {

        //           *
        //          / \
        //         /   \
        //        /     \
        //     na*       *pa/pb
        //      /       
        //   nb*

        Vec2 p1( 1.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0, -1.0);
        Vec2 p4( 1.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-2.0, -2.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_1(
                                                  &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mDst, a1);
        auto* b2 = org_b1->mNext->mDst;
        EXPECT_EQ(b2->mP, p3);
        EXPECT_EQ(b1, b2);
        EXPECT_EQ(b1->mPrev->mSrc, org_b1);
        EXPECT_EQ(b1->mNext->mDst->mP, p6);
    }


    for (long i = 0; i < 2; i++) {

        //           *
        //          / \
        //         /   \
        //        /     \
        //     nb*       *pa/pb
        //      /       
        //   na*

        Vec2 p1( 1.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-2.0, -2.0);
        Vec2 p4( 1.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0, -1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_1(
                                                   &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                   &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);

        EXPECT_EQ(org_b1->mNext->mDst, b1);
        auto* a2 = org_a1->mNext->mDst;
        EXPECT_EQ(a2->mP, p6);
        EXPECT_EQ(a1, a2);
        EXPECT_EQ(a1->mPrev->mSrc, org_a1);
        EXPECT_EQ(a1->mNext->mDst->mP, p3);
    }

    for (long i = 0; i < 2; i++) {

        //           *pa/pb
        //          /
        //         *curA/curB
        //        /
        //     nb*
        //      /
        //   na*

        Vec2 p1( 1.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-2.0, -2.0);
        Vec2 p4( 1.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0, -1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_1(
                                                   &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                   &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mNext->mDst, b1);
        auto* a2 = org_a1->mNext->mDst;
        EXPECT_EQ(a2->mP, p6);
        EXPECT_EQ(a1, a2);
        EXPECT_EQ(a1->mPrev->mSrc, org_a1);
        EXPECT_EQ(a1->mNext->mDst->mP, p3);
    }


    for (long i = 0; i < 2; i++) {

        //  pa/pb na/nb  
        //    *     *
        //     \   /
        //      \ /
        //       *
        //     ca/cb

        Vec2 p1(-1.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0,  1.0);
        Vec2 p4(-1.0,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0,  1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_1(
                                                  &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
    }


    // AR-BR sweeping up

    for (long i = 0; i < 2; i++) {

        //  pa/pb na/nb  
        //    *     *
        //     \   /
        //      \ /
        //       *
        //     ca/cb

        Vec2 p1(-1.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0,  1.0);
        Vec2 p4(-1.0,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0,  1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_1(
                                                &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
    }

    for (long i = 0; i < 2; i++) {

        //  *pa
        //   \
        //    *pb   *na/nb
        //     \   /
        //      \ /
        //       *
        //     ca/cb

        Vec2 p1(-2.0,  2.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0,  1.0);
        Vec2 p4(-1.0,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0,  1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_1(
                                                &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);

    }

    for (long i = 0; i < 2; i++) {

        //            *nb
        //  pa/pb    /
        //    *     *na
        //     \   /
        //      \ /
        //       *
        //     ca/cb

        Vec2 p1(-1.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0,  1.0);
        Vec2 p4(-1.0,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 2.0,  2.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_1(
                                                &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);

        EXPECT_EQ(org_a1->mNext->mDst, a1);
        auto* b2 = org_b1->mNext->mDst;
        EXPECT_EQ(b2->mP, p3);
        EXPECT_EQ(b2->mPrev->mSrc, org_b1);
        EXPECT_EQ(b1, b2);
        EXPECT_EQ(b1->mPrev->mSrc, org_b1);
        EXPECT_EQ(b1->mNext->mDst->mP, p6);

    }


    for (long i = 0; i < 2; i++) {

        //            *na
        //  pa/pb    /
        //    *     *nb
        //     \   /
        //      \ /
        //       *
        //     ca/cb

        Vec2 p1(-1.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 2.0,  2.0);
        Vec2 p4(-1.0,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0,  1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_1(
                                               &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);

        EXPECT_EQ(org_b1->mNext->mDst, b1);
        auto* a2 = org_a1->mNext->mDst;
        EXPECT_EQ(a2->mP, p6);
        EXPECT_EQ(a2->mPrev->mSrc, org_a1);
        EXPECT_EQ(a1, a2);
        EXPECT_EQ(a1->mPrev->mSrc, org_a1);
        EXPECT_EQ(a1->mNext->mDst->mP, p3);

    }


    for (long i = 0; i < 2; i++) {

        //            *na
        //           /
        //          *nb
        //         /
        //        /
        //       *ca/cb
        //      /
        //     /
        //    *pa/pb

        Vec2 p1(-1.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 2.0,  2.0);
        Vec2 p4(-1.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0,  1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_1(
                                               &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);

        EXPECT_EQ(org_b1->mNext->mDst, b1);
        auto* a2 = org_a1->mNext->mDst;
        EXPECT_EQ(a2->mP, p6);
        EXPECT_EQ(a2->mPrev->mSrc, org_a1);
        EXPECT_EQ(a1, a2);
        EXPECT_EQ(a1->mPrev->mSrc, org_a1);
        EXPECT_EQ(a1->mNext->mDst->mP, p3);

    }


    for (long i = 0; i < 2; i++) {

        //            *nb
        //           /
        //          *na
        //         /
        //        /
        //       *ca/cb
        //      /
        //     /
        //    *pa/pb

        Vec2 p1(-1.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0,  1.0);
        Vec2 p4(-1.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 2.0,  2.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_1(
                                                &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);

        EXPECT_EQ(org_a1->mNext->mDst, a1);
        auto* b2 = org_b1->mNext->mDst;
        EXPECT_EQ(b2->mP, p3);
        EXPECT_EQ(b2->mPrev->mSrc, org_b1);
        EXPECT_EQ(b1, b2);
        EXPECT_EQ(b1->mPrev->mSrc, org_b1);
        EXPECT_EQ(b1->mNext->mDst->mP, p6);

    }
}


/**  @brief CV_1_2
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test18) {

    for (long i = 0; i < 2; i++) {

        //        *pa/pb
        //       /
        //      *curA/curB
        //     / \
        //  nb*   *na

        Vec2 p1( 1.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.5);

        Vec2 p4( 1.0,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0, -1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_2(
                                                   &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                   &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);

        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
    }



    for (long i = 0; i < 2; i++) {

        //     na 
        //      * nb
        //      | *
        //      |/
        //      *curA/curB
        //     / 
        //    *
        //  pa/pb

        Vec2 p1(-1.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0,  2.0);

        Vec2 p4(-1.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0,  1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_2(
                                               &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
    }

}


/**  @brief CV_1_3
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test19) {


    // AL-BL sweeping down

    for (long i = 0; i < 2; i++) {

        //        *pa/pb
        //       /
        //      *curA/curB
        //     / \
        //  na*   \
        //         *nb

        Vec2 p1( 1.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0, -1.0);

        Vec2 p4( 1.0,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.5);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_3(
                                                  &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(b1, org_b1);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, true);
    }

    // AR-BR sweeping up

    for (long i = 0; i < 2; i++) {

        //     nb
        //      * na
        //      | *
        //      |/
        //      *curA/curB
        //     / 
        //    *
        //  pa/pb

        Vec2 p1(-1.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0,  1.0);

        Vec2 p4(-1.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  2.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_1_3(
                                                &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(b1, org_b1);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, true);
    }

}

/**  @brief CV_2_1
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test20) {

    // AL-BR sweeping down

    for (long i = 0; i < 2; i++) {

        // panb *
        //      |
        //      * curA/curB
        //      |
        // napb *

        Vec2 p1( 0.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0, -1.0);

        Vec2 p4( 0.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_2_1(
                                                 &a1, &b1, true, false, true);
        }
        else {
            mFinder.handleCoincidence(
                                                 &a1, &b1, true, false, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mPrev->mSrc);
    }

    for (long i = 0; i < 2; i++) {

        //   nb *
        //      |
        //   pa *
        //      |
        //      * curA/curB
        //      |
        // napb *

        Vec2 p1( 0.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0, -1.0);

        Vec2 p4( 0.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  2.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_2_1(
                                                  &a1, &b1, true, false, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, true, false, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mPrev->mSrc);
    }

    for (long i = 0; i < 2; i++) {

        // panb *
        //      |
        //      * curA/curB
        //      |
        //   pb *
        //      |
        //   na *

        Vec2 p1( 0.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0, -2.0);

        Vec2 p4( 0.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_2_1(
                                                 &a1, &b1, true, false, true);
        }
        else {
            mFinder.handleCoincidence(
                                                 &a1, &b1, true, false, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        auto* a2 = org_a1->mNext->mDst;
        EXPECT_EQ(a2->mP, p4);
        auto* a3 = a2->mNext->mDst;
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mPrev->mSrc);
    }

    for (long i = 0; i < 2; i++) {

        // panb *
        //      |
        //      * curA/curB
        //      |
        //   na *  b2  <= a1/b1
        //      |
        //   pb *  b3

        Vec2 p1( 0.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0, -1.0);

        Vec2 p4( 0.0, -2.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_2_1(
                                               &a1, &b1, true, false, true);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, true, false, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        auto* b2 = org_b1->mPrev->mSrc;
        EXPECT_EQ(b2->mP, p3);
        auto* b3 = b2->mPrev->mSrc;
        EXPECT_EQ(b3->mP, p4);
        EXPECT_EQ(b3->mNext->mDst, b2);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mPrev->mSrc);
    }

    // AR-BL sweeping down

    for (long i = 0; i < 2; i++) {

        // napb *
        //      |
        //      * curA/curB
        //      |
        // panb *

        Vec2 p1( 0.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0,  1.0);

        Vec2 p4( 0.0,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_2_1(
                                                 &a1, &b1, false, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                 &a1, &b1, false, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(a1, org_a1->mPrev->mSrc);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
    }

    for (long i = 0; i < 2; i++) {

        //   na *
        //      |
        //   pb *
        //      |
        //      * curA/curB
        //      |
        // panb *

        Vec2 p1( 0.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0,  2.0);

        Vec2 p4( 0.0,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_2_1(
                                                &a1, &b1, false, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, false, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(a1, org_a1->mPrev->mSrc);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
    }

    for (long i = 0; i < 2; i++) {

        //   na *
        //      |
        //   pb *
        //      |
        //      * curA/curB
        //      |
        //   pa * <- a1, b1=b2
        //      |
        //   nb * <- b3

        Vec2 p1( 0.0, -1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0,  2.0);

        Vec2 p4( 0.0,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -2.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_2_1(
                                                &a1, &b1, false, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, false, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(a1, org_a1->mPrev->mSrc);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        auto* b2 = org_b1->mNext->mDst;
        EXPECT_EQ(b1, b2);
        EXPECT_EQ(b1->mP, p1);
        auto* b3 = b2->mNext->mDst;
        EXPECT_EQ(b3->mP, p6);
        EXPECT_EQ(b3->mPrev->mSrc, b2);        
    }


    for (long i = 0; i < 2; i++) {

        //   na *
        //      |
        //   pb *
        //      |
        //      * curA/curB
        //      |
        //   nb * <- b1, a1=a2
        //      |
        //   pa * <- a3

        Vec2 p1( 0.0, -2.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0,  2.0);

        Vec2 p4( 0.0,  1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_2_1(
                                                  &a1, &b1, false, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, false, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);

        EXPECT_EQ(a1, org_a1->mPrev->mSrc);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        auto* a2 = org_a1->mPrev->mSrc;
        EXPECT_EQ(a1, a2);
        EXPECT_EQ(a1->mP, p6);
        auto* a3 = a2->mPrev->mSrc;
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(a3->mNext->mDst, a2);        

    }
}


/**  @brief CV_2_2
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test21) {

    // AL-BR sweeping down

    for (long i = 0; i < 2; i++) {

        //  panb*
        //      |
        //      |
        //      * curA/curB
        //     / \
        //    /   \
        // pb*     *na

        Vec2 p1( 0.0,  1.0);
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.0);

        Vec2 p4(-1.0, -1.0);
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0);

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_2_2(
                                                  &a1, &b1, true, false, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, true, false, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);

        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mPrev->mSrc);
    }


    // AR-BL sweeping down

    for (long i = 0; i < 2; i++) {

        // na*     *pb
        //    \   /
        //     \ /
        //      * curA/curB
        //      |
        //      |
        //  panb*

        Vec2 p1( 0.0, -1.0);// pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.0);// na

        Vec2 p4( 1.0,  1.0);// pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0);// nb


        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_2_2(
                                                 &a1, &b1, false, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                 &a1, &b1, false, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mPrev->mSrc);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
    }



    for (long i = 0; i < 2; i++) {

        // na*     *pb
        //    \   /
        //     \ /
        //      * curA/curB
        //      |
        //      |
        //    pa* <= a1,b1=b2 
        //      |
        //    nb* <= b3

        Vec2 p1( 0.0, -1.0);// pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.0);// na

        Vec2 p4( 1.0,  1.0);// pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -2.0);// nb


        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_2_2(
                                                  &a1, &b1, false, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, false, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mPrev->mSrc);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        auto* b2 = org_b1->mNext->mDst;
        EXPECT_EQ(b1, b2);
        EXPECT_EQ(b2->mP, p1);
        auto* b3 = b2->mNext->mDst;
        EXPECT_EQ(b3->mP, p6);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

    }

    for (long i = 0; i < 2; i++) {

        // na*     *pb
        //    \   /
        //     \ /
        //      * curA/curB
        //      |
        //      |
        //    nb* <= b1,a1=a2 
        //      |
        //    pa* <= a3

        Vec2 p1( 0.0, -2.0);// pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.0);// na

        Vec2 p4( 1.0,  1.0);// pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0);// nb


        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_2_2(
                                                  &a1, &b1, false, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, false, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mPrev->mSrc);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        auto* a2 = org_a1->mPrev->mSrc;
        EXPECT_EQ(a1, a2);
        EXPECT_EQ(a2->mP, p6);
        auto* a3 = a2->mPrev->mSrc;
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(a3->mNext->mDst, a2);

    }
}


/**  @brief CV_3_1
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test22) {

    // AL-BR sweeping down

    for (long i = 0; i < 2; i++) {

        // nb*     *pa
        //    \   /
        //     \ /
        //      * curA/curB
        //      |
        //      |
        //  pbna*

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0, -1.0); // na

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_3_1(
                                                  &a1, &b1, true, false, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, true, false, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mPrev->mSrc);
    }


    for (long i = 0; i < 2; i++) {

        // nb*     *pa
        //    \   /
        //     \ /
        //      * curA/curB
        //      |
        //      |
        //    na*
        //      |
        //    pb*

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0, -1.0); // na

        Vec2 p4( 0.0, -2.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_3_1(
                                              &a1, &b1, true, false, true);
        }
        else {
            mFinder.handleCoincidence(
                                              &a1, &b1, true, false, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mPrev->mSrc);
        auto* b2 = org_b1->mPrev->mSrc;
        EXPECT_EQ(b2->mP, p3);
        auto* b3 = b2->mPrev->mSrc;
        EXPECT_EQ(b3->mP, p4);
        EXPECT_EQ(b3->mNext->mDst, b2);

    }

    for (long i = 0; i < 2; i++) {

        // nb*     *pa
        //    \   /
        //     \ /
        //      * curA/curB
        //      |
        //      |
        //    pb*
        //      |
        //    na*

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0, -2.0); // na

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_3_1(
                                               &a1, &b1, true, false, true);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, true, false, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mPrev->mSrc);
        auto* a2 = org_a1->mNext->mDst;
        EXPECT_EQ(a2->mP, p4);
        auto* a3 = a2->mNext->mDst;
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(a3->mPrev->mSrc, a2);

    }

    // AR-BL  sweeping down

    for (long i = 0; i < 2; i++) {

        //      *na
        //      |
        //      *pb
        //      |
        //      *curA/curB
        //     / \
        //    /   \
        // pa*     \
        //          *nb

        Vec2 p1( 1.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0,  2.0); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.5); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_3_1(
                                                 &a1, &b1, false, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                 &a1, &b1, false, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mPrev->mSrc);
        EXPECT_EQ(b1, org_b1);

    }


    for (long i = 0; i < 2; i++) {

        //      *na
        //      |
        //      *pb
        //      |
        //      *curA/curB
        //     / \
        //    /   \
        //   /     *nb
        //pa*

        Vec2 p1( 1.0, -2.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0,  2.0); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_3_1(
                                                 &a1, &b1, false, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                 &a1, &b1, false, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mNext->mDst);

    }

}


/**  @brief CV_4_1
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test23) {

    // AL-BL sweeping down

    for (long i = 0; i < 2; i++) {

        //   *pa *pb
        //   |  /
        //   | /
        //   |/
        //   *curA/curB
        //    \
        //     \
        //  nbna*

        Vec2 p1( 0.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.0); // na

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_4_1(
                                                 &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                 &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, true);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);

    }


    for (long i = 0; i < 2; i++) {

        //   *pa *pb
        //   |  /
        //   | /
        //   |/
        //   *curA/curB
        //    \
        //     \
        //    nb* a1,b1
        //       \
        //      na* a2

        Vec2 p1( 0.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 2.0, -2.0); // na

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_4_1(
                                                 &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                 &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, true);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        EXPECT_EQ(a1->mP, p6);
        auto* a2 = a1->mNext->mDst;
        EXPECT_EQ(a2->mP, p3);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

    }

    for (long i = 0; i < 2; i++) {

        //   *pa *pb
        //   |  /
        //   | /
        //   |/
        //   *curA/curB
        //    \
        //     \
        //    na* a1,b1
        //       \
        //      nb* b2

        Vec2 p1( 0.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.0); // na

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 2.0, -2.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_4_1(
                                                  &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, true);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        EXPECT_EQ(b1->mP, p3);
        auto* b2 = b1->mNext->mDst;
        EXPECT_EQ(b2->mP, p6);
        EXPECT_EQ(b2->mPrev->mSrc, b1);

    }

    // AR-BR sweeping up

    for (long i = 0; i < 2; i++) {

        //    *nanb
        //     \
        //      \
        //       *curA/curB
        //      /|
        //     / |
        //  pb*  *pa

        Vec2 p1( 0.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.0); // na

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_4_1(
                                                &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, true);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);

    }


    for (long i = 0; i < 2; i++) {

        //   *nb b2
        //    \
        //     *na a1,b1
        //      \
        //       *curA/curB
        //      /|
        //     / |
        //  pb*  *pa

        Vec2 p1( 0.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.0); // na

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-2.0,  2.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_4_1(
                                               &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, true);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        EXPECT_EQ(b1->mP, p3);
        auto* b2 = b1->mNext->mDst;
        EXPECT_EQ(b2->mP, p6);
        EXPECT_EQ(b2->mPrev->mSrc, b1);

    }


    for (long i = 0; i < 2; i++) {

        //   *na a2
        //    \
        //     *nb a1,b1
        //      \
        //       *curA/curB
        //      /|
        //     / |
        //  pb*  *pa

        Vec2 p1( 0.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-2.0,  2.0); // na

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_4_1(
                                                &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, true);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        EXPECT_EQ(a1->mP, p6);
        auto* a2 = a1->mNext->mDst;
        EXPECT_EQ(a2->mP, p3);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

    }

}


/**  @brief CV_4_2
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test24) {


    // AL-BL sweeping down

    for (long i = 0; i < 2; i++) {

        //   *pb *pa
        //   |  /
        //   | /
        //   |/
        //   *curA/curB
        //    \
        //     \
        //  nbna*

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.0); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_4_2(
                                                 &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                 &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);

    }


    for (long i = 0; i < 2; i++) {

        //   *pb *pa
        //   |  /
        //   | /
        //   |/
        //   *curA/curB
        //    \
        //     \
        //    nb* a1,b1
        //       \
        //      na* a2

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 2.0, -2.0); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_4_2(
                                                 &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                 &a1, &b1, true, true, true);
        }
        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        EXPECT_EQ(a1->mP, p6);
        auto* a2 = a1->mNext->mDst;
        EXPECT_EQ(a2->mP, p3);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

    }

    for (long i = 0; i < 2; i++) {

        //   *pb *pa
        //   |  /
        //   | /
        //   |/
        //   *curA/curB
        //    \
        //     \
        //    na* a1,b1
        //       \
        //      nb* b2

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.0); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 2.0, -2.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_4_2(
                                                   &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                   &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        EXPECT_EQ(b1->mP, p3);
        auto* b2 = b1->mNext->mDst;
        EXPECT_EQ(b2->mP, p6);
        EXPECT_EQ(b2->mPrev->mSrc, b1);

    }

    // AR-BR sweeping up

    for (long i = 0; i < 2; i++) {

        //    *nanb
        //     \
        //      \
        //       *curA/curB
        //      /|
        //     / |
        //  pa*  *pb

        Vec2 p1(-1.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.0); // na

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_4_2(
                                               &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);

    }


    for (long i = 0; i < 2; i++) {

        //   *nb b2
        //    \
        //     *na a1,b1
        //      \
        //       *curA/curB
        //      /|
        //     / |
        //  pa*  *pb

        Vec2 p1(-2.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.0); // na

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-2.0,  2.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_4_2(
                                              &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                              &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        EXPECT_EQ(b1->mP, p3);
        auto* b2 = b1->mNext->mDst;
        EXPECT_EQ(b2->mP, p6);
        EXPECT_EQ(b2->mPrev->mSrc, b1);

    }


    for (long i = 0; i < 2; i++) {

        //   *na a2
        //    \
        //     *nb a1,b1
        //      \
        //       *curA/curB
        //      /|
        //     / |
        //  pa*  *pb

        Vec2 p1(-1.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-2.0,  2.0); // na

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_4_2(
                                               &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        EXPECT_EQ(a1->mP, p6);
        auto* a2 = a1->mNext->mDst;
        EXPECT_EQ(a2->mP, p3);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

    }
}


/**  @brief CV_5
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test25) {


    // AL-BR sweeping down

    for (long i = 0; i < 2; i++) {

        //  *nb   *pa
        //   \   /
        //    \ /
        //     *curA/curB
        //    / \
        //   /   \
        //  *pb   *na

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.0); // na

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_5(
                                               &a1, &b1, true, false, true);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, true, false, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_a1->mIsTouchingPoint, true);
        EXPECT_EQ(org_b1->mIsTouchingPoint, true);
        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mPrev->mSrc);

    }


    for (long i = 0; i < 2; i++) {

        //  *nb   *pa
        //   \   /
        //    \ /
        //     *curA/curB
        //    / \
        //   /   \
        //  /     *na
        // *pb

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.0); // na

        Vec2 p4(-2.0, -2.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_5(
                                                &a1, &b1, true, false, true);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, true, false, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_a1->mIsTouchingPoint, true);
        EXPECT_EQ(org_b1->mIsTouchingPoint, true);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1);

    }


    // AR-BL sweeping down

    for (long i = 0; i < 2; i++) {

        //  *na   *pb
        //   \   /
        //    \ /
        //     *curA/curB
        //    / \
        //   /   \
        //  *pa   *nb


        Vec2 p1(-1.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.0); // na

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_5(
                                                &a1, &b1, false, true,  true);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, false, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_a1->mIsTouchingPoint, true);
        EXPECT_EQ(org_b1->mIsTouchingPoint, true);
        EXPECT_EQ(b1, org_b1);
        EXPECT_EQ(a1, org_a1->mPrev->mSrc);

    }


    for (long i = 0; i < 2; i++) {

        //  *na   *pb
        //   \   /
        //    \ /
        //     *curA/curB
        //    / \
        //   /   \
        //  /     *nb
        // *pa

        Vec2 p1(-2.0, -2.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.0); // na

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_5(
                                               &a1, &b1, false, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, false, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_a1->mIsTouchingPoint, true);
        EXPECT_EQ(org_b1->mIsTouchingPoint, true);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
        EXPECT_EQ(a1, org_a1);

    }
}


/**  @brief CV_7
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test26) {


    // AL-BL sweeping down

    for (long i = 0; i < 2; i++) {

        //  *pa *pb
        //  |  /
        //  | /
        //  |/
        //  *curA/curB
        //  |\
        //  | \
        //  |  \
        //  *nb \
        //       *na

        Vec2 p1( 0.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.5); // na

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_7(
                                                  &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, true);
        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mNext->mDst);

    }


    for (long i = 0; i < 2; i++) {

        //  *pa *pb
        //  |  /
        //  | /
        //  |/
        //  *curA/curB
        //  |\
        //  | \
        //  |  \
        //  |   *na
        //  *nb    

        Vec2 p1( 0.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.0); // na

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.5); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_7(
                                                  &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, true);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1);

    }

    // AR-BR sweeping up

    for (long i = 0; i < 2; i++) {

        // *na
        //  \   *nb
        //   \  |
        //    \ |
        //     \|
        //      *curA/curB
        //     /|
        //    / |
        //   /  |
        //  *pb |
        //      *pa

        Vec2 p1( 0.0, -1.5); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.5); // na

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_7(
                                            &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                            &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, true);
        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mNext->mDst);

    }

    for (long i = 0; i < 2; i++) {

        //      *nb
        //  *na |
        //   \  |
        //    \ |
        //     \|
        //      *curA/curB
        //     /|
        //    / |
        //   /  |
        //  *pb |
        //      *pa

        Vec2 p1( 0.0, -1.5); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.0); // na

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.5); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_7(
                                              &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                              &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, false);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, true);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1);

    }
}


/**  @brief CV_8
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test27) {

    // AL-BL sweeping down

    for (long i = 0; i < 2; i++) {

        //  *pb *pa
        //  |  /
        //  | /
        //  |/
        //  *curA/curB
        //  |\
        //  | \
        //  |  \
        //  *na \
        //       *nb


        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0, -1.0); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.5); // nb


        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_8(
                                                  &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1);
    }


    for (long i = 0; i < 2; i++) {

        //  *pb *pa
        //  |  /
        //  | /
        //  |/
        //  *curA/curB
        //  |\
        //  | \
        //  |  \
        //  |   *nb
        //  *na

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0, -1.5); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb


        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_8(
                                                  &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                  &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
    }

    // AR-BR sweeping up

    for (long i = 0; i < 2; i++) {

        // *nb
        //  \   *na
        //   \  |
        //    \ |
        //     \|
        //      *curA/curB
        //     /|
        //    / |
        //   /  |
        //  *pa |
        //      *pb

        Vec2 p1(-1.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0,  1.0); // na

        Vec2 p4( 0.0, -1.5); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.5); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_8(
                                               &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1);

    }

    for (long i = 0; i < 2; i++) {

        //      *na
        //  *nb |
        //   \  |
        //    \ |
        //     \|
        //      *curA/curB
        //     /|
        //    / |
        //   /  |
        //  *pa |
        //      *pb

        Vec2 p1(-1.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 0.0,  1.5); // na

        Vec2 p4( 0.0, -1.5); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_8(
                                                &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mNext->mDst);

    }
}


/**  @brief CV_9
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test28) {

    // AL-BL sweeping down

    for (long i = 0; i < 2; i++) {

        //   pb  pa                                                 
        //   *  *
        //   | /
        //   |/
        //   *
        //   |\
        //   | \
        //   *  \
        //   nb  *na           

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.5); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0); // nb


        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_9(&a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(&a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
    }


    for (long i = 0; i < 2; i++) {

        //   pb  pa                                                 
        //   *  *
        //   | /
        //   |/
        //   *
        //   |\
        //   | \
        //   |  *na
        //   *nb 

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.0); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.5); // nb


        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_9(&a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(&a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1);
    }


    // AR-BR sweeping up

    for (long i = 0; i < 2; i++) {

        //      *nb                                                 
        // na*  |
        //    \ |
        //     \|
        //      *
        //     /|
        //    / |
        //   *  |
        //   pa *pb

        Vec2 p1(-1.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.0); // na

        Vec2 p4( 0.0, -1.5); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.5); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_9(&a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(&a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1);

    }


    for (long i = 0; i < 2; i++) {

        //  *na                                                 
        //   \  *nb
        //    \ |
        //     \|
        //      *
        //     /|
        //    / |
        //   *  |
        //   pa *pb

        Vec2 p1(-1.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.5); // na

        Vec2 p4( 0.0, -1.5); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb

        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_9(&a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(&a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mNext->mDst);

    }
}


/**  @brief CV_10
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test29) {

    // AL-BL sweeping down

    for (long i = 0; i < 2; i++) {

        //   pb  pa                                                 
        //   *  *
        //   | /
        //   |/
        //   *
        //   |\
        //   | \
        //   *  \
        //   nb  *na           

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.5); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0); // nb


        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_9(
                                               &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
    }


    for (long i = 0; i < 2; i++) {

        //   pb  pa                                                 
        //   *  *
        //   | /
        //   |/
        //   *
        //   |\
        //   | \
        //   |  *na 
        //   *nb

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3( 1.0, -1.0); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.5); // nb


        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_9(
                                                 &a1, &b1, true, true, true);
        }
        else {
            mFinder.handleCoincidence(
                                                 &a1, &b1, true, true, true);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1);
    }


    // AR-BR sweeping up

    for (long i = 0; i < 2; i++) {

        //  na 
        //   *   nb
        //    \  *
        //     \ |
        //      \|
        //       *cur
        //      /|
        //     / |
        //    *  |
        //   pa  *pb

        Vec2 p1(-1.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.5); // na

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb


        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_9(
                                               &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                               &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1);
        EXPECT_EQ(b1, org_b1->mNext->mDst);
    }


    for (long i = 0; i < 2; i++) {

        //   na  *nb
        //    *  |
        //     \ |
        //      \|
        //       *cur
        //      /|
        //     / |
        //    *  |
        //   pa  *pb

        Vec2 p1(-1.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0);
        Vec2 p3(-1.0,  1.0); // na

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.5); // nb


        CHNode *org_a1, *org_b1;
        createPartialChainsCV(p1, p2, p3, p4, p5, p6, &org_a1, &org_b1);
        CHNode *a1 = org_a1;
        CHNode *b1 = org_b1;

        if (i ==0) {
            mFinder.handleCV_9(
                                                &a1, &b1, false, false, false);
        }
        else {
            mFinder.handleCoincidence(
                                                &a1, &b1, false, false, false);
        }

        EXPECT_EQ(org_a1->mPeer, org_b1);
        EXPECT_EQ(org_b1->mPeer, org_a1);
        EXPECT_EQ(org_a1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_b1->mType, CHNode::NT_COINCIDENT_VERTEX);
        EXPECT_EQ(org_a1->mNext->mSrcInterior, true);
        EXPECT_EQ(org_a1->mPrev->mDstInterior, true);
        EXPECT_EQ(org_b1->mNext->mSrcInterior, false);
        EXPECT_EQ(org_b1->mPrev->mDstInterior, false);
        EXPECT_EQ(a1, org_a1->mNext->mDst);
        EXPECT_EQ(b1, org_b1);

    }
}


/**  @brief OE_1_1
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test30) {

    // AL-BL sweeping down

    for (long i = 0; i < 4; i++) {

        // ep**nb
        //   ||
        //   ||
        //   |*cur <= org_n 
        //   ||
        //   ||
        // en**pb

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb


        CHNode *org_n;

        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_e->mSrc;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_1_1(
                             &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, true, false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_1_1(
                             &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false, false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);
    }


    for (long i = 0; i < 4; i++) {

        // ep**nb
        //   ||
        //   ||
        //   |*cur
        //   ||
        //   ||
        // en*|
        //    *pb

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4( 0.0, -1.5); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb


        CHNode *org_n;

        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_e->mSrc;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_1_1(
                                &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                &a_param,  &b_param, true, true, false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_1_1(
                                &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                &b_param,  &a_param, false, false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);
    }


    for (long i = 0; i < 4; i++) {

        // ep**nb
        //   ||
        //   ||
        //   |*cur
        //   ||
        //   ||
        //   |*pb
        // en*

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.5); // en

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb


        CHNode *org_n;

        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_e->mSrc;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_1_1(
                              &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true, true, false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_1_1(
                              &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false, false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);
    }

    // AR-BL sweeping down

    for (long i = 0; i < 4; i++) {

        // en**pb
        //   ||
        //   ||
        //   |*cur <= org_n
        //   ||
        //   ||
        // ep**nb

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;

        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_1_1(
                             &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_1_1(
                             &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false, true, false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        // en**pb
        //   ||
        //   ||
        //   |*cur <= org_n
        //   ||
        //   ||
        // ep*|
        //    *nb

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.5); // nb

        CHNode *org_n;
        CHEdge *org_e;

        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_1_1(
                                &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_1_1(
                                &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                &b_param,  &a_param, false, true, false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);


    }


    for (long i = 0; i < 4; i++) {

        // en**pb
        //   ||
        //   ||
        //   |*cur <= org_n
        //   ||
        //   ||
        //   |*nb
        // ep*

        Vec2 p1( 0.0, -1.5); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;

        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_1_1(
                                &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_1_1(
                                &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                &b_param,  &a_param, false, true, false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);


    }
}


/**  @brief OE_1_2
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test31) {

    // AL-BR sweeping down

    for (long i = 0; i < 4; i++) {

        //nb*   *ep
        //   \  |
        //    \ |
        //     *|cur <= org_n 
        //     ||
        //     ||
        //   pb**en

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb


        CHNode *org_n;

        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_e->mSrc;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_1_2(
                              &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true, true, false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_1_2(
                              &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false, false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }



    for (long i = 0; i < 4; i++) {

        //nb*   *ep
        //   \  |
        //    \ |
        //     *|cur <= org_n 
        //     ||
        //     ||
        //   en*|
        //      *pb

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4( 0.0, -1.5); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb


        CHNode *org_n;

        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_e->mSrc;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_1_2(
                                &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                &a_param,  &b_param, true, true, false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_1_2(
                                &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                &b_param,  &a_param, false, false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //nb*   *ep
        //   \  |
        //    \ |
        //     *|cur <= org_n 
        //     ||
        //     ||
        //     |*pb
        //   en*

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.5); // en

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb


        CHNode *org_n;

        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_e->mSrc;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_1_2(
                               &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                               &a_param,  &b_param, true, true, false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_1_2(
                               &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                               &b_param,  &a_param, false, false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }

    // AR-BL sweeping down

    for (long i = 0; i < 4; i++) {

        // en**pb
        //   ||
        //   ||
        //   |*cur <= org_n
        //   | \
        //   |  \
        // ep*   *nb

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;

        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_1_2(
                               &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                               &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_1_2(
                               &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                               &b_param,  &a_param, false, true, false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        // en**pb
        //   ||
        //   ||
        //   |*cur <= org_n
        //   | \
        //   |  \
        // ep*   \
        //        *nb

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.5); // nb

        CHNode *org_n;
        CHEdge *org_e;

        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_1_2(
                              &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_1_2(
                              &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false, true, false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        // en**pb
        //   ||
        //   ||
        //   |*cur <= org_n
        //   | \
        //   |  \
        //   |   *nb
        // ep*

        Vec2 p1( 0.0, -1.5); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;

        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_1_2(
                              &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_1_2(
                              &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false, true, false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }
}


/**  @brief OE_2_1
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test32) {

    // AL-BR sweeping down

    for (long i = 0; i < 4; i++) {

        //   nb**ep
        //     ||
        //     ||
        //     *|cur <= org_n 
        //    / |
        //   /  |
        //pb*   *en

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_2_1(
                             &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, true, false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_2_1(
                             &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false, false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   nb**ep
        //     ||
        //     ||
        //     *|cur <= org_n 
        //    / |
        //   /  |
        //pb*   |
        //      *en

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.5); // en

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_2_1(
                           &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                           &a_param,  &b_param, true, true, false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_2_1(
                           &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                           &b_param,  &a_param, false, false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   nb**ep
        //     ||
        //     ||
        //     *|cur <= org_n 
        //    / |
        //   /  *en
        //pb*

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4(-1.0, -1.5); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_2_1(
                              &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true, true, false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_2_1(
                              &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false, false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }

    // AR-BL sweeping down

    for (long i = 0; i < 4; i++) {

        //   en*  *pb
        //     | /
        //     |/
        //     * cur <= org_n 
        //     ||
        //     ||
        //   ep**nb

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_2_1(
                             &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_2_1(
                             &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false, true, false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   en*  *pb
        //     | /
        //     |/
        //     * cur <= org_n 
        //     ||
        //     ||
        //   ep*|
        //      *nb

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.5); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_2_1(
                              &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_2_1(
                              &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false, true, false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   en*  *pb
        //     | /
        //     |/
        //     * cur <= org_n 
        //     ||
        //     ||
        //     |*nb
        //   ep*

        Vec2 p1( 0.0, -1.5); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_2_1(
                              &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_2_1(
                              &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false, true, false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }
}


/**  @brief OE_3_1
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test33) {

    // AL-BL sweeping down

    for (long i = 0; i < 4; i++) {

        //   pb**ep
        //     ||
        //     ||
        //     *|cur <= org_n 
        //     ||
        //     ||
        //   nb**en

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_3_1(
                          &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                          &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_3_1(
                          &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                          &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   pb**ep
        //     ||
        //     ||
        //     *|cur <= org_n 
        //     ||
        //     ||
        //   nb*|
        //      *en

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.5); // en

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_3_1(
                                &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_3_1(
                                &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   pb**ep
        //     ||
        //     ||
        //     *|cur <= org_n 
        //     ||
        //     ||
        //     |*en
        //   nb*

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.5); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_3_1(
                              &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_3_1(
                              &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    // AR-BR sweeping up

    for (long i = 0; i < 4; i++) {

        //   nb**en
        //     ||
        //     ||
        //     *|cur <= org_n 
        //     ||
        //     ||
        //   pb**ep

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_3_1(
                           &a_param, &b_param, true, false, false, false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                           &a_param,  &b_param, true, false, false, false);
        }
        else if (i == 2) {
            mFinder.handleOE_3_1(
                           &b_param, &a_param, false,  false, false, false);
        }
        else {
            mFinder.handleVertexOnEdge(
                           &b_param,  &a_param, false, false, false,false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a1);
        EXPECT_EQ(b_param, b1);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //      *en
        //   nb*|
        //     ||
        //     ||
        //     *|cur <= org_n 
        //     ||
        //     ||
        //   pb**ep

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.5); // en

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_3_1(
                             &a_param, &b_param, true, false, false, false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, false, false, false);
        }
        else if (i == 2) {
            mFinder.handleOE_3_1(
                             &b_param, &a_param, false,  false, false, false);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false, false, false,false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b1);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   nb*
        //     |*en
        //     ||
        //     ||
        //     *|cur <= org_n 
        //     ||
        //     ||
        //   pb**ep

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.5); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_3_1(
                         &a_param, &b_param, true, false, false, false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                         &a_param,  &b_param, true, false, false, false);
        }
        else if (i == 2) {
            mFinder.handleOE_3_1(
                         &b_param, &a_param, false,  false, false, false);
        }
        else {
            mFinder.handleVertexOnEdge(
                         &b_param,  &a_param, false, false, false,false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a1);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }
}


/**  @brief OE_3_2
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test34) {

    // AL-BL sweeping down

    for (long i = 0; i < 4; i++) {

        //   ep**pb
        //     ||
        //     ||
        //     |*cur <= org_n 
        //     | \
        //     |  \
        //   en*   *nb

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_3_2(
                              &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_3_2(
                              &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   ep**pb
        //     ||
        //     ||
        //     |*cur <= org_n 
        //     | \
        //     |  \
        //     |   *nb
        //   en*
        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.5); // en

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_3_2(
                             &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_3_2(
                             &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    // AR-BR sweeping up

    for (long i = 0; i < 4; i++) {

        //nb*   *en
        //   \  |
        //    \ |
        //     *| cur <= org_n 
        //     ||
        //     ||
        //   pb**ep

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_3_2(
                          &a_param, &b_param, true, false, false, false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                          &a_param,  &b_param, true, false, false, false);
        }
        else if (i == 2) {
            mFinder.handleOE_3_2(
                          &b_param, &a_param, false, false, false, false);
        }
        else {
            mFinder.handleVertexOnEdge(
                          &b_param,  &a_param,false, false, false, false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a1);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //      *en
        //nb*   |
        //   \  |
        //    \ |
        //     *| cur <= org_n 
        //     ||
        //     ||
        //   pb**ep

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.5); // en

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_3_2(
                      &a_param, &b_param, true, false, false, false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                      &a_param,  &b_param, true, false, false, false);
        }
        else if (i == 2) {
            mFinder.handleOE_3_2(
                      &b_param, &a_param, false, false, false, false);
        }
        else {
            mFinder.handleVertexOnEdge(
                      &b_param,  &a_param,false, false, false, false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b1);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }
}


/**  @brief OE_4_1
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test35) {

    // AL-BL sweeping down

    for (long i = 0; i < 4; i++) {

        //   ep*   *pb
        //     |  /
        //     | /
        //     |*cur <= org_n 
        //     ||
        //     ||
        //   en**nb

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_4_1(
                             &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_4_1(
                             &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   ep*   *pb
        //     |  /
        //     | /
        //     |*cur <= org_n 
        //     ||
        //     ||
        //   en*|
        //      *nb

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.5); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_4_1(
                            &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                            &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_4_1(
                            &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                            &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   ep*   *pb
        //     |  /
        //     | /
        //     |*cur <= org_n 
        //     ||
        //     ||
        //     |*nb
        //   en*

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.5); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0, -1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_4_1(
                               &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                               &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_4_1(
                               &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                               &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    // AR-BR sweeping up

    for (long i = 0; i < 4; i++) {

        //   nb**en
        //     ||
        //     ||
        //     *|cur <= org_n 
        //    / |
        //   /  |
        //pb*   *ep

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_4_1(
                            &a_param, &b_param, true, false, false, false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                            &a_param,  &b_param, true, false, false, false);
        }
        else if (i == 2) {
            mFinder.handleOE_4_1(
                            &b_param, &a_param, false, false, false, false);
        }
        else {
            mFinder.handleVertexOnEdge(
                            &b_param,  &a_param, false, false, false,false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a1);
        EXPECT_EQ(b_param, b1);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //      *en
        //   nb*|
        //     ||
        //     ||
        //     *|cur <= org_n 
        //    / |
        //   /  |
        //pb*   *ep

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.5); // en

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_4_1(
                               &a_param, &b_param, true, false, false, false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                               &a_param,  &b_param, true, false, false, false);
        }
        else if (i == 2) {
            mFinder.handleOE_4_1(
                               &b_param, &a_param, false, false, false, false);
        }
        else {
            mFinder.handleVertexOnEdge(
                               &b_param,  &a_param, false, false, false,false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b1);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   nb* 
        //     |*en
        //     ||
        //     ||
        //     *|cur <= org_n 
        //    / |
        //   /  |
        //pb*   *ep

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 0.0,  1.5); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_4_1(
                              &a_param, &b_param, true, false, false, false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true, false, false, false);
        }
        else if (i == 2) {
            mFinder.handleOE_4_1(
                              &b_param, &a_param, false, false, false, false);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false, false, false,false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a1);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }

}


/**  @brief OE_5_1
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test36) {

    // AL-BL sweeping down

    for (long i = 0; i < 4; i++) {

        //   ep*   *pb
        //      \ /
        //       *cur <= org_n 
        //      / \
        //   nb*   *en

        Vec2 p1(-1.0,  1.0); // ep
        Vec2 p3( 1.0, -1.0); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0, -1.0); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_5_1(
                               &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                               &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_5_1(
                               &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                               &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   ep*   *pb
        //      \ /
        //       *cur <= org_n 
        //      / \
        //     /   *en
        //  nb*
        Vec2 p1(-1.0,  1.0); // ep
        Vec2 p3( 1.0, -1.0); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0, -1.5); // nb


        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_5_1(
                             &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_5_1(
                             &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }

    // AR-BR sweeping up

    for (long i = 0; i < 4; i++) {

        //   en*   *nb
        //      \ /
        //       *cur <= org_n 
        //      / \
        //   pb*   *ep


        Vec2 p1( 1.0, -1.0); // ep
        Vec2 p3(-1.0,  1.0); // en

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0,  1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_5_1(
                            &a_param, &b_param, true, false, false, false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                            &a_param,  &b_param, true,false, false, false);
        }
        else if (i == 2) {
            mFinder.handleOE_5_1(
                            &b_param, &a_param, false,false, false, false);
        }
        else {
            mFinder.handleVertexOnEdge(
                            &b_param,  &a_param, false,false, false, false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b1);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //          *nb
        //   en*   /
        //      \ /
        //       *cur <= org_n 
        //      / \
        //   pb*   *ep


        Vec2 p1( 1.0, -1.0); // ep
        Vec2 p3(-1.0,  1.0); // en

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0,  1.5); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_5_1(
                              &a_param, &b_param, true, false, false, false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true,false, false, false);
        }
        else if (i == 2) {
            mFinder.handleOE_5_1(
                              &b_param, &a_param, false,false, false, false);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false,false, false, false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a1);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    // AR-BL sweeping down

    for (long i = 0; i < 4; i++) {

        //   pb*   *en
        //      \ /
        //       *cur <= org_n 
        //      / \
        //   ep*   *nb

        Vec2 p1(-1.0, -1.0); // ep
        Vec2 p3( 1.0,  1.0); // en

        Vec2 p4(-1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_5_1(
                             &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_5_1(
                             &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false,true, false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);
        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);
        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);
        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   pb*   *en
        //      \ /
        //       *cur <= org_n 
        //      / \
        //     /   *nb
        //  ep*

        Vec2 p1(-1.0, -1.5); // ep
        Vec2 p3( 1.0,  1.0); // en

        Vec2 p4(-1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_5_1(
                               &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                               &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_5_1(
                               &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                               &b_param,  &a_param, false,true, false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);
        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);
        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);
        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    // AL-BR sweeping down

    for (long i = 0; i < 4; i++) {

        //   nb*   *ep
        //      \ /
        //       *cur <= org_n 
        //      / \
        //   en*   *pb

        Vec2 p1( 1.0,  1.0); // ep
        Vec2 p3(-1.0, -1.0); // en

        Vec2 p4( 1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_5_1(
                             &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, true, false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_5_1(
                             &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false,false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);
        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);
        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);
        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   nb*   *ep
        //      \ /
        //       *cur <= org_n 
        //      / \
        //     /   *pb
        //  en*

        Vec2 p1( 1.0,  1.0); // ep
        Vec2 p3(-1.0, -1.5); // en

        Vec2 p4( 1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_5_1(
                                 &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                 &a_param,  &b_param, true, true, false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_5_1(
                                 &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                 &b_param,  &a_param, false,false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);
        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);
        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);
        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }
}


/**  @brief OE_5_2
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test37) {

    // AL-BL sweeping down

    for (long i = 0; i < 4; i++) {

        //   ep*  *pb
        //     | /
        //     |*cur <= org_n 
        //     | \
        //   en*  *nb

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_5_2(
                               &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                               &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_5_2(
                               &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                               &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //   ep*  *pb
        //     | /
        //     |*cur <= org_n 
        //     | \
        //     |  *nb
        //   en*
        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.5); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_5_2(
                              &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_5_2(
                              &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }

    // AR-BR sweeping up

    for (long i = 0; i < 4; i++) {

        // nb*  *en
        //    \ |
        //     *|cur <= org_n
        //    / |
        // pb*  *ep

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_5_2(
                             &a_param, &b_param, true, false,false,false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, false,false,false);
        }
        else if (i == 2) {
            mFinder.handleOE_5_2(
                             &b_param, &a_param, false,false,false,false);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false,false,false,false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a1);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //      *en
        // nb*  |
        //    \ |
        //     *|cur <= org_n
        //    / |
        // pb*  *ep

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.5); // en

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_5_2(
                               &a_param, &b_param, true, false,false,false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                               &a_param,  &b_param, true, false,false,false);
        }
        else if (i == 2) {
            mFinder.handleOE_5_2(
                               &b_param, &a_param, false,false,false,false);
        }
        else {
            mFinder.handleVertexOnEdge(
                               &b_param,  &a_param, false,false,false,false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b1);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }
}


/**  @brief OE_6_1
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test38) {

    // AL-BL sweeping down

    for (long i = 0; i < 4; i++) {

        //pb*  *ep
        //   \ |
        //    \|
        //     *cur <= org_n 
        //     |\
        //     | \
        //   en*  *nb

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4(-1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_6_1(
                                &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_6_1(
                                &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //pb*  *ep
        //   \ |
        //    \|
        //     *cur <= org_n 
        //     |\
        //     | \
        //     |  *nb
        //   en*

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.5); // en

        Vec2 p4(-1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_6_1(
                                 &a_param, &b_param, true, true, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                 &a_param,  &b_param, true, true, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_6_1(
                                 &b_param, &a_param, false, true, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                 &b_param,  &a_param, false, true, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    // AR-BR sweeping up

    for (long i = 0; i < 4; i++) {

        //nb*  *en
        //   \ |
        //    \|
        //     *cur <= org_n 
        //     |\
        //     | \
        //   ep*  *pb

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_6_1(
                               &a_param, &b_param, true, false,false,false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                               &a_param,  &b_param, true, false,false,false);
        }
        else if (i == 2) {
            mFinder.handleOE_6_1(
                               &b_param, &a_param, false, false,false,false);
        }
        else {
            mFinder.handleVertexOnEdge(
                               &b_param,  &a_param, false, false,false,false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a1);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //     *en
        //nb*  |
        //   \ |
        //    \|
        //     *cur <= org_n 
        //     |\
        //     | \
        //   ep*  *pb

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.5); // en

        Vec2 p4( 1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_6_1(
                              &a_param, &b_param, true, false,false,false);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                              &a_param,  &b_param, true, false,false,false);
        }
        else if (i == 2) {
            mFinder.handleOE_6_1(
                              &b_param, &a_param, false, false,false,false);
        }
        else {
            mFinder.handleVertexOnEdge(
                              &b_param,  &a_param, false, false,false,false);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b1);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    // AL-BR sweeping down

    for (long i = 0; i < 4; i++) {

        //    ep*  *nb
        //      | /
        //      |/
        //      *cur <= org_n 
        //     /|
        //    / |
        // pb*  *en

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0,  1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_6_1(
                                &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                &a_param,  &b_param, true, true,false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_6_1(
                                &b_param, &a_param, false,false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                &b_param,  &a_param, false,false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //    ep*  *nb
        //      | /
        //      |/
        //      *cur <= org_n 
        //     /|
        //    / |
        //   /  *en
        //pb*

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4(-1.0, -1.5); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0,  1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_6_1(
                                 &a_param, &b_param, true, true, false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                 &a_param,  &b_param, true, true,false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_6_1(
                                 &b_param, &a_param, false,false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                 &b_param,  &a_param, false,false, true, true);
        }

        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    // AR-BL sweeping down

    for (long i = 0; i < 4; i++) {

        //    en*  *pb
        //      | /
        //      |/
        //      *cur <= org_n 
        //     /|
        //    / |
        // nb*  *ep

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_6_1(
                                 &a_param, &b_param, true,false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                 &a_param,  &b_param, true,false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_6_1(
                                 &b_param, &a_param, false,true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                 &b_param,  &a_param, false,true,false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //    en*  *pb
        //      | /
        //      |/
        //      *cur <= org_n 
        //     /|
        //    / |
        //   /  *ep
        //nb*

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0, -1.5); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_6_1(
                                 &a_param, &b_param, true,false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                 &a_param,  &b_param, true,false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_6_1(
                                 &b_param, &a_param, false,true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                 &b_param,  &a_param, false,true,false, true);
        }

        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }
}


/**  @brief OE_6_2
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test39) {

    // AL-BR sweeping down

    for (long i = 0; i < 4; i++) {

        // nb*  *ep
        //    \ |
        //     \|
        //      *cur <= org_n 
        //     /|
        //    / |
        // pb*  *en

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_6_2(
                                 &a_param, &b_param, true, true,false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                                 &a_param,  &b_param, true, true,false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_6_2(
                                 &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                                 &b_param,  &a_param, false,false,true, true);
        }
        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, true);
        EXPECT_EQ(b2->mIsTouchingPoint, true);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        // nb*  *ep
        //    \ |
        //     \|
        //      *cur <= org_n 
        //     /|
        //    / |
        //   /  *en
        //pb*

        Vec2 p1( 0.0,  1.0); // ep
        Vec2 p3( 0.0, -1.0); // en

        Vec2 p4(-1.0, -1.5); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6(-1.0,  1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_ep;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_6_2(
                               &a_param, &b_param, true, true,false, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                               &a_param,  &b_param, true, true,false, true);
        }
        else if (i == 2) {
            mFinder.handleOE_6_2(
                               &b_param, &a_param, false, false, true, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                               &b_param,  &a_param, false,false,true, true);
        }
        auto* a1 = org_ep;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_n->mNext->mDst;
        auto* b2 = org_n;
        auto* b3 = org_n->mPrev->mSrc;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_en, a3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b1->mPrev->mSrc, b2);
        EXPECT_EQ(b3->mNext->mDst, b2);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, true);
        EXPECT_EQ(b2->mIsTouchingPoint, true);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    // AR-BL sweeping down

    for (long i = 0; i < 4; i++) {

        //    en*  *pb
        //      | /
        //      |/
        //      *cur <= org_n 
        //      |\
        //      | \
        //    ep*  *nb

        Vec2 p1( 0.0, -1.0); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_6_2(
                             &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_6_2(
                             &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false,true, false,true);
        }
        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, true);
        EXPECT_EQ(b2->mIsTouchingPoint, true);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);

    }


    for (long i = 0; i < 4; i++) {

        //    en*  *pb
        //      | /
        //      |/
        //      *cur <= org_n 
        //      |\
        //      | \
        //      |  *nb
        //    ep*

        Vec2 p1( 0.0, -1.5); // ep
        Vec2 p3( 0.0,  1.0); // en

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p5( 0.0,  0.0);
        Vec2 p6( 1.0, -1.0); // nb

        CHNode *org_n;
        CHEdge *org_e;
        createPartialChainsOE(p1, p3,p4, p5, p6, &org_e, &org_n);
        CHNode *org_ep = org_e->mSrc;
        CHNode *org_en = org_e->mDst;

        CHNode *a_param = org_en;
        CHNode *b_param = org_n;

        if (i == 0) {
            mFinder.handleOE_6_2(
                             &a_param, &b_param, true, false, true, true);
        }
        else if (i == 1) {
            mFinder.handleVertexOnEdge(
                             &a_param,  &b_param, true, false, true, true);
        }
        else if (i == 2) {
            mFinder.handleOE_6_2(
                             &b_param, &a_param, false, true, false, true);
        }
        else {
            mFinder.handleVertexOnEdge(
                             &b_param,  &a_param, false,true, false,true);
        }
        auto* a1 = org_en;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_n->mPrev->mSrc;
        auto* b2 = org_n;
        auto* b3 = org_n->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_ep, a3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b1->mNext->mDst, b2);
        EXPECT_EQ(b3->mPrev->mSrc, b2);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p5);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p5);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, true);
        EXPECT_EQ(b2->mIsTouchingPoint, true);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_SPLIT_POINT);
        EXPECT_EQ(b2->mType, CHNode::NT_ORIGINAL_VERTEX);
    }
}


/**  @brief IS_3
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test40) {

    // AL-BL sweeping down
    {
        // pb*     *pa
        //    \   /
        //     \ /
        //      *
        //     / \
        //    /   \
        // na*     *nb

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0); // int
        Vec2 p3(-1.0, -1.0); // na

        Vec2 p4(-1.0,  1.0); // pb
        Vec2 p6( 1.0, -1.0); // nb

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_eap;
        CHNode *b_param = org_ebp;

        mFinder.handleCrossing(
                                         &a_param, &b_param, true, true, true);

        auto* a1 = org_eap;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_ebp;
        auto* b2 = b1->mNext->mDst;
        auto* b3 = b2->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ean, a3);
        EXPECT_EQ(org_ebn, b3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b3->mPrev->mSrc, b2);
        EXPECT_EQ(b2->mPrev->mSrc, b1);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);

    }

    {
        // pb*     *pa
        //    \   /
        //     \ /
        //      *
        //     / \
        //    /   \
        //   /     *nb
        //na*

        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0); // int
        Vec2 p3(-2.0, -2.0); // na

        Vec2 p4(-1.0,  1.0); // pb
        Vec2 p6( 1.0, -1.0); // nb

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_eap;
        CHNode *b_param = org_ebp;

        mFinder.handleCrossing(
                                      &a_param, &b_param, true, true, true);

        auto* a1 = org_eap;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_ebp;
        auto* b2 = b1->mNext->mDst;
        auto* b3 = b2->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_ean, a3);
        EXPECT_EQ(org_ebn, b3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b3->mPrev->mSrc, b2);
        EXPECT_EQ(b2->mPrev->mSrc, b1);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);

    }

    // AL-BR sweeping down

    {
        //    pa*  *nb
        //      | /
        //      |/
        //      *
        //     /|
        //    / |
        // pb*  *na

        Vec2 p1( 0.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0); // int
        Vec2 p3( 0.0, -1.0); // na

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p6( 1.0,  1.0); // nb

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_eap;
        CHNode *b_param = org_ebn;

        mFinder.handleCrossing(
                                       &a_param, &b_param, true, false, true);

        auto* a1 = org_eap;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_ebn;
        auto* b2 = b1->mPrev->mSrc;
        auto* b3 = b2->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_ean, a3);
        EXPECT_EQ(org_ebp, b3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b3->mNext->mDst, b2);
        EXPECT_EQ(b2->mNext->mDst, b1);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);
    }


    {
        //    pa*  *nb
        //      | /
        //      |/
        //      *
        //     /|
        //    / |
        //   /  *na
        //pb*

        Vec2 p1( 0.0,  1.0); // pa
        Vec2 p2( 0.0,  0.0); // int
        Vec2 p3( 0.0, -1.0); // na

        Vec2 p4(-2.0, -2.0); // pb
        Vec2 p6( 1.0,  1.0); // nb

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_eap;
        CHNode *b_param = org_ebn;

        mFinder.handleCrossing(
                                      &a_param, &b_param, true, false, true);

        auto* a1 = org_eap;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_ebn;
        auto* b2 = b1->mPrev->mSrc;
        auto* b3 = b2->mPrev->mSrc;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ean, a3);
        EXPECT_EQ(org_ebp, b3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b3->mNext->mDst, b2);
        EXPECT_EQ(b2->mNext->mDst, b1);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);
    }



    // AR-BL sweeping down

    {
        //    na*  *pb
        //      | /
        //      |/
        //      *
        //     /|
        //    / |
        // nb*  *pa

        Vec2 p1( 0.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0); // int
        Vec2 p3( 0.0,  1.0); // na

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p6(-1.0, -1.0); // nb

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_ean;
        CHNode *b_param = org_ebp;

        mFinder.handleCrossing(
                                      &a_param, &b_param, false, true, true);

        auto* a1 = org_ean;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_ebp;
        auto* b2 = b1->mNext->mDst;
        auto* b3 = b2->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_eap, a3);
        EXPECT_EQ(org_ebn, b3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b3->mPrev->mSrc, b2);
        EXPECT_EQ(b2->mPrev->mSrc, b1);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);
    }


    {
        //    na*  *pb
        //      | /
        //      |/
        //      *
        //     /|
        //    / |
        //   /  *pa
        //nb*

        Vec2 p1( 0.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0); // int
        Vec2 p3( 0.0,  1.0); // na

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p6(-2.0, -2.0); // nb

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_ean;
        CHNode *b_param = org_ebp;

        mFinder.handleCrossing(
                                      &a_param, &b_param, false, true, true);

        auto* a1 = org_ean;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_ebp;
        auto* b2 = b1->mNext->mDst;
        auto* b3 = b2->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_eap, a3);
        EXPECT_EQ(org_ebn, b3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b3->mPrev->mSrc, b2);
        EXPECT_EQ(b2->mPrev->mSrc, b1);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);
    }


    // AR-BR sweeping up

    {
        // nb*     *na
        //    \   /
        //     \ /
        //      *
        //     / \
        //    /   \
        // pa*     *pb

        Vec2 p1(-1.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0); // int
        Vec2 p3( 1.0,  1.0); // na

        Vec2 p4( 1.0, -1.0); // pb
        Vec2 p6(-1.0,  1.0); // nb

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_eap;
        CHNode *b_param = org_ebp;

        mFinder.handleCrossing(
                                     &a_param, &b_param, false, false, false);

        auto* a1 = org_ean;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_ebn;
        auto* b2 = b1->mPrev->mSrc;
        auto* b3 = b2->mPrev->mSrc;

        EXPECT_EQ(a_param, a1);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_eap, a3);
        EXPECT_EQ(org_ebp, b3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b3->mNext->mDst, b2);
        EXPECT_EQ(b2->mNext->mDst, b1);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);
    }


    {
        //          *na
        // nb*     /
        //    \   /
        //     \ /
        //      *
        //     / \
        //    /   \
        // pa*     *pb

        Vec2 p1(-1.0, -1.0); // pa
        Vec2 p2( 0.0,  0.0); // int
        Vec2 p3( 2.0,  2.0); // na

        Vec2 p4( 1.0, -1.0); // pb
        Vec2 p6(-1.0,  1.0); // nb

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_eap;
        CHNode *b_param = org_ebp;

        mFinder.handleCrossing(
                                    &a_param, &b_param, false, false, false);

        auto* a1 = org_ean;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_ebn;
        auto* b2 = b1->mPrev->mSrc;
        auto* b3 = b2->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b1);

        EXPECT_EQ(org_eap, a3);
        EXPECT_EQ(org_ebp, b3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b3->mNext->mDst, b2);
        EXPECT_EQ(b2->mNext->mDst, b1);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, true);
        EXPECT_EQ(b2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mNext->mSrcInterior, false);
        EXPECT_EQ(a2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);
    }
}


/**  @brief IS_4
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test41) {

    // AL-BL sweeping down
    {
        // pa*     *pb
        //    \   /
        //     \ /
        //      *
        //     / \
        //    /   \
        // nb*     *na

        Vec2 p1(-1.0,  1.0); // pa
        Vec2 p3( 1.0, -1.0); // na
        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p6(-1.0, -1.0); // nb

        Vec2 p2( 0.0,  0.0); // int

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_eap;
        CHNode *b_param = org_ebp;

        mFinder.handleCrossing(
                                     &a_param, &b_param, true, true, true);

        auto* a1 = org_eap;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_ebp;
        auto* b2 = b1->mNext->mDst;
        auto* b3 = b2->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_ean, a3);
        EXPECT_EQ(org_ebn, b3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b3->mPrev->mSrc, b2);
        EXPECT_EQ(b2->mPrev->mSrc, b1);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);

    }

    {
        // pa*     *pb
        //    \   /
        //     \ /
        //      *
        //     / \
        //    /   \
        //   /     *na
        //nb*

        Vec2 p1(-1.0,  1.0); // pa
        Vec2 p3( 1.0, -1.0); // na

        Vec2 p4( 1.0,  1.0); // pb
        Vec2 p6(-2.0, -2.0); // nb

        Vec2 p2( 0.0,  0.0); // int

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_eap;
        CHNode *b_param = org_ebp;

        mFinder.handleCrossing(
                                      &a_param, &b_param, true, true, true);

        auto* a1 = org_eap;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_ebp;
        auto* b2 = b1->mNext->mDst;
        auto* b3 = b2->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ean, a3);
        EXPECT_EQ(org_ebn, b3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b3->mPrev->mSrc, b2);
        EXPECT_EQ(b2->mPrev->mSrc, b1);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);

    }

    // AR-BL sweeping down

    {
        //    pb*  *na
        //      | /
        //      |/
        //      *
        //     /|
        //    / |
        // pa*  *nb

        Vec2 p1(-1.0, -1.0); // pa
        Vec2 p3( 1.0,  1.0); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p6( 0.0, -1.0); // nb

        Vec2 p2( 0.0,  0.0); // int

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_ean;
        CHNode *b_param = org_ebp;

        mFinder.handleCrossing(
                                    &a_param, &b_param, false, true, true);


        auto* a1 = org_ean;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_ebp;
        auto* b2 = b1->mNext->mDst;
        auto* b3 = b2->mNext->mDst;

        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_eap, a3);
        EXPECT_EQ(org_ebn, b3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b3->mPrev->mSrc, b2);
        EXPECT_EQ(b2->mPrev->mSrc, b1);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);

    }

    {
        //    pb*  *na
        //      | /
        //      |/
        //      *
        //     /|
        //    / |
        //   /  *nb
        //pa*

        Vec2 p1(-2.0, -2.0); // pa
        Vec2 p3( 1.0,  1.0); // na

        Vec2 p4( 0.0,  1.0); // pb
        Vec2 p6( 0.0, -1.0); // nb

        Vec2 p2( 0.0,  0.0); // int

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_ean;
        CHNode *b_param = org_ebp;

        mFinder.handleCrossing(
                                  &a_param, &b_param, false, true, true);

        auto* a1 = org_ean;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_ebp;
        auto* b2 = b1->mNext->mDst;
        auto* b3 = b2->mNext->mDst;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_eap, a3);
        EXPECT_EQ(org_ebn, b3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b3->mPrev->mSrc, b2);
        EXPECT_EQ(b2->mPrev->mSrc, b1);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p4);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p6);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);

    }


    // AL-BR sweeping down

    {
        //    nb*  *pa
        //      | /
        //      |/
        //      *
        //     /|
        //    / |
        // na*  *pb
        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p3(-1.0, -1.0); // na

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p6( 0.0,  1.0); // nb

        Vec2 p2( 0.0,  0.0); // int

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_eap;
        CHNode *b_param = org_ebn;

        mFinder.handleCrossing(
                                  &a_param, &b_param, true, false, true);

        auto* a1 = org_eap;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_ebn;
        auto* b2 = b1->mPrev->mSrc;
        auto* b3 = b2->mPrev->mSrc;


        EXPECT_EQ(a_param, a3);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_ean, a3);
        EXPECT_EQ(org_ebp, b3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b3->mNext->mDst, b2);
        EXPECT_EQ(b2->mNext->mDst, b1);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);

    }


    {
        //    nb*  *pa
        //      | /
        //      |/
        //      *
        //     /|
        //    / |
        //   /  *pb
        //na*


        Vec2 p1( 1.0,  1.0); // pa
        Vec2 p3(-2.0, -2.0); // na

        Vec2 p4( 0.0, -1.0); // pb
        Vec2 p6( 0.0,  1.0); // nb

        Vec2 p2( 0.0,  0.0); // int

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_eap;
        CHNode *b_param = org_ebn;

        mFinder.handleCrossing(
                                 &a_param, &b_param, true, false, true);

        auto* a1 = org_eap;
        auto* a2 = a1->mNext->mDst;
        auto* a3 = a2->mNext->mDst;

        auto* b1 = org_ebn;
        auto* b2 = b1->mPrev->mSrc;
        auto* b3 = b2->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b3);

        EXPECT_EQ(org_ean, a3);
        EXPECT_EQ(org_ebp, b3);

        EXPECT_EQ(a3->mPrev->mSrc, a2);
        EXPECT_EQ(a2->mPrev->mSrc, a1);

        EXPECT_EQ(b3->mNext->mDst, b2);
        EXPECT_EQ(b2->mNext->mDst, b1);

        EXPECT_EQ(a1->mP, p1);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p3);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);
    }

    // AR-BR sweeping up

    {
        // na*     *nb
        //    \   /
        //     \ /
        //      *
        //     / \
        //    /   \
        // pb*     *pa

        Vec2 p1( 1.0, -1.0); // pa
        Vec2 p3(-1.0,  1.0); // na

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p6( 1.0,  1.0); // nb

        Vec2 p2( 0.0,  0.0); // int

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_eap;
        CHNode *b_param = org_ebp;

        mFinder.handleCrossing(
                                &a_param, &b_param, false, false, false);

        auto* a1 = org_ean;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_ebn;
        auto* b2 = b1->mPrev->mSrc;
        auto* b3 = b2->mPrev->mSrc;

        EXPECT_EQ(a_param, a2);
        EXPECT_EQ(b_param, b1);

        EXPECT_EQ(org_eap, a3);
        EXPECT_EQ(org_ebp, b3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b3->mNext->mDst, b2);
        EXPECT_EQ(b2->mNext->mDst, b1);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);
    }

    {
        //          *nb
        // na*     /
        //    \   /
        //     \ /
        //      *
        //     / \
        //    /   \
        // pb*     *pa

        Vec2 p1( 1.0, -1.0); // pa
        Vec2 p3(-1.0,  1.0); // na

        Vec2 p4(-1.0, -1.0); // pb
        Vec2 p6( 2.0,  2.0); // nb

        Vec2 p2( 0.0,  0.0); // int

        CHEdge *org_ea;
        CHEdge *org_eb;

        createPartialChainsIS(p1, p3, p4, p6, &org_ea, &org_eb);
        CHNode *org_eap = org_ea->mSrc;
        CHNode *org_ean = org_ea->mDst;
        CHNode *org_ebp = org_eb->mSrc;
        CHNode *org_ebn = org_eb->mDst;

        CHNode *a_param = org_eap;
        CHNode *b_param = org_ebp;

        mFinder.handleCrossing(
                                &a_param, &b_param, false, false, false);

        auto* a1 = org_ean;
        auto* a2 = a1->mPrev->mSrc;
        auto* a3 = a2->mPrev->mSrc;

        auto* b1 = org_ebn;
        auto* b2 = b1->mPrev->mSrc;
        auto* b3 = b2->mPrev->mSrc;

        EXPECT_EQ(a_param, a1);
        EXPECT_EQ(b_param, b2);

        EXPECT_EQ(org_eap, a3);
        EXPECT_EQ(org_ebp, b3);

        EXPECT_EQ(a3->mNext->mDst, a2);
        EXPECT_EQ(a2->mNext->mDst, a1);

        EXPECT_EQ(b3->mNext->mDst, b2);
        EXPECT_EQ(b2->mNext->mDst, b1);

        EXPECT_EQ(a1->mP, p3);
        EXPECT_EQ(a2->mP, p2);
        EXPECT_EQ(a3->mP, p1);
        EXPECT_EQ(b1->mP, p6);
        EXPECT_EQ(b2->mP, p2);
        EXPECT_EQ(b3->mP, p4);

        EXPECT_EQ(a2->mPeer, b2);
        EXPECT_EQ(b2->mPeer, a2);
        EXPECT_EQ(a2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mIsTouchingPoint, false);
        EXPECT_EQ(b2->mNext->mSrcInterior, false);
        EXPECT_EQ(b2->mPrev->mDstInterior, true);
        EXPECT_EQ(a2->mNext->mSrcInterior, true);
        EXPECT_EQ(a2->mPrev->mDstInterior, false);
        EXPECT_EQ(a2->mType, CHNode::NT_INTERSECTION);
        EXPECT_EQ(b2->mType, CHNode::NT_INTERSECTION);
    }

}

static void createHalfHull(
           double* pts, long len, CHNode** top, CHNode** bottom, bool thisIsA)
{
    CHNode* prevN;
    CHEdge* e;
    CHNode* firstN;

    for (long i = 0; i < len; i++) {
        Vec2 v(pts[2*i],pts[2*i+1]);
//        cerr << "create: " << v.x() << "," << v.y() << "\n";
        CHNode* n = new CHNode(v, thisIsA, i, -1, CHNode::NT_ORIGINAL_VERTEX);
        if (i == 0) {
            firstN    = n;
             (*top)    = n;
            (*bottom) = n;
        }
        else {
            e = new CHEdge();
            prevN->mNext = e;
            e->mSrc      = prevN;
            e->mDst      = n;
            n->mPrev     = e;
            if ((*top)->mP.y() < v.y() ||
                ((*top)->mP.y() == v.y() && (*top)->mP.x() > v.x()) ) {
                (*top) = n;
            }
            if ((*bottom)->mP.y() > v.y() ||
                ((*bottom)->mP.y() == v.y() && (*bottom)->mP.x() < v.x()) ) {
                (*bottom) = n;
            }
        }
        prevN = n;
    }
    e = new CHEdge();
    prevN->mNext  = e;
    e->mSrc       = prevN;
    e->mDst       = firstN;
    firstN->mPrev = e;

    (*top)->mTop = true;
    (*bottom)->mBottom = true;
//  cerr << "top: " << (*top)->mP.x() << "," << (*top)->mP.y() << "\n";
//  cerr << "bottom: " << (*bottom)->mP.x() << "," << (*bottom)->mP.y() <<"\n";
}

enum testVertexType {
    TT_NONE,
    TT_ORIGINAL_VERTEX,
    TT_ORIGINAL_VERTEX_ON_EDGE,
    TT_ORIGINAL_VERTEX_ON_EDGE_NEXT_INTERIOR,
    TT_ORIGINAL_VERTEX_ON_EDGE_PREV_INTERIOR,
    TT_ORIGINAL_VERTEX_ON_EDGE_BOTH_INTERIOR,
    TT_INTSEC_NEXT_INTERIOR,
    TT_INTSEC_PREV_INTERIOR,
    TT_COINCIDENT_VERTEX,
    TT_COINCIDENT_VERTEX_NEXT_INTERIOR,
    TT_COINCIDENT_VERTEX_PREV_INTERIOR,
    TT_COINCIDENT_VERTEX_BOTH_INTERIOR,
    TT_SPLIT_POINT,
    TT_SPLIT_POINT_NEXT_INTERIOR,
    TT_SPLIT_POINT_PREV_INTERIOR,
    TT_SPLIT_POINT_BOTH_INTERIOR,
    TT_TOUCHING_POINT_ORIGINAL_VERTEX,
    TT_TOUCHING_POINT_COINCIDENT_VERTEX,
    TT_TOUCHING_POINT_SPLIT_POINT,
    TT_END
};


class eElem {
  public:
    eElem(double x, double y, enum testVertexType type):
        mX(x),mY(y),mType(type){;}
    ~eElem(){;}
    double              mX;
    double              mY;
    enum testVertexType mType;
};


static bool inspectNode(CHNode* n, eElem& e)
{
    if (e.mX != 0.0 || e.mY != 0.0){
        if ( fabs(n->mP.x() - e.mX) > EPSILON_LINEAR ||
             fabs(n->mP.y() - e.mY) > EPSILON_LINEAR    ) {
            return false;
        }
    }

    switch(e.mType) {

      case TT_ORIGINAL_VERTEX:

        // Original Vertex.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer==nullptr                    &&
             !(n->mNext->mSrcInterior)            &&
             !(n->mPrev->mDstInterior)            &&
             n->mType==CHNode::NT_ORIGINAL_VERTEX    ) {

            return true;
        }
        else {
            return false;
        }
        break;

      case TT_ORIGINAL_VERTEX_ON_EDGE:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             !(n->mNext->mSrcInterior)            &&
             !(n->mPrev->mDstInterior)            &&
             n->mType==CHNode::NT_ORIGINAL_VERTEX    ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_ORIGINAL_VERTEX_ON_EDGE_NEXT_INTERIOR:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             n->mNext->mSrcInterior               &&
             !(n->mPrev->mDstInterior)            &&
             n->mType==CHNode::NT_ORIGINAL_VERTEX    ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_ORIGINAL_VERTEX_ON_EDGE_PREV_INTERIOR:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             !(n->mNext->mSrcInterior)            &&
             n->mPrev->mDstInterior               &&
             n->mType==CHNode::NT_ORIGINAL_VERTEX    ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_ORIGINAL_VERTEX_ON_EDGE_BOTH_INTERIOR:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             n->mNext->mSrcInterior               &&
             n->mPrev->mDstInterior               &&
             n->mType==CHNode::NT_ORIGINAL_VERTEX    ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_SPLIT_POINT:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             !(n->mNext->mSrcInterior)            &&
             !(n->mPrev->mDstInterior)            &&
             n->mType==CHNode::NT_SPLIT_POINT        ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_SPLIT_POINT_NEXT_INTERIOR:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             n->mNext->mSrcInterior               &&
             !(n->mPrev->mDstInterior)            &&
             n->mType==CHNode::NT_SPLIT_POINT        ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_SPLIT_POINT_PREV_INTERIOR:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             !(n->mNext->mSrcInterior)            &&
             n->mPrev->mDstInterior               &&
             n->mType==CHNode::NT_SPLIT_POINT        ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_SPLIT_POINT_BOTH_INTERIOR:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             n->mNext->mSrcInterior               &&
             n->mPrev->mDstInterior               &&
             n->mType==CHNode::NT_SPLIT_POINT        ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_INTSEC_NEXT_INTERIOR:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             n->mNext->mSrcInterior               &&
             !(n->mPrev->mDstInterior)            &&
             n->mType==CHNode::NT_INTERSECTION        ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_INTSEC_PREV_INTERIOR:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             !(n->mNext->mSrcInterior)            &&
             n->mPrev->mDstInterior               &&
             n->mType==CHNode::NT_INTERSECTION        ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_COINCIDENT_VERTEX:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             !(n->mNext->mSrcInterior)            &&
             !(n->mPrev->mDstInterior)            &&
             n->mType==CHNode::NT_COINCIDENT_VERTEX  ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_COINCIDENT_VERTEX_NEXT_INTERIOR:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             n->mNext->mSrcInterior               &&
             !(n->mPrev->mDstInterior)            &&
             n->mType==CHNode::NT_COINCIDENT_VERTEX  ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_COINCIDENT_VERTEX_PREV_INTERIOR:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             !(n->mNext->mSrcInterior)            &&
             n->mPrev->mDstInterior               &&
             n->mType==CHNode::NT_COINCIDENT_VERTEX  ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_COINCIDENT_VERTEX_BOTH_INTERIOR:

        // Original Vertex on Edge.
        if ( !(n->mIsTouchingPoint)               &&
             n->mPeer!=nullptr                    &&
             n->mNext->mSrcInterior               &&
             n->mPrev->mDstInterior               &&
             n->mType==CHNode::NT_COINCIDENT_VERTEX  ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_TOUCHING_POINT_ORIGINAL_VERTEX:

        // Original Vertex on Edge.
        if ( n->mIsTouchingPoint               &&
             n->mPeer!=nullptr                 &&
             !(n->mNext->mSrcInterior)         &&
             !(n->mPrev->mDstInterior)         &&
             n->mType==CHNode::NT_ORIGINAL_VERTEX  ) {
            return true;
        }
        else {
            return false;
        }
        break;

      case TT_TOUCHING_POINT_COINCIDENT_VERTEX:

        // Original Vertex on Edge.
        if ( n->mIsTouchingPoint               &&
             n->mPeer!=nullptr                 &&
             !(n->mNext->mSrcInterior)         &&
             !(n->mPrev->mDstInterior)         &&
             n->mType==CHNode::NT_COINCIDENT_VERTEX  ) {
            return true;
        }
        else {
            return false;
        }
        break;


      case TT_TOUCHING_POINT_SPLIT_POINT:

        // Original Vertex on Edge.
        if ( n->mIsTouchingPoint               &&
             n->mPeer!=nullptr                 &&
             !(n->mNext->mSrcInterior)         &&
             !(n->mPrev->mDstInterior)         &&
             n->mType==CHNode::NT_SPLIT_POINT ) {
            return true;
        }
        else {
            return false;
        }
        break;

      default:
        return false;
    }
}


static bool inspectHullsLL(
    CHNode*        topA,
    CHNode*        topB,
    vector<eElem>& expA,
    vector<eElem>& expB,
    bool           print
) {
    if (print) 
        cerr << "Hull A:\n";
    CHNode* cA;
    long i = 0;
    bool res;
    for (cA = topA; cA->mBottom != true; cA = cA->mNext->mDst) {
        if (print) 
            cA->debugPrint(cerr);
        res = inspectNode(cA, expA[i]);
        if (!res && !print) {return false;}
        i++;
    }
    if (print) 
        cA->debugPrint(cerr);    
    res = inspectNode(cA, expA[i]);
    if (!res && !print) {return false;}

    if (print) 
        cerr << "Hull B:\n";
    CHNode* cB;
    i = 0;
    for (cB = topB; cB->mBottom != true; cB = cB->mNext->mDst) {
        if (print) 
            cB->debugPrint(cerr);
        res = inspectNode(cB, expB[i]);
        if (!res && !print) {return false;}
        i++;
    }
    if (print) 
        cB->debugPrint(cerr);
    res = inspectNode(cB, expB[i]);
    if (!res && !print) {return false;}

    return true;
}


static bool inspectHullsRR(
    CHNode*        bottomA,
    CHNode*        bottomB,
    vector<eElem>& expA,
    vector<eElem>& expB,
    bool           print
) {
    if (print) 
        cerr << "Hull A:\n";
    CHNode* cA;
    long i = 0;
    bool res;
    for (cA = bottomA; cA->mTop != true; cA = cA->mNext->mDst) {
        if (print) 
            cA->debugPrint(cerr);
        res = inspectNode(cA, expA[i]);
        if (!res && !print) {return false;}
        i++;
    }
    if (print) 
        cA->debugPrint(cerr);    
    res = inspectNode(cA, expA[i]);
    if (!res && !print) {return false;}

    if (print) 
        cerr << "Hull B:\n";
    CHNode* cB;
    i = 0;
    for (cB = bottomB; cB->mTop != true; cB = cB->mNext->mDst) {
        if (print) 
            cB->debugPrint(cerr);
        res = inspectNode(cB, expB[i]);
        if (!res && !print) {return false;}
        i++;
    }
    if (print) 
        cB->debugPrint(cerr);
    res = inspectNode(cB, expB[i]);
    if (!res && !print) {return false;}

    return true;
}


static bool inspectHullsLR(
    CHNode*        topA,
    CHNode*        topB,
    vector<eElem>& expA,
    vector<eElem>& expB,
    bool           print
) {
    if (print) 
        cerr << "Hull A:\n";
    CHNode* cA;
    long i = 0;
    bool res;
    for (cA = topA; cA->mBottom != true; cA = cA->mNext->mDst) {
        if (print) 
            cA->debugPrint(cerr);
        res = inspectNode(cA, expA[i]);
        if (!res && !print) {return false;}
        i++;
    }
    if (print) 
        cA->debugPrint(cerr);    
    res = inspectNode(cA, expA[i]);
    if (!res && !print) {return false;}

    if (print) 
        cerr << "Hull B:\n";
    CHNode* cB;
    i = 0;
    for (cB = topB; cB->mBottom != true; cB = cB->mPrev->mSrc) {
        if (print) 
            cB->debugPrint(cerr);
        res = inspectNode(cB, expB[i]);
        if (!res && !print) {return false;}
        i++;
    }
    if (print) 
        cB->debugPrint(cerr);
    res = inspectNode(cB, expB[i]);
    if (!res && !print) {return false;}

    return true;
}


static bool inspectHullsRL(
    CHNode*        topA,
    CHNode*        topB,
    vector<eElem>& expA,
    vector<eElem>& expB,
    bool           print
) {
    if (print) 
        cerr << "Hull A:\n";
    CHNode* cA;
    long i = 0;
    bool res;
    for (cA = topA; cA->mBottom != true; cA = cA->mPrev->mSrc) {
        if (print) 
            cA->debugPrint(cerr);
        res = inspectNode(cA, expA[i]);
        if (!res && !print) {return false;}
        i++;
    }
    if (print) 
        cA->debugPrint(cerr);    
    res = inspectNode(cA, expA[i]);
    if (!res && !print) {return false;}

    if (print) 
        cerr << "Hull B:\n";
    CHNode* cB;
    i = 0;
    for (cB = topB; cB->mBottom != true; cB = cB->mNext->mDst) {
        if (print) 
            cB->debugPrint(cerr);
        res = inspectNode(cB, expB[i]);
        if (!res && !print) {return false;}
        i++;
    }
    if (print) 
        cB->debugPrint(cerr);
    res = inspectNode(cB, expB[i]);
    if (!res && !print) {return false;}

    return true;
}


/**  @brief IntersectionFinderConvexPolygon2D::sweepALBL()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [17.0, 13.0,  3.0,  2.0,  2.0,  3.0, 15.0, 19.0, 22.0]
 * points_ay = [35.0, 36.0, 30.0, 26.0, 11.0,  9.0,  4.0,  4.0,  7.0]
 * points_bx = [15.0, 10.0,  8.0,  4.0,  2.0,  1.0,  1.0,  2.0,  6.0, 11.0, 21.0, 24.0]
 * points_by = [35.0, 35.0, 34.0, 30.0, 26.0, 21.0, 15.0, 11.0,  7.0,  5.0,  5.0,  6.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test42) {

    double ptsA[] = {
        17.0,    35.0,
        13.0,    36.0,
         3.0,    30.0,
         2.0,    26.0,
         2.0,    11.0,
         3.0,     9.0,
        15.0,     4.0,
        19.0,     4.0,
        22.0,     7.0
    };

    double ptsB[] = {
        15.0,    35.0,
        10.0,    35.0,
         8.0,    34.0,
         4.0,    30.0,
         2.0,    26.0,
         1.0,    21.0,
         1.0,    15.0,
         2.0,    11.0,
         6.0,     7.0,
        11.0,     5.0,
        21.0,     5.0,
        24.0,     6.0
    };

    vector<eElem> expA;
    expA.emplace_back(13.0,    36.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back( 0.0,     0.0, TT_INTSEC_PREV_INTERIOR);
    expA.emplace_back( 3.0,    30.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back( 2.0,    26.0, TT_COINCIDENT_VERTEX_NEXT_INTERIOR);
    expA.emplace_back( 2.0,    11.0, TT_COINCIDENT_VERTEX_PREV_INTERIOR);
    expA.emplace_back( 3.0,     9.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back( 0.0,     0.0, TT_INTSEC_NEXT_INTERIOR);
    expA.emplace_back( 0.0,     0.0, TT_INTSEC_PREV_INTERIOR);
    expA.emplace_back(15.0,     4.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(19.0,     4.0, TT_ORIGINAL_VERTEX);

    vector<eElem> expB;
    expB.emplace_back(10.0,    35.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 8.0,    34.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 0.0,     0.0, TT_INTSEC_NEXT_INTERIOR);
    expB.emplace_back( 4.0,    30.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 2.0,    26.0, TT_COINCIDENT_VERTEX_PREV_INTERIOR);
    expB.emplace_back( 1.0,    21.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 1.0,    15.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 2.0,    11.0, TT_COINCIDENT_VERTEX_NEXT_INTERIOR);
    expB.emplace_back( 0.0,     0.0, TT_INTSEC_PREV_INTERIOR);
    expB.emplace_back( 6.0,     7.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(11.0,     5.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 0.0,     0.0, TT_INTSEC_NEXT_INTERIOR);
    expB.emplace_back(21.0,     5.0, TT_ORIGINAL_VERTEX);

    CHNode *topA;
    CHNode *bottomA;
    CHNode *topB;
    CHNode *bottomB;

    createHalfHull(
               ptsA, sizeof(ptsA)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
               ptsB, sizeof(ptsB)/sizeof(double)/2, &topB, &bottomB, false);

    mFinder.sweepALBL(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);

    auto res = inspectHullsLL(topA, topB, expA, expB, false);
    EXPECT_EQ(res, true);

    createHalfHull(
              ptsB, sizeof(ptsB)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
              ptsA, sizeof(ptsA)/sizeof(double)/2, &topB, &bottomB, false);

    mFinder.sweepALBL(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);

    res = inspectHullsLL(topA, topB, expB, expA, false);
    EXPECT_EQ(res, true);

}


/**  @brief IntersectionFinderConvexPolygon2D::sweepALBL()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [11.0,  8.0,  8.0,  8.0,  8.0, 10.0]
 * points_ay = [28.0, 25.0, 24.0, 22.0, 20.0, 18.0]
 * points_bx = [ 9.0,  8.0,  8.0,  8.0,  8.0,  9.0]
 * points_by = [28.0, 27.0, 23.0, 21.0, 18.0, 17.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test43) {

    double ptsA[] = {
        11.0,    28.0,
         8.0,    25.0,
         8.0,    24.0,
         8.0,    22.0,
         8.0,    20.0,
        10.0,    18.0
    };

    double ptsB[] = {
         9.0,    28.0,
         8.0,    27.0,
         8.0,    23.0,
         8.0,    21.0,
         8.0,    18.0,
         9.0,    17.0
    };

    vector<eElem> expA;
    expA.emplace_back(11.0,    28.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back( 8.0,    25.0, TT_ORIGINAL_VERTEX_ON_EDGE_PREV_INTERIOR);
    expA.emplace_back( 8.0,    24.0, TT_ORIGINAL_VERTEX_ON_EDGE);
    expA.emplace_back( 8.0,    23.0, TT_SPLIT_POINT);
    expA.emplace_back( 8.0,    22.0, TT_ORIGINAL_VERTEX_ON_EDGE);
    expA.emplace_back( 8.0,    21.0, TT_SPLIT_POINT);
    expA.emplace_back( 8.0,    20.0, TT_ORIGINAL_VERTEX_ON_EDGE_NEXT_INTERIOR);
    expA.emplace_back(10.0,    18.0, TT_ORIGINAL_VERTEX);

    vector<eElem> expB;
    expB.emplace_back( 9.0,    28.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 8.0,    27.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 8.0,    25.0, TT_SPLIT_POINT);
    expB.emplace_back( 8.0,    24.0, TT_SPLIT_POINT);
    expB.emplace_back( 8.0,    23.0, TT_ORIGINAL_VERTEX_ON_EDGE);
    expB.emplace_back( 8.0,    22.0, TT_SPLIT_POINT);
    expB.emplace_back( 8.0,    21.0, TT_ORIGINAL_VERTEX_ON_EDGE);
    expB.emplace_back( 8.0,    20.0, TT_SPLIT_POINT);
    expB.emplace_back( 8.0,    18.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 9.0,    17.0, TT_ORIGINAL_VERTEX);

    CHNode *topA;
    CHNode *bottomA;
    CHNode *topB;
    CHNode *bottomB;

    createHalfHull(
             ptsA, sizeof(ptsA)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
             ptsB, sizeof(ptsB)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepALBL(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    auto res = inspectHullsLL(topA, topB, expA, expB, false);
    EXPECT_EQ(res, true);

    createHalfHull(
             ptsB, sizeof(ptsB)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
             ptsA, sizeof(ptsA)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepALBL(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    res = inspectHullsLL(topA, topB, expB, expA, false);
    EXPECT_EQ(res, true);

}


/**  @brief IntersectionFinderConvexPolygon2D::sweepALBL()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [21.0, 17.0, 15.0, 21.0]
 * points_ay = [28.0, 27.0, 18.0, 15.0]
 * points_bx = [19.0, 15.0, 15.0, 19.0]
 * points_by = [29.0, 25.0, 18.0, 14.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test44) {

    double ptsA[] = {
        21.0,    28.0,
        17.0,    27.0,
        15.0,    18.0,
        21.0,    15.0
    };

    double ptsB[] = {
        19.0,    29.0,
        15.0,    25.0,
        15.0,    18.0,
        19.0,    14.0
    };

    vector<eElem> expA;
    expA.emplace_back(21.0,    28.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(17.0,    27.0, TT_ORIGINAL_VERTEX_ON_EDGE_BOTH_INTERIOR);
    expA.emplace_back(15.0,    18.0, TT_COINCIDENT_VERTEX_BOTH_INTERIOR);
    expA.emplace_back(21.0,    15.0, TT_ORIGINAL_VERTEX);

    vector<eElem> expB;
    expB.emplace_back(19.0,    29.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(17.0,    27.0, TT_SPLIT_POINT);
    expB.emplace_back(15.0,    25.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(15.0,    18.0, TT_COINCIDENT_VERTEX);
    expB.emplace_back(19.0,    14.0, TT_ORIGINAL_VERTEX);

    CHNode *topA;
    CHNode *bottomA;
    CHNode *topB;
    CHNode *bottomB;

    createHalfHull(
             ptsA, sizeof(ptsA)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
             ptsB, sizeof(ptsB)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepALBL(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    auto res = inspectHullsLL(topA, topB, expA, expB, false);
    EXPECT_EQ(res, true);

    createHalfHull(
              ptsB, sizeof(ptsB)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
              ptsA, sizeof(ptsA)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepALBL(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    res = inspectHullsLL(topA, topB, expB, expA, false);
    EXPECT_EQ(res, true);

}


/**  @brief IntersectionFinderConvexPolygon2D::sweepALBR()& 
 *          IntersectionFinderConvexPolygon2D::sweepARBL()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [14.0, 12.0, 10.0, 10.0, 14.0, 18.0]
 * points_ay = [24.0, 22.0, 18.0, 10.0,  6.0,  4.0]
 * points_bx = [ 2.0,  6.0, 10.0, 10.0,  8.0,  6.0]
 * points_by = [ 4.0,  6.0, 10.0, 18.0, 22.0, 24.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test45) {

    double ptsA[] = {
        14.0,    24.0,
        12.0,    22.0,
        10.0,    18.0,
        10.0,    10.0,
        14.0,     6.0,
        18.0,     4.0,
    };

    double ptsB[] = {
         2.0,     4.0,
         6.0,     6.0,
        10.0,    10.0,
        10.0,    18.0,
         8.0,    22.0,
         6.0,    24.0
    };

    vector<eElem> expA;
    expA.emplace_back(14.0,    24.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(12.0,    22.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(10.0,    18.0, TT_COINCIDENT_VERTEX);
    expA.emplace_back(10.0,    10.0, TT_COINCIDENT_VERTEX);
    expA.emplace_back(14.0,     6.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(18.0,     4.0, TT_ORIGINAL_VERTEX);

    vector<eElem> expB;
    expB.emplace_back( 6.0,    24.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 8.0,    22.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(10.0,    18.0, TT_COINCIDENT_VERTEX);
    expB.emplace_back(10.0,    10.0, TT_COINCIDENT_VERTEX);
    expB.emplace_back( 6.0,     6.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 2.0,     4.0, TT_ORIGINAL_VERTEX);

    CHNode *topA;
    CHNode *bottomA;
    CHNode *topB;
    CHNode *bottomB;

    createHalfHull(
           ptsA, sizeof(ptsA)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
           ptsB, sizeof(ptsB)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepALBR(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    auto res = inspectHullsLR(topA, topB, expA, expB, false);
    EXPECT_EQ(res, true);

    createHalfHull(
             ptsB, sizeof(ptsB)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
             ptsA, sizeof(ptsA)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepARBL(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);

    res = inspectHullsRL(topA, topB, expB, expA, false);
    EXPECT_EQ(res, true);

}


/** @brief IntersectionFinderConvexPolygon2D::sweepALBR()　& 
 *         IntersectionFinderConvexPolygon2D::sweepARBL()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [14.0, 12.0, 10.0, 10.0, 14.0, 18.0]
 * points_ay = [24.0, 22.0, 19.0,  9.0,  6.0,  4.0]
 * points_bx = [2.0,  6.0, 10.0, 10.0,  8.0,  6.0]
 * points_by = [4.0,  6.0, 10.0, 18.0, 22.0, 24.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test46) {

    double ptsA[] = {
        14.0,    24.0,
        12.0,    22.0,
        10.0,    19.0,
        10.0,     9.0,
        14.0,     6.0,
        18.0,     4.0,
    };

    double ptsB[] = {
         2.0,     4.0,
         6.0,     6.0,
        10.0,    10.0,
        10.0,    18.0,
         8.0,    22.0,
         6.0,    24.0
    };

    vector<eElem> expA;
    expA.emplace_back(14.0,    24.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(12.0,    22.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(10.0,    19.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(10.0,    18.0, TT_SPLIT_POINT);
    expA.emplace_back(10.0,    10.0, TT_SPLIT_POINT);
    expA.emplace_back(10.0,     9.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(14.0,     6.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(18.0,     4.0, TT_ORIGINAL_VERTEX);

    vector<eElem> expB;
    expB.emplace_back( 6.0,    24.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 8.0,    22.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(10.0,    18.0, TT_ORIGINAL_VERTEX_ON_EDGE);
    expB.emplace_back(10.0,    10.0, TT_ORIGINAL_VERTEX_ON_EDGE);
    expB.emplace_back( 6.0,     6.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 2.0,     4.0, TT_ORIGINAL_VERTEX);

    CHNode *topA;
    CHNode *bottomA;
    CHNode *topB;
    CHNode *bottomB;

    createHalfHull(
             ptsA, sizeof(ptsA)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
             ptsB, sizeof(ptsB)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepALBR(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    auto res = inspectHullsLR(topA, topB, expA, expB, false);
    EXPECT_EQ(res, true);

    createHalfHull(
               ptsB, sizeof(ptsB)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
               ptsA, sizeof(ptsA)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepARBL(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);

    res = inspectHullsRL(topA, topB, expB, expA, false);
    EXPECT_EQ(res, true);

}

/**  @brief IntersectionFinderConvexPolygon2D::sweepALBR()& 
 *          IntersectionFinderConvexPolygon2D::sweepARBL()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [14.0, 12.0, 10.0, 10.0, 14.0, 18.0]
 * points_ay = [24.0, 22.0, 17.0, 11.0,  6.0,  4.0]
 * points_bx = [ 2.0,  6.0, 10.0, 10.0,  8.0,  6.0]
 * points_by = [ 4.0,  6.0, 10.0, 18.0, 22.0, 24.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test47) {

    double ptsA[] = {
        14.0,    24.0,
        12.0,    22.0,
        10.0,    17.0,
        10.0,    11.0,
        14.0,     6.0,
        18.0,     4.0,
    };

    double ptsB[] = {
         2.0,     4.0,
         6.0,     6.0,
        10.0,    10.0,
        10.0,    18.0,
         8.0,    22.0,
         6.0,    24.0
    };

    vector<eElem> expA;
    expA.emplace_back(14.0,    24.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(12.0,    22.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(10.0,    17.0, TT_ORIGINAL_VERTEX_ON_EDGE);
    expA.emplace_back(10.0,    11.0, TT_ORIGINAL_VERTEX_ON_EDGE);
    expA.emplace_back(14.0,     6.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(18.0,     4.0, TT_ORIGINAL_VERTEX);

    vector<eElem> expB;
    expB.emplace_back( 6.0,    24.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 8.0,    22.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(10.0,    18.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(10.0,    17.0, TT_SPLIT_POINT);
    expB.emplace_back(10.0,    11.0, TT_SPLIT_POINT);
    expB.emplace_back(10.0,    10.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 6.0,     6.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 2.0,     4.0, TT_ORIGINAL_VERTEX);

    CHNode *topA;
    CHNode *bottomA;
    CHNode *topB;
    CHNode *bottomB;

    createHalfHull(
           ptsA, sizeof(ptsA)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
           ptsB, sizeof(ptsB)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepALBR(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    auto res = inspectHullsLR(topA, topB, expA, expB, false);
    EXPECT_EQ(res, true);

    createHalfHull(
              ptsB, sizeof(ptsB)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
              ptsA, sizeof(ptsA)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepARBL(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    res = inspectHullsRL(topA, topB, expB, expA, false);
    EXPECT_EQ(res, true);

}

/**  @brief IntersectionFinderConvexPolygon2D::sweepALBR()& 
 *          IntersectionFinderConvexPolygon2D::sweepARBL()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [14.0, 12.0, 10.0, 10.0, 14.0, 18.0]
 * points_ay = [24.0, 22.0, 19.0, 11.0,  6.0,  4.0]
 * points_bx = [ 2.0,  6.0, 10.0, 10.0,  8.0,  6.0]
 * points_by = [ 4.0,  6.0, 10.0, 18.0, 22.0, 24.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test48) {

    double ptsA[] = {
        14.0,    24.0,
        12.0,    22.0,
        10.0,    19.0,
        10.0,    11.0,
        14.0,     6.0,
        18.0,     4.0,
    };

    double ptsB[] = {
         2.0,     4.0,
         6.0,     6.0,
        10.0,    10.0,
        10.0,    18.0,
         8.0,    22.0,
         6.0,    24.0
    };

    vector<eElem> expA;
    expA.emplace_back(14.0,    24.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(12.0,    22.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(10.0,    19.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(10.0,    18.0, TT_SPLIT_POINT);
    expA.emplace_back(10.0,    11.0, TT_ORIGINAL_VERTEX_ON_EDGE);
    expA.emplace_back(14.0,     6.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(18.0,     4.0, TT_ORIGINAL_VERTEX);

    vector<eElem> expB;
    expB.emplace_back( 6.0,    24.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 8.0,    22.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(10.0,    18.0, TT_ORIGINAL_VERTEX_ON_EDGE);
    expB.emplace_back(10.0,    11.0, TT_SPLIT_POINT);
    expB.emplace_back(10.0,    10.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 6.0,     6.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 2.0,     4.0, TT_ORIGINAL_VERTEX);

    CHNode *topA;
    CHNode *bottomA;
    CHNode *topB;
    CHNode *bottomB;

    createHalfHull(
             ptsA, sizeof(ptsA)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
             ptsB, sizeof(ptsB)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepALBR(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    auto res = inspectHullsLR(topA, topB, expA, expB, false);
    EXPECT_EQ(res, true);

    createHalfHull(
             ptsB, sizeof(ptsB)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
             ptsA, sizeof(ptsA)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepARBL(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);

    res = inspectHullsRL(topA, topB, expB, expA, false);
    EXPECT_EQ(res, true);

}

/**  @brief IntersectionFinderConvexPolygon2D::sweepALBR() &
 *          IntersectionFinderConvexPolygon2D::sweepARBL()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [14.0, 12.0, 10.0, 14.0, 18.0]
 * points_ay = [24.0, 22.0, 14.0,  6.0,  4.0]
 * points_bx = [ 2.0,  6.0, 10.0, 10.0,  8.0,  6.0]
 * points_by = [ 4.0,  6.0, 10.0, 18.0, 22.0, 24.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test49) {

    double ptsA[] = {
        14.0,    24.0,
        12.0,    22.0,
        10.0,    14.0,
        14.0,     6.0,
        18.0,     4.0,
    };

    double ptsB[] = {
         2.0,     4.0,
         6.0,     6.0,
        10.0,    10.0,
        10.0,    18.0,
         8.0,    22.0,
         6.0,    24.0
    };

    vector<eElem> expA;
    expA.emplace_back(14.0,    24.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(12.0,    22.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(10.0,    14.0, TT_TOUCHING_POINT_ORIGINAL_VERTEX);
    expA.emplace_back(14.0,     6.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(18.0,     4.0, TT_ORIGINAL_VERTEX);

    vector<eElem> expB;
    expB.emplace_back( 6.0,    24.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 8.0,    22.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(10.0,    18.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(10.0,    14.0, TT_TOUCHING_POINT_SPLIT_POINT);
    expB.emplace_back(10.0,    10.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 6.0,     6.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 2.0,     4.0, TT_ORIGINAL_VERTEX);

    CHNode *topA;
    CHNode *bottomA;
    CHNode *topB;
    CHNode *bottomB;

    createHalfHull(
             ptsA, sizeof(ptsA)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
             ptsB, sizeof(ptsB)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepALBR(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    auto res = inspectHullsLR(topA, topB, expA, expB, false);
    EXPECT_EQ(res, true);

    createHalfHull(
           ptsB, sizeof(ptsB)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
           ptsA, sizeof(ptsA)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepARBL(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    res = inspectHullsRL(topA, topB, expB, expA, false);
    EXPECT_EQ(res, true);

}


/**  @brief IntersectionFinderConvexPolygon2D::sweepALBR() &
 *          IntersectionFinderConvexPolygon2D::sweepARBL()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [14.0, 12.0, 10.0, 14.0, 18.0]
 * points_ay = [24.0, 22.0, 14.0,  6.0,  4.0]
 * points_bx = [ 2.0,  6.0, 10.0,  8.0,  6.0]
 * points_by = [ 4.0,  6.0, 14.0, 22.0, 24.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test50) {

    double ptsA[] = {
        14.0,    24.0,
        12.0,    22.0,
        10.0,    14.0,
        14.0,     6.0,
        18.0,     4.0,
    };

    double ptsB[] = {
         2.0,     4.0,
         6.0,     6.0,
        10.0,    14.0,
         8.0,    22.0,
         6.0,    24.0
    };

    vector<eElem> expA;
    expA.emplace_back(14.0,    24.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(12.0,    22.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(10.0,    14.0, TT_TOUCHING_POINT_COINCIDENT_VERTEX);
    expA.emplace_back(14.0,     6.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(18.0,     4.0, TT_ORIGINAL_VERTEX);

    vector<eElem> expB;
    expB.emplace_back( 6.0,    24.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 8.0,    22.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(10.0,    14.0, TT_TOUCHING_POINT_COINCIDENT_VERTEX);
    expB.emplace_back( 6.0,     6.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 2.0,     4.0, TT_ORIGINAL_VERTEX);

    CHNode *topA;
    CHNode *bottomA;
    CHNode *topB;
    CHNode *bottomB;

    createHalfHull(
            ptsA, sizeof(ptsA)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
            ptsB, sizeof(ptsB)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepALBR(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    auto res = inspectHullsLR(topA, topB, expA, expB, false);
    EXPECT_EQ(res, true);

    createHalfHull(
               ptsB, sizeof(ptsB)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
               ptsA, sizeof(ptsA)/sizeof(double)/2, &topB, &bottomB, false);
    mFinder.sweepARBL(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    res = inspectHullsRL(topA, topB, expB, expA, false);
    EXPECT_EQ(res, true);

}


/**  @brief IntersectionFinderConvexPolygon2D::sweepARBR()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [ 4.0, 10.0, 14.0, 16.0, 13.0, 10.0,  5.0]
 * points_ay = [ 4.0,  6.0, 10.0, 17.0, 26.0, 29.0, 30.0]
 * points_bx = [ 6.0, 9.0, 15.0, 16.0, 15.0,  7.0]
 * points_by = [ 5.0, 7.0, 13.0, 20.0, 24.0, 32.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test51) {

    double ptsA[] = {
         4.0,     4.0,
        10.0,     6.0,
        14.0,    10.0,
        16.0,    17.0,
        13.0,    26.0,
        10.0,    29.0,
         5.0,    30.0
    };

    double ptsB[] = {
         6.0,     5.0,
         9.0,     7.0,
        15.0,    13.0,
        16.0,    20.0,
        15.0,    24.0,
         7.0,    32.0
    };

    vector<eElem> expA;
    expA.emplace_back( 4.0,     4.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(10.0,     6.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back(14.0,    10.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back( 0.0,     0.0, TT_INTSEC_NEXT_INTERIOR);
    expA.emplace_back( 0.0,     0.0, TT_INTSEC_PREV_INTERIOR);
    expA.emplace_back(16.0,    17.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back( 0.0,     0.0, TT_INTSEC_NEXT_INTERIOR);
    expA.emplace_back(13.0,    26.0, TT_ORIGINAL_VERTEX_ON_EDGE_PREV_INTERIOR);
    expA.emplace_back(10.0,    29.0, TT_ORIGINAL_VERTEX_ON_EDGE_NEXT_INTERIOR);
    expA.emplace_back( 5.0,    30.0, TT_ORIGINAL_VERTEX);

    vector<eElem> expB;
    expB.emplace_back( 6.0,     5.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 9.0,     7.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 0.0,     0.0, TT_INTSEC_PREV_INTERIOR);
    expB.emplace_back(15.0,    13.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 0.0,     0.0, TT_INTSEC_NEXT_INTERIOR);
    expB.emplace_back( 0.0,     0.0, TT_INTSEC_PREV_INTERIOR);
    expB.emplace_back(16.0,    20.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(15.0,    24.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back(13.0,    26.0, TT_SPLIT_POINT);
    expB.emplace_back(10.0,    29.0, TT_SPLIT_POINT);
    expB.emplace_back( 7.0,    32.0, TT_ORIGINAL_VERTEX);

    CHNode *topA;
    CHNode *bottomA;
    CHNode *topB;
    CHNode *bottomB;

    createHalfHull(
           ptsA, sizeof(ptsA)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
           ptsB, sizeof(ptsB)/sizeof(double)/2, &topB, &bottomB, false);

    mFinder.sweepARBR(bottomA, bottomB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    auto res = inspectHullsRR(bottomA, bottomB, expA, expB, false);
    EXPECT_EQ(res, true);

    createHalfHull(
               ptsB, sizeof(ptsB)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
               ptsA, sizeof(ptsA)/sizeof(double)/2, &topB, &bottomB, false);

    mFinder.sweepARBR(bottomA, bottomB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);
    res = inspectHullsRR(bottomA, bottomB, expB, expA, false);
    EXPECT_EQ(res, true);

}


/**  @brief IntersectionFinderConvexPolygon2D::sweepARBR()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [ 6.0,  9.0,  9.0,  9.0,  7.0]
 * points_ay = [ 9.0, 12.0, 17.0, 20.0, 22.0]
 * points_bx = [ 6.0,  9.0,  9.0,  9.0,  6.0]
 * points_by = [11.0, 14.0, 17.0, 22.0, 25.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test52) {

    double ptsA[] = {
         6.0,     9.0,
         9.0,    12.0,
         9.0,    17.0,
         9.0,    20.0,
         7.0,    22.0,
    };

    double ptsB[] = {
         6.0,    11.0,
         9.0,    14.0,
         9.0,    17.0,
         9.0,    22.0,
         6.0,    25.0
    };

    vector<eElem> expA;
    expA.emplace_back( 6.0,     9.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back( 9.0,    12.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back( 9.0,    14.0, TT_SPLIT_POINT);
    expA.emplace_back( 9.0,    17.0, TT_COINCIDENT_VERTEX);
    expA.emplace_back( 9.0,    20.0, TT_ORIGINAL_VERTEX_ON_EDGE_NEXT_INTERIOR);
    expA.emplace_back( 7.0,    22.0, TT_ORIGINAL_VERTEX);

    vector<eElem> expB;
    expB.emplace_back( 6.0,    11.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 9.0,    14.0, TT_ORIGINAL_VERTEX_ON_EDGE_PREV_INTERIOR);
    expB.emplace_back( 9.0,    17.0, TT_COINCIDENT_VERTEX);
    expB.emplace_back( 9.0,    20.0, TT_SPLIT_POINT);
    expB.emplace_back( 9.0,    22.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 6.0,    25.0, TT_ORIGINAL_VERTEX);

    CHNode *topA;
    CHNode *bottomA;
    CHNode *topB;
    CHNode *bottomB;

    createHalfHull(
               ptsA, sizeof(ptsA)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
               ptsB, sizeof(ptsB)/sizeof(double)/2, &topB, &bottomB, false);

    mFinder.sweepARBR(bottomA, bottomB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);

    auto res = inspectHullsRR(bottomA, bottomB, expA, expB, false);
    EXPECT_EQ(res, true);

    createHalfHull(
            ptsB, sizeof(ptsB)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
            ptsA, sizeof(ptsA)/sizeof(double)/2, &topB, &bottomB, false);

    mFinder.sweepARBR(bottomA, bottomB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);

    res = inspectHullsRR(bottomA, bottomB, expB, expA, false);
    EXPECT_EQ(res, true);

}


/**  @brief IntersectionFinderConvexPolygon2D::sweepALBL()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [ 3.0,  5.0,  3.0,  1.0]
 * points_ay = [14.0, 16.0, 18.0, 16.0]
 * points_bx = [ 3.0,  5.0,  3.0,  1.0]
 * points_by = [18.0, 20.0, 22.0, 20.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test53) {

    double ptsA[] = {
         3.0,    14.0,
         5.0,    16.0,
         3.0,    18.0,
         1.0,    16.0
    };

    double ptsB[] = {
         3.0,    18.0,
         5.0,    20.0,
         3.0,    22.0,
         1.0,    20.0
    };

    vector<eElem> expA;
    expA.emplace_back( 3.0,    14.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back( 5.0,    16.0, TT_ORIGINAL_VERTEX);
    expA.emplace_back( 3.0,    18.0, TT_ORIGINAL_VERTEX);

    vector<eElem> expB;
    expB.emplace_back( 3.0,    18.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 5.0,    20.0, TT_ORIGINAL_VERTEX);
    expB.emplace_back( 3.0,    22.0, TT_ORIGINAL_VERTEX);

    CHNode *topA;
    CHNode *bottomA;
    CHNode *topB;
    CHNode *bottomB;

    createHalfHull(
          ptsA, sizeof(ptsA)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
          ptsB, sizeof(ptsB)/sizeof(double)/2, &topB, &bottomB, false);

    mFinder.sweepARBR(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);

    auto res = inspectHullsRR(bottomA, bottomB, expA, expB, false);
    EXPECT_EQ(res, true);

    createHalfHull(
            ptsB, sizeof(ptsB)/sizeof(double)/2, &topA, &bottomA, true);
    createHalfHull(
            ptsA, sizeof(ptsA)/sizeof(double)/2, &topB, &bottomB, false);

    mFinder.sweepARBR(topA, topB, sizeof(ptsA)/sizeof(double)/2, sizeof(ptsB)/sizeof(double)/2);

    res = inspectHullsRR(bottomA, bottomB, expB, expA, false);
    EXPECT_EQ(res, true);
}


/**  @brief IntersectionFinderConvexPolygon2D::
 *          findIntersectionOfPolygonAndPolygon()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [2.0, 7.0, 7.0, 2.0, 2.0]
 * points_ay = [2.0, 2.0, 7.0, 7.0, 2.0]
 * points_bx = [4.0, 9.0, 7.0, 4.0, 1.0, 4.0]
 * points_by = [2.0, 7.0, 8.0, 9.0, 4.0, 2.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test54) {

    vector<Vec2> Ain;
    vector<Vec2> Bin;

    Ain.emplace_back(2.0, 2.0);
    Ain.emplace_back(7.0, 2.0);
    Ain.emplace_back(7.0, 7.0);
    Ain.emplace_back(2.0, 7.0);

    Bin.emplace_back(4.0, 2.0);
    Bin.emplace_back(9.0, 7.0);
    Bin.emplace_back(7.0, 8.0);
    Bin.emplace_back(4.0, 9.0);
    Bin.emplace_back(1.0, 4.0);

    IntersectionFinderConvexPolygon2D::InputElem ie_a1(Ain[0], 1);
    IntersectionFinderConvexPolygon2D::InputElem ie_a2(Ain[1], 2);
    IntersectionFinderConvexPolygon2D::InputElem ie_a3(Ain[2], 3);
    IntersectionFinderConvexPolygon2D::InputElem ie_a4(Ain[3], 4);

    IntersectionFinderConvexPolygon2D::InputElem ie_b1(Bin[0], 11);
    IntersectionFinderConvexPolygon2D::InputElem ie_b2(Bin[1], 12);
    IntersectionFinderConvexPolygon2D::InputElem ie_b3(Bin[2], 13);
    IntersectionFinderConvexPolygon2D::InputElem ie_b4(Bin[3], 14);
    IntersectionFinderConvexPolygon2D::InputElem ie_b5(Bin[4], 15);

    vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
    input_a01.push_back(&ie_a1);
    input_a01.push_back(&ie_a2);
    input_a01.push_back(&ie_a3);
    input_a01.push_back(&ie_a4);

    vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
    input_b01.push_back(&ie_b1);
    input_b01.push_back(&ie_b2);
    input_b01.push_back(&ie_b3);
    input_b01.push_back(&ie_b4);
    input_b01.push_back(&ie_b5);


    vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

    IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               Vec2(2.0,5.66667), 4, 1,  14, 15, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               Vec2(2.0,3.33333), 4, 1,  15, 11, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               Vec2(4.0, 2.0),    1, 2,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

    IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               Vec2(7.0, 5.0),    2, 3,  11,  12, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_05(
               Vec2(7.0, 7.0),    3,-1, -1, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

    IntersectionFinderConvexPolygon2D::OutputElem oe_06(
               Vec2(2.8, 7.0),    3, 4,  14,  15, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    expected_01.push_back(oe_01);
    expected_01.push_back(oe_02);
    expected_01.push_back(oe_03);
    expected_01.push_back(oe_04);
    expected_01.push_back(oe_05);
    expected_01.push_back(oe_06);

    vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;
    auto res = mFinder.
            findIntersection_polygon_polygon(input_a01, input_b01, output_01);
    EXPECT_EQ(res, true);
//    for (auto& v : output_01) {
//        cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//        cerr << "mIndexA: " << v.mIndexA << "\n";
//        cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//        cerr << "mIndexB: " << v.mIndexB << "\n";
//        cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//    }
    EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


    vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

    IntersectionFinderConvexPolygon2D::OutputElem oe_07(
               Vec2(2.0,5.66667), 14, 15, 4, 1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_08(
               Vec2(2.0,3.33333), 15, 11, 4, 1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_09(
               Vec2(4.0, 2.0),   11, -1,  1, 2,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_10(
               Vec2(7.0, 5.0),   11, 12,  2, 3,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_11(
               Vec2(7.0, 7.0),  -1, -1,  3,-1,
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

    IntersectionFinderConvexPolygon2D::OutputElem oe_12(
               Vec2(2.8, 7.0),  14, 15,  3, 4,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    expected_02.push_back(oe_07);
    expected_02.push_back(oe_08);
    expected_02.push_back(oe_09);
    expected_02.push_back(oe_10);
    expected_02.push_back(oe_11);
    expected_02.push_back(oe_12);

    vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;
    res = mFinder.
            findIntersection_polygon_polygon(input_b01, input_a01, output_02);
    EXPECT_EQ(res, true);
//    for (auto& v : output_02) {
//        cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//        cerr << "mIndexA: " << v.mIndexA << "\n";
//        cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//        cerr << "mIndexB: " << v.mIndexB << "\n";
//        cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//    }
    EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);


}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [2.0, 7.0, 7.0, 2.0, 2.0]
 * points_ay = [2.0, 2.0, 7.0, 7.0, 2.0]
 * points_bx = [4.0, 9.0, 7.0, 4.0, 1.0, 4.0]
 * points_by = [2.0, 7.0, 8.0, 9.0, 4.0, 2.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
*/
TEST_F(IntersectionFinderConvexPolygon2DTests, Test55) {

    vector<Vec2> Ain;
    vector<Vec2> Bin;

    Ain.emplace_back(2.0, 2.0);
    Ain.emplace_back(7.0, 2.0);
    Ain.emplace_back(7.0, 7.0);
    Ain.emplace_back(2.0, 7.0);

    Bin.emplace_back(4.0, 2.0);
    Bin.emplace_back(9.0, 7.0);
    Bin.emplace_back(7.0, 8.0);
    Bin.emplace_back(4.0, 9.0);
    Bin.emplace_back(1.0, 4.0);

    vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

    IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               Vec2(2.0,5.66667), 4, 1,  14, 15, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               Vec2(2.0,3.33333), 4, 1,  15, 11, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               Vec2(4.0, 2.0),    1, 2,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

    IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               Vec2(7.0, 5.0),    2, 3,  11,  12, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_05(
               Vec2(7.0, 7.0),    3,-1, -1, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

    IntersectionFinderConvexPolygon2D::OutputElem oe_06(
               Vec2(2.8, 7.0),    3, 4,  14,  15, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    expected_01.push_back(oe_01);
    expected_01.push_back(oe_02);
    expected_01.push_back(oe_03);
    expected_01.push_back(oe_04);
    expected_01.push_back(oe_05);
    expected_01.push_back(oe_06);

    IntersectionFinderConvexPolygon2D::InputElem ie_a1(Ain[0], 1);
    IntersectionFinderConvexPolygon2D::InputElem ie_a2(Ain[1], 2);
    IntersectionFinderConvexPolygon2D::InputElem ie_a3(Ain[2], 3);
    IntersectionFinderConvexPolygon2D::InputElem ie_a4(Ain[3], 4);

    IntersectionFinderConvexPolygon2D::InputElem ie_b1(Bin[0], 11);
    IntersectionFinderConvexPolygon2D::InputElem ie_b2(Bin[1], 12);
    IntersectionFinderConvexPolygon2D::InputElem ie_b3(Bin[2], 13);
    IntersectionFinderConvexPolygon2D::InputElem ie_b4(Bin[3], 14);
    IntersectionFinderConvexPolygon2D::InputElem ie_b5(Bin[4], 15);

    vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
    input_a01.push_back(&ie_a1);
    input_a01.push_back(&ie_a2);
    input_a01.push_back(&ie_a3);
    input_a01.push_back(&ie_a4);

    vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
    input_b01.push_back(&ie_b1);
    input_b01.push_back(&ie_b2);
    input_b01.push_back(&ie_b3);
    input_b01.push_back(&ie_b4);
    input_b01.push_back(&ie_b5);

    vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;
    auto res = mFinder.
            findIntersection_polygon_polygon(input_a01, input_b01, output_01);
    EXPECT_EQ(res, true);
//    for (auto& v : output_01) {
//        cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//        cerr << "mIndexA: " << v.mIndexA << "\n";
//        cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//        cerr << "mIndexB: " << v.mIndexB << "\n";
//        cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//    }
    EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


    vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

    IntersectionFinderConvexPolygon2D::OutputElem oe_07(
               Vec2(2.0,5.66667), 14, 15, 4, 1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_08(
               Vec2(2.0,3.33333), 15, 11, 4, 1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_09(
               Vec2(4.0, 2.0),    11, -1, 1, 2,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_10(
               Vec2(7.0, 5.0),    11, 12, 2, 3,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    IntersectionFinderConvexPolygon2D::OutputElem oe_11(
               Vec2(7.0, 7.0),    -1, -1, 3,-1,
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

    IntersectionFinderConvexPolygon2D::OutputElem oe_12(
               Vec2(2.8, 7.0),   14, 15,  3, 4,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

    expected_02.push_back(oe_07);
    expected_02.push_back(oe_08);
    expected_02.push_back(oe_09);
    expected_02.push_back(oe_10);
    expected_02.push_back(oe_11);
    expected_02.push_back(oe_12);

    vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;
    res = mFinder.
            findIntersection_polygon_polygon(input_b01, input_a01, output_02);
    EXPECT_EQ(res, true);
//    for (auto& v : output_01) {
//        cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//        cerr << "mIndexA: " << v.mIndexA << "\n";
//        cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//        cerr << "mIndexB: " << v.mIndexB << "\n";
//        cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//    }
    EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [ 3.0,  5.0,  3.0,  1.0,  3.0]
 * points_ay = [14.0, 16.0, 18.0, 16.0, 14.0]
 * points_bx = [ 3.0,  5.0,  3.0,  1.0,  3.0]
 * points_by = [18.0, 20.0, 22.0, 20.0, 18.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test56) {

    vector<Vec2> Ain;
    vector<Vec2> Bin;

    Ain.emplace_back(3.0, 14.0);
    Ain.emplace_back(5.0, 16.0);
    Ain.emplace_back(3.0, 18.0);
    Ain.emplace_back(1.0, 16.0);

    Bin.emplace_back(3.0, 18.0);
    Bin.emplace_back(5.0, 20.0);
    Bin.emplace_back(3.0, 22.0);
    Bin.emplace_back(1.0, 20.0);

    vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

    IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               Vec2(3.0, 18.0), 3, -1,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

    expected_01.push_back(oe_01);

    IntersectionFinderConvexPolygon2D::InputElem ie_a1(Ain[0], 1);
    IntersectionFinderConvexPolygon2D::InputElem ie_a2(Ain[1], 2);
    IntersectionFinderConvexPolygon2D::InputElem ie_a3(Ain[2], 3);
    IntersectionFinderConvexPolygon2D::InputElem ie_a4(Ain[3], 4);

    IntersectionFinderConvexPolygon2D::InputElem ie_b1(Bin[0], 11);
    IntersectionFinderConvexPolygon2D::InputElem ie_b2(Bin[1], 12);
    IntersectionFinderConvexPolygon2D::InputElem ie_b3(Bin[2], 13);
    IntersectionFinderConvexPolygon2D::InputElem ie_b4(Bin[3], 14);

    vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
    input_a01.push_back(&ie_a1);
    input_a01.push_back(&ie_a2);
    input_a01.push_back(&ie_a3);
    input_a01.push_back(&ie_a4);

    vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
    input_b01.push_back(&ie_b1);
    input_b01.push_back(&ie_b2);
    input_b01.push_back(&ie_b3);
    input_b01.push_back(&ie_b4);

    vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;
    auto res = mFinder.
            findIntersection_polygon_polygon(input_a01, input_b01, output_01);
    EXPECT_EQ(res, true);
//    for (auto& v : output_01) {
//        cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//        cerr << "mIndexA: " << v.mIndexA << "\n";
//        cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//        cerr << "mIndexB: " << v.mIndexB << "\n";
//        cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//    }
    EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);



    vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

    IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               Vec2(3.0, 18.0), 11, -1,  3, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

    expected_02.push_back(oe_02);

    vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;
    res = mFinder.
            findIntersection_polygon_polygon(input_b01, input_a01, output_02);
    EXPECT_EQ(res, true);
//    for (auto& v : output_02) {
//        cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//        cerr << "mIndexA: " << v.mIndexA << "\n";
//        cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//        cerr << "mIndexB: " << v.mIndexB << "\n";
//        cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//    }
    EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test57) {

    vector<Vec2> Ain;
    vector<Vec2> Bin;
    Ain.emplace_back(3.0, 14.0);

    Bin.emplace_back(3.0, 14.0);

    IntersectionFinderConvexPolygon2D::InputElem ie_a1(Ain[0], 1);

    IntersectionFinderConvexPolygon2D::InputElem ie_b1(Bin[0], 11);

    vector<IntersectionFinderConvexPolygon2D::InputElem*> input_a01;
    input_a01.push_back(&ie_a1);

    vector<IntersectionFinderConvexPolygon2D::InputElem*> input_b01;
    input_b01.push_back(&ie_b1);

    vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

    IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               Vec2(3.0, 14.0),  1, -1,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

    expected_01.push_back(oe_01);

    vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

    auto res = mFinder.
             findIntersection_polygon_polygon(input_a01, input_b01, output_01);
    EXPECT_EQ(res, true);
//    for (auto& v : output_01) {
//        cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//        cerr << "mIndexA: " << v.mIndexA << "\n";
//        cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//        cerr << "mIndexB: " << v.mIndexB << "\n";
//        cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//    }
    EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);

    vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

    IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               Vec2(3.0, 14.0),  11, -1, 1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

    expected_02.push_back(oe_02);

    vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

    res = mFinder.
             findIntersection_polygon_polygon(input_b01, input_a01, output_02);
    EXPECT_EQ(res, true);
//    for (auto& v : output_02) {
//        cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//        cerr << "mIndexA: " << v.mIndexA << "\n";
//        cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//        cerr << "mIndexB: " << v.mIndexB << "\n";
//        cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//    }
    EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test58) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

	Vec2 a1(-3.0, 0.0);
	Vec2 a2( 2.0, 0.0);
        Vec2 b1( 0.0, 0.0);
        Vec2 b2( 4.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               a2,  2, -1,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               b1,  1,  2,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               a2,  11, 12, 2, -1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               b1,  11, -1, 1,  2,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        expected_02.push_back(oe_03);
        expected_02.push_back(oe_04);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 */

TEST_F(IntersectionFinderConvexPolygon2DTests, Test59) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-3.0, 0.0);
        Vec2 a2( 4.0, 0.0);
        Vec2 b1( 0.0, 0.0);
        Vec2 b2( 2.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               b1,  1,  2,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               b2,  1,  2,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               b1,  11, -1, 1,  2,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               b2,  12, -1, 1,  2,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        expected_02.push_back(oe_03);
        expected_02.push_back(oe_04);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);
    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test60) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-3.0, 0.0);
        Vec2 a2( 1.0, 0.0);

        Vec2 b1( 1.0, 0.0);
        Vec2 b2( 2.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               a2,  2, -1,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_01.push_back(oe_01);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               a2,  11, -1, 2, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_02.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);
    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test61) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);


        Vec2 a1(-1.0,-1.0);
        Vec2 a2( 1.0, 1.0);

        Vec2 b1( 1.0,-1.0);
        Vec2 b2(-1.0, 1.0);

        Vec2 c1( 0.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;
        c1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        c1  = rotatedVec(c1.x(), c1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,  1, 2,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        expected_01.push_back(oe_01);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c1,  11, 12, 1, 2,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        expected_02.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test62) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1(-1.0,-1.0);
        Vec2 a2( 1.0, 1.0);

        Vec2 b1(-1.0, 0.0);
        Vec2 b2( 1.0, 2.0);

        Vec2 c1( 0.0, 0.0);

        a1 += trans;
        a2 += trans;
        b1 += trans;
        b2 += trans;
        c1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        c1  = rotatedVec(c1.x(), c1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, false);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }


        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, false);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }
    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test63) {


    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 2.0, 0.0);
        Vec2 a3( 1.0, 2.0);

        Vec2 b1( 0.5, 0.5);

        Vec2 c1( 0.5, 0.5);


        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        c1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        c1  = rotatedVec(c1.x(), c1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,  -1, -1,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        expected_01.push_back(oe_01);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c1,  11, -1, -1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        expected_02.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test64) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);


        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 2.0, 0.0);
        Vec2 a3( 1.0, 2.0);

        Vec2 b1( 0.5, 1.0);

        Vec2 c1( 0.5, 1.0);


        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        c1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        c1  = rotatedVec(c1.x(), c1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,  1, 3,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        expected_01.push_back(oe_01);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c1,  11, -1, 1, 3,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        expected_02.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test65) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 2.0, 0.0);
        Vec2 a3( 1.0, 2.0);
        Vec2 b1(-0.5, 1.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, false);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;
        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, false);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 2.0, 1.0, 0.0]
 * points_ay = [0.0, 0.0, 2.0, 0.0]
 * points_bx = [0.5, 1.5]
 * points_by = [0.5, 0.5]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test66) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

	Vec2 a1( 0.0, 0.0);
        Vec2 a2( 2.0, 0.0);
        Vec2 a3( 1.0, 2.0);

        Vec2 b1( 0.5, 0.5);
        Vec2 b2( 1.5, 0.5);

        Vec2 c1( 0.5, 0.5);
        Vec2 c2( 1.5, 0.5);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        c1 += trans;
        c2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,  -1, -1,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,  -1, -1,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c1,  11, -1, -1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c2,  12, -1, -1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        expected_02.push_back(oe_03);
        expected_02.push_back(oe_04);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 2.0, 1.0, 0.0]
 * points_ay = [0.0, 0.0, 2.0, 0.0]
 * points_bx = [ 1.0, 1.0]
 * points_by = [-0.5, 0.5]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test67) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 2.0, 0.0);
        Vec2 a3( 1.0, 2.0);

        Vec2 b1( 1.0,-0.5);
        Vec2 b2( 1.0, 0.5);

        Vec2 c1( 1.0, 0.0);
        Vec2 c2( 1.0, 0.5);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        c1 += trans;
        c2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1,  2,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,  -1, -1,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);



        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c1,  11, 12,  1,  2,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c2,  12, -1, -1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        expected_02.push_back(oe_03);
        expected_02.push_back(oe_04);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 2.0, 1.0, 0.0]
 * points_ay = [0.0, 0.0, 2.0, 0.0]
 * points_bx = [-1.0, 3.0]
 * points_by = [1.0, 1.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test68) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 2.0, 0.0);
        Vec2 a3( 1.0, 2.0);

        Vec2 b1(-1.0, 1.0);
        Vec2 b2( 3.0, 1.0);

        Vec2 c1( 0.5, 1.0);
        Vec2 c2( 1.5, 1.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        c1 += trans;
        c2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1,  3,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,  2,   3,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c1,  11, 12,  1,  3,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c2,  11, 12, 2,   3,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        expected_02.push_back(oe_03);
        expected_02.push_back(oe_04);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}



/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 2.0, 1.0, 0.0]
 * points_ay = [0.0, 0.0, 2.0, 0.0]
 * points_bx = [0.5, 1.5]
 * points_by = [1.0, 1.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test69) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 2.0, 0.0);
        Vec2 a3( 1.0, 2.0);

        Vec2 b1( 0.5, 1.0);
        Vec2 b2( 1.5, 1.0);

        Vec2 c1( 0.5, 1.0);
        Vec2 c2( 1.5, 1.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        c1 += trans;
        c2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1,  3,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,  2,   3,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);



        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c1,  11, -1,  1,  3,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c2,  12, -1, 2,   3,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        expected_02.push_back(oe_03);
        expected_02.push_back(oe_04);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);
    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 2.0, 1.0, 0.0]
 * points_ay = [0.0, 0.0, 2.0, 0.0]
 * points_bx = [1.0, 3.0]
 * points_by = [0.0, 0.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test70) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 2.0, 0.0);
        Vec2 a3( 1.0, 2.0);

        Vec2 b1( 1.0, 0.0);
        Vec2 b2( 3.0, 0.0);

        Vec2 c1( 1.0, 0.0);
        Vec2 c2( 2.0, 0.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        c1 += trans;
        c2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1,  2,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,  2,   -1,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c1,  11, -1,  1,  2,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c2,  11, 12, 2,   -1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        expected_02.push_back(oe_03);
        expected_02.push_back(oe_04);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);
    }
}



/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 2.0, 1.0, 0.0]
 * points_ay = [0.0, 0.0, 2.0, 0.0]
 * points_bx = [-1.0, 3.0]
 * points_by = [0.0, 0.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test71) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 2.0, 0.0);
        Vec2 a3( 1.0, 2.0);

        Vec2 b1(-1.0, 0.0);
        Vec2 b2( 3.0, 0.0);

        Vec2 c1( 0.0, 0.0);
        Vec2 c2( 2.0, 0.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        c1 += trans;
        c2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1, -1,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,   2, -1,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c1,  11, 12,  1, -1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c2,  11, 12,  2, -1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        expected_02.push_back(oe_03);
        expected_02.push_back(oe_04);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}



/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 4.0, 2.0, 0.0]
 * points_ay = [0.0, 0.0, 4.0, 0.0]
 * points_bx = [0.0, 2.0, 4.0, 0.0]
 * points_by = [2.0, -2.0, 2.0, 2.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test72) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 4.0, 0.0);
        Vec2 a3( 2.0, 4.0);

        Vec2 b1( 0.0, 2.0);
        Vec2 b2( 2.0,-2.0);
        Vec2 b3( 4.0, 2.0);

        Vec2 c1( 1.0, 0.0);
        Vec2 c2( 3.0, 0.0);
        Vec2 c3( 3.5, 1.0);
        Vec2 c4( 3.0, 2.0);
        Vec2 c5( 1.0, 2.0);
        Vec2 c6( 0.5, 1.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;
        c2 += trans;
        c3 += trans;
        c4 += trans;
        c5 += trans;
        c6 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);
        c3  = rotatedVec(c3.x(), c3.y(), i);
        c4  = rotatedVec(c4.x(), c4.y(), i);
        c5  = rotatedVec(c5.x(), c5.y(), i);
        c6  = rotatedVec(c6.x(), c6.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1,  2,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,   1,  2,  12, 13, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c3,   2,  3,  12, 13, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c4,   2,  3,  13, 11, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_05(
               c5,   3,  1,  13, 11, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_06(
               c6,   3,  1,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);
        expected_01.push_back(oe_03);
        expected_01.push_back(oe_04);
        expected_01.push_back(oe_05);
        expected_01.push_back(oe_06);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);



        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_07(
               c1,  11, 12,  1,  2,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_08(
               c2,  12, 13,  1,  2,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_09(
               c3,  12, 13,  2,  3,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_10(
               c4,  13, 11,  2,  3,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_11(
               c5,  13, 11,  3,  1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_12(
               c6,  11, 12,  3,  1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        expected_02.push_back(oe_07);
        expected_02.push_back(oe_08);
        expected_02.push_back(oe_09);
        expected_02.push_back(oe_10);
        expected_02.push_back(oe_11);
        expected_02.push_back(oe_12);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 4.0, 2.0, 0.0]
 * points_ay = [0.0, 0.0, 4.0, 0.0]
 * points_bx = [0.0, 2.0, 4.0, 0.0]
 * points_by = [6.0, 2.0, 6.0, 6.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test73) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 4.0, 0.0);
        Vec2 a3( 2.0, 4.0);

        Vec2 b1( 0.0, 6.0);
        Vec2 b2( 2.0, 2.0);
        Vec2 b3( 4.0, 6.0);

        Vec2 c1( 1.5, 3.0);
        Vec2 c2( 2.0, 2.0);
        Vec2 c3( 2.5, 3.0);
        Vec2 c4( 2.0, 4.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;
        c2 += trans;
        c3 += trans;
        c4 += trans;


        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);
        c3  = rotatedVec(c3.x(), c3.y(), i);
        c4  = rotatedVec(c4.x(), c4.y(), i);


        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1,  3,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,  -1, -1,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c3,   2,  3,  12, 13, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c4,   3, -1,  -1, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);
        expected_01.push_back(oe_03);
        expected_01.push_back(oe_04);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);



        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_05(
               c1,  11, 12,  1,  3,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_06(
               c2,  12, -1, -1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        IntersectionFinderConvexPolygon2D::OutputElem oe_07(
               c3,  12, 13,  2,  3,
               IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_08(
               c4,  -1, -1,  3, -1,
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        expected_02.push_back(oe_05);
        expected_02.push_back(oe_06);
        expected_02.push_back(oe_07);
        expected_02.push_back(oe_08);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);
    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [1.0, 2.0, 4.0, 1.0]
 * points_ay = [3.0, 0.0, 3.0, 3.0]
 * points_bx = [0.0, 3.0, 2.0, 0.0]
 * points_by = [3.0, 3.0, 6.0, 3.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test74) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 1.0, 3.0);
        Vec2 a2( 2.0, 0.0);
        Vec2 a3( 4.0, 3.0);

        Vec2 b1( 0.0, 3.0);
        Vec2 b2( 3.0, 3.0);
        Vec2 b3( 2.0, 6.0);

        Vec2 c1( 1.0, 3.0);
        Vec2 c2( 3.0, 3.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;
        c2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1, -1,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,   1,  3,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);

//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c1,  11, 12,  1, -1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c2,  12, -1,  1,  3,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        expected_02.push_back(oe_03);
        expected_02.push_back(oe_04);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);

//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 2.0, 3.0, 0.0]
 * points_ay = [3.0, 0.0, 3.0, 3.0]
 * points_bx = [0.0, 3.0, 2.0, 0.0]
 * points_by = [3.0, 3.0, 6.0, 3.0]
 *
 * points_ax = [7.5737, 9.66694, 10.5722, 7.5737]
 * points_ay = [4.13994, 1.20424, 4.23417, 4.13994]
 * points_bx = [7.5737, 10.5722, 9.47848, 7.5737]
 * points_by = [4.13994, 4.23417, 7.20128, 4.13994]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test75) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);


        Vec2 a1( 0.0, 3.0);
        Vec2 a2( 2.0, 0.0);
        Vec2 a3( 3.0, 3.0);

        Vec2 b1( 0.0, 3.0);
        Vec2 b2( 3.0, 3.0);
        Vec2 b3( 2.0, 6.0);

        Vec2 c1( 0.0, 3.0);
        Vec2 c2( 3.0, 3.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;
        c2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1, -1,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,   3, -1,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);
//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c1,  11, -1,  1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c2,  12, -1,  3, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_02.push_back(oe_03);
        expected_02.push_back(oe_04);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);
//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [1.0, 2.0, 3.0, 1.0]
 * points_ay = [3.0, 0.0, 3.0, 3.0]
 * points_bx = [0.0, 3.0, 2.0, 0.0]
 * points_by = [3.0, 3.0, 6.0, 3.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test76) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 1.0, 3.0);
        Vec2 a2( 2.0, 0.0);
        Vec2 a3( 3.0, 3.0);

        Vec2 b1( 0.0, 3.0);
        Vec2 b2( 3.0, 3.0);
        Vec2 b3( 2.0, 6.0);

        Vec2 c1( 1.0, 3.0);
        Vec2 c2( 3.0, 3.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;
        c2 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1, -1,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,   3, -1,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);
//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c1,  11, 12,  1, -1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c2,  12, -1,  3, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_02.push_back(oe_03);
        expected_02.push_back(oe_04);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);
//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 3.0, 2.0, 0.0]
 * points_ay = [0.0, 0.0, 3.0, 0.0]
 * points_bx = [0.0, 3.0, 1.0, 0.0]
 * points_by = [3.0, 3.0, 6.0, 3.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test77) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 3.0, 0.0);
        Vec2 a3( 2.0, 3.0);

        Vec2 b1( 0.0, 3.0);
        Vec2 b2( 3.0, 3.0);
        Vec2 b3( 1.0, 6.0);

        Vec2 c1( 2.0, 3.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   3, -1,  11, 12, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        expected_01.push_back(oe_01);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);
//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c1,  11, 12,   3, -1,
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        expected_02.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);
//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 3.0, 3.0, 0.0]
 * points_ay = [0.0, 0.0, 3.0, 0.0]
 * points_bx = [0.0, 3.0, 1.0, 0.0]
 * points_by = [3.0, 3.0, 6.0, 3.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test78) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 3.0, 0.0);
        Vec2 a3( 3.0, 3.0);

        Vec2 b1( 0.0, 3.0);
        Vec2 b2( 3.0, 3.0);
        Vec2 b3( 1.0, 6.0);

        Vec2 c1( 3.0, 3.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   3, -1,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_01.push_back(oe_01);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);
//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c1,  12, -1,  3, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_02.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);
//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 3.0, 0.0, 0.0]
 * points_ay = [0.0, 0.0, 3.0, 0.0]
 * points_bx = [0.0, 3.0, 1.0, 0.0]
 * points_by = [3.0, 3.0, 6.0, 3.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test79) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 3.0, 0.0);
        Vec2 a3( 0.0, 3.0);

        Vec2 b1( 0.0, 3.0);
        Vec2 b2( 3.0, 3.0);
        Vec2 b3( 1.0, 6.0);

        Vec2 c1( 0.0, 3.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   3, -1,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_01.push_back(oe_01);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);
//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c1,  11, -1,  3, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_02.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);
//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 5.0, 2.0, 0.0]
 * points_ay = [0.0, 0.0, 3.0, 0.0]
 * points_bx = [2.0, 5.0, 1.0, 2.0]
 * points_by = [3.0, 5.0, 5.0, 3.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test80) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 5.0, 0.0);
        Vec2 a3( 2.0, 3.0);

        Vec2 b1( 2.0, 3.0);
        Vec2 b2( 5.0, 5.0);
        Vec2 b3( 1.0, 5.0);

        Vec2 c1( 2.0, 3.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   3, -1,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_01.push_back(oe_01);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);
//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c1, 11, -1,  3, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_02.push_back(oe_02);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);
//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);
    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [1.0, 3.0, 2.0, 1.0]
 * points_ay = [1.0, 1.0, 3.0, 1.0]
 * points_bx = [0.0, 4.0, 2.0, 0.0]
 * points_by = [0.0, 0.0, 4.0, 0.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test81) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 1.0, 1.0);
        Vec2 a2( 3.0, 1.0);
        Vec2 a3( 2.0, 3.0);

        Vec2 b1( 0.0, 0.0);
        Vec2 b2( 4.0, 0.0);
        Vec2 b3( 2.0, 4.0);

        Vec2 c1( 1.0, 1.0);
        Vec2 c2( 3.0, 1.0);
        Vec2 c3( 2.0, 3.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;
        c2 += trans;
        c3 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);
        c3  = rotatedVec(c3.x(), c3.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1, -1,  -1, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,   2, -1,  -1, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c3,   3, -1,  -1, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);
        expected_01.push_back(oe_03);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);
//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c1,   -1, -1, 1, -1,
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_05(
               c2,   -1, -1, 2, -1,
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_06(
               c3,   -1, -1, 3, -1,
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        expected_02.push_back(oe_04);
        expected_02.push_back(oe_05);
        expected_02.push_back(oe_06);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);
//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);
    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 4.0, 2.0, 0.0]
 * points_ay = [0.0, 0.0, 4.0, 0.0]
 * points_bx = [1.0, 3.0, 2.0, 1.0]
 * points_by = [1.0, 1.0, 3.0, 1.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test82) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 4.0, 0.0);
        Vec2 a3( 2.0, 4.0);

        Vec2 b1( 1.0, 1.0);
        Vec2 b2( 3.0, 1.0);
        Vec2 b3( 2.0, 3.0);

        Vec2 c1( 1.0, 1.0);
        Vec2 c2( 3.0, 1.0);
        Vec2 c3( 2.0, 3.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;
        c2 += trans;
        c3 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);
        c3  = rotatedVec(c3.x(), c3.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   -1, -1,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,   -1, -1,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c3,   -1, -1,  13, -1, 
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);
        expected_01.push_back(oe_03);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);
//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c1,  11, -1,  -1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        IntersectionFinderConvexPolygon2D::OutputElem oe_05(
               c2,  12, -1,  -1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        IntersectionFinderConvexPolygon2D::OutputElem oe_06(
               c3,  13, -1,  -1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        expected_02.push_back(oe_04);
        expected_02.push_back(oe_05);
        expected_02.push_back(oe_06);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);
//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);
    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 4.0, 2.0, 0.0]
 * points_ay = [0.0, 0.0, 4.0, 0.0]
 * points_bx = [0.0, 4.0, 2.0, 0.0]
 * points_by = [0.0, 0.0, 2.0, 0.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test83) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 4.0, 0.0);
        Vec2 a3( 2.0, 4.0);

        Vec2 b1( 0.0, 0.0);
        Vec2 b2( 4.0, 0.0);
        Vec2 b3( 2.0, 2.0);

        Vec2 c1( 0.0, 0.0);
        Vec2 c2( 4.0, 0.0);
        Vec2 c3( 2.0, 2.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;
        c2 += trans;
        c3 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);
        c3  = rotatedVec(c3.x(), c3.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1, -1,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,   2, -1,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c3,   -1, -1,  13, -1, 
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);
        expected_01.push_back(oe_03);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);
//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c1,  11, -1,  1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_05(
               c2,  12, -1,  2, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_06(
               c3,  13, -1,  -1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        expected_02.push_back(oe_04);
        expected_02.push_back(oe_05);
        expected_02.push_back(oe_06);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);
//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}



/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 4.0, 2.0, 0.0]
 * points_ay = [0.0, 0.0, 4.0, 0.0]
 * points_bx = [1.0, 3.0, 2.0, 1.0]
 * points_by = [0.0, 0.0, 2.0, 0.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test84) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 4.0, 0.0);
        Vec2 a3( 2.0, 4.0);

        Vec2 b1( 1.0, 0.0);
        Vec2 b2( 3.0, 0.0);
        Vec2 b3( 2.0, 2.0);

        Vec2 c1( 1.0, 0.0);
        Vec2 c2( 3.0, 0.0);
        Vec2 c3( 2.0, 2.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;
        c2 += trans;
        c3 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);
        c3  = rotatedVec(c3.x(), c3.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,   1, 2,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,   1, 2,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c3,   -1, -1,  13, -1, 
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);
        expected_01.push_back(oe_03);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);
//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c1,  11, -1,  1, 2,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_05(
               c2,  12, -1,  1, 2,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE);

        IntersectionFinderConvexPolygon2D::OutputElem oe_06(
               c3,  13, -1,  -1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        expected_02.push_back(oe_04);
        expected_02.push_back(oe_05);
        expected_02.push_back(oe_06);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);
//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);

    }
}


/**  @brief IntersectionFinderConvexPolygon2D::findIntersection()
 * import numpy as np
 * import matplotlib.pyplot as plt
 *
 * points_ax = [0.0, 4.0, 2.0, 0.0]
 * points_ay = [0.0, 0.0, 4.0, 0.0]
 * points_bx = [1.0, 3.0, 2.0, 1.0]
 * points_by = [1.0, 1.0, 4.0, 1.0]
 *
 * fig, ax = plt.subplots()
 * ax.plot(points_ax, points_ay, color='r')
 * ax.plot(points_bx, points_by, color='b')
 * plt.show()
 */
TEST_F(IntersectionFinderConvexPolygon2DTests, Test85) {

    for (double i = 0.0; i < 2*pi; i += (pi/100.0)) {

        Vec2 trans((rand()%100)/10.0, (rand()%100)/10.0);

        Vec2 a1( 0.0, 0.0);
        Vec2 a2( 4.0, 0.0);
        Vec2 a3( 2.0, 4.0);

        Vec2 b1( 1.0, 1.0);
        Vec2 b2( 3.0, 1.0);
        Vec2 b3( 2.0, 4.0);

        Vec2 c1( 1.0, 1.0);
        Vec2 c2( 3.0, 1.0);
        Vec2 c3( 2.0, 4.0);

        a1 += trans;
        a2 += trans;
        a3 += trans;
        b1 += trans;
        b2 += trans;
        b3 += trans;

        c1 += trans;
        c2 += trans;
        c3 += trans;

        a1  = rotatedVec(a1.x(), a1.y(), i);
        a2  = rotatedVec(a2.x(), a2.y(), i);
        a3  = rotatedVec(a3.x(), a3.y(), i);
        b1  = rotatedVec(b1.x(), b1.y(), i);
        b2  = rotatedVec(b2.x(), b2.y(), i);
        b3  = rotatedVec(b3.x(), b3.y(), i);

        c1  = rotatedVec(c1.x(), c1.y(), i);
        c2  = rotatedVec(c2.x(), c2.y(), i);
        c3  = rotatedVec(c3.x(), c3.y(), i);

        IntersectionFinderConvexPolygon2D::InputElem ie_a1(a1, 1);
        IntersectionFinderConvexPolygon2D::InputElem ie_a2(a2, 2);
        IntersectionFinderConvexPolygon2D::InputElem ie_a3(a3, 3);
        IntersectionFinderConvexPolygon2D::InputElem ie_b1(b1, 11);
        IntersectionFinderConvexPolygon2D::InputElem ie_b2(b2, 12);
        IntersectionFinderConvexPolygon2D::InputElem ie_b3(b3, 13);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_a01;
        input_a01.push_back(ie_a1);
        input_a01.push_back(ie_a2);
        input_a01.push_back(ie_a3);

        vector<IntersectionFinderConvexPolygon2D::InputElem> input_b01;
        input_b01.push_back(ie_b1);
        input_b01.push_back(ie_b2);
        input_b01.push_back(ie_b3);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_01;

        IntersectionFinderConvexPolygon2D::OutputElem oe_01(
               c1,  -1, -1,  11, -1, 
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_02(
               c2,  -1, -1,  12, -1, 
               IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX);

        IntersectionFinderConvexPolygon2D::OutputElem oe_03(
               c3,   3, -1,  13, -1, 
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_01.push_back(oe_01);
        expected_01.push_back(oe_02);
        expected_01.push_back(oe_03);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_01;

        auto res = mFinder.
                          findIntersection(input_a01, input_b01, output_01);
        EXPECT_EQ(res, true);
//        for (auto& v : output_01) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_01, expected_01), true);


        vector<IntersectionFinderConvexPolygon2D::OutputElem> expected_02;

        IntersectionFinderConvexPolygon2D::OutputElem oe_04(
               c1, 11, -1,  -1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        IntersectionFinderConvexPolygon2D::OutputElem oe_05(
               c2, 12, -1,  -1, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR);

        IntersectionFinderConvexPolygon2D::OutputElem oe_06(
               c3, 13, -1,  3, -1,
               IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX);

        expected_02.push_back(oe_04);
        expected_02.push_back(oe_05);
        expected_02.push_back(oe_06);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> output_02;

        res = mFinder.
                          findIntersection(input_b01, input_a01, output_02);
        EXPECT_EQ(res, true);
//        for (auto& v : output_02) {
//            cerr << "v:(" << v.mP.x() << "," << v.mP.y() << ")\n";
//            cerr << "mIndexA: " << v.mIndexA << "\n";
//            cerr << "mIndexAaux: " << v.mIndexAaux << "\n";
//            cerr << "mIndexB: " << v.mIndexB << "\n";
//            cerr << "mIndexBaux: " << v.mIndexBaux << "\n";
//            cerr << "mType: " << v.mType << "\n";
//        }

        EXPECT_EQ(compareResultsIntsec(output_02, expected_02), true);
    }
}


} // namespace Makena

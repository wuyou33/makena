#include "gtest/gtest.h"
#include <vector>
#include "primitives.hpp"
#include "manifold.hpp"
#include "quaternion.hpp"
#include "jacobian_constraint.hpp"
#include "contact_manager.hpp"

#include <iostream>
#include <fstream>
#include <sstream>


namespace Makena {


class ContactManagerTest : public ::testing::Test {

  protected:  

    ContactManagerTest(){;}
    virtual ~ContactManagerTest(){;}
    virtual void SetUp() {;};
    virtual void TearDown() {;};

};


static VertexIt findVertex(
    const Vec3&      featurePoint,
    ConvexRigidBody& body
) {
    auto& ch = body.ConvexHull();
    for (auto vit  = ch.vertices().first; 
              vit != ch.vertices().second; vit++) {
        if (featurePoint == (*vit)->pLCS()) {
            return vit;
        }
    }
    return ch.vertices().second;
}


static FaceIt findFace(
    vector<Vec3>&    featurePoints,
    ConvexRigidBody& body
) {
    auto& ch = body.ConvexHull();
    for (auto fit = ch.faces().first; fit != ch.faces().second; fit++) {
        vector<long> foundCnts;
        for (long i = 0; i < featurePoints.size(); i++) {
            foundCnts.push_back(0);
        }
        for (auto& he : (*fit)->halfEdges()) {
            auto vit = (*he)->src();
            auto p = (*vit)->pLCS();
            for (long i = 0; i < foundCnts.size(); i++) {
                if (featurePoints[i]== p) {
                    foundCnts[i] = foundCnts[i] + 1;
                }
            }
        }
        long foundCnt = 0;
        for (long i = 0; i < featurePoints.size(); i++) {
            if (foundCnts[i] == 1) {
                foundCnt++;
            }
        }
        if (foundCnt==featurePoints.size()) {
            return fit;
        }
    }
    return ch.faces().second;
}


static EdgeIt findEdge(
    vector<Vec3>&    featurePoints,
    ConvexRigidBody& body
) {
    auto& ch = body.ConvexHull();
    for (auto eit = ch.edges().first; eit != ch.edges().second; eit++) {
        auto heit = (*eit)->he1();
        auto vit1 = (*heit)->src();
        auto vit2 = (*heit)->dst();
        auto p1   = (*vit1)->pLCS();
        auto p2   = (*vit2)->pLCS();
        if ( (featurePoints[0] == p1 && featurePoints[1] == p2) ||
             (featurePoints[1] == p1 && featurePoints[0] == p2)   ) {
            return eit;
        }
    }
    return ch.edges().second;
}


static VertexIt findVertex(
    vector<Vec3>&    featurePoints,
    ConvexRigidBody& body
) {
    auto& ch = body.ConvexHull();
    for (auto vit  = ch.vertices().first; 
              vit != ch.vertices().second; vit++) {
        if (featurePoints[0] == (*vit)->pLCS()) {
            return vit;
        }
    }
    return ch.vertices().second;
}



static double rand100()
{
    return ((rand()%10000)+1)/100.0;
}


static Vec3 randVec3D100()
{
    return Vec3(rand100(), rand100(), rand100());
}


static Quaternion randomRotQ()
{
    Vec3 v(0.0, 0.0, 0.0);
    while (v == Vec3(0.0, 0.0, 0.0)) {
        v = Vec3( (rand()%100)/100.0 - 0.5, 
                  (rand()%100)/100.0 - 0.5, 
                  (rand()%100)/100.0 - 0.5 );
    }
    
    Quaternion q(v, (rand()%100)/50.0 * M_PI);
    q.normalize();
    return q;
}


static void generateTestRigidBodies(
    vector<Vec3>&    polygon_1,
    const Vec3&      anchor_1,
    vector<Vec3>&    polygon_2,
    const Vec3&      anchor_2,
    Vec3             tiltAxis,
    const double&    rotationRad,
    ConvexRigidBody& body1,
    ConvexRigidBody& body2,
    const double&    zTwist,
    bool             swapBody12
) {

    enum predicate pred;

    Manifold m1;
    m1.findConvexHull(polygon_1, pred);
    auto m1Serial = m1.exportData();

    Vec3 CoM1 = anchor_1 * -1.0;

    Manifold m2;
    m2.findConvexHull(polygon_2, pred);
    auto m2Serial = m2.exportData();

    Vec3 CoM2 = anchor_2 * -1.0;

    Vec3 zAxis(0.0, 0.0, 1.0);
    Quaternion twistQ (zAxis, zTwist);

    tiltAxis.normalize();

    Quaternion tiltQ (tiltAxis, rotationRad);
    tiltQ = tiltQ * twistQ;

    CoM2 = tiltQ.rotate(CoM2);
    Mat3x3 I ( 1.0, 0.0, 0.0, 
               0.0, 1.0, 0.0, 
               0.0, 0.0, 1.0  );

    Quaternion Qunit(1.0, 0.0, 0.0, 0.0);

    tiltAxis.scale(1.0/10.0);

    if (!swapBody12) {

        body1.setObjectAttributes(m1Serial, 1.0, I, 1.0, 1.0, false);
        body1.setGeomConfig(
                        CoM1, Qunit, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0));
        body1.updateGeomConfig(0.0);
        body1.commitGeomConfig();
        body2.setObjectAttributes(m2Serial, 1.0, I, 1.0, 1.0, false);
        body2.setGeomConfig(
                        CoM2, tiltQ, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0));
        body2.updateGeomConfig(0.0);
        body2.commitGeomConfig();
    }
    else {
        body2.setObjectAttributes(m1Serial, 1.0, I, 1.0, 1.0, false);
        body2.setGeomConfig(
                        CoM1, Qunit, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0));
        body2.updateGeomConfig(0.0);
        body2.commitGeomConfig();

        body1.setObjectAttributes(m2Serial, 1.0, I, 1.0, 1.0, false);
        body1.setGeomConfig(
                        CoM2, tiltQ, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0));
        body1.updateGeomConfig(0.0);
        body1.commitGeomConfig();
    }
}


static void setActiveFeatures(
    vector<Vec3>&    featurePoints1tmp,
    vector<Vec3>&    featurePoints2tmp,
    ConvexRigidBody& body1,
    ConvexRigidBody& body2,
    ContactPairInfo& info,
    bool             swapBody12
) {
    vector<Vec3>& featurePoints1 = 
                             swapBody12?featurePoints2tmp:featurePoints1tmp;
    vector<Vec3>& featurePoints2 = 
                             swapBody12?featurePoints1tmp:featurePoints2tmp;

    info.mContactGenerationIrregular = false;

    if (featurePoints1.size() >= 3) {
        info.mType1 = ContactPairInfo::FT_FACE;
        info.mFit1  = findFace(featurePoints1, body1);

        if (featurePoints2.size() >= 3) {
            info.mType2 = ContactPairInfo::FT_FACE;
            info.mFit2  = findFace(featurePoints2, body2);
            info.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_1;
        }
        else if (featurePoints2.size() == 2) {
            info.mType2 = ContactPairInfo::FT_EDGE;
            info.mEit2  = findEdge(featurePoints2, body2);
            info.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_1;
        }
        else if (featurePoints2.size() == 1) {
            info.mType2 = ContactPairInfo::FT_VERTEX;
            info.mVit2  = findVertex(featurePoints2, body2);
            info.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_1;
        }
    }
    else if (featurePoints1.size() == 2) {
        info.mType1 = ContactPairInfo::FT_EDGE;
        info.mEit1  = findEdge(featurePoints1, body1);

        if (featurePoints2.size() >= 3) {
            info.mType2 = ContactPairInfo::FT_FACE;
            info.mFit2  = findFace(featurePoints2, body2);
            info.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;
        }
        else if (featurePoints2.size() == 2) {
            info.mType2 = ContactPairInfo::FT_EDGE;
            info.mEit2  = findEdge(featurePoints2, body2);
            info.mContactNormalType = ContactPairInfo::NT_EDGE_CROSS_EDGE;
                                        
        }
        else if (featurePoints2.size() == 1) {
            info.mType2 = ContactPairInfo::FT_VERTEX;
            info.mVit2  = findVertex(featurePoints2, body2);
            info.mContactNormalType = ContactPairInfo::NT_EDGE_NORMAL_1; 
        }
    }
    else if (featurePoints1.size() == 1) {
        info.mType1 = ContactPairInfo::FT_VERTEX;
        info.mVit1  = findVertex(featurePoints1, body1);

        if (featurePoints2.size() >= 3) {
            info.mType2 = ContactPairInfo::FT_FACE;
            info.mFit2  = findFace(featurePoints2, body2);
            info.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;
        }
        else if (featurePoints2.size() == 2) {
            info.mType2 = ContactPairInfo::FT_EDGE;
            info.mEit2  = findEdge(featurePoints2, body2);
            info.mContactNormalType = ContactPairInfo::NT_EDGE_NORMAL_2; 
        }
        else if (featurePoints2.size() == 1) {
            info.mType2 = ContactPairInfo::FT_VERTEX;
            info.mVit2  = findVertex(featurePoints2, body2);
            info.mContactNormalType = ContactPairInfo::NT_VERTEX_VERTEX_AVG;
        }
    }
}


// Move the object further to simulate the next tims step.
static void updateGeomConfig(
    const Vec3&       CoM1,
    const Quaternion& q1,
    const Vec3&       Vlin1,
    const Vec3&       Vang1,
    ConvexRigidBody&  body1
) {
    const Vec3       comNew1  = CoM1  + body1.CoM();
    Quaternion       qNewTmp  = q1 * body1.Q();
    qNewTmp.normalize();
    const Quaternion qNew1 = qNewTmp;
    const Vec3       vlinNew1 = Vlin1 + body1.Vlin();
    const Vec3       vangNew1 = Vang1 + body1.Vang();

    body1.setGeomConfig(comNew1, qNew1, vlinNew1, vangNew1);
}


static void rotateTestPatternRandomly(
    ConvexRigidBody&  body1,
    ConvexRigidBody&  body2,
    ContactPairInfo&  info,
    const Quaternion& q,
    const Vec3&       v
) {
    auto Q1cur = body1.Q();
    auto C1cur = body1.CoM();
    auto V1lincur = body1.Vlin();
    auto V1angcur = body1.Vang();
    auto Q2cur = body2.Q();
    auto C2cur = body2.CoM();
    auto V2lincur = body2.Vlin();
    auto V2angcur = body2.Vang();
    auto Q1new = q * Q1cur;
    auto Q2new = q * Q2cur;
    auto C1new = q.rotate(C1cur + v);
    auto C2new = q.rotate(C2cur + v);
    auto V1linnew = q.rotate(V1lincur);
    auto V2linnew = q.rotate(V2lincur);
    auto V1angnew = q.rotate(V1angcur);
    auto V2angnew = q.rotate(V2angcur);

    body1.setGeomConfig(C1new, Q1new, V1linnew, V1angnew);
    body1.updateGeomConfig(0.0);
    body1.commitGeomConfig();

    body2.setGeomConfig(C2new, Q2new, V2linnew, V2angnew);
    body2.updateGeomConfig(0.0);
    body2.commitGeomConfig();
}




/**  @brief getPerpendicularVectors
 */
TEST_F(ContactManagerTest, Test01) {

    for (long i = 0; i < 1000; i++) {
        Vec3 v_01 = randVec3D100();
        if (v_01.squaredNorm2() <= EPSILON_SQUARED) {
            continue;
        }
        v_01.normalize();

        ContactManager mgr01 (
                   0.0,   // KCorr
                   0.0,   // velocityThreshold
                   0.0,   // maxNumPointsPerContact
                   0.0,   // epsilonZero
                   0.0,   // epsilonAngle
                   0.0,   // epsilonZeroGJK
                   0.0,   // scalingGJK
                   0,     // maxNumIter
                   0,     // maxNumCycles
                   std::cerr );

        Vec3 v_02, v_03;
        mgr01.getPerpendicularVectors(v_01, v_02, v_03);

        EXPECT_EQ(fabs(v_01.dot(v_02)) <= EPSILON_SQUARED, true);
        EXPECT_EQ(fabs(v_01.dot(v_03)) <= EPSILON_SQUARED, true);
        EXPECT_EQ(fabs(v_02.dot(v_03)) <= EPSILON_SQUARED, true);
        EXPECT_EQ(fabs(v_02.squaredNorm2() - 1.0) <= EPSILON_SQUARED, true);
        EXPECT_EQ(fabs(v_03.squaredNorm2() - 1.0) <= EPSILON_SQUARED, true);
        Vec3 cr_12 = v_01.cross(v_02);
        EXPECT_EQ(cr_12, v_03);
        Vec3 cr_23 = v_02.cross(v_03);
        EXPECT_EQ(cr_23, v_01);
        Vec3 cr_31 = v_03.cross(v_01);
        EXPECT_EQ(cr_31, v_02);
    }

}


/**  @break findRelativeVelocityOfBody1RelativeToBody2()
 */
TEST_F(ContactManagerTest, Test02) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 0.0,  0.0,  1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 1.0);
        Vec3       Vang1 (-0.1, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, -1.0);
        Vec3       Vang2 (-0.1, 0.0, 0.0);

        Vec3          tiltAxis(-1.0, 0.0, 0.0);
        double       rotationRad = 0.0;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            0.0,
            false
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            false
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, body2);

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);

        auto vit1_01 = findVertex(Vec3( 0.0, 0.0,  1.0),body1);
        auto vit2_01 = findVertex(Vec3( 0.0, 0.0, -1.0),body2);

        BDVertex bd11(vit1_01, vit2_01);
        info.mLastVertices.clear();
        info.mLastVertices.push_back(bd11);

        Vec3 exp_01(0.0, 0.2, 2.0);
        exp_01 = rm.rotate(exp_01);


        ContactManager mgr01 (
                   0.0,   // KCorr
                   0.0,   // velocityThreshold
                   0.0,   // maxNumPointsPerContact
                   0.0,   // epsilonZero
                   0.0,   // epsilonAngle
                   0.0,   // epsilonZeroGJK
                   0.0,   // scalingGJK
                   0,     // maxNumIter
                   0,     // maxNumCycles
                   std::cerr );

        Vec3 v01 = mgr01.findRelativeVelocityOfBody1RelativeToBody2(
                                                   body1, body2, info, true);

        //EXPECT_EQ(v01, exp_01);
    }
}


/**  @brief decomposeQ()
 */
TEST_F(ContactManagerTest, Test03) {

    vector<Vec3> polygon_1;
    polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
    polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));

    Vec3         anchor_1(0.0,  0.0,  1.0);

    vector<Vec3> polygon_2;
    polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_2.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_2.push_back(Vec3(-1.0,  1.0, -1.0));

    Vec3         anchor_2(0.0,  0.0, -1.0);

    ConvexRigidBody body1(1);
    ConvexRigidBody body2(2);
    ContactPairInfo info;

    generateTestRigidBodies(
        polygon_1,
        anchor_1,
        polygon_2,
        anchor_2,
        Vec3(1.0, 0.0, 0.0),
        0.0,
        body1,
        body2,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3( 1.0, 1.0,  1.0),body1);
    auto vit1_02 = findVertex(Vec3(-1.0, 1.0,  1.0),body1);
    auto vit1_03 = findVertex(Vec3( 1.0,-1.0,  1.0),body1);
    auto vit1_04 = findVertex(Vec3(-1.0,-1.0,  1.0),body1);
//    auto vit1_05 = findVertex(Vec3( 1.0, 1.0, -1.0),body1);
//    auto vit1_06 = findVertex(Vec3(-1.0, 1.0, -1.0),body1);
//    auto vit1_07 = findVertex(Vec3( 1.0,-1.0, -1.0),body1);
//    auto vit1_08 = findVertex(Vec3(-1.0,-1.0, -1.0),body1);

    auto vit2_01 = findVertex(Vec3( 1.0, 1.0,  1.0),body2);
    auto vit2_02 = findVertex(Vec3(-1.0, 1.0,  1.0),body2);
    auto vit2_03 = findVertex(Vec3( 1.0,-1.0,  1.0),body2);
    auto vit2_04 = findVertex(Vec3(-1.0,-1.0,  1.0),body2);
//    auto vit2_05 = findVertex(Vec3( 1.0, 1.0, -1.0),body2);
//    auto vit2_06 = findVertex(Vec3(-1.0, 1.0, -1.0),body2);
//    auto vit2_07 = findVertex(Vec3( 1.0,-1.0, -1.0),body2);
//    auto vit2_08 = findVertex(Vec3(-1.0,-1.0, -1.0),body2);

    BDVertex bd11(vit1_01, vit2_01);
    BDVertex bd12(vit1_01, vit2_02);
    BDVertex bd13(vit1_01, vit2_03);
    BDVertex bd14(vit1_01, vit2_04);
    BDVertex bd21(vit1_02, vit2_01);
    BDVertex bd22(vit1_02, vit2_02);
    BDVertex bd23(vit1_02, vit2_03);
    BDVertex bd24(vit1_02, vit2_04);
    BDVertex bd31(vit1_03, vit2_01);
    BDVertex bd32(vit1_03, vit2_02);
    BDVertex bd33(vit1_03, vit2_03);
    BDVertex bd34(vit1_03, vit2_04);
    BDVertex bd41(vit1_04, vit2_01);
    BDVertex bd42(vit1_04, vit2_02);
    BDVertex bd43(vit1_04, vit2_03);
    BDVertex bd44(vit1_04, vit2_04);

    vector<BDVertex> Q;
    vector<VertexIt> vertices1;
    vector<VertexIt> vertices2;


    ContactManager mgr01 (
                   0.0,   // KCorr
                   0.0,   // velocityThreshold
                   0.0,   // maxNumPointsPerContact
                   0.0,   // epsilonZero
                   0.0,   // epsilonAngle
                   0.0,   // epsilonZeroGJK
                   0.0,   // scalingGJK
                   0,     // maxNumIter
                   0,     // maxNumCycles
                   std::cerr );

    Q.push_back(bd11);
    mgr01.decomposeQ(Q, vertices1, vertices2);
    EXPECT_EQ(vertices1.size(), 1);
    EXPECT_EQ(vertices2.size(), 1);
    EXPECT_EQ(vertices1[0], vit1_01);
    EXPECT_EQ(vertices2[0], vit2_01);

    Q.clear();
    vertices1.clear();
    vertices2.clear();
    Q.push_back(bd11);
    Q.push_back(bd12);
    mgr01.decomposeQ(Q, vertices1, vertices2);
    EXPECT_EQ(vertices1.size(), 1);
    EXPECT_EQ(vertices2.size(), 2);
    EXPECT_EQ(vertices1[0], vit1_01);
    EXPECT_EQ(vertices2[0], vit2_01);
    EXPECT_EQ(vertices2[1], vit2_02);

    Q.clear();
    vertices1.clear();
    vertices2.clear();
    Q.push_back(bd11);
    Q.push_back(bd12);
    Q.push_back(bd13);
    mgr01.decomposeQ(Q, vertices1, vertices2);
    EXPECT_EQ(vertices1.size(), 1);
    EXPECT_EQ(vertices2.size(), 3);
    EXPECT_EQ(vertices1[0], vit1_01);
    EXPECT_EQ(vertices2[0], vit2_01);
    EXPECT_EQ(vertices2[1], vit2_02);
    EXPECT_EQ(vertices2[2], vit2_03);

    Q.clear();
    vertices1.clear();
    vertices2.clear();
    Q.push_back(bd11);
    Q.push_back(bd12);
    Q.push_back(bd13);
    Q.push_back(bd14);
    mgr01.decomposeQ(Q, vertices1, vertices2);
    EXPECT_EQ(vertices1.size(), 1);
    EXPECT_EQ(vertices2.size(), 4);
    EXPECT_EQ(vertices1[0], vit1_01);
    EXPECT_EQ(vertices2[0], vit2_01);
    EXPECT_EQ(vertices2[1], vit2_02);
    EXPECT_EQ(vertices2[2], vit2_03);
    EXPECT_EQ(vertices2[3], vit2_04);

    Q.clear();
    vertices1.clear();
    vertices2.clear();
    Q.push_back(bd11);
    Q.push_back(bd12);
    Q.push_back(bd22);
    mgr01.decomposeQ(Q, vertices1, vertices2);
    EXPECT_EQ(vertices1.size(), 2);
    EXPECT_EQ(vertices2.size(), 2);
    EXPECT_EQ(vertices1[0], vit1_01);
    EXPECT_EQ(vertices1[1], vit1_02);
    EXPECT_EQ(vertices2[0], vit2_01);
    EXPECT_EQ(vertices2[1], vit2_02);

    Q.clear();
    vertices1.clear();
    vertices2.clear();
    Q.push_back(bd11);
    Q.push_back(bd21);
    Q.push_back(bd31);
    Q.push_back(bd41);
    mgr01.decomposeQ(Q, vertices1, vertices2);
    EXPECT_EQ(vertices1.size(), 4);
    EXPECT_EQ(vertices2.size(), 1);
    EXPECT_EQ(vertices1[0], vit1_01);
    EXPECT_EQ(vertices1[1], vit1_02);
    EXPECT_EQ(vertices1[2], vit1_03);
    EXPECT_EQ(vertices1[3], vit1_04);
    EXPECT_EQ(vertices2[0], vit2_01);

    Q.clear();
    vertices1.clear();
    vertices2.clear();
    Q.push_back(bd11);
    Q.push_back(bd21);
    Q.push_back(bd32);
    Q.push_back(bd41);

    mgr01.decomposeQ(Q, vertices1, vertices2);
    EXPECT_EQ(vertices1.size(), 4);
    EXPECT_EQ(vertices2.size(), 2);
    EXPECT_EQ(vertices1[0], vit1_01);
    EXPECT_EQ(vertices1[1], vit1_02);
    EXPECT_EQ(vertices2[0], vit2_01);
    EXPECT_EQ(vertices2[1], vit2_02);
    EXPECT_EQ(vertices2[2], vit2_03);
    EXPECT_EQ(vertices2[3], vit2_04);

    Q.clear();
    vertices1.clear();
    vertices2.clear();
    Q.push_back(bd11);
    Q.push_back(bd22);
    Q.push_back(bd33);
    Q.push_back(bd44);

    mgr01.decomposeQ(Q, vertices1, vertices2);
    EXPECT_EQ(vertices1.size(), 4);
    EXPECT_EQ(vertices2.size(), 4);
    EXPECT_EQ(vertices1[0], vit1_01);
    EXPECT_EQ(vertices1[1], vit1_02);
    EXPECT_EQ(vertices1[2], vit1_03);
    EXPECT_EQ(vertices1[3], vit1_04);
    EXPECT_EQ(vertices2[0], vit2_01);
    EXPECT_EQ(vertices2[1], vit2_02);
    EXPECT_EQ(vertices2[2], vit2_03);
    EXPECT_EQ(vertices2[3], vit2_04);

}


//
}// Makena

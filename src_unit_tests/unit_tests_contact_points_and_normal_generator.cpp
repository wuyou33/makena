#include "gtest/gtest.h"
#include <vector>
#include "primitives.hpp"
#include "manifold.hpp"
#include "quaternion.hpp"
#include "contact_points_and_normal_generator.hpp"

#include <iostream>
#include <fstream>
#include <sstream>


namespace Makena {


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


static VertexIt findVertex(
    const Vec3&      p,
    ConvexRigidBody& body
) {
    auto& ch = body.ConvexHull();
    for (auto vit  = ch.vertices().first; 
              vit != ch.vertices().second; vit++) {
        if (p == (*vit)->pLCS()) {
            return vit;
        }
    }
    return ch.vertices().second;
}


class ContactPointsAndNormalGeneratorTest : public ::testing::Test {

  protected:  

    ContactPointsAndNormalGeneratorTest(){;};
    virtual ~ContactPointsAndNormalGeneratorTest(){;};
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

static double rand100()
{
    return ((rand()%10000)+1)/100.0;
}


static Vec3 randVec3D100()
{
    return Vec3(rand100(), rand100(), rand100());
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
        body2.setGeomConfig(CoM2, tiltQ, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0));
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
        body1.setGeomConfig(CoM2, tiltQ, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0));
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
    vector<Vec3>& featurePoints1 = swapBody12?featurePoints2tmp:featurePoints1tmp;
    vector<Vec3>& featurePoints2 = swapBody12?featurePoints1tmp:featurePoints2tmp;

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



// Move the object further to simulate the next tims step.
static void setVelocities(
    const Vec3&       Vlin1,
    const Vec3&       Vang1,
    ConvexRigidBody&  body1
) {
    const Vec3       comNew1  = body1.CoM();
    const Quaternion qNewTmp  = body1.Q();
    const Vec3       vlinNew1 = Vlin1;
    const Vec3       vangNew1 = Vang1;

    body1.setGeomConfig(comNew1, qNewTmp, vlinNew1, vangNew1);
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


class ExpStructContactPointsAndNormalGenerator {

  public:

    ExpStructContactPointsAndNormalGenerator(
        enum ContactPairInfo::FeatureType t1, 
        vector<Vec3>& pts1,
        vector<Vec3>& intsec1,
        enum ContactPairInfo::FeatureType t2,
        vector<Vec3>& pts2,
        vector<Vec3>& intsec2
    ):  mType1(t1),
        mPts1(pts1),
        mIntsec1(intsec1),
        mType2(t2),
        mPts2(pts2),
        mIntsec2(intsec2){;}
  
    enum ContactPairInfo::FeatureType mType1;
    vector<Vec3> mPts1;
    vector<Vec3> mIntsec1;
    enum ContactPairInfo::FeatureType mType2;
    vector<Vec3> mPts2;
    vector<Vec3> mIntsec2;
};


static bool areSamePoints(vector<Vec3>& vec1, vector<Vec3>& vec2)
{
    if (vec1.size() != vec2.size()) {
        return false;
    }
    vector<long> cnts;
    for (long i = 0; i < vec1.size(); i++) {
        cnts.push_back(0);
    }
    for (long i = 0; i < vec1.size(); i++) {
        for (long j = 0; j < vec1.size(); j++) {
            if (vec1[i] == vec2[j]) {
                cnts[i] = cnts[i] + 1;
            }
        }
    }
    for (auto& c : cnts) {
        if (c!=1) {
            return false;
        }
    }
    return true;
}


static void showResult(
    ConvexRigidBody&         body1,
    ConvexRigidBody&         body2,
    ContactPairInfo&         info
) {
    cerr << "Body 1\n";
    if (info.mType1 == ContactPairInfo::FT_FACE) {
        cerr << "Face:\n";
        for (auto he : (*(info.mFit1))->halfEdges()) {
            auto vit = (*he)->src();            
            cerr << (*vit)->pLCS() << "\n"; 
        }
    }
    else if (info.mType1 == ContactPairInfo::FT_EDGE) {
        cerr << "Edge:\n";
        auto he = (*(info.mEit1))->he1();
        auto vit1 = (*he)->src();
        auto vit2 = (*he)->dst();
        cerr << (*vit1)->pLCS() << "\n";
        cerr << (*vit2)->pLCS() << "\n";
    }
    else {
        cerr << "Vertex: ";
        cerr << (*(info.mVit1))->pLCS() << "\n";
    }

    cerr << "Body 2\n";
    if (info.mType2 == ContactPairInfo::FT_FACE) {
        cerr << "Face:\n";
        for (auto he : (*(info.mFit2))->halfEdges()) {
            auto vit = (*he)->src();            
            cerr << (*vit)->pLCS() << "\n"; 
        }
    }
    else if (info.mType2 == ContactPairInfo::FT_EDGE) {
        cerr << "Edge:\n";
        auto he = (*(info.mEit2))->he1();
        auto vit1 = (*he)->src();
        auto vit2 = (*he)->dst();
        cerr << (*vit1)->pLCS() << "\n";
        cerr << (*vit2)->pLCS() << "\n";
    }
    else {
        cerr << "Vertex: ";
        cerr << (*(info.mVit2))->pLCS() << "\n";
    }

    cerr << "Intersection1:\n";
    for (auto& p : info.mIntsect1) {
        cerr << p << "\n";
    }
    cerr << "Intersection2:\n";
    for (auto& p : info.mIntsect2) {
        cerr << p << "\n";
    }

    cerr << "Contact Normal Type: ";
    switch (info.mContactNormalType) {
      case ContactPairInfo::NT_NONE:
        cerr << "NT_NONE";
        break;
      case ContactPairInfo::NT_FACE_NORMAL_1:
        cerr << "NT_FACE_NORMAL_1";
        break;
      case ContactPairInfo::NT_FACE_NORMAL_2:
        cerr << "NT_FACE_NORMAL_2";
        break;
      case ContactPairInfo::NT_EDGE_NORMAL_1:
        cerr << "NT_EDGE_NORMAL_1";
        break;
      case ContactPairInfo::NT_EDGE_NORMAL_2:
        cerr << "NT_EDGE_NORMAL_2";
        break;
      case ContactPairInfo::NT_EDGE_CROSS_EDGE:
        cerr << "NT_EDGE_CROSS_EDGE";
        break;
      case ContactPairInfo::NT_EDGE_EDGE_AVG:
        cerr << "NT_EDGE_EDGE_AVG";
        break;
      case ContactPairInfo::NT_VERTEX_VERTEX_AVG:
        cerr << "NT_VERTEX_VERTEX_AVG";
        break;
      case ContactPairInfo::NT_RELATIVE_VELOCITY:
        cerr << "NT_RELATIVE_VELOCITY";
        break;
      default:
        cerr <<  "<ERROR>";
        break;
    }
    cerr << "\n";
    cerr << "Contact Normal: " << info.mContactNormal1To2 << "\n";
}


static bool checkResult(
    ConvexRigidBody&         body1,
    ConvexRigidBody&         body2,
    ContactPairInfo&         info,
    ExpStructContactPointsAndNormalGenerator& exp,
    bool                     swapBody12
) {
    if (!swapBody12) {
        if (info.mType1 != exp.mType1 ||info.mType2 != exp.mType2) {
            return false;        
        }
    }
    else {
        //cerr << "swapped\n";
        if (info.mType2 != exp.mType1 ||info.mType1 != exp.mType2) {
            return false;        
        }
    }

    if (exp.mType1 == ContactPairInfo::FT_FACE) {
        auto fit = findFace(exp.mPts1, swapBody12?body2:body1);
        if (fit!=(swapBody12?info.mFit2:info.mFit1)) {
            return false;
        }
    }
    else if (exp.mType1 == ContactPairInfo::FT_EDGE) {
        auto eit = findEdge(exp.mPts1, swapBody12?body2:body1);
        if (eit!=(swapBody12?info.mEit2:info.mEit1)) {
            return false;
        }
    }
    else {
        auto vit = findVertex(exp.mPts1, swapBody12?body2:body1);
        if (vit!=(swapBody12?info.mVit2:info.mVit1)) {
            return false;
        }
    }

    if (exp.mType2 == ContactPairInfo::FT_FACE) {
        auto fit = findFace(exp.mPts2, swapBody12?body1:body2);
        if (fit!=(swapBody12?info.mFit1:info.mFit2)) {
            return false;
        }
    }
    else if (exp.mType2 == ContactPairInfo::FT_EDGE) {
        auto eit = findEdge(exp.mPts2, swapBody12?body1:body2);
        if (eit!=(swapBody12?info.mEit1:info.mEit2)) {
            return false;
        }
    }
    else {
        auto vit = findVertex(exp.mPts2, swapBody12?body1:body2);
        if (vit!=(swapBody12?info.mVit1:info.mVit2)) {
            return false;
        }
    }

    bool intsec1_OK = areSamePoints(swapBody12?info.mIntsect2:
                                               info.mIntsect1, exp.mIntsec1);
    bool intsec2_OK = areSamePoints(swapBody12?info.mIntsect1:
                                               info.mIntsect2, exp.mIntsec2);
    //cerr << "intsecs: " << intsec1_OK << "," << intsec2_OK << "\n";
    return intsec1_OK && intsec2_OK;
}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test01) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3(1.0, 1.0, 1.0), body1_01);
    auto vit2_01 = findVertex(Vec3(1.0, 1.0, 1.0), body2_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_VERTEX);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_VERTEX);
    EXPECT_EQ(info_01.mVit1,  vit1_01);
    EXPECT_EQ(info_01.mVit2,  vit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a1-b2 edge
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test02) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints2_01.push_back(Vec3( 1.0, -1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3(1.0, 1.0, 1.0), body1_01);
    auto vit2_01 = findVertex(Vec3(1.0, 1.0, 1.0), body2_01);
    auto vit2_02 = findVertex(Vec3(1.0,-1.0, 1.0), body2_01);

    auto eit2_01 = findEdge(featurePoints2_01, body2_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_01, vit2_02));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_VERTEX);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_EDGE);
    EXPECT_EQ(info_01.mVit1,  vit1_01);
    EXPECT_EQ(info_01.mEit2,  eit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a1-b2 face
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test03) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints2_01.push_back(Vec3(-1.0, -1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3(1.0, 1.0, 1.0), body1_01);
    auto vit2_01 = findVertex(Vec3(1.0, 1.0, 1.0), body2_01);
    auto vit2_02 = findVertex(Vec3(-1.0,-1.0, 1.0), body2_01);
    auto fit2_01 = findFace(featurePoints2_01, body2_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_01, vit2_02));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_VERTEX);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_FACE);
    EXPECT_EQ(info_01.mVit1,  vit1_01);
    EXPECT_EQ(info_01.mFit2,  fit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a2-b1 edge
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test04) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints1_01.push_back(Vec3( 1.0, -1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3(1.0, 1.0, 1.0), body1_01);
    auto vit1_02 = findVertex(Vec3(1.0,-1.0, 1.0), body1_01);
    auto vit2_01 = findVertex(Vec3(1.0, 1.0, 1.0), body2_01);

    auto eit1_01 = findEdge(featurePoints1_01, body1_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_02, vit2_01));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_EDGE);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_VERTEX);
    EXPECT_EQ(info_01.mEit1,  eit1_01);
    EXPECT_EQ(info_01.mVit2,  vit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a2-b1 face
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test05) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints1_01.push_back(Vec3(-1.0, -1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body1_01);
    auto vit1_02 = findVertex(Vec3(-1.0,-1.0, 1.0), body1_01);
    auto vit2_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body2_01);

    auto fit1_01 = findFace(featurePoints1_01, body1_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_02, vit2_01));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_FACE);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_VERTEX);
    EXPECT_EQ(info_01.mFit1,  fit1_01);
    EXPECT_EQ(info_01.mVit2,  vit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a2-b2 edge-edge
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test06) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints1_01.push_back(Vec3( 1.0, -1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints2_01.push_back(Vec3( 1.0, -1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3(1.0, 1.0, 1.0), body1_01);
    auto vit1_02 = findVertex(Vec3(1.0,-1.0, 1.0), body1_01);
    auto vit2_01 = findVertex(Vec3(1.0, 1.0, 1.0), body2_01);
    auto vit2_02 = findVertex(Vec3(1.0,-1.0, 1.0), body2_01);

    auto eit1_01 = findEdge(featurePoints1_01, body1_01);
    auto eit2_01 = findEdge(featurePoints2_01, body2_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_02, vit2_02));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_EDGE);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_EDGE);
    EXPECT_EQ(info_01.mEit1,  eit1_01);
    EXPECT_EQ(info_01.mEit2,  eit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a2-b2 face-face
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test07) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints1_01.push_back(Vec3(-1.0, -1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints2_01.push_back(Vec3(-1.0, -1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3(1.0, 1.0, 1.0), body1_01);
    auto vit1_02 = findVertex(Vec3(-1.0,-1.0, 1.0), body1_01);
    auto vit2_01 = findVertex(Vec3(1.0, 1.0, 1.0), body2_01);
    auto vit2_02 = findVertex(Vec3(-1.0,-1.0, 1.0), body2_01);

    auto fit1_01 = findFace(featurePoints1_01, body1_01);
    auto fit2_01 = findFace(featurePoints2_01, body2_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_02, vit2_02));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_FACE);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_FACE);
    EXPECT_EQ(info_01.mFit1,  fit1_01);
    EXPECT_EQ(info_01.mFit2,  fit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a1-b2, a1-b3
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test08) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints2_01.push_back(Vec3(-1.0, -1.0,  1.0));
    featurePoints2_01.push_back(Vec3( 1.0, -1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3(1.0, 1.0, 1.0), body1_01);

    auto vit2_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body2_01);
    auto vit2_02 = findVertex(Vec3(-1.0,-1.0, 1.0), body2_01);
    auto vit2_03 = findVertex(Vec3( 1.0,-1.0, 1.0), body2_01);

    auto fit2_01 = findFace(featurePoints2_01, body2_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_01, vit2_02));
    Q_01.push_back(BDVertex(vit1_01, vit2_03));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_VERTEX);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_FACE);
    EXPECT_EQ(info_01.mVit1,  vit1_01);
    EXPECT_EQ(info_01.mFit2,  fit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a2-b1, a3-b1
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test09) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints1_01.push_back(Vec3( 1.0, -1.0,  1.0));
    featurePoints1_01.push_back(Vec3(-1.0, -1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body1_01);
    auto vit1_02 = findVertex(Vec3( 1.0,-1.0, 1.0), body1_01);
    auto vit1_03 = findVertex(Vec3(-1.0,-1.0, 1.0), body1_01);

    auto vit2_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body2_01);

    auto fit1_01 = findFace(featurePoints1_01, body1_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_02, vit2_01));
    Q_01.push_back(BDVertex(vit1_03, vit2_01));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_FACE);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_VERTEX);
    EXPECT_EQ(info_01.mFit1,  fit1_01);
    EXPECT_EQ(info_01.mVit2,  vit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a1-b2, a2-b2 edge-edge
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test10) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints1_01.push_back(Vec3( 1.0, -1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints2_01.push_back(Vec3( 1.0, -1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body1_01);
    auto vit1_02 = findVertex(Vec3( 1.0,-1.0, 1.0), body1_01);

    auto vit2_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body2_01);
    auto vit2_02 = findVertex(Vec3( 1.0,-1.0, 1.0), body2_01);

    auto eit1_01 = findEdge(featurePoints1_01, body1_01);
    auto eit2_01 = findEdge(featurePoints2_01, body2_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_01, vit2_02));
    Q_01.push_back(BDVertex(vit1_02, vit2_02));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_EDGE);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_EDGE);
    EXPECT_EQ(info_01.mEit1,  eit1_01);
    EXPECT_EQ(info_01.mEit2,  eit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a1-b2, a2-b2 face-face
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test11) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints1_01.push_back(Vec3(-1.0, -1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints2_01.push_back(Vec3(-1.0, -1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body1_01);
    auto vit1_02 = findVertex(Vec3(-1.0,-1.0, 1.0), body1_01);

    auto vit2_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body2_01);
    auto vit2_02 = findVertex(Vec3(-1.0,-1.0, 1.0), body2_01);

    auto fit1_01 = findFace(featurePoints1_01, body1_01);
    auto fit2_01 = findFace(featurePoints2_01, body2_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_01, vit2_02));
    Q_01.push_back(BDVertex(vit1_02, vit2_02));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_FACE);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_FACE);
    EXPECT_EQ(info_01.mFit1,  fit1_01);
    EXPECT_EQ(info_01.mFit2,  fit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a2-b1, a2-b2 edge-edge
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test12) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints1_01.push_back(Vec3( 1.0, -1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints2_01.push_back(Vec3( 1.0, -1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body1_01);
    auto vit1_02 = findVertex(Vec3( 1.0,-1.0, 1.0), body1_01);

    auto vit2_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body2_01);
    auto vit2_02 = findVertex(Vec3( 1.0,-1.0, 1.0), body2_01);

    auto eit1_01 = findEdge(featurePoints1_01, body1_01);
    auto eit2_01 = findEdge(featurePoints2_01, body2_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_02, vit2_01));
    Q_01.push_back(BDVertex(vit1_02, vit2_02));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_EDGE);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_EDGE);
    EXPECT_EQ(info_01.mEit1,  eit1_01);
    EXPECT_EQ(info_01.mEit2,  eit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a2-b1, a2-b2 face-face
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test13) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints1_01.push_back(Vec3(-1.0, -1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints2_01.push_back(Vec3(-1.0, -1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body1_01);
    auto vit1_02 = findVertex(Vec3(-1.0,-1.0, 1.0), body1_01);

    auto vit2_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body2_01);
    auto vit2_02 = findVertex(Vec3(-1.0,-1.0, 1.0), body2_01);

    auto fit1_01 = findFace(featurePoints1_01, body1_01);
    auto fit2_01 = findFace(featurePoints2_01, body2_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_02, vit2_01));
    Q_01.push_back(BDVertex(vit1_02, vit2_02));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_FACE);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_FACE);
    EXPECT_EQ(info_01.mFit1,  fit1_01);
    EXPECT_EQ(info_01.mFit2,  fit2_01);

}


/**  @brief findActiveFeaturesFromBinaryDilation()
 *          a1-b1, a2-b2, a3-b3
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test14) {

    ConvexRigidBody body1_01(1);
    ConvexRigidBody body2_01(2);
    ContactPairInfo info_01;

    vector<Vec3> polygon_01;
    polygon_01.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_01.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_01.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_01.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints1_01;
    featurePoints1_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints1_01.push_back(Vec3( 1.0, -1.0,  1.0));
    featurePoints1_01.push_back(Vec3(-1.0, -1.0,  1.0));

    vector<Vec3> polygon_02;
    polygon_02.push_back(Vec3( 1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0,  1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0,  1.0));
    polygon_02.push_back(Vec3( 1.0,  1.0, -1.0));
    polygon_02.push_back(Vec3( 1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0, -1.0, -1.0));
    polygon_02.push_back(Vec3(-1.0,  1.0, -1.0));
    vector<Vec3> featurePoints2_01;
    featurePoints2_01.push_back(Vec3( 1.0,  1.0,  1.0));
    featurePoints2_01.push_back(Vec3( 1.0, -1.0,  1.0));
    featurePoints2_01.push_back(Vec3(-1.0, -1.0,  1.0));

    generateTestRigidBodies(
        polygon_01,
        Vec3(0.0, 0.0, 0.0),
        polygon_02,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        0.0,
        body1_01,
        body2_01,
        0.0,
        false
    );

    auto vit1_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body1_01);
    auto vit1_02 = findVertex(Vec3( 1.0,-1.0, 1.0), body1_01);
    auto vit1_03 = findVertex(Vec3(-1.0,-1.0, 1.0), body1_01);

    auto vit2_01 = findVertex(Vec3( 1.0, 1.0, 1.0), body2_01);
    auto vit2_02 = findVertex(Vec3( 1.0,-1.0, 1.0), body2_01);
    auto vit2_03 = findVertex(Vec3(-1.0,-1.0, 1.0), body2_01);

    auto fit1_01 = findFace(featurePoints1_01, body1_01);
    auto fit2_01 = findFace(featurePoints2_01, body2_01);

    vector<BDVertex> Q_01;
    Q_01.push_back(BDVertex(vit1_01, vit2_01));
    Q_01.push_back(BDVertex(vit1_02, vit2_02));
    Q_01.push_back(BDVertex(vit1_03, vit2_03));

    double epsilonZero  = 1.0e-10;
    double epsilonAngle = 0.1;

    ContactPointsAndNormalGenerator 
        gen_01(body1_01, body2_01, info_01, Q_01, 
            epsilonZero, epsilonAngle, true, 6, std::cerr);
    gen_01.findActiveFeaturesFromBinaryDilation();

    EXPECT_EQ(info_01.mType1, ContactPairInfo::FT_FACE);
    EXPECT_EQ(info_01.mType2, ContactPairInfo::FT_FACE);
    EXPECT_EQ(info_01.mFit1,  fit1_01);
    EXPECT_EQ(info_01.mFit2,  fit2_01);


}


/**  @brief adjustFeatures_FACE_EDGE() no penetration
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test15) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = 0.0;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1, expectedIntsec_2_1 );


        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }

}


/**  @brief adjustFeatures_FACE_EDGE() tilting no penetration
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test16) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double     yawTiltAxis = (rand100()/50.0 - 1.0)*0.2;
        Vec3       tiltAxis(-1.0 * cos(yawTiltAxis), sin(yawTiltAxis), 0.0);
        double     rotationRad = (rand100()/50.0  - 1.0) * 0.9 * atan(1.0/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1, expectedIntsec_2_1 );


        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }

}


/**  @brief adjustFeatures_FACE_EDGE() penetration
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test17) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(1.0/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }

}


/**  @brief adjustFeatures_FACE_EDGE() penetration overshoot
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test18) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(1.0/2.0) + (M_PI/2.0)*rand100()/100.0;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }

}


/**  @brief adjustFeatures_FACE_EDGE() penetration
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test19) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(1.0, 0.0, 0.0);
        double     rotationRad = atan(1.0/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }

}


/**  @brief adjustFeatures_FACE_EDGE() penetration overshoot
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test20) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(1.0, 0.0, 0.0);
        double     rotationRad = atan(1.0/2.0) + (M_PI/2.0)*rand100()/100.0;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);
    }
}


/**  @brief adjustFeatures_FACE_VERTEX() no penetration
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test21) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double     yawTiltAxis = 2.0 * M_PI * rand100()/100.0 - M_PI;
        Vec3       tiltAxis(cos(yawTiltAxis), sin(yawTiltAxis), 0.0);
        double     rotationRad = atan(1.0/sqrt(8.0)) * 0.9 *(rand100()/50.0 - 1.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_VERTEX, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_FACE_VERTEX() no penetration
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test22) {

    for (long i = 0; i < 20; i++) {


        bool swap12 = ((i%2)==0);
        long type = long((rand100()/4.0))%4;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));
        if (type==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else if (type==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        vector<Vec3> expectedIntsec_2_1;


        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis1( 1.0, 1.0, 0.0);
        Vec3       tiltAxis2(-1.0, 1.0, 0.0);
        Vec3       tiltAxis3(-1.0,-1.0, 0.0);
        Vec3       tiltAxis4( 1.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type==0) {
            tiltAxis = tiltAxis1;
        }
        else if (type==1) {
            tiltAxis = tiltAxis2;
        }
        else if (type==2) {
            tiltAxis = tiltAxis3;
        }
        else {
            tiltAxis = tiltAxis4;
        }

        double     rotationRad = atan(1.0/sqrt(8.0));

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }

}


/**  @brief adjustFeatures_FACE_VERTEX() no penetration
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test23) {

    for (long i = 0; i < 20; i++) {


        bool swap12 = ((i%2)==0);
        long type = long((rand100()/4.0))%4;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));
        if (type==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else if (type==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        vector<Vec3> expectedIntsec_2_1;


        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis1( 1.0, 1.0, 0.0);
        Vec3       tiltAxis2(-1.0, 1.0, 0.0);
        Vec3       tiltAxis3(-1.0,-1.0, 0.0);
        Vec3       tiltAxis4( 1.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type==0) {
            tiltAxis = tiltAxis1;
        }
        else if (type==1) {
            tiltAxis = tiltAxis2;
        }
        else if (type==2) {
            tiltAxis = tiltAxis3;
        }
        else {
            tiltAxis = tiltAxis4;
        }

        double     rotationRad = atan(1.0/sqrt(8.0)) + 
                                  M_PI * (rand100()/100.0) /4.0;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }

}


/**  @brief adjustFeatures_FACE_VERTEX() no penetration
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test24) {

    for (long i = 0; i < 20; i++) {


        bool swap12 = ((i%2)==0);
        long type = long((rand100()/4.0))%4;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));
        if (type==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else if (type==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        vector<Vec3> expectedIntsec_2_1;


        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis1( 1.0, 0.0, 0.0);
        Vec3       tiltAxis2( 0.0, 1.0, 0.0);
        Vec3       tiltAxis3(-1.0, 0.0, 0.0);
        Vec3       tiltAxis4( 0.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type==0) {
            tiltAxis = tiltAxis1;
        }
        else if (type==1) {
            tiltAxis = tiltAxis2;
        }
        else if (type==2) {
            tiltAxis = tiltAxis3;
        }
        else {
            tiltAxis = tiltAxis4;
        }

        double     rotationRad = atan(1.0/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }

}


/**  @brief adjustFeatures_FACE_VERTEX() no penetration
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test25) {

    for (long i = 0; i < 20; i++) {


        bool swap12 = ((i%2)==0);
        long type = long((rand100()/4.0))%4;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));
        if (type==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else if (type==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        vector<Vec3> expectedIntsec_2_1;


        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis1( 1.0, 0.0, 0.0);
        Vec3       tiltAxis2( 0.0, 1.0, 0.0);
        Vec3       tiltAxis3(-1.0, 0.0, 0.0);
        Vec3       tiltAxis4( 0.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type==0) {
            tiltAxis = tiltAxis1;
        }
        else if (type==1) {
            tiltAxis = tiltAxis2;
        }
        else if (type==2) {
            tiltAxis = tiltAxis3;
        }
        else {
            tiltAxis = tiltAxis4;
        }

        double     rotationRad = atan(1.0/sqrt(2.0)) + 
                                  M_PI * (rand100()/100.0) /4.0;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_CROSSING() 
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test26) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 0.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 0.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 0.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double     yawTiltAxis = 2.0 * M_PI * rand100()/100.0 - M_PI;
        Vec3       tiltAxis(cos(yawTiltAxis), sin(yawTiltAxis), 0.0);
        double     rotationRad = atan(1.0/2.0) * 0.9 * (rand100()/50.0 - 1.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_CROSSING() 
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test27) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type   = long(rand100()/2.0)%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 0.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 0.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 0.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        if (type==0) {       
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(1.0, 0.0, 0.0);
        double     rotationRad = (type==0)?atan(1.0/2.0):(-1.0*atan(1.0/2.0));

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_CROSSING() 
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test28) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type   = long(rand100()/2.0)%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 0.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 0.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 0.0, -2.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        if (type==0) {       
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(1.0, 0.0, 0.0);
        double     rotationRad = (type==0)?
                                     (atan(1.0/2.0) + M_PI*rand100()/400.0):
                                     (-1.0*(atan(1.0/2.0)+ M_PI*rand100()/400.0));

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_CROSSING() 
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test29) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        long type   = long(rand100()/4.0)%4;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 0.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  2.0,  1.0));
        activePoints_1.push_back(Vec3( 0.0, -2.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  2.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 0.0, -2.0,  1.0));
        if (type==0) {       
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type==1) {
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else if (type==2) {
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        if (type==0) {       
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else if (type==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else if (type==2) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       qAxis1, qAxis2;
        if (type==0) {
            qAxis1 = Vec3(-1.0, 0.0, 0.0);
            qAxis2 = Vec3( 0.0, 1.0, 0.0);
        }
        else if (type==1) {
            qAxis1 = Vec3(-1.0, 0.0, 0.0);
            qAxis2 = Vec3( 0.0,-1.0, 0.0);
        }
        else if (type==2) {
            qAxis1 = Vec3( 1.0, 0.0, 0.0);
            qAxis2 = Vec3( 0.0,-1.0, 0.0);
        }
        else {
            qAxis1 = Vec3( 1.0, 0.0, 0.0);
            qAxis2 = Vec3( 0.0, 1.0, 0.0);
        }
        
        double rotationRadTmp = atan(1.0/2.0);
        Quaternion q01(qAxis1, rotationRadTmp);
        Quaternion q02(qAxis2, rotationRadTmp);
        Quaternion q03 = q02 * q01;
        q03.normalize();
        double rotationRad = 2.0*acos(q03.s()) + M_PI * (rand100()/100.0)*0.25 ;
        auto tiltAxis = Vec3(q03.x(),q03.y(),q03.z());

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_PARALLEL()
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test30) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));


        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));


        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = atan(1.0/2.0)*2.0 + M_PI*0.5*(rand100()/100.0);
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_PARALLEL()
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test31) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));


        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));


        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = atan(1.0/2.0)*2.0 + M_PI*0.5*(rand100()/100.0);
        Vec3 tiltAxis(1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_PARALLEL() thin
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test32) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));


        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  0.2,  0.0));
        polygon_2.push_back(Vec3(-2.0,  0.2,  0.0));
        polygon_2.push_back(Vec3( 2.0, -0.2,  0.0));
        polygon_2.push_back(Vec3(-2.0, -0.2,  0.0));
        polygon_2.push_back(Vec3( 2.0,  0.2,  1.0));
        polygon_2.push_back(Vec3(-2.0,  0.2,  1.0));
        polygon_2.push_back(Vec3( 2.0, -0.2,  1.0));
        polygon_2.push_back(Vec3(-2.0, -0.2,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.2,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.2,  0.0));


        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad =
           std::min(atan(1.0/2.0) + atan(1.0/0.2)+rand100()/100.0, M_PI*0.95);
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_PARALLEL() thin
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test33) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));


        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  0.2,  0.0));
        polygon_2.push_back(Vec3(-2.0,  0.2,  0.0));
        polygon_2.push_back(Vec3( 2.0, -0.2,  0.0));
        polygon_2.push_back(Vec3(-2.0, -0.2,  0.0));
        polygon_2.push_back(Vec3( 2.0,  0.2,  1.0));
        polygon_2.push_back(Vec3(-2.0,  0.2,  1.0));
        polygon_2.push_back(Vec3( 2.0, -0.2,  1.0));
        polygon_2.push_back(Vec3(-2.0, -0.2,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0, -0.2,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0, -0.2,  0.0));


        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad =
           std::min(atan(1.0/2.0) + atan(1.0/0.2)+rand100()/100.0, M_PI*0.95);
        Vec3 tiltAxis(1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

//        showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_PARALLEL() thin-thin non penetrating
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test34) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  0.2,  0.0));
        polygon_1.push_back(Vec3(-2.0,  0.2,  0.0));
        polygon_1.push_back(Vec3( 2.0, -0.2,  0.0));
        polygon_1.push_back(Vec3(-2.0, -0.2,  0.0));
        polygon_1.push_back(Vec3( 2.0,  0.2, -1.0));
        polygon_1.push_back(Vec3(-2.0,  0.2, -1.0));
        polygon_1.push_back(Vec3( 2.0, -0.2, -1.0)); 
       polygon_1.push_back(Vec3(-2.0, -0.2, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  0.2,  0.0));
        polygon_2.push_back(Vec3(-2.0,  0.2,  0.0));
        polygon_2.push_back(Vec3( 2.0, -0.2,  0.0));
        polygon_2.push_back(Vec3(-2.0, -0.2,  0.0));
        polygon_2.push_back(Vec3( 2.0,  0.2,  1.0));
        polygon_2.push_back(Vec3(-2.0,  0.2,  1.0));
        polygon_2.push_back(Vec3( 2.0, -0.2,  1.0));
        polygon_2.push_back(Vec3(-2.0, -0.2,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad =atan(1.0/0.2)*2.0 *0.95 * (rand100()/50.0 -1.0);

        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_PARALLEL() thin-thin penetrating
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test35) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  0.2,  0.0));
        polygon_1.push_back(Vec3(-2.0,  0.2,  0.0));
        polygon_1.push_back(Vec3( 2.0, -0.2,  0.0));
        polygon_1.push_back(Vec3(-2.0, -0.2,  0.0));
        polygon_1.push_back(Vec3( 2.0,  0.2, -1.0));
        polygon_1.push_back(Vec3(-2.0,  0.2, -1.0));
        polygon_1.push_back(Vec3( 2.0, -0.2, -1.0)); 
       polygon_1.push_back(Vec3(-2.0, -0.2, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.2,  0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.2,  0.0));


        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  0.2,  0.0));
        polygon_2.push_back(Vec3(-2.0,  0.2,  0.0));
        polygon_2.push_back(Vec3( 2.0, -0.2,  0.0));
        polygon_2.push_back(Vec3(-2.0, -0.2,  0.0));
        polygon_2.push_back(Vec3( 2.0,  0.2,  1.0));
        polygon_2.push_back(Vec3(-2.0,  0.2,  1.0));
        polygon_2.push_back(Vec3( 2.0, -0.2,  1.0));
        polygon_2.push_back(Vec3(-2.0, -0.2,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.2,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.2,  0.0));


        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad =std::min(atan(1.0/0.2)*2.0+rand100()/100.0, M_PI*0.99);

        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_PARALLEL() thin-thin penetrating
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test36) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  0.2,  0.0));
        polygon_1.push_back(Vec3(-2.0,  0.2,  0.0));
        polygon_1.push_back(Vec3( 2.0, -0.2,  0.0));
        polygon_1.push_back(Vec3(-2.0, -0.2,  0.0));
        polygon_1.push_back(Vec3( 2.0,  0.2, -1.0));
        polygon_1.push_back(Vec3(-2.0,  0.2, -1.0));
        polygon_1.push_back(Vec3( 2.0, -0.2, -1.0)); 
       polygon_1.push_back(Vec3(-2.0, -0.2, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -0.2,  0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -0.2,  0.0));


        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  0.2,  0.0));
        polygon_2.push_back(Vec3(-2.0,  0.2,  0.0));
        polygon_2.push_back(Vec3( 2.0, -0.2,  0.0));
        polygon_2.push_back(Vec3(-2.0, -0.2,  0.0));
        polygon_2.push_back(Vec3( 2.0,  0.2,  1.0));
        polygon_2.push_back(Vec3(-2.0,  0.2,  1.0));
        polygon_2.push_back(Vec3( 2.0, -0.2,  1.0));
        polygon_2.push_back(Vec3(-2.0, -0.2,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0, -0.2,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0, -0.2,  0.0));


        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad =std::min(atan(1.0/0.2)*2.0+rand100()/100.0, M_PI*0.99);

        Vec3 tiltAxis(1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_PARALLEL() fat-fat non-penetrating
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test37) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,    0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,    0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  200.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  200.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -200.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -200.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  200.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  200.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -200.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -200.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,    0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,    0.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,    0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,    0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  200.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  200.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -200.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -200.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  200.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  200.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -200.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -200.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,    0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,    0.0, -1.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;

        Vec3 tiltAxis(1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_PARALLEL() fat-fat penetrating
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test38) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,    0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,    0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  200.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  200.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -200.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -200.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  200.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  200.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -200.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -200.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,      0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,      0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,    200.0,  0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,    200.0,  0.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,    0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,    0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  200.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  200.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -200.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -200.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  200.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  200.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -200.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -200.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,      0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,      0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,    200.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,    200.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = std::min(atan(1.0/200.0)*2.0 + M_PI*rand100()/100.0, M_PI*0.99);

        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_EDGE_PARALLEL() fat-fat penetrating
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test39) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,    0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,    0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  200.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  200.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -200.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -200.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  200.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  200.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -200.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -200.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,      0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,      0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,   -200.0,  0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,   -200.0,  0.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,    0.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,    0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  200.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  200.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -200.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -200.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  200.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  200.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -200.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -200.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,      0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,      0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,   -200.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,   -200.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = std::min(atan(1.0/200.0)*2.0 + M_PI*rand100()/100.0, M_PI*0.99);

        Vec3 tiltAxis(1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);

    }
}


/**  @brief adjustFeatures_EDGE_VERTEX no peoetration
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test40) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,      0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,      0.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = atan(1.0/2.0) * 0.99 * rand100()/100.0;
        double yawTiltAxis = 2.0 * M_PI * rand100()/100.0 - M_PI;
        Vec3 tiltAxis(cos(yawTiltAxis), sin(yawTiltAxis), 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_EDGE,   expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_VERTEX, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief adjustFeatures_EDGE_VERTEX edge-face 
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test41) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,      0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,      0.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = atan(1.0/2.0);
        Vec3 tiltAxis(0.0, 1.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief adjustFeatures_EDGE_VERTEX face-face 
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test42) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = atan(1.0/2.0)*2.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief adjustFeatures_EDGE_VERTEX face-edge
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test43) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = M_PI*0.25;

        double rotationRad = atan(1.0/2.0) + atan(1.0/sqrt(8.0));
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief adjustFeatures_EDGE_VERTEX edge-edge
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test44) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = M_PI*0.25;

        double rotationRad = atan(1.0/sqrt(8.0));
        Vec3 tiltAxis( 0.0, 1.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1, expectedIntsec_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief adjustFeatures_VERTEX_VERTEX no penetration
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test45) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0,  1.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = M_PI*0.25;

        double rotationRad = 0.99 * atan(1.0/sqrt(8.0)) * rand100()/100.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_VERTEX, expectedPoints_1_1, expectedIntsec_1_1,
          ContactPairInfo::FT_VERTEX, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief adjustFeatures_VERTEX_VERTEX face-face
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test46) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  1.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0,  1.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = atan(1.0/2.0) * 2.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
          ContactPairInfo::FT_FACE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief adjustFeatures_VERTEX_VERTEX edge-edge
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test47) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  0.5));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0, 0.5));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0,  0.5));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  0.5);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -0.5));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -0.5));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -0.5));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -0.5);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = atan(0.5/sqrt(8.0)) * 2.0;
        Vec3 tiltAxis(-1.0, 1.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_EDGE, expectedPoints_1_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief adjustFeatures_VERTEX_VERTEX face-edge
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test48) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  0.5));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0, 0.5));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0,  0.5));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));

        vector<Vec3> expectedIntsec_1_1;

        Vec3       anchor_1(0.0,  0.0,  0.5);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -0.5));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -0.5));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -0.5));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;

        Vec3       anchor_2(0.0,  0.0, -0.5);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = M_PI*0.25;

        double rotationRad = atan(0.5/sqrt(8.0)) + atan(0.5/2.0);
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_FACE, expectedPoints_1_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2_1, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.adjustFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief process_FACE_FACE
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test49) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  1.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0, 1.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0, 1.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0, 1.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0, 1.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0,  2.0, 1.0));
        expectedPoints_1.push_back(Vec3(-2.0,  2.0, 1.0));
        expectedPoints_1.push_back(Vec3( 2.0, -2.0, 1.0));
        expectedPoints_1.push_back(Vec3(-2.0, -2.0, 1.0));

        vector<Vec3> expectedIntsec_1;
        expectedIntsec_1.push_back(Vec3( 2.0,  2.0, 1.0));
        expectedIntsec_1.push_back(Vec3(-2.0,  2.0, 1.0));
        expectedIntsec_1.push_back(Vec3( 2.0, -2.0, 1.0));
        expectedIntsec_1.push_back(Vec3(-2.0, -2.0, 1.0));


        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  2.0,-1.0));
        activePoints_2.push_back(Vec3(-2.0,  2.0,-1.0));
        activePoints_2.push_back(Vec3( 2.0, -2.0,-1.0));
        activePoints_2.push_back(Vec3(-2.0, -2.0,-1.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 2.0,  2.0,-1.0));
        expectedPoints_2.push_back(Vec3(-2.0,  2.0,-1.0));
        expectedPoints_2.push_back(Vec3( 2.0, -2.0,-1.0));
        expectedPoints_2.push_back(Vec3(-2.0, -2.0,-1.0));

        vector<Vec3> expectedIntsec_2;
        expectedIntsec_2.push_back(Vec3( 2.0,  2.0,-1.0));
        expectedIntsec_2.push_back(Vec3(-2.0,  2.0,-1.0));
        expectedIntsec_2.push_back(Vec3( 2.0, -2.0,-1.0));
        expectedIntsec_2.push_back(Vec3(-2.0, -2.0,-1.0));

        Vec3       anchor_2(0.0,  0.0,-1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_FACE, expectedPoints_1, expectedIntsec_1,
          ContactPairInfo::FT_FACE, expectedPoints_2, expectedIntsec_2 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_FACE_NORMAL_1);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3(0.0, 0.0,-1.0)));
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3(0.0, 0.0, 1.0)));
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_FACE_FACE
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test50) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0, -3.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -3.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -3.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -3.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -5.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -5.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -5.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -5.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  2.0, -3.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0, -3.0));
        activePoints_1.push_back(Vec3( 2.0, -2.0, -3.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0, -3.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0,  2.0, -3.0));
        expectedPoints_1.push_back(Vec3(-2.0,  2.0, -3.0));
        expectedPoints_1.push_back(Vec3( 2.0, -2.0, -3.0));
        expectedPoints_1.push_back(Vec3(-2.0, -2.0, -3.0));

        vector<Vec3> expectedIntsec_1;
        expectedIntsec_1.push_back(Vec3( 2.0,  2.0, -3.0));
        expectedIntsec_1.push_back(Vec3(-2.0,  2.0, -3.0));
        expectedIntsec_1.push_back(Vec3( 2.0, -2.0, -3.0));
        expectedIntsec_1.push_back(Vec3(-2.0, -2.0, -3.0));


        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  6.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  6.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  6.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  6.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  2.0, 2.0));
        activePoints_2.push_back(Vec3(-2.0,  2.0, 2.0));
        activePoints_2.push_back(Vec3( 2.0, -2.0, 6.0));
        activePoints_2.push_back(Vec3(-2.0, -2.0, 6.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 2.0,  2.0, 2.0));
        expectedPoints_2.push_back(Vec3(-2.0,  2.0, 2.0));
        expectedPoints_2.push_back(Vec3( 2.0, -2.0, 6.0));
        expectedPoints_2.push_back(Vec3(-2.0, -2.0, 6.0));

        vector<Vec3> expectedIntsec_2;
        expectedIntsec_2.push_back(Vec3( 2.0,  2.0, 2.0));
        expectedIntsec_2.push_back(Vec3(-2.0,  2.0, 2.0));
        expectedIntsec_2.push_back(Vec3( 2.0, -2.0, 6.0));
        expectedIntsec_2.push_back(Vec3(-2.0, -2.0, 6.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_FACE, expectedPoints_1, expectedIntsec_1,
          ContactPairInfo::FT_FACE, expectedPoints_2, expectedIntsec_2 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        if (swap12) {
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_FACE_NORMAL_2);
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3(0.0, 0.0,-1.0)));
        }
        else {
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_FACE_NORMAL_1);
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3(0.0, 0.0, 1.0)));
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_FACE_FACE relative velocity
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test51) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0, -2.0, -3.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -3.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -7.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -7.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -7.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -7.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0, -2.0, -3.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0, -3.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0, -7.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0, -7.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0, -2.0, -3.0));
        expectedPoints_1.push_back(Vec3(-2.0, -2.0, -3.0));
        expectedPoints_1.push_back(Vec3( 2.0,  2.0, -7.0));
        expectedPoints_1.push_back(Vec3(-2.0,  2.0, -7.0));

        vector<Vec3> expectedIntsec_1;
        expectedIntsec_1.push_back(Vec3( 2.0, -2.0, -3.0));
        expectedIntsec_1.push_back(Vec3(-2.0, -2.0, -3.0));
        expectedIntsec_1.push_back(Vec3( 2.0,  2.0, -7.0));
        expectedIntsec_1.push_back(Vec3(-2.0,  2.0, -7.0));


        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 1.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  6.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  6.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  6.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  6.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  2.0, 2.0));
        activePoints_2.push_back(Vec3(-2.0,  2.0, 2.0));
        activePoints_2.push_back(Vec3( 2.0, -2.0, 6.0));
        activePoints_2.push_back(Vec3(-2.0, -2.0, 6.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 2.0,  2.0, 2.0));
        expectedPoints_2.push_back(Vec3(-2.0,  2.0, 2.0));
        expectedPoints_2.push_back(Vec3( 2.0, -2.0, 6.0));
        expectedPoints_2.push_back(Vec3(-2.0, -2.0, 6.0));

        vector<Vec3> expectedIntsec_2;
        expectedIntsec_2.push_back(Vec3( 2.0,  2.0, 2.0));
        expectedIntsec_2.push_back(Vec3(-2.0,  2.0, 2.0));
        expectedIntsec_2.push_back(Vec3( 2.0, -2.0, 6.0));
        expectedIntsec_2.push_back(Vec3(-2.0, -2.0, 6.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0,-1.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_FACE, expectedPoints_1, expectedIntsec_1,
          ContactPairInfo::FT_FACE, expectedPoints_2, expectedIntsec_2 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3(0.0, 0.0,-1.0)));
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3(0.0, 0.0, 1.0)));
        }
        EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_RELATIVE_VELOCITY);
        //showResult(body1, body2, info);
    }
}


/**  @brief process_FACE_FACE projection
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test52) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0, -2.0, -3.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -3.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -7.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -7.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -7.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -7.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0, -2.0, -3.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0, -3.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0, -7.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0, -7.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0, -2.0, -3.0));
        expectedPoints_1.push_back(Vec3(-2.0, -2.0, -3.0));
        expectedPoints_1.push_back(Vec3( 2.0,  2.0, -7.0));
        expectedPoints_1.push_back(Vec3(-2.0,  2.0, -7.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 2.0, -2.5, -2.5));
        expectedIntsec_1_1.push_back(Vec3(-2.0, -2.5, -2.5));
        expectedIntsec_1_1.push_back(Vec3( 2.0, -6.5,  1.5));
        expectedIntsec_1_1.push_back(Vec3(-2.0, -6.5,  1.5));

        vector<Vec3> expectedIntsec_1_2;
        expectedIntsec_1_2.push_back(Vec3( 2.0, -2.0, -3.0));
        expectedIntsec_1_2.push_back(Vec3(-2.0, -2.0, -3.0));
        expectedIntsec_1_2.push_back(Vec3( 2.0,  2.0, -7.0));
        expectedIntsec_1_2.push_back(Vec3(-2.0,  2.0, -7.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  6.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  6.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  6.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  6.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  2.0, 2.0));
        activePoints_2.push_back(Vec3(-2.0,  2.0, 2.0));
        activePoints_2.push_back(Vec3( 2.0, -2.0, 6.0));
        activePoints_2.push_back(Vec3(-2.0, -2.0, 6.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 2.0,  2.0, 2.0));
        expectedPoints_2.push_back(Vec3(-2.0,  2.0, 2.0));
        expectedPoints_2.push_back(Vec3( 2.0, -2.0, 6.0));
        expectedPoints_2.push_back(Vec3(-2.0, -2.0, 6.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedIntsec_2_1.push_back(Vec3(-2.0,  2.0,  2.0));
        expectedIntsec_2_1.push_back(Vec3( 2.0, -2.0,  6.0));
        expectedIntsec_2_1.push_back(Vec3(-2.0, -2.0,  6.0));

        vector<Vec3> expectedIntsec_2_2;
        expectedIntsec_2_2.push_back(Vec3( 2.0,  2.5,  1.5));
        expectedIntsec_2_2.push_back(Vec3(-2.0,  2.5,  1.5));
        expectedIntsec_2_2.push_back(Vec3( 2.0,  6.5, -2.5));
        expectedIntsec_2_2.push_back(Vec3(-2.0,  6.5, -2.5));


        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_FACE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_FACE, expectedPoints_2, expectedIntsec_2_1 );
        ExpStructContactPointsAndNormalGenerator exp2(
          ContactPairInfo::FT_FACE, expectedPoints_1, expectedIntsec_1_2,
          ContactPairInfo::FT_FACE, expectedPoints_2, expectedIntsec_2_2 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12)||
                  checkResult(body1, body2, info, exp2, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  -1.0*sqrt(0.5),  -1.0*sqrt(0.5))));
            if (info.mContactNormalType == ContactPairInfo::NT_FACE_NORMAL_2) {
                EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
            }
            else if (info.mContactNormalType == ContactPairInfo::NT_FACE_NORMAL_1) {
                EXPECT_EQ(checkResult(body1, body2, info, exp2, swap12), true);
            }
            else {
                EXPECT_EQ(false, true);
            }
        }

        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,   1.0*sqrt(0.5),  1.0*sqrt(0.5))));
            if (info.mContactNormalType == ContactPairInfo::NT_FACE_NORMAL_2) {
                EXPECT_EQ(checkResult(body1, body2, info, exp2, swap12), true);
            }
            else if (info.mContactNormalType == ContactPairInfo::NT_FACE_NORMAL_1) {
                EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
            }
            else {
                EXPECT_EQ(false, true);
            }
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_FACE_EDGE face normal
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test53) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0, -2.0,  2.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0,  2.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0,  2.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0, -2.0,  2.0));
        expectedPoints_1.push_back(Vec3(-2.0, -2.0,  2.0));
        expectedPoints_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedPoints_1.push_back(Vec3(-2.0,  2.0,  2.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 2.0,  0.0, 2.0));
        expectedIntsec_1_1.push_back(Vec3(-2.0,  0.0, 2.0));

        Vec3       anchor_1(0.0,  0.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  0.0, -2.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -2.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 2.0,  0.0, -2.0));
        expectedPoints_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 2.0,  0.0, -2.0));
        expectedIntsec_2_1.push_back(Vec3(-2.0,  0.0, -2.0));

        Vec3       anchor_2(0.0,  0.0, -2.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_FACE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0,-1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_FACE_NORMAL_2);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_FACE_NORMAL_1);
        }
//        showResult(body1, body2, info);
    }
}


/**  @brief process_FACE_EDGE edge normal
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test54) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);


        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -6.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -6.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -6.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -6.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0, -2.0, -2.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0, -2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0, -6.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0, -6.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0, -2.0, -2.0));
        expectedPoints_1.push_back(Vec3(-2.0, -2.0, -2.0));
        expectedPoints_1.push_back(Vec3( 2.0,  2.0, -6.0));
        expectedPoints_1.push_back(Vec3(-2.0,  2.0, -6.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 2.0,  0.0, -4.0));
        expectedIntsec_1_1.push_back(Vec3(-2.0,  0.0, -4.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 1.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  6.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  6.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  6.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  6.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  4.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  4.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  4.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  4.0));
        polygon_2.push_back(Vec3( 2.0,  0.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0,  2.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedPoints_2.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 2.0,  0.0, 2.0));
        expectedIntsec_2_1.push_back(Vec3(-2.0,  0.0, 2.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0,-1.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_FACE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0,-1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_NORMAL_1);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_NORMAL_2);
        }
//        showResult(body1, body2, info);
    }
}


/**  @brief process_FACE_EDGE relative velocity
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test55) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -6.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -6.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -6.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -6.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0, -2.0, -2.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0, -2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0, -6.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0, -6.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0, -2.0, -2.0));
        expectedPoints_1.push_back(Vec3(-2.0, -2.0, -2.0));
        expectedPoints_1.push_back(Vec3( 2.0,  2.0, -6.0));
        expectedPoints_1.push_back(Vec3(-2.0,  2.0, -6.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 2.0,  0.0, -4.0));
        expectedIntsec_1_1.push_back(Vec3(-2.0,  0.0, -4.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 1.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  4.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  4.0));
        polygon_2.push_back(Vec3( 2.0,  0.0,  4.0));
        polygon_2.push_back(Vec3(-2.0,  0.0,  4.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  0.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0,  2.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedPoints_2.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 2.0,  0.0, 2.0));
        expectedIntsec_2_1.push_back(Vec3(-2.0,  0.0, 2.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0,-1.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_FACE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0,-1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_RELATIVE_VELOCITY);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_RELATIVE_VELOCITY);
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_FACE_EDGE projection
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test56) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -6.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -6.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -6.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -6.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0, -2.0, -2.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0, -2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0, -6.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0, -6.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0, -2.0, -2.0));
        expectedPoints_1.push_back(Vec3(-2.0, -2.0, -2.0));
        expectedPoints_1.push_back(Vec3( 2.0,  2.0, -6.0));
        expectedPoints_1.push_back(Vec3(-2.0,  2.0, -6.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 2.0, -3.0, -1.0));
        expectedIntsec_1_1.push_back(Vec3(-2.0, -3.0, -1.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  4.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  4.0));
        polygon_2.push_back(Vec3( 2.0,  0.0,  4.0));
        polygon_2.push_back(Vec3(-2.0,  0.0,  4.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  0.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0,  2.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedPoints_2.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 2.0,  0.0, 2.0));
        expectedIntsec_2_1.push_back(Vec3(-2.0,  0.0, 2.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_FACE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  -1.0*sqrt(0.5), -1.0*sqrt(0.5))));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_FACE_NORMAL_2);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  sqrt(0.5), sqrt(0.5))));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_FACE_NORMAL_1);
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_FACE_VERTEX
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test57) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -6.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -6.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -6.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -6.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0, -2.0, -2.0));
        activePoints_1.push_back(Vec3(-2.0, -2.0, -2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0, -2.0));
        activePoints_1.push_back(Vec3(-2.0,  2.0, -2.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0, -2.0, -2.0));
        expectedPoints_1.push_back(Vec3(-2.0, -2.0, -2.0));
        expectedPoints_1.push_back(Vec3( 2.0,  2.0, -2.0));
        expectedPoints_1.push_back(Vec3(-2.0,  2.0, -2.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 0.0,  0.0, -2.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  4.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  4.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  4.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  4.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 0.0,  0.0,  0.0));


        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0,  0.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 0.0,  0.0,  0.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 0.0,  0.0, 0.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_FACE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_VERTEX, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, -1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_FACE_NORMAL_2);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_FACE_NORMAL_1);
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_EDGE_EDGE CROSSING
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test58) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0)); 
        polygon_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedPoints_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 0.0,  0.0,  2.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -2.0)); 
        polygon_2.push_back(Vec3( 0.0,  2.0, -2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -2.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 0.0,  2.0, -2.0));
        expectedPoints_2.push_back(Vec3( 0.0, -2.0, -2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 0.0,  0.0,-2.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_EDGE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, -1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_CROSS_EDGE);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_CROSS_EDGE);
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_EDGE_EDGE PARALLEL 1
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test59) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0)); 
        polygon_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedPoints_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedIntsec_1_1.push_back(Vec3(-2.0,  0.0,  2.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  0.0, -2.0)); 
        polygon_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -2.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 2.0,  0.0, -2.0));
        expectedPoints_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 2.0,  0.0,-2.0));
        expectedIntsec_2_1.push_back(Vec3(-2.0,  0.0,-2.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_EDGE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, -1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_EDGE_EDGE PARALLEL 2
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test60) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0)); 
        polygon_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedPoints_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 1.0,  0.0,  2.0));
        expectedIntsec_1_1.push_back(Vec3(-1.0,  0.0,  2.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 1.0,  0.0, -2.0)); 
        polygon_2.push_back(Vec3(-1.0,  0.0, -2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  0.0, -2.0));
        activePoints_2.push_back(Vec3(-1.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 1.0,  0.0, -2.0));
        expectedPoints_2.push_back(Vec3(-1.0,  0.0, -2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 1.0,  0.0,-2.0));
        expectedIntsec_2_1.push_back(Vec3(-1.0,  0.0,-2.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_EDGE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, -1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_EDGE_EDGE PARALLEL 3
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test61) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0)); 
        polygon_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedPoints_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedIntsec_1_1.push_back(Vec3(-1.0,  0.0,  2.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  0.0, -2.0)); 
        polygon_2.push_back(Vec3(-1.0,  0.0, -2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -2.0));
        activePoints_2.push_back(Vec3(-1.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 2.0,  0.0, -2.0));
        expectedPoints_2.push_back(Vec3(-1.0,  0.0, -2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 2.0,  0.0,-2.0));
        expectedIntsec_2_1.push_back(Vec3(-1.0,  0.0,-2.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_EDGE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, -1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_EDGE_EDGE PARALLEL 4
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test62) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0)); 
        polygon_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedPoints_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 1.0,  0.0,  2.0));
        expectedIntsec_1_1.push_back(Vec3(-2.0,  0.0,  2.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 1.0,  0.0, -2.0)); 
        polygon_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  0.0, -2.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 1.0,  0.0, -2.0));
        expectedPoints_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 1.0,  0.0,-2.0));
        expectedIntsec_2_1.push_back(Vec3(-2.0,  0.0,-2.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_EDGE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, -1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_EDGE_EDGE PARALLEL 5
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test63) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0)); 
        polygon_1.push_back(Vec3(-1.0,  0.0,  2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3(-1.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedPoints_1.push_back(Vec3(-1.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 1.0,  0.0,  2.0));
        expectedIntsec_1_1.push_back(Vec3(-1.0,  0.0,  2.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 1.0,  0.0, -2.0)); 
        polygon_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  0.0, -2.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 1.0,  0.0, -2.0));
        expectedPoints_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 1.0,  0.0,-2.0));
        expectedIntsec_2_1.push_back(Vec3(-1.0,  0.0,-2.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_EDGE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, -1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_EDGE_EDGE PARALLEL 6
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test64) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0)); 
        polygon_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedPoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 0.0,  0.0,  2.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0)); 
        polygon_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 0.0,  0.0,-2.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_EDGE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, -1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_EDGE_EDGE PARALLEL 7
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test65) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0)); 
        polygon_1.push_back(Vec3( 1.0,  0.0,  2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3( 1.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedPoints_1.push_back(Vec3( 1.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 1.0,  0.0,  2.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-1.0,  0.0, -2.0)); 
        polygon_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3(-1.0,  0.0, -2.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3(-1.0,  0.0, -2.0));
        expectedPoints_2.push_back(Vec3(-2.0,  0.0, -2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3(-1.0,  0.0,-2.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_EDGE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_EDGE, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, -1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_EDGE_VERTEX
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test66) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0)); 
        polygon_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 2.0,  0.0,  2.0));
        expectedPoints_1.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 0.0,  0.0,  2.0));

        Vec3       anchor_1(0.0,  0.0,  0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0)); 

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 0.0,  0.0,-2.0));

        Vec3       anchor_2(0.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_EDGE, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_VERTEX, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, -1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_NORMAL_2);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_EDGE_NORMAL_1);
        }
        //showResult(body1, body2, info);
    }
}


/**  @brief process_VERTEX_VERTEX
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test67) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 0.0,  0.0,  2.0)); 

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1;
        expectedPoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedIntsec_1_1;
        expectedIntsec_1_1.push_back(Vec3( 0.0,  0.0,  2.0));

        Vec3       anchor_1(0.0,  0.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 1.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0)); 
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0)); 

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2;
        expectedPoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedIntsec_2_1;
        expectedIntsec_2_1.push_back(Vec3( 0.0,  0.0,-2.0));

        Vec3       anchor_2(0.0,  0.0, -2.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0,-1.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double rotationRad = 0.0;
        Vec3 tiltAxis(-1.0, 0.0, 0.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactPointsAndNormalGenerator exp1(
          ContactPairInfo::FT_VERTEX, expectedPoints_1, expectedIntsec_1_1,
          ContactPairInfo::FT_VERTEX, expectedPoints_2, expectedIntsec_2_1 );

        generateTestRigidBodies(
            polygon_1,
            anchor_1,
            polygon_2,
            anchor_2,
            tiltAxis,
            rotationRad,
            body1,
            body2,
            YawAngle,
            swap12
        );

        setActiveFeatures(
            activePoints_1,
            activePoints_2,
            body1,
            body2,
            info,
            swap12
        );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
                   EPSILON_SQUARED, EPSILON_SQUARED, true, 6, std::cerr);
        //gen_01.setLogLevel(ContactPointsAndNormalGenerator::ALL);
        gen_01.generateContatctPointsAndNormalFromActiveFeatures();

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, -1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_VERTEX_VERTEX_AVG);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, randQ.rotate(Vec3( 0.0,  0.0, 1.0)));
            EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_VERTEX_VERTEX_AVG);
        }
        //showResult(body1, body2, info);
    }
}



/**  @brief reduceContactPoints()
 */
TEST_F(ContactPointsAndNormalGeneratorTest, Test68) {

    for (long i = 0; i < 20; i++) {

        long numPoints = 4*(i+1);

        vector<IntersectionFinderConvexPolygon2D::OutputElem> points_01;


        double randomRad = 2.0*M_PI*(rand100()/100.0) - M_PI;

        for (long j = 0; j < numPoints; j++) {

            double x =       cos( (double(j)/double(numPoints))*2.0*M_PI);
            double y = 4.0 * sin( (double(j)/double(numPoints))*2.0*M_PI);

            IntersectionFinderConvexPolygon2D::OutputElem 
                elem_01(
                    Vec2(x*cos(randomRad) - y*sin(randomRad),
                         x*sin(randomRad) + y*cos(randomRad) ),
                    0, 0, 0, 0, 
                    IntersectionFinderConvexPolygon2D::IT_NONE
                );
            points_01.push_back(elem_01);
        }

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;
        ContactPointsAndNormalGenerator 
            gen_01(body1, body2, info, info.mLastVertices,
            EPSILON_SQUARED, EPSILON_SQUARED, true, 8, std::cerr
        );

        vector<bool> res_01 = gen_01.reduceContactPoints(points_01);
        EXPECT_EQ(res_01.size(), points_01.size());
        for (long i = 0; i < points_01.size(); i++) {
            if (i%(numPoints/4)==0) {
                EXPECT_EQ(res_01[i], true);
            }
            else if ( (i-1)%(numPoints/2)!=0 && 
                      (i-1)%(numPoints/4)==0   ) {
                EXPECT_EQ(res_01[i], true);
            }
            else if ( (i+1)%(numPoints/2)!=0 && 
                      (i+1)%(numPoints/4)==0   ) {
                EXPECT_EQ(res_01[i], true);
            }
            else {
                EXPECT_EQ(res_01[i], false);
            }
        }
    }
}


} // namespace Makena

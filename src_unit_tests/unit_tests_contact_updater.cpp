#include "gtest/gtest.h"
#include <vector>
#include "primitives.hpp"
#include "manifold.hpp"
#include "quaternion.hpp"
#include "jacobian_constraint.hpp"
#include "contact_updater.hpp"

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


class ContactUpdaterTest : public ::testing::Test {

  protected:  

    ContactUpdaterTest(){;}
    virtual ~ContactUpdaterTest(){;}
    virtual void SetUp() {;};
    virtual void TearDown() {;};

};


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


//static double rand90()
//{
//    return ((rand()%9000)+1)/101.0;
//}


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
        body1.setGeomConfig(CoM1, Qunit, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0));
        body1.updateGeomConfig(0.0);
        body1.commitGeomConfig();

        body2.setObjectAttributes(m2Serial, 1.0, I, 1.0, 1.0, false);
        body2.setGeomConfig(CoM2, tiltQ, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0));
        body2.updateGeomConfig(0.0);
        body2.commitGeomConfig();
    }
    else {
        body2.setObjectAttributes(m1Serial, 1.0, I, 1.0, 1.0, false);
        body2.setGeomConfig(CoM1, Qunit, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0));
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


class ExpStructContactUpdater {

  public:

    ExpStructContactUpdater(
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
}


static bool checkResult(
    ConvexRigidBody&         body1,
    ConvexRigidBody&         body2,
    ContactPairInfo&         info,
    ExpStructContactUpdater& exp,
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


/**  @brief updateContactNormal() FACE-FACE #1
 */
TEST_F(ContactUpdaterTest, Test01) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3( 1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3         normal_1(0.0, 0.0, 1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0, -1.0));
        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3( 1.0, -1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0, -1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0, 0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3         tiltAxis(1.0, 0.0, 0.0);
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

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.updateContactNormal();

        normal_1 = rm.rotate(normal_1);
        normal_2 = rm.rotate(normal_2);
        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, normal_2);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, normal_1);
        }

        EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_FACE_NORMAL_1);
    }
}


/**  @brief updateContactNormal() FACE-FACE #2
 */
TEST_F(ContactUpdaterTest, Test02) {

    for (long i = 0; i < 20; i++) {


        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3( 1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3         normal_1(0.0, 0.0, 1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0, -1.0));
        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3( 1.0, -1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0, -1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0, -1.0*cos(M_PI*0.25), -1.0*sin(M_PI*0.25));

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3         tiltAxis(-1.0, 0.0, 0.0);
        double       rotationRad = M_PI * 0.25;

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

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, body2);

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.updateContactNormal();

        normal_1 = rm.rotate(normal_1);
        normal_2 = rm.rotate(normal_2);

        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, normal_2);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, normal_1);
        }
        EXPECT_EQ(info.mContactNormalType, ContactPairInfo::NT_FACE_NORMAL_1);
    }
}


/**  @brief updateContactNormal() FACE-EDGE #1
 */
TEST_F(ContactUpdaterTest, Test03) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3( 1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3         normal_1(0.0, 0.0, 1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  0.0, -1.0));
        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  0.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3         tiltAxis(-1.0, 0.0, 0.0);
        double       rotationRad = M_PI * 0.25;

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

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, body2);

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.updateContactNormal();

        normal_1 = rm.rotate(normal_1);
        normal_2 = rm.rotate(normal_2);

        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, normal_1*-1.0);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_FACE_NORMAL_2);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, normal_1);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_FACE_NORMAL_1);
        }
    }
}


/**  @brief updateContactNormal() FACE-EDGE #2
 */
TEST_F(ContactUpdaterTest, Test04) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3( 1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3         normal_1(0.0, 0.0, 1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  0.0, -1.0));
        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  0.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3         tiltAxis(-1.0, 1.0, 0.0);
        double       rotationRad = M_PI * 0.25 * (2.0*rand100()/100.0 - 1.0);

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

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, body2);

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.updateContactNormal();

        normal_1 = rm.rotate(normal_1);
        normal_2 = rm.rotate(normal_2);

        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, normal_1*-1.0);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_FACE_NORMAL_2);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, normal_1);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_FACE_NORMAL_1);
        }
    }
}


/**  @brief updateContactNormal() FACE-VERTEX 
 */
TEST_F(ContactUpdaterTest, Test05) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3( 1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3         normal_1(0.0, 0.0, 1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3         tiltAxis(-1.0, 1.0, 0.0);
        double       rotationRad = M_PI * 0.25 * (2.0*rand100()/100.0 - 1.0);

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

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, body2);

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.updateContactNormal();

        normal_1 = rm.rotate(normal_1);
        normal_2 = rm.rotate(normal_2);

        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, normal_1*-1.0);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_FACE_NORMAL_2);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, normal_1);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_FACE_NORMAL_1);
        }
    }
}


/**  @brief updateContactNormal() EDGE-EDGE parallel
 */
TEST_F(ContactUpdaterTest, Test06) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  0.0,  1.0));
        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  0.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3         normal_1(0.0, 0.0, 1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  0.0, -1.0));
        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  0.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3         tiltAxis( 0.0, 1.0, 0.0);
        double       rotationRad = M_PI * 0.25 * (2.0*rand100()/100.0 - 1.0);

        Vec3         normal_2(-1.0*sin(rotationRad),  
                               0.0, 
                              -1.0*cos(rotationRad));

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

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, body2);

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);
        normal_1 = rm.rotate(normal_1);
        normal_2 = rm.rotate(normal_2);

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.updateContactNormal();

        Vec3 normal_01 = normal_1 - normal_2;
        normal_01.normalize();

        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, normal_01*-1.0);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, normal_01);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
    }
}


/**  @brief updateContactNormal() EDGE-EDGE parallel
 */
TEST_F(ContactUpdaterTest, Test07) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  0.0,  1.0));
        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  0.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3         normal_1(0.0, 0.0, 1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  0.0, -1.0));
        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  0.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3         tiltAxis(-1.0, 0.0, 0.0);
        double       rotationRad = M_PI * 0.25 * (2.0*rand100()/100.0 - 1.0);

        Vec3         normal_2( 0.0,
                              -1.0*sin(rotationRad),  
                              -1.0*cos(rotationRad));

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

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, body2);

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);
        normal_1 = rm.rotate(normal_1);
        normal_2 = rm.rotate(normal_2);

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.updateContactNormal();

        Vec3 normal_01 = normal_1 - normal_2;
        normal_01.normalize();

        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, normal_01*-1.0);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, normal_01);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_EDGE_EDGE_AVG);
        }
    }
}


/**  @brief updateContactNormal() EDGE-EDGE crossing
 */
TEST_F(ContactUpdaterTest, Test08) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  0.0,  1.0));
        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  0.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3         normal_1(0.0, 0.0, 1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  0.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  0.0, -1.0));
        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  0.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  0.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3         tiltAxis( 0.0, 0.0, 1.0);
        double       rotationRad = 0.0;
        while (rotationRad <= 0.01) {
            rotationRad = M_PI * 0.5 * (2.0*rand100()/100.0 - 1.0);
        }

        Vec3         normal_2( 0.0, 0.0, -1.0);

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

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, body2);

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);
        normal_1 = rm.rotate(normal_1);
        normal_2 = rm.rotate(normal_2);

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.updateContactNormal();

        Vec3 normal_01 = normal_1 - normal_2;
        normal_01.normalize();

        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, normal_1*-1.0);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_EDGE_CROSS_EDGE);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, normal_1);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_EDGE_CROSS_EDGE);
        }
    }
}


/**  @brief updateContactNormal() EDGE-EDGE crossing
 */
TEST_F(ContactUpdaterTest, Test09) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  0.0,  1.0));
        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  0.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3         normal_1(0.0, 0.0, 1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0, -1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3         tiltAxis(-1.0, 0.0, 0.0);
        double       rotationRad = M_PI * 0.25 * (2.0*rand100()/100.0 - 1.0);

        Vec3         normal_2( 0.0, -1.0*sin(rotationRad), -1.0*cos(rotationRad));

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

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, body2);

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);
        normal_1 = rm.rotate(normal_1);
        normal_2 = rm.rotate(normal_2);

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.updateContactNormal();

        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, normal_2);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_EDGE_CROSS_EDGE);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, normal_2*-1.0);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_EDGE_CROSS_EDGE);
        }
    }
}


/**  @brief updateContactNormal() EDGE-EDGE crossing
 */
TEST_F(ContactUpdaterTest, Test10) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  0.0,  1.0));
        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  0.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3         dir_1(1.0, 0.0, 0.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0, -1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3         tiltAxis( 0.0, 1.0, 0.0);
        double       rotationRad = M_PI * 0.25 * (2.0*rand100()/100.0 - 1.0);

        Vec3         dir_2(1.0, 1.0, 0.0);
        Quaternion q_2(tiltAxis, rotationRad);
        dir_2 = q_2.rotate(dir_2);

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

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, body2);

        auto cr = dir_1.cross(dir_2);
        if (cr.z() < 0.0) {
            cr.scale(-1.0);
        }
        cr.normalize();

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);
        cr = rm.rotate(cr);

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0,  6, std::cerr);
        upd.updateContactNormal();

        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, cr*-1.0);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_EDGE_CROSS_EDGE);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, cr);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_EDGE_CROSS_EDGE);
        }
    }
}


/**  @brief updateContactNormal() EDGE-VERTEX
 */
TEST_F(ContactUpdaterTest, Test11) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  0.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  0.0,  1.0));
        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  0.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  0.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3         normal_1(0.0, 0.0, 1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

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
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3          tiltAxis = randVec3D100();
        tiltAxis.normalize();
        double       rotationRad = M_PI * 0.5 * (2.0*rand100()/100.0 - 1.0);

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

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, body2);

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);
        normal_1 = rm.rotate(normal_1);

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.updateContactNormal();

        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, normal_1*-1.0);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_EDGE_NORMAL_2);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, normal_1);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_EDGE_NORMAL_1);
        }
    }
}


/**  @brief updateContactNormal() VERTEX-VERTEX
 */
TEST_F(ContactUpdaterTest, Test12) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 0.0,  0.0,  1.0));
        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3         normal_1(0.0, 0.0, 1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0, 0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3          tiltAxis = randVec3D100();
        tiltAxis.normalize();
        double       rotationRad = M_PI * 0.5 * (2.0*rand100()/100.0 - 1.0);

        Quaternion q_01(tiltAxis, rotationRad);
        normal_2 = q_01.rotate(normal_2);

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

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, body2);

        auto rm = randomRotQ(); 
        auto rd = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, rm, rd);
        normal_1 = rm.rotate(normal_1);
        normal_2 = rm.rotate(normal_2);

        auto normal_01 = normal_1 - normal_2;
        normal_01.normalize();

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.updateContactNormal();

        if (swap12) {
            EXPECT_EQ(info.mContactNormal1To2, normal_01*-1.0);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_VERTEX_VERTEX_AVG);
        }
        else {
            EXPECT_EQ(info.mContactNormal1To2, normal_01);
            EXPECT_EQ(info.mContactNormalType, 
                      ContactPairInfo::NT_VERTEX_VERTEX_AVG);
        }
    }
}


/**  @brief makeInitialBDVertices() FACE-FACE;
 */
TEST_F(ContactUpdaterTest, Test13) {

    for (long i = 0; i < 20; i++) {

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3( 1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0, -1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0, -1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3( 1.0, -1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0, -1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0, 0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3          tiltAxis(1.0, 0.0, 0.0);
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

        auto vit1_01 = findVertex(Vec3( 1.0, 1.0,  1.0),body1);
        auto vit1_02 = findVertex(Vec3(-1.0, 1.0,  1.0),body1);
        auto vit1_03 = findVertex(Vec3( 1.0,-1.0,  1.0),body1);
        auto vit1_04 = findVertex(Vec3(-1.0,-1.0,  1.0),body1);

        auto vit2_01 = findVertex(Vec3( 1.0, 1.0, -1.0),body2);
        auto vit2_02 = findVertex(Vec3(-1.0, 1.0, -1.0),body2);
        auto vit2_03 = findVertex(Vec3( 1.0,-1.0, -1.0),body2);
        auto vit2_04 = findVertex(Vec3(-1.0,-1.0, -1.0),body2);


        Vec3 p1_01    = (*vit1_01)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p1_02    = (*vit1_02)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p1_03    = (*vit1_03)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p1_04    = (*vit1_04)->pGCS(body1.Qmat(), body1.CoM());

        Vec3 p2_01    = (*vit2_01)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_02    = (*vit2_02)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_03    = (*vit2_03)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_04    = (*vit2_04)->pGCS(body2.Qmat(), body2.CoM());

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.makeInitialBDVertices();
        auto& res_01 = info.mLastVertices[0];

        EXPECT_EQ(res_01.v1()==vit1_01||res_01.v1()==vit1_02||res_01.v1()==vit1_03||res_01.v1()==vit1_04, true);
        EXPECT_EQ(res_01.v2()==vit2_01||res_01.v2()==vit2_02||res_01.v2()==vit2_03||res_01.v2()==vit2_04, true);

        EXPECT_EQ(res_01.p()==(p1_01 - p2_01)||
                  res_01.p()==(p1_02 - p2_01)||
                  res_01.p()==(p1_03 - p2_01)||
                  res_01.p()==(p1_04 - p2_01)||
                  res_01.p()==(p1_01 - p2_02)||
                  res_01.p()==(p1_02 - p2_02)||
                  res_01.p()==(p1_03 - p2_02)||
                  res_01.p()==(p1_04 - p2_02)||
                  res_01.p()==(p1_01 - p2_03)||
                  res_01.p()==(p1_02 - p2_03)||
                  res_01.p()==(p1_03 - p2_03)||
                  res_01.p()==(p1_04 - p2_03)||
                  res_01.p()==(p1_01 - p2_04)||
                  res_01.p()==(p1_02 - p2_04)||
                  res_01.p()==(p1_03 - p2_04)||
                  res_01.p()==(p1_04 - p2_04), true);
    }
}


/**  @brief makeInitialBDVertices() FACE-EDGE;
 */
TEST_F(ContactUpdaterTest, Test14) {

    for (long i = 0; i < 20; i++) {

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3( 1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0, -1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0, -1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0, 0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3          tiltAxis(1.0, 0.0, 0.0);
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

        auto vit1_01 = findVertex(Vec3( 1.0, 1.0,  1.0),body1);
        auto vit1_02 = findVertex(Vec3(-1.0, 1.0,  1.0),body1);
        auto vit1_03 = findVertex(Vec3( 1.0,-1.0,  1.0),body1);
        auto vit1_04 = findVertex(Vec3(-1.0,-1.0,  1.0),body1);

        auto vit2_01 = findVertex(Vec3( 1.0, 1.0, -1.0),body2);
        auto vit2_02 = findVertex(Vec3(-1.0, 1.0, -1.0),body2);

        Vec3 p1_01    = (*vit1_01)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p1_02    = (*vit1_02)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p1_03    = (*vit1_03)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p1_04    = (*vit1_04)->pGCS(body1.Qmat(), body1.CoM());

        Vec3 p2_01    = (*vit2_01)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_02    = (*vit2_02)->pGCS(body2.Qmat(), body2.CoM());

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.makeInitialBDVertices();
        auto& res_01 = info.mLastVertices[0];

        EXPECT_EQ(res_01.v1()==vit1_01||res_01.v1()==vit1_02||res_01.v1()==vit1_03||res_01.v1()==vit1_04, true);
        EXPECT_EQ(res_01.v2()==vit2_01||res_01.v2()==vit2_02, true);
        EXPECT_EQ(res_01.p()==(p1_01 - p2_01)||
                  res_01.p()==(p1_02 - p2_01)||
                  res_01.p()==(p1_03 - p2_01)||
                  res_01.p()==(p1_04 - p2_01)||
                  res_01.p()==(p1_01 - p2_02)||
                  res_01.p()==(p1_02 - p2_02)||
                  res_01.p()==(p1_03 - p2_02)||
                  res_01.p()==(p1_04 - p2_02), true);
        

    }
}


/**  @brief makeInitialBDVertices() EDGE-FACE;
 */
TEST_F(ContactUpdaterTest, Test15) {

    for (long i = 0; i < 20; i++) {

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0, -1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3( 1.0, -1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0, -1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0, 0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3          tiltAxis(1.0, 0.0, 0.0);
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

        auto vit1_01 = findVertex(Vec3( 1.0, 1.0,  1.0),body1);
        auto vit1_02 = findVertex(Vec3(-1.0, 1.0,  1.0),body1);
        auto vit2_01 = findVertex(Vec3( 1.0, 1.0, -1.0),body2);
        auto vit2_02 = findVertex(Vec3(-1.0, 1.0, -1.0),body2);
        auto vit2_03 = findVertex(Vec3( 1.0,-1.0, -1.0),body2);
        auto vit2_04 = findVertex(Vec3(-1.0,-1.0, -1.0),body2);

        Vec3 p1_01    = (*vit1_01)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p1_02    = (*vit1_02)->pGCS(body1.Qmat(), body1.CoM());

        Vec3 p2_01    = (*vit2_01)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_02    = (*vit2_02)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_03    = (*vit2_03)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_04    = (*vit2_04)->pGCS(body2.Qmat(), body2.CoM());

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.makeInitialBDVertices();
        auto& res_01 = info.mLastVertices[0];

        EXPECT_EQ(res_01.v1()==vit1_01||res_01.v1()==vit1_02, true);
        EXPECT_EQ(res_01.v2()==vit2_01||res_01.v2()==vit2_02||res_01.v2()==vit2_03||res_01.v2()==vit2_04, true);
        EXPECT_EQ(res_01.p()==(p1_01 - p2_01)||
                  res_01.p()==(p1_01 - p2_02)||
                  res_01.p()==(p1_01 - p2_03)||
                  res_01.p()==(p1_01 - p2_04)||
                  res_01.p()==(p1_02 - p2_01)||
                  res_01.p()==(p1_02 - p2_02)||
                  res_01.p()==(p1_02 - p2_03)||
                  res_01.p()==(p1_02 - p2_04), true);
        
    }
}


/**  @brief makeInitialBDVertices() FACE-VERTEX;
 */
TEST_F(ContactUpdaterTest, Test16) {

    for (long i = 0; i < 20; i++) {

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3( 1.0, -1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0, -1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0, -1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0, 0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3          tiltAxis(1.0, 0.0, 0.0);
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

        auto vit1_01 = findVertex(Vec3( 1.0, 1.0,  1.0),body1);
        auto vit1_02 = findVertex(Vec3(-1.0, 1.0,  1.0),body1);
        auto vit1_03 = findVertex(Vec3( 1.0,-1.0,  1.0),body1);
        auto vit1_04 = findVertex(Vec3(-1.0,-1.0,  1.0),body1);
        auto vit2_01 = findVertex(Vec3( 1.0, 1.0, -1.0),body2);

        Vec3 p1_01    = (*vit1_01)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p1_02    = (*vit1_02)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p1_03    = (*vit1_03)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p1_04    = (*vit1_04)->pGCS(body1.Qmat(), body1.CoM());

        Vec3 p2_01    = (*vit2_01)->pGCS(body2.Qmat(), body2.CoM());

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.makeInitialBDVertices();
        auto& res_01 = info.mLastVertices[0];

        EXPECT_EQ(res_01.v1()==vit1_01||res_01.v1()==vit1_02||res_01.v1()==vit1_03||res_01.v1()==vit1_04, true);
        EXPECT_EQ(res_01.v2()==vit2_01, true);
        EXPECT_EQ(res_01.p()==(p1_01 - p2_01)||
                  res_01.p()==(p1_02 - p2_01)||
                  res_01.p()==(p1_03 - p2_01)||
                  res_01.p()==(p1_04 - p2_01)  , true);
        
    }
}


/**  @brief makeInitialBDVertices() VERTEX-FACE;
 */
TEST_F(ContactUpdaterTest, Test17) {

    for (long i = 0; i < 20; i++) {

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0, -1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3( 1.0, -1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0, -1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0, 0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3          tiltAxis(1.0, 0.0, 0.0);
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

        auto vit1_01 = findVertex(Vec3( 1.0, 1.0,  1.0),body1);

        auto vit2_01 = findVertex(Vec3( 1.0, 1.0, -1.0),body2);
        auto vit2_02 = findVertex(Vec3(-1.0, 1.0, -1.0),body2);
        auto vit2_03 = findVertex(Vec3( 1.0,-1.0, -1.0),body2);
        auto vit2_04 = findVertex(Vec3(-1.0,-1.0, -1.0),body2);

        Vec3 p1_01    = (*vit1_01)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p2_01    = (*vit2_01)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_02    = (*vit2_02)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_03    = (*vit2_03)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_04    = (*vit2_04)->pGCS(body2.Qmat(), body2.CoM());

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0,  6, std::cerr);
        upd.makeInitialBDVertices();
        auto& res_01 = info.mLastVertices[0];

        EXPECT_EQ(res_01.v1()==vit1_01, true);
        EXPECT_EQ(res_01.v2()==vit2_01||res_01.v2()==vit2_02||res_01.v2()==vit2_03||res_01.v2()==vit2_04, true);

        EXPECT_EQ(res_01.p()==(p1_01 - p2_01)||
                  res_01.p()==(p1_01 - p2_02)||
                  res_01.p()==(p1_01 - p2_03)||
                  res_01.p()==(p1_01 - p2_04)  , true);
        
    }
}


/**  @brief makeInitialBDVertices() EDGE-EDGE
 */
TEST_F(ContactUpdaterTest, Test18) {

    for (long i = 0; i < 20; i++) {

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0, -1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0, 0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3          tiltAxis(1.0, 0.0, 0.0);
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

        auto vit1_01 = findVertex(Vec3( 1.0, 1.0,  1.0),body1);
        auto vit1_02 = findVertex(Vec3(-1.0, 1.0,  1.0),body1);
        auto vit2_01 = findVertex(Vec3( 1.0, 1.0, -1.0),body2);
        auto vit2_02 = findVertex(Vec3(-1.0, 1.0, -1.0),body2);

        Vec3 p1_01    = (*vit1_01)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p1_02    = (*vit1_02)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_01    = (*vit2_01)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_02    = (*vit2_02)->pGCS(body2.Qmat(), body2.CoM());

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.makeInitialBDVertices();
        auto& res_01 = info.mLastVertices[0];

        EXPECT_EQ(res_01.v1()==vit1_01||res_01.v1()==vit1_02, true);
        EXPECT_EQ(res_01.v2()==vit2_01||res_01.v2()==vit2_02, true);
        EXPECT_EQ(res_01.p()==(p1_01 - p2_01)||
                  res_01.p()==(p1_01 - p2_02)||
                  res_01.p()==(p1_02 - p2_01)||
                  res_01.p()==(p1_02 - p2_02)  , true);



    }
}


/**  @brief makeInitialBDVertices() EDGE-VERTEX
 */
TEST_F(ContactUpdaterTest, Test19) {

    for (long i = 0; i < 20; i++) {

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));
        activePoints_1.push_back(Vec3(-1.0,  1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0, -1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0, 0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3          tiltAxis(1.0, 0.0, 0.0);
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

        auto vit1_01 = findVertex(Vec3( 1.0, 1.0,  1.0),body1);
        auto vit1_02 = findVertex(Vec3(-1.0, 1.0,  1.0),body1);
        auto vit2_01 = findVertex(Vec3( 1.0, 1.0, -1.0),body2);

        Vec3 p1_01    = (*vit1_01)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p1_02    = (*vit1_02)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_01    = (*vit2_01)->pGCS(body2.Qmat(), body2.CoM());

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.makeInitialBDVertices();
        auto& res_01 = info.mLastVertices[0];

        EXPECT_EQ(res_01.v1()==vit1_01||res_01.v1()==vit1_02, true);
        EXPECT_EQ(res_01.v2()==vit2_01, true);
        EXPECT_EQ(res_01.p()==(p1_01 - p2_01)||
                  res_01.p()==(p1_02 - p2_01)  , true);

    }
}


/**  @brief makeInitialBDVertices() VERTEX-EDGE
 */
TEST_F(ContactUpdaterTest, Test20) {

    for (long i = 0; i < 20; i++) {

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0, -1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));
        activePoints_2.push_back(Vec3(-1.0,  1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0, 0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3          tiltAxis(1.0, 0.0, 0.0);
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

        auto vit1_01 = findVertex(Vec3( 1.0, 1.0,  1.0),body1);
        auto vit2_01 = findVertex(Vec3( 1.0, 1.0, -1.0),body2);
        auto vit2_02 = findVertex(Vec3(-1.0, 1.0, -1.0),body2);


        Vec3 p1_01    = (*vit1_01)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p2_01    = (*vit2_01)->pGCS(body2.Qmat(), body2.CoM());
        Vec3 p2_02    = (*vit2_02)->pGCS(body2.Qmat(), body2.CoM());

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.makeInitialBDVertices();
        auto& res_01 = info.mLastVertices[0];

        EXPECT_EQ(res_01.v1()==vit1_01, true);
        EXPECT_EQ(res_01.v2()==vit2_01||res_01.v2()==vit2_02, true);
        EXPECT_EQ(res_01.p()==(p1_01 - p2_01)||
                  res_01.p()==(p1_01 - p2_02)  , true);
    }
}


/**  @brief makeInitialBDVertices() VERTEX-VERTEX
 */
TEST_F(ContactUpdaterTest, Test21) {

    for (long i = 0; i < 20; i++) {

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0, -1.0));
        polygon_1.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_1.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_1.push_back(Vec3(-1.0,  1.0,  1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 1.0,  1.0,  1.0));

        Vec3         anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0, -1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0, -1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0, -1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 1.0,  1.0, -1.0));

        Vec3         anchor_2(0.0,  0.0, -1.0);

        Vec3         normal_2(0.0, 0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);

        Vec3          tiltAxis(1.0, 0.0, 0.0);
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

        auto vit1_01 = findVertex(Vec3( 1.0, 1.0,  1.0),body1);
        auto vit2_01 = findVertex(Vec3( 1.0, 1.0, -1.0),body2);
        Vec3 p1_01    = (*vit1_01)->pGCS(body1.Qmat(), body1.CoM());
        Vec3 p2_01    = (*vit2_01)->pGCS(body2.Qmat(), body2.CoM());

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);
        upd.makeInitialBDVertices();
        auto& res_01 = info.mLastVertices[0];

        EXPECT_EQ(res_01.v1()==vit1_01, true);
        EXPECT_EQ(res_01.v2()==vit2_01, true);
        EXPECT_EQ(res_01.p(), p1_01 - p2_01);

    }
}


/**  @brief decomposeQ()
 */
TEST_F(ContactUpdaterTest, Test22) {

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

    ContactUpdater upd(
          body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                              EPSILON_SQUARED, 1.0, 6,  std::cerr);

    Q.push_back(bd11);
    upd.decomposeQ(Q, vertices1, vertices2);
    EXPECT_EQ(vertices1.size(), 1);
    EXPECT_EQ(vertices2.size(), 1);
    EXPECT_EQ(vertices1[0], vit1_01);
    EXPECT_EQ(vertices2[0], vit2_01);

    Q.clear();
    vertices1.clear();
    vertices2.clear();
    Q.push_back(bd11);
    Q.push_back(bd12);
    upd.decomposeQ(Q, vertices1, vertices2);
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
    upd.decomposeQ(Q, vertices1, vertices2);
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
    upd.decomposeQ(Q, vertices1, vertices2);
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
    upd.decomposeQ(Q, vertices1, vertices2);
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
    upd.decomposeQ(Q, vertices1, vertices2);
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

    upd.decomposeQ(Q, vertices1, vertices2);
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

    upd.decomposeQ(Q, vertices1, vertices2);
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



/**  @brief  findRelativeVelocityOfBody1RelativeToBody2();
 */
TEST_F(ContactUpdaterTest, Test23) {

    for (long i = 0; i < 20; i++) {

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

        ContactUpdater upd(
              body1, body2, info, EPSILON_SQUARED, EPSILON_SQUARED, 
                                  EPSILON_SQUARED, 1.0, 6,  std::cerr);

        auto v01 = upd.findRelativeVelocityOfBody1RelativeToBody2();
        EXPECT_EQ(v01, exp_01);
    }
}


}// Makena

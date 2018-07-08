#include "gtest/gtest.h"
#include <vector>
#include "primitives.hpp"
#include "manifold.hpp"
#include "quaternion.hpp"
#include "jacobian_constraint.hpp"
#include "contact_updater_face_edge.hpp"

#include <iostream>
#include <fstream>
#include <sstream>


namespace Makena {


class ContactUpdater_FACE_EDGETest : public ::testing::Test {

  protected:  

    ContactUpdater_FACE_EDGETest(){;}
    virtual ~ContactUpdater_FACE_EDGETest(){;}
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

/*
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

/*
static double rand90()
{
    return ((rand()%9000)+1)/101.0;
}
*/
static Vec3 randVec3D100()
{
    return Vec3(rand100(), rand100(), rand100());
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


static bool isAnyOf(
    enum predicate pTest,
    enum predicate p1,
    enum predicate p2,
    enum predicate p3
) {
    return pTest==p1 || pTest==p2 || pTest==p3;
}


static bool isAnyOf(
    enum predicate pTest,
    enum predicate p1,
    enum predicate p2
) {
    return pTest==p1 || pTest==p2;
}


static bool isAnyOf(
    enum predicate pTest,
    enum predicate p1
) {
    return pTest==p1;
}
*/

class ExpStructFaceEdge {
  public:

    ExpStructFaceEdge(
        enum ContactPairInfo::FeatureType t1, 
        double p11x,
        double p11y,
        double p12x,
        double p12y,
        enum ContactPairInfo::FeatureType t2,
        double p21x,
        double p21y,
        double p22x,
        double p22y
    ):  mType1(t1),
        p11(p11x, p11y,  0.5),
        p12(p12x, p12y,  0.5),
        mType2(t2),
        p21(p21x, p21y,  -0.5),
        p22(p22x, p22y,  -0.5){;}
  
    enum ContactPairInfo::FeatureType mType1;
    Vec3 p11;
    Vec3 p12;
    enum ContactPairInfo::FeatureType mType2;
    Vec3 p21;
    Vec3 p22;
};


static void generateTestPattern(
    vector<Vec2>&    polygon2D_1,
    const Vec2&      anchor2D_1,
    vector<Vec2>&    polygon2D_2,
    const Vec2&      anchor2D_2,
    const Vec2&      tiltAxisXY,
    const double&    rotationRad,
    ConvexRigidBody& body1,
    ConvexRigidBody& body2,
    ContactPairInfo& info,
    bool             edgeOn2
) {
    // CoM is (0.0, 0.0, -0.5)
    vector<Vec3> polygon3D_1;
    for (auto v2 : polygon2D_1) {

        polygon3D_1.push_back(Vec3(v2.x(), v2.y(), 0.5));
    }
    polygon3D_1.push_back(Vec3(0.0, 0.0, -0.5));

    Manifold m1;
    enum predicate pred;
    m1.findConvexHull(polygon3D_1, pred);
    auto m1Serial = m1.exportData();
    Vec3 CoM1(0.0, 0.0, -0.5);

    // CoM is (0.0, 0.0, +0.5)
    vector<Vec3> polygon3D_2;
    for (auto v2 : polygon2D_2) {
        polygon3D_2.push_back(Vec3(v2.x(), v2.y(),-0.5));
    }
    polygon3D_2.push_back(Vec3(-1.0, -1.0, 0.5));
    polygon3D_2.push_back(Vec3(-1.0,  1.0, 0.5));
    polygon3D_2.push_back(Vec3( 1.0, -1.0, 0.5));
    polygon3D_2.push_back(Vec3( 1.0,  1.0, 0.5));

    Manifold m2;
    m2.findConvexHull(polygon3D_2, pred);
    auto m2Serial = m2.exportData();
    Vec3 CoM2(0.0, 0.0, 0.5);

    Vec3 anchor3D_1(anchor2D_1.x(), anchor2D_1.y(),  0.5);
    Vec3 anchor3D_2(anchor2D_2.x(), anchor2D_2.y(), -0.5);

    // Offsets to bring the anchor points to the origin.
    CoM1 = anchor3D_1 * -1.0;
    CoM2 = anchor3D_2 * -1.0;

    // Tilt body 2 around origin.
    Vec3       tilt3D(tiltAxisXY.x(), tiltAxisXY.y(), 0.0);
    tilt3D.normalize();
    Quaternion tiltQ (tilt3D, rotationRad);
    CoM2 = tiltQ.rotate(CoM2);
    Mat3x3 I ( 1.0, 0.0, 0.0, 
               0.0, 1.0, 0.0, 
               0.0, 0.0, 1.0  );

    Quaternion Qunit(1.0, 0.0, 0.0, 0.0);

    if (edgeOn2) {

        body1.setObjectAttributes(m1Serial, 1.0, I, 1.0, 1.0, false);
        body1.setGeomConfig(
                        CoM1, Qunit, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0));
        body1.updateGeomConfig(0.0);
        body1.commitGeomConfig();

        body2.setObjectAttributes(m2Serial, 1.0, I, 1.0, 1.0, false);
        tilt3D.scale(1.0/10.0);
        body2.setGeomConfig(CoM2, tiltQ, Vec3(0.0, 0.0, 0.0), tilt3D);
        body2.updateGeomConfig(0.0);
        body2.commitGeomConfig();

        info.mType1 = ContactPairInfo::FT_FACE;
        auto& ch1 = body1.ConvexHull();
        for (auto fit = ch1.faces().first; fit != ch1.faces().second; fit++) {
            if ((*fit)->nLCS() == Vec3(0.0, 0.0, 1.0)) {
                info.mFit1  = fit;
            }
        }

        info.mType2 = ContactPairInfo::FT_EDGE;
        auto& ch2 = body2.ConvexHull();
        for (auto eit = ch2.edges().first; eit != ch2.edges().second; eit++) {
            auto heit = (*eit)->he1();
            auto vit1 = (*heit)->src();
            auto vit2 = (*heit)->dst();
            auto p1   = (*vit1)->pLCS();
            auto p2   = (*vit2)->pLCS();
            if ( (p1 == polygon3D_2[0] && p2 == polygon3D_2[1])||
                 (p1 == polygon3D_2[1] && p2 == polygon3D_2[0])  ) {
                info.mEit2  = eit;
            }
        }
        info.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_1;
        info.mContactNormal1To2 = Vec3(0.0, 0.0, 1.0);

    }
    else {

        body2.setObjectAttributes(m1Serial, 1.0, I, 1.0, 1.0, false);
        body2.setGeomConfig(
                        CoM1, Qunit, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0));
        body2.updateGeomConfig(0.0);
        body2.commitGeomConfig();

        body1.setObjectAttributes(m2Serial, 1.0, I, 1.0, 1.0, false);
        tilt3D.scale(1.0/10.0);
        body1.setGeomConfig(CoM2, tiltQ, Vec3(0.0, 0.0, 0.0), tilt3D);
        body1.updateGeomConfig(0.0);
        body1.commitGeomConfig();

        info.mType2 = ContactPairInfo::FT_FACE;
        auto& ch2 = body2.ConvexHull();
        for (auto fit = ch2.faces().first; fit != ch2.faces().second; fit++) {
            if ((*fit)->nLCS() == Vec3(0.0, 0.0, 1.0)) {
                info.mFit2  = fit;
            }
        }

        info.mType1 = ContactPairInfo::FT_EDGE;
        auto& ch1 = body1.ConvexHull();
        for (auto eit = ch1.edges().first; eit != ch1.edges().second; eit++) {
            auto heit = (*eit)->he1();
            auto vit1 = (*heit)->src();
            auto vit2 = (*heit)->dst();
            auto p1   = (*vit1)->pLCS();
            auto p2   = (*vit2)->pLCS();
            if ( (p1 == polygon3D_2[0] && p2 == polygon3D_2[1])||
                 (p1 == polygon3D_2[1] && p2 == polygon3D_2[0])  ) {
                info.mEit1  = eit;
            }
        }
        info.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;
        info.mContactNormal1To2 = Vec3(0.0, 0.0,1.0);

    }
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
    auto V1cur = body1.Vang();
    auto Q2cur = body2.Q();
    auto C2cur = body2.CoM();
    auto V2cur = body2.Vang();
    auto Q1new = q * Q1cur;
    auto Q2new = q * Q2cur;
    auto C1new = q.rotate(C1cur + v);
    auto C2new = q.rotate(C2cur + v);
    auto V1new = q.rotate(V1cur);
    auto V2new = q.rotate(V2cur);

/*
    cerr << "Q1new: " << Q1new.s() << "," <<
                         Q1new.x() << "," <<
                         Q1new.y() << "," <<
                         Q1new.z() << "\n";

    cerr << "Q2new: " << Q2new.s() << "," <<
                         Q2new.x() << "," <<
                         Q2new.y() << "," <<
                         Q2new.z() << "\n";
*/


    body1.setGeomConfig(C1new, Q1new, Vec3(0.0, 0.0, 0.0), V1new);
    body1.updateGeomConfig(0.0);
    body1.commitGeomConfig();

    body2.setGeomConfig(C2new, Q2new, Vec3(0.0, 0.0, 0.0), V2new);
    body2.updateGeomConfig(0.0);
    body2.commitGeomConfig();

    info.mContactNormal1To2 = q.rotate(info.mContactNormal1To2);
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
        cerr << (*(info.mVit1))->pLCS() << "\n";
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
    ConvexRigidBody& body1,
    ConvexRigidBody& body2,
    ContactPairInfo& info,
    ExpStructFaceEdge&       exp,
    bool             edgeOn2
) {
    if (edgeOn2) {
        if (info.mType1 != exp.mType1 ||info.mType2 != exp.mType2) {
            cerr << "return1 : " << info.mType1 << "," << info.mType2 << "\n";
            return false;        
        }
    }
    else {
        if (info.mType1 != exp.mType2 ||info.mType2 != exp.mType1) {
            cerr << "return2 : " << info.mType1 << "," << info.mType2 << "\n";
            return false;        
        }
    }

    auto& ch1 = body1.ConvexHull();
    if (info.mType1 == ContactPairInfo::FT_VERTEX) {
        bool found = false;
        for (auto vit =  ch1.vertices().first; 
                  vit != ch1.vertices().second; vit++) {

            auto p = (*vit)->pLCS();
            if ( ((p == exp.p11 && edgeOn2)||(p == exp.p21 && !edgeOn2))
                  && info.mVit1 == vit ) {
                found = true;
                break;
            }
        }
        if (!found) {
            return false;
        }
    }
    else if (info.mType1 == ContactPairInfo::FT_EDGE) {
        bool found = false;
        for (auto eit = ch1.edges().first; eit != ch1.edges().second; eit++) {
            auto heit = (*eit)->he1();
            auto v1   = (*heit)->src();
            auto v2   = (*heit)->dst();
            auto p1 = (*v1)->pLCS();
            auto p2 = (*v2)->pLCS();
            if (edgeOn2) {
                if ( (p1 == exp.p11 && p2 == exp.p12) ||
                     (p1 == exp.p12 && p2 == exp.p11)    ) {
                    if (info.mEit1 == eit) {
                        found = true;
                        break;
                    }
                }
            }
            else {
                if ( (p1 == exp.p21 && p2 == exp.p22) ||
                     (p1 == exp.p22 && p2 == exp.p21)    ) {
                    if (info.mEit1 == eit) {
                        found = true;
                        break;
                    }
                }
            }
        }
        if (!found) {
            return false;
        }
    }
    else if (info.mType1 == ContactPairInfo::FT_FACE) {
        bool found = false;
        for (auto fit = ch1.faces().first; fit != ch1.faces().second; fit++) {
            auto n = (*fit)->nLCS();
            if (n == Vec3(0.0, 0.0, 1.0)) {
                if (info.mFit1 == fit) {
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            return false;
        }
    }

    auto& ch2 = body2.ConvexHull();
    if (info.mType2 == ContactPairInfo::FT_VERTEX) {
        bool found = false;
        for (auto vit =  ch2.vertices().first; 
                  vit != ch2.vertices().second; vit++) {

            auto p = (*vit)->pLCS();
            if (edgeOn2) {
                if (p == exp.p21 && info.mVit2 == vit) {
                    found = true;
                    break;
                }
            }
            else {
                if (p == exp.p11 && info.mVit2 == vit) {
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            return false;
        }
    }
    else if (info.mType2 == ContactPairInfo::FT_EDGE) {
        bool found = false;
        for (auto eit = ch2.edges().first; eit != ch2.edges().second; eit++) {
            auto heit = (*eit)->he1();
            auto v1   = (*heit)->src();
            auto v2   = (*heit)->dst();
            auto p1 = (*v1)->pLCS();
            auto p2 = (*v2)->pLCS();
            if (edgeOn2) {
                if ( (p1 == exp.p21 && p2 == exp.p22) ||
                     (p1 == exp.p22 && p2 == exp.p21)    ) {
                    if (info.mEit2 == eit) {
                        found = true;
                        break;
                    }
                }
            }
            else {
                if ( (p1 == exp.p11 && p2 == exp.p12) ||
                     (p1 == exp.p12 && p2 == exp.p11)    ) {
                    if (info.mEit2 == eit) {
                        found = true;
                        break;
                    }
                }
            }
        }
        if (!found) {
            return false;
        }
    }
    else if (info.mType2 == ContactPairInfo::FT_FACE) {

        bool found = false;
        for (auto fit = ch2.faces().first; fit != ch2.faces().second; fit++) {
            auto n = (*fit)->nLCS();
            if (n == Vec3(0.0, 0.0,  1.0)) {
                if (info.mFit2 == fit) {
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            return false;
        }
    }
    return true;
}



/**  @brief update() finds intersection and the faces are not tilted.
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test01) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0, 1.0));
        polygon2D_1.push_back(Vec2( 1.0,-1.0));
        polygon2D_1.push_back(Vec2(-1.0,-1.0));
        polygon2D_1.push_back(Vec2(-1.0, 1.0));
        Vec2 anchor2D_1( 0.0, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2( 0.5, 0.0));
        polygon2D_2.push_back(Vec2(-0.5, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0, -1.0);
//        double          rotationRad = (M_PI /180.0) * 0.000000000000000001;
        double          rotationRad = 0.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_FACE,  0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_EDGE,  0.5,  0.0, -0.5,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief update() finds intersection and the faces are not tilted.
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test02) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( 0.0, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(1.0, 0.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_FACE,  0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_EDGE,  1.0,  0.0, -1.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);
    }
}


/**  @brief update() finds intersection and the faces are tilted.
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test03) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( 0.0, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0, -1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_FACE,    0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_VERTEX, -1.0,  0.0,  0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);
    }
}


/**  @brief update() finds intersection and the faces are not tilted.
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test04) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -1.5, 0.5);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0, -1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,   -1.0,  1.0, -2.0,  0.0,
            ContactPairInfo::FT_EDGE,   -1.0,  0.0,  1.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);
    }
}


/**  @brief update() finds intersection and the faces are not tilted.
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test05) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -1.5, 0.5);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0, 1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_FACE,    0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_VERTEX,  1.0,  0.0,  0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);
    }
}


/**  @brief update() finds intersection and the faces are not tilted.
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test06) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -1.5, 0.5);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0, 1.0);
        double          rotationRad = (M_PI /180.0) * 0.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_FACE,    0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_EDGE,    1.0,  0.0, -1.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);
    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test07) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1(  0.0, 1.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-0.5, 0.0));
        polygon2D_2.push_back(Vec2( 0.5, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(1.0, 0.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,    1.0,  1.0, -1.0,  1.0,
            ContactPairInfo::FT_EDGE,    0.5,  0.0, -0.5,  0.0 );

        ExpStructFaceEdge exp2(
            ContactPairInfo::FT_FACE,    0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_EDGE,    0.5,  0.0, -0.5,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        if (!res01) { 
            // EXPECT_EQ(res01, false);
            EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2)||
                      checkResult(body1, body2, info, exp2, edgeOn2), true);
        }
    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test08) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1(  0.0, 1.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-0.5, 0.0));
        polygon2D_2.push_back(Vec2( 0.5, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0, -1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,    1.0,  1.0, -1.0,  1.0,
            ContactPairInfo::FT_VERTEX, -0.5,  0.0,  0.0,  0.0 );

        ExpStructFaceEdge exp2(
            ContactPairInfo::FT_FACE,    0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_VERTEX, -0.5,  0.0,  0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();

        //EXPECT_EQ(res01, false);
        if (res01==false) {
            EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2)||
                      checkResult(body1, body2, info, exp2, edgeOn2), true);
        }
    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test09) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1(  0.0, 1.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-0.5, 0.0));
        polygon2D_2.push_back(Vec2( 0.5, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0,  1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,    1.0,  1.0, -1.0,  1.0,
            ContactPairInfo::FT_VERTEX,  0.5,  0.0,  0.0,  0.0 );

        ExpStructFaceEdge exp2(
            ContactPairInfo::FT_FACE,    0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_VERTEX,  0.5,  0.0,  0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2)||
                  checkResult(body1, body2, info, exp2, edgeOn2), true);
    }
}



/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test10) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1(  0.0, 1.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(1.0,  0.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,    1.0,  1.0, -1.0,  1.0,
            ContactPairInfo::FT_EDGE,    1.0,  0.0, -1.0,  0.0 );

        ExpStructFaceEdge exp2(
            ContactPairInfo::FT_FACE,    0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_EDGE,    1.0,  0.0, -1.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2)||
                  checkResult(body1, body2, info, exp2, edgeOn2), true);
    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test11) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1(  0.0, 1.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0, -1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,    1.0,  1.0, -1.0,  1.0,
            ContactPairInfo::FT_VERTEX, -1.0,  0.0,  0.0,  0.0 );

        ExpStructFaceEdge exp2(
            ContactPairInfo::FT_FACE,    0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_VERTEX, -1.0,  0.0,  0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2)||
                  checkResult(body1, body2, info, exp2, edgeOn2), true);
    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test12) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1(  0.0, 1.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0,  1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,    1.0,  1.0, -1.0,  1.0,
            ContactPairInfo::FT_VERTEX,  1.0,  0.0,  0.0,  0.0 );

        ExpStructFaceEdge exp2(
            ContactPairInfo::FT_FACE,    0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_VERTEX,  1.0,  0.0,  0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2)||
                  checkResult(body1, body2, info, exp2, edgeOn2), true);
    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test13) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1(  1.0, 1.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(1.0,  0.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,  1.0,  1.0, -1.0,  1.0,
            ContactPairInfo::FT_EDGE,  1.0,  0.0, -1.0,  0.0 );

        ExpStructFaceEdge exp2(
            ContactPairInfo::FT_FACE,    0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_EDGE,    1.0,  0.0, -1.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2)||
                  checkResult(body1, body2, info, exp2, edgeOn2), true);
    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test14) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -1.0, 1.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0, -1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_VERTEX, -1.0,  1.0, 0.0,  0.0,
            ContactPairInfo::FT_EDGE,   1.0,  0.0, -1.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        //EXPECT_EQ(res01, false);
        if (res01==false) {
            EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);
        }
    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test15) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -1.0, 1.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0, 1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,    -1.0,  1.0, 1.0,  1.0,
            ContactPairInfo::FT_VERTEX,   1.0,  0.0, 0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}



/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test16) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -1.0, 1.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        Vec2 anchor2D_2( 1.0, 0.0);

        Vec2            tiltAxisXY(1.0, 0.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_VERTEX,  -1.0,  1.0, 0.0,  0.0,
            ContactPairInfo::FT_VERTEX,   1.0,  0.0, 0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);


        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test17) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -1.0, 1.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        Vec2 anchor2D_2( 1.0, 0.0);

        Vec2            tiltAxisXY(0.0, -1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_VERTEX,  -1.0,  1.0, 0.0,  0.0,
            ContactPairInfo::FT_VERTEX,   1.0,  0.0, 0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test18) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -1.0, 1.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-1.0, 0.0));
        polygon2D_2.push_back(Vec2( 1.0, 0.0));
        Vec2 anchor2D_2( 1.0, 0.0);

        Vec2            tiltAxisXY(0.0,  1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_VERTEX,  -1.0,  1.0, 0.0,  0.0,
            ContactPairInfo::FT_VERTEX,   1.0,  0.0, 0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test19) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( 0.0, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-2.0, 0.0));
        polygon2D_2.push_back(Vec2( 2.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.5);

        Vec2            tiltAxisXY(1.0,  0.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_FACE,   0.0,  0.0, 0.0,  0.0,
            ContactPairInfo::FT_EDGE,  -2.0,  0.0, 2.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test20) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( 0.0, 0.5);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-2.0, 0.0));
        polygon2D_2.push_back(Vec2( 2.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0,  -1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,  -1.0,  1.0, -2.0,  0.0,
            ContactPairInfo::FT_EDGE,  -2.0,  0.0,  2.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test21) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( 0.0, 0.5);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-2.0, 0.0));
        polygon2D_2.push_back(Vec2( 2.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(0.0,   1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,   1.0,  1.0,  2.0,  0.0,
            ContactPairInfo::FT_EDGE,  -2.0,  0.0,  2.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test22) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( 0.0, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-2.0, 0.0));
        polygon2D_2.push_back(Vec2( 2.0, 0.0));
        Vec2 anchor2D_2( 0.0, 0.0);

        Vec2            tiltAxisXY(1.0,   0.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_FACE,   0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_EDGE,  -2.0,  0.0,  2.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test23) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -2.0, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-2.0, 0.0));
        polygon2D_2.push_back(Vec2( 2.0, 0.0));
        Vec2 anchor2D_2( -2.0, 0.0);

        Vec2            tiltAxisXY(0.0,  -1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_VERTEX,  -2.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_VERTEX,  -2.0,  0.0,  0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test24) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -2.0, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-2.0, 0.0));
        polygon2D_2.push_back(Vec2( 2.0, 0.0));
        Vec2 anchor2D_2( -2.0, 0.0);

        Vec2            tiltAxisXY(0.0,   1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_FACE,     0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_VERTEX,   2.0,  0.0,  0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        //showResult(body1, body2, info);
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}




/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test25) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -1.5, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2( 0.0,  1.0));
        polygon2D_2.push_back(Vec2( 0.0, -1.0));
        Vec2 anchor2D_2(  0.0, 0.0);

        Vec2            tiltAxisXY(0.0,   1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_FACE,     0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_EDGE,     0.0,  1.0,  0.0, -1.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test26) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -1.5, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2( 0.0,  1.0));
        polygon2D_2.push_back(Vec2( 0.0, -1.0));
        Vec2 anchor2D_2(  0.0, 0.0);

        Vec2            tiltAxisXY(1.0,   0.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,    -2.0,  0.0, -1.0, -1.0,
            ContactPairInfo::FT_EDGE,     0.0,  1.0,  0.0, -1.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test27) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -1.5, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2( 0.0,  1.0));
        polygon2D_2.push_back(Vec2( 0.0, -1.0));
        Vec2 anchor2D_2(  0.0, 0.0);

        Vec2            tiltAxisXY(-1.0,   0.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_EDGE,    -2.0,  0.0,  -1.0, 1.0,
            ContactPairInfo::FT_EDGE,     0.0,  1.0,  0.0, -1.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test28) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -2.0, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-0.5,  0.0));
        polygon2D_2.push_back(Vec2( 0.5,  0.0));
        Vec2 anchor2D_2( -0.5, 0.0);

        Vec2            tiltAxisXY( 1.0,   0.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_FACE,     0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_EDGE,    -0.5,  0.0,  0.5,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test29) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -2.0, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-0.5,  0.0));
        polygon2D_2.push_back(Vec2( 0.5,  0.0));
        Vec2 anchor2D_2( -0.5, 0.0);

        Vec2            tiltAxisXY( 0.0,  -1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_VERTEX,    -2.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_VERTEX,    -0.5,  0.0,  0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}


/**  @brief update() finds intersection
 */
TEST_F(ContactUpdater_FACE_EDGETest, Test30) {

    for (long i = 0; i < 1000; i++) {

        vector<Vec2> polygon2D_1;
        polygon2D_1.push_back(Vec2( 1.0,  1.0));
        polygon2D_1.push_back(Vec2(-1.0,  1.0));
        polygon2D_1.push_back(Vec2(-2.0,  0.0));
        polygon2D_1.push_back(Vec2(-1.0, -1.0));
        polygon2D_1.push_back(Vec2( 1.0, -1.0));
        polygon2D_1.push_back(Vec2( 2.0,  0.0));
        Vec2 anchor2D_1( -2.0, 0.0);

        vector<Vec2> polygon2D_2;
        polygon2D_2.push_back(Vec2(-0.5,  0.0));
        polygon2D_2.push_back(Vec2( 0.5,  0.0));
        Vec2 anchor2D_2( -0.5, 0.0);

        Vec2            tiltAxisXY( 0.0,  1.0);
        double          rotationRad = (M_PI /180.0) * 10.0;

        ExpStructFaceEdge exp1(
            ContactPairInfo::FT_FACE,       0.0,  0.0,  0.0,  0.0,
            ContactPairInfo::FT_VERTEX,     0.5,  0.0,  0.0,  0.0 );

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        bool edgeOn2 = i%2==0;

        generateTestPattern(
            polygon2D_1,
            anchor2D_1,
            polygon2D_2,
            anchor2D_2,
            tiltAxisXY,
            rotationRad,
            body1,
            body2,
            info,
            edgeOn2
        );

        rotateTestPatternRandomly(
                     body1, body2, info, randomRotQ(), randVec3D100()*0.05);

        ContactUpdater_FACE_EDGE updater(body1, body2, info, 
                                  EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res01 = updater.update();
        EXPECT_EQ(res01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, edgeOn2), true);

    }
}





}// Makena


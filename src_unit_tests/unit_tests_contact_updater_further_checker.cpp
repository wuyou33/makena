#include "gtest/gtest.h"
#include <vector>
#include "primitives.hpp"
#include "manifold.hpp"
#include "quaternion.hpp"
#include "contact_updater_further_checker.hpp"

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


class ContactUpdaterFurtherCheckerTest : public ::testing::Test {

  protected:  
    ContactUpdaterFurtherCheckerTest(){;}
    virtual ~ContactUpdaterFurtherCheckerTest(){;}
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


static Vec2 randVec2D100()
{
    return Vec2(rand100(), rand100());
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


class ExpStructContactUpdaterFurtherChecker {

  public:

    ExpStructContactUpdaterFurtherChecker(
        enum ContactPairInfo::FeatureType t1, 
        vector<Vec3>& pts1,
        enum ContactPairInfo::FeatureType t2,
        vector<Vec3>& pts2
    ):  mType1(t1),
        mPts1(pts1),
        mType2(t2),
        mPts2(pts2){;}
  
    enum ContactPairInfo::FeatureType mType1;
    vector<Vec3> mPts1;
    enum ContactPairInfo::FeatureType mType2;
    vector<Vec3> mPts2;
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
    ExpStructContactUpdaterFurtherChecker& exp,
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

    return true;
}


static void showResults(
    vector<enum 
           Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER
    >& predMain,
    vector<enum 
           Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER
    >& predAux
) {
    cerr << "predMain:\n";
    for (auto& e : predMain) {
        switch (e) {
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_NONE:
            cerr << "WC_NONE";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_ON_APEX:
            cerr << "WC_ON_APEX";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_ON_EDGE1:
            cerr << "WC_ON_EDGE1";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_ON_EDGE2:
            cerr << "WC_ON_EDGE2";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INSIDE:
            cerr << "WC_INSIDE";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_OUTSIDE:
            cerr << "WC_OUTSIDE";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INTSEC_ON_APEX:
            cerr << "WC_INTSEC_ON_APEX";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INTSEC_ON_EDGE1:
            cerr << "WC_INTSEC_ON_EDGE1";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INTSEC_ON_EDGE2:
            cerr << "WC_INTSEC_ON_EDGE2";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INTSEC_ON_EDGE12:
            cerr << "WC_INTSEC_ON_EDGE12";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_PARSE_RESULT_BEGIN:
            cerr << "WC_PARSE_RESULT_BEGIN";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_COMPLETELY_INSIDE:
            cerr << "WC_COMPLETELY_INSIDE";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_EDGE1_INSIDE_VERTICES:
            cerr << "WC_EDGE1_INSIDE_VERTICES";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_EDGE1_TOUCHING:
            cerr << "WC_EDGE1_TOUCHING";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_EDGE2_INSIDE_VERTICES:
            cerr << "WC_EDGE2_INSIDE_VERTICES";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_EDGE2_TOUCHING:
            cerr << "WC_EDGE2_TOUCHING";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_APEX_ONE_EDGE_CUT_ACROSS:
            cerr << "WC_APEX_ONE_EDGE_CUT_ACROSS";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_APEX_ONE_EDGE_TOUCHING:
            cerr << "WC_APEX_ONE_EDGE_TOUCHING";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_APEX_INSIDE_VERTICES:
            cerr << "WC_APEX_INSIDE_VERTICES";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_APEX_COINCIDENT_VERTEX:
            cerr << "WC_APEX_COINCIDENT_VERTEX";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_END:
            cerr << "WC_END";
            break;
        }
        cerr << "\n";
    }

    cerr << "predAux:\n";
    for (auto& e : predAux) {
        switch (e) {
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_NONE:
            cerr << "WC_NONE";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_ON_APEX:
            cerr << "WC_ON_APEX";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_ON_EDGE1:
            cerr << "WC_ON_EDGE1";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_ON_EDGE2:
            cerr << "WC_ON_EDGE2";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INSIDE:
            cerr << "WC_INSIDE";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_OUTSIDE:
            cerr << "WC_OUTSIDE";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INTSEC_ON_APEX:
            cerr << "WC_INTSEC_ON_APEX";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INTSEC_ON_EDGE1:
            cerr << "WC_INTSEC_ON_EDGE1";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INTSEC_ON_EDGE2:
            cerr << "WC_INTSEC_ON_EDGE2";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INTSEC_ON_EDGE12:
            cerr << "WC_INTSEC_ON_EDGE12";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_PARSE_RESULT_BEGIN:
            cerr << "WC_PARSE_RESULT_BEGIN";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_COMPLETELY_INSIDE:
            cerr << "WC_COMPLETELY_INSIDE";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_EDGE1_INSIDE_VERTICES:
            cerr << "WC_EDGE1_INSIDE_VERTICES";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_EDGE1_TOUCHING:
            cerr << "WC_EDGE1_TOUCHING";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_EDGE2_INSIDE_VERTICES:
            cerr << "WC_EDGE2_INSIDE_VERTICES";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_EDGE2_TOUCHING:
            cerr << "WC_EDGE2_TOUCHING";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_APEX_ONE_EDGE_CUT_ACROSS:
            cerr << "WC_APEX_ONE_EDGE_CUT_ACROSS";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_APEX_ONE_EDGE_TOUCHING:
            cerr << "WC_APEX_ONE_EDGE_TOUCHING";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_APEX_INSIDE_VERTICES:
            cerr << "WC_APEX_INSIDE_VERTICES";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_APEX_COINCIDENT_VERTEX:
            cerr << "WC_APEX_COINCIDENT_VERTEX";
            break;
          case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_END:
            cerr << "WC_END";
            break;
        }
        cerr << "\n";
    }
}


static void showResult(
    enum Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER res
) {
    switch (res) {
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_NONE:
        cerr << "WC_NONE";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_ON_APEX:
        cerr << "WC_ON_APEX";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_ON_EDGE1:
        cerr << "WC_ON_EDGE1";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_ON_EDGE2:
        cerr << "WC_ON_EDGE2";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INSIDE:
        cerr << "WC_INSIDE";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_OUTSIDE:
        cerr << "WC_OUTSIDE";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INTSEC_ON_APEX:
        cerr << "WC_INTSEC_ON_APEX";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INTSEC_ON_EDGE1:
        cerr << "WC_INTSEC_ON_EDGE1";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INTSEC_ON_EDGE2:
        cerr << "WC_INTSEC_ON_EDGE2";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_INTSEC_ON_EDGE12:
        cerr << "WC_INTSEC_ON_EDGE12";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_PARSE_RESULT_BEGIN:
        cerr << "WC_PARSE_RESULT_BEGIN";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_COMPLETELY_INSIDE:
        cerr << "WC_COMPLETELY_INSIDE";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_EDGE1_INSIDE_VERTICES:
        cerr << "WC_EDGE1_INSIDE_VERTICES";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_EDGE1_TOUCHING:
        cerr << "WC_EDGE1_TOUCHING";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_EDGE2_INSIDE_VERTICES:
        cerr << "WC_EDGE2_INSIDE_VERTICES";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_EDGE2_TOUCHING:
        cerr << "WC_EDGE2_TOUCHING";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_APEX_ONE_EDGE_CUT_ACROSS:
        cerr << "WC_APEX_ONE_EDGE_CUT_ACROSS";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_APEX_ONE_EDGE_TOUCHING:
        cerr << "WC_APEX_ONE_EDGE_TOUCHING";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_APEX_INSIDE_VERTICES:
        cerr << "WC_APEX_INSIDE_VERTICES";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_APEX_COINCIDENT_VERTEX:
        cerr << "WC_APEX_COINCIDENT_VERTEX";
        break;
      case Makena::ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER::WC_END:
        cerr << "WC_END";
        break;
    }
    cerr << "\n";
}


/**  @brief checkFeatures_FACE_EDGE() no penetration
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test01) {

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
        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();

        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);

        EXPECT_EQ(res_01, false);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() touching f-f #1.
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test02) {

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() touching f-f #2.
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test03) {

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() shallow penetration f-f #1.
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test04) {

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

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(1.0/2.0) * 1.01;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() shallow penetration f-f #2.
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test05) {

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

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis( 1.0, 0.0, 0.0);
        double     rotationRad = atan(1.0/2.0) * 1.01;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() deep penetration f-f #1.
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test06) {

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

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(1.0/2.0) + M_PI * 0.49;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() deep penetration f-f #2.
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test07) {

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

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis( 1.0, 0.0, 0.0);
        double     rotationRad = atan(1.0/2.0) + M_PI * 0.49;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() total immersion
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test08) {

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

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad =
                       M_PI + atan(1.0/2.0) * (rand100()/50.0 - 1.0) * 0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() too deep #1
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test09) {

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

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(1.0/2.0) + M_PI * 0.51;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() too deep #2
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test10) {

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

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis( 1.0, 0.0, 0.0);
        double     rotationRad = atan(1.0/2.0) + M_PI * 0.51;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() no penetration obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test11) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));

        Vec3       anchor_2(0.0,  0.0, -0.01);

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() touching f-f #1. obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test12) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(0.01/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() touching f-f #2. obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test13) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(1.0, 0.0, 0.0);
        double     rotationRad = atan(0.01/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() shallow penetration f-f #1. obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test14) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(0.01/2.0) * 1.01;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() shallow penetration f-f #2.obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test15) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis( 1.0, 0.0, 0.0);
        double     rotationRad = atan(0.01/2.0) * 1.01;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() deep penetration f-f #1.obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test16) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(0.01/2.0) + M_PI * 0.49;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() deep penetration f-f #2.obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test17) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis( 1.0, 0.0, 0.0);
        double     rotationRad = atan(0.01/2.0) + M_PI * 0.49;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() total immersion obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test18) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad =
                       M_PI + atan(0.01/2.0) * (rand100()/50.0 - 1.0) * 0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() too deep #1 obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test19) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(0.01/2.0) + M_PI * 0.51;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() too deep #2 obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test20) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis( 1.0, 0.0, 0.0);
        double     rotationRad = atan(0.01/2.0) + M_PI * 0.51;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() no penetration accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test21) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() touching f-f #1. accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test22) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(100.0/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() touching f-f #2. accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test23) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(1.0, 0.0, 0.0);
        double     rotationRad = atan(100.0/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() shallow penetration f-f #1. accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test24) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(100.0/2.0) * 1.01;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() shallow penetration f-f #2.accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test25) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis( 1.0, 0.0, 0.0);
        double     rotationRad = atan(100.0/2.0) * 1.01;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() deep penetration f-f #1.accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test26) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(100.0/2.0) + 
                                 atan(2.0/100.0)*(rand100()/100.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() deep penetration f-f #2.accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test27) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis( 1.0, 0.0, 0.0);
        double     rotationRad = atan(100.0/2.0) + 
                                 atan(2.0/100.0)*(rand100()/100.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() total immersion accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test28) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad =
                       M_PI + atan(100.0/2.0) * (rand100()/50.0 - 1.0) * 0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() too deep #1 accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test29) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis(-1.0, 0.0, 0.0);
        double     rotationRad = atan(100.0/2.0) + M_PI * 0.501;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_EDGE() too deep #2 accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test30) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis( 1.0, 0.0, 0.0);
        double     rotationRad = atan(100.0/2.0) + M_PI * 0.501;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);

        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() no penetration
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test31) {

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

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double     yawTiltAxis = 2.0 * M_PI * (rand100()/100.0) - M_PI;
        Vec3       tiltAxis( cos(yawTiltAxis), sin(yawTiltAxis), 0.0);
        double     rotationRad = atan(1.0/sqrt(8.0)) * (rand100()/100.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_VERTEX, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() touching FACE-EDGE
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test32) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0, 1.0, 0.0);
        Vec3       tiltAxis02(-1.0, 1.0, 0.0);
        Vec3       tiltAxis03(-1.0,-1.0, 0.0);
        Vec3       tiltAxis04( 1.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }
        double     rotationRad = atan(1.0/sqrt(8.0));

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() FACE-EDGE deep penetration
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test33) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0, 1.0, 0.0);
        Vec3       tiltAxis02(-1.0, 1.0, 0.0);
        Vec3       tiltAxis03(-1.0,-1.0, 0.0);
        Vec3       tiltAxis04( 1.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }
        double     rotationRad = atan(1.0/sqrt(8.0)) + 
                                 (M_PI*0.5 - atan(1.0/sqrt(8.0)))*0.99 *
                                 (rand100()/100.0);
        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() FACE-FACE touching
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test34) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0,  0.0,  0.0);
        Vec3       tiltAxis02( 0.0,  1.0,  0.0);
        Vec3       tiltAxis03(-1.0,  0.0,  0.0);
        Vec3       tiltAxis04( 0.0, -1.0,  0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }
        double     rotationRad = atan(1.0/2.0);
        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() FACE-FACE deep penetration
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test35) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0,  0.0,  0.0);
        Vec3       tiltAxis02( 0.0,  1.0,  0.0);
        Vec3       tiltAxis03(-1.0,  0.0,  0.0);
        Vec3       tiltAxis04( 0.0, -1.0,  0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }

        double     rotationRad = atan(1.0/2.0) + M_PI*0.5*(rand100()/100.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() too deep penetration to run GJK.
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test36) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0,  0.0,  0.0);
        Vec3       tiltAxis02( 0.0,  1.0,  0.0);
        Vec3       tiltAxis03(-1.0,  0.0,  0.0);
        Vec3       tiltAxis04( 0.0, -1.0,  0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }

        double     rotationRad = atan(1.0/2.0) + M_PI*0.5 + 
                                 (M_PI - ((atan(1.0/2.0) + M_PI*0.5)))*
                                 rand100()/100.0;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() too deep penetration to run GJK #2
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test37) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0, 1.0, 0.0);
        Vec3       tiltAxis02(-1.0, 1.0, 0.0);
        Vec3       tiltAxis03(-1.0,-1.0, 0.0);
        Vec3       tiltAxis04( 1.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }
        double     rotationRad = (atan(1.0/sqrt(8.0)) + M_PI*0.5) +
                                 (M_PI - (atan(1.0/sqrt(8.0)) + M_PI*0.5)) *
                                 ((rand100()+1.0)/101.0);
        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() no penetration obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test38) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -0.01));

        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -0.01));

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double     yawTiltAxis = 2.0 * M_PI * (rand100()/100.0) - M_PI;
        Vec3       tiltAxis( cos(yawTiltAxis), sin(yawTiltAxis), 0.0);
        double     rotationRad = atan(0.01/sqrt(8.0)) * (rand100()/100.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_VERTEX, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() touching FACE-EDGE obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test39) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -0.01));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0, 1.0, 0.0);
        Vec3       tiltAxis02(-1.0, 1.0, 0.0);
        Vec3       tiltAxis03(-1.0,-1.0, 0.0);
        Vec3       tiltAxis04( 1.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }
        double     rotationRad = atan(0.01/sqrt(8.0));

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() FACE-EDGE deep penetration obtuse ang
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test40) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -0.01));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0, 1.0, 0.0);
        Vec3       tiltAxis02(-1.0, 1.0, 0.0);
        Vec3       tiltAxis03(-1.0,-1.0, 0.0);
        Vec3       tiltAxis04( 1.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }
        double     rotationRad = atan(0.01/sqrt(8.0)) + 
                                 (M_PI*0.5 - atan(0.01/sqrt(8.0)))*0.99 *
                                 (rand100()/100.0);
        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() FACE-FACE touching obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test41) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -0.01));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0,  0.0,  0.0);
        Vec3       tiltAxis02( 0.0,  1.0,  0.0);
        Vec3       tiltAxis03(-1.0,  0.0,  0.0);
        Vec3       tiltAxis04( 0.0, -1.0,  0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }
        double     rotationRad = atan(0.01/2.0);
        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() FACE-FACE deep penetration obtuse ang
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test42) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -0.01));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0,  0.0,  0.0);
        Vec3       tiltAxis02( 0.0,  1.0,  0.0);
        Vec3       tiltAxis03(-1.0,  0.0,  0.0);
        Vec3       tiltAxis04( 0.0, -1.0,  0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }

        double     rotationRad = atan(0.01/2.0) + M_PI*0.5*(rand100()/100.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() too deep penetration to run GJK. obtus
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test43) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -0.01));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0,  0.0,  0.0);
        Vec3       tiltAxis02( 0.0,  1.0,  0.0);
        Vec3       tiltAxis03(-1.0,  0.0,  0.0);
        Vec3       tiltAxis04( 0.0, -1.0,  0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }

        double     rotationRad = atan(0.01/2.0) + M_PI*0.5 + 
                                 (M_PI - ((atan(0.01/2.0) + M_PI*0.5)))*
                                 rand100()/100.0;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() too deep penetration to run GJK #2
 *          obtuse angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test44) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -0.01));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0, 1.0, 0.0);
        Vec3       tiltAxis02(-1.0, 1.0, 0.0);
        Vec3       tiltAxis03(-1.0,-1.0, 0.0);
        Vec3       tiltAxis04( 1.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }
        double     rotationRad = (atan(0.01/sqrt(8.0)) + M_PI*0.5) +
                                 (M_PI - (atan(0.01/sqrt(8.0)) + M_PI*0.5)) *
                                 ((rand100()+1.0)/101.0);
        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() no penetration accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test45) {

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -100.0));

        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -100.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        double     yawTiltAxis = 2.0 * M_PI * (rand100()/100.0) - M_PI;
        Vec3       tiltAxis( cos(yawTiltAxis), sin(yawTiltAxis), 0.0);
        double     rotationRad =atan(100.0/sqrt(8.0)) * (rand100()/100.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_VERTEX, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() touching FACE-EDGE accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test46) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0, 1.0, 0.0);
        Vec3       tiltAxis02(-1.0, 1.0, 0.0);
        Vec3       tiltAxis03(-1.0,-1.0, 0.0);
        Vec3       tiltAxis04( 1.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }
        double     rotationRad = atan(100.0/sqrt(8.0));

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() FACE-EDGE deep penetration accute ang
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test47) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0, 1.0, 0.0);
        Vec3       tiltAxis02(-1.0, 1.0, 0.0);
        Vec3       tiltAxis03(-1.0,-1.0, 0.0);
        Vec3       tiltAxis04( 1.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }
        double     rotationRad = atan(100.0/sqrt(8.0)) + 
                                 (M_PI*0.5 - atan(100.0/sqrt(8.0)))*0.99 *
                                 (rand100()/100.0);
        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() FACE-FACE touching accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test48) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0,  0.0,  0.0);
        Vec3       tiltAxis02( 0.0,  1.0,  0.0);
        Vec3       tiltAxis03(-1.0,  0.0,  0.0);
        Vec3       tiltAxis04( 0.0, -1.0,  0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }
        double     rotationRad = atan(100.0/2.0);
        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() FACE-FACE deep penetration accute ang
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test49) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0,  0.0,  0.0);
        Vec3       tiltAxis02( 0.0,  1.0,  0.0);
        Vec3       tiltAxis03(-1.0,  0.0,  0.0);
        Vec3       tiltAxis04( 0.0, -1.0,  0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }

        double     rotationRad = atan(100.0/2.0) + 
                                 atan(1.0/100.0)*(rand100()/100.0) * 0.9;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() too deep penetration to run GJK. accut
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test50) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0,  0.0,  0.0);
        Vec3       tiltAxis02( 0.0,  1.0,  0.0);
        Vec3       tiltAxis03(-1.0,  0.0,  0.0);
        Vec3       tiltAxis04( 0.0, -1.0,  0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }

        double     rotationRad = atan(100.0/2.0) + M_PI*0.5 + 
                                 (M_PI - ((atan(100.0/2.0) + M_PI*0.5)))*
                                 rand100()/100.0;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_FACE_VERTEX() too deep penetration to run GJK #2
 *          accute angle
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test51) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%4;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else if (type01==1) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else if (type01==2) {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3       tiltAxis01( 1.0, 1.0, 0.0);
        Vec3       tiltAxis02(-1.0, 1.0, 0.0);
        Vec3       tiltAxis03(-1.0,-1.0, 0.0);
        Vec3       tiltAxis04( 1.0,-1.0, 0.0);
        Vec3       tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else if (type01==1) {
            tiltAxis = tiltAxis02;
        }
        else if (type01==2) {
            tiltAxis = tiltAxis03;
        }
        else {
            tiltAxis = tiltAxis04;
        }
        double     rotationRad = (atan(100.0/sqrt(8.0)) + M_PI*0.5) +
                                 (M_PI - (atan(100.0/sqrt(8.0)) + M_PI*0.5)) *
                                 ((rand100()+1.0)/101.0);
        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE,   expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() no penetration
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test52) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -1.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -1.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -1.0));


        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 0.0, 1.0, 0.0);
        Vec3 tiltAxis;
        if (type01 == 0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double     rotationRad = atan(1.0/2.0) * (rand100()/100.0 - 1.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE,   expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() touching edge-edge
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test53) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

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
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -1.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -1.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -1.0));


        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);

        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(1.0/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() touching edge-edge
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test54) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -1.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -1.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -1.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01( 0.0, 1.0, 0.0);
        Vec3 tiltAxis02( 0.0,-1.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(1.0/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() edge-edge deep penetration
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test55) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

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
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -1.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -1.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -1.0));


        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(1.0/2.0) + M_PI*0.5*
                             (rand100()/100.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() edge-edge deep penetration
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test56) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -1.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -1.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -1.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01( 0.0, 1.0, 0.0);
        Vec3 tiltAxis02( 0.0,-1.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(1.0/2.0) + M_PI*0.5*
                             (rand100()/100.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() too deep
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test57) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -1.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -1.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -1.0));

        vector<Vec3> expectedPoints_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = (atan(1.0/2.0) + M_PI*0.5) +
                             (M_PI - ((atan(1.0/2.0) + M_PI*0.5))) *
                             (rand100()/100.0)*0.99;

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() edge-edge deep penetration
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test58) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

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

        Vec3       anchor_1(0.0,  0.0,  1.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -1.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -1.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -1.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -1.0));

        vector<Vec3> expectedPoints_2_1;

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01( 0.0, 1.0, 0.0);
        Vec3 tiltAxis02( 0.0,-1.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = (atan(1.0/2.0) + M_PI*0.5) +
                             (M_PI - ((atan(1.0/2.0) + M_PI*0.5))) *
                             (rand100()/100.0)*0.99;

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() no penetration obtuse-obtuse
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test59) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -0.01));
        polygon_2.push_back(Vec3( 0.0, -2.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -0.01));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -0.01));


        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 0.0, 1.0, 0.0);
        Vec3 tiltAxis;
        if (type01 == 0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double     rotationRad = atan(0.01/2.0) * (rand100()/100.0 - 1.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE,   expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() touching edge-edge 
 *          obtuse-obtuse
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test60) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -0.01));
        polygon_2.push_back(Vec3( 0.0, -2.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -0.01));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -0.01));


        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);

        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(0.01/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() touching edge-edge
 *          obtuse-obtuse
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test61) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -0.01));
        polygon_2.push_back(Vec3( 0.0, -2.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -0.01));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -0.01));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01( 0.0, 1.0, 0.0);
        Vec3 tiltAxis02( 0.0,-1.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(0.01/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() edge-edge deep penetration
 *          obtuse-obtuse
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test62) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -0.01));
        polygon_2.push_back(Vec3( 0.0, -2.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -0.01));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -0.01));


        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(0.01/2.0) + M_PI*0.5*
                             (rand100()/100.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() edge-edge deep penetration
 *          obtuse-obtuse
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test63) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -0.01));
        polygon_2.push_back(Vec3( 0.0, -2.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -0.01));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -0.01));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -0.01));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01( 0.0, 1.0, 0.0);
        Vec3 tiltAxis02( 0.0,-1.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(0.01/2.0) + M_PI*0.5*
                             (rand100()/100.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() too deep
 *          obtuse-obtuse
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test64) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -0.01));
        polygon_2.push_back(Vec3( 0.0, -2.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -0.01));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -0.01));

        vector<Vec3> expectedPoints_2_1;

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = (atan(0.01/2.0) + M_PI*0.5) +
                             (M_PI - ((atan(0.01/2.0) + M_PI*0.5))) *
                             (rand100()/100.0)*0.99;

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() edge-edge deep penetration
 *          obtuse-obtuse
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test65) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -0.01));
        polygon_2.push_back(Vec3( 0.0, -2.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -0.01));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -0.01));

        vector<Vec3> expectedPoints_2_1;

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01( 0.0, 1.0, 0.0);
        Vec3 tiltAxis02( 0.0,-1.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = (atan(0.01/2.0) + M_PI*0.5) +
                             (M_PI - ((atan(0.01/2.0) + M_PI*0.5))) *
                             (rand100()/100.0)*0.99;

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() no penetration
 *  accute-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test66) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  100.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  100.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  100.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  100.0));

        Vec3       anchor_1(0.0,  0.0,  100.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -100.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -100.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -100.0));


        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 0.0, 1.0, 0.0);
        Vec3 tiltAxis;
        if (type01 == 0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double     rotationRad =atan(100.0/2.0)*(rand100()/100.0-1.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE,   expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() touching edge-edge 
 *          accute-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test67) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  100.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  100.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  100.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  100.0));
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  100.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -100.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -100.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -100.0));


        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);

        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(100.0/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() touching edge-edge
 *          accute-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test68) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  100.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  100.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  100.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  100.0));

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -100.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -100.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01( 0.0, 1.0, 0.0);
        Vec3 tiltAxis02( 0.0,-1.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(100.0/2.0);

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() edge-edge deep penetration
 *          accute-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test69) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  100.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  100.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  100.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  100.0));
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  100.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -100.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -100.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -100.0));


        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(100.0/2.0) + atan(2.0/100.0)*
                             (rand100()/100.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() edge-edge deep penetration
 *          accute-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test70) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  100.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  100.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  100.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  100.0));

        Vec3       anchor_1(0.0,  0.0,  100.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -100.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -100.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -100.0));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01( 0.0, 1.0, 0.0);
        Vec3 tiltAxis02( 0.0,-1.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(100.0/2.0) + atan(2.0/100.0)*
                             (rand100()/100.0)*0.99;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() too deep
 *          accute-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test71) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  100.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  100.0));

        vector<Vec3> expectedPoints_1_1;

        Vec3       anchor_1(0.0,  0.0,  100.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -100.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -100.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -100.0));

        vector<Vec3> expectedPoints_2_1;

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(100.0/2.0) + M_PI * 0.5 + 
                             atan(2.0/100.0)*(rand100()/100.0)*0.99;

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_CROSSING() edge-edge too deep
 *          accute-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test72) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  100.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  100.0));

        vector<Vec3> expectedPoints_1_1;

        Vec3       anchor_1(0.0,  0.0,  100.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  2.0, -100.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  2.0, -100.0));
        activePoints_2.push_back(Vec3( 0.0, -2.0, -100.0));

        vector<Vec3> expectedPoints_2_1;

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01( 0.0, 1.0, 0.0);
        Vec3 tiltAxis02( 0.0,-1.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(100.0/2.0) + M_PI * 0.5 + 
                             atan(2.0/100.0)*(rand100()/100.0)*0.99;

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

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() edge-edge no penetration
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test73) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

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

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(1.0/2.0)*2.0 * (rand100()/100.0)*0.99;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() face-face touching
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test74) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

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
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

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
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(1.0/2.0)*2.0;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() face-face deep penetration
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test75) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

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
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

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
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(1.0/2.0)*2.0 + 
                             M_PI*0.5 * (rand100()/100.0)*0.99 ;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() too deep.
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test76) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

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
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

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
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = (atan(1.0/2.0)*2.0 + M_PI*0.5) + 
                             (M_PI - (atan(1.0/2.0)*2.0 + M_PI*0.5)) *
                             (rand100()/100.0) * 0.99 ;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}






















/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() edge-edge no penetration
 *          obtuse-obtuse
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test77) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(0.01/2.0)*2.0 * (rand100()/100.0)*0.99;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() face-face touching
 *          obtuse-obtuse
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test78) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(0.01/2.0)*2.0;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() face-face deep penetration
 *          obtuse-obtuse
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test79) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(0.01/2.0)*2.0 + 
                             M_PI*0.5 * (rand100()/100.0)*0.99 ;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() too deep.
 *          obtuse-obtuse
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test80) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3(-2.0,  0.0, -0.01));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -0.01));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -0.01));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -0.01));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -0.01));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -0.01);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = (atan(0.01/2.0)*2.0 + M_PI*0.5) + 
                             (M_PI - (atan(0.01/2.0)*2.0 + M_PI*0.5)) *
                             (rand100()/100.0) * 0.99 ;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}





























/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() edge-edge no penetration
 *          accute-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test81) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  100.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  100.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  100.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  100.0));

        Vec3       anchor_1(0.0,  0.0,  100.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(100.0/2.0)*2.0 * (rand100()/100.0)*0.99;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() face-face touching
 *          accute-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test82) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  100.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  100.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  100.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  100.0));
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  100.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(100.0/2.0)*2.0;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() face-face deep penetration
 *          accute-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test83) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  100.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  100.0));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  100.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  100.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  100.0));
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  100.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(100.0/2.0)*2.0 +
                             atan(2.0/100.0)*2.0 * (rand100()/100.0)*0.99 ;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        auto randVec = randVec3D100()*0.05;
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        // showResult(body1, body2, info);
    }
}



















/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() edge-edge no penetration
 *          obtuse-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test84) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = 
                 (atan(100.0/2.0)+atan(0.01/2.0)) * (rand100()/100.0)*0.99;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() face-face touching
 *          obtuse-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test85) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = atan(100.0/2.0) +  atan(0.01/2.0);

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() face-face deep penetration
 *          obtuse-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test86) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = (atan(100.0/2.0) + atan(0.01/2.0)) + 
                             (M_PI - (atan(100.0/2.0) + atan(0.01/2.0)))*
                             (rand100()/100.0)*0.99;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_EDGE_EDGE_PARALLEL() too deep.
 *          obtuse-accute
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test87) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);
        long type01 = long(rand100())%2;

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3(-2.0,  0.0,  0.01));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -1.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -1.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -1.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 2.0,  0.0,  0.01));
        activePoints_1.push_back(Vec3(-2.0,  0.0,  0.01));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0,  0.01));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0,  0.01));
        if (type01==0) {
            expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_1(0.0,  0.0,  0.01);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3(-2.0,  0.0, -100.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  1.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  1.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 2.0,  0.0, -100.0));
        activePoints_2.push_back(Vec3(-2.0,  0.0, -100.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0, -100.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  0.0, -100.0));
        if (type01==0) {
            expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));
        }
        else {
            expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));
            expectedPoints_2_1.push_back(Vec3(-2.0, -2.0,  0.0));
        }

        Vec3       anchor_2(0.0,  0.0, -100.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis01(-1.0, 0.0, 0.0);
        Vec3 tiltAxis02( 1.0, 0.0, 0.0);
        Vec3 tiltAxis;
        if (type01==0) {
            tiltAxis = tiltAxis01;
        }
        else {
            tiltAxis = tiltAxis02;
        }

        double rotationRad = M_PI;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);
        info.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;

        auto randQ =  randomRotQ();           
        rotateTestPatternRandomly(
                     body1, body2, info, randQ, randVec3D100()*0.05);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        //EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


static Vec2 rotate(const Vec2& v, double a)
{
    double c = cos(a);
    double s = sin(a);
    return Vec2(v.x() * c - v.y() * s, 
                v.x() * s + v.y() * c );
}

static void roatateTestPatternrandomly(
    vector<Vec2>& fan,
    Vec2&         apex,
    Vec2&         side1,
    Vec2&         side2,
    double        ang,
    Vec2&         trans
) {
    for (auto& v : fan) {
        v = rotate(v + trans, ang);
    }
    apex = rotate(apex + trans, ang);
    side1 = rotate(side1, ang);
    side2 = rotate(side2, ang);
}


static bool checkResult(
    vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
    >& v1,

    vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
    >& v2

) {
    if (v1.size()!=v2.size()) {
        return false;
    }
    for (long i = 0 ; i < v1.size() ; i++) {
        if (v1[i]!=v2[i]) {
            return false;
        }
    }
    return true;
}


/**  @brief classifyPointsAgainstWedge() outside #1
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test88) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 2.0, 1.0));
        fan.push_back(Vec2( 2.0, 2.0));
        fan.push_back(Vec2(-2.0, 1.0));
        fan.push_back(Vec2(-2.0, 0.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);

        std::reverse(fan.begin(), fan.end());
        std::reverse(predMainExp01.begin(), predMainExp01.end());
        std::reverse(predAuxExp01.begin(),  predAuxExp01.end());

        auto res_02 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_02, false);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() outside #2
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test89) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( -5.0,-4.0));
        fan.push_back(Vec2( -4.0,-3.0));
        fan.push_back(Vec2( -2.0,-1.0));
        fan.push_back(Vec2(  0.0, 1.0));
        fan.push_back(Vec2(  2.0,-1.0));
        fan.push_back(Vec2(  4.0,-3.0));
        fan.push_back(Vec2(  5.0,-4.0));
        fan.push_back(Vec2(  5.0, 5.0));
        fan.push_back(Vec2( -5.0, 5.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

        std::reverse(fan.begin(), fan.end());
        std::reverse(predMainExp01.begin(), predMainExp01.end());
        std::reverse(predAuxExp01.begin(),  predAuxExp01.end());

        auto res_02 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_02, false);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() on apex
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test90) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( -1.0, 0.0));
        fan.push_back(Vec2(  1.0, 0.0));
        fan.push_back(Vec2(  1.0, 1.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INTSEC_ON_APEX);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() touching apex
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test91) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(  0.0, 0.0));
        fan.push_back(Vec2(  1.0, 0.0));
        fan.push_back(Vec2(  1.0, 1.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_APEX);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() colinear to side 1
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test92) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( -1.0,  1.0));
        fan.push_back(Vec2(  1.0, -1.0));
        fan.push_back(Vec2(  1.0,  1.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE1);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_EDGE1_TOUCHING);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() colinear to side 2
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test93) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(  1.0,  1.0));
        fan.push_back(Vec2( -1.0, -1.0));
        fan.push_back(Vec2( -1.0,  1.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE2);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_EDGE2_TOUCHING);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() colinear to side 1
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test94) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(  1.0, -1.0));
        fan.push_back(Vec2(  2.0, -2.0));
        fan.push_back(Vec2(  1.0,  1.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE1);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE1);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_EDGE1_TOUCHING);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() colinear to side 2
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test95) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( -1.0, -1.0));
        fan.push_back(Vec2( -2.0, -2.0));
        fan.push_back(Vec2( -1.0,  1.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE2);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE2);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_EDGE2_TOUCHING);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() touching side 1 from outside
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test96) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(  1.0, -1.0));
        fan.push_back(Vec2(  2.0, -1.0));
        fan.push_back(Vec2(  1.0,  1.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE1);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() touching side 1 from inside
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test97) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(  1.0, -1.0));
        fan.push_back(Vec2(  0.0, -1.0));
        fan.push_back(Vec2(  0.0, -2.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE1);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() touching side 1 and side 2 inside
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test98) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(  1.0, -1.0));
        fan.push_back(Vec2( -1.0, -1.0));
        fan.push_back(Vec2(  0.0, -2.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE1);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE2);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() touching side 1 and intersects side 2
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test99) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(  1.0, -1.0));
        fan.push_back(Vec2( -2.0, -1.0));
        fan.push_back(Vec2( -2.0,  1.0));
        fan.push_back(Vec2(  1.0,  1.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE1);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INTSEC_ON_EDGE2);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() touching side 2 from outside
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test100) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( -1.0, -1.0));
        fan.push_back(Vec2( -2.0, -1.0));
        fan.push_back(Vec2( -1.0,  1.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE2);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() touching side 2 from inside
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test101) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( -1.0, -1.0));
        fan.push_back(Vec2(  0.0, -1.0));
        fan.push_back(Vec2(  0.0, -2.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE2);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() touching side 2 and intersects side 1
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test102) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( -1.0, -1.0));
        fan.push_back(Vec2(  2.0, -1.0));
        fan.push_back(Vec2(  2.0,  1.0));
        fan.push_back(Vec2( -1.0,  1.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_ON_EDGE2);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INTSEC_ON_EDGE1);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() intersects apex from outside into inside
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test103) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(  0.5,  1.0));
        fan.push_back(Vec2( -0.5, -1.0));
        fan.push_back(Vec2(  2.0,  0.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INTSEC_ON_APEX);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INTSEC_ON_EDGE1);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() intersects side 1
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test104) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(  2.0, -1.0));
        fan.push_back(Vec2(  0.0, -1.0));
        fan.push_back(Vec2(  1.0,  0.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INTSEC_ON_EDGE1);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INTSEC_ON_EDGE1);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() intersects side 2
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test105) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( -2.0, -1.0));
        fan.push_back(Vec2(  0.0, -1.0));
        fan.push_back(Vec2( -1.0,  0.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INTSEC_ON_EDGE2);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INTSEC_ON_EDGE2);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief classifyPointsAgainstWedge() intersects both side 1 and side 2
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test106) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( -2.0, -1.0));
        fan.push_back(Vec2(  2.0, -1.0));
        fan.push_back(Vec2(  0.0, 10.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMainExp01, predAuxExp01;
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predMainExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_OUTSIDE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_INTSEC_ON_EDGE12);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        predAuxExp01.push_back(Makena::ContactUpdaterFurtherChecker::WC_NONE);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan,
            apex,
            side1,
            side1.perp(),
            side2,
            (side2.perp())*-1.0,
            predMain,
            predAux  
        );

        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(predMain, predMainExp01), true);
        EXPECT_EQ(checkResult(predAux,  predAuxExp01),  true);
        
        //showResults(predMain, predAux);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() outside
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test107) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 2.0, 1.0));
        fan.push_back(Vec2( 2.0, 2.0));
        fan.push_back(Vec2(-2.0, 1.0));
        fan.push_back(Vec2(-2.0, 0.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, false);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_NONE);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() inside
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test108) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 1.0,-1.0));
        fan.push_back(Vec2(-1.0,-1.0));
        fan.push_back(Vec2( 0.0,-2.0));
        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_COMPLETELY_INSIDE);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() tricky 1
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test109) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 1.0, 0.0));
        fan.push_back(Vec2(-1.0, 0.0));
        fan.push_back(Vec2(-1.0,-2.0));
        fan.push_back(Vec2( 1.0,-2.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_APEX_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() edge touching apex
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test110) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 1.0, 0.0));
        fan.push_back(Vec2(-1.0, 0.0));
        fan.push_back(Vec2(-1.0, 2.0));
        fan.push_back(Vec2( 1.0, 2.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_APEX_ONE_EDGE_TOUCHING);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() vertex touching apex
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test111) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 0.0, 0.0));
        fan.push_back(Vec2(-1.0, 2.0));
        fan.push_back(Vec2( 1.0, 2.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_APEX_COINCIDENT_VERTEX);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test112) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 1.0,  0.0));
        fan.push_back(Vec2( 0.0, -1.0));
        fan.push_back(Vec2( 2.0, -1.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE1_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test113) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(-1.0,  0.0));
        fan.push_back(Vec2( 0.0, -1.0));
        fan.push_back(Vec2(-2.0, -1.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE2_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test114) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 1.0,  0.0));
        fan.push_back(Vec2( 0.0, -1.0));
        fan.push_back(Vec2( 1.0, -2.0));
        fan.push_back(Vec2( 4.0, -2.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE1_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test115) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(-1.0,  0.0));
        fan.push_back(Vec2( 0.0, -1.0));
        fan.push_back(Vec2(-1.0, -2.0));
        fan.push_back(Vec2(-4.0, -2.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE2_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test116) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 1.0,  0.0));
        fan.push_back(Vec2( 1.0, -1.0));
        fan.push_back(Vec2( 4.0, -2.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE1_TOUCHING);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test117) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(-1.0,  0.0));
        fan.push_back(Vec2(-1.0, -1.0));
        fan.push_back(Vec2(-4.0, -2.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE2_TOUCHING);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test118) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 2.0,  0.0));
        fan.push_back(Vec2( 1.0, -1.0));
        fan.push_back(Vec2( 1.0, -2.0));
        fan.push_back(Vec2( 5.0, -2.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE1_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test119) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(-2.0,  0.0));
        fan.push_back(Vec2(-1.0, -1.0));
        fan.push_back(Vec2(-1.0, -2.0));
        fan.push_back(Vec2(-5.0, -2.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE2_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test120) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 2.0,  0.0));
        fan.push_back(Vec2( 1.0, -1.0));
        fan.push_back(Vec2( 2.0, -2.0));
        fan.push_back(Vec2( 5.0, -2.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE1_TOUCHING);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test121) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(-2.0,  0.0));
        fan.push_back(Vec2(-1.0, -1.0));
        fan.push_back(Vec2(-2.0, -2.0));
        fan.push_back(Vec2(-5.0, -2.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE2_TOUCHING);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test122) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 0.0,  1.0));
        fan.push_back(Vec2( 0.0, -1.0));
        fan.push_back(Vec2( 5.0, -1.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE1_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test123) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 0.0,  1.0));
        fan.push_back(Vec2( 0.0, -1.0));
        fan.push_back(Vec2(-5.0, -1.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE2_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test124) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2( 1.0,  1.0));
        fan.push_back(Vec2( 0.0,  0.0));
        fan.push_back(Vec2( 0.0, -1.0));
        fan.push_back(Vec2( 5.0, -1.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE1_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test125) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(-1.0,  1.0));
        fan.push_back(Vec2( 0.0,  0.0));
        fan.push_back(Vec2( 0.0, -1.0));
        fan.push_back(Vec2(-5.0, -1.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_EDGE2_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test126) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(-1.0,  1.0));
        fan.push_back(Vec2( 0.0,  0.0));
        fan.push_back(Vec2( 1.0,  1.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_APEX_COINCIDENT_VERTEX);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test127) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(-1.0,  1.0));
        fan.push_back(Vec2( 0.0, -1.0));
        fan.push_back(Vec2( 1.0,  1.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_APEX_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test128) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(-5.0,  1.0));
        fan.push_back(Vec2(-2.0, -3.0));
        fan.push_back(Vec2( 2.0, -3.0));
        fan.push_back(Vec2( 5.0,  1.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_APEX_INSIDE_VERTICES);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test129) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(-5.0,  1.0));
        fan.push_back(Vec2(-5.0, -3.0));
        fan.push_back(Vec2( 5.0, -3.0));
        fan.push_back(Vec2( 5.0,  1.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_APEX_ONE_EDGE_CUT_ACROSS);

    }
}


/**  @brief parseClassifiersPointsAgainstWedge() 
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test130) {

    for (long i = 0; i < 100; i++) {

        // TEST PARAM BEGIN
        vector<Vec2> fan;
        fan.push_back(Vec2(-5.0,  2.0));
        fan.push_back(Vec2(-5.0,  0.2));
        fan.push_back(Vec2( 5.0, -0.2));
        fan.push_back(Vec2( 5.0,  2.0));

        Vec2 apex(0.0, 0.0);
        Vec2 side1(1.0, -1.0);
        Vec2 side2(-1.0,-1.0);
        // TEST PARAM END

        vector<
           enum 
           Makena::ContactUpdaterFurtherChecker::
           WEDGE_CLASSIFIER
        > predMain, predAux;

        ConvexRigidBody body1(1);
        ConvexRigidBody body2(2);
        ContactPairInfo info;

        auto randAng = 2.0 * M_PI * rand100()/100.0;
        auto randV = randVec2D100()*0.05;

        roatateTestPatternrandomly(fan, apex, side1, side2, randAng, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);

        auto res_01 = upd_01.classifyPointsAgainstWedge(
            fan, apex, side1, side1.perp(), side2, (side2.perp())*-1.0, predMain, predAux);
        EXPECT_EQ(res_01, true);
        //showResults(predMain, predAux);
        auto res_02 = upd_01.parseClassifiersPointsAgainstWedge(predMain, predAux);
        //showResult(res_02);
        EXPECT_EQ(res_02, Makena::ContactUpdaterFurtherChecker::WC_APEX_ONE_EDGE_TOUCHING);

    }
}



/**  @brief checkFeatures_VERTEX_EDGE() no penetration
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test131) {

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
        activePoints_1.push_back(Vec3(-2.0,  2.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  2.0));

        Vec3       anchor_1(0.0,  2.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -10.0));
        polygon_2.push_back(Vec3( 1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  10.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -10.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -10.0));

        Vec3       anchor_2(0.0,  0.0, -10.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis(-1.0, 0.0, 0.0);
        double rotationRad = 0.0;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE,   expectedPoints_1_1,
            ContactPairInfo::FT_VERTEX, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE() completely inside
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test132) {

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
        activePoints_1.push_back(Vec3(-2.0,  2.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  2.0));

        Vec3       anchor_1(0.0,  2.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -10.0));
        polygon_2.push_back(Vec3( 1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  10.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -10.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -10.0));

        Vec3       anchor_2(0.0,  0.0, -10.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( 1.0, 0.0, 0.0);
        double rotationRad = M_PI * 0.75;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE,   expectedPoints_1_1,
            ContactPairInfo::FT_VERTEX, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE() completely inside boundary case.
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test133) {

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
        activePoints_1.push_back(Vec3(-2.0,  2.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  2.0));

        Vec3       anchor_1(0.0,  2.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -10.0));
        polygon_2.push_back(Vec3( 1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  10.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -10.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -10.0));

        Vec3       anchor_2(0.0,  0.0, -10.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( 1.0, 0.0, 0.0);
        double rotationRad = M_PI * 0.5 + atan(1.0/20.0);

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE,   expectedPoints_1_1,
            ContactPairInfo::FT_VERTEX, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE() completely inside boundary case.
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test134) {

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
        activePoints_1.push_back(Vec3(-2.0,  2.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  2.0));

        Vec3       anchor_1(0.0,  2.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -10.0));
        polygon_2.push_back(Vec3( 1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  10.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -10.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -10.0));

        Vec3       anchor_2(0.0,  0.0, -10.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis(-1.0, 0.0, 0.0);
        double rotationRad = M_PI + atan(1.0/20.0);

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE,   expectedPoints_1_1,
            ContactPairInfo::FT_VERTEX, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, true);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE() face-face
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test135) {

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
        activePoints_1.push_back(Vec3(-2.0,  2.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, -2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, -2.0));

        Vec3       anchor_1(0.0,  2.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -10.0));
        polygon_2.push_back(Vec3( 1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  10.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -10.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -10.0));
        expectedPoints_2_1.push_back(Vec3( 1.0,  1.0,  10.0));
        expectedPoints_2_1.push_back(Vec3(-1.0,  1.0,  10.0));

        Vec3       anchor_2(0.0,  0.0, -10.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis(-1.0, 0.0, 0.0);
        double rotationRad = M_PI;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE()  face-face
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test136) {

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
        activePoints_1.push_back(Vec3(-2.0,  2.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  2.0));

        Vec3       anchor_1(0.0,  2.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -10.0));
        polygon_2.push_back(Vec3( 1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  10.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -10.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -10.0));
        expectedPoints_2_1.push_back(Vec3( 1.0, -1.0,  10.0));
        expectedPoints_2_1.push_back(Vec3(-1.0, -1.0,  10.0));

        Vec3       anchor_2(0.0,  0.0, -10.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( 1.0, 0.0, 0.0);
        double rotationRad = M_PI * 0.5;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE() face-face
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test137) {

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
        activePoints_1.push_back(Vec3(-2.0,  2.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, -2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, -2.0));

        Vec3       anchor_1(0.0,  2.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -10.0));
        polygon_2.push_back(Vec3( 1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  10.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -10.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -10.0));
        expectedPoints_2_1.push_back(Vec3( 1.0,  1.0,  10.0));
        expectedPoints_2_1.push_back(Vec3(-1.0,  1.0,  10.0));

        Vec3       anchor_2(0.0,  0.0, -10.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis(-1.0, 0.0, 0.0);
        double rotationRad = M_PI - atan(1.0/20.0);

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE()  face-face
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test138) {

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
        activePoints_1.push_back(Vec3(-2.0,  2.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  2.0));

        Vec3       anchor_1(0.0,  2.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -10.0));
        polygon_2.push_back(Vec3( 1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  10.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  10.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  10.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -10.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -10.0));
        expectedPoints_2_1.push_back(Vec3( 1.0, -1.0,  10.0));
        expectedPoints_2_1.push_back(Vec3(-1.0, -1.0,  10.0));

        Vec3       anchor_2(0.0,  0.0, -10.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( 1.0, 0.0, 0.0);
        double rotationRad = M_PI * 0.5 - atan(1.0/20.0);

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);
        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE()  face-edge
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test139) {

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
        activePoints_1.push_back(Vec3(-2.0,  2.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0,  2.0));

        Vec3       anchor_1(0.0,  2.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -10.0));
        polygon_2.push_back(Vec3( 1.0,  0.0,  10.0));
        polygon_2.push_back(Vec3(-1.0,  0.0,  10.0));
        polygon_2.push_back(Vec3( 0.0, -1.0,  10.0));
        polygon_2.push_back(Vec3( 0.0,  1.0,  10.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -10.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -10.0));
        expectedPoints_2_1.push_back(Vec3( 0.0, -1.0,  10.0));

        Vec3       anchor_2(0.0,  0.0, -10.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( 1.0, 0.0, 0.0);
        double rotationRad = M_PI * 0.5 - atan(1.0/20.0);

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE()  face-edge
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test140) {

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
        activePoints_1.push_back(Vec3(-2.0,  2.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, -2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, -2.0));

        Vec3       anchor_1(0.0,  2.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -10.0));
        polygon_2.push_back(Vec3( 1.0,  0.0,  10.0));
        polygon_2.push_back(Vec3(-1.0,  0.0,  10.0));
        polygon_2.push_back(Vec3( 0.0, -1.0,  10.0));
        polygon_2.push_back(Vec3( 0.0,  1.0,  10.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -10.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -10.0));
        expectedPoints_2_1.push_back(Vec3( 0.0,  1.0,  10.0));

        Vec3       anchor_2(0.0,  0.0, -10.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis(-1.0, 0.0, 0.0);
        double rotationRad = M_PI - atan(1.0/20.0);

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE()  face-edge
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test141) {

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
        activePoints_1.push_back(Vec3(-2.0,  2.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  2.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0,  2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, -2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, -2.0));

        Vec3       anchor_1(0.0,  2.0,  2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -10.0));
        polygon_2.push_back(Vec3( 1.0,  0.0,  10.0));
        polygon_2.push_back(Vec3(-1.0,  0.0,  10.0));
        polygon_2.push_back(Vec3( 0.0, -1.0,  10.0));
        polygon_2.push_back(Vec3( 0.0,  1.0,  10.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -10.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -10.0));
        expectedPoints_2_1.push_back(Vec3( 0.0,  1.0,  10.0));

        Vec3       anchor_2(0.0,  0.0, -10.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis(-1.0, 0.0, 0.0);
        double rotationRad = M_PI;

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
                   
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        //EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE()  edge-edge coincident
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test142) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3(-2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0, 2.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 1.0,  0.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  0.0,  1.0));
        polygon_2.push_back(Vec3( 0.0, -1.0,  1.0));
        polygon_2.push_back(Vec3( 0.0,  1.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 1.0,  0.0,  1.0));

        Vec3       anchor_2(0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( 0.0, 1.0, 0.0);
        double rotationRad = atan(2.0/1.0);

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


        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );


        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE()  edge-edge coincident overshoot
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test143) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3(-2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, 0.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, 0.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3(-5.0,  0.0,  0.0));
        polygon_2.push_back(Vec3( 0.0, -1.0,  1.0));
        polygon_2.push_back(Vec3( 0.0,  0.0,  0.0));
        polygon_2.push_back(Vec3( 0.0,  2.0, -1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3(-5.0,  0.0, 0.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3(-5.0,  0.0,  0.0));
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0,  0.0));
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0, -1.0));

        Vec3       anchor_2(-5.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( 0.0, 1.0, 0.0);
        double rotationRad = 0.1;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE()  edge-edge coincident overshoot
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test144) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3(-2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0, -2.0, 0.0));
        expectedPoints_1_1.push_back(Vec3( 2.0, -2.0, 0.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3(-5.0,  0.0,  0.0));
        polygon_2.push_back(Vec3( 0.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 0.0,  0.0,  0.0));
        polygon_2.push_back(Vec3( 0.0, -2.0, -1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3(-5.0,  0.0, 0.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3(-5.0,  0.0,  0.0));
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0,  0.0));
        expectedPoints_2_1.push_back(Vec3( 0.0, -2.0, -1.0));

        Vec3       anchor_2(-5.0,  0.0, 0.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( 0.0, 1.0, 0.0);
        double rotationRad = 0.1;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE()  edge-edge coincident overshoot
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test145) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3(-2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0, 2.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 1.0,  1.0,  1.0));
        expectedPoints_2_1.push_back(Vec3( 1.0, -1.0,  1.0));

        Vec3       anchor_2( 0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( 0.0, 1.0, 0.0);
        double rotationRad = atan(2.0/1.0);

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_EDGE()  edge-edge coincident overshoot
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test146) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 2.0,  0.0,  2.0));
        polygon_1.push_back(Vec3(-2.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3(-2.0,  0.0,  2.0));
        activePoints_1.push_back(Vec3( 2.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 2.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  0.0, 2.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -1.0));
        polygon_2.push_back(Vec3( 1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0,  1.0,  1.0));
        polygon_2.push_back(Vec3( 1.0, -1.0,  1.0));
        polygon_2.push_back(Vec3(-1.0, -1.0,  1.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -1.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -1.0));
        expectedPoints_2_1.push_back(Vec3( 1.0,  1.0,  1.0));
        expectedPoints_2_1.push_back(Vec3( 1.0, -1.0,  1.0));

        Vec3       anchor_2( 0.0,  0.0, -1.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( 0.0, 1.0, 0.0);
        double rotationRad = atan(2.0/1.0) * 1.1;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}



/**  @brief checkFeatures_VERTEX_VERTEX()  no penetration
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test147) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0, 2.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -2.0));

        Vec3       anchor_2( 0.0,  0.0, -2.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( -1.0, 0.0, 0.0);
        double rotationRad = 0.0;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_VERTEX, expectedPoints_1_1,
            ContactPairInfo::FT_VERTEX, expectedPoints_2_1 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_VERTEX()  touching face-face
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test148) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, 0.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2( 0.0,  0.0, -2.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( -1.0, 0.0, 0.0);
        double rotationRad = atan(2.0/2.0)*2.0;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_VERTEX()  touching face-face overshoot
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test149) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);


        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, 0.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2( 0.0,  0.0, -2.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( -1.0, 0.0, 0.0);
        double rotationRad = atan(2.0/2.0)*2.0 * 1.1;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_VERTEX()  touching edge-edge
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test150) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);


        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, 0.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));

        Vec3       anchor_2( 0.0,  0.0, -2.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( -1.0, 1.0, 0.0);
        double rotationRad = atan(2.0/sqrt(8.0))*2.0 * 1.0;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_EDGE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ(checkResult(body1, body2, info, exp1, swap12), true);
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_VERTEX()  touching edge-edge
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test151) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, 0.0));

        vector<Vec3> expectedPoints_1_2;
        expectedPoints_1_2.push_back(Vec3( 0.0,  0.0, 2.0));
        expectedPoints_1_2.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_1_2.push_back(Vec3( 2.0, -2.0, 0.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3( 2.0, -2.0,  0.0));

        vector<Vec3> expectedPoints_2_2;
        expectedPoints_2_2.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2_2.push_back(Vec3( 2.0,  2.0,  0.0));
        expectedPoints_2_2.push_back(Vec3(-2.0,  2.0,  0.0));

        Vec3       anchor_2( 0.0,  0.0, -2.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis( -1.0, 1.0, 0.0);
        double rotationRad = atan(2.0/sqrt(8.0))*2.0 * 1.1;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

        ExpStructContactUpdaterFurtherChecker exp2(
            ContactPairInfo::FT_FACE, expectedPoints_1_2,
            ContactPairInfo::FT_FACE, expectedPoints_2_2 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ( checkResult(body1, body2, info, exp1, swap12)||
                   checkResult(body1, body2, info, exp2, swap12),  true );
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_VERTEX()  touching edge-edge
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test152) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, 0.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0));
        polygon_2.push_back(Vec3( 2.0,  0.0,  0.0));
        polygon_2.push_back(Vec3( 0.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  0.0,  0.0));
        polygon_2.push_back(Vec3( 0.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 0.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  0.0,  2.0));
        polygon_2.push_back(Vec3( 0.0, -2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0,  0.0));

        Vec3       anchor_2( 0.0,  0.0, -2.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis(-1.0, 0.0, 0.0);
        double rotationRad = atan(2.0/2.0)*2.0;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_EDGE, expectedPoints_2_1 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ( checkResult(body1, body2, info, exp1, swap12),  true );
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_VERTEX()  touching edge-edge
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test153) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, 0.0));

        vector<Vec3> expectedPoints_1_2;
        expectedPoints_1_2.push_back(Vec3( 0.0,  0.0, 2.0));
        expectedPoints_1_2.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_1_2.push_back(Vec3(-2.0,  2.0, 0.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0));
        polygon_2.push_back(Vec3( 2.0,  0.0,  0.0));
        polygon_2.push_back(Vec3( 0.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  0.0,  0.0));
        polygon_2.push_back(Vec3( 0.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 0.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  0.0,  2.0));
        polygon_2.push_back(Vec3( 0.0, -2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0,  0.0));

        vector<Vec3> expectedPoints_2_2;
        expectedPoints_2_2.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2_2.push_back(Vec3( 0.0,  2.0,  0.0));
        expectedPoints_2_2.push_back(Vec3(-2.0,  0.0,  0.0));

        Vec3       anchor_2( 0.0,  0.0, -2.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis(-1.0, 0.0, 0.0);
        double rotationRad = atan(2.0/2.0)*2.0*1.1;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

        ExpStructContactUpdaterFurtherChecker exp2(
            ContactPairInfo::FT_FACE, expectedPoints_1_2,
            ContactPairInfo::FT_FACE, expectedPoints_2_2 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ( checkResult(body1, body2, info, exp1, swap12)||
                   checkResult(body1, body2, info, exp2, swap12),  true );
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_VERTEX()  complex intersection
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test154) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);


        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, 0.0));

        vector<Vec3> expectedPoints_1_2;
        expectedPoints_1_2.push_back(Vec3( 0.0,  0.0, 2.0));
        expectedPoints_1_2.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_1_2.push_back(Vec3(-2.0,  2.0, 0.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0));
        polygon_2.push_back(Vec3( 2.0,  0.0,  0.0));
        polygon_2.push_back(Vec3( 0.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  0.0,  0.0));
        polygon_2.push_back(Vec3( 0.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 0.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0,  0.0,  2.0));
        polygon_2.push_back(Vec3( 0.0, -2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  0.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2_1.push_back(Vec3( 0.0,  2.0,  0.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  0.0,  0.0));

        vector<Vec3> expectedPoints_2_2;
        expectedPoints_2_2.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2_2.push_back(Vec3( 0.0,  2.0,  0.0));
        expectedPoints_2_2.push_back(Vec3(-2.0,  0.0,  0.0));

        Vec3       anchor_2( 0.0,  0.0, -2.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis(-1.0, 0.0, 0.0);
        double rotationRad = atan(2.0/2.0)*2.0*1.1;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

        ExpStructContactUpdaterFurtherChecker exp2(
            ContactPairInfo::FT_FACE, expectedPoints_1_2,
            ContactPairInfo::FT_FACE, expectedPoints_2_2 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ( checkResult(body1, body2, info, exp1, swap12)||
                   checkResult(body1, body2, info, exp2, swap12),  true );
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_VERTEX()
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test155) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 3.0,  0.0,  0.0));
        polygon_1.push_back(Vec3(-3.0,  0.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 3.0,  0.0, -2.0));
        polygon_1.push_back(Vec3(-3.0,  0.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, 0.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 3.0,  0.0,  0.0));
        polygon_2.push_back(Vec3(-3.0,  0.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 3.0,  0.0,  2.0));
        polygon_2.push_back(Vec3(-3.0,  0.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0, 0.0));

        Vec3       anchor_2( 0.0,  0.0, -2.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.0;

        Vec3 tiltAxis(-1.0, 0.0, 0.0);
        double rotationRad = atan(2.0/2.0)*2.0;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ( checkResult(body1, body2, info, exp1, swap12),  true );
        //showResult(body1, body2, info);
    }
}


/**  @brief checkFeatures_VERTEX_VERTEX()
 */
TEST_F(ContactUpdaterFurtherCheckerTest, Test156) {

    for (long i = 0; i < 20; i++) {

        bool swap12 = ((i%2)==0);

        vector<Vec3> polygon_1;
        polygon_1.push_back(Vec3( 0.0,  0.0,  2.0));
        polygon_1.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_1.push_back(Vec3( 3.0,  0.0,  0.0));
        polygon_1.push_back(Vec3(-3.0,  0.0,  0.0));
        polygon_1.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_1.push_back(Vec3( 2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0,  2.0, -2.0));
        polygon_1.push_back(Vec3( 2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3(-2.0, -2.0, -2.0));
        polygon_1.push_back(Vec3( 3.0,  0.0, -2.0));
        polygon_1.push_back(Vec3(-3.0,  0.0, -2.0));

        vector<Vec3> activePoints_1;
        activePoints_1.push_back(Vec3( 0.0,  0.0,  2.0));

        vector<Vec3> expectedPoints_1_1;
        expectedPoints_1_1.push_back(Vec3( 0.0,  0.0, 2.0));
        expectedPoints_1_1.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_1_1.push_back(Vec3(-2.0,  2.0, 0.0));

        Vec3       anchor_1(0.0,  0.0, 2.0);

        Vec3       CoM1  (0.0, 0.0, 0.0);
        Quaternion q1    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin1 (0.0, 0.0, 0.0);
        Vec3       Vang1 (0.0, 0.0, 0.0);

        vector<Vec3> polygon_2;
        polygon_2.push_back(Vec3( 0.0,  0.0, -2.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  0.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  0.0));
        polygon_2.push_back(Vec3( 3.0,  0.0,  0.0));
        polygon_2.push_back(Vec3(-3.0,  0.0,  0.0));
        polygon_2.push_back(Vec3( 2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0,  2.0,  2.0));
        polygon_2.push_back(Vec3( 2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3(-2.0, -2.0,  2.0));
        polygon_2.push_back(Vec3( 3.0,  0.0,  2.0));
        polygon_2.push_back(Vec3(-3.0,  0.0,  2.0));

        vector<Vec3> activePoints_2;
        activePoints_2.push_back(Vec3( 0.0,  0.0, -2.0));

        vector<Vec3> expectedPoints_2_1;
        expectedPoints_2_1.push_back(Vec3( 0.0,  0.0, -2.0));
        expectedPoints_2_1.push_back(Vec3( 2.0,  2.0, 0.0));
        expectedPoints_2_1.push_back(Vec3(-2.0,  2.0, 0.0));

        Vec3       anchor_2( 0.0,  0.0, -2.0);

        Vec3       CoM2  (0.0, 0.0, 0.0);
        Quaternion q2    (1.0, 0.0, 0.0, 0.0);
        Vec3       Vlin2 (0.0, 0.0, 0.0);
        Vec3       Vang2 (0.0, 0.0, 0.0);
        double     YawAngle = 0.2*rand100()/100.0 - 0.1;

        Vec3 tiltAxis(-1.0, 0.0, 0.0);
        double rotationRad = atan(2.0/2.0)*2.0*1.5;

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

        ExpStructContactUpdaterFurtherChecker exp1(
            ContactPairInfo::FT_FACE, expectedPoints_1_1,
            ContactPairInfo::FT_FACE, expectedPoints_2_1 );

        updateGeomConfig(CoM1, q1, Vlin1, Vang1, swap12?body2:body1);
        updateGeomConfig(CoM2, q2, Vlin2, Vang2, swap12?body1:body2);

        auto randQ =  randomRotQ();           
        auto randV = randVec3D100()*0.05;
        rotateTestPatternRandomly(body1, body2, info, randQ, randV);

        ContactUpdaterFurtherChecker
            upd_01(body1,body2,info,EPSILON_SQUARED,EPSILON_SQUARED,std::cerr);
        auto res_01 = upd_01.checkFeatures();
        auto randQmat = randQ.rotationMatrix();
        randQmat.transposeInPlace();
        auto N = randQmat * (info.mContactNormal1To2);
        EXPECT_EQ(res_01, false);

        EXPECT_EQ( checkResult(body1, body2, info, exp1, swap12),  true );
        //showResult(body1, body2, info);
    }
}


} // namespace Makena

#include "gtest/gtest.h"
#include <vector>
#include "primitives.hpp"
#include "manifold.hpp"
#include "quaternion.hpp"
#include "intersection_finder.hpp"

#include <iostream>
#include <fstream>
#include <sstream>


namespace Makena {



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


static VertexIt findVertex(Manifold& m, const Vec3& p)
{
    auto vPair = m.vertices();
    for (auto vit = vPair.first ; vit != vPair.second; vit++) {
        if ((*vit)->pLCS() == p) {
            return vit;
        }
    }
    return vPair.second;
}


static EdgeIt findEdge(Manifold& m, const Vec3& p1, const Vec3& p2)
{
    auto ePair = m.edges();
    for (auto eit = ePair.first ; eit != ePair.second; eit++) {
        auto he = (*eit)->he1();
        auto src = (*he)->src();
        auto dst = (*he)->dst();
        auto sp  = (*src)->pLCS();
        auto dp  = (*dst)->pLCS();
        if ( (p1 == sp && p2 == dp) ||
             (p1 == dp && p2 == sp)   ) {
            return eit;
        }
    }
    return ePair.second;
}


static FaceIt findFace(Manifold& m, const Vec3& n)
{
    auto fPair = m.faces();
    for (auto fit = fPair.first ; fit != fPair.second; fit++) {
        if ((*fit)->nLCS() == n) {
            return fit;
        }
    }
    return fPair.second;
}


bool static compareAttributes(
    bool                            checkPoint,
    IntersectionFinder::Attributes &a1,
    IntersectionFinder::Attributes &a2
) {
    if (checkPoint) {
        if (a1.mP != a2.mP) {
            return false;
        }
    }

    if (a1.mPred != a2.mPred) {
        return false;
    }

    switch (a1.mPred) {


        case IF_VERTEX_VERTEX:
          if (a1.mVit1 == a2.mVit1 && a1.mVit2 == a2.mVit2) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_VERTEX_EDGE:
          if (a1.mVit1 == a2.mVit1 && a1.mEit2 == a2.mEit2) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_VERTEX_FACE:
          if (a1.mVit1 == a2.mVit1 && a1.mFit2 == a2.mFit2) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_VERTEX_INTERIOR:
          if (a1.mVit1 == a2.mVit1) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_EDGE_VERTEX:
          if (a1.mEit1 == a2.mEit1 && a1.mVit2 == a2.mVit2) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_EDGE_EDGE:
          if (a1.mEit1 == a2.mEit1 && a1.mEit2 == a2.mEit2) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_EDGE_FACE:
          if (a1.mEit1 == a2.mEit1 && a1.mFit2 == a2.mFit2) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_EDGE_INTERIOR:
          if (a1.mEit1 == a2.mEit1) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_FACE_VERTEX:
          if (a1.mFit1 == a2.mFit1 && a1.mVit2 == a2.mVit2) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_FACE_EDGE:
          if (a1.mFit1 == a2.mFit1 && a1.mEit2 == a2.mEit2) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_FACE_FACE:
          if (a1.mFit1 == a2.mFit1 && a1.mFit2 == a2.mFit2) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_FACE_INTERIOR:
          if (a1.mFit1 == a2.mFit1) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_INTERIOR_VERTEX:
          if (a1.mVit2 == a2.mVit2) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_INTERIOR_EDGE:
          if (a1.mEit2 == a2.mEit2) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_INTERIOR_FACE:
          if (a1.mFit2 == a2.mFit2) {
              return true;
          }
          else {
              return false;
          }
          break;

        case IF_INTERIOR_INTERIOR:
          return true;
          break; 

        default:
          return false;
          break; 
    }
    return false;
}


bool static compareAttributes(
    vector<IntersectionFinder::Attributes>& v1,
    vector<IntersectionFinder::Attributes>& v2
) {
    if (v1.size() != v2.size()) {
        return false;
    }
    vector<long> cnts;
    for (long i = 0; i < v1.size(); i++) {
        cnts.push_back(0);
    }
    for (long i = 0; i < v1.size() ; i++) {
        auto & a1 = v1[i];
        for (long j = 0; j < v1.size() ; j++) {
            auto & a2 = v2[j];
            if (compareAttributes(true, a1, a2)) {
                cnts[i] = cnts[i] + 1;
            }
        }
    }

    for (long i = 0; i < v1.size(); i++) {
        if (cnts[i] != 1) {
            return false;
        }
    }
    return true;
}


static bool checkResult(
    IntersectionFinder&             finder,
    VertexIt                        vit,
    IntersectionFinder::Attributes& attrExp
) {
    if ((*vit)->pLCS() != attrExp.mP) {
        return false;
    }

    IntersectionFinder::Attributes& attrActual=finder.vertexAttributes3D(vit);
    bool res = compareAttributes(true, attrExp, attrActual);
    return res;

}


static bool checkResult(
    bool                            checkPoint,
    IntersectionFinder::Attributes& attr1,
    IntersectionFinder::Attributes& attr2
) {
    bool res = compareAttributes(checkPoint, attr1, attr2);
    return res;
}


static bool checkResult(
    IntersectionFinder&             finder,
    const Vec3&                     ep1,
    const Vec3&                     ep2,
    EdgeIt                          eit,
    IntersectionFinder::Attributes& attrExp
) {

    auto he = (*eit)->he1();
    auto src = (*he)->src();
    auto dst = (*he)->dst();
    auto p1  = (*src)->pLCS(); 
    auto p2  = (*dst)->pLCS(); 

    if ((p1 == ep1 && p2 == ep2) || (p1 == ep2 && p2 == ep1)) {
        ;
    }
    else {
        return false;
    }
    IntersectionFinder::Attributes& attrActual=finder.edgeAttributes3D(eit);
    bool res = compareAttributes(false, attrExp, attrActual);

    return res;
}


static bool checkResult(
    IntersectionFinder&             finder,
    FaceIt                          fit,
    IntersectionFinder::Attributes& attrExp
) {
    IntersectionFinder::Attributes& attrActual=finder.faceAttributes3D(fit);
    bool res = compareAttributes(false, attrExp, attrActual);
    return res;
}


static bool checkResults3D(
    IntersectionFinder&                     finder,
    vector<IntersectionFinder::Attributes>& vatts,
    vector<pair<Vec3, Vec3> >&              eps,
    vector<IntersectionFinder::Attributes>& eatts,
    vector<IntersectionFinder::Attributes>& fatts
) {
    if (finder.dimension()!=3) {
        return false;
    }

    Manifold& mani = finder.hull3D();

    auto vPair = mani.vertices();
    vector<long> Vcnts;
    for (long i = 0; i < vatts.size(); i++) {
        Vcnts.push_back(0); 
    }
    long numVertices = 0;
    for (auto vit = vPair.first ; vit != vPair.second; vit++) {
        for (long i = 0; i < vatts.size(); i++) {
            if (checkResult(finder, vit, vatts[i])) {
                Vcnts[i] = Vcnts[i] + 1;
            }
        }
        numVertices++;
    }
    if (numVertices != vatts.size()) {
        return false;
    }
    for (long i = 0; i < vatts.size(); i++) { 
         if(Vcnts[i] != 1) {
             return false;
         }
    }

    vector<long> Ecnts;
    for (long i = 0; i < eatts.size(); i++) {
        Ecnts.push_back(0);
    }
    long numEdges = 0;
    auto ePair = mani.edges();
    for (long i = 0; i < eatts.size(); i++) {

        for (auto eit = ePair.first ; eit != ePair.second; eit++) {

            if (checkResult(finder,eps[i].first,eps[i].second,eit,eatts[i])){
                Ecnts[i] = Ecnts[i] + 1;
            }
        }
        numEdges++;
    }
    if (numEdges != eatts.size()) {
        return false;
    }
    for (long i = 0; i < eatts.size(); i++) { 
         if(Ecnts[i] != 1) {
             return false;
         }
    }

    vector<long> Fcnts;
    for (long i = 0; i < fatts.size(); i++) {
        Fcnts.push_back(0);
    }
    long numFaces = 0;
    auto fPair = mani.faces();
    for (long i = 0; i < fatts.size(); i++) {
        for (auto fit = fPair.first ; fit != fPair.second; fit++) {
            if (checkResult(finder, fit, fatts[i])) {
                Fcnts[i] = Fcnts[i] + 1;
            }
        }
        numFaces++;
    }
    if (numFaces != fatts.size()) {
        return false;
    }
    for (long i = 0; i < fatts.size(); i++) { 
         if(Fcnts[i] != 1) {
             return false;
         }
    }

    return true;

}


static bool checkResults2D(
    IntersectionFinder&                     finder,
    vector<IntersectionFinder::Attributes>& vatts,
    vector<pair<Vec3, Vec3> >&              eps,
    vector<IntersectionFinder::Attributes>& eatts,
    IntersectionFinder::Attributes&         fatts
) {

    if (finder.dimension()!=2) {
        return false;
    }
    if (finder.size2D() != vatts.size() || finder.size2D() != eatts.size()) {
        return false;        
    }

    bool found1 = false;
    for (long offset = 0; offset < finder.size2D(); offset++) {

        bool found2 = true;
        for (long i1 = 0; i1 < finder.size2D(); i1++) {
            long i2 = (i1 + offset)%finder.size2D();
            if (!checkResult(true, vatts[i1], finder.vertexAttributes2D(i2))){
                found2 = false;
                break;
            }
        }
        if (found2){
            found1 = true;
            break;
        }

    }
    if (!found1) {
        for (long offset = 0; offset < finder.size2D(); offset++) {

            bool found2 = true;
            for (long i1 = 0; i1 < finder.size2D(); i1++) {

                long i2 = (finder.size2D() - 1) - 
                          ((i1 + offset)%finder.size2D());
                if (!checkResult(
                           true, vatts[i1], finder.vertexAttributes2D(i2))){
                    found2 = false;
                    break;
                }
            }
            if (found2){
                found1 = true;
                break;
            }
        }

    }

    if (!found1) {
        return false;
    }

    found1 = false;
    for (long offset = 0; offset < finder.size2D(); offset++) {
        bool found2 = true;
        for (long i1 = 0; i1 < finder.size2D(); i1++) {

            long i2      = (i1 + offset)%finder.size2D();
            long i2next  = (i2 + 1)%finder.size2D();
            if (!checkResult(false, eatts[i1], finder.edgeAttributes2D(i2))){
                found2 = false;
                break;
            }
            auto& p1 = finder.vertexAttributes2D(i2).mP;
            auto& p2 = finder.vertexAttributes2D(i2next).mP;
            auto& ep1 = eps[i1].first;
            auto& ep2 = eps[i1].second;
            if ((p1 == ep1 && p2 == ep2 )||(p1 == ep2 && p2 == ep1 )) {
                ;
            }
            else {
                found2 = false;
                break;
            }
        }
        if (found2){
            found1 = true;
            break;
        }

    }
    if (!found1) {
        for (long offset = 0; offset < finder.size2D(); offset++) {
            bool found2 = true;
            for (long i1 = 0; i1 < finder.size2D(); i1++) {

                long i2 = (finder.size2D() - 1) - 
                          ((i1 + offset)%finder.size2D());
                long i2next  = (i2 + 1)%finder.size2D();
                if (!checkResult(
                           false, eatts[i1], finder.edgeAttributes2D(i2))){
                    found2 = false;
                    break;
                }

                auto& p1 = finder.vertexAttributes2D(i2).mP;
                auto& p2 = finder.vertexAttributes2D(i2next).mP;
                auto& ep1 = eps[i1].first;
                auto& ep2 = eps[i1].second;
                if ((p1 == ep1 && p2 == ep2 )||(p1 == ep2 && p2 == ep1 )) {
                    ;
                }
                else {
                    found2 = false;
                    break;
                }
            }
            if (found2){
                found1 = true;
                break;
            }
        }

    }

    if (!found1) {
        return false;
    }

    auto res = checkResult(false, fatts, finder.faceAttributes2D());
    return true;
}


static bool checkResults1D(
    IntersectionFinder&                     finder,
    vector<IntersectionFinder::Attributes>& vatts,
    pair<Vec3, Vec3>&                       ep,
    IntersectionFinder::Attributes&         eatts
) {

    if (finder.dimension()!=1) {
        return false;
    }

    auto& vExp1 = finder.vertexAttributes1D(0);
    auto& vExp2 = finder.vertexAttributes1D(1);

    if ( ( checkResult(true, vatts[0], vExp1) &&
           checkResult(true, vatts[1], vExp2)    ) ||
         ( checkResult(true, vatts[0], vExp2) &&
           checkResult(true, vatts[1], vExp1)    )    ) {
        ;
    }
    else {
        return false;
    }

    if ( ( ep.first == vExp1.mP && ep.second == vExp2.mP ) ||
         ( ep.first == vExp2.mP && ep.second == vExp1.mP )    ) {
        ;
    }
    else {
        return false;
    }

    if (!checkResult(false, eatts, finder.edgeAttributes1D())) {
        return false;
    }

    return true;
}



static bool checkResults0D(
    IntersectionFinder&             finder,
    IntersectionFinder::Attributes& vatt
) {

    if (finder.dimension()!=0) {
        return false;
    }

    auto& vExp1 = finder.vertexAttributes0D();

    auto res = checkResult(true, vatt, vExp1);

    return res;
}


class IntersectionFinderTest : public ::testing::Test {

  protected:  

    IntersectionFinderTest(){;}
    virtual ~IntersectionFinderTest(){;}
    virtual void SetUp() {;};
    virtual void TearDown() {;};

};


/**  @brief prepareVertices();
 *          prepareEdgesAndFaces()
 */
TEST_F(IntersectionFinderTest, Test01) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Quaternion Q01(1.0, 0.0, 0.0, 0.0), Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0),   CoM02 (0.0, 0.0, 0.0);

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        Vec3 sepAxis(1.0, 0.0, 0.0);

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Quaternion q01     = randomRotQ();
        Mat3x3     qMat01  = q01.rotationMatrix();
        Vec3       trans01 = randVec3D100();

        Qmat01  =  qMat01 * Qmat01;
        Qmat02  =  qMat01 * Qmat02;
        CoM01   =  qMat01 * CoM01;
        CoM02   =  qMat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  qMat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();


        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            EXPECT_EQ((*vit)->IFisActive(), true);
            bool found = false;
            for(auto v2:finder01.mActiveVertices1){if(vit==v2){found=true;}}
            EXPECT_EQ(found, true);
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            EXPECT_EQ((*vit)->IFisActive(), true);
            bool found = false;
            for(auto v2:finder01.mActiveVertices2){if(vit==v2){found=true;}}
            EXPECT_EQ(found, true);
        }

        finder01.prepareVertices(sepAxis);


        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            EXPECT_EQ((*vit)->IFisActive(), true);
            bool found = false;
            for(auto v2:finder01.mActiveVertices1){if(vit==v2){found=true;}}
            EXPECT_EQ(found, true);
        }

        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            EXPECT_EQ((*vit)->IFisActive(), true);
            bool found = false;
            for(auto v2:finder01.mActiveVertices2){if(vit==v2){found=true;}}
            EXPECT_EQ(found, true);
        }

        finder01.prepareEdgesAndFaces();

        auto ePair01 = m01.edges();
        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges1){if(eit==e2){found=true;}}
            EXPECT_EQ(found, true);
        }

        auto ePair02 = m02.edges();
        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges2){if(eit==e2){found=true;}}
            EXPECT_EQ(found, true);
        }

        auto fPair01 = m01.faces();
        for (auto fit = fPair01.first; fit != fPair01.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces1){if(fit==f2){found=true;}}
            EXPECT_EQ(found, true);
        }

        auto fPair02 = m02.faces();
        for (auto fit = fPair02.first; fit != fPair02.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces2){if(fit==f2){found=true;}}
            EXPECT_EQ(found, true);
        }
    }
}


/**  @brief prepareVertices() #2; Separating
 *          prepareEdgesAndFaces()
 */
TEST_F(IntersectionFinderTest, Test02) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(10.0, rand100() - 50.0, rand100() - 50.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();

        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            EXPECT_EQ((*vit)->IFisActive(), true);
            bool found = false;
            for(auto v2:finder01.mActiveVertices1){if(vit==v2){found=true;}}
            EXPECT_EQ(found, true);
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            EXPECT_EQ((*vit)->IFisActive(), true);
            bool found = false;
            for(auto v2:finder01.mActiveVertices2){if(vit==v2){found=true;}}
            EXPECT_EQ(found, true);
        }

        finder01.mActiveVertices1.clear();
        finder01.mActiveVertices2.clear();

        finder01.prepareVertices(sepAxis);
        finder01.prepareEdgesAndFaces();

        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            EXPECT_EQ((*vit)->IFisActive(), false);
            bool found = false;
            for(auto v2:finder01.mActiveVertices1){if(vit==v2){found=true;}}
            EXPECT_EQ(found, false);
        }

        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            EXPECT_EQ((*vit)->IFisActive(), false);
            bool found = false;
            for(auto v2:finder01.mActiveVertices2){if(vit==v2){found=true;}}
            EXPECT_EQ(found, false);
        }

        auto ePair01 = m01.edges();
        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges1){if(eit==e2){found=true;}}
            EXPECT_EQ(found, false);
        }

        auto ePair02 = m02.edges();
        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges2){if(eit==e2){found=true;}}
            EXPECT_EQ(found, false);
        }

        auto fPair01 = m01.faces();
        for (auto fit = fPair01.first; fit != fPair01.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces1){if(fit==f2){found=true;}}
            EXPECT_EQ(found, false);
        }

        auto fPair02 = m02.faces();
        for (auto fit = fPair02.first; fit != fPair02.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces2){if(fit==f2){found=true;}}
            EXPECT_EQ(found, false);
        }
    }
    
}


/**  @brief prepareVertices() #3; Touching
 *          prepareEdgesAndFaces()
 */
TEST_F(IntersectionFinderTest, Test03) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, rand100() - 50.0, rand100() - 50.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices(sepAxis);
        finder01.prepareEdgesAndFaces();

        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()==1.0) {
                EXPECT_EQ((*vit)->IFisActive(), true);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices1){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, true);
            }
            else {
                EXPECT_EQ((*vit)->IFisActive(), false);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices1){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, false);
            }
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0) {
                EXPECT_EQ((*vit)->IFisActive(), true);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices2){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, true);
            }
            else {
                EXPECT_EQ((*vit)->IFisActive(), false);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices2){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, false);
            }
        }

        auto ePair01 = m01.edges();
        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges1){if(eit==e2){found=true;}}
            auto he = (*eit)->he1();
            auto p1 = (*((*he)->src()))->pLCS();
            auto p2 = (*((*he)->dst()))->pLCS();
            if (p1.x()==1.0 ||p2.x()==1.0) {
                EXPECT_EQ(found, true);
            }
            else {
                EXPECT_EQ(found, false);
            }
        }

        auto ePair02 = m02.edges();
        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges2){if(eit==e2){found=true;}}
            auto he = (*eit)->he1();
            auto p1 = (*((*he)->src()))->pLCS();
            auto p2 = (*((*he)->dst()))->pLCS();
            if (p1.x()==-1.0 ||p2.x()==-1.0) {
                EXPECT_EQ(found, true);
            }
            else {
                EXPECT_EQ(found, false);
            }
        }

        auto fPair01 = m01.faces();
        for (auto fit = fPair01.first; fit != fPair01.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces1){if(fit==f2){found=true;}}
            EXPECT_EQ(found, (*fit)->nLCS()!=Vec3(-1.0, 0.0, 0.0));
        }

        auto fPair02 = m02.faces();
        for (auto fit = fPair02.first; fit != fPair02.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces2){if(fit==f2){found=true;}}
            EXPECT_EQ(found, (*fit)->nLCS()!=Vec3( 1.0, 0.0, 0.0));
        }
    }

}


/**  @brief prepareVertices() #4; Overlapping
 *          prepareEdgesAndFaces()
 */
TEST_F(IntersectionFinderTest, Test04) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(1.8, rand100() - 50.0, rand100() - 50.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices(sepAxis);
        finder01.prepareEdgesAndFaces();

        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()==1.0) {
                EXPECT_EQ((*vit)->IFisActive(), true);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices1){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, true);
            }
            else {
                EXPECT_EQ((*vit)->IFisActive(), false);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices1){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, false);
            }
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0) {
                EXPECT_EQ((*vit)->IFisActive(), true);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices2){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, true);
            }
            else {
                EXPECT_EQ((*vit)->IFisActive(), false);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices2){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, false);
            }
        }

        auto ePair01 = m01.edges();
        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges1){if(eit==e2){found=true;}}
            auto he = (*eit)->he1();
            auto p1 = (*((*he)->src()))->pLCS();
            auto p2 = (*((*he)->dst()))->pLCS();
            if (p1.x()==1.0 ||p2.x()==1.0) {
                EXPECT_EQ(found, true);
            }
            else {
                EXPECT_EQ(found, false);
            }
        }

        auto ePair02 = m02.edges();
        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges2){if(eit==e2){found=true;}}
            auto he = (*eit)->he1();
            auto p1 = (*((*he)->src()))->pLCS();
            auto p2 = (*((*he)->dst()))->pLCS();
            if (p1.x()==-1.0 ||p2.x()==-1.0) {
                EXPECT_EQ(found, true);
            }
            else {
                EXPECT_EQ(found, false);
            }
        }

        auto fPair01 = m01.faces();
        for (auto fit = fPair01.first; fit != fPair01.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces1){if(fit==f2){found=true;}}
            EXPECT_EQ(found, (*fit)->nLCS()!=Vec3(-1.0, 0.0, 0.0));
        }

        auto fPair02 = m02.faces();
        for (auto fit = fPair02.first; fit != fPair02.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces2){if(fit==f2){found=true;}}
            EXPECT_EQ(found, (*fit)->nLCS()!=Vec3( 1.0, 0.0, 0.0));
        }

    }
}


/**  @brief prepareVertices() #5; Overlapping
 *          prepareEdgesAndFaces()
 */
TEST_F(IntersectionFinderTest, Test05a) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM( 1.0 + 1.0*cos(M_PI*0.2), 
                            rand100() - 50.0, 
                            rand100() - 50.0    );

        Quaternion Q01(Vec3(0.0,-1.0, 0.0), M_PI*0.0);
        Quaternion Q02(Vec3(0.0, 1.0, 0.0), M_PI*0.2);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices(sepAxis);
        finder01.prepareEdgesAndFaces();

        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()== 1.0) {
                EXPECT_EQ((*vit)->IFisActive(), true);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices1){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, true);
            }
            else {
                EXPECT_EQ((*vit)->IFisActive(), false);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices1){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, false);
            }
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0 && (*vit)->pLCS().z()==-1.0) {
                EXPECT_EQ((*vit)->IFisActive(), true);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices2){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, true);
            }
            else {
                EXPECT_EQ((*vit)->IFisActive(), false);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices2){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, false);
            }
        }

        auto ePair01 = m01.edges();
        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges1){if(eit==e2){found=true;}}
            auto he = (*eit)->he1();
            auto p1 = (*((*he)->src()))->pLCS();
            auto p2 = (*((*he)->dst()))->pLCS();
            if (p1.x()==1.0||p2.x()==1.0) {
                EXPECT_EQ(found, true);
            }
            else {
                EXPECT_EQ(found, false);
            }
        }

        auto ePair02 = m02.edges();
        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges2){if(eit==e2){found=true;}}
            auto he = (*eit)->he1();
            auto p1 = (*((*he)->src()))->pLCS();
            auto p2 = (*((*he)->dst()))->pLCS();
            if ((p1.x()==-1.0 && p1.z()==-1.0)||(p2.x()==-1.0&&p2.z()==-1.0)) {
                EXPECT_EQ(found, true);
            }
            else {
                EXPECT_EQ(found, false);
            }
        }

        auto fPair01 = m01.faces();
        for (auto fit = fPair01.first; fit != fPair01.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces1){if(fit==f2){found=true;}}
            EXPECT_EQ(found, (*fit)->nLCS()!=Vec3(-1.0, 0.0, 0.0));
        }

        auto fPair02 = m02.faces();
        for (auto fit = fPair02.first; fit != fPair02.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces2){if(fit==f2){found=true;}}
            EXPECT_EQ(found, (*fit)->nLCS()!=Vec3( 1.0, 0.0, 0.0)&&
                             (*fit)->nLCS()!=Vec3( 0.0, 0.0, 1.0)   );
        }

    }
}


/**  @brief prepareVertices() #5; Overlapping
 *          prepareEdgesAndFaces()
 */
TEST_F(IntersectionFinderTest, Test05b) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM( 1.0 + 1.0*cos(M_PI*0.2), 
                            rand100() - 50.0, 
                            rand100() - 50.0    );

        Quaternion Q01(Vec3(0.0,-1.0, 0.0), M_PI*0.2);
        Quaternion Q02(Vec3(0.0, 1.0, 0.0), M_PI*0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices(sepAxis);
        finder01.prepareEdgesAndFaces();

        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()== 1.0 && (*vit)->pLCS().z()== -1.0) {
                EXPECT_EQ((*vit)->IFisActive(), true);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices1){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, true);
            }
            else {
                EXPECT_EQ((*vit)->IFisActive(), false);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices1){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, false);
            }
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0) {
                EXPECT_EQ((*vit)->IFisActive(), true);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices2){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, true);
            }
            else {
                EXPECT_EQ((*vit)->IFisActive(), false);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices2){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, false);
            }
        }

        auto ePair01 = m01.edges();
        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges1){if(eit==e2){found=true;}}
            auto he = (*eit)->he1();
            auto p1 = (*((*he)->src()))->pLCS();
            auto p2 = (*((*he)->dst()))->pLCS();
            if ((p1.x()== 1.0 && p1.z()==-1.0)||(p2.x()== 1.0&&p2.z()==-1.0)) {
                EXPECT_EQ(found, true);
            }
            else {
                EXPECT_EQ(found, false);
            }
        }

        auto ePair02 = m02.edges();
        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges2){if(eit==e2){found=true;}}
            auto he = (*eit)->he1();
            auto p1 = (*((*he)->src()))->pLCS();
            auto p2 = (*((*he)->dst()))->pLCS();
            if (p1.x()==-1.0||p2.x()==-1.0) {
                EXPECT_EQ(found, true);
            }
            else {
                EXPECT_EQ(found, false);
            }
        }

        auto fPair01 = m01.faces();
        for (auto fit = fPair01.first; fit != fPair01.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces1){if(fit==f2){found=true;}}
            EXPECT_EQ(found, (*fit)->nLCS()!=Vec3(-1.0, 0.0, 0.0)&&
                             (*fit)->nLCS()!=Vec3( 0.0, 0.0, 1.0)   );
        }

        auto fPair02 = m02.faces();
        for (auto fit = fPair02.first; fit != fPair02.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces2){if(fit==f2){found=true;}}
            EXPECT_EQ(found, (*fit)->nLCS()!=Vec3( 1.0, 0.0, 0.0));
        }

    }
}


/**  @brief prepareVertices() #5; Overlapping
 *          prepareEdgesAndFaces()
 */
TEST_F(IntersectionFinderTest, Test06) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(1.0, rand100() - 50.0, rand100() - 50.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(Vec3(0.0, 1.0, 0.0), M_PI*0.2);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices(sepAxis);
        finder01.prepareEdgesAndFaces();

        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()== 1.0) {
                EXPECT_EQ((*vit)->IFisActive(), true);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices1){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, true);
            }
            else {
                EXPECT_EQ((*vit)->IFisActive(), false);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices1){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, false);
            }
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0) {
                EXPECT_EQ((*vit)->IFisActive(), true);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices2){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, true);
            }
            else {
                EXPECT_EQ((*vit)->IFisActive(), false);
                bool fd = false;
                for(auto v2:finder01.mActiveVertices2){if(vit==v2){fd=true;}}
                EXPECT_EQ(fd, false);
            }
        }

        auto ePair01 = m01.edges();
        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges1){if(eit==e2){found=true;}}
            auto he = (*eit)->he1();
            auto p1 = (*((*he)->src()))->pLCS();
            auto p2 = (*((*he)->dst()))->pLCS();
            if (p1.x()==1.0||p2.x()==1.0) {
                EXPECT_EQ(found, true);
            }
            else {
                EXPECT_EQ(found, false);
            }
        }

        auto ePair02 = m02.edges();
        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {
            bool found = false;
            for(auto e2:finder01.mActiveEdges2){if(eit==e2){found=true;}}
            auto he = (*eit)->he1();
            auto p1 = (*((*he)->src()))->pLCS();
            auto p2 = (*((*he)->dst()))->pLCS();
            if (p1.x()==-1.0||p2.x()==-1.0) {
                EXPECT_EQ(found, true);
            }
            else {
                EXPECT_EQ(found, false);
            }
        }

        auto fPair01 = m01.faces();
        for (auto fit = fPair01.first; fit != fPair01.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces1){if(fit==f2){found=true;}}
            EXPECT_EQ(found, (*fit)->nLCS()!=Vec3(-1.0, 0.0, 0.0));
        }

        auto fPair02 = m02.faces();
        for (auto fit = fPair02.first; fit != fPair02.second; fit++) {
            bool found = false;
            for(auto f2:finder01.mActiveFaces2){if(fit==f2){found=true;}}
            EXPECT_EQ(found, (*fit)->nLCS()!=Vec3( 1.0, 0.0, 0.0));
        }
    }
}


/**  @brief findIntersectionsEdgesOn1FacesOn2()
 *          findIntersectionsFacesOn1EdgesOn2()
 *          findRemainingInteriorVertices()
 */
TEST_F(IntersectionFinderTest, Test07) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(4.0, rand100() - 50.0, rand100() - 50.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Vec3       v01;
        while (v01.squaredNorm2()<=EPSILON_SQUARED){v01=randVec3D100();}
        
        Quaternion Q02(v01, M_PI*(rand100()/50.0 - 1.0));
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();

        EXPECT_EQ(finder01.mMapVV.size(), 0);
        EXPECT_EQ(finder01.mMapVE.size(), 0);
        EXPECT_EQ(finder01.mMapVF.size(), 0);
        EXPECT_EQ(finder01.mMapEV.size(), 0);
        EXPECT_EQ(finder01.mMapEE.size(), 0);
        EXPECT_EQ(finder01.mMapEF.size(), 0);
        EXPECT_EQ(finder01.mMapFV.size(), 0);
        EXPECT_EQ(finder01.mMapFE.size(), 0);
        EXPECT_EQ(finder01.mMapVI.size(), 0);
        EXPECT_EQ(finder01.mMapIV.size(), 0);
    }
}


/**  @brief findIntersectionsEdgesOn1FacesOn2()
 *          findIntersectionsFacesOn1EdgesOn2()
 *          findRemainingInteriorVertices()
 */
TEST_F(IntersectionFinderTest, Test08) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(0.0, 0.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();

        EXPECT_EQ(finder01.mMapVV.size(), 8);
        EXPECT_EQ(finder01.mMapVV.size(), 8);
        EXPECT_EQ(finder01.mMapVE.size(), 0);
        EXPECT_EQ(finder01.mMapVF.size(), 0);
        EXPECT_EQ(finder01.mMapEV.size(), 0);
        EXPECT_EQ(finder01.mMapEE.size(), 0);
        EXPECT_EQ(finder01.mMapEF.size(), 0);
        EXPECT_EQ(finder01.mMapFV.size(), 0);
        EXPECT_EQ(finder01.mMapFE.size(), 0);
        EXPECT_EQ(finder01.mMapVI.size(), 0);
        EXPECT_EQ(finder01.mMapIV.size(), 0);

        map<long, long> Vmap01;
        map<long, long> Vmap02;
        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            Vmap01[(*vit)->id()] = 0;
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            Vmap02[(*vit)->id()] = 0;
        }

        auto Qmat01Inv = Qmat01.transpose();

        for (auto&e : finder01.mMapVV) {
            Vmap01[e.first.first ] = Vmap01[e.first.first ] + 1;
            Vmap02[e.first.second] = Vmap02[e.first.second] + 1;
            auto vit1 = m01.vertexIt(e.first.first);
            auto vit2 = m02.vertexIt(e.first.second);
            EXPECT_EQ((*vit1)->pLCS(), Qmat01Inv * (e.second - CoM01));
            EXPECT_EQ((*vit2)->pLCS(), Qmat01Inv * (e.second - CoM01));
        }

        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            EXPECT_EQ(Vmap01[(*vit)->id()], 1);
        }
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            EXPECT_EQ(Vmap02[(*vit)->id()], 1);
        }
    }
}


/**  @brief findIntersectionsEdgesOn1FacesOn2()
 *          findIntersectionsFacesOn1EdgesOn2()
 *          findRemainingInteriorVertices()
 */
TEST_F(IntersectionFinderTest, Test09) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 0.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();

        EXPECT_EQ(finder01.mMapVV.size(), 4);
        EXPECT_EQ(finder01.mMapVE.size(), 0);
        EXPECT_EQ(finder01.mMapVF.size(), 0);
        EXPECT_EQ(finder01.mMapEV.size(), 0);
        EXPECT_EQ(finder01.mMapEE.size(), 0);
        EXPECT_EQ(finder01.mMapEF.size(), 0);
        EXPECT_EQ(finder01.mMapFV.size(), 0);
        EXPECT_EQ(finder01.mMapFE.size(), 0);
        EXPECT_EQ(finder01.mMapVI.size(), 0);
        EXPECT_EQ(finder01.mMapIV.size(), 0);

        map<long, long> Vmap01;
        map<long, long> Vmap02;
        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            Vmap01[(*vit)->id()] = 0;
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            Vmap02[(*vit)->id()] = 0;
        }

        auto Qmat01Inv = Qmat01.transpose();
        auto Qmat02Inv = Qmat02.transpose();

        for (auto&e : finder01.mMapVV) {
            Vmap01[e.first.first ] = Vmap01[e.first.first ] + 1;
            Vmap02[e.first.second] = Vmap02[e.first.second] + 1;
            auto vit1 = m01.vertexIt(e.first.first);
            auto vit2 = m02.vertexIt(e.first.second);
            EXPECT_EQ((*vit1)->pLCS(), Qmat01Inv * (e.second - CoM01));
            EXPECT_EQ((*vit2)->pLCS(), Qmat02Inv * (e.second - CoM02));
        }

        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()==1.0) {
                EXPECT_EQ(Vmap01[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap01[(*vit)->id()], 0);
            }
        }
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0) {
                EXPECT_EQ(Vmap02[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap02[(*vit)->id()], 0);
            }
        }
    }
}


/**  @brief findIntersectionsEdgesOn1FacesOn2()
 *          findIntersectionsFacesOn1EdgesOn2()
 *          findRemainingInteriorVertices()
 */
TEST_F(IntersectionFinderTest, Test10) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 1.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();

        EXPECT_EQ(finder01.mMapVV.size(), 0);
        EXPECT_EQ(finder01.mMapVE.size(), 2);
        EXPECT_EQ(finder01.mMapVF.size(), 0);
        EXPECT_EQ(finder01.mMapEV.size(), 2);
        EXPECT_EQ(finder01.mMapEE.size(), 0);
        EXPECT_EQ(finder01.mMapEF.size(), 0);
        EXPECT_EQ(finder01.mMapFV.size(), 0);
        EXPECT_EQ(finder01.mMapFE.size(), 0);
        EXPECT_EQ(finder01.mMapVI.size(), 0);
        EXPECT_EQ(finder01.mMapIV.size(), 0);

        map<long, long> Vmap01;
        map<long, long> Vmap02;
        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            Vmap01[(*vit)->id()] = 0;
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            Vmap02[(*vit)->id()] = 0;
        }

        map<pair<long,long>,long> Emap01;
        map<pair<long,long>,long> Emap02;
        auto ePair01 = m01.edges();
        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {
            auto id1 = (*((*((*eit)->he1()))->src()))->id();
            auto id2 = (*((*((*eit)->he1()))->dst()))->id();
            if (id1 > id2) { std::swap(id1,id2); }
            Emap01[make_pair(id1, id2)] = 0;
        }

        auto ePair02 = m02.edges();
        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {
            auto id1 = (*((*((*eit)->he1()))->src()))->id();
            auto id2 = (*((*((*eit)->he1()))->dst()))->id();
            if (id1 > id2) { std::swap(id1,id2); }
            Emap02[make_pair(id1, id2)] = 0;
        }

        auto Qmat01Inv = Qmat01.transpose();
        auto Qmat02Inv = Qmat02.transpose();

        for (auto&e : finder01.mMapEV) {
            auto eid = e.first.first;
            auto vid = e.first.second;
            Emap01[eid] = Emap01[eid] + 1;
            Vmap02[vid] = Vmap02[vid] + 1;
            auto eit1 = m01.edgeIt(eid);
            auto vit2 = m02.vertexIt(vid);
            EXPECT_EQ((*vit2)->pLCS(), Qmat02Inv * (e.second - CoM02));
        }

        for (auto&e : finder01.mMapVE) {
            auto vid = e.first.first;
            auto eid = e.first.second;
            Vmap01[vid] = Vmap01[vid] + 1;
            Emap02[eid] = Emap02[eid] + 1;
            auto vit1 = m01.vertexIt(vid);
            auto eit2 = m02.edgeIt(eid);
            EXPECT_EQ((*vit1)->pLCS(), Qmat01Inv * (e.second - CoM01));
        }

        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()==1.0&&(*vit)->pLCS().y()==1.0) {
                EXPECT_EQ(Vmap01[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap01[(*vit)->id()], 0);
            }
        }

        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0&&(*vit)->pLCS().y()==-1.0) {
                EXPECT_EQ(Vmap02[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap02[(*vit)->id()], 0);
            }
        }

        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {

            auto src = (*((*eit)->he1()))->src();
            auto dst = (*((*eit)->he1()))->dst();

            auto id1 = (*src)->id();
            auto id2 = (*dst)->id();
            if (id1 > id2) { std::swap(id1,id2); }

            if ( ( ((*src)->pLCS().x()==1.0 && (*src)->pLCS().z()== 1.0)&&
                   ((*dst)->pLCS().x()==1.0 && (*dst)->pLCS().z()== 1.0)  )||
                 ( ((*src)->pLCS().x()==1.0 && (*src)->pLCS().z()== -1.0)&&
                   ((*dst)->pLCS().x()==1.0 && (*dst)->pLCS().z()== -1.0)  )) {
                EXPECT_EQ(Emap01[make_pair(id1, id2)], 1);
            }
            else {
                EXPECT_EQ(Emap01[make_pair(id1, id2)], 0);
            }
        }

        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {

            auto src = (*((*eit)->he1()))->src();
            auto dst = (*((*eit)->he1()))->dst();

            auto id1 = (*src)->id();
            auto id2 = (*dst)->id();
            if (id1 > id2) { std::swap(id1,id2); }

            if ( ( ((*src)->pLCS().x()==-1.0 && (*src)->pLCS().z()== 1.0)&&
                   ((*dst)->pLCS().x()==-1.0 && (*dst)->pLCS().z()== 1.0)  )||
                 ( ((*src)->pLCS().x()==-1.0 && (*src)->pLCS().z()== -1.0)&&
                   ((*dst)->pLCS().x()==-1.0 && (*dst)->pLCS().z()== -1.0) )){
                EXPECT_EQ(Emap02[make_pair(id1, id2)], 1);
            }
            else {
                EXPECT_EQ(Emap02[make_pair(id1, id2)], 0);
            }
        }
    }
}


/**  @brief findIntersectionsEdgesOn1FacesOn2()
 *          findIntersectionsFacesOn1EdgesOn2()
 *          findRemainingInteriorVertices()
 */
TEST_F(IntersectionFinderTest, Test11) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 1.0, 1.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();

        EXPECT_EQ(finder01.mMapVV.size(), 0);
        EXPECT_EQ(finder01.mMapVE.size(), 0);
        EXPECT_EQ(finder01.mMapVF.size(), 1);
        EXPECT_EQ(finder01.mMapEV.size(), 0);
        EXPECT_EQ(finder01.mMapEE.size(), 2);
        EXPECT_EQ(finder01.mMapEF.size(), 0);
        EXPECT_EQ(finder01.mMapFV.size(), 1);
        EXPECT_EQ(finder01.mMapFE.size(), 0);
        EXPECT_EQ(finder01.mMapVI.size(), 0);
        EXPECT_EQ(finder01.mMapIV.size(), 0);

        map<long, long> Vmap01;
        map<long, long> Vmap02;
        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            Vmap01[(*vit)->id()] = 0;
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            Vmap02[(*vit)->id()] = 0;
        }

        map<pair<long,long>,long> Emap01;
        map<pair<long,long>,long> Emap02;
        auto ePair01 = m01.edges();
        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {
            auto id1 = (*((*((*eit)->he1()))->src()))->id();
            auto id2 = (*((*((*eit)->he1()))->dst()))->id();
            if (id1 > id2) { std::swap(id1,id2); }
            Emap01[make_pair(id1, id2)] = 0;
        }

        auto ePair02 = m02.edges();
        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {
            auto id1 = (*((*((*eit)->he1()))->src()))->id();
            auto id2 = (*((*((*eit)->he1()))->dst()))->id();
            if (id1 > id2) { std::swap(id1,id2); }
            Emap02[make_pair(id1, id2)] = 0;
        }

        map<long,long> Fmap01;
        map<long,long> Fmap02;
        auto fPair01 = m01.faces();
        for (auto fit = fPair01.first; fit != fPair01.second; fit++) {
            Fmap01[(*fit)->id()] = 0;
        }

        auto fPair02 = m02.faces();
        for (auto fit = fPair02.first; fit != fPair02.second; fit++) {
            Fmap02[(*fit)->id()] = 0;
        }

        auto Qmat01Inv = Qmat01.transpose();
        auto Qmat02Inv = Qmat02.transpose();

        for (auto&e : finder01.mMapEE) {
            auto eid1 = e.first.first;
            auto eid2 = e.first.second;
            Emap01[eid1] = Emap01[eid1] + 1;
            Emap02[eid2] = Emap02[eid2] + 1;
            auto eit1 = m01.edgeIt(eid1);
            auto eit2 = m02.edgeIt(eid2);
            EXPECT_EQ(
                ( Vec3( 1.0,  0.0,  1.0) == Qmat01Inv * (e.second - CoM01) &&
                  Vec3(-1.0, -1.0,  0.0) == Qmat02Inv * (e.second - CoM02) )||
                ( Vec3( 1.0,  1.0,  0.0) == Qmat01Inv * (e.second - CoM01) &&
                  Vec3(-1.0,  0.0, -1.0) == Qmat02Inv * (e.second - CoM02) ), 
                true);
        }

        for (auto&e : finder01.mMapVF) {
            auto vid = e.first.first;
            auto fid = e.first.second;
            Vmap01[vid] = Vmap01[vid] + 1;
            Fmap02[fid] = Fmap02[fid] + 1;
            auto vit1 = m01.vertexIt(vid);
            auto fit2 = m02.faceIt(fid);
            EXPECT_EQ((*vit1)->pLCS(), Qmat01Inv * (e.second - CoM01));
        }

        for (auto&e : finder01.mMapFV) {
            auto fid = e.first.first;
            auto vid = e.first.second;
            Fmap01[fid] = Vmap01[fid] + 1;
            Vmap02[vid] = Fmap02[vid] + 1;
            auto fit1 = m01.faceIt(fid);
            auto vit2 = m02.vertexIt(vid);
            EXPECT_EQ((*vit2)->pLCS(), Qmat02Inv * (e.second - CoM02));
        }

        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()==1.0&&(*vit)->pLCS().y()==1.0&&
                (*vit)->pLCS().z()==1.0 ) {
                EXPECT_EQ(Vmap01[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap01[(*vit)->id()], 0);
            }
        }

        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0&&(*vit)->pLCS().y()==-1.0&&
                (*vit)->pLCS().z()==-1.0 ) {
                EXPECT_EQ(Vmap02[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap02[(*vit)->id()], 0);
            }
        }

        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {

            auto src = (*((*eit)->he1()))->src();
            auto dst = (*((*eit)->he1()))->dst();

            auto id1 = (*src)->id();
            auto id2 = (*dst)->id();
            if (id1 > id2) { std::swap(id1,id2); }

            if ( ( ((*src)->pLCS().x()==1.0 && (*src)->pLCS().z()== 1.0)&&
                   ((*dst)->pLCS().x()==1.0 && (*dst)->pLCS().z()== 1.0)  )||
                 ( ((*src)->pLCS().x()==1.0 && (*src)->pLCS().y()== 1.0)&&
                   ((*dst)->pLCS().x()==1.0 && (*dst)->pLCS().y()== 1.0)  )){
                EXPECT_EQ(Emap01[make_pair(id1, id2)], 1);
            }
            else {
                EXPECT_EQ(Emap01[make_pair(id1, id2)], 0);
            }
        }

        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {

            auto src = (*((*eit)->he1()))->src();
            auto dst = (*((*eit)->he1()))->dst();

            auto id1 = (*src)->id();
            auto id2 = (*dst)->id();
            if (id1 > id2) { std::swap(id1,id2); }

            if ( ( ((*src)->pLCS().x()==-1.0 && (*src)->pLCS().z()== -1.0)&&
                   ((*dst)->pLCS().x()==-1.0 && (*dst)->pLCS().z()== -1.0)  )||
                 ( ((*src)->pLCS().x()==-1.0 && (*src)->pLCS().y()== -1.0)&&
                   ((*dst)->pLCS().x()==-1.0 && (*dst)->pLCS().y()== -1.0) )){
                EXPECT_EQ(Emap02[make_pair(id1, id2)], 1);
            }
            else {
                EXPECT_EQ(Emap02[make_pair(id1, id2)], 0);
            }
        }

        for (auto fit = fPair01.first; fit != fPair01.second; fit++) {
            if ((*fit)->nLCS() == Vec3(1.0, 0.0, 0.0)) {
                EXPECT_EQ(Fmap01[(*fit)->id()], 1);
            }
            else {
                EXPECT_EQ(Fmap01[(*fit)->id()], 0);
            }
        }

        for (auto fit = fPair02.first; fit != fPair02.second; fit++) {
            if ((*fit)->nLCS() == Vec3(-1.0, 0.0, 0.0)) {
                EXPECT_EQ(Fmap02[(*fit)->id()], 1);
            }
            else {
                EXPECT_EQ(Fmap02[(*fit)->id()], 0);
            }
        }
    }
}


/**  @brief findIntersectionsEdgesOn1FacesOn2()
 *          findIntersectionsFacesOn1EdgesOn2()
 *          findRemainingInteriorVertices()
 */
TEST_F(IntersectionFinderTest, Test12) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(1.0, 1.0, 1.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();

        EXPECT_EQ(finder01.mMapVV.size(), 0);
        EXPECT_EQ(finder01.mMapVE.size(), 0);
        EXPECT_EQ(finder01.mMapVF.size(), 0);
        EXPECT_EQ(finder01.mMapEV.size(), 0);
        EXPECT_EQ(finder01.mMapEE.size(), 0);
        EXPECT_EQ(finder01.mMapEF.size(), 3);
        EXPECT_EQ(finder01.mMapFV.size(), 0);
        EXPECT_EQ(finder01.mMapFE.size(), 3);
        EXPECT_EQ(finder01.mMapVI.size(), 1);
        EXPECT_EQ(finder01.mMapIV.size(), 1);

        map<long, long> Vmap01;
        map<long, long> Vmap02;
        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            Vmap01[(*vit)->id()] = 0;
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            Vmap02[(*vit)->id()] = 0;
        }

        map<pair<long,long>,long> Emap01;
        map<pair<long,long>,long> Emap02;
        auto ePair01 = m01.edges();
        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {
            auto id1 = (*((*((*eit)->he1()))->src()))->id();
            auto id2 = (*((*((*eit)->he1()))->dst()))->id();
            if (id1 > id2) { std::swap(id1,id2); }
            Emap01[make_pair(id1, id2)] = 0;
        }

        auto ePair02 = m02.edges();
        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {
            auto id1 = (*((*((*eit)->he1()))->src()))->id();
            auto id2 = (*((*((*eit)->he1()))->dst()))->id();
            if (id1 > id2) { std::swap(id1,id2); }
            Emap02[make_pair(id1, id2)] = 0;
        }

        map<long,long> Fmap01;
        map<long,long> Fmap02;
        auto fPair01 = m01.faces();
        for (auto fit = fPair01.first; fit != fPair01.second; fit++) {
            Fmap01[(*fit)->id()] = 0;
        }

        auto fPair02 = m02.faces();
        for (auto fit = fPair02.first; fit != fPair02.second; fit++) {
            Fmap02[(*fit)->id()] = 0;
        }

        auto Qmat01Inv = Qmat01.transpose();
        auto Qmat02Inv = Qmat02.transpose();

        for (auto&e : finder01.mMapEF) {
            auto eid1 = e.first.first;
            auto fid2 = e.first.second;
            Emap01[eid1] = Emap01[eid1] + 1;
            Fmap02[fid2] = Fmap02[fid2] + 1;
            auto eit1 = m01.edgeIt(eid1);
            auto fit2 = m02.faceIt(fid2);
            auto he = (*eit1)->he1();
            auto p1 = (*((*he)->src()))->pLCS();
            auto p2 = (*((*he)->dst()))->pLCS();
            
            EXPECT_EQ(p1*0.5 + p2*0.5, Qmat01Inv * (e.second - CoM01));

        }

        for (auto&e : finder01.mMapFE) {
            auto fid1 = e.first.first;
            auto eid2 = e.first.second;
            Fmap01[fid1] = Fmap01[fid1] + 1;
            Emap02[eid2] = Emap02[eid2] + 1;
            auto fit1 = m01.faceIt(fid1);
            auto eit2 = m02.edgeIt(eid2);
            auto he = (*eit2)->he1();
            auto p1 = (*((*he)->src()))->pLCS();
            auto p2 = (*((*he)->dst()))->pLCS();
            
            EXPECT_EQ(p1*0.5 + p2*0.5, Qmat02Inv * (e.second - CoM02));

        }

        for (auto&e : finder01.mMapVI) {
            auto vid1 = e.first;
            Vmap01[vid1] = Vmap01[vid1] + 1;
            auto vit1 = m01.vertexIt(vid1);
            auto p1 = (*vit1)->pLCS();
            EXPECT_EQ(p1, Qmat01Inv * (e.second - CoM01));
        }

        for (auto&e : finder01.mMapIV) {
            auto vid2 = e.first;
            Vmap02[vid2] = Vmap02[vid2] + 1;
            auto vit2 = m01.vertexIt(vid2);
            auto p2 = (*vit2)->pLCS();
            EXPECT_EQ(p2, Qmat02Inv * (e.second - CoM02));
        }

        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()==1.0&&(*vit)->pLCS().y()==1.0&&
                (*vit)->pLCS().z()==1.0 ) {
                EXPECT_EQ(Vmap01[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap01[(*vit)->id()], 0);
            }
        }

        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0&&(*vit)->pLCS().y()==-1.0&&
                (*vit)->pLCS().z()==-1.0 ) {
                EXPECT_EQ(Vmap02[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap02[(*vit)->id()], 0);
            }
        }

        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {

            auto src = (*((*eit)->he1()))->src();
            auto dst = (*((*eit)->he1()))->dst();

            auto id1 = (*src)->id();
            auto id2 = (*dst)->id();
            if (id1 > id2) { std::swap(id1,id2); }

            if ( ( ((*src)->pLCS().x()==1.0 && (*src)->pLCS().z()== 1.0)&&
                   ((*dst)->pLCS().x()==1.0 && (*dst)->pLCS().z()== 1.0)  )||
                 ( ((*src)->pLCS().x()==1.0 && (*src)->pLCS().y()== 1.0)&&
                   ((*dst)->pLCS().x()==1.0 && (*dst)->pLCS().y()== 1.0)  )||
                 ( ((*src)->pLCS().y()==1.0 && (*src)->pLCS().z()== 1.0)&&
                ((*dst)->pLCS().y()==1.0 && (*dst)->pLCS().z()== 1.0)  )  ) {

                EXPECT_EQ(Emap01[make_pair(id1, id2)], 1);
            }
            else {
                EXPECT_EQ(Emap01[make_pair(id1, id2)], 0);
            }
        }

        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {

            auto src = (*((*eit)->he1()))->src();
            auto dst = (*((*eit)->he1()))->dst();

            auto id1 = (*src)->id();
            auto id2 = (*dst)->id();
            if (id1 > id2) { std::swap(id1,id2); }

            if ( ( ((*src)->pLCS().x()== -1.0 && (*src)->pLCS().z()== -1.0)&&
                  ((*dst)->pLCS().x()== -1.0 && (*dst)->pLCS().z()== -1.0)  )||
                 ( ((*src)->pLCS().x()== -1.0 && (*src)->pLCS().y()== -1.0)&&
                  ((*dst)->pLCS().x()== -1.0 && (*dst)->pLCS().y()== -1.0)  )||
                 ( ((*src)->pLCS().y()== -1.0 && (*src)->pLCS().z()== -1.0)&&
               ((*dst)->pLCS().y()== -1.0 && (*dst)->pLCS().z()== -1.0)  )  ) {

                EXPECT_EQ(Emap02[make_pair(id1, id2)], 1);
            }
            else {
                EXPECT_EQ(Emap02[make_pair(id1, id2)], 0);
            }
        }

        for (auto fit = fPair01.first; fit != fPair01.second; fit++) {
            if ((*fit)->nLCS() == Vec3(1.0, 0.0, 0.0)||
                (*fit)->nLCS() == Vec3(0.0, 1.0, 0.0)||
                (*fit)->nLCS() == Vec3(0.0, 0.0, 1.0)  ) {
                EXPECT_EQ(Fmap01[(*fit)->id()], 1);
            }
            else {
                EXPECT_EQ(Fmap01[(*fit)->id()], 0);
            }
        }

        for (auto fit = fPair02.first; fit != fPair02.second; fit++) {
            if ((*fit)->nLCS() == Vec3(-1.0,  0.0,  0.0)||
                (*fit)->nLCS() == Vec3( 0.0, -1.0,  0.0)||
                (*fit)->nLCS() == Vec3( 0.0,  0.0, -1.0)  ) {
                EXPECT_EQ(Fmap02[(*fit)->id()], 1);
            }
            else {
                EXPECT_EQ(Fmap02[(*fit)->id()], 0);
            }
        }
    }
}


/**  @brief findIntersectionsEdgesOn1FacesOn2()
 *          findIntersectionsFacesOn1EdgesOn2()
 *          findRemainingInteriorVertices()
 */
TEST_F(IntersectionFinderTest, Test13) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(1.0, 0.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();

        EXPECT_EQ(finder01.mMapVV.size(), 0);
        EXPECT_EQ(finder01.mMapVE.size(), 4);
        EXPECT_EQ(finder01.mMapVF.size(), 0);
        EXPECT_EQ(finder01.mMapEV.size(), 4);
        EXPECT_EQ(finder01.mMapEE.size(), 0);
        EXPECT_EQ(finder01.mMapEF.size(), 0);
        EXPECT_EQ(finder01.mMapFV.size(), 0);
        EXPECT_EQ(finder01.mMapFE.size(), 0);
        EXPECT_EQ(finder01.mMapVI.size(), 0);
        EXPECT_EQ(finder01.mMapIV.size(), 0);


        map<long, long> Vmap01;
        map<long, long> Vmap02;
        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            Vmap01[(*vit)->id()] = 0;
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            Vmap02[(*vit)->id()] = 0;
        }

        map<pair<long,long>,long> Emap01;
        map<pair<long,long>,long> Emap02;
        auto ePair01 = m01.edges();
        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {
            auto id1 = (*((*((*eit)->he1()))->src()))->id();
            auto id2 = (*((*((*eit)->he1()))->dst()))->id();
            if (id1 > id2) { std::swap(id1,id2); }
            Emap01[make_pair(id1, id2)] = 0;
        }

        auto ePair02 = m02.edges();
        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {
            auto id1 = (*((*((*eit)->he1()))->src()))->id();
            auto id2 = (*((*((*eit)->he1()))->dst()))->id();
            if (id1 > id2) { std::swap(id1,id2); }
            Emap02[make_pair(id1, id2)] = 0;
        }

        auto Qmat01Inv = Qmat01.transpose();
        auto Qmat02Inv = Qmat02.transpose();

        for (auto&e : finder01.mMapEV) {
            auto eid = e.first.first;
            auto vid = e.first.second;
            Emap01[eid] = Emap01[eid] + 1;
            Vmap02[vid] = Vmap02[vid] + 1;
            auto eit1 = m01.edgeIt(eid);
            auto vit2 = m02.vertexIt(vid);
            EXPECT_EQ((*vit2)->pLCS(), Qmat02Inv * (e.second-CoM02));
        }

        for (auto&e : finder01.mMapVE) {
            auto vid = e.first.first;
            auto eid = e.first.second;
            Vmap01[vid] = Vmap01[vid] + 1;
            Emap02[eid] = Emap02[eid] + 1;
            auto vit1 = m01.vertexIt(vid);
            auto eit2 = m02.edgeIt(eid);
            EXPECT_EQ((*vit1)->pLCS(), Qmat01Inv * (e.second - CoM01));
        }

        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()==1.0) {
                EXPECT_EQ(Vmap01[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap01[(*vit)->id()], 0);
            }
        }

        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0) {
                EXPECT_EQ(Vmap02[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap02[(*vit)->id()], 0);
            }
        }

        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {

            auto src = (*((*eit)->he1()))->src();
            auto dst = (*((*eit)->he1()))->dst();

            auto id1 = (*src)->id();
            auto id2 = (*dst)->id();
            if (id1 > id2) { std::swap(id1,id2); }

            if ( ( ((*src)->pLCS().y()==1.0 && (*src)->pLCS().z()==  1.0)&&
                   ((*dst)->pLCS().y()==1.0 && (*dst)->pLCS().z()==  1.0)  )||
                 ( ((*src)->pLCS().y()==1.0 && (*src)->pLCS().z()== -1.0)&&
                   ((*dst)->pLCS().y()==1.0 && (*dst)->pLCS().z()== -1.0) )||
                 ( ((*src)->pLCS().y()==-1.0 && (*src)->pLCS().z()==  1.0)&&
                   ((*dst)->pLCS().y()==-1.0 && (*dst)->pLCS().z()==  1.0)  )||
                 ( ((*src)->pLCS().y()==-1.0 && (*src)->pLCS().z()== -1.0)&&
                   ((*dst)->pLCS().y()==-1.0 && (*dst)->pLCS().z()== -1.0)  )){
                EXPECT_EQ(Emap01[make_pair(id1, id2)], 1);
            }
            else {
                EXPECT_EQ(Emap01[make_pair(id1, id2)], 0);
            }
        }

        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {

            auto src = (*((*eit)->he1()))->src();
            auto dst = (*((*eit)->he1()))->dst();

            auto id1 = (*src)->id();
            auto id2 = (*dst)->id();
            if (id1 > id2) { std::swap(id1,id2); }

            if ( ( ((*src)->pLCS().y()==1.0 && (*src)->pLCS().z()==  1.0)&&
                   ((*dst)->pLCS().y()==1.0 && (*dst)->pLCS().z()==  1.0)  )||
                 ( ((*src)->pLCS().y()==1.0 && (*src)->pLCS().z()== -1.0)&&
                   ((*dst)->pLCS().y()==1.0 && (*dst)->pLCS().z()== -1.0) )||
                 ( ((*src)->pLCS().y()==-1.0 && (*src)->pLCS().z()==  1.0)&&
                   ((*dst)->pLCS().y()==-1.0 && (*dst)->pLCS().z()==  1.0)  )||
                 ( ((*src)->pLCS().y()==-1.0 && (*src)->pLCS().z()== -1.0)&&
                   ((*dst)->pLCS().y()==-1.0 && (*dst)->pLCS().z()== -1.0)  )){
                EXPECT_EQ(Emap02[make_pair(id1, id2)], 1);
            }
            else {
                EXPECT_EQ(Emap02[make_pair(id1, id2)], 0);
            }
        }
    }
}


/**  @brief findIntersectionsEdgesOn1FacesOn2()
 *          findIntersectionsFacesOn1EdgesOn2()
 *          findRemainingInteriorVertices()
 */
TEST_F(IntersectionFinderTest, Test14) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 2.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();

        EXPECT_EQ(finder01.mMapVV.size(), 2);
        EXPECT_EQ(finder01.mMapVE.size(), 0);
        EXPECT_EQ(finder01.mMapVF.size(), 0);
        EXPECT_EQ(finder01.mMapEV.size(), 0);
        EXPECT_EQ(finder01.mMapEE.size(), 0);
        EXPECT_EQ(finder01.mMapEF.size(), 0);
        EXPECT_EQ(finder01.mMapFV.size(), 0);
        EXPECT_EQ(finder01.mMapFE.size(), 0);
        EXPECT_EQ(finder01.mMapVI.size(), 0);
        EXPECT_EQ(finder01.mMapIV.size(), 0);

        map<long, long> Vmap01;
        map<long, long> Vmap02;
        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            Vmap01[(*vit)->id()] = 0;
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            Vmap02[(*vit)->id()] = 0;
        }

        auto Qmat01Inv = Qmat01.transpose();
        auto Qmat02Inv = Qmat02.transpose();

        for (auto&e : finder01.mMapVV) {
            Vmap01[e.first.first ] = Vmap01[e.first.first ] + 1;
            Vmap02[e.first.second] = Vmap02[e.first.second] + 1;
            auto vit1 = m01.vertexIt(e.first.first);
            auto vit2 = m02.vertexIt(e.first.second);
            EXPECT_EQ((*vit1)->pLCS(), Qmat01Inv * (e.second - CoM01));
            EXPECT_EQ((*vit2)->pLCS(), Qmat02Inv * (e.second - CoM02));
        }

        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()==1.0 && (*vit)->pLCS().y()==1.0 ) {
                EXPECT_EQ(Vmap01[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap01[(*vit)->id()], 0);
            }
        }
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0 && (*vit)->pLCS().y()==-1.0 ) {
                EXPECT_EQ(Vmap02[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap02[(*vit)->id()], 0);
            }
        }
    }
}


/**  @brief findIntersectionsEdgesOn1FacesOn2()
 *          findIntersectionsFacesOn1EdgesOn2()
 *          findRemainingInteriorVertices()
 */
TEST_F(IntersectionFinderTest, Test15) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 2.0, 1.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();

        EXPECT_EQ(finder01.mMapVV.size(), 0);
        EXPECT_EQ(finder01.mMapVE.size(), 1);
        EXPECT_EQ(finder01.mMapVF.size(), 0);
        EXPECT_EQ(finder01.mMapEV.size(), 1);
        EXPECT_EQ(finder01.mMapEE.size(), 0);
        EXPECT_EQ(finder01.mMapEF.size(), 0);
        EXPECT_EQ(finder01.mMapFV.size(), 0);
        EXPECT_EQ(finder01.mMapFE.size(), 0);
        EXPECT_EQ(finder01.mMapVI.size(), 0);
        EXPECT_EQ(finder01.mMapIV.size(), 0);


        map<long, long> Vmap01;
        map<long, long> Vmap02;
        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            Vmap01[(*vit)->id()] = 0;
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            Vmap02[(*vit)->id()] = 0;
        }

        map<pair<long,long>,long> Emap01;
        map<pair<long,long>,long> Emap02;
        auto ePair01 = m01.edges();
        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {
            auto id1 = (*((*((*eit)->he1()))->src()))->id();
            auto id2 = (*((*((*eit)->he1()))->dst()))->id();
            if (id1 > id2) { std::swap(id1,id2); }
            Emap01[make_pair(id1, id2)] = 0;
        }

        auto ePair02 = m02.edges();
        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {
            auto id1 = (*((*((*eit)->he1()))->src()))->id();
            auto id2 = (*((*((*eit)->he1()))->dst()))->id();
            if (id1 > id2) { std::swap(id1,id2); }
            Emap02[make_pair(id1, id2)] = 0;
        }

        auto Qmat01Inv = Qmat01.transpose();
        auto Qmat02Inv = Qmat02.transpose();

        for (auto&e : finder01.mMapEV) {
            auto eid = e.first.first;
            auto vid = e.first.second;
            Emap01[eid] = Emap01[eid] + 1;
            Vmap02[vid] = Vmap02[vid] + 1;
            auto eit1 = m01.edgeIt(eid);
            auto vit2 = m02.vertexIt(vid);
            EXPECT_EQ((*vit2)->pLCS(), Qmat02Inv * (e.second-CoM02));
        }

        for (auto&e : finder01.mMapVE) {
            auto vid = e.first.first;
            auto eid = e.first.second;
            Vmap01[vid] = Vmap01[vid] + 1;
            Emap02[eid] = Emap02[eid] + 1;
            auto vit1 = m01.vertexIt(vid);
            auto eit2 = m02.edgeIt(eid);
            EXPECT_EQ((*vit1)->pLCS(), Qmat01Inv * (e.second - CoM01));
        }

        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()==1.0 && (*vit)->pLCS().y()==1.0 &&
                (*vit)->pLCS().z()==1.0                               ) {
                EXPECT_EQ(Vmap01[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap01[(*vit)->id()], 0);
            }
        }

        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0 && (*vit)->pLCS().y()==-1.0 &&
                (*vit)->pLCS().z()==-1.0                               ) {
                EXPECT_EQ(Vmap02[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap02[(*vit)->id()], 0);
            }
        }

        for (auto eit = ePair01.first; eit != ePair01.second; eit++) {

            auto src = (*((*eit)->he1()))->src();
            auto dst = (*((*eit)->he1()))->dst();

            auto id1 = (*src)->id();
            auto id2 = (*dst)->id();
            if (id1 > id2) { std::swap(id1,id2); }

            if ( ( ((*src)->pLCS().x()==1.0 && (*src)->pLCS().y()==  1.0)&&
                   ((*dst)->pLCS().x()==1.0 && (*dst)->pLCS().y()==  1.0)  )) {
                EXPECT_EQ(Emap01[make_pair(id1, id2)], 1);
            }
            else {
                EXPECT_EQ(Emap01[make_pair(id1, id2)], 0);
            }
        }

        for (auto eit = ePair02.first; eit != ePair02.second; eit++) {

            auto src = (*((*eit)->he1()))->src();
            auto dst = (*((*eit)->he1()))->dst();

            auto id1 = (*src)->id();
            auto id2 = (*dst)->id();
            if (id1 > id2) { std::swap(id1,id2); }

            if ( ( ((*src)->pLCS().x()==-1.0 && (*src)->pLCS().y()== -1.0)&&
                   ((*dst)->pLCS().x()==-1.0 && (*dst)->pLCS().y()== -1.0) )){
                EXPECT_EQ(Emap02[make_pair(id1, id2)], 1);
            }
            else {
                EXPECT_EQ(Emap02[make_pair(id1, id2)], 0);
            }
        }
    }
}


/**  @brief findIntersectionsEdgesOn1FacesOn2()
 *          findIntersectionsFacesOn1EdgesOn2()
 *          findRemainingInteriorVertices()
 */
TEST_F(IntersectionFinderTest, Test16) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 2.0, 2.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();

        EXPECT_EQ(finder01.mMapVV.size(), 1);
        EXPECT_EQ(finder01.mMapVE.size(), 0);
        EXPECT_EQ(finder01.mMapVF.size(), 0);
        EXPECT_EQ(finder01.mMapEV.size(), 0);
        EXPECT_EQ(finder01.mMapEE.size(), 0);
        EXPECT_EQ(finder01.mMapEF.size(), 0);
        EXPECT_EQ(finder01.mMapFV.size(), 0);
        EXPECT_EQ(finder01.mMapFE.size(), 0);
        EXPECT_EQ(finder01.mMapVI.size(), 0);
        EXPECT_EQ(finder01.mMapIV.size(), 0);

        map<long, long> Vmap01;
        map<long, long> Vmap02;
        auto vPair01 = m01.vertices();
        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            Vmap01[(*vit)->id()] = 0;
        }

        auto vPair02 = m02.vertices();
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            Vmap02[(*vit)->id()] = 0;
        }

        auto Qmat01Inv = Qmat01.transpose();
        auto Qmat02Inv = Qmat02.transpose();

        for (auto&e : finder01.mMapVV) {
            Vmap01[e.first.first ] = Vmap01[e.first.first ] + 1;
            Vmap02[e.first.second] = Vmap02[e.first.second] + 1;
            auto vit1 = m01.vertexIt(e.first.first);
            auto vit2 = m02.vertexIt(e.first.second);
            EXPECT_EQ((*vit1)->pLCS(), Qmat01Inv * (e.second - CoM01));
            EXPECT_EQ((*vit2)->pLCS(), Qmat02Inv * (e.second - CoM02));
        }

        for (auto vit = vPair01.first; vit != vPair01.second; vit++) {
            if ((*vit)->pLCS().x()==1.0 && (*vit)->pLCS().y()==1.0 &&
                (*vit)->pLCS().z()==1.0                               ) {
                EXPECT_EQ(Vmap01[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap01[(*vit)->id()], 0);
            }
        }
        for (auto vit = vPair02.first; vit != vPair02.second; vit++) {
            if ((*vit)->pLCS().x()==-1.0 && (*vit)->pLCS().y()==-1.0 &&
                (*vit)->pLCS().z()==-1.0                               ) {
                EXPECT_EQ(Vmap02[(*vit)->id()], 1);
            }
            else {
                EXPECT_EQ(Vmap02[(*vit)->id()], 0);
            }
        }
    }
}

/**  @brief constructVertexAttributes() and findDimension();
 *          from Test07
 */
TEST_F(IntersectionFinderTest, Test17) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(4.0, rand100() - 50.0, rand100() - 50.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Vec3       v01;
        while (v01.squaredNorm2()<=EPSILON_SQUARED){v01=randVec3D100();}
        
        Quaternion Q02(v01, M_PI*(rand100()/50.0 - 1.0));
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();
        finder01.constructVertexAttributes();
        finder01.findDimension();

        EXPECT_EQ(finder01.mVertexAttributes.size(), 0);
        EXPECT_EQ(finder01.mDimension, -1);
    }
}


/**  @brief constructVertexAttributes() and findDimension();
 *          from test08
 */
TEST_F(IntersectionFinderTest, Test18) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(0.0, 0.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();
        finder01.constructVertexAttributes();
        finder01.findDimension();

        IntersectionFinder::Attributes att_01;
        att_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        att_01.mPred = IF_VERTEX_VERTEX;
        att_01.mVit1 = findVertex(m01, Vec3(1.0, 1.0, 1.0));
        att_01.mVit2 = findVertex(m02, Vec3(1.0, 1.0, 1.0));

        IntersectionFinder::Attributes att_02;
        att_02.mP    = rotQmat01 * Vec3(1.0, 1.0, -1.0) + trans01;
        att_02.mPred = IF_VERTEX_VERTEX;
        att_02.mVit1 = findVertex(m01, Vec3(1.0, 1.0,-1.0));
        att_02.mVit2 = findVertex(m02, Vec3(1.0, 1.0,-1.0));

        IntersectionFinder::Attributes att_03;
        att_03.mP    = rotQmat01 * Vec3(1.0, -1.0, 1.0) + trans01;
        att_03.mPred = IF_VERTEX_VERTEX;
        att_03.mVit1 = findVertex(m01, Vec3(1.0, -1.0, 1.0));
        att_03.mVit2 = findVertex(m02, Vec3(1.0, -1.0, 1.0));

        IntersectionFinder::Attributes att_04;
        att_04.mP    = rotQmat01 * Vec3(1.0, -1.0, -1.0) + trans01;
        att_04.mPred = IF_VERTEX_VERTEX;
        att_04.mVit1 = findVertex(m01, Vec3(1.0, -1.0, -1.0));
        att_04.mVit2 = findVertex(m02, Vec3(1.0, -1.0, -1.0));

        IntersectionFinder::Attributes att_05;
        att_05.mP    = rotQmat01 * Vec3(-1.0, 1.0, 1.0) + trans01;
        att_05.mPred = IF_VERTEX_VERTEX;
        att_05.mVit1 = findVertex(m01, Vec3(-1.0, 1.0, 1.0));
        att_05.mVit2 = findVertex(m02, Vec3(-1.0, 1.0, 1.0));

        IntersectionFinder::Attributes att_06;
        att_06.mP    = rotQmat01 * Vec3(-1.0, 1.0, -1.0) + trans01;
        att_06.mPred = IF_VERTEX_VERTEX;
        att_06.mVit1 = findVertex(m01, Vec3(-1.0, 1.0,-1.0));
        att_06.mVit2 = findVertex(m02, Vec3(-1.0, 1.0,-1.0));

        IntersectionFinder::Attributes att_07;
        att_07.mP    = rotQmat01 * Vec3(-1.0, -1.0, 1.0) + trans01;
        att_07.mPred = IF_VERTEX_VERTEX;
        att_07.mVit1 = findVertex(m01, Vec3(-1.0, -1.0, 1.0));
        att_07.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, 1.0));

        IntersectionFinder::Attributes att_08;
        att_08.mP    = rotQmat01 * Vec3(-1.0, -1.0, -1.0) + trans01;
        att_08.mPred = IF_VERTEX_VERTEX;
        att_08.mVit1 = findVertex(m01, Vec3(-1.0, -1.0, -1.0));
        att_08.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> atts_01;
        atts_01.push_back(att_01);
        atts_01.push_back(att_02);
        atts_01.push_back(att_03);
        atts_01.push_back(att_04);
        atts_01.push_back(att_05);
        atts_01.push_back(att_06);
        atts_01.push_back(att_07);
        atts_01.push_back(att_08);

        EXPECT_EQ(
             compareAttributes(finder01.mVertexAttributes,atts_01),true);
        EXPECT_EQ(finder01.mDimension, 3);        

        finder01.processThreeDimensionalIntersection();

    }
}


/**  @brief constructVertexAttributes() and findDimension()
 *    from test09
 */  
TEST_F(IntersectionFinderTest, Test19) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 0.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();
        finder01.constructVertexAttributes();
        finder01.findDimension();

        IntersectionFinder::Attributes att_01;
        att_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        att_01.mPred = IF_VERTEX_VERTEX;
        att_01.mVit1 = findVertex(m01, Vec3(1.0, 1.0, 1.0));
        att_01.mVit2 = findVertex(m02, Vec3(-1.0, 1.0, 1.0));

        IntersectionFinder::Attributes att_02;
        att_02.mP    = rotQmat01 * Vec3(1.0, 1.0, -1.0) + trans01;
        att_02.mPred = IF_VERTEX_VERTEX;
        att_02.mVit1 = findVertex(m01, Vec3(1.0, 1.0,-1.0));
        att_02.mVit2 = findVertex(m02, Vec3(-1.0, 1.0,-1.0));

        IntersectionFinder::Attributes att_03;
        att_03.mP    = rotQmat01 * Vec3(1.0, -1.0, 1.0) + trans01;
        att_03.mPred = IF_VERTEX_VERTEX;
        att_03.mVit1 = findVertex(m01, Vec3(1.0, -1.0, 1.0));
        att_03.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, 1.0));

        IntersectionFinder::Attributes att_04;
        att_04.mP    = rotQmat01 * Vec3(1.0, -1.0, -1.0) + trans01;
        att_04.mPred = IF_VERTEX_VERTEX;
        att_04.mVit1 = findVertex(m01, Vec3(1.0, -1.0, -1.0));
        att_04.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> atts_01;
        atts_01.push_back(att_01);
        atts_01.push_back(att_02);
        atts_01.push_back(att_03);
        atts_01.push_back(att_04);

        EXPECT_EQ(
             compareAttributes(finder01.mVertexAttributes,atts_01),true);
        EXPECT_EQ(finder01.mDimension, 2);
    }
}


/**  @brief  constructVertexAttributes() and findDimension()
 *          from test10
 */
TEST_F(IntersectionFinderTest, Test20) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 1.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();
        finder01.constructVertexAttributes();
        finder01.findDimension();

        IntersectionFinder::Attributes att_01;
        att_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        att_01.mPred = IF_VERTEX_EDGE;
        att_01.mVit1 = findVertex(m01, Vec3( 1.0,  1.0, 1.0));
        att_01.mEit2 = findEdge  (m02, Vec3(-1.0,  1.0, 1.0),
                                       Vec3(-1.0, -1.0, 1.0));

        IntersectionFinder::Attributes att_02;
        att_02.mP    = rotQmat01 * Vec3(1.0, 1.0, -1.0) + trans01;
        att_02.mPred = IF_VERTEX_EDGE;
        att_02.mVit1 = findVertex(m01, Vec3( 1.0,  1.0, -1.0));
        att_02.mEit2 = findEdge  (m02, Vec3(-1.0,  1.0, -1.0),
                                       Vec3(-1.0, -1.0, -1.0));

        IntersectionFinder::Attributes att_03;
        att_03.mP    = rotQmat01 * Vec3( 1.0, 0.0, 1.0) + trans01;
        att_03.mPred = IF_EDGE_VERTEX;
        att_03.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, 1.0),
                                       Vec3( 1.0, -1.0, 1.0));
        att_03.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, 1.0));

        IntersectionFinder::Attributes att_04;
        att_04.mP    = rotQmat01 * Vec3( 1.0, 0.0, -1.0) + trans01;
        att_04.mPred = IF_EDGE_VERTEX;
        att_04.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, -1.0),
                                       Vec3( 1.0, -1.0, -1.0));
        att_04.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> atts_01;
        atts_01.push_back(att_01);
        atts_01.push_back(att_02);
        atts_01.push_back(att_03);
        atts_01.push_back(att_04);

        EXPECT_EQ(
             compareAttributes(finder01.mVertexAttributes,atts_01),true);
        EXPECT_EQ(finder01.mDimension, 2);
    }
}


/**  @brief constructVertexAttributes() and findDimension()
 *          from test11
 */
TEST_F(IntersectionFinderTest, Test21) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 1.0, 1.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();
        finder01.constructVertexAttributes();
        finder01.findDimension();

        IntersectionFinder::Attributes att_01;
        att_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        att_01.mPred = IF_VERTEX_FACE;
        att_01.mVit1 = findVertex(m01, Vec3( 1.0,  1.0, 1.0));
        att_01.mFit2 = findFace  (m02, Vec3(-1.0,  0.0, 0.0));

        IntersectionFinder::Attributes att_02;
        att_02.mP    = rotQmat01 * Vec3( 1.0,  0.0,  0.0) + trans01;
        att_02.mPred = IF_FACE_VERTEX;
        att_02.mFit1 = findFace  (m01, Vec3( 1.0,  0.0,  0.0));
        att_02.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        IntersectionFinder::Attributes att_03;
        att_03.mP    = rotQmat01 * Vec3( 1.0, 0.0, 1.0) + trans01;
        att_03.mPred = IF_EDGE_EDGE;
        att_03.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, 1.0),
                                       Vec3( 1.0, -1.0, 1.0));
        att_03.mEit2 = findEdge  (m02, Vec3(-1.0, -1.0, 1.0),
                                       Vec3(-1.0, -1.0,-1.0));

        IntersectionFinder::Attributes att_04;
        att_04.mP    = rotQmat01 * Vec3( 1.0, 1.0, 0.0) + trans01;
        att_04.mPred = IF_EDGE_EDGE;
        att_04.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0,  1.0),
                                       Vec3( 1.0,  1.0, -1.0));
        att_04.mEit2 = findEdge  (m02, Vec3(-1.0,  1.0, -1.0),
                                       Vec3(-1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> atts_01;
        atts_01.push_back(att_01);
        atts_01.push_back(att_02);
        atts_01.push_back(att_03);
        atts_01.push_back(att_04);

        EXPECT_EQ(
             compareAttributes(atts_01,finder01.mVertexAttributes),true);
        EXPECT_EQ(finder01.mDimension, 2);
    }
}


/**  @brief constructVertexAttributes() and findDimension()
 *          from test12
 */
TEST_F(IntersectionFinderTest, Test22) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(1.0, 1.0, 1.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();
        finder01.constructVertexAttributes();
        finder01.findDimension();

        IntersectionFinder::Attributes att_01;
        att_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 0.0) + trans01;
        att_01.mPred = IF_EDGE_FACE;
        att_01.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, 1.0),
                                       Vec3( 1.0,  1.0,-1.0));
        att_01.mFit2 = findFace  (m02, Vec3( 0.0,  0.0, -1.0));

        IntersectionFinder::Attributes att_02;
        att_02.mP    = rotQmat01 * Vec3(1.0, 0.0, 1.0) + trans01;
        att_02.mPred = IF_EDGE_FACE;
        att_02.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, 1.0),
                                       Vec3( 1.0, -1.0, 1.0));
        att_02.mFit2 = findFace  (m02, Vec3( 0.0, -1.0, 0.0));

        IntersectionFinder::Attributes att_03;
        att_03.mP    = rotQmat01 * Vec3(0.0, 1.0, 1.0) + trans01;
        att_03.mPred = IF_EDGE_FACE;
        att_03.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, 1.0),
                                       Vec3(-1.0,  1.0, 1.0));
        att_03.mFit2 = findFace  (m02, Vec3(-1.0, 0.0, 0.0));

        IntersectionFinder::Attributes att_04;
        att_04.mP    = rotQmat01 * Vec3(1.0, 0.0, 0.0) + trans01;
        att_04.mPred = IF_FACE_EDGE;
        att_04.mFit1 = findFace  (m01, Vec3( 1.0,  0.0, 0.0));
        att_04.mEit2 = findEdge  (m02, Vec3( 1.0, -1.0,-1.0),
                                       Vec3(-1.0, -1.0,-1.0));

        IntersectionFinder::Attributes att_05;
        att_05.mP    = rotQmat01 * Vec3(0.0, 1.0, 0.0) + trans01;
        att_05.mPred = IF_FACE_EDGE;
        att_05.mFit1 = findFace  (m01, Vec3( 0.0,  1.0, 0.0));
        att_05.mEit2 = findEdge  (m02, Vec3(-1.0,  1.0,-1.0),
                                       Vec3(-1.0, -1.0,-1.0));

        IntersectionFinder::Attributes att_06;
        att_06.mP    = rotQmat01 * Vec3(0.0, 0.0, 1.0) + trans01;
        att_06.mPred = IF_FACE_EDGE;
        att_06.mFit1 = findFace  (m01, Vec3( 0.0,  0.0, 1.0));
        att_06.mEit2 = findEdge  (m02, Vec3(-1.0, -1.0, 1.0),
                                       Vec3(-1.0, -1.0,-1.0));

        IntersectionFinder::Attributes att_07;
        att_07.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        att_07.mPred = IF_VERTEX_INTERIOR;
        att_07.mVit1 = findVertex(m01, Vec3( 1.0, 1.0, 1.0));

        IntersectionFinder::Attributes att_08;
        att_08.mP    = rotQmat01 * Vec3(0.0, 0.0, 0.0) + trans01;
        att_08.mPred = IF_INTERIOR_VERTEX;
        att_08.mVit2 = findVertex(m02, Vec3( -1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> atts_01;
        atts_01.push_back(att_01);
        atts_01.push_back(att_02);
        atts_01.push_back(att_03);
        atts_01.push_back(att_04);
        atts_01.push_back(att_05);
        atts_01.push_back(att_06);
        atts_01.push_back(att_07);
        atts_01.push_back(att_08);

        EXPECT_EQ(
             compareAttributes(atts_01,finder01.mVertexAttributes),true);
        EXPECT_EQ(finder01.mDimension, 3);
    }
}


/**  @brief constructVertexAttributes() and findDimension()
 *          from test13
 */
TEST_F(IntersectionFinderTest, Test23) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(1.0, 0.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();
        finder01.constructVertexAttributes();
        finder01.findDimension();

        IntersectionFinder::Attributes att_01;
        att_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        att_01.mPred = IF_VERTEX_EDGE;
        att_01.mVit1 = findVertex(m01, Vec3( 1.0,  1.0,  1.0));
        att_01.mEit2 = findEdge  (m02, Vec3( 1.0,  1.0,  1.0),
                                       Vec3(-1.0,  1.0,  1.0));

        IntersectionFinder::Attributes att_02;
        att_02.mP    = rotQmat01 * Vec3(1.0, 1.0, -1.0) + trans01;
        att_02.mPred = IF_VERTEX_EDGE;
        att_02.mVit1 = findVertex(m01, Vec3( 1.0,  1.0, -1.0));
        att_02.mEit2 = findEdge  (m02, Vec3( 1.0,  1.0, -1.0),
                                       Vec3(-1.0,  1.0, -1.0));

        IntersectionFinder::Attributes att_03;
        att_03.mP    = rotQmat01 * Vec3(1.0, -1.0, -1.0) + trans01;
        att_03.mPred = IF_VERTEX_EDGE;
        att_03.mVit1 = findVertex(m01, Vec3( 1.0,  -1.0, -1.0));
        att_03.mEit2 = findEdge  (m02, Vec3( 1.0,  -1.0, -1.0),
                                       Vec3(-1.0,  -1.0, -1.0));

        IntersectionFinder::Attributes att_04;
        att_04.mP    = rotQmat01 * Vec3(1.0, -1.0,  1.0) + trans01;
        att_04.mPred = IF_VERTEX_EDGE;
        att_04.mVit1 = findVertex(m01, Vec3( 1.0,  -1.0,  1.0));
        att_04.mEit2 = findEdge  (m02, Vec3( 1.0,  -1.0,  1.0),
                                       Vec3(-1.0,  -1.0,  1.0));

        IntersectionFinder::Attributes att_05;
        att_05.mP    = rotQmat01 * Vec3(0.0, 1.0, 1.0) + trans01;
        att_05.mPred = IF_EDGE_VERTEX;
        att_05.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0,  1.0),
                                       Vec3(-1.0,  1.0,  1.0));
        att_05.mVit2 = findVertex(m02, Vec3(-1.0,  1.0,  1.0));

        IntersectionFinder::Attributes att_06;
        att_06.mP    = rotQmat01 * Vec3(0.0, 1.0, -1.0) + trans01;
        att_06.mPred = IF_EDGE_VERTEX;
        att_06.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, -1.0),
                                       Vec3(-1.0,  1.0, -1.0));
        att_06.mVit2 = findVertex(m02, Vec3(-1.0,  1.0, -1.0));

        IntersectionFinder::Attributes att_07;
        att_07.mP    = rotQmat01 * Vec3(0.0, -1.0, -1.0) + trans01;
        att_07.mPred = IF_EDGE_VERTEX;
        att_07.mEit1 = findEdge  (m01, Vec3( 1.0, -1.0, -1.0),
                                       Vec3(-1.0, -1.0, -1.0));
        att_07.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        IntersectionFinder::Attributes att_08;
        att_08.mP    = rotQmat01 * Vec3(0.0, -1.0, 1.0) + trans01;
        att_08.mPred = IF_EDGE_VERTEX;
        att_08.mEit1 = findEdge  (m01, Vec3( 1.0, -1.0, 1.0),
                                       Vec3(-1.0, -1.0, 1.0));
        att_08.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, 1.0));

        vector<IntersectionFinder::Attributes> atts_01;
        atts_01.push_back(att_01);
        atts_01.push_back(att_02);
        atts_01.push_back(att_03);
        atts_01.push_back(att_04);
        atts_01.push_back(att_05);
        atts_01.push_back(att_06);
        atts_01.push_back(att_07);
        atts_01.push_back(att_08);

        EXPECT_EQ(
             compareAttributes(atts_01,finder01.mVertexAttributes),true);
        EXPECT_EQ(finder01.mDimension, 3);
    }
}


/**  @brief constructVertexAttributes() and findDimension()
 *          from test14
 */
TEST_F(IntersectionFinderTest, Test24) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 2.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();
        finder01.constructVertexAttributes();
        finder01.findDimension();

        IntersectionFinder::Attributes att_01;
        att_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        att_01.mPred = IF_VERTEX_VERTEX;
        att_01.mVit1 = findVertex(m01, Vec3( 1.0,  1.0,  1.0));
        att_01.mVit2 = findVertex(m02, Vec3(-1.0, -1.0,  1.0));

        IntersectionFinder::Attributes att_02;
        att_02.mP    = rotQmat01 * Vec3(1.0, 1.0, -1.0) + trans01;
        att_02.mPred = IF_VERTEX_VERTEX;
        att_02.mVit1 = findVertex(m01, Vec3( 1.0,  1.0, -1.0));
        att_02.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> atts_01;
        atts_01.push_back(att_01);
        atts_01.push_back(att_02);

        EXPECT_EQ(
             compareAttributes(atts_01,finder01.mVertexAttributes),true);
        EXPECT_EQ(finder01.mDimension, 1);
    }
}


/**  @brief constructVertexAttributes() and findDimension()
 *          from test15
 */
TEST_F(IntersectionFinderTest, Test25) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 2.0, 1.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();
        finder01.constructVertexAttributes();
        finder01.findDimension();

        IntersectionFinder::Attributes att_01;
        att_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        att_01.mPred = IF_VERTEX_EDGE;
        att_01.mVit1 = findVertex(m01, Vec3( 1.0,  1.0, 1.0));
        att_01.mEit2 = findEdge  (m02, Vec3(-1.0, -1.0,  1.0),
                                       Vec3(-1.0, -1.0, -1.0));

        IntersectionFinder::Attributes att_02;
        att_02.mP    = rotQmat01 * Vec3(1.0, 1.0,  0.0) + trans01;
        att_02.mPred = IF_EDGE_VERTEX;
        att_02.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0,  1.0),
                                       Vec3( 1.0,  1.0, -1.0));
        att_02.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> atts_01;
        atts_01.push_back(att_01);
        atts_01.push_back(att_02);

        EXPECT_EQ(
             compareAttributes(atts_01,finder01.mVertexAttributes),true);
        EXPECT_EQ(finder01.mDimension, 1);
    }
}


/**  @brief constructVertexAttributes() and findDimension()
 *          from test16
 */
TEST_F(IntersectionFinderTest, Test26) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 2.0, 2.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.prepareVertices();
        finder01.prepareEdgesAndFaces();
        finder01.findIntersectionsEdgesOn1FacesOn2();
        finder01.findIntersectionsFacesOn1EdgesOn2();
        finder01.findRemainingInteriorVertices();
        finder01.constructVertexAttributes();
        finder01.findDimension();

        IntersectionFinder::Attributes att_01;
        att_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        att_01.mPred = IF_VERTEX_VERTEX;
        att_01.mVit1 = findVertex(m01, Vec3( 1.0,  1.0,  1.0));
        att_01.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> atts_01;
        atts_01.push_back(att_01);

        EXPECT_EQ(
             compareAttributes(atts_01,finder01.mVertexAttributes),true);
        EXPECT_EQ(finder01.mDimension, 0);
    }
}


/**  @brief find();
 *          from test08
 */
TEST_F(IntersectionFinderTest, Test27) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(0.0, 0.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.find(sepAxis);

        IntersectionFinder::Attributes vatt_01;
        vatt_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        vatt_01.mPred = IF_VERTEX_VERTEX;
        vatt_01.mVit1 = findVertex(m01, Vec3(1.0, 1.0, 1.0));
        vatt_01.mVit2 = findVertex(m02, Vec3(1.0, 1.0, 1.0));

        IntersectionFinder::Attributes vatt_02;
        vatt_02.mP    = rotQmat01 * Vec3(1.0, 1.0, -1.0) + trans01;
        vatt_02.mPred = IF_VERTEX_VERTEX;
        vatt_02.mVit1 = findVertex(m01, Vec3(1.0, 1.0,-1.0));
        vatt_02.mVit2 = findVertex(m02, Vec3(1.0, 1.0,-1.0));

        IntersectionFinder::Attributes vatt_03;
        vatt_03.mP    = rotQmat01 * Vec3(1.0, -1.0, 1.0) + trans01;
        vatt_03.mPred = IF_VERTEX_VERTEX;
        vatt_03.mVit1 = findVertex(m01, Vec3(1.0, -1.0, 1.0));
        vatt_03.mVit2 = findVertex(m02, Vec3(1.0, -1.0, 1.0));

        IntersectionFinder::Attributes vatt_04;
        vatt_04.mP    = rotQmat01 * Vec3(1.0, -1.0, -1.0) + trans01;
        vatt_04.mPred = IF_VERTEX_VERTEX;
        vatt_04.mVit1 = findVertex(m01, Vec3(1.0, -1.0, -1.0));
        vatt_04.mVit2 = findVertex(m02, Vec3(1.0, -1.0, -1.0));

        IntersectionFinder::Attributes vatt_05;
        vatt_05.mP    = rotQmat01 * Vec3(-1.0, 1.0, 1.0) + trans01;
        vatt_05.mPred = IF_VERTEX_VERTEX;
        vatt_05.mVit1 = findVertex(m01, Vec3(-1.0, 1.0, 1.0));
        vatt_05.mVit2 = findVertex(m02, Vec3(-1.0, 1.0, 1.0));

        IntersectionFinder::Attributes vatt_06;
        vatt_06.mP    = rotQmat01 * Vec3(-1.0, 1.0, -1.0) + trans01;
        vatt_06.mPred = IF_VERTEX_VERTEX;
        vatt_06.mVit1 = findVertex(m01, Vec3(-1.0, 1.0,-1.0));
        vatt_06.mVit2 = findVertex(m02, Vec3(-1.0, 1.0,-1.0));

        IntersectionFinder::Attributes vatt_07;
        vatt_07.mP    = rotQmat01 * Vec3(-1.0, -1.0, 1.0) + trans01;
        vatt_07.mPred = IF_VERTEX_VERTEX;
        vatt_07.mVit1 = findVertex(m01, Vec3(-1.0, -1.0, 1.0));
        vatt_07.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, 1.0));

        IntersectionFinder::Attributes vatt_08;
        vatt_08.mP    = rotQmat01 * Vec3(-1.0, -1.0, -1.0) + trans01;
        vatt_08.mPred = IF_VERTEX_VERTEX;
        vatt_08.mVit1 = findVertex(m01, Vec3(-1.0, -1.0, -1.0));
        vatt_08.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> vatts_01;
        vatts_01.push_back(vatt_01);
        vatts_01.push_back(vatt_02);
        vatts_01.push_back(vatt_03);
        vatts_01.push_back(vatt_04);
        vatts_01.push_back(vatt_05);
        vatts_01.push_back(vatt_06);
        vatts_01.push_back(vatt_07);
        vatts_01.push_back(vatt_08);

        vector<pair<Vec3, Vec3> > ePoints;
        vector<IntersectionFinder::Attributes> eatts_01;

        auto eps_01 = make_pair(rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01,
                                rotQmat01 * Vec3(-1.0, 1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_01;
        eatt_01.mPred = IF_EDGE_EDGE;
        eatt_01.mEit1 = findEdge(m01, Vec3( 1.0,  1.0,  1.0),
                                      Vec3(-1.0,  1.0,  1.0));
        eatt_01.mEit2 = findEdge(m02, Vec3( 1.0,  1.0,  1.0),
                                      Vec3(-1.0,  1.0,  1.0));

        auto eps_02 = make_pair(rotQmat01 * Vec3(-1.0, 1.0, 1.0) + trans01,
                                rotQmat01 * Vec3(-1.0,-1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_02;
        eatt_02.mPred = IF_EDGE_EDGE;
        eatt_02.mEit1 = findEdge(m01, Vec3(-1.0,  1.0,  1.0),
                                      Vec3(-1.0, -1.0,  1.0));
        eatt_02.mEit2 = findEdge(m02, Vec3(-1.0,  1.0,  1.0),
                                      Vec3(-1.0, -1.0,  1.0));

        auto eps_03 = make_pair(rotQmat01 * Vec3(-1.0,-1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0,-1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_03;
        eatt_03.mPred = IF_EDGE_EDGE;
        eatt_03.mEit1 = findEdge(m01, Vec3(-1.0, -1.0,  1.0),
                                      Vec3( 1.0, -1.0,  1.0));
        eatt_03.mEit2 = findEdge(m02, Vec3(-1.0, -1.0,  1.0),
                                      Vec3( 1.0, -1.0,  1.0));

        auto eps_04 = make_pair(rotQmat01 * Vec3( 1.0,-1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_04;
        eatt_04.mPred = IF_EDGE_EDGE;
        eatt_04.mEit1 = findEdge(m01, Vec3( 1.0, -1.0,  1.0),
                                      Vec3( 1.0,  1.0,  1.0));
        eatt_04.mEit2 = findEdge(m02, Vec3( 1.0, -1.0,  1.0),
                                      Vec3( 1.0,  1.0,  1.0));

        auto eps_05 = make_pair(rotQmat01 * Vec3( 1.0, 1.0,  1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0, -1.0) + trans01 );

        IntersectionFinder::Attributes eatt_05;
        eatt_05.mPred = IF_EDGE_EDGE;
        eatt_05.mEit1 = findEdge(m01, Vec3( 1.0,  1.0,  1.0),
                                      Vec3( 1.0,  1.0, -1.0));
        eatt_05.mEit2 = findEdge(m02, Vec3( 1.0,  1.0,  1.0),
                                      Vec3( 1.0,  1.0, -1.0));

        auto eps_06 = make_pair(rotQmat01 * Vec3( 1.0,-1.0,  1.0) + trans01,
                                rotQmat01 * Vec3( 1.0,-1.0, -1.0) + trans01 );

        IntersectionFinder::Attributes eatt_06;
        eatt_06.mPred = IF_EDGE_EDGE;
        eatt_06.mEit1 = findEdge(m01, Vec3( 1.0, -1.0,  1.0),
                                      Vec3( 1.0, -1.0, -1.0));
        eatt_06.mEit2 = findEdge(m02, Vec3( 1.0, -1.0,  1.0),
                                      Vec3( 1.0, -1.0, -1.0));

        auto eps_07 = make_pair(rotQmat01 * Vec3(-1.0,-1.0,  1.0) + trans01,
                                rotQmat01 * Vec3(-1.0,-1.0, -1.0) + trans01 );

        IntersectionFinder::Attributes eatt_07;
        eatt_07.mPred = IF_EDGE_EDGE;
        eatt_07.mEit1 = findEdge(m01, Vec3(-1.0, -1.0,  1.0),
                                      Vec3(-1.0, -1.0, -1.0));
        eatt_07.mEit2 = findEdge(m02, Vec3(-1.0, -1.0,  1.0),
                                      Vec3(-1.0, -1.0, -1.0));

        auto eps_08 = make_pair(rotQmat01 * Vec3(-1.0, 1.0,  1.0) + trans01,
                                rotQmat01 * Vec3(-1.0, 1.0, -1.0) + trans01 );

        IntersectionFinder::Attributes eatt_08;
        eatt_08.mPred = IF_EDGE_EDGE;
        eatt_08.mEit1 = findEdge(m01, Vec3(-1.0,  1.0,  1.0),
                                      Vec3(-1.0,  1.0, -1.0));
        eatt_08.mEit2 = findEdge(m02, Vec3(-1.0,  1.0,  1.0),
                                      Vec3(-1.0,  1.0, -1.0));

        auto eps_09 = make_pair(rotQmat01 * Vec3( 1.0, 1.0, -1.0) + trans01,
                                rotQmat01 * Vec3(-1.0, 1.0, -1.0) + trans01 );

        IntersectionFinder::Attributes eatt_09;
        eatt_09.mPred = IF_EDGE_EDGE;
        eatt_09.mEit1 = findEdge(m01, Vec3( 1.0,  1.0, -1.0),
                                      Vec3(-1.0,  1.0, -1.0));
        eatt_09.mEit2 = findEdge(m02, Vec3( 1.0,  1.0, -1.0),
                                      Vec3(-1.0,  1.0, -1.0));

        auto eps_10 = make_pair(rotQmat01 * Vec3(-1.0, 1.0, -1.0) + trans01,
                                rotQmat01 * Vec3(-1.0,-1.0, -1.0) + trans01 );

        IntersectionFinder::Attributes eatt_10;
        eatt_10.mPred = IF_EDGE_EDGE;
        eatt_10.mEit1 = findEdge(m01, Vec3(-1.0,  1.0,  -1.0),
                                      Vec3(-1.0, -1.0,  -1.0));
        eatt_10.mEit2 = findEdge(m02, Vec3(-1.0,  1.0,  -1.0),
                                      Vec3(-1.0, -1.0,  -1.0));

        auto eps_11 = make_pair(rotQmat01 * Vec3(-1.0,-1.0, -1.0) + trans01,
                                rotQmat01 * Vec3( 1.0,-1.0, -1.0) + trans01 );

        IntersectionFinder::Attributes eatt_11;
        eatt_11.mPred = IF_EDGE_EDGE;
        eatt_11.mEit1 = findEdge(m01, Vec3(-1.0, -1.0,  -1.0),
                                      Vec3( 1.0, -1.0,  -1.0));
        eatt_11.mEit2 = findEdge(m02, Vec3(-1.0, -1.0,  -1.0),
                                      Vec3( 1.0, -1.0,  -1.0));

        auto eps_12 = make_pair(rotQmat01 * Vec3( 1.0,-1.0, -1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0, -1.0) + trans01 );

        IntersectionFinder::Attributes eatt_12;
        eatt_12.mPred = IF_EDGE_EDGE;
        eatt_12.mEit1 = findEdge(m01, Vec3( 1.0, -1.0,  -1.0),
                                      Vec3( 1.0,  1.0,  -1.0));
        eatt_12.mEit2 = findEdge(m02, Vec3( 1.0, -1.0,  -1.0),
                                      Vec3( 1.0,  1.0,  -1.0));

        ePoints.push_back(eps_01);
        ePoints.push_back(eps_02);
        ePoints.push_back(eps_03);
        ePoints.push_back(eps_04);
        ePoints.push_back(eps_05);
        ePoints.push_back(eps_06);
        ePoints.push_back(eps_07);
        ePoints.push_back(eps_08);
        ePoints.push_back(eps_09);
        ePoints.push_back(eps_10);
        ePoints.push_back(eps_11);
        ePoints.push_back(eps_12);

        eatts_01.push_back(eatt_01);
        eatts_01.push_back(eatt_02);
        eatts_01.push_back(eatt_03);
        eatts_01.push_back(eatt_04);
        eatts_01.push_back(eatt_05);
        eatts_01.push_back(eatt_06);
        eatts_01.push_back(eatt_07);
        eatts_01.push_back(eatt_08);
        eatts_01.push_back(eatt_09);
        eatts_01.push_back(eatt_10);
        eatts_01.push_back(eatt_11);
        eatts_01.push_back(eatt_12);

        IntersectionFinder::Attributes fatt_01;
        fatt_01.mP    = rotQmat01 * Vec3( 1.0, 0.0, 0.0);
        fatt_01.mPred = IF_FACE_FACE;
        fatt_01.mFit1 = findFace(m01, Vec3( 1.0, 0.0, 0.0));
        fatt_01.mFit2 = findFace(m02, Vec3( 1.0, 0.0, 0.0));

        IntersectionFinder::Attributes fatt_02;
        fatt_02.mP    = rotQmat01 * Vec3( 0.0, 1.0, 0.0);
        fatt_02.mPred = IF_FACE_FACE;
        fatt_02.mFit1 = findFace(m01, Vec3( 0.0, 1.0, 0.0));
        fatt_02.mFit2 = findFace(m02, Vec3( 0.0, 1.0, 0.0));

        IntersectionFinder::Attributes fatt_03;
        fatt_03.mP    = rotQmat01 * Vec3( 0.0, 0.0, 1.0);
        fatt_03.mPred = IF_FACE_FACE;
        fatt_03.mFit1 = findFace(m01, Vec3( 0.0, 0.0, 1.0));
        fatt_03.mFit2 = findFace(m02, Vec3( 0.0, 0.0, 1.0));

        IntersectionFinder::Attributes fatt_04;
        fatt_04.mP    = rotQmat01 * Vec3(-1.0, 0.0, 0.0);
        fatt_04.mPred = IF_FACE_FACE;
        fatt_04.mFit1 = findFace(m01, Vec3(-1.0, 0.0, 0.0));
        fatt_04.mFit2 = findFace(m02, Vec3(-1.0, 0.0, 0.0));

        IntersectionFinder::Attributes fatt_05;
        fatt_05.mP    = rotQmat01 * Vec3( 0.0, -1.0, 0.0);
        fatt_05.mPred = IF_FACE_FACE;
        fatt_05.mFit1 = findFace(m01, Vec3( 0.0, -1.0, 0.0));
        fatt_05.mFit2 = findFace(m02, Vec3( 0.0, -1.0, 0.0));

        IntersectionFinder::Attributes fatt_06;
        fatt_06.mP    = rotQmat01 * Vec3( 0.0, 0.0, -1.0);
        fatt_06.mPred = IF_FACE_FACE;
        fatt_06.mFit1 = findFace(m01, Vec3( 0.0, 0.0, -1.0));
        fatt_06.mFit2 = findFace(m02, Vec3( 0.0, 0.0, -1.0));

        vector<IntersectionFinder::Attributes> fatts_01;
        fatts_01.push_back(fatt_01);
        fatts_01.push_back(fatt_02);
        fatts_01.push_back(fatt_03);
        fatts_01.push_back(fatt_04);
        fatts_01.push_back(fatt_05);
        fatts_01.push_back(fatt_06);


        auto res = checkResults3D(
                       finder01, vatts_01, ePoints, eatts_01, fatts_01);
        EXPECT_EQ(res, true);

    }
}






/**  @brief constructVertexAttributes() and findDimension()
 *    from test09
 */  
TEST_F(IntersectionFinderTest, Test28) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 0.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;


        double epsilonZero    = EPSILON_LINEAR*1.0;
        double epsilonZeroPCA = EPSILON_LINEAR*1.0;
        double epsilonAngle   = EPSILON_LINEAR*1.0;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.find(sepAxis);

        IntersectionFinder::Attributes vatt_01;
        vatt_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        vatt_01.mPred = IF_VERTEX_VERTEX;
        vatt_01.mVit1 = findVertex(m01, Vec3(1.0, 1.0, 1.0));
        vatt_01.mVit2 = findVertex(m02, Vec3(-1.0, 1.0, 1.0));

        IntersectionFinder::Attributes vatt_02;
        vatt_02.mP    = rotQmat01 * Vec3(1.0, 1.0, -1.0) + trans01;
        vatt_02.mPred = IF_VERTEX_VERTEX;
        vatt_02.mVit1 = findVertex(m01, Vec3(1.0, 1.0,-1.0));
        vatt_02.mVit2 = findVertex(m02, Vec3(-1.0, 1.0,-1.0));

        IntersectionFinder::Attributes vatt_03;
        vatt_03.mP    = rotQmat01 * Vec3(1.0, -1.0, -1.0) + trans01;
        vatt_03.mPred = IF_VERTEX_VERTEX;
        vatt_03.mVit1 = findVertex(m01, Vec3(1.0, -1.0, -1.0));
        vatt_03.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        IntersectionFinder::Attributes vatt_04;
        vatt_04.mP    = rotQmat01 * Vec3(1.0, -1.0, 1.0) + trans01;
        vatt_04.mPred = IF_VERTEX_VERTEX;
        vatt_04.mVit1 = findVertex(m01, Vec3(1.0, -1.0, 1.0));
        vatt_04.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, 1.0));

        vector<IntersectionFinder::Attributes> vatts_01;
        vatts_01.push_back(vatt_01);
        vatts_01.push_back(vatt_02);
        vatts_01.push_back(vatt_03);
        vatts_01.push_back(vatt_04);

        vector<pair<Vec3, Vec3> > ePoints;
        vector<IntersectionFinder::Attributes> eatts_01;

        auto eps_01 = make_pair(rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0,-1.0) + trans01 );

        IntersectionFinder::Attributes eatt_01;
        eatt_01.mPred = IF_EDGE_EDGE;
        eatt_01.mEit1 = findEdge(m01, Vec3( 1.0,  1.0,  1.0),
                                      Vec3( 1.0,  1.0, -1.0));
        eatt_01.mEit2 = findEdge(m02, Vec3(-1.0,  1.0,  1.0),
                                      Vec3(-1.0,  1.0, -1.0));

        auto eps_02 = make_pair(rotQmat01 * Vec3( 1.0, 1.0,-1.0) + trans01,
                                rotQmat01 * Vec3( 1.0,-1.0,-1.0) + trans01 );

        IntersectionFinder::Attributes eatt_02;
        eatt_02.mPred = IF_EDGE_EDGE;
        eatt_02.mEit1 = findEdge(m01, Vec3( 1.0,  1.0, -1.0),
                                      Vec3( 1.0, -1.0, -1.0));
        eatt_02.mEit2 = findEdge(m02, Vec3(-1.0,  1.0, -1.0),
                                      Vec3(-1.0, -1.0, -1.0));

        auto eps_03 = make_pair(rotQmat01 * Vec3( 1.0,-1.0,-1.0) + trans01,
                                rotQmat01 * Vec3( 1.0,-1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_03;
        eatt_03.mPred = IF_EDGE_EDGE;
        eatt_03.mEit1 = findEdge(m01, Vec3( 1.0, -1.0, -1.0),
                                      Vec3( 1.0, -1.0,  1.0));
        eatt_03.mEit2 = findEdge(m02, Vec3(-1.0, -1.0, -1.0),
                                      Vec3(-1.0, -1.0,  1.0));

        auto eps_04 = make_pair(rotQmat01 * Vec3( 1.0,-1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_04;
        eatt_04.mPred = IF_EDGE_EDGE;
        eatt_04.mEit1 = findEdge(m01, Vec3( 1.0, -1.0,  1.0),
                                      Vec3( 1.0,  1.0,  1.0));
        eatt_04.mEit2 = findEdge(m02, Vec3(-1.0, -1.0,  1.0),
                                      Vec3(-1.0,  1.0,  1.0));

        ePoints.push_back(eps_01);
        ePoints.push_back(eps_02);
        ePoints.push_back(eps_03);
        ePoints.push_back(eps_04);

        eatts_01.push_back(eatt_01);
        eatts_01.push_back(eatt_02);
        eatts_01.push_back(eatt_03);
        eatts_01.push_back(eatt_04);

        IntersectionFinder::Attributes fatt_01;
        fatt_01.mP    = rotQmat01 * Vec3( 1.0, 0.0, 0.0);
        fatt_01.mPred = IF_FACE_FACE;
        fatt_01.mFit1 = findFace(m01, Vec3( 1.0, 0.0, 0.0));
        fatt_01.mFit2 = findFace(m02, Vec3(-1.0, 0.0, 0.0));


        auto res = checkResults2D(
                               finder01, vatts_01, ePoints, eatts_01, fatt_01);
        EXPECT_EQ(res, true);
    }
}


/**  @brief  constructVertexAttributes() and findDimension()
 *          from test10
 */
TEST_F(IntersectionFinderTest, Test29) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 1.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.find(sepAxis);

        IntersectionFinder::Attributes vatt_01;
        vatt_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        vatt_01.mPred = IF_VERTEX_EDGE;
        vatt_01.mVit1 = findVertex(m01, Vec3( 1.0,  1.0, 1.0));
        vatt_01.mEit2 = findEdge  (m02, Vec3(-1.0,  1.0, 1.0),
                                       Vec3(-1.0, -1.0, 1.0));

        IntersectionFinder::Attributes vatt_02;
        vatt_02.mP    = rotQmat01 * Vec3(1.0, 1.0, -1.0) + trans01;
        vatt_02.mPred = IF_VERTEX_EDGE;
        vatt_02.mVit1 = findVertex(m01, Vec3( 1.0,  1.0, -1.0));
        vatt_02.mEit2 = findEdge  (m02, Vec3(-1.0,  1.0, -1.0),
                                       Vec3(-1.0, -1.0, -1.0));

        IntersectionFinder::Attributes vatt_03;
        vatt_03.mP    = rotQmat01 * Vec3( 1.0, 0.0, -1.0) + trans01;
        vatt_03.mPred = IF_EDGE_VERTEX;
        vatt_03.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, -1.0),
                                       Vec3( 1.0, -1.0, -1.0));
        vatt_03.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        IntersectionFinder::Attributes vatt_04;
        vatt_04.mP    = rotQmat01 * Vec3( 1.0, 0.0, 1.0) + trans01;
        vatt_04.mPred = IF_EDGE_VERTEX;
        vatt_04.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, 1.0),
                                       Vec3( 1.0, -1.0, 1.0));
        vatt_04.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, 1.0));

        vector<IntersectionFinder::Attributes> vatts_01;
        vatts_01.push_back(vatt_01);
        vatts_01.push_back(vatt_02);
        vatts_01.push_back(vatt_03);
        vatts_01.push_back(vatt_04);

        vector<pair<Vec3, Vec3> > ePoints;
        vector<IntersectionFinder::Attributes> eatts_01;

        auto eps_01 = make_pair(rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0,-1.0) + trans01 );

        IntersectionFinder::Attributes eatt_01;
        eatt_01.mPred = IF_EDGE_FACE;
        eatt_01.mEit1 = findEdge(m01, Vec3( 1.0,  1.0,  1.0),
                                      Vec3( 1.0,  1.0, -1.0));
        eatt_01.mFit2 = findFace(m02, Vec3(-1.0,  0.0,  0.0));

        auto eps_02 = make_pair(rotQmat01 * Vec3( 1.0, 1.0,-1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 0.0,-1.0) + trans01 );

        IntersectionFinder::Attributes eatt_02;
        eatt_02.mPred = IF_EDGE_EDGE;
        eatt_02.mEit1 = findEdge(m01, Vec3( 1.0,  1.0, -1.0),
                                      Vec3( 1.0, -1.0, -1.0));
        eatt_02.mEit2 = findEdge(m02, Vec3(-1.0,  1.0, -1.0),
                                      Vec3(-1.0, -1.0, -1.0));

        auto eps_03 = make_pair(rotQmat01 * Vec3( 1.0, 0.0,-1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 0.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_03;
        eatt_03.mPred = IF_FACE_EDGE;
        eatt_03.mFit1 = findFace(m01, Vec3( 1.0,  0.0,  0.0));

        eatt_03.mEit2 = findEdge(m02, Vec3(-1.0, -1.0, -1.0),
                                      Vec3(-1.0, -1.0,  1.0));

        auto eps_04 = make_pair(rotQmat01 * Vec3( 1.0, 0.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_04;
        eatt_04.mPred = IF_EDGE_EDGE;
        eatt_04.mEit1 = findEdge(m01, Vec3( 1.0,  1.0,  1.0),
                                      Vec3( 1.0, -1.0,  1.0));
        eatt_04.mEit2 = findEdge(m02, Vec3(-1.0,  1.0,  1.0),
                                      Vec3(-1.0, -1.0,  1.0));

        ePoints.push_back(eps_01);
        ePoints.push_back(eps_02);
        ePoints.push_back(eps_03);
        ePoints.push_back(eps_04);

        eatts_01.push_back(eatt_01);
        eatts_01.push_back(eatt_02);
        eatts_01.push_back(eatt_03);
        eatts_01.push_back(eatt_04);

        IntersectionFinder::Attributes fatt_01;
        fatt_01.mP    = rotQmat01 * Vec3( 1.0, 0.0, 0.0);
        fatt_01.mPred = IF_FACE_FACE;
        fatt_01.mFit1 = findFace(m01, Vec3( 1.0, 0.0, 0.0));
        fatt_01.mFit2 = findFace(m02, Vec3(-1.0, 0.0, 0.0));


        auto res = checkResults2D(
                          finder01, vatts_01, ePoints, eatts_01, fatt_01);
        EXPECT_EQ(res, true);



    }
}


/**  @brief find()
 *          from test11
 */
TEST_F(IntersectionFinderTest, Test30) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 1.0, 1.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.find(sepAxis);

        IntersectionFinder::Attributes vatt_01;
        vatt_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        vatt_01.mPred = IF_VERTEX_FACE;
        vatt_01.mVit1 = findVertex(m01, Vec3( 1.0,  1.0, 1.0));
        vatt_01.mFit2 = findFace  (m02, Vec3(-1.0,  0.0, 0.0));


        IntersectionFinder::Attributes vatt_02;
        vatt_02.mP    = rotQmat01 * Vec3( 1.0, 0.0, 1.0) + trans01;
        vatt_02.mPred = IF_EDGE_EDGE;
        vatt_02.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, 1.0),
                                       Vec3( 1.0, -1.0, 1.0));
        vatt_02.mEit2 = findEdge  (m02, Vec3(-1.0, -1.0, 1.0),
                                       Vec3(-1.0, -1.0,-1.0));

        IntersectionFinder::Attributes vatt_03;
        vatt_03.mP    = rotQmat01 * Vec3( 1.0,  0.0,  0.0) + trans01;
        vatt_03.mPred = IF_FACE_VERTEX;
        vatt_03.mFit1 = findFace  (m01, Vec3( 1.0,  0.0,  0.0));
        vatt_03.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        IntersectionFinder::Attributes vatt_04;
        vatt_04.mP    = rotQmat01 * Vec3( 1.0, 1.0, 0.0) + trans01;
        vatt_04.mPred = IF_EDGE_EDGE;
        vatt_04.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0,  1.0),
                                       Vec3( 1.0,  1.0, -1.0));
        vatt_04.mEit2 = findEdge  (m02, Vec3(-1.0,  1.0, -1.0),
                                       Vec3(-1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> vatts_01;
        vatts_01.push_back(vatt_01);
        vatts_01.push_back(vatt_02);
        vatts_01.push_back(vatt_03);
        vatts_01.push_back(vatt_04);


        vector<pair<Vec3, Vec3> > ePoints;
        vector<IntersectionFinder::Attributes> eatts_01;

        auto eps_01 = make_pair(rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 0.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_01;
        eatt_01.mPred = IF_EDGE_FACE;
        eatt_01.mEit1 = findEdge(m01, Vec3( 1.0,  1.0,  1.0),
                                      Vec3( 1.0, -1.0,  1.0));
        eatt_01.mFit2 = findFace(m02, Vec3(-1.0,  0.0,  0.0));

        auto eps_02 = make_pair(rotQmat01 * Vec3( 1.0, 0.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 0.0, 0.0) + trans01 );

        IntersectionFinder::Attributes eatt_02;
        eatt_02.mPred = IF_FACE_EDGE;
        eatt_02.mFit1 = findFace(m01, Vec3( 1.0,  0.0,  0.0));
        eatt_02.mEit2 = findEdge(m02, Vec3(-1.0, -1.0,  1.0),
                                      Vec3(-1.0, -1.0, -1.0));

        auto eps_03 = make_pair(rotQmat01 * Vec3( 1.0, 0.0, 0.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0, 0.0) + trans01 );

        IntersectionFinder::Attributes eatt_03;
        eatt_03.mPred = IF_FACE_EDGE;
        eatt_03.mFit1 = findFace(m01, Vec3( 1.0,  0.0,  0.0));
        eatt_03.mEit2 = findEdge(m02, Vec3(-1.0,  1.0, -1.0),
                                      Vec3(-1.0, -1.0, -1.0));

        auto eps_04 = make_pair(rotQmat01 * Vec3( 1.0, 1.0, 0.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_04;
        eatt_04.mPred = IF_EDGE_FACE;
        eatt_04.mEit1 = findEdge(m01, Vec3( 1.0,  1.0,  1.0),
                                      Vec3( 1.0,  1.0, -1.0));
        eatt_04.mFit2 = findFace(m02, Vec3(-1.0,  0.0,  0.0));

        ePoints.push_back(eps_01);
        ePoints.push_back(eps_02);
        ePoints.push_back(eps_03);
        ePoints.push_back(eps_04);

        eatts_01.push_back(eatt_01);
        eatts_01.push_back(eatt_02);
        eatts_01.push_back(eatt_03);
        eatts_01.push_back(eatt_04);

        IntersectionFinder::Attributes fatt_01;
        fatt_01.mP    = rotQmat01 * Vec3( 1.0, 0.0, 0.0);
        fatt_01.mPred = IF_FACE_FACE;
        fatt_01.mFit1 = findFace(m01, Vec3( 1.0, 0.0, 0.0));
        fatt_01.mFit2 = findFace(m02, Vec3(-1.0, 0.0, 0.0));

        auto res = checkResults2D(
                          finder01, vatts_01, ePoints, eatts_01, fatt_01);
        EXPECT_EQ(res, true);

    }
}



/**  @brief find()
 *          from test12
 */
TEST_F(IntersectionFinderTest, Test31) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(1.0, 1.0, 1.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.find();

        IntersectionFinder::Attributes vatt_01;
        vatt_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 0.0) + trans01;
        vatt_01.mPred = IF_EDGE_FACE;
        vatt_01.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, 1.0),
                                       Vec3( 1.0,  1.0,-1.0));
        vatt_01.mFit2 = findFace  (m02, Vec3( 0.0,  0.0, -1.0));

        IntersectionFinder::Attributes vatt_02;
        vatt_02.mP    = rotQmat01 * Vec3(1.0, 0.0, 1.0) + trans01;
        vatt_02.mPred = IF_EDGE_FACE;
        vatt_02.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, 1.0),
                                       Vec3( 1.0, -1.0, 1.0));
        vatt_02.mFit2 = findFace  (m02, Vec3( 0.0, -1.0, 0.0));

        IntersectionFinder::Attributes vatt_03;
        vatt_03.mP    = rotQmat01 * Vec3(0.0, 1.0, 1.0) + trans01;
        vatt_03.mPred = IF_EDGE_FACE;
        vatt_03.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, 1.0),
                                       Vec3(-1.0,  1.0, 1.0));
        vatt_03.mFit2 = findFace  (m02, Vec3(-1.0, 0.0, 0.0));

        IntersectionFinder::Attributes vatt_04;
        vatt_04.mP    = rotQmat01 * Vec3(1.0, 0.0, 0.0) + trans01;
        vatt_04.mPred = IF_FACE_EDGE;
        vatt_04.mFit1 = findFace  (m01, Vec3( 1.0,  0.0, 0.0));
        vatt_04.mEit2 = findEdge  (m02, Vec3( 1.0, -1.0,-1.0),
                                       Vec3(-1.0, -1.0,-1.0));

        IntersectionFinder::Attributes vatt_05;
        vatt_05.mP    = rotQmat01 * Vec3(0.0, 1.0, 0.0) + trans01;
        vatt_05.mPred = IF_FACE_EDGE;
        vatt_05.mFit1 = findFace  (m01, Vec3( 0.0,  1.0, 0.0));
        vatt_05.mEit2 = findEdge  (m02, Vec3(-1.0,  1.0,-1.0),
                                       Vec3(-1.0, -1.0,-1.0));

        IntersectionFinder::Attributes vatt_06;
        vatt_06.mP    = rotQmat01 * Vec3(0.0, 0.0, 1.0) + trans01;
        vatt_06.mPred = IF_FACE_EDGE;
        vatt_06.mFit1 = findFace  (m01, Vec3( 0.0,  0.0, 1.0));
        vatt_06.mEit2 = findEdge  (m02, Vec3(-1.0, -1.0, 1.0),
                                       Vec3(-1.0, -1.0,-1.0));

        IntersectionFinder::Attributes vatt_07;
        vatt_07.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        vatt_07.mPred = IF_VERTEX_INTERIOR;
        vatt_07.mVit1 = findVertex(m01, Vec3( 1.0, 1.0, 1.0));

        IntersectionFinder::Attributes vatt_08;
        vatt_08.mP    = rotQmat01 * Vec3(0.0, 0.0, 0.0) + trans01;
        vatt_08.mPred = IF_INTERIOR_VERTEX;
        vatt_08.mVit2 = findVertex(m02, Vec3( -1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> vatts_01;
        vatts_01.push_back(vatt_01);
        vatts_01.push_back(vatt_02);
        vatts_01.push_back(vatt_03);
        vatts_01.push_back(vatt_04);
        vatts_01.push_back(vatt_05);
        vatts_01.push_back(vatt_06);
        vatts_01.push_back(vatt_07);
        vatts_01.push_back(vatt_08);

        vector<pair<Vec3, Vec3> > ePoints;
        vector<IntersectionFinder::Attributes> eatts_01;

        auto eps_01 = make_pair(rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 0.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_01;
        eatt_01.mPred = IF_EDGE_INTERIOR;
        eatt_01.mEit1 = findEdge(m01, Vec3( 1.0,  1.0,  1.0),
                                      Vec3( 1.0, -1.0,  1.0));

        auto eps_02 = make_pair(rotQmat01 * Vec3( 1.0, 0.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 0.0, 0.0) + trans01 );

        IntersectionFinder::Attributes eatt_02;
        eatt_02.mPred = IF_FACE_FACE;
        eatt_02.mFit1 = findFace(m01, Vec3( 1.0,  0.0,  0.0));
        eatt_02.mFit2 = findFace(m02, Vec3( 0.0, -1.0,  0.0));

        auto eps_03 = make_pair(rotQmat01 * Vec3( 1.0, 0.0, 0.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0, 0.0) + trans01 );

        IntersectionFinder::Attributes eatt_03;
        eatt_03.mPred = IF_FACE_FACE;
        eatt_03.mFit1 = findFace(m01, Vec3( 1.0,  0.0,  0.0));
        eatt_03.mFit2 = findFace(m02, Vec3( 0.0,  0.0, -1.0));

        auto eps_04 = make_pair(rotQmat01 * Vec3( 1.0, 1.0, 0.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_04;
        eatt_04.mPred = IF_EDGE_INTERIOR;
        eatt_04.mEit1 = findEdge(m01, Vec3( 1.0,  1.0,  1.0),
                                      Vec3( 1.0,  1.0, -1.0));

        auto eps_05 = make_pair(rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 0.0, 1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_05;
        eatt_05.mPred = IF_EDGE_INTERIOR;
        eatt_05.mEit1 = findEdge(m01, Vec3( 1.0,  1.0,  1.0),
                                      Vec3(-1.0,  1.0,  1.0));

        auto eps_06 = make_pair(rotQmat01 * Vec3( 1.0, 0.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 0.0, 0.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_06;
        eatt_06.mPred = IF_FACE_FACE;
        eatt_06.mFit1 = findFace(m01, Vec3( 0.0,  0.0,  1.0));
        eatt_06.mFit2 = findFace(m02, Vec3( 0.0, -1.0,  0.0));

        auto eps_07 = make_pair(rotQmat01 * Vec3( 1.0, 0.0, 0.0) + trans01,
                                rotQmat01 * Vec3( 0.0, 0.0, 0.0) + trans01 );

        IntersectionFinder::Attributes eatt_07;
        eatt_07.mPred = IF_INTERIOR_EDGE;
        eatt_07.mEit2 = findEdge(m02, Vec3( 1.0, -1.0, -1.0),
                                      Vec3(-1.0, -1.0, -1.0));

        auto eps_08 = make_pair(rotQmat01 * Vec3( 1.0, 1.0, 0.0) + trans01,
                                rotQmat01 * Vec3( 0.0, 1.0, 0.0) + trans01 );

        IntersectionFinder::Attributes eatt_08;

        eatt_08.mPred = IF_FACE_FACE;
        eatt_08.mFit1 = findFace(m01, Vec3( 0.0,  1.0,  0.0));
        eatt_08.mFit2 = findFace(m02, Vec3( 0.0,  0.0, -1.0));

        auto eps_09 = make_pair(rotQmat01 * Vec3( 0.0, 1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 0.0, 0.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_09;

        eatt_09.mPred = IF_FACE_FACE;
        eatt_09.mFit1 = findFace(m01, Vec3( 0.0,  0.0,  1.0));
        eatt_09.mFit2 = findFace(m02, Vec3(-1.0,  0.0,  0.0));

        auto eps_10 = make_pair(rotQmat01 * Vec3( 0.0, 0.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 0.0, 0.0, 0.0) + trans01 );

        IntersectionFinder::Attributes eatt_10;

        eatt_10.mPred = IF_INTERIOR_EDGE;
        eatt_10.mEit2 = findEdge(m02, Vec3(-1.0, -1.0,  1.0),
                                      Vec3(-1.0, -1.0, -1.0));

        auto eps_11 = make_pair(rotQmat01 * Vec3( 0.0, 0.0, 0.0) + trans01,
                                rotQmat01 * Vec3( 0.0, 1.0, 0.0) + trans01 );

        IntersectionFinder::Attributes eatt_11;

        eatt_11.mPred = IF_INTERIOR_EDGE;
        eatt_11.mEit2 = findEdge(m02, Vec3(-1.0,  1.0, -1.0),
                                      Vec3(-1.0, -1.0, -1.0));

        auto eps_12 = make_pair(rotQmat01 * Vec3( 0.0, 1.0, 0.0) + trans01,
                                rotQmat01 * Vec3( 0.0, 1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_12;

        eatt_12.mPred = IF_FACE_FACE;
        eatt_12.mFit1 = findFace(m01, Vec3( 0.0,  1.0,  0.0));
        eatt_12.mFit2 = findFace(m02, Vec3(-1.0,  0.0,  0.0));

        ePoints.push_back(eps_01);
        ePoints.push_back(eps_02);
        ePoints.push_back(eps_03);
        ePoints.push_back(eps_04);
        ePoints.push_back(eps_05);
        ePoints.push_back(eps_06);
        ePoints.push_back(eps_07);
        ePoints.push_back(eps_08);
        ePoints.push_back(eps_09);
        ePoints.push_back(eps_10);
        ePoints.push_back(eps_11);
        ePoints.push_back(eps_12);

        eatts_01.push_back(eatt_01);
        eatts_01.push_back(eatt_02);
        eatts_01.push_back(eatt_03);
        eatts_01.push_back(eatt_04);
        eatts_01.push_back(eatt_05);
        eatts_01.push_back(eatt_06);
        eatts_01.push_back(eatt_07);
        eatts_01.push_back(eatt_08);
        eatts_01.push_back(eatt_09);
        eatts_01.push_back(eatt_10);
        eatts_01.push_back(eatt_11);
        eatts_01.push_back(eatt_12);

        IntersectionFinder::Attributes fatt_01;
        fatt_01.mP    = rotQmat01 * Vec3( 1.0, 0.0, 0.0);
        fatt_01.mPred = IF_FACE_INTERIOR;
        fatt_01.mFit1 = findFace(m01, Vec3( 1.0, 0.0, 0.0));

        IntersectionFinder::Attributes fatt_02;
        fatt_02.mP    = rotQmat01 * Vec3( 0.0, 1.0, 0.0);
        fatt_02.mPred = IF_FACE_INTERIOR;
        fatt_02.mFit1 = findFace(m01, Vec3( 0.0, 1.0, 0.0));

        IntersectionFinder::Attributes fatt_03;
        fatt_03.mP    = rotQmat01 * Vec3( 0.0, 0.0, 1.0);
        fatt_03.mPred = IF_FACE_INTERIOR;
        fatt_03.mFit1 = findFace(m01, Vec3( 0.0, 0.0, 1.0));

        IntersectionFinder::Attributes fatt_04;
        fatt_04.mP    = rotQmat01 * Vec3(-1.0, 0.0, 0.0);
        fatt_04.mPred = IF_INTERIOR_FACE;
        fatt_04.mFit2 = findFace(m02, Vec3(-1.0, 0.0, 0.0));

        IntersectionFinder::Attributes fatt_05;
        fatt_05.mP    = rotQmat01 * Vec3( 0.0, -1.0, 0.0);
        fatt_05.mPred = IF_INTERIOR_FACE;
        fatt_05.mFit2 = findFace(m02, Vec3(0.0, -1.0, 0.0));

        IntersectionFinder::Attributes fatt_06;
        fatt_06.mP    = rotQmat01 * Vec3( 0.0, 0.0, -1.0);
        fatt_06.mPred = IF_INTERIOR_FACE;
        fatt_06.mFit2 = findFace(m02, Vec3(0.0, 0.0, -1.0));

        vector<IntersectionFinder::Attributes> fatts_01;
        fatts_01.push_back(fatt_01);
        fatts_01.push_back(fatt_02);
        fatts_01.push_back(fatt_03);
        fatts_01.push_back(fatt_04);
        fatts_01.push_back(fatt_05);
        fatts_01.push_back(fatt_06);

        auto res = checkResults3D(
                          finder01, vatts_01, ePoints, eatts_01, fatts_01);
        EXPECT_EQ(res, true);

    }
}


/**  @brief find()
 *          from test13
 */
TEST_F(IntersectionFinderTest, Test32) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(1.0, 0.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.find();

        IntersectionFinder::Attributes vatt_01;
        vatt_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        vatt_01.mPred = IF_VERTEX_EDGE;
        vatt_01.mVit1 = findVertex(m01, Vec3( 1.0,  1.0,  1.0));
        vatt_01.mEit2 = findEdge  (m02, Vec3( 1.0,  1.0,  1.0),
                                       Vec3(-1.0,  1.0,  1.0));

        IntersectionFinder::Attributes vatt_02;
        vatt_02.mP    = rotQmat01 * Vec3(1.0, 1.0, -1.0) + trans01;
        vatt_02.mPred = IF_VERTEX_EDGE;
        vatt_02.mVit1 = findVertex(m01, Vec3( 1.0,  1.0, -1.0));
        vatt_02.mEit2 = findEdge  (m02, Vec3( 1.0,  1.0, -1.0),
                                       Vec3(-1.0,  1.0, -1.0));

        IntersectionFinder::Attributes vatt_03;
        vatt_03.mP    = rotQmat01 * Vec3(1.0, -1.0, -1.0) + trans01;
        vatt_03.mPred = IF_VERTEX_EDGE;
        vatt_03.mVit1 = findVertex(m01, Vec3( 1.0,  -1.0, -1.0));
        vatt_03.mEit2 = findEdge  (m02, Vec3( 1.0,  -1.0, -1.0),
                                       Vec3(-1.0,  -1.0, -1.0));

        IntersectionFinder::Attributes vatt_04;
        vatt_04.mP    = rotQmat01 * Vec3(1.0, -1.0,  1.0) + trans01;
        vatt_04.mPred = IF_VERTEX_EDGE;
        vatt_04.mVit1 = findVertex(m01, Vec3( 1.0,  -1.0,  1.0));
        vatt_04.mEit2 = findEdge  (m02, Vec3( 1.0,  -1.0,  1.0),
                                       Vec3(-1.0,  -1.0,  1.0));

        IntersectionFinder::Attributes vatt_05;
        vatt_05.mP    = rotQmat01 * Vec3(0.0, 1.0, 1.0) + trans01;
        vatt_05.mPred = IF_EDGE_VERTEX;
        vatt_05.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0,  1.0),
                                       Vec3(-1.0,  1.0,  1.0));
        vatt_05.mVit2 = findVertex(m02, Vec3(-1.0,  1.0,  1.0));

        IntersectionFinder::Attributes vatt_06;
        vatt_06.mP    = rotQmat01 * Vec3(0.0, 1.0, -1.0) + trans01;
        vatt_06.mPred = IF_EDGE_VERTEX;
        vatt_06.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0, -1.0),
                                       Vec3(-1.0,  1.0, -1.0));
        vatt_06.mVit2 = findVertex(m02, Vec3(-1.0,  1.0, -1.0));

        IntersectionFinder::Attributes vatt_07;
        vatt_07.mP    = rotQmat01 * Vec3(0.0, -1.0, -1.0) + trans01;
        vatt_07.mPred = IF_EDGE_VERTEX;
        vatt_07.mEit1 = findEdge  (m01, Vec3( 1.0, -1.0, -1.0),
                                       Vec3(-1.0, -1.0, -1.0));
        vatt_07.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        IntersectionFinder::Attributes vatt_08;
        vatt_08.mP    = rotQmat01 * Vec3(0.0, -1.0, 1.0) + trans01;
        vatt_08.mPred = IF_EDGE_VERTEX;
        vatt_08.mEit1 = findEdge  (m01, Vec3( 1.0, -1.0, 1.0),
                                       Vec3(-1.0, -1.0, 1.0));
        vatt_08.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, 1.0));

        vector<IntersectionFinder::Attributes> vatts_01;
        vatts_01.push_back(vatt_01);
        vatts_01.push_back(vatt_02);
        vatts_01.push_back(vatt_03);
        vatts_01.push_back(vatt_04);
        vatts_01.push_back(vatt_05);
        vatts_01.push_back(vatt_06);
        vatts_01.push_back(vatt_07);
        vatts_01.push_back(vatt_08);

        vector<pair<Vec3, Vec3> > ePoints;
        vector<IntersectionFinder::Attributes> eatts_01;

        auto eps_01 = make_pair(rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0,-1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_01;
        eatt_01.mPred = IF_EDGE_FACE;
        eatt_01.mEit1 = findEdge(m01, Vec3( 1.0,  1.0,  1.0),
                                      Vec3( 1.0, -1.0,  1.0));
        eatt_01.mFit2 = findFace(m02, Vec3( 0.0,  0.0,  1.0));

        auto eps_02 = make_pair(rotQmat01 * Vec3( 1.0,-1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0,-1.0,-1.0) + trans01 );

        IntersectionFinder::Attributes eatt_02;
        eatt_02.mPred = IF_EDGE_FACE;
        eatt_02.mEit1 = findEdge(m01, Vec3( 1.0, -1.0,  1.0),
                                      Vec3( 1.0, -1.0, -1.0));
        eatt_02.mFit2 = findFace(m02, Vec3( 0.0, -1.0,  0.0));

        auto eps_03 = make_pair(rotQmat01 * Vec3( 1.0,-1.0,-1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0,-1.0) + trans01 );

        IntersectionFinder::Attributes eatt_03;
        eatt_03.mPred = IF_EDGE_FACE;
        eatt_03.mEit1 = findEdge(m01, Vec3( 1.0, -1.0, -1.0),
                                      Vec3( 1.0,  1.0, -1.0));
        eatt_03.mFit2 = findFace(m02, Vec3( 0.0,  0.0, -1.0));

        auto eps_04 = make_pair(rotQmat01 * Vec3( 1.0, 1.0,-1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_04;
        eatt_04.mPred = IF_EDGE_FACE;
        eatt_04.mEit1 = findEdge(m01, Vec3( 1.0,  1.0, -1.0),
                                      Vec3( 1.0,  1.0,  1.0));
        eatt_04.mFit2 = findFace(m02, Vec3( 0.0,  1.0,  0.0));

        auto eps_05 = make_pair(rotQmat01 * Vec3( 0.0, 1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_05;
        eatt_05.mPred = IF_EDGE_EDGE;
        eatt_05.mEit1 = findEdge(m01, Vec3( 1.0,  1.0,  1.0),
                                      Vec3(-1.0,  1.0,  1.0));
        eatt_05.mEit2 = findEdge(m02, Vec3( 1.0,  1.0,  1.0),
                                      Vec3(-1.0,  1.0,  1.0));

        auto eps_06 = make_pair(rotQmat01 * Vec3( 0.0,-1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 1.0,-1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_06;
        eatt_06.mPred = IF_EDGE_EDGE;
        eatt_06.mEit1 = findEdge(m01, Vec3( 1.0, -1.0,  1.0),
                                      Vec3(-1.0, -1.0,  1.0));
        eatt_06.mEit2 = findEdge(m02, Vec3( 1.0, -1.0,  1.0),
                                      Vec3(-1.0, -1.0,  1.0));

        auto eps_07 = make_pair(rotQmat01 * Vec3( 0.0,-1.0,-1.0) + trans01,
                                rotQmat01 * Vec3( 1.0,-1.0,-1.0) + trans01 );

        IntersectionFinder::Attributes eatt_07;
        eatt_07.mPred = IF_EDGE_EDGE;
        eatt_07.mEit1 = findEdge(m01, Vec3( 1.0, -1.0, -1.0),
                                      Vec3(-1.0, -1.0, -1.0));
        eatt_07.mEit2 = findEdge(m02, Vec3( 1.0, -1.0, -1.0),
                                      Vec3(-1.0, -1.0, -1.0));

        auto eps_08 = make_pair(rotQmat01 * Vec3( 0.0, 1.0,-1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0,-1.0) + trans01 );

        IntersectionFinder::Attributes eatt_08;
        eatt_08.mPred = IF_EDGE_EDGE;
        eatt_08.mEit1 = findEdge(m01, Vec3( 1.0,  1.0, -1.0),
                                      Vec3(-1.0,  1.0, -1.0));
        eatt_08.mEit2 = findEdge(m02, Vec3( 1.0,  1.0, -1.0),
                                      Vec3(-1.0,  1.0, -1.0));

        auto eps_09 = make_pair(rotQmat01 * Vec3( 0.0, 1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 0.0,-1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_09;
        eatt_09.mPred = IF_FACE_EDGE;
        eatt_09.mFit1 = findFace(m01, Vec3( 0.0,  0.0,  1.0));
        eatt_09.mEit2 = findEdge(m02, Vec3(-1.0,  1.0,  1.0),
                                      Vec3(-1.0, -1.0,  1.0));

        auto eps_10 = make_pair(rotQmat01 * Vec3( 0.0,-1.0, 1.0) + trans01,
                                rotQmat01 * Vec3( 0.0,-1.0,-1.0) + trans01 );

        IntersectionFinder::Attributes eatt_10;
        eatt_10.mPred = IF_FACE_EDGE;
        eatt_10.mFit1 = findFace(m01, Vec3( 0.0, -1.0,  0.0));
        eatt_10.mEit2 = findEdge(m02, Vec3(-1.0, -1.0,  1.0),
                                      Vec3(-1.0, -1.0, -1.0));

        auto eps_11 = make_pair(rotQmat01 * Vec3( 0.0,-1.0,-1.0) + trans01,
                                rotQmat01 * Vec3( 0.0, 1.0,-1.0) + trans01 );

        IntersectionFinder::Attributes eatt_11;
        eatt_11.mPred = IF_FACE_EDGE;
        eatt_11.mFit1 = findFace(m01, Vec3( 0.0,  0.0, -1.0));
        eatt_11.mEit2 = findEdge(m02, Vec3(-1.0, -1.0, -1.0),
                                      Vec3(-1.0,  1.0, -1.0));

        auto eps_12 = make_pair(rotQmat01 * Vec3( 0.0, 1.0,-1.0) + trans01,
                                rotQmat01 * Vec3( 0.0, 1.0, 1.0) + trans01 );

        IntersectionFinder::Attributes eatt_12;
        eatt_12.mPred = IF_FACE_EDGE;
        eatt_12.mFit1 = findFace(m01, Vec3( 0.0,  1.0,  0.0));
        eatt_12.mEit2 = findEdge(m02, Vec3(-1.0,  1.0, -1.0),
                                      Vec3(-1.0,  1.0,  1.0));


        ePoints.push_back(eps_01);
        ePoints.push_back(eps_02);
        ePoints.push_back(eps_03);
        ePoints.push_back(eps_04);
        ePoints.push_back(eps_05);
        ePoints.push_back(eps_06);
        ePoints.push_back(eps_07);
        ePoints.push_back(eps_08);
        ePoints.push_back(eps_09);
        ePoints.push_back(eps_10);
        ePoints.push_back(eps_11);
        ePoints.push_back(eps_12);

        eatts_01.push_back(eatt_01);
        eatts_01.push_back(eatt_02);
        eatts_01.push_back(eatt_03);
        eatts_01.push_back(eatt_04);
        eatts_01.push_back(eatt_05);
        eatts_01.push_back(eatt_06);
        eatts_01.push_back(eatt_07);
        eatts_01.push_back(eatt_08);
        eatts_01.push_back(eatt_09);
        eatts_01.push_back(eatt_10);
        eatts_01.push_back(eatt_11);
        eatts_01.push_back(eatt_12);

        IntersectionFinder::Attributes fatt_01;
        fatt_01.mP    = rotQmat01 * Vec3( 1.0, 0.0, 0.0);
        fatt_01.mPred = IF_FACE_INTERIOR;
        fatt_01.mFit1 = findFace(m01, Vec3( 1.0, 0.0, 0.0));

        IntersectionFinder::Attributes fatt_02;
        fatt_02.mP    = rotQmat01 * Vec3( 0.0, 1.0, 0.0);
        fatt_02.mPred = IF_FACE_FACE;
        fatt_02.mFit1 = findFace(m01, Vec3( 0.0, 1.0, 0.0));
        fatt_02.mFit2 = findFace(m02, Vec3( 0.0, 1.0, 0.0));

        IntersectionFinder::Attributes fatt_03;
        fatt_03.mP    = rotQmat01 * Vec3( 0.0, 0.0, 1.0);
        fatt_03.mPred = IF_FACE_FACE;
        fatt_03.mFit1 = findFace(m01, Vec3( 0.0, 0.0, 1.0));
        fatt_03.mFit2 = findFace(m02, Vec3( 0.0, 0.0, 1.0));

        IntersectionFinder::Attributes fatt_04;
        fatt_04.mP    = rotQmat01 * Vec3(-1.0, 0.0, 0.0);
        fatt_04.mPred = IF_INTERIOR_FACE;
        fatt_04.mFit2 = findFace(m02, Vec3(-1.0, 0.0, 0.0));

        IntersectionFinder::Attributes fatt_05;
        fatt_05.mP    = rotQmat01 * Vec3( 0.0, -1.0, 0.0);
        fatt_05.mPred = IF_FACE_FACE;
        fatt_05.mFit1 = findFace(m01, Vec3(0.0, -1.0, 0.0));
        fatt_05.mFit2 = findFace(m02, Vec3(0.0, -1.0, 0.0));

        IntersectionFinder::Attributes fatt_06;
        fatt_06.mP    = rotQmat01 * Vec3( 0.0, 0.0, -1.0);
        fatt_06.mPred = IF_FACE_FACE;
        fatt_06.mFit1 = findFace(m01, Vec3(0.0, 0.0, -1.0));
        fatt_06.mFit2 = findFace(m02, Vec3(0.0, 0.0, -1.0));

        vector<IntersectionFinder::Attributes> fatts_01;
        fatts_01.push_back(fatt_01);
        fatts_01.push_back(fatt_02);
        fatts_01.push_back(fatt_03);
        fatts_01.push_back(fatt_04);
        fatts_01.push_back(fatt_05);
        fatts_01.push_back(fatt_06);

        auto res = checkResults3D(
                          finder01, vatts_01, ePoints, eatts_01, fatts_01);
        EXPECT_EQ(res, true);
    }
}


/**  @brief find()
 *          from test14
 */
TEST_F(IntersectionFinderTest, Test33) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 2.0, 0.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.find(sepAxis);

        IntersectionFinder::Attributes vatt_01;
        vatt_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        vatt_01.mPred = IF_VERTEX_VERTEX;
        vatt_01.mVit1 = findVertex(m01, Vec3( 1.0,  1.0,  1.0));
        vatt_01.mVit2 = findVertex(m02, Vec3(-1.0, -1.0,  1.0));

        IntersectionFinder::Attributes vatt_02;
        vatt_02.mP    = rotQmat01 * Vec3(1.0, 1.0, -1.0) + trans01;
        vatt_02.mPred = IF_VERTEX_VERTEX;
        vatt_02.mVit1 = findVertex(m01, Vec3( 1.0,  1.0, -1.0));
        vatt_02.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> vatts_01;
        vatts_01.push_back(vatt_01);
        vatts_01.push_back(vatt_02);

        auto eps_01 = make_pair(rotQmat01 * Vec3( 1.0, 1.0,  1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0, -1.0) + trans01 );

        IntersectionFinder::Attributes eatt_01;
        eatt_01.mPred = IF_EDGE_EDGE;
        eatt_01.mEit1 = findEdge(m01, Vec3(  1.0,  1.0,  1.0),
                                      Vec3(  1.0,  1.0, -1.0) );
        eatt_01.mEit2 = findEdge(m02, Vec3( -1.0, -1.0,  1.0),
                                      Vec3( -1.0, -1.0, -1.0) );
        

        EXPECT_EQ(checkResults1D(finder01, vatts_01, eps_01, eatt_01), true);

    }
}


/**  @brief find()
 *          from test15
 */
TEST_F(IntersectionFinderTest, Test34) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 2.0, 1.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.find(sepAxis);

        IntersectionFinder::Attributes vatt_01;
        vatt_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        vatt_01.mPred = IF_VERTEX_EDGE;
        vatt_01.mVit1 = findVertex(m01, Vec3( 1.0,  1.0, 1.0));
        vatt_01.mEit2 = findEdge  (m02, Vec3(-1.0, -1.0,  1.0),
                                       Vec3(-1.0, -1.0, -1.0));

        IntersectionFinder::Attributes vatt_02;
        vatt_02.mP    = rotQmat01 * Vec3(1.0, 1.0,  0.0) + trans01;
        vatt_02.mPred = IF_EDGE_VERTEX;
        vatt_02.mEit1 = findEdge  (m01, Vec3( 1.0,  1.0,  1.0),
                                       Vec3( 1.0,  1.0, -1.0));
        vatt_02.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        vector<IntersectionFinder::Attributes> vatts_01;
        vatts_01.push_back(vatt_01);
        vatts_01.push_back(vatt_02);

        auto eps_01 = make_pair(rotQmat01 * Vec3( 1.0, 1.0,  1.0) + trans01,
                                rotQmat01 * Vec3( 1.0, 1.0,  0.0) + trans01 );

        IntersectionFinder::Attributes eatt_01;
        eatt_01.mPred = IF_EDGE_EDGE;
        eatt_01.mEit1 = findEdge(m01, Vec3(  1.0,  1.0,  1.0),
                                      Vec3(  1.0,  1.0, -1.0) );
        eatt_01.mEit2 = findEdge(m02, Vec3( -1.0, -1.0,  1.0),
                                      Vec3( -1.0, -1.0, -1.0) );
        

        EXPECT_EQ(checkResults1D(finder01, vatts_01, eps_01, eatt_01), true);



    }
}


/**  @brief find()
 *          from test16
 */
TEST_F(IntersectionFinderTest, Test35) {

    for (long i = 0; i < 100; i++) {

        Manifold m01,    m02;

        Vec3 p01( -1.0,  -1.0,   1.0);// frontLowerLeft
        Vec3 p02( -1.0,   1.0,   1.0);// frontUpperLeft
        Vec3 p03(  1.0,   1.0,   1.0);// frontUpperRight
        Vec3 p04(  1.0,  -1.0,   1.0);// frontLowerRight
        Vec3 p05( -1.0,  -1.0,  -1.0);// backLowerLeft
        Vec3 p06( -1.0,   1.0,  -1.0);// backUpperLeft
        Vec3 p07(  1.0,   1.0,  -1.0);// backUpperRight
        Vec3 p08(  1.0,  -1.0,  -1.0);// backLowerRight

        m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
        m02.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);

        Vec3       sepAxis(1.0, 0.0, 0.0);
        Vec3       dispCoM(2.0, 2.0, 2.0);

        Quaternion Q01(1.0, 0.0, 0.0, 0.0);
        Quaternion Q02(1.0, 0.0, 0.0, 0.0);
        Mat3x3     Qmat01 = Q01.rotationMatrix();
        Mat3x3     Qmat02 = Q02.rotationMatrix();
        Vec3       CoM01 (0.0, 0.0, 0.0);
        Vec3       CoM02 (0.0, 0.0, 0.0);

        CoM02 += dispCoM;

        Quaternion rotQ01    = randomRotQ();
        Mat3x3     rotQmat01 = rotQ01.rotationMatrix();
        Vec3       trans01   = randVec3D100();

        Qmat01  =  rotQmat01 * Qmat01;
        Qmat02  =  rotQmat01 * Qmat02;
        CoM01   =  rotQmat01 * CoM01;
        CoM02   =  rotQmat01 * CoM02;
        CoM01   += trans01;
        CoM02   += trans01;
        sepAxis =  rotQmat01 * sepAxis;

        double epsilonZero    = EPSILON_LINEAR;
        double epsilonZeroPCA = EPSILON_LINEAR*1000.0;
        double epsilonAngle   = EPSILON_LINEAR;

        IntersectionFinder finder01( m01, Qmat01, CoM01, m02, Qmat02, CoM02,
                                    epsilonZero, epsilonZeroPCA,epsilonAngle );

        finder01.find(sepAxis);

        IntersectionFinder::Attributes vatt_01;
        vatt_01.mP    = rotQmat01 * Vec3(1.0, 1.0, 1.0) + trans01;
        vatt_01.mPred = IF_VERTEX_VERTEX;
        vatt_01.mVit1 = findVertex(m01, Vec3( 1.0,  1.0,  1.0));
        vatt_01.mVit2 = findVertex(m02, Vec3(-1.0, -1.0, -1.0));

        EXPECT_EQ(checkResults0D(finder01, vatt_01), true);
    }
}

}// Makena

#include "contact_updater_edge_edge.hpp"

/**
 * @file contact_updater_edge_edge.cpp
 *
 * @brief
 */
namespace Makena {


ContactUpdater_EDGE_EDGE::ContactUpdater_EDGE_EDGE(
    ConvexRigidBody&         body1,
    ConvexRigidBody&         body2,
    ContactPairInfo&         info,
    const double&            epsilonZero,
    const double&            epsilonAngle,
    std::ostream&            logStream
)
    :Loggable(logStream),
     mBody1(body1),
     mBody2(body2),
     mInfo(info),
     mEpsilonZero(epsilonZero),
     mEpsilonAngle(epsilonAngle){;}


ContactUpdater_EDGE_EDGE::~ContactUpdater_EDGE_EDGE(){;}


bool ContactUpdater_EDGE_EDGE:: update()
{
    if (mInfo.mContactNormalType == ContactPairInfo::NT_EDGE_CROSS_EDGE) {
        return update_CROSS();
    }
    else {
        return update_PARALLEL();
    }
}


bool ContactUpdater_EDGE_EDGE:: update_CROSS()
{
    rotateWorld(mInfo.mContactNormal1To2);

    bool intersectionFound = findIntersection2D();

    if (intersectionFound) {
        dispatchAndUpdate();
        return false;
    }
    else {
        return true;
    }
}


bool ContactUpdater_EDGE_EDGE:: update_PARALLEL()
{
    auto heit1 = ((*mInfo.mEit1))->he1();
    auto vit11 = (*heit1)->src();
    auto vit12 = (*heit1)->dst();
    Vec3 p11   = (*vit11)->pGCS(mBody1.Qmat(), mBody1.CoM());
    Vec3 p12   = (*vit12)->pGCS(mBody1.Qmat(), mBody1.CoM());
    Vec3 v1    = p12 - p11;
    v1.normalize();
    auto heit2  = ((*mInfo.mEit2))->he1();
    auto vit21 = (*heit2)->src();
    auto vit22 = (*heit2)->dst();
    Vec3 p21    = (*vit21)->pGCS(mBody2.Qmat(), mBody2.CoM());
    Vec3 p22    = (*vit22)->pGCS(mBody2.Qmat(), mBody2.CoM());
    Vec3 v2     = p22 - p21;
    v2.normalize();
    if (v1.dot(v2) < 0.0) {
        v2.scale(-1.0);
    }
    Vec3 zDir  = v1 + v2; 
    zDir.normalize();

    // Contact Normal => X
    // Edge Direction => Z
    rotateWorld(zDir, mInfo.mContactNormal1To2);

    p11 = (*vit11)->pGCS(mRotMat1, mCom1);
    p12 = (*vit12)->pGCS(mRotMat1, mCom1);
    p21 = (*vit21)->pGCS(mRotMat2, mCom2);
    p22 = (*vit22)->pGCS(mRotMat2, mCom2);

    if (p11.z() > p12.z()) {
        std::swap(vit11, vit12);
        std::swap(p11,   p12);
    }

    if (p21.z() > p22.z()) {
        std::swap(vit21, vit22);
        std::swap(p21,   p22);
    }

    // p11 < p12 < p21 << p22 || p21 < p22 < p11 << p12
    if ( (p12.z() + mEpsilonZero) < p21.z() ||
         (p22.z() + mEpsilonZero) < p11.z()    ) {
        return true;// No intersection
    }

    Vec2 v1_2d(p12.z() - p11.z(), p12.x() - p11.x());
    Vec2 v2_2d(p22.z() - p21.z(), p22.x() - p21.x());
    v1_2d.normalize();
    v2_2d.normalize();
    Vec2 v2_2d_perp = v2_2d.perp();
    auto dot = v1_2d.dot(v2_2d_perp);
    if (fabs(dot) <= mEpsilonAngle) {
        ; // No need to demote.
    }
    else if (dot > 0.0) {
        // p12, p22 are penetrating
        if (fabs(p12.z()-p22.z())<=mEpsilonAngle) {
            mInfo.mType1 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit1  = vit12;
            mInfo.mType2 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit2  = vit22;
        }
        else if (p12.z() < p22.z()) {
            mInfo.mType1 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit1  = vit12;
        }
        else {
            mInfo.mType2 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit2  = vit22;
        }
    }
    else {
        // p11, p21 are penetrating
        if (fabs(p11.z()-p21.z())<=mEpsilonAngle) {
            mInfo.mType1 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit1  = vit11;
            mInfo.mType2 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit2  = vit21;
        }
        else if (p11.z() > p21.z()) {
            mInfo.mType1 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit1  = vit11;
        }
        else {
            mInfo.mType2 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit2  = vit21;
        }
    }
    return false;
}



void ContactUpdater_EDGE_EDGE::rotateWorld(const Vec3& zDir, const Vec3& xDir)
{
    const Vec3 yDir = zDir.cross(xDir);
    Mat3x3 rMat(xDir, yDir, zDir);
    rMat.transposeInPlace();
    mRotMat1 = rMat * mBody1.Qmat();
    mCom1    = rMat * mBody1.CoM();
    mRotMat2 = rMat * mBody2.Qmat();
    mCom2    = rMat * mBody2.CoM();
}





void ContactUpdater_EDGE_EDGE::rotateWorld(const Vec3& zDir)
{
    Mat3x3 rMat = rotMatAlignZDirToZAxis(zDir);

    mRotMat1 = rMat * mBody1.Qmat();
    mCom1    = rMat * mBody1.CoM();
    mRotMat2 = rMat * mBody2.Qmat();
    mCom2    = rMat * mBody2.CoM();

}


Mat3x3 ContactUpdater_EDGE_EDGE::rotMatAlignZDirToZAxis(const Vec3& zDir)
{
    // Find the new coordinate axes with zDir for new z-axis.
    Vec3 r3 = zDir;
    double ax = fabs(r3.x());
    double ay = fabs(r3.y());
    double az = fabs(r3.z());
    Vec3 r1;
    if (ax <= ay && ax <= az) {
        Vec3 xUnit(1.0, 0.0, 0.0);
        r1 = r3.cross(xUnit);
    }
    else if (ay <= ax && ay <= az) {
        Vec3 yUnit(0.0, 1.0, 0.0);
        r1 = r3.cross(yUnit);
    }
    else {
        Vec3 zUnit(0.0, 0.0, 1.0);
        r1 = r3.cross(zUnit);
    }

    r1.normalize();
    // r3 x r1 = r2 => r1 x r2 = r3
    const Vec3 r2 = r3.cross(r1);


    // Find the rotation matrix to rotate a point in GCS into the new CS.
    Mat3x3 rMat(r1, r2, r3);
    rMat.transposeInPlace();

    return rMat;
}

bool ContactUpdater_EDGE_EDGE::findIntersection2D()
{
    vector< IntSec2D::InputElem> inputElem1;
    vector< IntSec2D::InputElem> inputElem2;
    mIntsec.clear();

    auto heit1 = ((*mInfo.mEit1))->he1();
    mVit11     = (*heit1)->src();
    mVit12     = (*heit1)->dst();
    Vec3 p11   = ((*mVit11)->pGCS(mRotMat1, mCom1));
    Vec3 p12   = ((*mVit12)->pGCS(mRotMat1, mCom1));
    Vec2 p11_2d(p11.x(), p11.y());
    Vec2 p12_2d(p12.x(), p12.y());
    IntSec2D::InputElem ie11(p11_2d, 1);
    inputElem1.push_back(ie11);
    IntSec2D::InputElem ie12(p12_2d, 2);
    inputElem1.push_back(ie12);

    auto heit2 = ((*mInfo.mEit2))->he1();
    mVit21     = (*heit2)->src();
    mVit22     = (*heit2)->dst();
    Vec3 p21   = ((*mVit21)->pGCS(mRotMat2, mCom2));
    Vec3 p22   = ((*mVit22)->pGCS(mRotMat2, mCom2));
    Vec2 p21_2d(p21.x(), p21.y());
    Vec2 p22_2d(p22.x(), p22.y());

    IntSec2D::InputElem ie21(p21_2d, 1);
    inputElem2.push_back(ie21);
    IntSec2D::InputElem ie22(p22_2d, 2);
    inputElem2.push_back(ie22);

    IntSec2D IntsecFinder2D(
                      mEpsilonZero, mEpsilonZero, mEpsilonAngle, mLogStream);
    return IntsecFinder2D.findIntersection(inputElem1, inputElem2, mIntsec);
}

void ContactUpdater_EDGE_EDGE::dispatchAndUpdate()
{

    if (mIntsec.size()==1) {
        auto& oe = mIntsec[0];
        switch(oe.mType) {
          case IntSec2D::IT_EDGE_VERTEX:
            processIntsec_EDGE_VERTEX();
            break;
          case IntSec2D::IT_VERTEX_EDGE:
            processIntsec_VERTEX_EDGE();
            break;
          case IntSec2D::IT_VERTEX_VERTEX:
            processIntsec_VERTEX_VERTEX();
            break;
          default:
            break;
        }
    }
    // (mIntsec.size()==2)
    // This is a parallel EDGE-EDGE situation 
    // and no need to update mInfo.
}

void ContactUpdater_EDGE_EDGE::processIntsec_EDGE_VERTEX()
{
    auto& oe  = mIntsec[0];
    mInfo.mType2 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit2  = (oe.mIndexB==1)?mVit21:mVit22;
}


void ContactUpdater_EDGE_EDGE::processIntsec_VERTEX_EDGE()
{
    auto& oe  = mIntsec[0];
    mInfo.mType1 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit1  = (oe.mIndexA==1)?mVit11:mVit12;
}


void ContactUpdater_EDGE_EDGE::processIntsec_VERTEX_VERTEX()
{
    auto& oe  = mIntsec[0];
    mInfo.mType1 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit1  = (oe.mIndexA==1)?mVit11:mVit12;
    mInfo.mType2 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit2  = (oe.mIndexB==1)?mVit21:mVit22;
}



}// namespace Makena


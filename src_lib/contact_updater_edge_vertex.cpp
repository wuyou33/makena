#include "contact_updater_edge_vertex.hpp"

/**
 * @file contact_updater_edge_vertex.cpp
 *
 * @brief
 */
namespace Makena {


ContactUpdater_EDGE_VERTEX::ContactUpdater_EDGE_VERTEX(
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


ContactUpdater_EDGE_VERTEX::~ContactUpdater_EDGE_VERTEX(){;}


bool ContactUpdater_EDGE_VERTEX::update()
{
    if (mInfo.mType1 == ContactPairInfo::FT_EDGE) {
        auto eit   = mInfo.mEit1;
        auto heit  = (*eit)->he1();
        auto vit11 = (*heit)->src();
        auto vit12 = (*heit)->dst();
        auto p11   = (*vit11)->pGCS(mBody1.Qmat(), mBody1.CoM());
        auto p12   = (*vit12)->pGCS(mBody1.Qmat(), mBody1.CoM());

        auto vit2  = mInfo.mVit2;
        auto p2    = (*vit2)->pGCS(mBody2.Qmat(), mBody2.CoM());
        auto v11_2 = p2 - p11;
        if (v11_2.squaredNorm2()<=mEpsilonZero) {
            mInfo.mType1 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit1  = vit11;
            return false;
        }

        auto v12_2 = p2 - p12;
        if (v12_2.squaredNorm2()<=mEpsilonZero) {
            mInfo.mType1 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit1  = vit12;
            return false;
        }

        auto v11_12 = p12 - p11;
        double sqDist11_12 = v11_12.squaredNorm2();
        double sqDist11_2  = v11_2.squaredNorm2();
        v11_12.normalize();
        v11_2.normalize();
        auto cr = v11_12.cross(v11_2);
        if (v11_12.dot(v11_2) >= 0.0           &&
            cr.squaredNorm2() <= mEpsilonAngle &&
            sqDist11_2        <= sqDist11_12      ) {
            return false;
        }
        else {
            return true;
        }
    }
    else {

        auto eit   = mInfo.mEit2;
        auto heit  = (*eit)->he1();
        auto vit21 = (*heit)->src();
        auto vit22 = (*heit)->dst();
        auto p21   = (*vit21)->pGCS(mBody2.Qmat(), mBody2.CoM());
        auto p22   = (*vit22)->pGCS(mBody2.Qmat(), mBody2.CoM());

        auto vit1  = mInfo.mVit1;
        auto p1    = (*vit1)->pGCS(mBody1.Qmat(), mBody1.CoM());

        auto v21_1 = p1 - p21;
        if (v21_1.squaredNorm2()<=mEpsilonZero) {
            mInfo.mType2 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit2  = vit21;
            return false;
        }

        auto v22_1 = p1 - p22;
        if (v22_1.squaredNorm2()<=mEpsilonZero) {
            mInfo.mType2 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit2  = vit22;
            return false;
        }

        auto v21_22 = p22 - p21;
        double sqDist21_22 = v21_22.squaredNorm2();
        double sqDist21_1  = v21_1.squaredNorm2();

        v21_22.normalize();
        v21_1.normalize();
        auto cr = v21_22.cross(v21_1);
        if (v21_22.dot(v21_1) >= 0.0           &&
            cr.squaredNorm2() <= mEpsilonAngle &&
            sqDist21_1        <= sqDist21_22      ) {
            return false;
        }
        else {
            return true;
        }
    }
}




}// namespace Makena


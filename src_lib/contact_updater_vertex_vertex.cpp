#include "contact_updater_vertex_vertex.hpp"

/**
 * @file contact_updater_vertex_vertex.cpp
 *
 * @brief
 */
namespace Makena {


ContactUpdater_VERTEX_VERTEX::ContactUpdater_VERTEX_VERTEX(
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


ContactUpdater_VERTEX_VERTEX::~ContactUpdater_VERTEX_VERTEX(){;}


bool ContactUpdater_VERTEX_VERTEX::update()
{
    auto vit1  = mInfo.mVit1;
    auto p1    = (*vit1)->pGCS(mBody1.Qmat(), mBody1.CoM());

    auto vit2  = mInfo.mVit2;
    auto p2    = (*vit2)->pGCS(mBody2.Qmat(), mBody2.CoM());

    auto v12 = p2 - p1;
    if (v12.squaredNorm2() <= mEpsilonZero) {
        return false;
    }
    else {
        return true;
    }
}



}// namespace Makena


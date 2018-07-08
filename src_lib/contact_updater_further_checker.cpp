#include "convex_hull_2d.hpp"
#include "contact_updater_further_checker.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file contact_updater_further_checker.cpp
 *
 * @brief
 */
namespace Makena {

using namespace std;

/** @
 *
 *
 *
 */
ContactUpdaterFurtherChecker::ContactUpdaterFurtherChecker(
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
     mEpsilonAngle(epsilonAngle)
     {;}


ContactUpdaterFurtherChecker::~ContactUpdaterFurtherChecker(){;}


bool ContactUpdaterFurtherChecker::checkFeatures()
{
    if ( mInfo.mType1 == ContactPairInfo::FT_FACE ) {
        if ( mInfo.mType2 == ContactPairInfo::FT_FACE ) {
            return false;
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE ) {
            return checkFeatures_FACE_EDGE();
        }
        else {
            return checkFeatures_FACE_VERTEX();
        }

    }
    if ( mInfo.mType1 == ContactPairInfo::FT_EDGE ) {

        if ( mInfo.mType2 == ContactPairInfo::FT_FACE ) {
            return checkFeatures_EDGE_FACE();
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE ) {
            return checkFeatures_EDGE_EDGE();
        }
        else {
            return checkFeatures_EDGE_VERTEX();
        }

    }
    else {
        if ( mInfo.mType2 == ContactPairInfo::FT_FACE ) {
            return checkFeatures_VERTEX_FACE();
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE ) {
            return checkFeatures_VERTEX_EDGE();
        }
        else {
            return checkFeatures_VERTEX_VERTEX();
        }

    }
}


bool ContactUpdaterFurtherChecker::checkFeatures_FACE_EDGE()
{
    auto body2_he1  = (*mInfo.mEit2)->he1();
    auto body2_he2  = (*mInfo.mEit2)->he2();
    auto body2_vit1 = (*body2_he1)->src();
    auto body2_vit2 = (*body2_he1)->dst();
    const auto body2_p1   = (*body2_vit1)->pGCS(mBody2.Qmat(), mBody2.CoM());
    const auto body2_p2   = (*body2_vit2)->pGCS(mBody2.Qmat(), mBody2.CoM());
    auto body2_v12  = body2_p2 - body2_p1;
    body2_v12.normalize();
    auto body2_fit1 = (*body2_he1)->face();
    auto body2_fit2 = (*body2_he2)->face();

    const auto body2_n1 = (*body2_fit1)->nGCS(mBody2.Qmat());
    const auto body2_n2 = (*body2_fit2)->nGCS(mBody2.Qmat());

    const auto body2_d1 = body2_n1.cross(body2_v12); // from v12 inward face 1
    const auto body2_d2 = body2_v12.cross(body2_n2); // from v12 inward face 2

    const auto body1_n  = (*mInfo.mFit1)->nGCS(mBody1.Qmat());

    const auto body2_d1_dot_body1_n = body2_d1.dot(body1_n);
    const auto body2_d2_dot_body1_n = body2_d2.dot(body1_n);

    const auto body2_n1_dot_body1_n = body2_n1.dot(body1_n);
    const auto body2_n2_dot_body1_n = body2_n2.dot(body1_n);

    if (body2_d1_dot_body1_n >= mEpsilonAngle && 
        body2_d2_dot_body1_n >= mEpsilonAngle   ) {

        return false; // Good.

    }
    else if (body2_d1_dot_body1_n <= 0.0 && 
             body2_d2_dot_body1_n <= 0.0   ) {
        return true; // Il-conditioned. Need to run GJK.

    }
    else if (body2_n1_dot_body1_n >= -1.0 * mEpsilonAngle && 
             body2_n2_dot_body1_n >= -1.0 * mEpsilonAngle   ) {
        return true; // Il-conditioned. Need to run GJK.

    }
    else if (body2_d1_dot_body1_n < body2_d2_dot_body1_n && 
             body2_n1_dot_body1_n < 0.0                     ) {

        // Face 1 is penetrating.
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2  = body2_fit1;
        return false;

    }
    else if (body2_n2_dot_body1_n < 0.0) {

        // Face 2 is penetrating.
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2  = body2_fit2;
        return false;

    }    
    else {
        return true; // Il-conditioned. Need to run GJK.

    }
}


bool ContactUpdaterFurtherChecker::checkFeatures_EDGE_FACE()
{
    auto body1_he1  = (*mInfo.mEit1)->he1();
    auto body1_he2  = (*mInfo.mEit1)->he2();
    auto body1_vit1 = (*body1_he1)->src();
    auto body1_vit2 = (*body1_he1)->dst();
    const auto body1_p1   = (*body1_vit1)->pGCS(mBody1.Qmat(), mBody1.CoM());
    const auto body1_p2   = (*body1_vit2)->pGCS(mBody1.Qmat(), mBody1.CoM());
    auto body1_v12  = body1_p2 - body1_p1;
    body1_v12.normalize();
    auto body1_fit1 = (*body1_he1)->face();
    auto body1_fit2 = (*body1_he2)->face();

    const auto body1_n1 = (*body1_fit1)->nGCS(mBody1.Qmat());
    const auto body1_n2 = (*body1_fit2)->nGCS(mBody1.Qmat());
    const auto body1_d1 = body1_n1.cross(body1_v12);
    const auto body1_d2 = body1_v12.cross(body1_n2);

    const auto body2_n  = (*mInfo.mFit2)->nGCS(mBody2.Qmat());

    const auto body1_d1_dot_body2_n = body1_d1.dot(body2_n);
    const auto body1_d2_dot_body2_n = body1_d2.dot(body2_n);

    const auto body1_n1_dot_body2_n = body1_n1.dot(body2_n);
    const auto body1_n2_dot_body2_n = body1_n2.dot(body2_n);

    if (body1_d1_dot_body2_n >= mEpsilonAngle && 
        body1_d2_dot_body2_n >= mEpsilonAngle   ) {

        return false; // Good.

    }
    else if (body1_d1_dot_body2_n <= 0.0 && 
             body1_d2_dot_body2_n <= 0.0   ) {
        return true; // Il-conditioned. Need to run GJK.
    }
    else if (body1_n1_dot_body2_n >= -1.0 * mEpsilonAngle && 
             body1_n2_dot_body2_n >= -1.0 * mEpsilonAngle   ) {
        return true; // Il-conditioned. Need to run GJK.

    }
    else if (body1_d1_dot_body2_n < body1_d2_dot_body2_n &&
             body1_n1_dot_body2_n < 0.0                     ) {

        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1  = body1_fit1;
        return false;

    }
    else if (body1_n2_dot_body2_n < 0.0) {

        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1  = body1_fit2;
        return false;

    }    
    else {
        return true; // Il-conditioned. Need to run GJK.

    }
}


bool ContactUpdaterFurtherChecker::checkFeatures_FACE_VERTEX()
{
    const auto body1_n = (*mInfo.mFit1)->nGCS(mBody1.Qmat());

    const auto body2_p = (*mInfo.mVit2)->pGCS(mBody2.Qmat(), mBody2.CoM());

    //           *previt
    //          /
    //         /
    //        v
    //       *-------------->* nextit
    //   body2_p   body2_he
    //
    //
    double maxDepth = 1.0;
    HalfEdgeIt maxHEit;
    long numFans = 0;
    long numInteriorFans = 0;
    for (auto body2_he : (*mInfo.mVit2)->halfEdges()) {

        if ((*body2_he)->src()==mInfo.mVit2) {
            numFans++;
            VertexIt vit = (*body2_he)->dst();
            auto edgeDir = (*vit)->pGCS(mBody2.Qmat(),mBody2.CoM()) - body2_p;
            edgeDir.normalize();
            const auto body1_n_dot_d = body1_n.dot(edgeDir);
 
            if (body1_n_dot_d <= mEpsilonAngle) {
                numInteriorFans++;
                if (body1_n_dot_d <= maxDepth) {
                    maxDepth = body1_n_dot_d;
                    maxHEit  = body2_he;
                }
            }
        }
    }

    if (numFans == numInteriorFans) {
        return true; // Total submersion. Need to run GJK.
    }

    if (maxDepth < 1.0) {
        // Check the neighbor edegs around the vertex.
        auto hePrev  = (*maxHEit)->prev();
        auto vPrev   = (*hePrev)->src();
        auto edgeDirPrev =(*vPrev)->pGCS(mBody2.Qmat(),mBody2.CoM()) - body2_p;
        edgeDirPrev.normalize();
        const auto body1_n_dot_dPrev = body1_n.dot(edgeDirPrev);
        auto fPrev   = (*maxHEit)->face();
        const auto body1_n_dot_nPrev = 
                       body1_n.dot((*fPrev)->nGCS(mBody2.Qmat()));
        auto heBuddy = (*maxHEit)->buddy();
        auto heNext  = (*heBuddy)->next();
        auto vNext   = (*heNext)->dst();
        auto edgeDirNext =(*vNext)->pGCS(mBody2.Qmat(),mBody2.CoM()) - body2_p;
        edgeDirNext.normalize();
        const auto body1_n_dot_dNext = body1_n.dot(edgeDirNext);
        auto fNext   = (*heBuddy)->face();
        const auto body1_n_dot_nNext = 
                       body1_n.dot((*fNext)->nGCS(mBody2.Qmat()));
        if (body1_n_dot_dPrev < body1_n_dot_dNext && 
            body1_n_dot_nPrev < mEpsilonAngle        ) {
            if (body1_n_dot_dPrev <= mEpsilonAngle) {

                // Face penetration
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit2  = fPrev;
            }
            else {
                // Edge penetration.
                mInfo.mType2 = ContactPairInfo::FT_EDGE;
                mInfo.mEit2  = (*maxHEit)->edge();
            }
        }
        else if (body1_n_dot_nNext < mEpsilonAngle) {
            if (body1_n_dot_dNext <= mEpsilonAngle) {
                // Face penetration
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit2  = fNext;
            }
            else {
                // Edge penetration.
                mInfo.mType2 = ContactPairInfo::FT_EDGE;
                mInfo.mEit2  = (*maxHEit)->edge();
            }
        }
        else {
            return true;
        }
    }
    return false;
}



bool ContactUpdaterFurtherChecker::checkFeatures_VERTEX_FACE()
{
    const auto body2_n = (*mInfo.mFit2)->nGCS(mBody2.Qmat());

    const auto body1_p = (*mInfo.mVit1)->pGCS(mBody1.Qmat(), mBody1.CoM());

    //           *previt
    //          /
    //         /
    //        v
    //       *-------------->* nextit
    //   body1_p   body1_he
    //
    //
    double maxDepth = 1.0;
    HalfEdgeIt maxHEit;
    long numFans = 0;
    long numInteriorFans = 0;
    for (auto body1_he : (*mInfo.mVit1)->halfEdges()) {

        if ((*body1_he)->src()==mInfo.mVit1) {
            numFans++;
            VertexIt vit = (*body1_he)->dst();
            auto edgeDir = (*vit)->pGCS(mBody1.Qmat(),mBody1.CoM()) - body1_p;
            edgeDir.normalize();
            const auto body2_n_dot_d = body2_n.dot(edgeDir);
 
            if (body2_n_dot_d <= mEpsilonAngle) {
                numInteriorFans++;
                if (body2_n_dot_d <= maxDepth) {
                    maxDepth = body2_n_dot_d;
                    maxHEit  = body1_he;
                }
            }
        }
    }

    if (numFans == numInteriorFans) {
        return true; // Total submersion. Need to run GJK.
    }

    if (maxDepth < 1.0) {
        // Check the neighbor edegs around the vertex.
        auto hePrev  = (*maxHEit)->prev();
        auto vPrev   = (*hePrev)->src();
        auto edgeDirPrev =(*vPrev)->pGCS(mBody1.Qmat(),mBody1.CoM()) - body1_p;
        edgeDirPrev.normalize();
        const auto body2_n_dot_dPrev = body2_n.dot(edgeDirPrev);
        auto fPrev   = (*maxHEit)->face();
        const auto body2_n_dot_nPrev = 
                       body2_n.dot((*fPrev)->nGCS(mBody1.Qmat()));
        auto heBuddy = (*maxHEit)->buddy();
        auto heNext  = (*heBuddy)->next();
        auto vNext   = (*heNext)->dst();
        auto edgeDirNext =(*vNext)->pGCS(mBody1.Qmat(),mBody1.CoM()) - body1_p;
        edgeDirNext.normalize();
        const auto body2_n_dot_dNext = body2_n.dot(edgeDirNext);
        auto fNext   = (*heBuddy)->face();
        const auto body2_n_dot_nNext = 
                       body2_n.dot((*fNext)->nGCS(mBody1.Qmat()));

        if (body2_n_dot_dPrev < body2_n_dot_dNext &&
            body2_n_dot_nPrev < mEpsilonAngle        ) {
            if (body2_n_dot_dPrev <= mEpsilonAngle) {
                // Face penetration
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = fPrev;
            }
            else {
                // Edge penetration.
                mInfo.mType1 = ContactPairInfo::FT_EDGE;
                mInfo.mEit1  = (*maxHEit)->edge();
            }
        }
        else if (body2_n_dot_nNext < mEpsilonAngle) {
            if (body2_n_dot_dNext <= mEpsilonAngle) {
                // Face penetration
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = fNext;
            }
            else {
                // Edge penetration.
                mInfo.mType1 = ContactPairInfo::FT_EDGE;
                mInfo.mEit1  = (*maxHEit)->edge();
            }
        }
        else {
            return true;
        }
    }
    return false;
}


bool ContactUpdaterFurtherChecker::checkFeatures_EDGE_EDGE()
{
    if (mInfo.mContactNormalType == ContactPairInfo::NT_EDGE_EDGE_AVG) {
        return checkFeatures_EDGE_EDGE_PARALLEL();
    }
    else { // mInfo.mContactNormalType == ContactPairInfo::NT_EDGE_CROSS_EDGE)
        return checkFeatures_EDGE_EDGE_CROSSING();
    }
}


bool ContactUpdaterFurtherChecker::checkFeatures_EDGE_EDGE_CROSSING()
{
    auto body1_he1  = (*mInfo.mEit1)->he1();
    auto body1_he2  = (*mInfo.mEit1)->he2();

    auto body1_fit1 = (*body1_he1)->face();
    auto body1_fit2 = (*body1_he2)->face();

    const auto body1_n1   = (*body1_fit1)->nGCS(mBody1.Qmat());
    const auto body1_n2   = (*body1_fit2)->nGCS(mBody1.Qmat());

    auto body1_vit1 = (*body1_he1)->src();
    auto body1_vit2 = (*body1_he1)->dst();
    const auto body1_p1   = (*body1_vit1)->pGCS(mBody1.Qmat(), mBody1.CoM());
    const auto body1_p2   = (*body1_vit2)->pGCS(mBody1.Qmat(), mBody1.CoM());
    auto body1_v12  = body1_p2 - body1_p1;
    body1_v12.normalize();

    const auto body1_d1 = body1_n1.cross(body1_v12);
    const auto body1_d2 = body1_v12.cross(body1_n2);

    const auto body1_en = (*mInfo.mEit1)->nGCS(mBody1.Qmat());

    auto body2_he1  = (*mInfo.mEit2)->he1();
    auto body2_he2  = (*mInfo.mEit2)->he2();

    auto body2_fit1 = (*body2_he1)->face();
    auto body2_fit2 = (*body2_he2)->face();

    const auto body2_n1   = (*body2_fit1)->nGCS(mBody2.Qmat());
    const auto body2_n2   = (*body2_fit2)->nGCS(mBody2.Qmat());

    auto body2_vit1 = (*body2_he1)->src();
    auto body2_vit2 = (*body2_he1)->dst();
    const auto body2_p1   = (*body2_vit1)->pGCS(mBody2.Qmat(), mBody2.CoM());
    const auto body2_p2   = (*body2_vit2)->pGCS(mBody2.Qmat(), mBody2.CoM());
    auto body2_v12  = body2_p2 - body2_p1;
    body2_v12.normalize();

    const auto body2_d1 = body2_n1.cross(body2_v12);
    const auto body2_d2 = body2_v12.cross(body2_n2);

    const auto body2_en = (*mInfo.mEit2)->nGCS(mBody2.Qmat());

    auto N1to2 = body1_v12.cross(body2_v12);
    N1to2.normalize();
    const auto directionReference = (*mInfo.mEit1)->nGCS(mBody1.Qmat()) - 
                                    (*mInfo.mEit2)->nGCS(mBody2.Qmat());
    if (N1to2.dot(directionReference) < 0.0) {
        N1to2.scale(-1.0);
    }

    // Check if either of the incident faces on body 1 penetrates into 
    // the edge on body 2.
    const auto pen_body1_d1 = N1to2.dot(body1_d1) * -1.0;
    const auto pen_body1_d2 = N1to2.dot(body1_d2) * -1.0;
    const auto pen_body2_d1 = N1to2.dot(body2_d1);
    const auto pen_body2_d2 = N1to2.dot(body2_d2);

    if ( pen_body1_d1 >= mEpsilonAngle &&
         pen_body1_d2 >= mEpsilonAngle &&
         pen_body2_d1 >= mEpsilonAngle &&
         pen_body2_d2 >= mEpsilonAngle    ) {
        return false; // Good. No penetration.
    }

    if ( (pen_body1_d1 <= mEpsilonAngle &&
          pen_body1_d2 <= mEpsilonAngle   ) ||
         (pen_body2_d1 <= mEpsilonAngle &&
          pen_body2_d2 <= mEpsilonAngle   )   ) {
        return true; // Il-conditioned. Need to run GJK.
    }

    if ( pen_body1_d1 <= mEpsilonAngle                 &&
         body1_n1.dot(body2_en) <= -1.0*mEpsilonAngle     ) {
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1 = body1_fit1;
    }
    else if ( pen_body1_d2 <= mEpsilonAngle                &&
              body1_n2.dot(body2_en) <= -1.0*mEpsilonAngle    ) {
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1 = body1_fit2;
    }
    else if (pen_body1_d1 <= mEpsilonAngle || pen_body1_d2 <= mEpsilonAngle ) {
        // The edge is penetrating but the normal is wrong. Run GJK.
        return true;
    }

    if ( pen_body2_d1 <= mEpsilonAngle               &&
         body2_n1.dot(body1_en) <= -1.0*mEpsilonAngle  ) {
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2 = body2_fit1;
    }
    else if ( pen_body2_d2 <= mEpsilonAngle                &&
              body2_n2.dot(body1_en) <= -1.0*mEpsilonAngle    ) {
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2 = body2_fit2;
    }
    else if (pen_body2_d1 <= mEpsilonAngle || pen_body2_d2 <= mEpsilonAngle ) {
        // The edge is penetrating but the normal is wrong. Run GJK.
        return true;
    }

    return false;
}


bool ContactUpdaterFurtherChecker::checkFeatures_EDGE_EDGE_PARALLEL()
{
    auto body1_he1  = (*mInfo.mEit1)->he1();
    auto body1_he2  = (*mInfo.mEit1)->he2();
    auto body1_fit1 = (*body1_he1)->face();
    auto body1_fit2 = (*body1_he2)->face();
    const auto body1_n1   = (*body1_fit1)->nGCS(mBody1.Qmat());
    const auto body1_n2   = (*body1_fit2)->nGCS(mBody1.Qmat());
    auto body1_vit1 = (*body1_he1)->src();
    auto body1_vit2 = (*body1_he1)->dst();
    const auto body1_p1   = (*body1_vit1)->pGCS(mBody1.Qmat(), mBody1.CoM());
    const auto body1_p2   = (*body1_vit2)->pGCS(mBody1.Qmat(), mBody1.CoM());
    auto body1_v12  = body1_p2 - body1_p1;
    body1_v12.normalize();
    const auto body1_d1   = body1_n1.cross(body1_v12);
    const auto body1_d2   = body1_v12.cross(body1_n2);

    auto body2_he1  = (*mInfo.mEit2)->he1();
    auto body2_he2  = (*mInfo.mEit2)->he2();
    auto body2_fit1 = (*body2_he1)->face();
    auto body2_fit2 = (*body2_he2)->face();
    const auto body2_n1   = (*body2_fit1)->nGCS(mBody2.Qmat());
    const auto body2_n2   = (*body2_fit2)->nGCS(mBody2.Qmat());
    auto body2_vit1 = (*body2_he1)->src();
    auto body2_vit2 = (*body2_he1)->dst();
    const auto body2_p1   = (*body2_vit1)->pGCS(mBody2.Qmat(), mBody2.CoM());
    const auto body2_p2   = (*body2_vit2)->pGCS(mBody2.Qmat(), mBody2.CoM());
    auto body2_v12  = body2_p2 - body2_p1;
    body2_v12.normalize();
    const auto body2_d1   = body2_n1.cross(body2_v12);
    const auto body2_d2   = body2_v12.cross(body2_n2);

    auto cr_11_21 = body1_d1.cross(body2_d1);
    auto cr_12_21 = body1_d2.cross(body2_d1);
    auto cr_11_22 = body1_d1.cross(body2_d2);
    auto cr_12_22 = body1_d2.cross(body2_d2);

    if (body1_v12.dot(body2_v12) < 0.0) {
        if ( body1_d1.dot(body2_d1)  >  0.0           && 
             cr_11_21.squaredNorm2() <= mEpsilonAngle   ) {
            //         .....d12
            //         ...../
            //  d11/d21<===*
            //         .....\
            //         .....d22
            if (body1_n1.dot(body2_n1) <= -1.0*mEpsilonAngle) {
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit1;
                mInfo.mFit2  = body2_fit1;
            }
            else {
                return true;
            }
        }
        else if ( body1_d2.dot(body2_d2)  >  0.0           && 
                  cr_12_22.squaredNorm2() <= mEpsilonAngle    ) {
            //         .....d21
            //         ...../
            //  d12/d22<===*
            //         .....\
            //         .....d11
            if (body1_n2.dot(body2_n2) <= -1.0*mEpsilonAngle) {
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit2;
                mInfo.mFit2  = body2_fit2;
            }
            else {
                return true;
            }
        }
        else if ( body1_d1.dot(body2_d2)  >  0.0           && 
                  cr_11_22.squaredNorm2() <= mEpsilonAngle    ) {
            if (isCCW(body1_d1, body1_d2, body2_d1, body1_v12)) {
                //         .d21.d12 
                //         ...\./
                // d11/d22 <===*
                if (body1_n1.dot(body2_n1) <= -1.0*mEpsilonAngle) {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit1;
                    mInfo.mFit2  = body2_fit1;
                }
                else {
                    return true;
                }
            }
            else {
                //         .d12.d21
                //         ...\./
                // d11/d22 <===*
                if (body1_n2.dot(body2_n2) <= -1.0*mEpsilonAngle) {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit2;
                    mInfo.mFit2  = body2_fit2;
                }
                else {
                    return true;
                }
            }
        }
        else if ( body1_d2.dot(body2_d1)  >                0.0 && 
                  cr_12_21.squaredNorm2() <= mEpsilonAngle        ) {
            if (isCCW(body1_d1, body1_d2, body2_d2, body1_v12)) {
                // d12/d21 <===*
                //         .../.\
                //         .d22.d11
                if (body1_n2.dot(body2_n2) <= -1.0*mEpsilonAngle) {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit2;
                    mInfo.mFit2  = body2_fit2;
                }
                else {
                    return true;
                }
            }
            else {
                // d12/d21 <===*
                //         .../.\
                //         .d11.d22
                if (body1_n1.dot(body2_n1) <= -1.0*mEpsilonAngle) {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit1;
                    mInfo.mFit2  = body2_fit1;
                }
                else {
                    return true;
                }
            }
        }
        else {
            bool ccw_11_12_21 = isCCW(body1_d1, body1_d2, body2_d1, body1_v12);
            bool ccw_11_12_22 = isCCW(body1_d1, body1_d2, body2_d2, body1_v12);
            bool ccw_11_21_22 = isCCW(body1_d1, body2_d1, body2_d2, body1_v12);
            bool ccw_12_21_22 = isCCW(body1_d2, body2_d1, body2_d2, body1_v12);

            if ( !ccw_11_12_21 && !ccw_11_12_22 && 
                  ccw_11_21_22 &&  ccw_12_21_22    ) {
                // d11->d21->d22->d12  ( Separate)
                //
                //           d22 d21
                //             \ /
                //              *
                //             / \
                //          d12   d11
                ;
            }
            else if (  ccw_11_12_21 &&  ccw_11_12_22 && 
                       ccw_11_21_22 &&  ccw_12_21_22    ) {
                // d11->d12->d21->d22  ( Body 2 is completerly inside Body 1)
                //
                //            _-*-_
                //          _- / \ -_
                //      d12-  /   \  -d11
                //          d21   d22
                if (body1_d2.dot(body2_d1) > body2_d2.dot(body1_d1) ) {
                    if (body1_n2.dot(body2_n2) <= -1.0*mEpsilonAngle) {
                        mInfo.mType1 = ContactPairInfo::FT_FACE;
                        mInfo.mType2 = ContactPairInfo::FT_FACE;
                        mInfo.mFit1  = body1_fit2;
                        mInfo.mFit2  = body2_fit2;
                    }
                    else {
                        return true;
                    }
                }
                else {
                    if (body1_n1.dot(body2_n1) <= -1.0*mEpsilonAngle) {
                        mInfo.mType1 = ContactPairInfo::FT_FACE;
                        mInfo.mType2 = ContactPairInfo::FT_FACE;
                        mInfo.mFit1  = body1_fit1;
                        mInfo.mFit2  = body2_fit1;
                    }
                    else {
                        return true;
                    }
                }
            }
            else if ( !ccw_11_12_21 &&  ccw_11_12_22 && 
                       ccw_11_21_22 && !ccw_12_21_22    ) {
                // d11->d21->d12->d22  ( d22 is penetrating into d12)
                //
                //            _-*-_
                //          _- / \ -_
                //      d21-  /   \  -d11
                //          d12   d22
                if (body1_n2.dot(body2_n2) <= -1.0*mEpsilonAngle) {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit2;
                    mInfo.mFit2  = body2_fit2;
                }
                else {
                    return true;
                }
            }
            else if (  ccw_11_12_21 && !ccw_11_12_22 && 
                      !ccw_11_21_22 &&  ccw_12_21_22    ) {
                // d11->d22->d12->d21  (d11 is penetrating into d21)
                //
                //            _-*-_
                //          _- / \ -_
                //      d12-  /   \  -d22
                //          d21   d11
                if (body1_n1.dot(body2_n1) <= -1.0*mEpsilonAngle) {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit1;
                    mInfo.mFit2  = body2_fit1;
                }
                else {
                    return true;
                }
            }
            else if ( !ccw_11_12_21 && !ccw_11_12_22 && 
                      !ccw_11_21_22 && !ccw_12_21_22    ) {
                // d11->d22->d21->d12  ( Body 1 is completerly inside Body 2)
                //
                //            _-*-_
                //          _- / \ -_
                //      d21-  /   \  -d22
                //          d12   d11
                if (body2_d1.dot(body1_d2) > body1_d1.dot(body2_d2) ) {
                    if (body1_n1.dot(body2_n1) <= -1.0*mEpsilonAngle) {
                        mInfo.mType1 = ContactPairInfo::FT_FACE;
                        mInfo.mType2 = ContactPairInfo::FT_FACE;
                        mInfo.mFit1  = body1_fit1;
                        mInfo.mFit2  = body2_fit1;
                    }
                    else {
                        return true;
                    }
                }
                else {
                    if (body1_n2.dot(body2_n2) <= -1.0*mEpsilonAngle) {
                        mInfo.mType1 = ContactPairInfo::FT_FACE;
                        mInfo.mType2 = ContactPairInfo::FT_FACE;
                        mInfo.mFit1  = body1_fit2;
                        mInfo.mFit2  = body2_fit2;
                    }
                    else {
                        return true;
                    }
                }
            }
            else {
                // d11->d12->d22->d21  (Forbidden)
                log(ERROR, __FILE__, __LINE__, "Forbidden predicates");
                return true;
            }
        }
    }
    else {
        if ( body1_d1.dot(body2_d2)  >  0.0           && 
             cr_11_22.squaredNorm2() <= mEpsilonAngle    ) {
            //         .....d12
            //         ...../
            //  d11/d22<===*
            //         .....\
            //         .....d21
            if (body1_n1.dot(body2_n2) <= -1.0*mEpsilonAngle) {
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit1;
                mInfo.mFit2  = body2_fit2;
            }
            else {
                return true;
            }
        }
        else if ( body1_d2.dot(body2_d1)  >  0.0           && 
                  cr_12_21.squaredNorm2() <= mEpsilonAngle    ) {
            //         .....d22
            //         ...../
            //  d12/d21<===*
            //         .....\
            //         .....d11
            if (body1_n2.dot(body2_n1) <= -1.0*mEpsilonAngle) {
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit2;
                mInfo.mFit2  = body2_fit1;
            }
            else {
                return true;
            }
        }
        else if ( body1_d1.dot(body2_d1)  >  0.0           && 
                  cr_11_21.squaredNorm2() <= mEpsilonAngle    ) {
            if (isCCW(body1_d1, body1_d2, body2_d2, body1_v12)) {
                //         .d22.d12 
                //         ...\./
                // d11/d21 <===*
                if (body1_n1.dot(body2_n2) <= -1.0*mEpsilonAngle) {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit1;
                    mInfo.mFit2  = body2_fit2;
                }
                else {
                    return true;
                }
            }
            else {
                //         .d12.d22
                //         ...\./
                // d11/d21 <===*
                if (body1_n2.dot(body2_n1) <= -1.0*mEpsilonAngle) {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit2;
                    mInfo.mFit2  = body2_fit1;
                }
                else {
                    return true;
                }
            }
        }
        else if ( body1_d2.dot(body2_d2)  >  0.0           && 
                  cr_12_22.squaredNorm2() <= mEpsilonAngle    ) {
            if (isCCW(body1_d1, body1_d2, body2_d1, body1_v12)) {
                // d12/d22 <===*
                //         .../.\
                //         .d21.d11
                if (body1_n2.dot(body2_n1) <= -1.0*mEpsilonAngle) {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit2;
                    mInfo.mFit2  = body2_fit1;
                }
                else {
                    return true;
                }
            }
            else {
                // d12/d22 <===*
                //         .../.\
                //         .d11.d21
                if (body1_n1.dot(body2_n2) <= -1.0*mEpsilonAngle) {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit1;
                    mInfo.mFit2  = body2_fit2;
                }
                else {
                    return true;
                }
            }
        }
        else {
            bool ccw_11_12_22 = isCCW(body1_d1, body1_d2, body2_d2, body1_v12);
            bool ccw_11_12_21 = isCCW(body1_d1, body1_d2, body2_d1, body1_v12);
            bool ccw_11_22_21 = isCCW(body1_d1, body2_d2, body2_d1, body1_v12);
            bool ccw_12_22_21 = isCCW(body1_d2, body2_d2, body2_d1, body1_v12);

            if ( !ccw_11_12_22 && !ccw_11_12_21 && 
                  ccw_11_22_21 &&  ccw_12_22_21    ) {
                // d11->d22->d21->d12  ( Separate)
                //
                //           d21 d22
                //             \ /
                //              *
                //             / \
                //          d12   d11
                ;
            }
            else if (  ccw_11_12_22 &&  ccw_11_12_21 && 
                       ccw_11_22_21 &&  ccw_12_22_21    ) {
                // d11->d12->d22->d21  ( Body 2 is completerly inside Body 1)
                //
                //            _-*-_
                //          _- / \ -_
                //      d12-  /   \  -d11
                //          d22   d21
                if (body1_d2.dot(body2_d2) > body2_d1.dot(body1_d1) ) {
                    if (body1_n2.dot(body2_n1) <= -1.0*mEpsilonAngle) {
                        mInfo.mType1 = ContactPairInfo::FT_FACE;
                        mInfo.mType2 = ContactPairInfo::FT_FACE;
                        mInfo.mFit1  = body1_fit2;
                        mInfo.mFit2  = body2_fit1;
                    }
                    else {
                        return true;
                    }
                }
                else {
                    if (body1_n1.dot(body2_n2) <= -1.0*mEpsilonAngle) {
                        mInfo.mType1 = ContactPairInfo::FT_FACE;
                        mInfo.mType2 = ContactPairInfo::FT_FACE;
                        mInfo.mFit1  = body1_fit1;
                        mInfo.mFit2  = body2_fit2;
                    }    
                    else {
                        return true;
                    }
                }
            }
            else if ( !ccw_11_12_22 &&  ccw_11_12_21 && 
                       ccw_11_22_21 && !ccw_12_22_21    ) {
                // d11->d22->d12->d21  ( d21 is penetrating into d12)
                //
                //            _-*-_
                //          _- / \ -_
                //      d22-  /   \  -d11
                //          d12   d21
                if (body1_n2.dot(body2_n1) <= -1.0*mEpsilonAngle) {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit2;
                    mInfo.mFit2  = body2_fit1;
                }
                else {
                    return true;
                }
            }
            else if (  ccw_11_12_22 && !ccw_11_12_21 && 
                      !ccw_11_22_21 &&  ccw_12_22_21    ) {
                // d11->d21->d12->d22  (d11 is penetrating into d22)
                //
                //            _-*-_
                //          _- / \ -_
                //      d12-  /   \  -d21
                //          d22   d11
                if (body1_n1.dot(body2_n2) <= -1.0*mEpsilonAngle) {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit1;
                    mInfo.mFit2  = body2_fit2;
                }
                else {
                    return true;
                }
            }
            else if ( !ccw_11_12_22 && !ccw_11_12_21 && 
                      !ccw_11_22_21 && !ccw_12_22_21    ) {
                // d11->d21->d22->d12  ( Body 1 is completerly inside Body 2)
                //
                //            _-*-_
                //          _- / \ -_
                //      d22-  /   \  -d21
                //          d12   d11
                if (body2_d2.dot(body1_d2) > body1_d1.dot(body2_d1) ) {
                    if (body1_n1.dot(body2_n2) <= -1.0*mEpsilonAngle) {
                        mInfo.mType1 = ContactPairInfo::FT_FACE;
                        mInfo.mType2 = ContactPairInfo::FT_FACE;
                        mInfo.mFit1  = body1_fit1;
                        mInfo.mFit2  = body2_fit2;
                    }
                    else {
                        return true;
                    }
                }
                else {
                    if (body1_n2.dot(body2_n1) <= -1.0*mEpsilonAngle) {
                        mInfo.mType1 = ContactPairInfo::FT_FACE;
                        mInfo.mType2 = ContactPairInfo::FT_FACE;
                        mInfo.mFit1  = body1_fit2;
                        mInfo.mFit2  = body2_fit1;
                    }
                    else {
                        return true;
                    }
                }
            }
            else {
                // d11->d12->d21->d22  (Forbidden)
                log(ERROR, __FILE__, __LINE__, "Forbidden predicates");
                return true;
            }
        }
    }
    return false;
}


bool ContactUpdaterFurtherChecker::isCCW(
    const Vec3& v1,
    const Vec3& v2,
    const Vec3& v3,
    const Vec3& axis
) {
    auto cr13 = v1.cross(v3);
    auto cr21 = v2.cross(v1);
    auto cr32 = v3.cross(v2);

    if (cr13.squaredNorm2() <= mEpsilonAngle) {
        if (v1.dot(v3) < 0.0 ) {
            return axis.dot(v1.cross(v2))>= 0.0 &&
                   axis.dot(v2.cross(v3))>= 0.0 ;
        }
        else {
            return true;
        }
    }
    else if (cr21.squaredNorm2() <= mEpsilonAngle) {
        if (v2.dot(v1) <  0.0 ) { 
            return axis.dot(v2.cross(v3))>= 0.0 &&
                   axis.dot(v3.cross(v1))>= 0.0 ;
        }
        else {
            return true;
        }
    }
    else if (cr32.squaredNorm2() <= mEpsilonAngle) {
        if (v3.dot(v2) <  0.0 ) { 
            return axis.dot(v3.cross(v1))>= 0.0 &&
                   axis.dot(v1.cross(v2))>= 0.0 ;
        }
        else {
            return true;
        }
    }
    else if (axis.dot(cr13) >= 0.0) {
        return axis.dot(v1.cross(v2))>= 0.0 &&
               axis.dot(v2.cross(v3))>= 0.0 ;
    }
    else if (axis.dot(cr21) >= 0.0) {
        return axis.dot(v2.cross(v3))>= 0.0 &&
               axis.dot(v3.cross(v1))>= 0.0 ;
    }
    else {
        return true;
    }
}



bool ContactUpdaterFurtherChecker::isRayInsideCone(
    const Vec3& dir,
    VertexIt    vitApex,
    bool        coneOnBody1 
) {
    for (auto hePillar : (*vitApex)->halfEdges()) {
        if ((*hePillar)->src()==vitApex) {
            auto fit = (*hePillar)->face();
            auto fn  = (*fit)->nGCS(coneOnBody1?mBody1.Qmat():mBody2.Qmat());
            if (dir.dot(fn) >= mEpsilonAngle) {
                return false;
            }
        }
    }
    return true;
}


void ContactUpdaterFurtherChecker::rotateWorld(
    const Vec3& yDir,
    const Vec3& zDir,
    Mat3x3&     rotMat1,
    Vec3&       com1,
    Mat3x3&     rotMat2,
    Vec3&       com2
) {
    Vec3 xDir = yDir.cross(zDir);
    xDir.normalize();
    Mat3x3 rMat(xDir, yDir, zDir);
    rMat.transposeInPlace();

    rotMat1 = rMat * mBody1.Qmat();
    com1    = rMat * mBody1.CoM();
    rotMat2 = rMat * mBody2.Qmat();
    com2    = rMat * mBody2.CoM();
}


void ContactUpdaterFurtherChecker::rotateWorld(
    const Vec3& zDir,
    Mat3x3&     rotMat1,
    Vec3&       com1,
    Mat3x3&     rotMat2,
    Vec3&       com2
) {
    Mat3x3 rMat = rotMatAlignZDirToZAxis(zDir);

    rotMat1 = rMat * mBody1.Qmat();
    com1    = rMat * mBody1.CoM();
    rotMat2 = rMat * mBody2.Qmat();
    com2    = rMat * mBody2.CoM();
}


Mat3x3 ContactUpdaterFurtherChecker::rotMatAlignZDirToZAxis(const Vec3& zDir)
{
    // Find the new coordinate axes with zDir for new z-axis.
    Vec3 r3 = zDir;
    double ax = fabs(r3.x());
    double ay = fabs(r3.y());
    double az = fabs(r3.z());

    Vec3 r1;
    if (ax <= ay && ay <= az) {
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



bool ContactUpdaterFurtherChecker::findIntersection(
    const Vec2& p1,
    const Vec2& p2,
    const Vec2& p3,
    const Vec2& p4,
    Vec2&       intsec
) {
    const double denom = (p1.x() - p2.x())*(p3.y() - p4.y()) -
                         (p1.y() - p2.y())*(p3.x() - p4.x()) ;

    if (fabs(denom)<=mEpsilonZero) {
        return false;
    }

    const double x1y2_y1x2 = p1.x()*p2.y() - p1.y()*p2.x();
    const double x3y4_y3x4 = p3.x()*p4.y() - p3.y()*p4.x();

    intsec.setX((x1y2_y1x2*(p3.x()-p4.x())-(p1.x()-p2.x())*x3y4_y3x4)/denom);
    intsec.setY((x1y2_y1x2*(p3.y()-p4.y())-(p1.y()-p2.y())*x3y4_y3x4)/denom);

    return true;
}


/** @brief check for further penetrations for EDGE_VERTEX cases.
 *
 *  This rarely occurs in reality and the penetration check is complex.
 *  Let W be the wedge formed by the two incident faces of the active edge
 *  on Body1.
 *
 *                      
 *                                  _-*
 *          Active edge on Body 1 _-   \
 *         (mInfo.Eit1)        _-       \
 *                           *-          \
 *                          /.\           .
 *                         /...\           .
 *                        /..W..\
 *                       .........
 *                      ...........
 *
 *  Let C be the circular composition of surrounding edges around the apex
 *  vertex A (active vertex, mInfo.Vit2).
 *  Please note that those surrounding edges are not directly incident to
 *  the apex.
 *                             C
 *                           *---*-*
 *                          / . .   \
 *                         *...*.....*
 *                         |  .A.   /
 *                         \ .   . /
 *                          *-----*
 *
 *  The problem of finding the penetration of two Bodies if they are anchored
 *  at the edge and the vertex is reduced to an 2D intersection problem as
 *  follows.
 *  We project the features of W, A, and C onto the 2D plane that perpendicular
 *  to the active edge on Body 1.
 *  W becomes an open triangle with its apex at A, and C is a closed curve,
 *  which may not be convex, plane, or simple.
 *
 *              ________
 *             /        \
 *             \______   \C
 *              A*    \   \
 *              / \   /   /
 *             /   \  \  /
 *            /  W  \  \/
 *           .       .
 *          .         .
 *
 *  If W and C do not intersect or touch, then there will be no penetration
 *  between Body 1 and Body2 except for the active edge and the apex vertex.
 *  If they intersect or touch with each other, then we classify it into 
 *  one of the following types.
 *
        

 *  Type  : WC_APEX_COINCIDENT_VERTEX,
 *          A vertex of C (a pillar edge of the cone on Body 2) is coincident
 *          to A. Promote the active feature of Body 2 to the pillar edge.
 *
 *                         /
 *                        /
 *                  -----o*
 *                      / \
 *                     /   \
 *                    /     \
 *
 *
 *  Type  : WC_APEX_INSIDE_VERTICES,
 *          C intersects with both sides of W and there are one or multiple
 *          vertex of C inside W. Promote the active feature of Body 2 to
 *          a pillar edge inside W.
 *
 *                          *
 *                         / \
 *                        /   \    .
 *                    .  /     \  /
 *                     \/       \/
 *                     /\       /\
 *                    /  .--o--.  \
 *
 *  Type  : WC_APEX_ONE_EDGE_TOUCHING or WC_APEX_ONE_EDGE_CUT_ACROSS
 *          An edge of C (face of the cone on Body 2) touchies W at A or
 *          cuts across W on both sides.
 *          Prompte the active feature of Body 2 to the face.
 *
 *                    o--*--o         *
 *                      / \          / \
 *                     /   \     o--+---+--o
 *                    /     \      /     \
 *
 *
 *  Type  : WC_EDGE1_INSIDE_VERTICES or WC_EDGE2_INSIDE_VERTICES
 *          C intersects on one side of W and there are one or more vertices
 *          of C inside W.
 *          Promote the active features of Body 1 and 2 to the faces.
 *          Pick the face on Body 2 most aligned with the wedge face.
 *          
 *                        *
 *                       / \|
 *                      /   +
 *                     /    |\
 *                    /     o \
 *                   /       \ \
 *                  /         o \
 *                             \ \
 *                              o-+--
 *                                 \
 *                                  \
 *
 *  Type  : WC_EDGE1_TOUCHING or WC_EDGE2_TOUCHING,
 *          An edge of C touches W on a side of W.
 *          Promote the active features of Body 1 and 2 to the touching faces.
 *
 *                       *               *o           o 
 *                      / \o            / \\           \
 *                     /   \\          /   \\          *\
 *                    /     \\        /     \o        / \\
 *                   /       \o      /       \       /   \o
 *
 *          A vertex of C touching on a side of W.
 *          Promote the active feature of Body 1 to the face, and the active
 *          feature of Body 2 to the edge.
 *
 *                       * |
 *                      / \|
 *                     /   o--
 *                    /     \
 *                   /       \
 *                  /         \
 *
 */
bool ContactUpdaterFurtherChecker::checkFeatures_EDGE_VERTEX()
{

    auto body1_he1  = (*(mInfo.mEit1))->he1();
    auto body1_he2  = (*(mInfo.mEit1))->he2();
    auto body1_vit1 = (*body1_he1)->src();
    auto body1_vit2 = (*body1_he1)->dst();
    auto body1_p1   = (*body1_vit1)->pGCS(mBody1.Qmat(), mBody1.CoM());
    auto body1_p2   = (*body1_vit2)->pGCS(mBody1.Qmat(), mBody1.CoM());
    auto body1_v12  = body1_p2 - body1_p1;
    body1_v12.normalize();
    auto body1_n    = (*mInfo.mEit1)->nGCS(mBody1.Qmat());
    // Rotate the bodies such that the edge normal (body1_n) is aligned with
    // Y-axis and the edge direction (from he1's src to dst, body1_v12) with 
    // Z-axis.
    Mat3x3 rotMat1, rotMat2;
    Vec3   com1,    com2;
    rotateWorld(body1_n, body1_v12, rotMat1, com1, rotMat2, com2);

    body1_p1        = (*body1_vit1)->pGCS(rotMat1, com1);
    body1_p2        = (*body1_vit2)->pGCS(rotMat1, com1);
    body1_v12       = body1_p2 - body1_p1;
    body1_v12.normalize();
    auto body1_fit1 = (*body1_he1)->face();
    auto body1_fit2 = (*body1_he2)->face();
    auto body1_fn1  = (*body1_fit1)->nGCS(rotMat1);
    auto body1_fn2  = (*body1_fit2)->nGCS(rotMat1);

    auto body1_dn1  = body1_fn1.cross(body1_v12);
    auto body1_dn2  = body1_v12.cross(body1_fn2);
    // Project the features onto the new XY-plane.
    vector<Vec2>       body2_fan;
    vector<HalfEdgeIt> body2_hes;
    auto body2_p   = (*mInfo.mVit2)->pGCS(rotMat2, com2);
    for (auto body2_he : (*mInfo.mVit2)->halfEdges()) {
        if ((*body2_he)->src()==mInfo.mVit2) {
            auto body2_vitFan = (*body2_he)->dst();
            auto body2_vFan   = (*body2_vitFan)->pGCS(rotMat2, com2) -
                                body2_p;
            body2_vFan.normalize();
            Vec2 p_2D(body2_p.x() + body2_vFan.x(), body2_p.y() + body2_vFan.y());
            body2_fan.push_back(p_2D);
            body2_hes.push_back(body2_he);
        }            
    }
    Vec2 apex           (body1_p1.x(),         body1_p1.y());
    Vec2 body1_side1    (body1_dn1.x(),        body1_dn1.y());
    body1_side1.normalize();
    Vec2 body1_side1perp(-1.0 * body1_dn1.y(), body1_dn1.x());
    Vec2 body1_side2    (body1_dn2.x(),        body1_dn2.y());
    body1_side2.normalize();
    Vec2 body1_side2perp(body1_dn2.y(),        -1.0*body1_dn2.x());
    // Classify the vertex w.r.t the wedge.
    vector<enum WEDGE_CLASSIFIER> predMain;
    vector<enum WEDGE_CLASSIFIER> predAux;

    bool intersecting = classifyPointsAgainstWedge(
                            body2_fan,
                            apex,
                            body1_side1,
                            body1_side1perp,
                            body1_side2,
                            body1_side2perp,
                            predMain,
                            predAux
                        );
    if (!intersecting) {
        return false;
    }
    auto parseResult = parseClassifiersPointsAgainstWedge(predMain, predAux);
    return updateFeatures_EDGE_VERTEX(
        parseResult,
        predMain,
        predAux,
        body2_hes,
        apex,
        body2_fan,
        body1_side1perp,
        body1_side2perp,
        body1_fit1,
        body1_fit2,
        true
    );
}


bool ContactUpdaterFurtherChecker::classifyPointsAgainstWedge(
    vector<Vec2>&                  fan,
    const Vec2&                    apex,
    const Vec2&                    side1,
    const Vec2&                    side1perp,
    const Vec2&                    side2,
    const Vec2&                    side2perp,
    vector<enum WEDGE_CLASSIFIER>& predMain,
    vector<enum WEDGE_CLASSIFIER>& predAux
) {
    predMain.clear();
    predAux.clear();

    bool intersecting = false;
    // Classify the points on the fan.
    for (auto& f : fan) {
        auto v = f - apex;
        if (v.squaredNorm2() <= mEpsilonZero) {
            predMain.push_back(WC_ON_APEX);
            intersecting = true;
        }
        else if (fabs(v.dot(side1perp)) <= mEpsilonAngle && 
                 v.dot(side1) > 0.0                         ) {
            predMain.push_back(WC_ON_EDGE1);
            intersecting = true;
        }
        else if (fabs(v.dot(side2perp)) <= mEpsilonAngle &&
                 v.dot(side2) > 0.0                         ) {
            predMain.push_back(WC_ON_EDGE2);
            intersecting = true;
        }
        else if (v.dot(side1perp) < 0.0 && v.dot(side2perp) < 0.0 &&
                 v.dot(side1)     > 0.0     && v.dot(side2) > 0.0    ) {
            predMain.push_back(WC_INSIDE);
            intersecting = true;
        }
        else {
            predMain.push_back(WC_OUTSIDE);
        }
        predAux.push_back(WC_NONE);
    }

    // Further classify the intersections.
    auto numElems = predMain.size();
    for (long i = 0; i < numElems; i++) {
        auto i2 = (i+1)%numElems;
        auto pred1 = predMain[i];
        auto pred2 = predMain[i2];

        auto v12 = fan[i2] - fan[i];
        if (v12.squaredNorm2() <= mEpsilonZero) {
            if (pred1 == WC_ON_EDGE1 || pred2 == WC_ON_EDGE1) {
                predAux[i] = WC_EDGE1_TOUCHING;
            }
            else if (pred1 == WC_ON_EDGE2 || pred2 == WC_ON_EDGE2) {
                predAux[i] = WC_EDGE2_TOUCHING;
            }
        }
        else {
            v12.normalize();
            if ( (pred1 == WC_INSIDE  && pred2 == WC_OUTSIDE)||
                 (pred1 == WC_OUTSIDE && pred2 == WC_INSIDE )  ) {
                Vec2 intsec1, intsec2;
                if (findIntersection(
                                  fan[i], fan[i2], apex, apex + side1, intsec1)) {
                    auto vi1a = intsec1 - apex;
                    auto v1i = intsec1 - fan[i];
                    auto v2i = intsec1 - fan[i2];
    
                    if (vi1a.squaredNorm2() <= mEpsilonZero) {
                        predAux[i] = WC_INTSEC_ON_APEX;
                    }
                    else if (vi1a.dot(side1) > 0.0             &&
                             v12.dot(v1i) >= -1.0*mEpsilonZero &&
                             v12.dot(v2i) <= mEpsilonZero        ) {
                        predAux[i] = WC_INTSEC_ON_EDGE1;
                    }
                    else if (findIntersection(
                                  fan[i], fan[i2], apex, apex + side2, intsec2)) {
                        auto vi2a = intsec2 - apex;
                        v1i = intsec2 - fan[i];
                        v2i = intsec2 - fan[i2];
    
                        if (vi2a.squaredNorm2()<=mEpsilonZero) {
                            predAux[i] = WC_INTSEC_ON_APEX;
                        }
                        else if (vi2a.dot(side2) > 0.0             &&
                                 v12.dot(v1i) >= -1.0*mEpsilonZero &&
                                 v12.dot(v2i) <= mEpsilonZero         ) {
                            predAux[i] = WC_INTSEC_ON_EDGE2;
                        }
                    }
                }
                else if (findIntersection(
                                  fan[i], fan[i2], apex, apex + side2, intsec2)) {
                    auto v = intsec2 - apex;
                    if (v.squaredNorm2() <= mEpsilonZero) {
                        predAux[i] = WC_INTSEC_ON_APEX;
                    }
                    else if (v.dot(side2) > 0.0) {
                        predAux[i] = WC_INTSEC_ON_EDGE2;
                    }
                }
            }
            else if (pred1 == WC_OUTSIDE && pred2 == WC_OUTSIDE) {
                Vec2 intsec1;
                if (findIntersection(
                                  fan[i], fan[i2], apex, apex + side1, intsec1)) {
                    auto v = intsec1 - apex;
                    if (v.squaredNorm2()<=mEpsilonZero) {
                        predAux[i] = WC_INTSEC_ON_APEX;
                        intersecting = true;
                    }
                    else if (v.dot(side1) > 0.0) {
                        auto v1i = intsec1 - fan[i];
                        if (v12.dot(v1i) >= -1.0*mEpsilonZero) {
                            auto v2i = intsec1 - fan[i2];
                            if (v12.dot(v2i) <= mEpsilonZero) {
                                predAux[i] = WC_INTSEC_ON_EDGE12;
                                intersecting = true;
                            }
                        }
                    }
                }
            }
            else if (pred1 == WC_ON_EDGE1||pred2 == WC_ON_EDGE1) {
                if (predAux[i]==WC_NONE) {
                    if (pred1 == WC_ON_EDGE1 && pred2 == WC_ON_EDGE1) {
                        predAux[i] = WC_EDGE1_TOUCHING;
                    }
                    else if (pred1 == WC_ON_APEX || pred2 == WC_ON_APEX) {
                        predAux[i] = WC_EDGE1_TOUCHING;
                    }
                    else if (pred1 == WC_OUTSIDE || pred2 == WC_OUTSIDE) {
                        if (fabs(side1perp.dot(v12))<=mEpsilonAngle) {
                            predAux[i] = WC_EDGE1_TOUCHING;
                        }
                        else {
                            Vec2 intsec2;
                            if (findIntersection(
                                  fan[i], fan[i2], apex, apex + side2, intsec2)) {
                                auto v = intsec2 - apex;
                                if (v.dot(side2) > 0.0) {
                                   auto v1i = intsec2 - fan[i];
                                   if (v12.dot(v1i) >= -1.0*mEpsilonZero) {
                                       auto v2i = intsec2 - fan[i2];
                                       if (v12.dot(v2i) <= mEpsilonZero) {
                                           predAux[i] = WC_INTSEC_ON_EDGE2;
                                       }
                                   }
                               }
                           }
                        }
                    }
                }            
            }
            else if (pred1 == WC_ON_EDGE2||pred2 == WC_ON_EDGE2) {
                if (predAux[i]==WC_NONE) {
                    if (pred1 == WC_ON_EDGE2 && pred2 == WC_ON_EDGE2) {
                        predAux[i] = WC_EDGE2_TOUCHING;
                    }
                    else if (pred1 == WC_ON_APEX || pred2 == WC_ON_APEX) {
                        predAux[i] = WC_EDGE2_TOUCHING;
                    }
                    else if (pred1 == WC_OUTSIDE || pred2 == WC_OUTSIDE) {
                        if (fabs(side2perp.dot(v12))<=mEpsilonAngle) {
                            predAux[i] = WC_EDGE2_TOUCHING;
                        }
                        else {
                            Vec2 intsec1;
                            if (findIntersection(
                                  fan[i], fan[i2], apex, apex + side1, intsec1)) {
                                auto v = intsec1 - apex;
                                if (v.dot(side1) > 0.0) {
                                   auto v1i = intsec1 - fan[i];
                                   if (v12.dot(v1i) >= -1.0*mEpsilonZero) {
                                       auto v2i = intsec1 - fan[i2];
                                       if (v12.dot(v2i) <= mEpsilonZero) {
                                           predAux[i] = WC_INTSEC_ON_EDGE1;
                                       }
                                   }
                               }
                           }
                        }
                    }
                }            
            }
        }
    }
    return intersecting;
}


enum ContactUpdaterFurtherChecker::WEDGE_CLASSIFIER 
ContactUpdaterFurtherChecker::parseClassifiersPointsAgainstWedge(
    vector<enum WEDGE_CLASSIFIER>& predMain,
    vector<enum WEDGE_CLASSIFIER>& predAux
) {
    long cntON_APEX          = 0;
    long cntON_EDGE1         = 0;
    long cntON_EDGE2         = 0;
    long cntINSIDE           = 0;
    long cntOUTSIDE          = 0;
    long cntINTSEC_ON_APEX   = 0;
    long cntINTSEC_ON_EDGE1  = 0;
    long cntINTSEC_ON_EDGE2  = 0;
    long cntINTSEC_ON_EDGE12 = 0;
    long cntAux_NONE         = 0;

    const long numElems = predMain.size();
    // Scan the predicates and determine the penetration.
    for (long i = 0; i < numElems; i++) {
        switch (predMain[i]) {
          case WC_ON_APEX:
            cntON_APEX++;
            break;
          case WC_ON_EDGE1:
            cntON_EDGE1++;
            break;
          case WC_ON_EDGE2:
            cntON_EDGE2++;
            break;
          case WC_INSIDE:
            cntINSIDE++;
            break;
          case WC_OUTSIDE:
            cntOUTSIDE++;
            break;
          default:
            break;
        }
    }

    if (cntOUTSIDE==0) {
        return WC_COMPLETELY_INSIDE;
    }

    for (long i = 0; i < numElems; i++) {
        switch (predAux[i]) {
          case WC_NONE:
            cntAux_NONE++;
            break;
          case WC_INTSEC_ON_APEX:
            cntINTSEC_ON_APEX++;
            break;
          case WC_INTSEC_ON_EDGE12:
            cntINTSEC_ON_EDGE12++;
            break;
          case WC_INTSEC_ON_EDGE1:
            cntINTSEC_ON_EDGE1++;
            break;
          case WC_INTSEC_ON_EDGE2:
            cntINTSEC_ON_EDGE2++;
            break;
          default:
            break;
        }
    }

    if (cntOUTSIDE==0 && cntAux_NONE==numElems) {
        return WC_NONE;// Early out
    }

    if ( cntON_EDGE1 + cntINTSEC_ON_EDGE1 >  0 &&
         cntON_EDGE2 + cntINTSEC_ON_EDGE2 == 0 && 
         cntINTSEC_ON_EDGE12              == 0    ) {
        // Penetration on the face 1 on Body 1.
        if (cntINSIDE > 0) {
            return WC_EDGE1_INSIDE_VERTICES;
        }
        else if (cntON_EDGE1 > 0) {
            return WC_EDGE1_TOUCHING;
        }
    }
    else if ( cntON_EDGE2 + cntINTSEC_ON_EDGE2 >  0 &&
              cntON_EDGE1 + cntINTSEC_ON_EDGE1 == 0 && 
              cntINTSEC_ON_EDGE12              == 0    ) {
        // Penetration on the face 1 on Body 1.
        if (cntINSIDE > 0) {
            return WC_EDGE2_INSIDE_VERTICES;
        }
        else if (cntON_EDGE2 > 0) {
            return WC_EDGE2_TOUCHING;
        }
    }
    else {
        if (cntINSIDE > 0) {
            return WC_APEX_INSIDE_VERTICES;
        }
        else if (cntINTSEC_ON_EDGE12 > 0) {
            return WC_APEX_ONE_EDGE_CUT_ACROSS;
        }
        else if (cntINTSEC_ON_APEX > 0) {
            return WC_APEX_ONE_EDGE_TOUCHING;
        }
        else if (cntON_APEX > 0) {
            return WC_APEX_COINCIDENT_VERTEX;
        }
        else {
            log(ERROR, __FILE__, __LINE__, "Forbidden predicates");
            return WC_NONE;
        }
    }
    return WC_NONE;
}


bool ContactUpdaterFurtherChecker::updateFeatures_EDGE_VERTEX(
    const enum WEDGE_CLASSIFIER    parseResult,
    vector<enum WEDGE_CLASSIFIER>& predMain,
    vector<enum WEDGE_CLASSIFIER>& predAux,
    vector<HalfEdgeIt>&            hes,
    const Vec2&                    apex,
    vector<Vec2>&                  fan,
    const Vec2&                    side1perp,
    const Vec2&                    side2perp,
    FaceIt                         fit1,
    FaceIt                         fit2,
    const bool                     wedgeOnBody1
) {
    const long numElems = fan.size();
    switch (parseResult) {

        case WC_COMPLETELY_INSIDE:
            return true;
            break;

        case WC_EDGE1_INSIDE_VERTICES:
          {

              if (wedgeOnBody1) {
                  mInfo.mType1 = ContactPairInfo::FT_FACE;
                  mInfo.mFit1  = fit1;
              }
              else {
                  mInfo.mType2 = ContactPairInfo::FT_FACE;
                  mInfo.mFit2  = fit1;
              }
              // Find the deepest inside vertex under edge 1 and then find the
              // less steep incident face.
              double maxDepth = 0.0;
              long   maxIndex = 0;
              for (long i = 0; i < numElems; i++) {
                  if (predMain[i]==WC_INSIDE) {
                      auto dot = side1perp.dot(fan[i] - apex);
                      if (dot < maxDepth) {
                          maxDepth = dot;
                          maxIndex = i;
                      }
                  }
              }
              long nextIndex = (maxIndex + 1)%numElems;
              long prevIndex = (maxIndex + numElems - 1)%numElems;
              Vec2 vPrev     = fan[maxIndex]  - fan[prevIndex];
              Vec2 vNext     = fan[nextIndex] - fan[maxIndex];
              if (wedgeOnBody1) {
                  mInfo.mType2 = ContactPairInfo::FT_FACE;
              }
              else {
                  mInfo.mType1 = ContactPairInfo::FT_FACE;
              }

              if (vPrev.squaredNorm2() <= mEpsilonZero) {
                  if (wedgeOnBody1) {
                      mInfo.mFit2 = (*(hes[prevIndex]))->face();
                  }
                  else {
                      mInfo.mFit1 = (*(hes[prevIndex]))->face();
                  }
              }
              else if (vNext.squaredNorm2() <= mEpsilonZero) {
                  if (wedgeOnBody1) {
                      mInfo.mFit2 = (*(hes[maxIndex]))->face();
                  }
                  else {
                      mInfo.mFit1 = (*(hes[maxIndex]))->face();
                  }
              }
              else {
                  vPrev.normalize();
                  vNext.normalize();

                  if (side1perp.dot(vPrev) < side1perp.dot(vNext)) {
                      if (wedgeOnBody1) {
                          mInfo.mFit2 = (*(hes[prevIndex]))->face();
                      }
                      else {
                          mInfo.mFit1 = (*(hes[prevIndex]))->face();
                      }
                  }
                  else {
                      if (wedgeOnBody1) {
                          mInfo.mFit2 = (*(hes[maxIndex]))->face();
                      }
                      else {
                          mInfo.mFit1 = (*(hes[maxIndex]))->face();
                      }
                  }
              } 
          }
          break;

        case WC_EDGE1_TOUCHING:
          {
              if (wedgeOnBody1) {
                  mInfo.mType1 = ContactPairInfo::FT_FACE;
                  mInfo.mFit1  = fit1;
              }
              else {
                  mInfo.mType2 = ContactPairInfo::FT_FACE;
                  mInfo.mFit2  = fit1;
              }


              // Find the touching vertex and check if an incident face is
              // touching
              long touchingIndex = -1;
              for (long i = 0; i < numElems; i++) {
                  if (predMain[i]==WC_ON_EDGE1) {
                      touchingIndex = i;
                      break;
                  }
              }
              long nextIndex = (touchingIndex + 1)%numElems;
              long prevIndex = (touchingIndex + numElems - 1)%numElems;
              Vec2 vPrev     = fan[prevIndex] - fan[touchingIndex];
              Vec2 vNext     = fan[nextIndex] - fan[touchingIndex];
              if (vPrev.squaredNorm2() <= mEpsilonZero) {
                  if (wedgeOnBody1) {
                      mInfo.mType2 = ContactPairInfo::FT_FACE;
                      mInfo.mFit2  = (*(hes[prevIndex]))->face();
                  }
                  else {
                      mInfo.mType1 = ContactPairInfo::FT_FACE;
                      mInfo.mFit1  = (*(hes[prevIndex]))->face();
                  }
              }
              if (vNext.squaredNorm2() <= mEpsilonZero) {
                  if (wedgeOnBody1) {
                      mInfo.mType2 = ContactPairInfo::FT_FACE;
                      mInfo.mFit2  = (*(hes[touchingIndex]))->face();
                  }
                  else {
                      mInfo.mType1 = ContactPairInfo::FT_FACE;
                      mInfo.mFit1  = (*(hes[touchingIndex]))->face();
                  }
              }
              else {
                  vPrev.normalize();
                  vNext.normalize();
                  if (side1perp.dot(vPrev) <= mEpsilonAngle) {
                      if (wedgeOnBody1) {
                          mInfo.mType2 = ContactPairInfo::FT_FACE;
                          mInfo.mFit2  = (*(hes[prevIndex]))->face();
                      }
                      else {
                          mInfo.mType1 = ContactPairInfo::FT_FACE;
                          mInfo.mFit1  = (*(hes[prevIndex]))->face();
                      }
                  }
                  else if  (side1perp.dot(vNext) <= mEpsilonAngle) {
                      if (wedgeOnBody1) {
                          mInfo.mType2 = ContactPairInfo::FT_FACE;
                          mInfo.mFit2  = (*(hes[touchingIndex]))->face();
                      }
                      else {
                          mInfo.mType1 = ContactPairInfo::FT_FACE;
                          mInfo.mFit1  = (*(hes[touchingIndex]))->face();
                      }
                  }
                  else {
                      if (wedgeOnBody1) {
                          mInfo.mType2 = ContactPairInfo::FT_EDGE;
                          mInfo.mEit2  = (*(hes[touchingIndex]))->edge();
                      }
                      else {
                          mInfo.mType1 = ContactPairInfo::FT_EDGE;
                          mInfo.mEit1  = (*(hes[touchingIndex]))->edge();
                      }
                  }
              }
          }
          break;

        case WC_EDGE2_INSIDE_VERTICES:
          {
              if (wedgeOnBody1) {
                  mInfo.mType1 = ContactPairInfo::FT_FACE;
                  mInfo.mFit1  = fit2;
              }
              else {
                  mInfo.mType2 = ContactPairInfo::FT_FACE;
                  mInfo.mFit2  = fit2;
              }

              // Find the deepest inside vertex under edge 1 and then find the
              // less steep incident face.
              double maxDepth = 0.0;
              long   maxIndex = 0;
              for (long i = 0; i < numElems; i++) {
                  if (predMain[i]==WC_INSIDE) {
                      auto dot = side2perp.dot(fan[i] - apex);
                      if (dot < maxDepth) {
                          maxDepth = dot;
                          maxIndex = i;
                      }
                  }
              }
              long nextIndex = (maxIndex + 1)%numElems;
              long prevIndex = (maxIndex + numElems - 1)%numElems;
              Vec2 vPrev     = fan[maxIndex]  - fan[prevIndex];
              Vec2 vNext     = fan[nextIndex] - fan[maxIndex];
              if (wedgeOnBody1) {
                  mInfo.mType2 = ContactPairInfo::FT_FACE;
              }
              else {
                  mInfo.mType1 = ContactPairInfo::FT_FACE;
              }

              if (vPrev.squaredNorm2() <= mEpsilonZero) {
                  if (wedgeOnBody1) {
                      mInfo.mFit2 = (*(hes[prevIndex]))->face();
                  }
                  else {
                      mInfo.mFit1 = (*(hes[prevIndex]))->face();
                  }
              }
              else if (vNext.squaredNorm2() <= mEpsilonZero) {
                  if (wedgeOnBody1) {
                      mInfo.mFit2 = (*(hes[maxIndex]))->face();
                  }
                  else {
                      mInfo.mFit1 = (*(hes[maxIndex]))->face();
                  }
              }
              else {
                  vPrev.normalize();
                  vNext.normalize();

                  if (side2perp.dot(vPrev) < side2perp.dot(vNext)) {
                      if (wedgeOnBody1) {
                          mInfo.mFit2 = (*(hes[prevIndex]))->face();
                      }
                      else {
                          mInfo.mFit1 = (*(hes[prevIndex]))->face();
                      }
                  }
                  else {
                      if (wedgeOnBody1) {
                          mInfo.mFit2 = (*(hes[maxIndex]))->face();
                      }
                      else {
                          mInfo.mFit1 = (*(hes[maxIndex]))->face();
                      }
                  }
              }
          }
          break;

        case WC_EDGE2_TOUCHING:
          {

              if (wedgeOnBody1) {
                  mInfo.mType1 = ContactPairInfo::FT_FACE;
                  mInfo.mFit1  = fit2;
              }
              else {
                  mInfo.mType2 = ContactPairInfo::FT_FACE;
                  mInfo.mFit2  = fit2;
              }
              // Find the touching vertex and check if an incident face is
              // touching
              long touchingIndex = -1;
              for (long i = 0; i < numElems; i++) {
                  if (predMain[i]==WC_ON_EDGE2) {
                      touchingIndex = i;
                      break;
                  }
              }
              long nextIndex = (touchingIndex + 1)%numElems;
              long prevIndex = (touchingIndex + numElems - 1)%numElems;
              Vec2 vPrev     = fan[prevIndex] - fan[touchingIndex];
              Vec2 vNext     = fan[nextIndex] - fan[touchingIndex];
              if (vPrev.squaredNorm2() <= mEpsilonZero) {
                  if (wedgeOnBody1) {
                      mInfo.mType2 = ContactPairInfo::FT_FACE;
                      mInfo.mFit2  = (*(hes[prevIndex]))->face();
                  }
                  else {
                      mInfo.mType1 = ContactPairInfo::FT_FACE;
                      mInfo.mFit1  = (*(hes[prevIndex]))->face();
                  }
              }
              else if  (vNext.squaredNorm2() <= mEpsilonZero) {
                  if (wedgeOnBody1) {
                      mInfo.mType2 = ContactPairInfo::FT_FACE;
                      mInfo.mFit2  = (*(hes[touchingIndex]))->face();
                  }
                  else {
                      mInfo.mType1 = ContactPairInfo::FT_FACE;
                      mInfo.mFit1  = (*(hes[touchingIndex]))->face();
                  }
              }
              else {
                  vPrev.normalize();
                  vNext.normalize();
                  if (side2perp.dot(vPrev) <= mEpsilonAngle) {
                      if (wedgeOnBody1) {
                          mInfo.mType2 = ContactPairInfo::FT_FACE;
                          mInfo.mFit2  = (*(hes[prevIndex]))->face();
                      }
                      else {
                          mInfo.mType1 = ContactPairInfo::FT_FACE;
                          mInfo.mFit1  = (*(hes[prevIndex]))->face();
                      }
                  }
                  else if  (side2perp.dot(vNext) <= mEpsilonAngle) {
                      if (wedgeOnBody1) {
                          mInfo.mType2 = ContactPairInfo::FT_FACE;
                          mInfo.mFit2  = (*(hes[touchingIndex]))->face();
                      }
                      else {
                          mInfo.mType1 = ContactPairInfo::FT_FACE;
                          mInfo.mFit1  = (*(hes[touchingIndex]))->face();
                      }
                  }
                  else {
                      if (wedgeOnBody1) {
                          mInfo.mType2 = ContactPairInfo::FT_EDGE;
                          mInfo.mEit2  = (*(hes[touchingIndex]))->edge();
                      }
                      else {
                          mInfo.mType1 = ContactPairInfo::FT_EDGE;
                          mInfo.mEit1  = (*(hes[touchingIndex]))->edge();
                      }
                  }
              }
          }
          break;
    
        case WC_APEX_ONE_EDGE_CUT_ACROSS:
          {
              for (long i = 0; i < numElems; i++) {
                  if (predAux[i]==WC_INTSEC_ON_EDGE12) {
                      if (wedgeOnBody1) {
                          mInfo.mType2 = ContactPairInfo::FT_FACE;
                          mInfo.mFit2  = (*(hes[i]))->face();
                      }
                      else {
                          mInfo.mType1 = ContactPairInfo::FT_FACE;
                          mInfo.mFit1  = (*(hes[i]))->face();
                      }
                      break;
                  }
              }
          }
          break;

        case WC_APEX_ONE_EDGE_TOUCHING:
          {
              for (long i = 0; i < numElems; i++) {
                  if (predAux[i]==WC_INTSEC_ON_APEX) {

                      if (wedgeOnBody1) {
                          mInfo.mType2 = ContactPairInfo::FT_FACE;
                          mInfo.mFit2  = (*(hes[i]))->face();
                      }
                      else {
                          mInfo.mType1 = ContactPairInfo::FT_FACE;
                          mInfo.mFit1  = (*(hes[i]))->face();
                      }
                      break;
                  }
              }
          }
          break;

        case WC_APEX_INSIDE_VERTICES:
          {
              // Find the deepest inside vertex from the apex.
              // Check the incident faces. If both go up toward apex,
              // make it EDGE-EDGE (parallel). If one of them goes down.
              // makd it FACE-EDGE.
              double maxDepth = apex.y();
              long   maxIndex = 0;
              for (long i = 0; i < numElems; i++) {
                  if (predMain[i]==WC_INSIDE) {
                      if (fan[i].y() < maxDepth) {
                          maxDepth = fan[i].y();
                          maxIndex = i;
                      }
                  }
              }
              long nextIndex = (maxIndex + 1)%numElems;
              long prevIndex = (maxIndex + numElems - 1)%numElems;
              if (fan[maxIndex].y() > fan[prevIndex].y()) {

                  if (fan[maxIndex].x() > fan[prevIndex].x()) {
                      if (wedgeOnBody1) {
                          mInfo.mType1 = ContactPairInfo::FT_FACE;
                          mInfo.mFit1  = fit2;
                      }
                      else {
                          mInfo.mType2 = ContactPairInfo::FT_FACE;
                          mInfo.mFit2  = fit2;
                      }
                  }
                  else {
                      if (wedgeOnBody1) {
                          mInfo.mType1 = ContactPairInfo::FT_FACE;
                          mInfo.mFit1  = fit1;
                      }
                      else {
                          mInfo.mType2 = ContactPairInfo::FT_FACE;
                          mInfo.mFit2  = fit1;
                      }
                  }
                  if (wedgeOnBody1) {
                      mInfo.mType2 = ContactPairInfo::FT_FACE;
                      mInfo.mFit2  = (*(hes[prevIndex]))->face();
                  }
                  else {
                      mInfo.mType1 = ContactPairInfo::FT_FACE;
                      mInfo.mFit1  = (*(hes[prevIndex]))->face();
                  }
              }
              else if  (fan[maxIndex].y() > fan[nextIndex].y()) {

                  if (fan[maxIndex].x() > fan[nextIndex].x()) {
                      if (wedgeOnBody1) {
                          mInfo.mType1 = ContactPairInfo::FT_FACE;
                          mInfo.mFit1  = fit2;
                      }
                      else {
                          mInfo.mType2 = ContactPairInfo::FT_FACE;
                          mInfo.mFit2  = fit2;
                      }
                  }
                  else {
                      if (wedgeOnBody1) {
                          mInfo.mType1 = ContactPairInfo::FT_FACE;
                          mInfo.mFit1  = fit1;
                      }
                      else {
                          mInfo.mType2 = ContactPairInfo::FT_FACE;
                          mInfo.mFit2  = fit1;
                      }
                  }
                  if (wedgeOnBody1) {
                      mInfo.mType2 = ContactPairInfo::FT_FACE;
                      mInfo.mFit2  = (*(hes[maxIndex]))->face();
                  }
                  else {
                      mInfo.mType1 = ContactPairInfo::FT_FACE;
                      mInfo.mFit1  = (*(hes[maxIndex]))->face();
                  }
              }
              else {
                  if (wedgeOnBody1) {
                      mInfo.mType2 = ContactPairInfo::FT_EDGE;
                      mInfo.mEit2  = (*(hes[maxIndex]))->edge();
                  }
                  else {
                      mInfo.mType1 = ContactPairInfo::FT_EDGE;
                      mInfo.mEit1  = (*(hes[maxIndex]))->edge();
                  }
              }
          }
          break;

        case WC_APEX_COINCIDENT_VERTEX:
          {
            for (long i = 0; i < numElems; i++) {
                if (predMain[i]==WC_ON_APEX) {
                    if (wedgeOnBody1) {
                        mInfo.mType2 = ContactPairInfo::FT_EDGE;
                        mInfo.mEit2  = (*(hes[i]))->edge();
                    }
                    else {
                        mInfo.mType1 = ContactPairInfo::FT_EDGE;
                        mInfo.mEit1  = (*(hes[i]))->edge();
                    }
                    break;
                }
            }
          }
          break;

        default:
          break;
    }
    return false;
}


bool ContactUpdaterFurtherChecker::checkFeatures_VERTEX_EDGE()
{
    auto body2_he1  = (*(mInfo.mEit2))->he1();
    auto body2_he2  = (*(mInfo.mEit2))->he2();
    auto body2_vit1 = (*body2_he1)->src();
    auto body2_vit2 = (*body2_he1)->dst();
    auto body2_p1   = (*body2_vit1)->pGCS(mBody2.Qmat(), mBody2.CoM());
    auto body2_p2   = (*body2_vit2)->pGCS(mBody2.Qmat(), mBody2.CoM());
    auto body2_v12  = body2_p2 - body2_p1;
    body2_v12.normalize();
    auto body2_n    = (*mInfo.mEit2)->nGCS(mBody2.Qmat());

    // Rotate the bodies such that the edge normal (body1_n) is aligned with
    // Y-axis and the edge direction (from he1's src to dst, body1_v12) with 
    // Z-axis.
    Mat3x3 rotMat1, rotMat2;
    Vec3   com1,    com2;
    rotateWorld(body2_n, body2_v12, rotMat1, com1, rotMat2, com2);

    body2_p1        = (*body2_vit1)->pGCS(rotMat2, com2);
    body2_p2        = (*body2_vit2)->pGCS(rotMat2, com2);
    body2_v12       = body2_p2 - body2_p1;
    body2_v12.normalize();
    auto body2_fit1 = (*body2_he1)->face();
    auto body2_fit2 = (*body2_he2)->face();
    auto body2_fn1  = (*body2_fit1)->nGCS(rotMat2);
    auto body2_fn2  = (*body2_fit2)->nGCS(rotMat2);
    auto body2_dn1  = body2_fn1.cross(body2_v12);
    auto body2_dn2  = body2_v12.cross(body2_fn2);

    // Project the features onto the new XY-plane.
    vector<Vec2>       body1_fan;
    vector<HalfEdgeIt> body1_hes;
    auto body1_p   = (*mInfo.mVit1)->pGCS(rotMat1, com1);
    for (auto body1_he : (*mInfo.mVit1)->halfEdges()) {
        if ((*body1_he)->src()==mInfo.mVit1) {
            auto body1_vitFan = (*body1_he)->dst();
            auto body1_vFan   = (*body1_vitFan)->pGCS(rotMat1, com1) -
                                body1_p;
            body1_vFan.normalize();
            Vec2 p_2D(body1_p.x() + body1_vFan.x(), body1_p.y() + body1_vFan.y());
            body1_fan.push_back(p_2D);
            body1_hes.push_back(body1_he);
        }            
    }

    Vec2 apex           (body2_p1.x(),         body2_p1.y());
    Vec2 body2_side1    (body2_dn1.x(),        body2_dn1.y());
    body2_side1.normalize();
    Vec2 body2_side1perp(-1.0 * body2_dn1.y(), body2_dn1.x());
    Vec2 body2_side2    (body2_dn2.x(),        body2_dn2.y());
    body2_side2.normalize();
    Vec2 body2_side2perp(body2_dn2.y(),        -1.0*body2_dn2.x());

    // Classify the vertex w.r.t the wedge.
    vector<enum WEDGE_CLASSIFIER> predMain;
    vector<enum WEDGE_CLASSIFIER> predAux;

    bool intersecting = classifyPointsAgainstWedge(
                            body1_fan,
                            apex,
                            body2_side1,
                            body2_side1perp,
                            body2_side2,
                            body2_side2perp,
                            predMain,
                            predAux
                        );

    if (!intersecting) {
        return false;
    }

    auto parseResult = parseClassifiersPointsAgainstWedge(predMain, predAux);
    return updateFeatures_EDGE_VERTEX(
        parseResult,
        predMain,
        predAux,
        body1_hes,
        apex,
        body1_fan,
        body2_side1perp,
        body2_side2perp,
        body2_fit1,
        body2_fit2,
        false
    );

}


bool ContactUpdaterFurtherChecker::checkFeatures_VERTEX_VERTEX()
{
    vector<Vec3>       body1_interior, body2_interior;

    auto body1_vitApex = mInfo.mVit1;
    auto body1_pApex   = (*body1_vitApex)->pGCS(mBody1.Qmat(), mBody1.CoM());
    for (auto hePillar : (*body1_vitApex)->halfEdges()) {
        if ((*hePillar)->src()==body1_vitApex) {
            auto vit = (*hePillar)->dst();
            auto p   = (*vit)->pGCS(mBody1.Qmat(), mBody1.CoM());
            auto v   = p - body1_pApex;
            v.normalize();
            if (isRayInsideCone(v, mInfo.mVit2, false)) {
                body1_interior.push_back(v);
            }
        }
    }

    auto body2_vitApex = mInfo.mVit2;
    auto body2_pApex   = (*body2_vitApex)->pGCS(mBody2.Qmat(), mBody2.CoM());
    for (auto hePillar : (*body2_vitApex)->halfEdges()) {
        if ((*hePillar)->src()==body2_vitApex) {
            auto vit = (*hePillar)->dst();
            auto p   = (*vit)->pGCS(mBody2.Qmat(), mBody2.CoM());
            auto v   = p - body2_pApex;
            v.normalize();
            if (isRayInsideCone(v, mInfo.mVit1, true)) {
                body2_interior.push_back(v);
            }
        }
    }

    if (body1_interior.size()==0 && body2_interior.size()==0) {
        return false;
    }

    if ( (body1_interior.size()==((*body1_vitApex)->halfEdges()).size()/2)||
         (body2_interior.size()==((*body2_vitApex)->halfEdges()).size()/2)   ){
        return true; // One cone is completely inside the other.
    }

    Vec3 dAvg(0.0, 0.0, 0.0);
    for (auto& v : body1_interior) {
       dAvg += v;
    }
    for (auto& v : body2_interior) {
       dAvg += v;
    }
    dAvg.normalize();
    Mat3x3 rotMat1, rotMat2;
    Vec3   com1,    com2;
    rotateWorld(dAvg, rotMat1, com1, rotMat2, com2);

    vector<IntersectionFinderConvexPolygon2D::InputElem> fan1_2D, fan2_2D;
    vector<IntersectionFinderConvexPolygon2D::OutputElem> intsec_2D;
    vector<HalfEdgeIt> body1_hes, body2_hes;

    body1_pApex   = (*body1_vitApex)->pGCS(rotMat1, com1);
    long index = 0;
    for (auto hePillar : (*body1_vitApex)->halfEdges()) {
        if ((*hePillar)->src()==body1_vitApex) {
            auto vit = (*hePillar)->dst();
            auto v   = (*vit)->pGCS(rotMat1, com1) - body1_pApex;
            v.normalize();
            auto p   = body1_pApex + v;
            IntersectionFinderConvexPolygon2D::InputElem 
                                                e(Vec2(p.x(), p.y()), index++);
            fan1_2D.push_back(e);
            body1_hes.push_back(hePillar);
        }
    }
    long nElems1 = index;
    body2_pApex   = (*body2_vitApex)->pGCS(rotMat2, com2);
    index = 0;
    for (auto hePillar : (*body2_vitApex)->halfEdges()) {
        if ((*hePillar)->src()==body2_vitApex) {
            auto vit = (*hePillar)->dst();
            auto v   = (*vit)->pGCS(rotMat2, com2) - body2_pApex;
            v.normalize();
            auto p   = body2_pApex + v;
            IntersectionFinderConvexPolygon2D::InputElem 
                                                e(Vec2(p.x(), p.y()), index++);
            fan2_2D.push_back(e);
            body2_hes.push_back(hePillar);
        }
    }
    long nElems2 = index;

    vector<IntersectionFinderConvexPolygon2D::InputElem> 
        fan1_2D_rev(fan1_2D.rbegin(), fan1_2D.rend()),
        fan2_2D_rev(fan2_2D.rbegin(), fan2_2D.rend());

    IntersectionFinderConvexPolygon2D 
                 finder(mEpsilonZero, mEpsilonZero, mEpsilonAngle, mLogStream);
    finder.findIntersection(fan1_2D_rev, fan2_2D_rev, intsec_2D);
    if (intsec_2D.size()==0) {
        // No penetration.
        return false;
    }
    else if  (intsec_2D.size()==1) {

        // Touching at one point. edge-point or point-point.
        auto ia1 = intsec_2D[0].mIndexA;
        auto ia2 = intsec_2D[0].mIndexAaux;
        auto ib1 = intsec_2D[0].mIndexB;
        auto ib2 = intsec_2D[0].mIndexBaux;

        switch (intsec_2D[0].mType) {

          case IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX:
            mInfo.mType1 = ContactPairInfo::FT_FACE;
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            if (ia1==0||ia2==0) {
                if (ia1==1||ia2==1) {
                    mInfo.mFit1 = (*(body1_hes[0]))->face();
                }
                else {                
                    mInfo.mFit1 = (*(body1_hes[nElems1-1]))->face();
                }
            }
            else {
                mInfo.mFit1 = (*(body1_hes[std::min(ia1, ia2)]))->face();
            }           
            mInfo.mEit2 = (*(body2_hes[ib1]))->edge();
            break;

          case IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX:
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1 = (*(body1_hes[ia1]))->edge();
            mInfo.mEit2 = (*(body2_hes[ib1]))->edge();
            break;

          case IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE:
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mType2 = ContactPairInfo::FT_FACE;
            mInfo.mEit1 = (*(body1_hes[ia1]))->edge();
            if (ib1==0||ib2==0) {
                if (ib1==1||ib2==1) {
                    mInfo.mFit2 = (*(body2_hes[0]))->face();
                }
                else {                
                    mInfo.mFit2 = (*(body2_hes[nElems2-1]))->face();
                }
            }
            else {
                mInfo.mFit2 = (*(body2_hes[std::min(ib1, ib2)]))->face();
            }           
            break;

          case IntersectionFinderConvexPolygon2D::IT_EDGE_EDGE:
            mInfo.mType1 = ContactPairInfo::FT_FACE;
            mInfo.mType2 = ContactPairInfo::FT_FACE;
            if (ia1==0||ia2==0) {
                if (ia1==1||ia2==1) {
                    mInfo.mFit1 = (*(body1_hes[0]))->face();
                }
                else {                
                    mInfo.mFit1 = (*(body1_hes[nElems1-1]))->face();
                }
            }
            else {
                mInfo.mFit1 = (*(body1_hes[std::min(ia1, ia2)]))->face();
            }           
            if (ib1==0||ib2==0) {
                if (ib1==1||ib2==1) {
                    mInfo.mFit2 = (*(body2_hes[0]))->face();
                }
                else {                
                    mInfo.mFit2 = (*(body2_hes[nElems2-1]))->face();
                }
            }
            else {
                mInfo.mFit2 = (*(body2_hes[std::min(ib1, ib2)]))->face();
            }           
            break;

          default:
            log(ERROR, __FILE__, __LINE__, "Forbidden predicates");
            break;
        }
    }
    else if  (intsec_2D.size()==2) {

        // Touching in edge-edge.
        bool edgeOnBody1Found = false;
        bool edgeOnBody2Found = false;

        auto ia1 = intsec_2D[0].mIndexA;
        auto ia2 = intsec_2D[0].mIndexAaux;
        auto ib1 = intsec_2D[0].mIndexB;
        auto ib2 = intsec_2D[0].mIndexBaux;
        long body1Vid1 = ia1;
        long body2Vid1 = ib1;

        switch (intsec_2D[0].mType) {
          case IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX:
            mInfo.mType1 = ContactPairInfo::FT_FACE;
            if (ia1==0||ia2==0) {
                if (ia1==1||ia2==1) {
                    mInfo.mFit1 = (*(body1_hes[0]))->face();
                }
                else {                
                    mInfo.mFit1 = (*(body1_hes[nElems1-1]))->face();
                }
            }
            else {
                mInfo.mFit1 = (*(body1_hes[std::min(ia1, ia2)]))->face();
            }           
            edgeOnBody1Found = true;
            break;

          case IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX:
            break;

          case IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE:
            mInfo.mType2 = ContactPairInfo::FT_FACE;
            if (ib1==0||ib2==0) {
                if (ib1==1||ib2==1) {
                    mInfo.mFit2 = (*(body2_hes[0]))->face();
                }
                else {                
                    mInfo.mFit2 = (*(body2_hes[nElems2-1]))->face();
                }
            }
            else {
                mInfo.mFit2 = (*(body2_hes[std::min(ib1, ib2)]))->face();
            }           
            edgeOnBody2Found = true;
            break;

          default:
            log(ERROR, __FILE__, __LINE__, "Forbidden predicates");
            break;
        }

        ia1 = intsec_2D[1].mIndexA;
        ia2 = intsec_2D[1].mIndexAaux;
        ib1 = intsec_2D[1].mIndexB;
        ib2 = intsec_2D[1].mIndexBaux;
        long body1Vid2 = ia1;
        long body2Vid2 = ib1;

        switch (intsec_2D[1].mType) {
          case IntersectionFinderConvexPolygon2D::IT_EDGE_VERTEX:
            if (!edgeOnBody1Found) {
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                if (ia1==0||ia2==0) {
                    if (ia1==1||ia2==1) {
                        mInfo.mFit1 = (*(body1_hes[0]))->face();
                    }
                    else {                
                        mInfo.mFit1 = (*(body1_hes[nElems1-1]))->face();
                    }
                }
                else {
                    mInfo.mFit1 = (*(body1_hes[std::min(ia1, ia2)]))->face();
                }           
                edgeOnBody1Found = true;
            }
            break;

          case IntersectionFinderConvexPolygon2D::IT_VERTEX_VERTEX:
            break;

          case IntersectionFinderConvexPolygon2D::IT_VERTEX_EDGE:
            if (!edgeOnBody2Found) {
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                if (ib1==0||ib2==0) {
                    if (ib1==1||ib2==1) {
                        mInfo.mFit2 = (*(body2_hes[0]))->face();
                    }
                    else {                
                        mInfo.mFit2 = (*(body2_hes[nElems2-1]))->face();
                    }
                }
                else {
                    mInfo.mFit2 = (*(body2_hes[std::min(ib1, ib2)]))->face();
                }           
                edgeOnBody2Found = true;
            }
            break;

          default:
            log(ERROR, __FILE__, __LINE__, "Forbidden predicates");
            break;
        }

        if (!edgeOnBody1Found) {

            if ( (body1Vid1==0 && body1Vid2==1) ||
                 (body1Vid1==1 && body1Vid2==0)   ) {
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = (*(body1_hes[0]))->face();
            }
            else if ( (body1Vid1==0           && body1Vid2==(nElems1-1)) ||
                      (body1Vid1==(nElems1-1) && body1Vid2==0              ) ){
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = (*(body1_hes[nElems1-1]))->face();
            }
            else if (body1Vid1+1 == body1Vid2) {
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = (*(body1_hes[body1Vid1]))->face();
            }
            else  if (body1Vid2+1 == body1Vid1) {
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = (*(body1_hes[body1Vid2]))->face();

            }
            else {
                log(ERROR, __FILE__, __LINE__, "Inconsistent vertex IDs");
            }
        }

        if (!edgeOnBody2Found) {

            if ( (body2Vid1==0 && body2Vid2==1) ||
                 (body2Vid1==1 && body2Vid2==0)   ) {
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit2  = (*(body2_hes[0]))->face();
            }
            else if ( (body2Vid1==0           && body2Vid2==(nElems2-1)) ||
                      (body2Vid1==(nElems2-1) && body2Vid2==0              ) ){
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit2  = (*(body2_hes[nElems2-1]))->face();
            }
            else if (body2Vid1+1 == body2Vid2) {
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit2  = (*(body2_hes[body2Vid1]))->face();
            }
            else  if (body2Vid2+1 == body2Vid1) {
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit2  = (*(body2_hes[body2Vid2]))->face();
            }
            else {
                log(ERROR, __FILE__, __LINE__, "Inconsistent vertex IDs");
            }
        }
    }
    else {
        Vec2 vPenDir;
        auto found = 
                guessPenetrationDirectionFromIntersection(intsec_2D, vPenDir);

        if (!found) {
            guessPenetrationDirectionFromSpread(
                                         fan1_2D, fan2_2D, intsec_2D,vPenDir);
        }

        double extremeDot1,   extremeDot2;
        long   extremeIndex1, extremeIndex2;
        for (long i = 0 ; i < intsec_2D.size(); i++) {
            auto& e = intsec_2D[i];
            if (i==0) {
                extremeDot1   = vPenDir.dot(e.mP);
                extremeDot2   = vPenDir.dot(e.mP);
                extremeIndex1 = e.mIndexA;
                extremeIndex2 = e.mIndexB; 
            }
            else {
                auto dot = vPenDir.dot(e.mP);
                if (extremeDot1 < dot) {
                    extremeDot1   = dot;
                    extremeIndex1 = e.mIndexA;
                }
                if (extremeDot2 > dot) {
                    extremeDot2   = dot;
                    extremeIndex2 = e.mIndexB;
                }
            }
        }
        findAppropriateFacePair(extremeIndex1, extremeIndex2,
                                      fan1_2D, fan2_2D, body1_hes, body2_hes);
    }
    return false;
}

void ContactUpdaterFurtherChecker::findAppropriateFacePair(
    long                                                  indexBody1,
    long                                                  indexBody2,
    vector<IntersectionFinderConvexPolygon2D::InputElem>& fan1_2D,
    vector<IntersectionFinderConvexPolygon2D::InputElem>& fan2_2D,
    vector<HalfEdgeIt>&                                   body1_hes,
    vector<HalfEdgeIt>&                                   body2_hes
) {

    long numElems1      = fan1_2D.size();
    long numElems2      = fan2_2D.size();
    long indexBody1Prev = (indexBody1 + numElems1 - 1)%numElems1;
    long indexBody1Cur  = indexBody1;
    long indexBody1Next = (indexBody1 + 1)%numElems1;
    long indexBody2Prev = (indexBody2 + numElems2 - 1)%numElems2;
    long indexBody2Cur  = indexBody2;
    long indexBody2Next = (indexBody2 + 1)%numElems2;

    auto& pBody1Prev = fan1_2D[indexBody1Prev].mP;
    auto& pBody1Cur  = fan1_2D[indexBody1Cur].mP;
    auto& pBody1Next = fan1_2D[indexBody1Next].mP;

    auto& pBody2Prev = fan2_2D[indexBody2Prev].mP;
    auto& pBody2Cur  = fan2_2D[indexBody2Cur].mP;
    auto& pBody2Next = fan2_2D[indexBody2Next].mP;

    auto  he11 = body1_hes[indexBody1Prev];
    auto  he12 = body1_hes[indexBody1Cur];
    auto  he21 = body2_hes[indexBody2Prev];
    auto  he22 = body2_hes[indexBody2Cur];

    auto  fit11 = (*he11)->face();
    auto  fit12 = (*he12)->face();
    auto  fit21 = (*he21)->face();
    auto  fit22 = (*he22)->face();
    auto  fn11  = (*fit11)->nGCS(mBody1.Qmat());
    auto  fn12  = (*fit12)->nGCS(mBody1.Qmat());
    auto  fn21  = (*fit21)->nGCS(mBody2.Qmat());
    auto  fn22  = (*fit22)->nGCS(mBody2.Qmat());

    double area11(0.0),     area12(0.0),     area21(0.0),     area22(0.0);
    double distance11(0.0), distance12(0.0), distance21(0.0), distance22(0.0);
    double alignment11(0.0),alignment12(0.0),alignment21(0.0),alignment22(0.0);

    if (fn11.dot(fn21) <= -1.0*mEpsilonAngle) {
        findProjectedAreaAndDistance(pBody1Prev, pBody1Cur, 
                                     pBody2Prev, pBody2Cur,
                                     area11, distance11, alignment11);
    }
    if (fn11.dot(fn22) <= -1.0*mEpsilonAngle) {
        findProjectedAreaAndDistance(pBody1Prev, pBody1Cur, 
                                     pBody2Cur, pBody2Next,
                                     area12, distance12, alignment12);
    }
    if (fn12.dot(fn21) <= -1.0*mEpsilonAngle) {
        findProjectedAreaAndDistance(pBody1Cur, pBody1Next, 
                                     pBody2Prev, pBody2Cur,
                                     area21, distance21, alignment21);
    }
    if (fn12.dot(fn22) <= -1.0*mEpsilonAngle) {
        findProjectedAreaAndDistance(pBody1Cur, pBody1Next, 
                                     pBody2Cur, pBody2Next,
                                     area22, distance22, alignment22);
    }
    
    double score11 = area11 / ( (mEpsilonZero*10.0 + distance11 ) * 
                                (mEpsilonZero*10.0 + alignment11)   );
    double score12 = area12 / ( (mEpsilonZero*10.0 + distance12 ) * 
                                (mEpsilonZero*10.0 + alignment12)   );
    double score21 = area21 / ( (mEpsilonZero*10.0 + distance21 ) * 
                                (mEpsilonZero*10.0 + alignment21)   );
    double score22 = area22 / ( (mEpsilonZero*10.0 + distance22 ) * 
                                (mEpsilonZero*10.0 + alignment22)   );

    if (score11 == 0.0 && score12 == 0.0 && score21 == 0.0 && score22 == 0.0) {
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*(body1_hes[indexBody1]))->edge();
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = (*(body2_hes[indexBody2]))->edge();
    }
    else if (score11 >= score12 && score11 >= score21 && score11 >= score22) {
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1  = (*(body1_hes[indexBody1Prev]))->face();
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2  = (*(body2_hes[indexBody2Prev]))->face();
    }
    else if (score12 >= score11 && score12 >= score21 && score12 >= score22) {
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1  = (*(body1_hes[indexBody1Prev]))->face();
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2  = (*(body2_hes[indexBody2Cur]))->face();
    }
    else if (score21 >= score11 && score21 >= score12 && score21 >= score22) {
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1  = (*(body1_hes[indexBody1Cur]))->face();
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2  = (*(body2_hes[indexBody2Prev]))->face();
    }
    else {
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1  = (*(body1_hes[indexBody1Cur]))->face();
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2  = (*(body2_hes[indexBody2Cur]))->face();
    }
}


void ContactUpdaterFurtherChecker::guessPenetrationDirectionFromSpread(
    vector<IntersectionFinderConvexPolygon2D::InputElem>&  fan1_2D,
    vector<IntersectionFinderConvexPolygon2D::InputElem>&  fan2_2D,
    vector<IntersectionFinderConvexPolygon2D::OutputElem>& intsec_2D,
    Vec2&                                                 dir
) {
    vector<Vec2> fanPoints_2D, intsecPoints_2D;

    Vec2 mean1(0.0, 0.0);
    for (auto& e : fan1_2D) {
        fanPoints_2D.push_back(e.mP);
        mean1 += e.mP;
    }
    mean1.scale(1.0/fan1_2D.size());
    Vec2 mean2(0.0, 0.0);

    for (auto& e : fan2_2D) {
        fanPoints_2D.push_back(e.mP);
        mean2 += e.mP;
    }
    mean2.scale(1.0/fan2_2D.size());

    vector<long> fanPoints_2Dind = findConvexHull2D(fanPoints_2D);

    vector<Vec2> fanPoints_2DCH;
    for (auto i : fanPoints_2Dind) {
        fanPoints_2DCH.push_back(fanPoints_2D[i]);
    }

    for (auto& e : intsec_2D) {
        intsecPoints_2D.push_back(e.mP);
    }

    Vec2 spreadFan, spreadIntsec;
    Vec2 meanFan,   meanIntsec;
    Vec2 axis1Fan,  axis1Intsec, axis2Fan,  axis2Intsec;

    findPrincipalComponents(fanPoints_2DCH, 
                           spreadFan,    meanFan,    axis1Fan,    axis2Fan);
                    
    findPrincipalComponents(intsecPoints_2D,
                           spreadIntsec, meanIntsec, axis1Intsec, axis2Intsec);

    // Guess the penetration direction based on the spreads of 
    // the fans and the intersection.
    Vec2 dir1to2 = mean2 - mean1;
    dir1to2.normalize();
    if (dir1to2.dot(axis1Fan) < 0.0) {
        axis1Fan.scale(-1.0);
    }        
    if (dir1to2.dot(axis2Intsec) < 0.0) {
        axis2Intsec.scale(-1.0);
    }        

    if (axis1Fan.dot(axis2Intsec) <= COSINE_75_DEGREE) {
        dir = axis1Fan + axis2Intsec;
        dir.normalize();
    }
    else {
        // If the intersection's primary axis deviates too much from
        // fan's, then take the fan's.
        dir = axis1Fan;
    }
}


bool ContactUpdaterFurtherChecker::guessPenetrationDirectionFromIntersection(
    vector<IntersectionFinderConvexPolygon2D::OutputElem>& intsec,
    Vec2&                                                  dir
) {
    vector<long> vi;
    vector<long> iv;
    vector<long> others;
    for (long i = 0; i < intsec.size(); i++) {
        if (intsec[i].mType ==
                IntersectionFinderConvexPolygon2D::IT_VERTEX_INTERIOR) {
            vi.push_back(i);
        }
        else if  (intsec[i].mType ==
                IntersectionFinderConvexPolygon2D::IT_INTERIOR_VERTEX) {
            iv.push_back(i);
        }
        else {
            others.push_back(i);
        }
    }

    Vec2 body1_avg(0.0, 0.0);
    Vec2 body2_avg(0.0, 0.0);
    if (vi.size()>0) {
        if (iv.size()>0) {
            for (auto i : vi) {
                body2_avg += intsec[i].mP;
            }
            body2_avg.scale(1.0/vi.size());
            for (auto i : iv) {
                body1_avg += intsec[i].mP;
            }
            body1_avg.scale(1.0/iv.size());

        }
        else {
            for (auto i : vi) {
                body2_avg += intsec[i].mP;
            }
            body2_avg.scale(1.0/vi.size());
            for (auto i : others) {
                body1_avg += intsec[i].mP;
            }
            body1_avg.scale(1.0/others.size());
        }
    }
    else {
        if (iv.size()>0) {
            for (auto i : iv) {
                body1_avg += intsec[i].mP;
            }
            body1_avg.scale(1.0/iv.size());
            for (auto i : others) {
                body2_avg += intsec[i].mP;
            }
            body2_avg.scale(1.0/others.size());
        }
        else {
            return false;
        }
    }
    dir = body2_avg - body1_avg;
    return true;
}


void ContactUpdaterFurtherChecker::findProjectedAreaAndDistance(
    const Vec2& p11,
    const Vec2& p12,
    const Vec2& p21,
    const Vec2& p22,
    double&     area,
    double&     distance,
    double&     alignment
) {
    Vec2 v1 = p12 - p11;
    double v1Len = v1.norm2();
    v1.scale(1.0/v1Len);
    Vec2 v1perp(-1.0*v1.y(), v1.x());

    Vec2 v2 = p22 - p21;
    v2.normalize();

    alignment = fabs(v1perp.dot(v2));

    Vec2 vTest1 = p21 - p11;
    Vec2 vTest2 = p21 - p12;
    Vec2 vTest3 = p22 - p11;
    Vec2 vTest4 = p22 - p12;

    auto test1 = vTest1.dot(v1);
    auto test2 = vTest2.dot(v1);
    auto test3 = vTest3.dot(v1);
    auto test4 = vTest4.dot(v1);

    if (test1 >= 0.0 && test2 <= 0.0) {

        // P21's projection is inside (p11, p12).

        if (test3 >= 0.0 && test4 <= 0.0) {
            // P22's projection is inside (p11, p12).
            area = fabs(test1 - test3);
            auto dot1 = vTest1.dot(v1perp);
            auto dot2 = vTest3.dot(v1perp);
            if (dot1 * dot2 <= 0.0) {
                distance = 0.0;
            }              
            else {
                distance = std::min(fabs(dot1), fabs(dot2));
            }
        }
        else if (test3 <= 0.0) {
            // P22's projection is outside of p11
            area = test1;
            Vec2 intsec;
            findIntersection(p11, p11 + v1perp, p21, p22, intsec);
            auto vTest5 = intsec - p11;
            auto dot1 = vTest5.dot(v1perp);
            auto dot2 = vTest1.dot(v1perp);
            if (dot1 * dot2 <= 0.0) {
                distance = 0.0;
            }              
            else {
                distance = std::min(fabs(dot1), fabs(dot2));
            }
        }
        else {
            // P22's projection is outside of p12
            area = -1.0 * test2;
            Vec2 intsec;
            findIntersection(p12, p12 + v1perp, p21, p22, intsec);
            auto vTest5 = intsec - p12;
            auto dot1 = vTest5.dot(v1perp);
            auto dot2 = vTest2.dot(v1perp);
            if (dot1 * dot2 <= 0.0) {
                distance = 0.0;
            }              
            else {
                distance = std::min(fabs(dot1), fabs(dot2));
            }
        }

    }
    else if (test1 <= 0.0) {

        // P21's projection is outside p11.

        if (test3 >= 0.0 && test4 <= 0.0) {
            // P22's projection is inside (p11, p12).
            area = test3;
            Vec2 intsec;
            findIntersection(p11, p11 + v1perp, p21, p22, intsec);
            auto vTest5 = intsec - p11;
            auto dot1 = vTest5.dot(v1perp);
            auto dot2 = vTest3.dot(v1perp);
            if (dot1 * dot2 <= 0.0) {
                distance = 0.0;
            }              
            else {
                distance = std::min(fabs(dot1), fabs(dot2));
            }
        }
        else if (test3 <= 0.0) {
            // P22's projection is outside of p11
            area = 0.0;
            distance = 0.0;
        }
        else {
            // P22's projection is outside of p12
            area = v1Len;
            Vec2 intsec1, intsec2;
            findIntersection(p11, p11 + v1perp, p21, p22, intsec1);
            findIntersection(p12, p12 + v1perp, p21, p22, intsec2);
            auto vTest5 = intsec1 - p11;
            auto vTest6 = intsec2 - p12;
            auto dot1 = vTest5.dot(v1perp);
            auto dot2 = vTest6.dot(v1perp);
            if (dot1 * dot2 <= 0.0) {
                distance = 0.0;
            }              
            else {
                distance = std::min(fabs(dot1), fabs(dot2));
            }
        }

    }
    else {

        // P21's projection is outside p12.

        if (test3 >= 0.0 && test4 <= 0.0) {
            // P22's projection is inside (p11, p12).
            area = -1.0 * test4;
            Vec2 intsec;
            findIntersection(p12, p12 + v1perp, p21, p22, intsec);
            auto vTest5 = intsec - p12;
            auto dot1 = vTest5.dot(v1perp);
            auto dot2 = vTest4.dot(v1perp);
            if (dot1 * dot2 <= 0.0) {
                distance = 0.0;
            }              
            else {
                distance = std::min(fabs(dot1), fabs(dot2));
            }
        }
        else if (test3 <= 0.0) {
            // P22's projection is outside of p11
            area = v1Len;
            Vec2 intsec1, intsec2;
            findIntersection(p11, p11 + v1perp, p21, p22, intsec1);
            findIntersection(p12, p12 + v1perp, p21, p22, intsec2);
            auto vTest5 = intsec1 - p11;
            auto vTest6 = intsec2 - p12;
            auto dot1 = vTest5.dot(v1perp);
            auto dot2 = vTest6.dot(v1perp);
            if (dot1 * dot2 <= 0.0) {
                distance = 0.0;
            }              
            else {
                distance = std::min(fabs(dot1), fabs(dot2));
            }
        }
        else {
            // P22's projection is outside of p12
            area = 0.0;
            distance = 0.0;
        }
    }
}


}// namespace Makena

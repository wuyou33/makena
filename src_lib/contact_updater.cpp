#include "contact_updater.hpp"
#include "contact_updater_face_face.hpp"
#include "contact_updater_face_edge.hpp"
#include "contact_updater_face_vertex.hpp"
#include "contact_updater_edge_edge.hpp"
#include "contact_updater_edge_vertex.hpp"
#include "contact_updater_vertex_vertex.hpp"
#include "contact_updater_further_checker.hpp"
#include "contact_points_and_normal_generator.hpp"

#include "gjk_origin_finder.hpp"
#include "bd_boundary_simplex_finder.hpp"
#include "contact_points_and_normal_generator.hpp"

/**
 * @file contact_updater.cpp
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
ContactUpdater::ContactUpdater(
    ConvexRigidBody&  body1,
    ConvexRigidBody&  body2,
    ContactPairInfo&  info,
    const double&     epsilonZero,
    const double&     epsilonAngle,
    const double&     epsilonZeroGJK,
    const double&     scalingGJK,
    const double&     maxNumPointsPerContact,
    std::ostream&     logStream
)
    :Loggable(logStream),
     mBody1(body1),
     mBody2(body2),
     mInfo(info),
     mEpsilonZero(epsilonZero),
     mEpsilonAngle(epsilonAngle),
     mEpsilonZeroGJK(epsilonZeroGJK),
     mScalingGJK(scalingGJK),
     mMaxNumPointsPerContact(maxNumPointsPerContact)
     {;}


ContactUpdater::~ContactUpdater(){;}


bool ContactUpdater::update()
{
    updateContactNormal();
    // Add direction from Body1 to Body2 to the exit direction to perturb
    // the feature normal to avoid pathetic boundary simplex if two bodies
    // are cleanly parallel.
    Vec3 perturb = mBody2.CoM() - mBody1.CoM();
    perturb.normalize();
    perturb.scale(0.0001);

    if (mInfo.mContactGenerationIrregular) {
        auto noContact = rediscover(true, mInfo.mContactNormal1To2 + perturb);
        return noContact;
    }
    bool needRediscovery, featurePairUpdated;

    checkForUpdate(needRediscovery, featurePairUpdated);

    if (needRediscovery) {
 
        auto noContact = rediscover(true, mInfo.mContactNormal1To2 + perturb);
        return noContact;
    }
    if (featurePairUpdated) {
        updateContactNormal();
    }
    checkFurther(needRediscovery, featurePairUpdated);

    if (needRediscovery) {
        auto noContact = rediscover(true, mInfo.mContactNormal1To2 + perturb);
        return noContact;
    }
    if (featurePairUpdated) {
        updateContactNormal();
    }
    auto noContact = findContactPointPairs();
    if (noContact) {    
        return true;
    }

    return false;
}


Vec3 ContactUpdater::findRelativeVelocityOfBody1RelativeToBody2()
{
    vector<VertexIt> vertices1;
    vector<VertexIt> vertices2;

    decomposeQ(mInfo.mLastVertices, vertices1, vertices2);

    Vec3 r1(0.0, 0.0, 0.0);
    for (auto vit : vertices1) {
        r1 += ((*vit)->pLCS());
    }
    r1.scale(1.0 / (double)(vertices1.size()));
    r1 = mBody1.Qmat() * r1;
    Vec3 r2(0.0, 0.0, 0.0);
    for (auto vit : vertices2) {
        r2 += ((*vit)->pLCS());
    }
    r2.scale(1.0 / (double)(vertices2.size()));
    r2 = mBody2.Qmat() * r2;

    const Vec3 relV = (mBody1.Vlin() + mBody1.Vang().cross(r1)) - 
                      (mBody2.Vlin() + mBody2.Vang().cross(r2)) ;
    return relV;
}


void ContactUpdater::decomposeQ(
    vector<BDVertex>& Q,
    vector<VertexIt>& vertices1,
    vector<VertexIt>& vertices2
) {

    vertices1.clear();
    vertices2.clear();

    for (size_t i = 0; i < Q.size(); i++) {
        bool found = false;
        for (size_t j = 0; j < i; j++) {
            if (Q[j].v1() == Q[i].v1()) {
                found = true;
                break;
            }
        }
        if (!found) {
            vertices1.push_back(Q[i].v1());
        }

        found = false;
        for (size_t j = 0; j < i; j++) {
            if (Q[j].v2() == Q[i].v2()) {
                found = true;
                break;
            }
        }
        if (!found) {
            vertices2.push_back(Q[i].v2());
        }
    }
}


bool ContactUpdater::rediscover(
    bool        preferExternalExitDir,
    const Vec3& exitDir
) {
    makeInitialBDVertices();
    GJKOriginFinder originFinder (mBody1, mBody2, mInfo.mLastVertices, 
                    100, 8, mEpsilonZeroGJK, false, mScalingGJK, mLogStream);
    auto found = originFinder.findOrigin(false);
    if (!found) {
        return true;
    }
                  
    ContactDiscoverer discoverer(
        mBody1,
        mBody2,
        EPSILON_LINEAR,
        EPSILON_LINEAR, 
        0.1,
        0.05,
        0,
        true,
        std::cerr
    );

    mInfo = discoverer.discover();

    if (mInfo.mType1 ==ContactPairInfo::FT_NONE) {
        return true;
    }

    ContactPointsAndNormalGenerator pnGenerator(
        mBody1,
        mBody2,
        mInfo,
        mInfo.mLastVertices,
        mEpsilonZero,
        mEpsilonAngle,
        false, // Not using temp config
        mMaxNumPointsPerContact,
        mLogStream
    );

//    pnGenerator.adjustFeatures();
    pnGenerator.generateContatctPointsAndNormalFromActiveFeatures();
    return mInfo.mIntsect1.size()==0;
}


void ContactUpdater::makeInitialBDVertices()
{
    mInfo.mLastVertices.clear();

    VertexIt vit1, vit2;

    if ( mInfo.mType1 == ContactPairInfo::FT_FACE ) {
        auto heitit1 = (*(mInfo.mFit1))->halfEdges().begin();
        auto heit1   = *heitit1;
        vit1         = (*heit1)->src();
    }
    else if ( mInfo.mType1 == ContactPairInfo::FT_EDGE ) {
        auto he = (*(mInfo.mEit1))->he1();
        vit1    = (*he)->src();
    }
    else {
        vit1 = mInfo.mVit1;
    }

    if ( mInfo.mType2 == ContactPairInfo::FT_FACE ) {
        auto heitit2 = (*(mInfo.mFit2))->halfEdges().begin();
        auto heit2   = *heitit2;
        vit2         = (*heit2)->src();
    }
    else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE ) {
        auto he = (*(mInfo.mEit2))->he1();
        vit2    = (*he)->src();
    }
    else {
        vit2 = mInfo.mVit2;
    }

    BDVertex bd(vit1, vit2);
    bd.updateP( mBody1.ConvexHull(), mBody1.Qmat(), mBody1.CoM(),
                mBody2.ConvexHull(), mBody2.Qmat(), mBody2.CoM() );

    mInfo.mLastVertices.push_back(bd);

}


bool ContactUpdater::findContactPointPairs()
{
    ContactPointsAndNormalGenerator generator(
        mBody1,
        mBody2,
        mInfo,
        mInfo.mLastVertices,
        mEpsilonZero,
        mEpsilonAngle,
        false, // Not using temp config
        mMaxNumPointsPerContact,
        mLogStream
    );

    generator.generateContatctPointsAndNormalFromActiveFeatures();

    return mInfo.mIntsect1.size()==0;
}


void ContactUpdater::checkFurther(
    bool& needRediscovery,
    bool& featurePairUpdated
) {

    auto savedType1 = mInfo.mType1;
    auto savedType2 = mInfo.mType2;
    auto savedVit1  = mInfo.mVit1;
    auto savedVit2  = mInfo.mVit2;
    auto savedEit1  = mInfo.mEit1;
    auto savedEit2  = mInfo.mEit2;
    auto savedFit1  = mInfo.mFit1;
    auto savedFit2  = mInfo.mFit2;

    ContactUpdaterFurtherChecker checker(
                                     mBody1,
                                     mBody2,
                                     mInfo,
                                     mEpsilonZero,
                                     mEpsilonAngle,
                                     mLogStream
                                 );

    needRediscovery = checker.checkFeatures();

    if (needRediscovery) {
        return;
    }

    if ( savedType1 == ContactPairInfo::FT_FACE ) {
        if ( savedType2 == ContactPairInfo::FT_FACE ) {
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_FACE &&
                                   mInfo.mType2 == ContactPairInfo::FT_FACE &&
                                   mInfo.mFit1  == savedFit1                &&
                                   mInfo.mFit2  == savedFit2                 );
        }
        else if ( savedType2 == ContactPairInfo::FT_EDGE ) {
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_FACE &&
                                   mInfo.mType2 == ContactPairInfo::FT_EDGE &&
                                   mInfo.mFit1  == savedFit1                &&
                                   mInfo.mEit2  == savedEit2                 );
        }
        else {
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_FACE &&
                                   mInfo.mType2 == ContactPairInfo::FT_VERTEX&&
                                   mInfo.mFit1  == savedFit1                &&
                                   mInfo.mVit2  == savedVit2                 );
        }
    }
    if ( savedType1 == ContactPairInfo::FT_EDGE ) {
        if ( savedType2 == ContactPairInfo::FT_FACE ) {
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_EDGE &&
                                   mInfo.mType2 == ContactPairInfo::FT_FACE &&
                                   mInfo.mEit1  == savedEit1                &&
                                   mInfo.mFit2  == savedFit2                 );
        }
        else if ( savedType2 == ContactPairInfo::FT_EDGE ) {
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_EDGE &&
                                   mInfo.mType2 == ContactPairInfo::FT_EDGE &&
                                   mInfo.mEit1  == savedEit1                &&
                                   mInfo.mEit2  == savedEit2                 );
        }
        else {
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_EDGE &&
                                   mInfo.mType2 == ContactPairInfo::FT_VERTEX&&
                                   mInfo.mEit1  == savedEit1                &&
                                   mInfo.mVit2  == savedVit2                 );
        }
    }
    else {
        if ( savedType2 == ContactPairInfo::FT_FACE ) {
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_VERTEX&&
                                   mInfo.mType2 == ContactPairInfo::FT_FACE &&
                                   mInfo.mVit1  == savedVit1                &&
                                   mInfo.mFit2  == savedFit2                 );
        }
        else if ( savedType2 == ContactPairInfo::FT_EDGE ) {
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_VERTEX&&
                                   mInfo.mType2 == ContactPairInfo::FT_EDGE &&
                                   mInfo.mVit1  == savedVit1                &&
                                   mInfo.mEit2  == savedEit2                 );
        }
        else {
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_VERTEX&&
                                   mInfo.mType2 == ContactPairInfo::FT_VERTEX&&
                                   mInfo.mVit1  == savedVit1                &&
                                   mInfo.mVit2  == savedVit2                 );
        }
    }
}


void ContactUpdater::checkForUpdate(
    bool& needRediscovery,
    bool& featurePairUpdated
) {
    if ( mInfo.mType1 == ContactPairInfo::FT_FACE ) {
        if ( mInfo.mType2 == ContactPairInfo::FT_FACE ) {
            auto savedFit1 = mInfo.mFit1;
            auto savedFit2 = mInfo.mFit2;
            ContactUpdater_FACE_FACE updater(
                mBody1,
                mBody2,
                mInfo,
                mEpsilonZero,
                mEpsilonAngle,
                mLogStream
            );

            needRediscovery = updater.update();

            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_FACE &&
                                   mInfo.mType2 == ContactPairInfo::FT_FACE &&
                                   mInfo.mFit1  == savedFit1                &&
                                   mInfo.mFit2  == savedFit2                 );
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE ) {
            auto savedFit1 = mInfo.mFit1;
            auto savedEit2 = mInfo.mEit2;
            ContactUpdater_FACE_EDGE updater(
                mBody1,
                mBody2,
                mInfo,
                mEpsilonZero,
                mEpsilonAngle,
                mLogStream
            );
            needRediscovery = updater.update();
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_FACE &&
                                   mInfo.mType2 == ContactPairInfo::FT_EDGE &&
                                   mInfo.mFit1  == savedFit1                &&
                                   mInfo.mEit2  == savedEit2                 );
        }
        else {
            auto savedFit1 = mInfo.mFit1;
            auto savedVit2 = mInfo.mVit2;
            ContactUpdater_FACE_VERTEX updater(
                mBody1,
                mBody2,
                mInfo,
                mEpsilonZero,
                mEpsilonAngle,
                mLogStream
            );
            needRediscovery = updater.update();
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_FACE &&
                                   mInfo.mType2 == ContactPairInfo::FT_VERTEX&&
                                   mInfo.mFit1  == savedFit1                &&
                                   mInfo.mVit2  == savedVit2                 );
        }
    }
    else if ( mInfo.mType1 == ContactPairInfo::FT_EDGE ) {
        if ( mInfo.mType2 == ContactPairInfo::FT_FACE ) {
            auto savedEit1 = mInfo.mEit1;
            auto savedFit2 = mInfo.mFit2;
            ContactUpdater_FACE_EDGE updater(
                mBody1,
                mBody2,
                mInfo,
                mEpsilonZero,
                mEpsilonAngle,
                mLogStream
            );
            needRediscovery = updater.update();
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_EDGE &&
                                   mInfo.mType2 == ContactPairInfo::FT_FACE &&
                                   mInfo.mEit1  == savedEit1                &&
                                   mInfo.mFit2  == savedFit2                 );
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE ) {
            auto savedEit1 = mInfo.mEit1;
            auto savedEit2 = mInfo.mEit2;
            ContactUpdater_EDGE_EDGE updater(
                mBody1,
                mBody2,
                mInfo,
                mEpsilonZero,
                mEpsilonAngle,
                mLogStream
            );
            needRediscovery = updater.update();
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_EDGE &&
                                   mInfo.mType2 == ContactPairInfo::FT_EDGE &&
                                   mInfo.mEit1  == savedEit1                &&
                                   mInfo.mEit2  == savedEit2                 );
        }
        else {
            auto savedEit1 = mInfo.mEit1;
            auto savedVit2 = mInfo.mVit2;
            ContactUpdater_EDGE_VERTEX updater(
                mBody1,
                mBody2,
                mInfo,
                mEpsilonZero,
                mEpsilonAngle,
                mLogStream
            );
            needRediscovery = updater.update();
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_EDGE &&
                                   mInfo.mType2 == ContactPairInfo::FT_VERTEX&&
                                   mInfo.mEit1  == savedEit1                &&
                                   mInfo.mVit2  == savedVit2                 );
        }
    }
    else {
        if ( mInfo.mType2 == ContactPairInfo::FT_FACE ) {
            auto savedVit1 = mInfo.mVit1;
            auto savedFit2 = mInfo.mFit2;
            ContactUpdater_FACE_VERTEX updater(
                mBody1,
                mBody2,
                mInfo,
                mEpsilonZero,
                mEpsilonAngle,
                mLogStream
            );
            needRediscovery = updater.update();
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_VERTEX&&
                                   mInfo.mType2 == ContactPairInfo::FT_FACE &&
                                   mInfo.mVit1  == savedVit1                &&
                                   mInfo.mFit2  == savedFit2                 );
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE ) {
            auto savedVit1 = mInfo.mVit1;
            auto savedEit2 = mInfo.mEit2;
            ContactUpdater_EDGE_VERTEX updater(
                mBody1,
                mBody2,
                mInfo,
                mEpsilonZero,
                mEpsilonAngle,
                mLogStream
            );
            needRediscovery = updater.update();
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_VERTEX&&
                                   mInfo.mType2 == ContactPairInfo::FT_EDGE &&
                                   mInfo.mVit1  == savedVit1                &&
                                   mInfo.mEit2  == savedEit2                 );
        }
        else {
            auto savedVit1 = mInfo.mVit1;
            auto savedVit2 = mInfo.mVit2;
            ContactUpdater_VERTEX_VERTEX updater(
                mBody1,
                mBody2,
                mInfo,
                mEpsilonZero,
                mEpsilonAngle,
                mLogStream
            );
            needRediscovery = updater.update();
            featurePairUpdated =!( mInfo.mType1 == ContactPairInfo::FT_VERTEX&&
                                   mInfo.mType2 == ContactPairInfo::FT_VERTEX&&
                                   mInfo.mVit1  == savedVit1                &&
                                   mInfo.mVit2  == savedVit2                 );
        }
    }
}


void ContactUpdater::updateContactNormal()
{
    if ( mInfo.mType1 == ContactPairInfo::FT_FACE) {
        mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_1;
        mInfo.mContactNormal1To2 = (*mInfo.mFit1)->nGCS(mBody1.Qmat());
    }
    else if ( mInfo.mType1 == ContactPairInfo::FT_EDGE) {
        if ( mInfo.mType2 == ContactPairInfo::FT_FACE) {
            mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;
            mInfo.mContactNormal1To2 = 
                             (*mInfo.mFit2)->nGCS(mBody2.Qmat()) * -1.0;
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE) {
            auto vit1src = (*((*mInfo.mEit1)->he1()))->src();
            auto vit1dst = (*((*mInfo.mEit1)->he1()))->dst();
            auto vit2src = (*((*mInfo.mEit2)->he1()))->src();
            auto vit2dst = (*((*mInfo.mEit2)->he1()))->dst();

            auto d1      = (*vit1dst)->pGCS(mBody1.Qmat()) -
                           (*vit1src)->pGCS(mBody1.Qmat()) ;

            auto d2      = (*vit2dst)->pGCS(mBody2.Qmat()) -
                           (*vit2src)->pGCS(mBody2.Qmat()) ;

            auto n1      =  (*mInfo.mEit1)->nGCS(mBody1.Qmat());
            auto n2      =  (*mInfo.mEit2)->nGCS(mBody2.Qmat());
            auto nAvg    = n1 - n2;
            nAvg.normalize();
            auto cr = d1.cross(d2);
            if (cr.squaredNorm2() <= mEpsilonAngle) {
                // two edges are parallel. Take average of the two edge normals
                mInfo.mContactNormal1To2 = nAvg;
                mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;
            }
            else {
                cr.normalize();

                if (fabs(nAvg.dot(cr)) >= sqrt(2.0)/2.0) {

                    if (cr.dot(mBody2.CoM() - mBody1.CoM())< 0.0) {

                        cr.scale(-1.0);
                    }
                    mInfo.mContactNormal1To2 = cr;
                    mInfo.mContactNormalType = 
                                          ContactPairInfo::NT_EDGE_CROSS_EDGE;
                }
                else {

                    mInfo.mContactNormal1To2 = nAvg;
                    mInfo.mContactNormalType = 
                                          ContactPairInfo::NT_EDGE_EDGE_AVG;
                }
            }
        }
        else {
            mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_NORMAL_1;
            mInfo.mContactNormal1To2 = (*mInfo.mEit1)->nGCS(mBody1.Qmat());
        }
    }
    else {
        if ( mInfo.mType2 == ContactPairInfo::FT_FACE) {
            mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;
            mInfo.mContactNormal1To2 = 
                              (*mInfo.mFit2)->nGCS(mBody2.Qmat()) * -1.0;
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE) {
            mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_NORMAL_2;
            mInfo.mContactNormal1To2 = 
                              (*mInfo.mEit2)->nGCS(mBody2.Qmat()) * -1.0;
        }
        else {
            mInfo.mContactNormal1To2 =
                (*mInfo.mVit1)->nGCS(mBody1.Qmat()) + 
                (*mInfo.mVit2)->nGCS(mBody2.Qmat()) * -1.0;
            mInfo.mContactNormal1To2.normalize();
            mInfo.mContactNormalType = ContactPairInfo::NT_VERTEX_VERTEX_AVG;
        }
    }
}


#ifdef UNIT_TESTS
void ContactUpdater::showResult()
{
    cerr << "Body 1\n";
    if (mInfo.mType1 == ContactPairInfo::FT_FACE) {
        cerr << "Face:\n";
        for (auto he : (*(mInfo.mFit1))->halfEdges()) {
            auto vit = (*he)->src();
            cerr << (*vit)->pLCS() << "\n";
        }
    }
    else if (mInfo.mType1 == ContactPairInfo::FT_EDGE) {
        cerr << "Edge:\n";
        auto he = (*(mInfo.mEit1))->he1();
        auto vit1 = (*he)->src();
        auto vit2 = (*he)->dst();
        cerr << (*vit1)->pLCS() << "\n";
        cerr << (*vit2)->pLCS() << "\n";
    }
    else {
        cerr << "Vertex: ";
        cerr << (*(mInfo.mVit1))->pLCS() << "\n";
    }

    cerr << "Body 2\n";
    if (mInfo.mType2 == ContactPairInfo::FT_FACE) {
        cerr << "Face:\n";
        for (auto he : (*(mInfo.mFit2))->halfEdges()) {
            auto vit = (*he)->src();
            cerr << (*vit)->pLCS() << "\n";
        }
    }
    else if (mInfo.mType2 == ContactPairInfo::FT_EDGE) {
        cerr << "Edge:\n";
        auto he = (*(mInfo.mEit2))->he1();
        auto vit1 = (*he)->src();
        auto vit2 = (*he)->dst();
        cerr << (*vit1)->pLCS() << "\n";
        cerr << (*vit2)->pLCS() << "\n";
    }
    else {
        cerr << "Vertex: ";
        cerr << (*(mInfo.mVit2))->pLCS() << "\n";
    }

    cerr << "Intersection1:\n";
    for (auto& p : mInfo.mIntsect1) {
        cerr << p << "\n";
    }
    cerr << "Intersection2:\n";
    for (auto& p : mInfo.mIntsect2) {
        cerr << p << "\n";
    }
}
#endif


}// namespace Makena

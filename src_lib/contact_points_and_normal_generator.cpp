#include "contact_points_and_normal_generator.hpp"
#include "intersection_convex_polygon_2d.hpp"
#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file gjk_contact_generator.cpp
 *
 * @brief
 */
namespace Makena {

using namespace std;

ContactPointsAndNormalGenerator::ContactPointsAndNormalGenerator(
    ConvexRigidBody&         body1,
    ConvexRigidBody&         body2,
    ContactPairInfo&         info,
    vector<BDVertex>&        Q,
    const double&            epsilonZero,
    const double&            epsilonAngle,
    bool                     useGeomConfigTemp,
    const long               maxNumPointsPerContact,
    std::ostream&            logStream
)
    :Loggable(logStream),
     mBody1(body1),
     mBody2(body2),
     mInfo(info),
     mQ(Q),
     mEpsilonZero(epsilonZero),
     mEpsilonAngle(epsilonAngle),
     mUseGeomConfigTemp(useGeomConfigTemp),
     mMaxNumPointsPerContact(maxNumPointsPerContact),
     mQmat1(useGeomConfigTemp?(mBody1.QmatTemp()):(mBody1.Qmat())),
     mQmat2(useGeomConfigTemp?(mBody2.QmatTemp()):(mBody2.Qmat())),
     mCoM1(useGeomConfigTemp?(mBody1.CoMTemp()):(mBody1.CoM())),
     mCoM2(useGeomConfigTemp?(mBody2.CoMTemp()):(mBody2.CoM()))
     {;}


ContactPointsAndNormalGenerator::~ContactPointsAndNormalGenerator(){;}


void ContactPointsAndNormalGenerator::findActiveFeaturesFromBinaryDilation()
{
    vector<VertexIt> vits1;
    vector<VertexIt> vits2;

    findVerticesOfRigidBody(vits1, vits2);

    findFeatureOfRigidBody(mBody1.ConvexHull(), 
                  vits1, mInfo.mFit1, mInfo.mEit1, mInfo.mVit1, mInfo.mType1);

    findFeatureOfRigidBody(mBody2.ConvexHull(),
                  vits2, mInfo.mFit2, mInfo.mEit2, mInfo.mVit2, mInfo.mType2);
}


void ContactPointsAndNormalGenerator::
generateContatctPointsAndNormalFromActiveFeatures()
{
    mInfo.mIntsect1.clear();
    mInfo.mIntsect2.clear();
    mInfo.mIntsect1reduced.clear();
    mInfo.mIntsect2reduced.clear();
    mInfo.mContactGenerationIrregular = false;

    if ( mInfo.mType1 == ContactPairInfo::FT_FACE) {
        if ( mInfo.mType2 == ContactPairInfo::FT_FACE) {
            process_FACE_FACE();
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE) {

            process_FACE_EDGE();
        }
        else {
            process_FACE_VERTEX();
        }
    }
    else if ( mInfo.mType1 == ContactPairInfo::FT_EDGE) {
        if ( mInfo.mType2 == ContactPairInfo::FT_FACE) {
            process_EDGE_FACE();
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE) {
            process_EDGE_EDGE();
        }
        else {
            process_EDGE_VERTEX();
        }
    }
    else {
        if ( mInfo.mType2 == ContactPairInfo::FT_FACE) {
            process_VERTEX_FACE();
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE) {
            process_VERTEX_EDGE();
        }
        else {
            process_VERTEX_VERTEX();
        }
    }
}


Vec3 ContactPointsAndNormalGenerator::
findRelativeVelocityOfBody1RelativeToBody2()
{
    Vec3 r1(0.0, 0.0, 0.0);
    if ( mInfo.mType1 == ContactPairInfo::FT_FACE) {
        auto& hes = (*(mInfo.mFit1))->halfEdges();
        for (auto he : hes) {
            auto vit = (*he)->src();
            r1 += ((*vit)->pLCS());
        }
        r1.scale(1.0/hes.size());
    }
    else if ( mInfo.mType1 == ContactPairInfo::FT_EDGE) {
        auto he = (*(mInfo.mEit1))->he1();
        auto vit1 = (*he)->src();
        auto vit2 = (*he)->dst();
        r1 += ((*vit1)->pLCS());
        r1 += ((*vit2)->pLCS());
        r1.scale(0.5);
    }
    else {
        r1 = (*(mInfo.mVit1))->pLCS();
    }

    Vec3 r2(0.0, 0.0, 0.0);
    if ( mInfo.mType2 == ContactPairInfo::FT_FACE) {
        auto& hes = (*(mInfo.mFit2))->halfEdges();
        for (auto he : hes) {
            auto vit = (*he)->src();
            r2 += ((*vit)->pLCS());
        }
        r1.scale(1.0/hes.size());
    }
    else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE) {
        auto he = (*(mInfo.mEit2))->he1();
        auto vit1 = (*he)->src();
        auto vit2 = (*he)->dst();
        r2 += ((*vit1)->pLCS());
        r2 += ((*vit2)->pLCS());
        r2.scale(0.5);
    }
    else {
        r2 = (*(mInfo.mVit2))->pLCS();
    }

    const Vec3 relV = (mBody1.Vlin() + mBody1.Vang().cross(r1)) - 
                      (mBody2.Vlin() + mBody2.Vang().cross(r2)) ;
    return relV;
}



void ContactPointsAndNormalGenerator::process_FACE_FACE()
{
    // Attempt 1: Face Normal 1
    Vec3 normalDir1 = (*mInfo.mFit1)->nGCS(mQmat1);

    process_FACE_FACE_tryOneDirection(normalDir1);
    if (mInfo.mIntsect1.size() > 0) {
        mInfo.mContactNormal1To2 = normalDir1;
        mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_1;
        return;
    }

    mInfo.mContactGenerationIrregular = true;

    // Attempt 2: Face Normal 2
    Vec3 normalDir2 = (*mInfo.mFit2)->nGCS(mQmat2);
                         
    normalDir2.scale(-1.0);
    process_FACE_FACE_tryOneDirection(normalDir2);
    if (mInfo.mIntsect1.size() > 0) {
        mInfo.mContactNormal1To2 = normalDir2;
        mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;
        return;
    }



    // Attempt 3: Relative Velocity of Body 1 to 2.
    log(INFO, __FILE__, __LINE__,
        "Can't find intersection. Trying Relative Velocity");
    Vec3 relVel = findRelativeVelocityOfBody1RelativeToBody2();
    relVel.normalize();
    if (relVel.dot(normalDir1) >= mEpsilonAngle && 
        relVel.dot(normalDir2) >= mEpsilonAngle   ) {

        process_FACE_FACE_tryOneDirection(relVel);
        if (mInfo.mIntsect1.size() > 0) {
            mInfo.mContactNormal1To2 = relVel;
            mInfo.mContactNormalType = ContactPairInfo::NT_RELATIVE_VELOCITY;
            return;
        }
        // Back-up: Projection of Body 2 onto the plane of Body 1.
        log(INFO, __FILE__, __LINE__,
            "Can't find intersection. Using Projection");
        process_FACE_FACE_projection();    
    }
    else {
        log(INFO, __FILE__, __LINE__,
            "Can't find intersection. Using Projection");
        process_FACE_FACE_projection();    
    }
}



void ContactPointsAndNormalGenerator::process_FACE_FACE_tryOneDirection(
    const Vec3& normalDir
) {
    Mat3x3 rotMat1,  rotMat2;
    Vec3   com1,     com2;

    rotateWorld(normalDir, rotMat1, com1, rotMat2, com2);

    vector<IntersectionFinderConvexPolygon2D::InputElem> points1;
    for (auto& he : (*mInfo.mFit1)->halfEdges()) {
        auto p3 = (*((*he)->src()))->pGCS(rotMat1, com1);
        Vec2 p2(p3.x(), p3.y());
        points1.push_back(IntersectionFinderConvexPolygon2D::InputElem(p2, 0));
    }

    vector<IntersectionFinderConvexPolygon2D::InputElem> points2;
    for (auto& he : (*mInfo.mFit2)->halfEdges()) {
        auto p3 = (*((*he)->src()))->pGCS(rotMat2, com2);
        Vec2 p2(p3.x(), p3.y());
        points2.push_back(IntersectionFinderConvexPolygon2D::InputElem(p2, 0));
    }
    std::reverse(points2.begin(), points2.end());

    vector<IntersectionFinderConvexPolygon2D::OutputElem> intsec;
    IntersectionFinderConvexPolygon2D finder;

    auto found = finder.findIntersection(points1, points2, intsec);

    if (!found) {
        return;
    }
    from2DXYto3DPoints(intsec, rotMat1, com1, rotMat2, com2);

}


void ContactPointsAndNormalGenerator::process_FACE_FACE_projection()
{
    // Body1's normal
    Vec3 normalDir1 = (*mInfo.mFit1)->nGCS(mQmat1);

    Mat3x3 rotMat11,  rotMat12;
    Vec3   com11,     com12;
    rotateWorld(normalDir1, rotMat11, com11, rotMat12, com12);

    Vec3         rotatedPoint1;
    vector<Vec3> rotatedPoints2;
    auto he1 = *((*mInfo.mFit1)->halfEdges()).begin();
    rotatedPoint1 = ((*((*he1)->src()))->pGCS(rotMat11, com11));
    for (auto& he2 : (*mInfo.mFit2)->halfEdges()) {
        rotatedPoints2.push_back((*((*he2)->src()))->pGCS(rotMat12, com12));
    }
    double maxDepth1 = 0.0;
    for (auto& p : rotatedPoints2) {
        maxDepth1 = std::min(maxDepth1, rotatedPoint1.z() - p.z());
        p.setZ(rotatedPoint1.z());
    }

    // Body2's normal. We don't negate it yet to check the depth and compare
    // it with the other configuration below.
    Vec3 normalDir2 = (*mInfo.mFit2)->nGCS(mQmat2);

    Mat3x3 rotMat21,  rotMat22;
    Vec3   com21,     com22;
    rotateWorld(normalDir2, rotMat21, com21, rotMat22, com22);

    Vec3         rotatedPoint2;
    vector<Vec3> rotatedPoints1;
    auto he2 = *((*mInfo.mFit2)->halfEdges()).begin();
    rotatedPoint2 = ((*((*he2)->src()))->pGCS(rotMat22, com22));

    for (auto& he1 : (*mInfo.mFit1)->halfEdges()) {
        rotatedPoints1.push_back((*((*he1)->src()))->pGCS(rotMat21, com21));
    }
    double maxDepth2 = 0.0;
    for (auto& p : rotatedPoints1) {
        maxDepth2 = std::min(maxDepth2, rotatedPoint2.z() - p.z());
        p.setZ(rotatedPoint2.z());
    }

    if (maxDepth1 > maxDepth2) {
        // Depth 1 is shallower. Project Face 2 onto Plane of Face 1.
        Mat3x3 Qmat1inv = rotMat11.transpose();
        for (auto& p: rotatedPoints2) {
            // Those points are outside Face 1. (Projection of Face 2)
            mInfo.mIntsect1.push_back(Qmat1inv * (p - com11));
        }
        for (auto& he2 : (*mInfo.mFit2)->halfEdges()) {
            mInfo.mIntsect2.push_back((*((*he2)->src()))->pLCS());
        }

        mInfo.mContactNormal1To2 = normalDir1;
        mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_1;

    }
    else {
        // Depth 2 is shallower. Project Face 1 onto Plane of Face 2.
        Mat3x3 Qmat2inv = rotMat22.transpose();
        for (auto& p: rotatedPoints1) {
            // Those points are outside Face 2. (Projection of Face 1)
            mInfo.mIntsect2.push_back(Qmat2inv * (p - com22));
        }
        for (auto& he1 : (*mInfo.mFit1)->halfEdges()) {
            mInfo.mIntsect1.push_back((*((*he1)->src()))->pLCS());
        }
        mInfo.mContactNormal1To2 = normalDir2 * -1.0;
        mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;

    }
    mInfo.mIntsect1reduced = mInfo.mIntsect1;
    mInfo.mIntsect2reduced = mInfo.mIntsect2;
}



/** @brief
 * 
 *   Strategy: If one of the incident faces on Body 2is penetrating to the
 *             face of Body 1, then try face-face, otherwise, face-edge.
 */
void ContactPointsAndNormalGenerator::process_FACE_EDGE()
{

    // Attempt 1: Face Normal 1
    Vec3 normalDir1 = (*mInfo.mFit1)->nGCS(mQmat1);

    process_FACE_EDGE_tryOneDirection(normalDir1);
    if (mInfo.mIntsect1.size() > 0) {
        mInfo.mContactNormal1To2 = normalDir1;
        mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_1;
        return;
    }

    mInfo.mContactGenerationIrregular = true;

    // Attempt 2: Edge Normal 2
    Vec3 normalDir2 = (*mInfo.mEit2)->nGCS(mQmat2);
    normalDir2.scale(-1.0);

    process_FACE_EDGE_tryOneDirection(normalDir2);
    if (mInfo.mIntsect1.size() > 0) {
        mInfo.mContactNormal1To2 = normalDir2;
        mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_NORMAL_2;
        return;
    }



    // Attempt 3: Relative Velocity of Body 1 to 2.
    log(INFO, __FILE__, __LINE__,
        "Can't find intersection. Trying Relative Velocity");

    Vec3 relVel = findRelativeVelocityOfBody1RelativeToBody2();
    relVel.normalize();
    process_FACE_EDGE_tryOneDirection(relVel);

    if (mInfo.mIntsect1.size() > 0) {
        mInfo.mContactNormal1To2 = relVel;
        mInfo.mContactNormalType = ContactPairInfo::NT_RELATIVE_VELOCITY;
        return;
    }


    // Back-up: Projection of Body 2 onto the plane of Body 1.
    log(INFO, __FILE__, __LINE__,
        "Can't find intersection. Using Projection");
    process_FACE_EDGE_projection();    

}


void ContactPointsAndNormalGenerator::process_FACE_EDGE_tryOneDirection(
    const Vec3& normalDir
) {
    Mat3x3 rotMat1,  rotMat2;
    Vec3   com1,     com2;

    rotateWorld(normalDir, rotMat1, com1, rotMat2, com2);

    vector<IntersectionFinderConvexPolygon2D::InputElem> points1;
    for (auto& he : (*mInfo.mFit1)->halfEdges()) {
        auto p3 = (*((*he)->src()))->pGCS(rotMat1, com1);
        Vec2 p2(p3.x(), p3.y());
        points1.push_back(IntersectionFinderConvexPolygon2D::InputElem(p2, 0));

    }
    vector<IntersectionFinderConvexPolygon2D::InputElem> points2;
    auto he21  = (*(mInfo.mEit2))->he1();
    auto vit21 = (*he21)->src();
    auto vit22 = (*he21)->dst();
    auto p21   = (*vit21)->pGCS(rotMat2, com2);
    auto p22   = (*vit22)->pGCS(rotMat2, com2);
    Vec2 p21_2D(p21.x(), p21.y());
    Vec2 p22_2D(p22.x(), p22.y());
    points2.push_back(IntersectionFinderConvexPolygon2D::InputElem(p21_2D, 0));
    points2.push_back(IntersectionFinderConvexPolygon2D::InputElem(p22_2D, 1));
    vector<IntersectionFinderConvexPolygon2D::OutputElem> intsec;
    IntersectionFinderConvexPolygon2D finder;
    auto found = finder.findIntersection(points1, points2, intsec);
    if (!found) {
        return;
    }
    from2DXYto3DPoints(intsec, rotMat1, com1, rotMat2, com2);
}


void ContactPointsAndNormalGenerator::process_FACE_EDGE_projection()
{
    // Body1's normal
    Vec3 normalDir1 = (*mInfo.mFit1)->nGCS(mQmat1);

    Mat3x3 rotMat1,  rotMat2;
    Vec3   com1,     com2;
    rotateWorld(normalDir1, rotMat1, com1, rotMat2, com2);

    auto he1 = *((*mInfo.mFit1)->halfEdges()).begin();
    auto rotatedPoint1 = ((*((*he1)->src()))->pGCS(rotMat1, com1));

    auto he2   = (*(mInfo.mEit2))->he1();
    auto rotatedPoint21 = (*((*he2)->src()))->pGCS(rotMat2, com2);
    auto rotatedPoint22 = (*((*he2)->dst()))->pGCS(rotMat2, com2);
    rotatedPoint21.setZ(rotatedPoint1.z());
    rotatedPoint22.setZ(rotatedPoint1.z());

    Mat3x3 Qmat1inv = rotMat1.transpose();
    // Those points are outside Face 1. (Projection of Edge 2)
    mInfo.mIntsect1.push_back(Qmat1inv * (rotatedPoint21 - com1));
    mInfo.mIntsect1.push_back(Qmat1inv * (rotatedPoint22 - com1));

    mInfo.mIntsect2.push_back((*((*he2)->src()))->pLCS());
    mInfo.mIntsect2.push_back((*((*he2)->dst()))->pLCS());

    mInfo.mContactNormal1To2 = normalDir1;
    mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_1;

    mInfo.mIntsect1reduced = mInfo.mIntsect1;
    mInfo.mIntsect2reduced = mInfo.mIntsect2;

}


/** @brief
 * 
 *   Strategy: If one of the incident faces on Body 1 is penetrating to the
 *             face of Body 2, then try face-face, otherwise, edge-face.
 */
void ContactPointsAndNormalGenerator::process_EDGE_FACE()
{
    // Attempt 1: Face Normal 2
    Vec3 normalDir1 = (*mInfo.mFit2)->nGCS(mQmat2);
    normalDir1.scale(-1.0);
    process_EDGE_FACE_tryOneDirection(normalDir1);
    if (mInfo.mIntsect1.size() > 0) {
        mInfo.mContactNormal1To2 = normalDir1;
        mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;
        return;
    }

    mInfo.mContactGenerationIrregular = true;

    // Attempt 2: Edge Normal 1
    Vec3 normalDir2 = (*mInfo.mEit1)->nGCS(mQmat1);

    process_EDGE_FACE_tryOneDirection(normalDir2);
    if (mInfo.mIntsect1.size() > 0) {
        mInfo.mContactNormal1To2 = normalDir2;
        mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_NORMAL_1;
        return;
    }

    // Attempt 3: Relative Velocity
    log(INFO, __FILE__, __LINE__,
            "Can't find intersection. Trying Relative Velocity");
    Vec3 relVel = findRelativeVelocityOfBody1RelativeToBody2();
    relVel.normalize();
    process_EDGE_FACE_tryOneDirection(relVel);
    if (mInfo.mIntsect1.size() > 0) {
        mInfo.mContactNormal1To2 = relVel;
        mInfo.mContactNormalType = ContactPairInfo::NT_RELATIVE_VELOCITY;
        return;
    }

    // Back-up: Projection of Body 1 onto the plane of Body 2.
    log(INFO, __FILE__, __LINE__,
        "Can't find intersection. Using Projection");
    process_EDGE_FACE_projection();    
}


void ContactPointsAndNormalGenerator::process_EDGE_FACE_tryOneDirection(
    const Vec3& normalDir
) {

    Mat3x3 rotMat1,  rotMat2;
    Vec3   com1,     com2;

    rotateWorld(normalDir, rotMat1, com1, rotMat2, com2);
    vector<IntersectionFinderConvexPolygon2D::InputElem> points1;
    auto he11  = (*(mInfo.mEit1))->he1();
    auto vit11 = (*he11)->src();
    auto vit12 = (*he11)->dst();
    auto p11   = (*vit11)->pGCS(rotMat1, com1);
    auto p12   = (*vit12)->pGCS(rotMat1, com1);
    Vec2 p11_2D(p11.x(), p11.y());
    Vec2 p12_2D(p12.x(), p12.y());
    points1.push_back(IntersectionFinderConvexPolygon2D::InputElem(p11_2D, 0));
    points1.push_back(IntersectionFinderConvexPolygon2D::InputElem(p12_2D, 1));

    vector<IntersectionFinderConvexPolygon2D::InputElem> points2;
    for (auto& he : (*mInfo.mFit2)->halfEdges()) {
        auto p3 = (*((*he)->src()))->pGCS(rotMat2, com2);
        Vec2 p2(p3.x(), p3.y());
        points2.push_back(IntersectionFinderConvexPolygon2D::InputElem(p2, 0));
    }
    std::reverse(points2.begin(), points2.end());
    vector<IntersectionFinderConvexPolygon2D::OutputElem> intsec;
    IntersectionFinderConvexPolygon2D finder;
    auto found = finder.findIntersection(points1, points2, intsec);
    if (!found) {
        return;
    }
    from2DXYto3DPoints(intsec, rotMat1, com1, rotMat2, com2);
}


void ContactPointsAndNormalGenerator::process_EDGE_FACE_projection()
{
    // Body2's normal
    Vec3 normalDir2 = (*mInfo.mFit2)->nGCS(mQmat2);

    Mat3x3 rotMat1,  rotMat2;
    Vec3   com1,     com2;
    rotateWorld(normalDir2, rotMat1, com1, rotMat2, com2);

    auto he2 = *((*mInfo.mFit2)->halfEdges()).begin();
    auto rotatedPoint2 = ((*((*he2)->src()))->pGCS(rotMat2, com2));

    auto he1   = (*(mInfo.mEit1))->he1();
    auto rotatedPoint11 = (*((*he1)->src()))->pGCS(rotMat1, com1);
    auto rotatedPoint12 = (*((*he1)->dst()))->pGCS(rotMat1, com1);
    rotatedPoint11.setZ(rotatedPoint2.z());
    rotatedPoint12.setZ(rotatedPoint2.z());

    Mat3x3 Qmat2inv = rotMat2.transpose();
    // Those points are outside Face 2. (Projection of Edge 1)
    mInfo.mIntsect2.push_back(Qmat2inv * (rotatedPoint11 - com2));
    mInfo.mIntsect2.push_back(Qmat2inv * (rotatedPoint12 - com2));

    mInfo.mIntsect1.push_back((*((*he1)->src()))->pLCS());
    mInfo.mIntsect1.push_back((*((*he1)->dst()))->pLCS());

    mInfo.mContactNormal1To2 = normalDir2 * -1.0;
    mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;

    mInfo.mIntsect1reduced = mInfo.mIntsect1;
    mInfo.mIntsect2reduced = mInfo.mIntsect2;

}


void ContactPointsAndNormalGenerator::process_FACE_VERTEX()
{
    // Project the vertex on Body 2 onto the plane of the face on Body 1.

    // Body1's normal
    Vec3 normalDir1 = (*mInfo.mFit1)->nGCS(mQmat1);

    Mat3x3 rotMat1,  rotMat2;
    Vec3   com1,     com2;
    rotateWorld(normalDir1, rotMat1, com1, rotMat2, com2);

    auto he1 = *((*mInfo.mFit1)->halfEdges()).begin();
    auto rotatedPoint1 = ((*((*he1)->src()))->pGCS(rotMat1, com1));

    auto rotatedPoint2 = (*(mInfo.mVit2))->pGCS(rotMat2, com2);
    rotatedPoint2.setZ(rotatedPoint1.z());

    Mat3x3 Qmat1inv = rotMat1.transpose();
    mInfo.mIntsect1.push_back(Qmat1inv * (rotatedPoint2 - com1));
    mInfo.mIntsect2.push_back((*(mInfo.mVit2))->pLCS());

    mInfo.mContactNormal1To2 = normalDir1;
    mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_1;

    mInfo.mIntsect1reduced = mInfo.mIntsect1;
    mInfo.mIntsect2reduced = mInfo.mIntsect2;

}


void ContactPointsAndNormalGenerator::process_EDGE_EDGE()
{
    auto body1_he1      = (*mInfo.mEit1)->he1();
    auto body1_vit1     = (*body1_he1)->src();
    auto body1_vit2     = (*body1_he1)->dst();
    const auto body1_p1 = (*body1_vit1)->pGCS(mBody1.Qmat(), mBody1.CoM());
    const auto body1_p2 = (*body1_vit2)->pGCS(mBody1.Qmat(), mBody1.CoM());
    auto body1_v12  = body1_p2 - body1_p1;
    body1_v12.normalize();
    auto body2_he1      = (*mInfo.mEit2)->he1();
    auto body2_vit1     = (*body2_he1)->src();
    auto body2_vit2     = (*body2_he1)->dst();
    const auto body2_p1 = (*body2_vit1)->pGCS(mBody2.Qmat(), mBody2.CoM());
    const auto body2_p2 = (*body2_vit2)->pGCS(mBody2.Qmat(), mBody2.CoM());
    auto body2_v12  = body2_p2 - body2_p1;
    body2_v12.normalize();

    auto body1_en       = (*mInfo.mEit1)->nGCS(mBody1.Qmat());
    auto body2_en       = (*mInfo.mEit2)->nGCS(mBody2.Qmat());
    auto enAvg          = body1_en - body2_en;
    enAvg.normalize();

    auto cr = body1_v12.cross(body2_v12);
    if (cr.squaredNorm2() <= mEpsilonAngle) {
        process_EDGE_EDGE_PARALLEL(body1_v12);

    }
    else {
        cr.normalize();
        auto dot = cr.dot(enAvg);
        
        if ( fabs(dot) >= sqrt(2.0)/2.0 ) {
            if (dot < 0.0) {
                cr.scale(-1.0);
            }
            // the angle between the cross and the edge normal <= pi/2.0
            process_EDGE_EDGE_CROSSING(cr);

        }
        else {
            process_EDGE_EDGE_PARALLEL(body1_v12);

        }
    }
}


void ContactPointsAndNormalGenerator::process_EDGE_EDGE_CROSSING(
    const Vec3& zDir
) {

    // Find the closest point of two edges.
    Mat3x3 rotMat1,  rotMat2;
    Vec3   com1,     com2;
    rotateWorld(zDir, rotMat1, com1, rotMat2, com2);

    auto body1_he1      = (*mInfo.mEit1)->he1();
    auto body1_vit1     = (*body1_he1)->src();
    auto body1_vit2     = (*body1_he1)->dst();
    const auto body1_p1 = (*body1_vit1)->pGCS(rotMat1, com1);
    const auto body1_p2 = (*body1_vit2)->pGCS(rotMat1, com1);

    auto body2_he1      = (*mInfo.mEit2)->he1();
    auto body2_vit1     = (*body2_he1)->src();
    auto body2_vit2     = (*body2_he1)->dst();
    const auto body2_p1 = (*body2_vit1)->pGCS(rotMat2, com2);
    const auto body2_p2 = (*body2_vit2)->pGCS(rotMat2, com2);

    Vec2 p11(body1_p1.x(), body1_p1.y());
    Vec2 p12(body1_p2.x(), body1_p2.y());
    Vec2 v1 = p12 - p11;
    Vec2 p21(body2_p1.x(), body2_p1.y());
    Vec2 p22(body2_p2.x(), body2_p2.y());
    Vec2 v2 = p22 - p21;    

    Vec2 v12 = p21 - p11;
    Vec2 v2perp = v2.perp();
    double s = v12.dot(v2perp) / v1.dot(v2perp);
    Vec2 intsec = p11 + v1 * s;
    Vec3 p1(intsec.x(), intsec.y(), body1_p1.z());
    Vec3 p2(intsec.x(), intsec.y(), body2_p1.z());

    Mat3x3 Qmat1inv = rotMat1.transpose();
    Mat3x3 Qmat2inv = rotMat2.transpose();
    mInfo.mIntsect1.push_back(Qmat1inv * (p1 - com1));
    mInfo.mIntsect2.push_back(Qmat2inv * (p2 - com2));

    if (zDir.dot(mCoM2 - mCoM1)< 0.0) {
        mInfo.mContactNormal1To2 = zDir * -1.0;
    }
    else {                
        mInfo.mContactNormal1To2 = zDir;
    }
    mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_CROSS_EDGE;

    mInfo.mIntsect1reduced = mInfo.mIntsect1;
    mInfo.mIntsect2reduced = mInfo.mIntsect2;


}


Vec3 ContactPointsAndNormalGenerator::interpolate(
    const Vec3&  p1,
    const Vec3&  p2,
    const double z
) {
    const double zSpan = p2.z() - p1.z();
    if (fabs(zSpan)< mEpsilonZero) {
        return p1;
    }
    const double s = (z - p1.z())/ zSpan;

    return Vec3( s * (p2.x() - p1.x()) + p1.x(),
                 s * (p2.y() - p1.y()) + p1.y(),
                 z                               );
}


void ContactPointsAndNormalGenerator::process_EDGE_EDGE_PARALLEL(
    const Vec3& axisDir
) {

    Mat3x3 rotMat1,  rotMat2;
    Vec3   com1,     com2;
    rotateWorld(axisDir, rotMat1, com1, rotMat2, com2);

    // Align the edge along z-axis and find the overlap.
    auto body1_he1      = (*mInfo.mEit1)->he1();
    auto body1_vit1     = (*body1_he1)->src();
    auto body1_vit2     = (*body1_he1)->dst();
    auto body1_p1       = (*body1_vit1)->pGCS(rotMat1, com1);
    auto body1_p2       = (*body1_vit2)->pGCS(rotMat1, com1);
    auto body2_he1      = (*mInfo.mEit2)->he1();
    auto body2_vit1     = (*body2_he1)->src();
    auto body2_vit2     = (*body2_he1)->dst();
    auto body2_p1       = (*body2_vit1)->pGCS(rotMat2, com2);
    auto body2_p2       = (*body2_vit2)->pGCS(rotMat2, com2);

    double z11 = body1_p1.z();
    double z12 = body1_p2.z();
    if (z11 > z12) {
        swap(z11, z12);
        swap(body1_p1, body1_p2);
        swap(body1_vit1, body1_vit2);
    }
    double z21 = body2_p1.z();
    double z22 = body2_p2.z();
    if (z21 > z22) {
        swap(z21, z22);
        swap(body2_p1, body2_p2);
        swap(body2_vit1, body2_vit2);
    }

    // Classify the edge placements.
    if (z12 - mEpsilonZero <= z21) {
        // 11 12 21 22
        //  o--o  *--* 
        //     ^  ^
        mInfo.mIntsect1.push_back((*body1_vit2)->pLCS());
        mInfo.mIntsect2.push_back((*body2_vit1)->pLCS());
    }
    else if (z22 - mEpsilonZero <= z11) {
        // 21 22 11 12
        //  *--*  o--o
        //     ^  ^
        mInfo.mIntsect1.push_back((*body1_vit1)->pLCS());
        mInfo.mIntsect2.push_back((*body2_vit2)->pLCS());
    }
    else if (z11<=z21) {
        if (z22<=z12) {
            // 11 21 22 12 
            //  o--*--*--o
            //     ^  ^
            const auto v1 = interpolate(body1_p1, body1_p2, body2_p1.z());
            const auto v2 = interpolate(body1_p1, body1_p2, body2_p2.z());
            Mat3x3 Qmat1inv = rotMat1.transpose();
            mInfo.mIntsect1.push_back(Qmat1inv * (v1 - com1));
            mInfo.mIntsect1.push_back(Qmat1inv * (v2 - com1));
            mInfo.mIntsect2.push_back((*body2_vit1)->pLCS());
            mInfo.mIntsect2.push_back((*body2_vit2)->pLCS());
        }
        else {
            // 11 21 12 22
            //  o--*--o--*
            //     ^  ^
            const auto v1 = interpolate(body1_p1, body1_p2, body2_p1.z());
            const auto v2 = interpolate(body2_p1, body2_p2, body1_p2.z());
            Mat3x3 Qmat1inv = rotMat1.transpose();
            mInfo.mIntsect1.push_back(Qmat1inv * (v1 - com1));
            mInfo.mIntsect1.push_back((*body1_vit2)->pLCS());
            Mat3x3 Qmat2inv = rotMat2.transpose();
            mInfo.mIntsect2.push_back((*body2_vit1)->pLCS());
            mInfo.mIntsect2.push_back(Qmat2inv * (v2 - com2));
        }
    }
    else {
        if (z12<=z22) {
            // 21 11 12 22
            //  *--o--o--*
            //     ^  ^
            const auto v1 = interpolate(body2_p1, body2_p2, body1_p1.z());
            const auto v2 = interpolate(body2_p1, body2_p2, body1_p2.z());

            mInfo.mIntsect1.push_back((*body1_vit1)->pLCS());
            mInfo.mIntsect1.push_back((*body1_vit2)->pLCS());
            Mat3x3 Qmat2inv = rotMat2.transpose();
            mInfo.mIntsect2.push_back(Qmat2inv * (v1 - com2));
            mInfo.mIntsect2.push_back(Qmat2inv * (v2 - com2));
        }
        else {
            // 21 11 22 12
            //  *--o--*--o
            //     ^  ^
            const auto v1 = interpolate(body1_p1, body1_p2, body2_p2.z());
            const auto v2 = interpolate(body2_p1, body2_p2, body1_p1.z());

            Mat3x3 Qmat1inv = rotMat1.transpose();
            mInfo.mIntsect1.push_back((*body1_vit1)->pLCS());
            mInfo.mIntsect1.push_back(Qmat1inv * (v1 - com1));

            Mat3x3 Qmat2inv = rotMat2.transpose();
            mInfo.mIntsect2.push_back(Qmat2inv * (v2 - com2));
            mInfo.mIntsect2.push_back((*body2_vit2)->pLCS());
        }
    }

    mInfo.mContactNormal1To2 = (*mInfo.mEit1)->nGCS(mQmat1) -
                               (*mInfo.mEit2)->nGCS(mQmat2) ;
    mInfo.mContactNormal1To2.normalize();
    mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;


    mInfo.mIntsect1reduced = mInfo.mIntsect1;
    mInfo.mIntsect2reduced = mInfo.mIntsect2;


}


void ContactPointsAndNormalGenerator::process_EDGE_VERTEX()
{
    auto he    = (*mInfo.mEit1)->he1();
    auto pSrc  = (*((*he)->src()))->pGCS(mQmat1, mCoM1);
    auto pDst  = (*((*he)->dst()))->pGCS(mQmat1, mCoM1);
    auto pTest = (*(mInfo.mVit2) )->pGCS(mQmat2, mCoM2);

    auto pProj = findClosestPointOnEdge(pSrc, pDst, pTest);

    Mat3x3 Qmat1inv = mQmat1.transpose();

    mInfo.mIntsect1.push_back(Qmat1inv * (pProj - mCoM1));
    mInfo.mIntsect2.push_back((*mInfo.mVit2)->pLCS());

    mInfo.mContactNormal1To2 = (*mInfo.mEit1)->nGCS(mQmat1);
    mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_NORMAL_1;

    mInfo.mIntsect1reduced = mInfo.mIntsect1;
    mInfo.mIntsect2reduced = mInfo.mIntsect2;

}


void ContactPointsAndNormalGenerator::process_VERTEX_FACE()
{
    // Project the vertex on Body 1 onto the plane of the face on Body 2.

    // Body2's normal
    Vec3 normalDir2 = (*mInfo.mFit2)->nGCS(mQmat2);

    Mat3x3 rotMat1,  rotMat2;
    Vec3   com1,     com2;
    rotateWorld(normalDir2, rotMat1, com1, rotMat2, com2);

    auto he2 = *((*mInfo.mFit2)->halfEdges()).begin();
    auto rotatedPoint2 = ((*((*he2)->src()))->pGCS(rotMat2, com2));

    auto rotatedPoint1 = (*(mInfo.mVit1))->pGCS(rotMat1, com1);
    rotatedPoint1.setZ(rotatedPoint2.z());

    mInfo.mIntsect1.push_back((*(mInfo.mVit1))->pLCS());
    Mat3x3 Qmat2inv = rotMat2.transpose();
    mInfo.mIntsect2.push_back(Qmat2inv * (rotatedPoint1 - com2));

    mInfo.mContactNormal1To2 = normalDir2 * -1.0;
    mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;

    mInfo.mIntsect1reduced = mInfo.mIntsect1;
    mInfo.mIntsect2reduced = mInfo.mIntsect2;
}


void ContactPointsAndNormalGenerator::process_VERTEX_EDGE()
{
    auto he    = (*mInfo.mEit2)->he1();
    auto pSrc  = (*((*he)->src()))->pGCS(mQmat2, mCoM2);
    auto pDst  = (*((*he)->dst()))->pGCS(mQmat2, mCoM2);
    auto pTest = (*(mInfo.mVit1) )->pGCS(mQmat1, mCoM1);
    auto pProj = findClosestPointOnEdge(pSrc, pDst, pTest);
    Mat3x3 Qmat2inv = mQmat2.transpose();

    mInfo.mIntsect1.push_back((*mInfo.mVit1)->pLCS());
    mInfo.mIntsect2.push_back(Qmat2inv * (pProj - mCoM2));

    mInfo.mContactNormal1To2 = (*mInfo.mEit2)->nGCS(mQmat2) * -1.0;
    mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_NORMAL_2;

    mInfo.mIntsect1reduced = mInfo.mIntsect1;
    mInfo.mIntsect2reduced = mInfo.mIntsect2;
}


void ContactPointsAndNormalGenerator::process_VERTEX_VERTEX()
{
    mInfo.mIntsect1.push_back((*mInfo.mVit1)->pLCS());
    mInfo.mIntsect2.push_back((*mInfo.mVit2)->pLCS());

    mInfo.mContactNormal1To2 = (*mInfo.mVit1)->nGCS(mQmat1) - 
                               (*mInfo.mVit2)->nGCS(mQmat2) ;
    mInfo.mContactNormal1To2.normalize();
    mInfo.mContactNormalType = ContactPairInfo::NT_VERTEX_VERTEX_AVG;


    mInfo.mIntsect1reduced = mInfo.mIntsect1;
    mInfo.mIntsect2reduced = mInfo.mIntsect2;
}


vector<bool> ContactPointsAndNormalGenerator::reduceContactPoints(
    vector<IntersectionFinderConvexPolygon2D::OutputElem>& elems2D
) {
    vector<bool> flags;
    auto numElems = elems2D.size();
    if (numElems <= mMaxNumPointsPerContact) {
        for (long i = 0; i < numElems; i++) {        
            flags.push_back(true);
        }
    }
    else {
        // PCA and find 3 or 4 extreme points.
        vector<Vec2> points2D;
        for (long i = 0; i < numElems; i++) {
            points2D.push_back(elems2D[i].mP);
            flags.push_back(false);
        }

        Vec2 spread, mean, axis1, dummy;
        findPrincipalComponents(points2D, spread, mean, axis1, dummy);
        long   minIndex1, maxIndex1;
        double minDot,   maxDot;
        for (long i = 0; i < numElems; i++) {
            double dot = points2D[i].dot(axis1);
            if (i==0) {
                minIndex1 = 0;
                maxIndex1 = 0;
                minDot    = dot;
                maxDot    = dot;
            }
            else {
                if (minDot > dot) {
                    minIndex1 = i;
                    minDot    = dot;
                }
                if (maxDot < dot) {
                    maxIndex1 = i;
                    maxDot    = dot;
                }
            }
        }

        Vec2 axis2 = (points2D[maxIndex1] - points2D[minIndex1]);
        axis2.normalize();
        axis2 = axis2.perp();
        long minIndex2, maxIndex2;
        for (long i = 0; i < numElems; i++) {
            double dot = points2D[i].dot(axis2);
            if (i==0) {
                minIndex2 = 0;
                maxIndex2 = 0;
                minDot    = dot;
                maxDot    = dot;
            }
            else {
                if (minDot > dot) {
                    minIndex2 = i;
                    minDot    = dot;
                }
                if (maxDot < dot) {
                    maxIndex2 = i;
                    maxDot    = dot;
                }
            }
        }
        flags[minIndex1] = true;
        flags[maxIndex1] = true;
        flags[minIndex2] = true;
        flags[maxIndex2] = true;

        long assigned = 0;
        vector<double> angles1, angles2;
        for (long i = 0; i < numElems; i++) {

            if (flags[i]) {
                assigned++;
                angles1.push_back(-1.0);
                angles2.push_back(-1.0);
            }
            else {
                long iPrev = (i+numElems-1)%numElems;
                long iNext = (i+1)%numElems;

                Vec2 vPrev = points2D[i] -     points2D[iPrev];
                vPrev.normalize();
                Vec2 vNext = points2D[iNext] - points2D[i];
                vNext.normalize();
                double fdot = fabs(vPrev.dot(vNext.perp()));
                angles1.push_back(fdot);
                angles2.push_back(fdot);
            }
        }
        std::sort(angles1.begin(), angles1.end());

        auto numToBeFilled = mMaxNumPointsPerContact - assigned;

        if (numToBeFilled > 0) {
            auto threshold = angles1[numElems - numToBeFilled];
            for (long i = 0; i < numElems && numToBeFilled > 0; i++) {
                if (angles2[i]>=threshold) {
                    flags[i] = true;
                    numToBeFilled--;
                }
            }
        }
    }
    return flags;
}


void ContactPointsAndNormalGenerator::from2DXYto3DPoints(
    vector<IntersectionFinderConvexPolygon2D::OutputElem>&
                   points2D,                         
    Mat3x3&        rotMat1,
    Vec3&          com1,
    Mat3x3&        rotMat2,
    Vec3&          com2
) {
    Vec3 pBase, N;
    auto numElems = points2D.size();
    vector<bool> reducedPointsFlags = reduceContactPoints(points2D);

    if (mInfo.mType1 == ContactPairInfo::FT_FACE) {
        auto he = *((*mInfo.mFit1)->halfEdges().begin());
        pBase   = (*((*he)->src()))->pGCS(rotMat1, com1);
        N       = (*mInfo.mFit1)->nGCS(rotMat1);
    }
    else if (mInfo.mType1 == ContactPairInfo::FT_EDGE) {
        auto he = (*mInfo.mEit1)->he1();
        pBase   = (*((*he)->src()))->pGCS(rotMat1, com1);
        N       = (*mInfo.mEit1)->nGCS(rotMat1);
    }
    else {
        auto vit = mInfo.mVit1;
        pBase   = (*vit)->pGCS(rotMat1, com1);
        N       = (*vit)->nGCS(rotMat1);
    }


    vector<Vec3> projection3D1;
    orthoProj(points2D, N, pBase, projection3D1);

    // TODO:
    if (projection3D1.size()==0) {
        return;
    }
    Mat3x3 Qmat1inv = rotMat1.transpose();

    for (long i = 0; i < numElems; i++) {
        auto& pProj = projection3D1[i];
        auto  p = Qmat1inv * (pProj - com1);
        mInfo.mIntsect1.push_back(p);
        if (reducedPointsFlags[i]) {
            mInfo.mIntsect1reduced.push_back(p);
        }
    }

    if (mInfo.mType2 == ContactPairInfo::FT_FACE) {
        auto he = *((*mInfo.mFit2)->halfEdges().begin());
        pBase   = (*((*he)->src()))->pGCS(rotMat2, com2);
        N       = (*mInfo.mFit2)->nGCS(rotMat2);
    }
    else if (mInfo.mType2 == ContactPairInfo::FT_EDGE) {
        auto he = (*mInfo.mEit2)->he1();
        pBase   = (*((*he)->src()))->pGCS(rotMat2, com2);
        N       = (*mInfo.mEit2)->nGCS(rotMat2);
    }
    else {
        auto vit = mInfo.mVit2;
        pBase   = (*vit)->pGCS(rotMat2, com2);
        N       = (*vit)->nGCS(rotMat2);
    }


    vector<Vec3> projection3D2;
    orthoProj(points2D, N, pBase, projection3D2);
    // TODO:
    if (projection3D2.size()==0) {;
        mInfo.mIntsect1.clear();
        mInfo.mIntsect1reduced.clear();
        return;
    }


    Mat3x3 Qmat2inv = rotMat2.transpose();

    for (long i = 0; i < numElems; i++) {
        auto& pProj = projection3D2[i];
        auto  p = Qmat2inv * (pProj - com2);
        mInfo.mIntsect2.push_back(p);
        if (reducedPointsFlags[i]) {
            mInfo.mIntsect2reduced.push_back(p);
        }
    }
}


void ContactPointsAndNormalGenerator::orthoProj(
   vector<IntersectionFinderConvexPolygon2D::OutputElem>&
                 points2D,                         
   const Vec3&   N,
   const Vec3&   pBase,
   vector<Vec3>& points3D
) {
    if (fabs(N.z()) <= mEpsilonAngle) {
        log(ERROR, __FILE__, __LINE__,"Contact Normal illconditioned");
        // TODO: TEST
        return;
        double z = pBase.z();
        for (auto p2 : points2D) {
            points3D.push_back(Vec3(p2.mP.x(), p2.mP.y(), z));
        }
    }
    else {
        for (auto p2 : points2D) {
            double z = -1.0*((p2.mP.x() - pBase.x())*N.x() + 
                             (p2.mP.y() - pBase.y())*N.y()   ) / N.z() +
                       pBase.z();
            points3D.push_back(Vec3(p2.mP.x(), p2.mP.y(), z));
        }
    }
}


Vec3 ContactPointsAndNormalGenerator::findProjectedPointOnPlane(
    const Vec3& pBase,
    const Vec3& N,
    const Vec3& pTest
) {
    double s = N.dot(pTest - pBase);
    return pTest - (N * s);
}


Vec3 ContactPointsAndNormalGenerator::findClosestPointOnEdge(
    const Vec3&      p1,
    const Vec3&      p2,
    const Vec3&      b
) {
    Vec3  v12     = p2 - p1;
    Vec3  v1b     = b  - p1;
    double sqDist = v12.squaredNorm2();
    double t      = v1b.dot(v12) / sqDist;

    if (t >= 1.0 - mEpsilonZero) {
        return p2;
    }
    else if (t <= mEpsilonZero) {
      return p1;
    }
    else {
        return p1 + (v12 * t);
    }
}


void ContactPointsAndNormalGenerator::findClosestPointsBetweenTwoEdges(
    const Vec3& a1,
    const Vec3& a2,
    const Vec3& b1,
    const Vec3& b2,
    Vec3&       pA,
    Vec3&       pB
) {
    Vec3 va = a2 - a1;
    Vec3 vb = b2 - b1;

    Vec3 va_cr_vb = va.cross(vb);
    if (va_cr_vb.squaredNorm2() <= mEpsilonZero) {

        // v1 and v2 are parallel.
        Vec3   pa1 = findClosestPointOnEdge(b1, b2, a1);
        Vec3   pa2 = findClosestPointOnEdge(b1, b2, a2);
        Vec3   pb1 = findClosestPointOnEdge(a1, a2, b1);
        Vec3   pb2 = findClosestPointOnEdge(a1, a2, b2);

        Vec3   va1 = pa1 - a1;
        Vec3   va2 = pa2 - a2;
        Vec3   vb1 = b1  - pb1;
        Vec3   vb2 = b2  - pb2;

        double da1 = va1.squaredNorm2();
        double da2 = va2.squaredNorm2();
        double db1 = vb1.squaredNorm2();
        double db2 = vb2.squaredNorm2();

        double d = da1;
        pA = a1;
        pB = pa1;
        if (d > da2) {
            pA = a2; pB = pa2; d = da2;
        }
        if (d > db1) {
            pA = pb1; pB = b1; d = db1;
        }
        if (d > db2) {
            pA = pb2; pB = b2;
        }
        return;
    }

    Vec3   vab        = b1 - a1;
    double va_2       = va.squaredNorm2();
    double vb_2       = vb.squaredNorm2();
    double va_dot_vb  = va.dot(vb);
    double vab_dot_va = vab.dot(va);
    double vab_dot_vb = vab.dot(vb);

    double s  = (vab_dot_va * vb_2      - vab_dot_vb * va_dot_vb ) / 
                (va_2 * vb_2            - va_dot_vb * va_dot_vb  );

    double t  = (vab_dot_va * va_dot_vb - vab_dot_vb * va_2      ) /
                (va_2 * vb_2            - va_dot_vb * va_dot_vb  );

    s = std::min(s, 1.0);
    s = std::max(s, 0.0);
    t = std::min(t, 1.0);
    t = std::max(t, 0.0);

    pA = a1 + (va * s);
    pB = b1 + (vb * t);

}


void ContactPointsAndNormalGenerator::rotateWorld(
    Mat3x3&     rotMat1,
    Vec3&       com1,
    Mat3x3&     rotMat2,
    Vec3&       com2
) {
    rotateWorld(mInfo.mContactNormal1To2, rotMat1, com1, rotMat2, com2);
}


void ContactPointsAndNormalGenerator::rotateWorld(
    const Vec3& zDir,
    Mat3x3&     rotMat1,
    Vec3&       com1,
    Mat3x3&     rotMat2,
    Vec3&       com2
) {

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

    rotMat1 = rMat * mQmat1;
    com1    = rMat * mCoM1;
    rotMat2 = rMat * mQmat2;
    com2    = rMat * mCoM2;
}


void ContactPointsAndNormalGenerator::findContactNormal()
{
    if (mInfo.mType1==ContactPairInfo::FT_FACE) {
        // Use the face normal.
        mInfo.mContactNormal1To2 =(*mInfo.mFit1)->nGCS(mQmat1);
        mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_1;
    }
    else if (mInfo.mType2==ContactPairInfo::FT_FACE) {
        // Use the face normal.
        mInfo.mContactNormal1To2 = (*mInfo.mFit2)->nGCS(mQmat2) * -1.0;
        mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;
    }
    else if (mInfo.mType1==ContactPairInfo::FT_EDGE) {

        if (mInfo.mType2==ContactPairInfo::FT_EDGE) {

            // Use the cross product of edges if two edges are parallel, 
            // use the avg of edge normals.
            auto vit1src = (*((*mInfo.mEit1)->he1()))->src();
            auto vit1dst = (*((*mInfo.mEit1)->he1()))->dst();
            auto vit2src = (*((*mInfo.mEit2)->he1()))->src();
            auto vit2dst = (*((*mInfo.mEit2)->he1()))->dst();

            auto d1      = (*vit1dst)->pGCS(mQmat1) - (*vit1src)->pGCS(mQmat1);
            auto d2      = (*vit2dst)->pGCS(mQmat2) - (*vit2src)->pGCS(mQmat2);
            d1.normalize();
            d2.normalize();
            auto cr = d1.cross(d2);

            auto body1_en       = (*mInfo.mEit1)->nGCS(mBody1.Qmat());
            auto body2_en       = (*mInfo.mEit2)->nGCS(mBody2.Qmat());
            auto enAvg          = body1_en - body2_en;
            enAvg.normalize();

            if (cr.squaredNorm2() <= mEpsilonZero) {
                // two edges are parallel. Take average of the two edge normals
                mInfo.mContactNormal1To2 = enAvg;
                mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;
            }
            else {
                cr.normalize();
                auto dot = cr.dot(enAvg);
                if ( fabs(dot) >= sqrt(2.0)/2.0 ) {
                    if (cr.dot(mBody2.CoM() - mBody1.CoM())< 0.0) {
                        cr.scale(-1.0);
                    }
                    // the angle between the cross and the edge normal<=pi/2.0
                    mInfo.mContactNormal1To2 = cr;
                    mInfo.mContactNormalType = 
                                          ContactPairInfo::NT_EDGE_CROSS_EDGE;
                }
                else {
                    mInfo.mContactNormal1To2 = enAvg;
                    mInfo.mContactNormalType = 
                                          ContactPairInfo::NT_EDGE_EDGE_AVG;
                }
            }
        }
        else {
            // Use the edge normal.
            mInfo.mContactNormal1To2 = (*mInfo.mEit1)->nGCS(mQmat1);
            mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_NORMAL_1;
        }
    }
    else {

        if (mInfo.mType2==ContactPairInfo::FT_EDGE) {
            // Use the edge normal.
            mInfo.mContactNormal1To2 = (*mInfo.mEit2)->nGCS(mQmat2) * -1.0;
            mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_NORMAL_2;
        }
        else {
            // Use the avg of vertex normals
            mInfo.mContactNormal1To2 = (*mInfo.mVit1)->nGCS(mQmat1) -
                                       (*mInfo.mVit2)->nGCS(mQmat2) ;
            mInfo.mContactNormal1To2.normalize();
            mInfo.mContactNormalType = ContactPairInfo::NT_VERTEX_VERTEX_AVG;
        }
    }

    mInfo.mContactNormal1To2.normalize();
}


void ContactPointsAndNormalGenerator::adjustFeatures()
{
    if ( mInfo.mType1==ContactPairInfo::FT_FACE ) {
        if ( mInfo.mType2==ContactPairInfo::FT_FACE ) {
            ;
        }
        else if ( mInfo.mType2==ContactPairInfo::FT_EDGE ) {
            adjustFeatures_FACE_EDGE();
        }
        else {
            adjustFeatures_FACE_VERTEX();
        }
    }
    else if ( mInfo.mType1==ContactPairInfo::FT_EDGE ) {
        if ( mInfo.mType2==ContactPairInfo::FT_FACE ) {
            adjustFeatures_EDGE_FACE();    
        }
        else if ( mInfo.mType2==ContactPairInfo::FT_EDGE ) {
            adjustFeatures_EDGE_EDGE();
        }
        else {
            adjustFeatures_EDGE_VERTEX();
        }
    }
    else {
        if ( mInfo.mType2==ContactPairInfo::FT_FACE ) {
            adjustFeatures_VERTEX_FACE();
        }
        else if ( mInfo.mType2==ContactPairInfo::FT_EDGE ) {
            adjustFeatures_VERTEX_EDGE();
        }
        else {
            adjustFeatures_VERTEX_VERTEX();
        }
    }
}


void ContactPointsAndNormalGenerator::adjustFeatures_FACE_EDGE()
{
    auto body2_he1  = (*mInfo.mEit2)->he1();
    auto body2_he2  = (*mInfo.mEit2)->he2();
    auto body2_vit1 = (*body2_he1)->src();
    auto body2_vit2 = (*body2_he1)->dst();
    const auto body2_p1 = (*body2_vit1)->pGCS(mQmat2, mCoM2);
    const auto body2_p2 = (*body2_vit2)->pGCS(mQmat2, mCoM2);

    auto body2_v12  = body2_p2 - body2_p1;
    body2_v12.normalize();
    auto body2_fit1 = (*body2_he1)->face();
    auto body2_fit2 = (*body2_he2)->face();

    const auto body2_n1 = (*body2_fit1)->nGCS(mQmat2);
    const auto body2_n2 = (*body2_fit2)->nGCS(mQmat2);

    const auto body2_d1 = body2_n1.cross(body2_v12); // from v12 inward face 1
    const auto body2_d2 = body2_v12.cross(body2_n2); // from v12 inward face 2

    const auto body1_n  = (*mInfo.mFit1)->nGCS(mQmat1);

    const auto body2_d1_dot_body1_n = body2_d1.dot(body1_n);
    const auto body2_d2_dot_body1_n = body2_d2.dot(body1_n);

    if (body2_d1_dot_body1_n >= mEpsilonAngle && 
        body2_d2_dot_body1_n >= mEpsilonAngle   ) {
        // test face-edge
        return;
    }
    else if (body2_d1_dot_body1_n < body2_d2_dot_body1_n) {
        // test face-face 1
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2  = body2_fit1;
    }
    else {
        // test face-face 2
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2  = body2_fit2;
    }
}


void ContactPointsAndNormalGenerator::adjustFeatures_EDGE_FACE()
{
    auto body1_he1  = (*mInfo.mEit1)->he1();
    auto body1_he2  = (*mInfo.mEit1)->he2();
    auto body1_vit1 = (*body1_he1)->src();
    auto body1_vit2 = (*body1_he1)->dst();

    const auto body1_p1 = (*body1_vit1)->pGCS(mQmat1, mCoM1);
    const auto body1_p2 = (*body1_vit2)->pGCS(mQmat1, mCoM1);

    auto body1_v12  = body1_p2 - body1_p1;
    body1_v12.normalize();
    auto body1_fit1 = (*body1_he1)->face();
    auto body1_fit2 = (*body1_he2)->face();

    const auto body1_n1 = (*body1_fit1)->nGCS(mQmat1);
    const auto body1_n2 = (*body1_fit2)->nGCS(mQmat1);

    const auto body1_d1 = body1_n1.cross(body1_v12); // from v12 inward face 1
    const auto body1_d2 = body1_v12.cross(body1_n2); // from v12 inward face 2

    const auto body2_n  = (*mInfo.mFit2)->nGCS(mQmat2);

    const auto body1_d1_dot_body2_n = body1_d1.dot(body2_n);
    const auto body1_d2_dot_body2_n = body1_d2.dot(body2_n);

    if (body1_d1_dot_body2_n >= mEpsilonAngle && 
        body1_d2_dot_body2_n >= mEpsilonAngle   ) {
        // test edge-face
        return;
    }
    else if (body1_d1_dot_body2_n < body1_d2_dot_body2_n) {
        // test face1-face
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1  = body1_fit1;
    }
    else {
        // test face2-face
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1  = body1_fit2;
    }
}


void ContactPointsAndNormalGenerator::adjustFeatures_FACE_VERTEX()
{
    const auto body1_n = (*mInfo.mFit1)->nGCS(mQmat1);

    const auto body2_p = (*mInfo.mVit2)->pGCS(mQmat2, mCoM2);

    //           *previt
    //          /
    //         /
    //        v
    //       *-------------->* nextit
    //   body2_p   body2_he
    //
    //
    long penetratedEdgeCount = 0;
    HalfEdgeIt minItHalfEdge;
    double body1_n_dot_prevd;
    double body1_n_dot_firstd;
    bool first = true;
    VertexIt lastit;
    for (auto body2_he : (*mInfo.mVit2)->halfEdges()) {

        if ((*body2_he)->src()==mInfo.mVit2) {
            if (first) {
                VertexIt previt =(*((*((*body2_he)->buddy()))->next()))->dst();
                auto prevd = (*previt)->pGCS(mQmat2, mCoM2) - body2_p;
                prevd.normalize();
                body1_n_dot_prevd = body1_n.dot(prevd);
                lastit = previt;
            }

            VertexIt nextit = (*body2_he)->dst();
            auto nextd = (*nextit)->pGCS(mQmat2, mCoM2) - body2_p;
            nextd.normalize();

            const auto body1_n_dot_nextd = body1_n.dot(nextd);
 
            if (body1_n_dot_nextd <= mEpsilonAngle) {

                if (first) {
                    body1_n_dot_firstd = body1_n_dot_nextd;
                }

                if (penetratedEdgeCount == 0) {
                    minItHalfEdge = body2_he;
                }
                else if ( penetratedEdgeCount == 1             && 
                          body1_n_dot_prevd   >  mEpsilonAngle    ) {
                    if (nextit==lastit&&body1_n_dot_firstd <= mEpsilonAngle) {
                        minItHalfEdge = body2_he;                        
                    }
                    else {
                        break;
                    }
                }
                else if (penetratedEdgeCount>1) {
                    break;
                }
                penetratedEdgeCount++;
            }

            body1_n_dot_prevd = body1_n_dot_nextd;
            first = false;
        }
    }

    if (penetratedEdgeCount==1) { 
        // Edge penetration.
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = (*minItHalfEdge)->edge();
    }
    else if (penetratedEdgeCount==2) { 
        // Face penetration.
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2  = (*minItHalfEdge)->face();
    }
}


void ContactPointsAndNormalGenerator::adjustFeatures_VERTEX_FACE()
{
    const auto body2_n = (*mInfo.mFit2)->nGCS(mQmat2);

    const auto body1_p = (*mInfo.mVit1)->pGCS(mQmat1, mCoM1);

    //           *previt
    //          /
    //         /
    //        v
    //       *-------------->* nextit
    //   body1_p   body1_he
    //
    //
    long penetratedEdgeCount = 0;
    HalfEdgeIt minItHalfEdge;
    double body2_n_dot_prevd;
    double body2_n_dot_firstd;
    bool first = true;
    VertexIt lastit;
    for (auto body1_he : (*mInfo.mVit1)->halfEdges()) {

        if ((*body1_he)->src()==mInfo.mVit1) {
            if (first) {
                VertexIt previt =(*((*((*body1_he)->buddy()))->next()))->dst();
                auto prevd = (*previt)->pGCS(mQmat1, mCoM1) - body1_p;
                prevd.normalize();
                body2_n_dot_prevd = body2_n.dot(prevd);
                lastit = previt;
            }

            VertexIt nextit = (*body1_he)->dst();
            auto nextd = (*nextit)->pGCS(mQmat1, mCoM1) - body1_p;
            nextd.normalize();

            const auto body2_n_dot_nextd = body2_n.dot(nextd);
 
            if (body2_n_dot_nextd <= mEpsilonAngle) {

                if (first) {
                    body2_n_dot_firstd = body2_n_dot_nextd;
                }

                if (penetratedEdgeCount == 0) {
                    minItHalfEdge = body1_he;
                }
                else if ( penetratedEdgeCount == 1             && 
                          body2_n_dot_prevd   >  mEpsilonAngle    ) {
                    if (nextit==lastit&&body2_n_dot_firstd <= mEpsilonAngle) {
                        minItHalfEdge = body1_he;                        
                    }
                    else {
                        break;
                    }
                }
                else if (penetratedEdgeCount>1) {
                    break;
                }
                penetratedEdgeCount++;
            }

            body2_n_dot_prevd = body2_n_dot_nextd;
            first = false;
        }
    }

    if (penetratedEdgeCount==1) { 
        // Edge penetration.
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*minItHalfEdge)->edge();
    }
    else if (penetratedEdgeCount==2) { 
        // Face penetration.
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1  = (*minItHalfEdge)->face();
    }
}


void ContactPointsAndNormalGenerator::adjustFeatures_EDGE_EDGE()
{
    auto body1_he1      = (*mInfo.mEit1)->he1();
    auto body1_vit1     = (*body1_he1)->src();
    auto body1_vit2     = (*body1_he1)->dst();
    const auto body1_p1 = (*body1_vit1)->pGCS(mBody1.Qmat(), mBody1.CoM());
    const auto body1_p2 = (*body1_vit2)->pGCS(mBody1.Qmat(), mBody1.CoM());
    auto body1_v12  = body1_p2 - body1_p1;
    body1_v12.normalize();

    auto body2_he1      = (*mInfo.mEit2)->he1();
    auto body2_vit1     = (*body2_he1)->src();
    auto body2_vit2     = (*body2_he1)->dst();
    const auto body2_p1 = (*body2_vit1)->pGCS(mBody2.Qmat(), mBody2.CoM());
    const auto body2_p2 = (*body2_vit2)->pGCS(mBody2.Qmat(), mBody2.CoM());
    auto body2_v12  = body2_p2 - body2_p1;
    body2_v12.normalize();
    auto cr = body1_v12.cross(body2_v12);

    auto body1_en       = (*mInfo.mEit1)->nGCS(mBody1.Qmat());
    auto body2_en       = (*mInfo.mEit2)->nGCS(mBody2.Qmat());
    auto enAvg          = body1_en - body2_en;
    enAvg.normalize();
    if (cr.squaredNorm2() <= mEpsilonAngle) {
        adjustFeatures_EDGE_EDGE_PARALLEL();
    }
    else {
        cr.normalize();
        if (fabs(enAvg.dot(cr)) <= sqrt(2.0)/2.0) {
            adjustFeatures_EDGE_EDGE_PARALLEL();
        }
        else {
            adjustFeatures_EDGE_EDGE_CROSSING();
        }
    }
}


void ContactPointsAndNormalGenerator::adjustFeatures_EDGE_EDGE_CROSSING()
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

    if ( pen_body1_d1 <= mEpsilonAngle ) {
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1 = body1_fit1;
    }
    else if ( pen_body1_d2 <= mEpsilonAngle ) {    
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1 = body1_fit2;
    }

    if ( pen_body2_d1 <= mEpsilonAngle ) {
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2 = body2_fit1;
    }
    else if ( pen_body2_d2 <= mEpsilonAngle ) {    
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2 = body2_fit2;
    }
}


void ContactPointsAndNormalGenerator::adjustFeatures_EDGE_EDGE_PARALLEL()
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

    if (body1_v12.dot(body2_v12) < 0.0) {
        if ( body1_d1.dot(body2_d1)       >  0.0 && 
             fabs(body1_n1.dot(body2_d1)) <= mEpsilonAngle ) {
            //         .....d12
            //         ...../
            //  d11/d21<===*
            //         .....\
            //         .....d22
            mInfo.mType1 = ContactPairInfo::FT_FACE;
            mInfo.mType2 = ContactPairInfo::FT_FACE;
            mInfo.mFit1  = body1_fit1;
            mInfo.mFit2  = body2_fit1;
        }
        else if ( body1_d2.dot(body2_d2)       >  0.0 && 
                  fabs(body2_n2.dot(body1_d2)) <= mEpsilonAngle ) {
            //         .....d21
            //         ...../
            //  d12/d22<===*
            //         .....\
            //         .....d11
            mInfo.mType1 = ContactPairInfo::FT_FACE;
            mInfo.mType2 = ContactPairInfo::FT_FACE;
            mInfo.mFit1  = body1_fit2;
            mInfo.mFit2  = body2_fit2;
        }
        else if ( body1_d1.dot(body2_d2)       >  0.0 && 
                  fabs(body2_n2.dot(body1_d1)) <= mEpsilonAngle ) {
            if (isCCW(body1_d1, body1_d2, body2_d1, body1_v12)) {
                //         .d21.d12 
                //         ...\./
                // d11/d22 <===*
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit1;
                mInfo.mFit2  = body2_fit1;
            }
            else {
                //         .d12.d21
                //         ...\./
                // d11/d22 <===*
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit2;
                mInfo.mFit2  = body2_fit2;
            }
        }
        else if ( body1_d2.dot(body2_d1)       >  0.0 && 
                  fabs(body2_n1.dot(body1_d2)) <= mEpsilonAngle ) {
            if (isCCW(body1_d1, body1_d2, body2_d2, body1_v12)) {
                // d12/d21 <===*
                //         .../.\
                //         .d22.d11
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit2;
                mInfo.mFit2  = body2_fit2;
            }
            else {
                // d12/d21 <===*
                //         .../.\
                //         .d11.d22
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit1;
                mInfo.mFit2  = body2_fit1;
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
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit2;
                    mInfo.mFit2  = body2_fit2;
                }
                else {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit1;
                    mInfo.mFit2  = body2_fit1;
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
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit2;
                mInfo.mFit2  = body2_fit2;
            }
            else if (  ccw_11_12_21 && !ccw_11_12_22 && 
                      !ccw_11_21_22 &&  ccw_12_21_22    ) {
                // d11->d22->d12->d21  (d11 is penetrating into d21)
                //
                //            _-*-_
                //          _- / \ -_
                //      d12-  /   \  -d22
                //          d21   d11
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit1;
                mInfo.mFit2  = body2_fit1;
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
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit1;
                    mInfo.mFit2  = body2_fit1;
                }
                else {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit2;
                    mInfo.mFit2  = body2_fit2;
                }
            }
            else {
                // d11->d12->d22->d21  (Forbidden)
                log(ERROR, __FILE__, __LINE__, "Forbidden predicates");
            }
        }
    }
    else {

        if ( body1_d1.dot(body2_d2)       >  0.0 && 
             fabs(body1_n1.dot(body2_d2)) <= mEpsilonAngle ) {
            //         .....d12
            //         ...../
            //  d11/d22<===*
            //         .....\
            //         .....d21
            mInfo.mType1 = ContactPairInfo::FT_FACE;
            mInfo.mType2 = ContactPairInfo::FT_FACE;
            mInfo.mFit1  = body1_fit1;
            mInfo.mFit2  = body2_fit2;
        }
        else if ( body1_d2.dot(body2_d1)       >  0.0 && 
                  fabs(body2_n1.dot(body1_d2)) <= mEpsilonAngle ) {
            //         .....d22
            //         ...../
            //  d12/d21<===*
            //         .....\
            //         .....d11
            mInfo.mType1 = ContactPairInfo::FT_FACE;
            mInfo.mType2 = ContactPairInfo::FT_FACE;
            mInfo.mFit1  = body1_fit2;
            mInfo.mFit2  = body2_fit1;
        }
        else if ( body1_d1.dot(body2_d1)       >  0.0 && 
                  fabs(body2_n1.dot(body1_d1)) <= mEpsilonAngle ) {
            if (isCCW(body1_d1, body1_d2, body2_d2, body1_v12)) {
                //         .d22.d12 
                //         ...\./
                // d11/d21 <===*
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit1;
                mInfo.mFit2  = body2_fit2;
            }
            else {
                //         .d12.d22
                //         ...\./
                // d11/d21 <===*
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit2;
                mInfo.mFit2  = body2_fit1;
            }
        }
        else if ( body1_d2.dot(body2_d2)       >  0.0 && 
                  fabs(body2_n2.dot(body1_d2)) <= mEpsilonAngle ) {
            if (isCCW(body1_d1, body1_d2, body2_d1, body1_v12)) {
                // d12/d22 <===*
                //         .../.\
                //         .d21.d11
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit2;
                mInfo.mFit2  = body2_fit1;
            }
            else {
                // d12/d22 <===*
                //         .../.\
                //         .d11.d21
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit1;
                mInfo.mFit2  = body2_fit2;
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
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit2;
                    mInfo.mFit2  = body2_fit1;
                }
                else {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit1;
                    mInfo.mFit2  = body2_fit2;
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
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit2;
                mInfo.mFit2  = body2_fit1;
            }
            else if (  ccw_11_12_22 && !ccw_11_12_21 && 
                      !ccw_11_22_21 &&  ccw_12_22_21    ) {
                // d11->d21->d12->d22  (d11 is penetrating into d22)
                //
                //            _-*-_
                //          _- / \ -_
                //      d12-  /   \  -d21
                //          d22   d11
                mInfo.mType1 = ContactPairInfo::FT_FACE;
                mInfo.mType2 = ContactPairInfo::FT_FACE;
                mInfo.mFit1  = body1_fit1;
                mInfo.mFit2  = body2_fit2;
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
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit1;
                    mInfo.mFit2  = body2_fit2;
                }
                else {
                    mInfo.mType1 = ContactPairInfo::FT_FACE;
                    mInfo.mType2 = ContactPairInfo::FT_FACE;
                    mInfo.mFit1  = body1_fit2;
                    mInfo.mFit2  = body2_fit1;
                }
            }
            else {
                // d11->d12->d21->d22  (Forbidden)
                log(ERROR, __FILE__, __LINE__, "Forbidden predicates");
            }
        }
    }
}


bool ContactPointsAndNormalGenerator::isCCW(
    const Vec3& v1,
    const Vec3& v2,
    const Vec3& v3,
    const Vec3& axis
) {
    auto cr13 = v1.cross(v3);
    auto cr21 = v2.cross(v1);
    auto cr32 = v3.cross(v2);
    if (cr13.squaredNorm2() <= mEpsilonZero) {
        if (v1.dot(v3) < 0.0 ) {
            return axis.dot(v1.cross(v2))>= 0.0 &&
                   axis.dot(v2.cross(v3))>= 0.0 ;
        }
        else {
            return true;
        }
    }
    else if (cr21.squaredNorm2() <= mEpsilonZero) {
        if (v2.dot(v1) <  0.0 ) { 
            return axis.dot(v2.cross(v3))>= 0.0 &&
                   axis.dot(v3.cross(v1))>= 0.0 ;
        }
        else {
            return true;
        }
    }
    else if (cr32.squaredNorm2() <= mEpsilonZero) {
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


/** @brief Strategy: 
 *         We only check shallow penetration of incident features, and
 *         do not perform deep penetration checks on the followin basis.
 *         - Edge-Vertex pair occurs rarely.
 *         - If there is a deep penetration, there would be more than
 *           two boundary vertices (2- or 3-simplex) in the binary dilation.
 *         - Deep penetration check is too complex.
 */
void ContactPointsAndNormalGenerator::adjustFeatures_EDGE_VERTEX()
{
    auto he1  = (*mInfo.mEit1)->he1();
    auto he2  = (*mInfo.mEit1)->he2();

    auto fit1 = (*he1)->face();
    auto fit2 = (*he2)->face();

    auto n1 = (*fit1)->nGCS(mQmat1);
    n1.normalize();
    auto n2 = (*fit2)->nGCS(mQmat1);
    n2.normalize();
    auto d  = (*((*he1)->dst()))->pGCS(mQmat1) - 
              (*((*he1)->src()))->pGCS(mQmat1) ;
    d.normalize();

    auto d1 = n1.cross(d);
    auto d2 = d.cross(n2);

    auto p =  (*mInfo.mVit2)->pGCS(mQmat2);

    FaceIt FACE_FACE_fit1;
    FaceIt FACE_FACE_fit2;
    double FACE_FACE_minVal;
    bool   FACE_FACE_found = false;

    FaceIt FACE_EDGE_fit;
    EdgeIt FACE_EDGE_eit;
    double FACE_EDGE_minVal;
    bool   FACE_EDGE_found = false;

    EdgeIt EDGE_EDGE_eit;
    double EDGE_EDGE_minVal;
    bool   EDGE_EDGE_found = false;

    FaceIt EDGE_FACE_fit;
    double EDGE_FACE_minVal;
    bool   EDGE_FACE_found = false;

    for (auto he : (*mInfo.mVit2)->halfEdges()) {
        if ((*he)->src()==mInfo.mVit2) {
            auto adjit = (*he)->dst();
            auto fitTest = (*he)->face();
            auto nTest   = (*fitTest)->nGCS(mQmat2);
            nTest.normalize();

            auto dTest  = (p - (*adjit)->pGCS(mQmat2));
            dTest.normalize();                          

            auto hePrev = (*he)->prev();
            auto previt = (*hePrev)->src(); 
            auto dPrev = (p - (*previt)->pGCS(mQmat2));
            dPrev.normalize();

            auto cr1 = n1.cross(nTest).squaredNorm2();

            if ( (cr1 <= mEpsilonAngle                      ) && 
                 (dTest.dot(d1) < 0.0 || dPrev.dot(d1) < 0.0)    ) {
                if (!FACE_FACE_found) {
                    FACE_FACE_fit1   = fit1;
                    FACE_FACE_fit2   = fitTest;
                    FACE_FACE_minVal = cr1;
                    FACE_FACE_found  = true;
                }
                else if (cr1 < FACE_FACE_minVal) {
                    FACE_FACE_fit1   = fit1;
                    FACE_FACE_fit2   = fitTest;
                    FACE_FACE_minVal = cr1;
                }
            }
            auto cr2 = n2.cross(nTest).squaredNorm2();

            if ( (cr2 <= mEpsilonAngle                      ) && 
                 (dTest.dot(d2) < 0.0 || dPrev.dot(d2) < 0.0)    ) {
                if (!FACE_FACE_found) {
                    FACE_FACE_fit1   = fit2;
                    FACE_FACE_fit2   = fitTest;
                    FACE_FACE_minVal = cr2;
                    FACE_FACE_found  = true;
                }
                else if (cr1 < FACE_FACE_minVal) {
                    FACE_FACE_fit1   = fit2;
                    FACE_FACE_fit2   = fitTest;
                    FACE_FACE_minVal = cr2;
                }
            }

            auto dot1 = fabs(n1.dot(dTest));
            if (dot1 <= mEpsilonAngle && dTest.dot(d1) < 0.0 ) {

                if (!FACE_EDGE_found) {
                    FACE_EDGE_fit    = fit1;
                    FACE_EDGE_eit    = (*he)->edge();
                    FACE_EDGE_minVal = dot1;
                    FACE_EDGE_found  = true;
                }
                else if (dot1 < FACE_EDGE_minVal) {
                    FACE_EDGE_fit    = fit1;
                    FACE_EDGE_eit    = (*he)->edge();
                    FACE_EDGE_minVal = dot1;
                }
            }

            auto dot2 = fabs(n2.dot(dTest));
            if (dot2 <= mEpsilonAngle && dTest.dot(d2) < 0.0 ) {
                if (!FACE_EDGE_found) {
                    FACE_EDGE_fit    = fit2;
                    FACE_EDGE_eit    = (*he)->edge();
                    FACE_EDGE_minVal = dot2;
                    FACE_EDGE_found  = true;
                }
                else if (dot2 < FACE_EDGE_minVal) {
                    FACE_EDGE_fit    = fit2;
                    FACE_EDGE_eit    = (*he)->edge();
                    FACE_EDGE_minVal = dot2;
                }
            }

            auto cr3 = d.cross(dTest).squaredNorm2();
            if (cr3 <= mEpsilonAngle) {
                if (!EDGE_EDGE_found) {
                    EDGE_EDGE_eit    = (*he)->edge();
                    EDGE_EDGE_minVal = cr3;
                    EDGE_EDGE_found  = true;
                }
                else if (cr3 < EDGE_EDGE_minVal) {
                    EDGE_EDGE_eit    = (*he)->edge();
                    EDGE_EDGE_minVal = cr3;
                }
            }
            auto dot3 = fabs(d.dot(nTest));
            if (dot3 <= mEpsilonAngle &&
                (isCCW(dTest*-1.0, d,      dPrev*-1.0, nTest)||
                 isCCW(dTest*-1.0, d*-1.0, dPrev*-1.0, nTest)  )) {
                if (!EDGE_FACE_found) {
                    EDGE_FACE_fit    = fitTest;
                    EDGE_FACE_minVal = dot3;
                    EDGE_FACE_found  = true;
                }
                else if (dot3 < EDGE_FACE_minVal) {
                    EDGE_FACE_fit    = fitTest;
                    EDGE_FACE_minVal = dot3;
                }
            }
        }
    }

    if (EDGE_EDGE_found) {

        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = EDGE_EDGE_eit;

    }
    else if (FACE_FACE_found) {

        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit1  = FACE_FACE_fit1;
        mInfo.mFit2  = FACE_FACE_fit2;

    }
    else if (FACE_EDGE_found) {

        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mFit1  = FACE_EDGE_fit;
        mInfo.mEit2  = FACE_EDGE_eit;

    }
    else if (EDGE_FACE_found) {

        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2  = EDGE_FACE_fit;

    }
}


/** @brief this is a reversed version of adjustFeatures_EDGE_VERTEX()
 *         see comments there.
 */
void ContactPointsAndNormalGenerator::adjustFeatures_VERTEX_EDGE()
{
    auto he1  = (*mInfo.mEit2)->he1();
    auto he2  = (*mInfo.mEit2)->he2();

    auto fit1 = (*he1)->face();
    auto fit2 = (*he2)->face();

    auto n1   = (*fit1)->nGCS(mQmat2);
    n1.normalize();
    auto n2   = (*fit2)->nGCS(mQmat2);
    n2.normalize();

    auto d    = (*((*he1)->dst()))->pGCS(mQmat2) - 
                (*((*he1)->src()))->pGCS(mQmat2) ;
    d.normalize();

    auto d1 = n1.cross(d);
    auto d2 = d.cross(n2);

    auto p = (*mInfo.mVit1)->pGCS(mQmat1);

    FaceIt FACE_FACE_fit1;
    FaceIt FACE_FACE_fit2;
    double FACE_FACE_minVal;
    bool   FACE_FACE_found = false;

    EdgeIt EDGE_FACE_eit;
    FaceIt EDGE_FACE_fit;
    double EDGE_FACE_minVal;
    bool   EDGE_FACE_found = false;

    EdgeIt EDGE_EDGE_eit;
    double EDGE_EDGE_minVal;
    bool   EDGE_EDGE_found = false;

    FaceIt FACE_EDGE_fit;
    double FACE_EDGE_minVal;
    bool   FACE_EDGE_found = false;

    for (auto he : (*mInfo.mVit1)->halfEdges()) {
        if ((*he)->src()==mInfo.mVit1) {

            VertexIt adjit = (*he)->dst();
            auto fitTest = (*he)->face();
            auto nTest   = (*fitTest)->nGCS(mQmat1);
            nTest.normalize();

            auto dTest  = (p - (*adjit)->pGCS(mQmat1));
            dTest.normalize();                          

            auto hePrev = (*he)->prev();
            auto previt = (*hePrev)->src(); 
            auto dPrev = (p - (*previt)->pGCS(mQmat1));
            dPrev.normalize();

            auto cr1 = n1.cross(nTest).squaredNorm2();

            if ( (cr1 <= mEpsilonAngle                      ) && 
                 (dTest.dot(d1) < 0.0 || dPrev.dot(d1) < 0.0)    ) {
                if (!FACE_FACE_found) {
                    FACE_FACE_fit1   = fitTest;
                    FACE_FACE_fit2   = fit1;
                    FACE_FACE_minVal = cr1;
                    FACE_FACE_found  = true;
                }
                else if (cr1 < FACE_FACE_minVal) {
                    FACE_FACE_fit1   = fitTest;
                    FACE_FACE_fit2   = fit1;
                    FACE_FACE_minVal = cr1;
                }
            }
            auto cr2 = n2.cross(nTest).squaredNorm2();

            if ( (cr2 <= mEpsilonAngle                      ) && 
                 (dTest.dot(d2) < 0.0 || dPrev.dot(d2) < 0.0)    ) {
                if (!FACE_FACE_found) {
                    FACE_FACE_fit1   = fitTest;
                    FACE_FACE_fit2   = fit2;
                    FACE_FACE_minVal = cr2;
                    FACE_FACE_found  = true;
                }
                else if (cr1 < FACE_FACE_minVal) {
                    FACE_FACE_fit1   = fitTest;
                    FACE_FACE_fit2   = fit2;
                    FACE_FACE_minVal = cr2;
                }
            }

            auto dot1 = fabs(n1.dot(dTest));
            if (dot1 <= mEpsilonAngle && dTest.dot(d1) < 0.0 ) {
                if (!EDGE_FACE_found) {
                    EDGE_FACE_fit    = fit1;
                    EDGE_FACE_eit    = (*he)->edge();
                    EDGE_FACE_minVal = dot1;
                    EDGE_FACE_found  = true;
                }
                else if (dot1 < EDGE_FACE_minVal) {
                    EDGE_FACE_fit    = fit1;
                    EDGE_FACE_eit    = (*he)->edge();
                    EDGE_FACE_minVal = dot1;
                }
            }

            auto dot2 = fabs(n2.dot(dTest));
            if (dot2 <= mEpsilonAngle && dTest.dot(d2) < 0.0 ) {
                if (!EDGE_FACE_found) {
                    EDGE_FACE_fit    = fit2;
                    EDGE_FACE_eit    = (*he)->edge();
                    EDGE_FACE_minVal = dot2;
                    EDGE_FACE_found  = true;
                }
                else if (dot2 < EDGE_FACE_minVal) {
                    EDGE_FACE_fit    = fit2;
                    EDGE_FACE_eit    = (*he)->edge();
                    EDGE_FACE_minVal = dot2;
                }
            }

            auto cr3 = d.cross(dTest).squaredNorm2();
            if (cr3 <= mEpsilonAngle) {
                if (!EDGE_EDGE_found) {
                    EDGE_EDGE_eit    = (*he)->edge();
                    EDGE_EDGE_minVal = cr3;
                    EDGE_EDGE_found  = true;
                }
                else if (cr3 < EDGE_EDGE_minVal) {
                    EDGE_EDGE_eit    = (*he)->edge();
                    EDGE_EDGE_minVal = cr3;
                }
            }

            auto dot3 = fabs(d.dot(nTest));
            if (dot3 <= mEpsilonAngle) {
                if (!FACE_EDGE_found) {
                    FACE_EDGE_fit    = fitTest;
                    FACE_EDGE_minVal = dot3;
                    FACE_EDGE_found  = true;
                }
                else if (dot3 < FACE_EDGE_minVal) {
                    FACE_EDGE_fit    = fitTest;
                    FACE_EDGE_minVal = dot3;
                }
            }
        }
    }

    if (EDGE_EDGE_found) {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = EDGE_EDGE_eit;

    }
    else if (FACE_FACE_found) {

        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit1  = FACE_FACE_fit1;
        mInfo.mFit2  = FACE_FACE_fit2;

    }
    else if (EDGE_FACE_found) {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mEit1  = EDGE_FACE_eit;
        mInfo.mFit2  = EDGE_FACE_fit;
    }
    else if (FACE_EDGE_found) {

        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1  = FACE_EDGE_fit;

    }
}


/** @brief Strategy: 
 *         We only check shallow penetration of incident features, and
 *         do not perform deep penetration checks on the followin basis.
 *         - Vertex-Vertex pair occurs rarely.
 *         - If there is a deep penetration, there would be more than
 *           two boundary vertices (2- or 3-simplex) in the binary dilation.
 *         - Deep penetration check is too complex.
 */
void ContactPointsAndNormalGenerator::adjustFeatures_VERTEX_VERTEX()
{
    auto p1 = (*mInfo.mVit1)->pGCS(mQmat1);
    auto p2 = (*mInfo.mVit2)->pGCS(mQmat2);

    FaceIt FACE_FACE_fit1;
    FaceIt FACE_FACE_fit2;
    double FACE_FACE_minVal;
    bool   FACE_FACE_found = false;

    EdgeIt EDGE_FACE_eit;
    FaceIt EDGE_FACE_fit;
    double EDGE_FACE_minVal;
    bool   EDGE_FACE_found = false;

    EdgeIt FACE_EDGE_eit;
    FaceIt FACE_EDGE_fit;
    double FACE_EDGE_minVal;
    bool   FACE_EDGE_found = false;

    EdgeIt EDGE_EDGE_eit1;
    EdgeIt EDGE_EDGE_eit2;
    double EDGE_EDGE_minVal;
    bool   EDGE_EDGE_found = false;

    for (auto he1 : (*mInfo.mVit1)->halfEdges()) {
        if ((*he1)->src()==mInfo.mVit1) {
            VertexIt adjit1 = (*he1)->dst();
            auto fit1 = (*he1)->face();
            auto n1 = (*fit1)->nGCS(mQmat1);
            n1.normalize();
            auto d1 = (p1 - (*adjit1)->pGCS(mQmat1));
            d1.normalize();
            auto hePrev1 = (*he1)->prev();
            auto previt1 = (*hePrev1)->src(); 
            auto dPrev1 = (p1 - (*previt1)->pGCS(mQmat1));
            dPrev1.normalize();

            for (auto he2 : (*mInfo.mVit2)->halfEdges()) {
                VertexIt adjit2 = (*he2)->dst();
                if ((*he2)->src()==mInfo.mVit2) {
                    auto fit2 = (*he2)->face();
                    auto n2  = (*fit2)->nGCS(mQmat2);
                    n2.normalize();
                    auto d2 = (p2 - (*adjit2)->pGCS(mQmat2));
                    d2.normalize();

                    auto hePrev2 = (*he2)->prev();
                    auto previt2 = (*hePrev2)->src(); 
                    auto dPrev2  = (p2 - (*previt2)->pGCS(mQmat2));
                    dPrev2.normalize();

                    auto cr1 = n1.cross(n2).squaredNorm2();

                    auto cr_11 = d1.cross(d2).squaredNorm2();
                    auto cr_12 = d1.cross(dPrev2).squaredNorm2();
                    auto cr_21 = dPrev1.cross(d2).squaredNorm2();
                    auto cr_22 = dPrev1.cross(dPrev2).squaredNorm2();

                    if (cr1 <= mEpsilonAngle) {
                        bool found = false;
                        if (cr_11 <= mEpsilonAngle ) {
                            if (cr_22 <= mEpsilonAngle) {
                                found = true;
                            }
                            else if(isCCW(d1*-1.0,dPrev2*-1.0,dPrev1*-1.0,n1)){
                                found = true;
                            }
                        }
                        else if(cr_12 <= mEpsilonAngle) {
                            if (cr_21 <= mEpsilonAngle) {
                                found = true;
                            }
                            else if(isCCW(d1*-1.0,d2*-1.0,dPrev1*-1.0,n1)) {
                                found = true;
                            }
                        }                        
                        else if (cr_21 <= mEpsilonAngle ) { 
                            if(isCCW(d1*-1.0,dPrev2*-1.0,dPrev1*-1.0,n1)){
                                found = true;
                            }
                        }
                        else if (cr_22 <= mEpsilonAngle ) { 
                            if(isCCW(d1*-1.0,d2*-1.0,dPrev1*-1.0,n1)){
                                found = true;
                            }
                        }
                        else if 
                            (isCCW(d1*-1.0, d2* -1.0,    dPrev1*-1.0, n1)||
                             isCCW(d1*-1.0, dPrev2*-1.0, dPrev1*-1.0, n1)  ) {
                            found = true;
                        }
                        if (found) {
                            if (!FACE_FACE_found) {
                                FACE_FACE_fit1   = fit1;
                                FACE_FACE_fit2   = fit2;
                                FACE_FACE_minVal = cr1;
                                FACE_FACE_found  = true;
                            }
                            else if (cr1 < FACE_FACE_minVal) {
                                FACE_FACE_fit1   = fit1;
                                FACE_FACE_fit2   = fit2;
                                FACE_FACE_minVal = cr1;
                            }
                        }
                    }

                    auto dot1 = fabs(n1.dot(d2));
                    if (dot1 <= mEpsilonAngle &&
                        cr_11 > mEpsilonAngle &&
                        cr_12 > mEpsilonAngle &&
                        isCCW(d1*-1.0, d2*-1.0, dPrev1*-1.0, n1) ) {
                        if (!FACE_EDGE_found) {
                            FACE_EDGE_fit    = fit1;
                            FACE_EDGE_eit    = (*he2)->edge();
                            FACE_EDGE_minVal = dot1;
                            FACE_EDGE_found  = true;
                        }
                        else if (dot1 < FACE_EDGE_minVal) {
                            FACE_EDGE_fit    = fit1;
                            FACE_EDGE_eit    = (*he2)->edge();
                            FACE_EDGE_minVal = dot1;
                        }
                    }

                    auto dot2 = fabs(n2.dot(d1));
                    if (dot2 <= mEpsilonAngle &&
                        cr_11 > mEpsilonAngle &&
                        cr_21 > mEpsilonAngle &&
                        isCCW(d2*-1.0, d1*-1.0, dPrev2*-1.0, n2) ) {
                        if (!EDGE_FACE_found) {
                            EDGE_FACE_eit    = (*he1)->edge();
                            EDGE_FACE_fit    = fit2;
                            EDGE_FACE_minVal = dot2;
                            EDGE_FACE_found  = true;
                        }
                        else if (dot2 < EDGE_FACE_minVal) {
                            EDGE_FACE_eit    = (*he1)->edge();
                            EDGE_FACE_fit    = fit2;
                            EDGE_FACE_minVal = dot2;
                        }
                    }

                    if (cr_11 <= mEpsilonAngle && d1.dot(d2)>0.0) {
                        if (!EDGE_EDGE_found) {
                            EDGE_EDGE_eit1   = (*he1)->edge();
                            EDGE_EDGE_eit2   = (*he2)->edge();
                            EDGE_EDGE_minVal = cr_11;
                            EDGE_EDGE_found  = true;
                        }
                        else if (cr_11 < EDGE_EDGE_minVal) {
                            EDGE_EDGE_eit1   = (*he1)->edge();
                            EDGE_EDGE_eit2   = (*he2)->edge();
                            EDGE_EDGE_minVal = cr_11;
                        }
                    }
                }
            }
        }
    }

    if (FACE_FACE_found) {

        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit1 = FACE_FACE_fit1;
        mInfo.mFit2 = FACE_FACE_fit2;

    }
    else if (EDGE_EDGE_found) {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1 = EDGE_EDGE_eit1;
        mInfo.mEit2 = EDGE_EDGE_eit2;

    }
    else if (FACE_EDGE_found) {

        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mFit1 = FACE_EDGE_fit;
        mInfo.mEit2 = FACE_EDGE_eit;

    }
    else if (EDGE_FACE_found) {
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mEit1 = EDGE_FACE_eit;
        mInfo.mFit2 = EDGE_FACE_fit;

    }
}


void ContactPointsAndNormalGenerator::findVerticesOfRigidBody(
    vector<VertexIt>& vits1,
    vector<VertexIt>& vits2
) {

    for (auto& q : mQ) {

        bool found1 = false;
        for (auto& v1 : vits1) {
            if (q.v1() ==v1) {
                found1 = true;
                break;
            }
        }
        if (!found1) {
            vits1.push_back(q.v1());
        }

        bool found2 = false;
        for (auto& v2 : vits2) {
            if (q.v2() ==v2) {
                found2 = true;
                break;
            }
        }
        if (!found2) {
            vits2.push_back(q.v2());
        }
    }

}



void ContactPointsAndNormalGenerator::findFeatureOfRigidBody(
    Manifold&                          ch,
    vector<VertexIt>&                  vits,
    FaceIt&                            fit,
    EdgeIt&                            eit,
    VertexIt&                          vit,
    enum ContactPairInfo::FeatureType& type
) {
    auto  noFace   = ch.faces().second;
    auto  noEdge   = ch.edges().second;
    auto  noVertex = ch.vertices().second;

    fit  = noFace;
    eit  = noEdge;
    vit  = noVertex;
    type = ContactPairInfo::FT_NONE;

    if (vits.size() >= 3) {

        fit = ch.findFace(vits[0], vits[1], vits[2]);
        if (fit == noFace) {
            log(WARNING, __FILE__, __LINE__,
                                       "Can't find a face from three points");
            fit = ch.findFace(vits[0],vits[1]);
            if (fit == noFace) {
                fit = ch.findFace(vits[1],vits[2]);
                if (fit == noFace) {
                    fit = ch.findFace(vits[2],vits[0]);
                }
            }
        }

        if (fit == noFace) {
            log(ERROR, __FILE__, __LINE__,
                               "Can't find a face from pairs of three points");
            eit = ch.findEdge(vits[0],vits[1]);
            if (eit == noEdge) {
                eit = ch.findEdge(vits[1],vits[2]);
                if (eit == noEdge) {
                    eit = ch.findEdge(vits[2],vits[0]);
                }
            }

            if (eit == noEdge) {
                log(ERROR, __FILE__, __LINE__,
                              "Can't find an edge from pairs of three points");
                vit = vits[0];
                type = ContactPairInfo::FT_VERTEX;
            }
            else {
                type = ContactPairInfo::FT_EDGE;
            }
        }
        else{
            type = ContactPairInfo::FT_FACE;
        }

    }
    else if (vits.size() == 2) {
        eit = ch.findEdge(vits[0],vits[1]);
        if (eit == noEdge) {
            fit = ch.findFace(vits[0],vits[1]);
            if (fit == noFace) {
                log(ERROR, __FILE__, __LINE__,
                                 "Can't find an edge of a face from a pair.");

                vit = vits[0];
                type = ContactPairInfo::FT_VERTEX;
            }
            else {
                type = ContactPairInfo::FT_FACE;
            }
        }
        else {
            type = ContactPairInfo::FT_EDGE;
        }
    }
    else {
        vit = vits[0];
        type = ContactPairInfo::FT_VERTEX;
    }

}


#ifdef UNIT_TESTS


void ContactPointsAndNormalGenerator::updateGeomConfigInternally()
{
    if (mUseGeomConfigTemp) {
        mQmat1 = mBody1.QmatTemp();
        mQmat2 = mBody2.QmatTemp();
        mCoM1  = mBody1.CoMTemp();
        mCoM2  = mBody2.CoMTemp();
    }
    else {
        mQmat1 = mBody1.Qmat();
        mQmat2 = mBody2.Qmat();
        mCoM1  = mBody1.CoM();
        mCoM2  = mBody2.CoM();
    }
}


void ContactPointsAndNormalGenerator::renderObj1ActiveFace(
    Vec3&         color,
    float         alpha,
    vector<Vec3>& vertices,
    vector<Vec3>& colors,
    vector<float>&alphas,
    vector<Vec3>& normals
){
    if (mInfo.mType1==ContactPairInfo::FT_FACE) {
        mBody1.ConvexHull().makeOpenGLVerticesColorsNormalsForTriangles(
            mInfo.mFit1,
            color,
            alpha,
            vertices,
            colors,
            alphas,
            normals, 
            true,
            mQmat1,
            mCoM1
        );

    }
}


void ContactPointsAndNormalGenerator::renderObj2ActiveFace(
    Vec3&         color,
    float         alpha,
    vector<Vec3>& vertices,
    vector<Vec3>& colors,
    vector<float>&alphas,
    vector<Vec3>& normals
){
    if (mInfo.mType2==ContactPairInfo::FT_FACE) {
        mBody2.ConvexHull().makeOpenGLVerticesColorsNormalsForTriangles(
            mInfo.mFit2,
            color,
            alpha,
            vertices,
            colors,
            alphas,
            normals, 
            true,
            mQmat2,
            mCoM2
        );
    }
}


void ContactPointsAndNormalGenerator::renderObj1ActiveEdge(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
){
    if (mInfo.mType1==ContactPairInfo::FT_EDGE) {

        auto vit1 = (*((*mInfo.mEit1)->he1()))->src();
        auto vit2 = (*((*mInfo.mEit1)->he1()))->dst();
        auto p1   = (*vit1)->pGCS(mQmat1, mCoM1);
        auto p2   = (*vit2)->pGCS(mQmat1, mCoM1);

        vector<Vec3> points;
        for (long i =0; i <= 100; i++) {
            Vec3 p = p1 * (i/100.0) + p2* (1.0 - i/100.0);
            points.push_back(p);
        }

        Manifold::makeOpenGLVerticesColorsForPoints(
            points,
            color,
            vertices,
            colors
        );

    }
}


void ContactPointsAndNormalGenerator::renderObj2ActiveEdge(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
){
    if (mInfo.mType2==ContactPairInfo::FT_EDGE) {

        auto vit1 = (*((*mInfo.mEit2)->he1()))->src();
        auto vit2 = (*((*mInfo.mEit2)->he1()))->dst();
        auto p1 = (*vit1)->pGCS(mQmat2, mCoM2);
        auto p2 = (*vit2)->pGCS(mQmat2, mCoM2);
        vector<Vec3> points;
        for (long i =0; i <= 100; i++) {
            Vec3 p = p1 * (i/100.0) + p2* (1.0 - i/100.0);
            points.push_back(p);
        }
        Manifold::makeOpenGLVerticesColorsForPoints(
            points,
            color,
            vertices,
            colors
        );
    }
}


void ContactPointsAndNormalGenerator::renderObj1ActiveVertex(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
){
    if (mInfo.mType1==ContactPairInfo::FT_VERTEX) {

        auto p = (*mInfo.mVit1)->pGCS(mQmat1, mCoM1);

        vector<Vec3> points;
        points.push_back(p);

        Manifold::makeOpenGLVerticesColorsForPoints(
            points,
            color,
            vertices,
            colors
        );
    }
}


void ContactPointsAndNormalGenerator::renderObj2ActiveVertex(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
){
    if (mInfo.mType2==ContactPairInfo::FT_VERTEX) {

        auto p = (*mInfo.mVit2)->pGCS(mQmat2, mCoM2);

        vector<Vec3> points;
        points.push_back(p);

        Manifold::makeOpenGLVerticesColorsForPoints(
            points,
            color,
            vertices,
            colors
        );

    }
}


void ContactPointsAndNormalGenerator::renderContactPoints(
    Vec3&         color1,
    Vec3&         color2,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
) {

    vector<Vec3> points1;
    for (auto& pLCS : mInfo.mIntsect1) {
        auto pGCS = (mQmat1 * pLCS) + mCoM1;
        points1.push_back(pGCS);
    }

    Manifold::makeOpenGLVerticesColorsForPoints(
        points1,
        color1,
        vertices,
        colors
    );

    vector<Vec3> points2;
    for (auto& pLCS : mInfo.mIntsect2) {
        auto pGCS = (mQmat2 * pLCS) + mCoM2;
        points2.push_back(pGCS);
    }

    Manifold::makeOpenGLVerticesColorsForPoints(
        points2,
        color2,
        vertices,
        colors
    );
}


string ContactPointsAndNormalGenerator::FeatureTypeString(ContactPairInfo::FeatureType t)
{
    switch (t) {
      case ContactPairInfo::FT_NONE:
        return "FT_NONE";;
        break;
      case ContactPairInfo::FT_VERTEX:
        return "FT_VERTEX";
        break;
      case ContactPairInfo::FT_EDGE:
        return "FT_EDGE";
        break;
      case ContactPairInfo::FT_FACE:
        return "FT_FACE";
        break;
    };
    return "";
}

#endif

}// namespace Makena


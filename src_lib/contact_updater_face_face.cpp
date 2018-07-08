#include "contact_updater_face_face.hpp"

/**
 * @file contact_updater_face_face.cpp
 *
 * @brief
 */
namespace Makena {


long ContactUpdater_FACE_FACE::edgeIndexA(IntSec2D::OutputElem& e)
{
    if (e.mIndexA==0&&e.mIndexAaux==mLastIndexA) {
        return e.mIndexAaux;
    }
    else {
        return e.mIndexA;  
    }
}


long ContactUpdater_FACE_FACE::edgeIndexB(IntSec2D::OutputElem& e)
{
    if (e.mIndexB==0&&e.mIndexBaux==mLastIndexB) {
        return e.mIndexBaux;
    }
    else {
        return e.mIndexB;  
    }
}


bool ContactUpdater_FACE_FACE::isVertexIncidentToEdge(VertexIt v, HalfEdgeIt h)
{
    return (*h)->src()==v ||(*h)->dst()==v;
}


bool ContactUpdater_FACE_FACE::findHalfEdgeBody1(
    const long  index1,
    const long  index2,
    HalfEdgeIt& heOut
) {
    auto numEdges = mHalfEdges1.size();
    auto heit1 = mHalfEdges1[(index1+numEdges-1)%numEdges];
    auto heit2 = mHalfEdges1[index1];
    auto vit1  = mVertices1[index1];
    auto vit2  = mVertices1[index2];

    if ((*heit1)->dst()!= vit1) {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
    if ((*heit2)->src()!= vit1) {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
    if ((*heit1)->src()== vit2) {
        heOut = heit1;
        return true;
    }
    if ((*heit2)->dst()== vit2) {
        heOut = heit2;
        return true;
    }
    else {
        return false;
    }
}


bool ContactUpdater_FACE_FACE::findHalfEdgeBody2(
    const long  index1,
    const long  index2,
    HalfEdgeIt& heOut
) {
    auto numEdges = mHalfEdges2.size();
    auto heit1 = mHalfEdges2[(index1+numEdges-1)%numEdges];
    auto heit2 = mHalfEdges2[index1];
    auto vit1  = mVertices2[index1];
    auto vit2  = mVertices2[index2];

    if ((*heit1)->dst()!= vit1) {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
    if ((*heit2)->src()!= vit1) {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
    if ((*heit1)->src()== vit2) {
        heOut = heit1;
        return true;
    }
    if ((*heit2)->dst()== vit2) {
        heOut = heit2;
        return true;
    }
    else {
        return false;
    }
}


bool ContactUpdater_FACE_FACE::areParallel(
    EdgeIt     eit1,
    const long body1,
    EdgeIt     eit2,
    const long body2
) {
    auto he1   = (*eit1)->he1();
    auto he2   = (*eit2)->he1();
    auto vit11 = (*he1)->src();
    auto vit12 = (*he1)->dst();
    auto vit21 = (*he2)->src();
    auto vit22 = (*he2)->dst();
    Vec3 p11, p12, p21, p22;

    if (body1==1) {
        p11 = (*vit11)->pGCS(mRotMat1, mCom1);
        p12 = (*vit12)->pGCS(mRotMat1, mCom1);
    }
    else {
        p11 = (*vit11)->pGCS(mRotMat2, mCom2);
        p12 = (*vit12)->pGCS(mRotMat2, mCom2);
    }

    if (body2==1) {
        p21 = (*vit21)->pGCS(mRotMat1, mCom1);
        p22 = (*vit22)->pGCS(mRotMat1, mCom1);
    }
    else {
        p21 = (*vit21)->pGCS(mRotMat2, mCom2);
        p22 = (*vit22)->pGCS(mRotMat2, mCom2);
    }

    auto v1 = p12 - p11;
    auto v2 = p22 - p21;
    if (v1.squaredNorm2() <= mEpsilonZero ||
        v2.squaredNorm2() <= mEpsilonZero   ) {
        return true;
    }
    v1.normalize();
    v2.normalize();
    auto cr = v1.cross(v2);

    return cr.squaredNorm2() <= mEpsilonAngle;

}


ContactUpdater_FACE_FACE::ContactUpdater_FACE_FACE(
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


ContactUpdater_FACE_FACE::~ContactUpdater_FACE_FACE(){;}


bool ContactUpdater_FACE_FACE:: update()
{

    bool isTilted, abort;

    checkForTilt(isTilted, abort);

    if (abort) {
        return true;
    }

    rotateWorld(mInfo.mContactNormal1To2);

    bool intersectionFound = findIntersection2D();

    if (intersectionFound) {

        if (isTilted) {

            long index1, index2;
            Vec2 tiltDir2D = findTilt();

            findExtremeFeature2D(tiltDir2D, index1, index2);

            dispatchAndUpdate(index1, index2);

        }
        else {
            checkAndUpdateEdgeCases();

        }

        return  false;        
    }
    else {
        return true;
    }
}


void ContactUpdater_FACE_FACE::checkAndUpdateEdgeCases()
{
    if (mIntsec.size()==1) {
        auto& oe = mIntsec[0];
        switch (oe.mType) {        
          case IntSec2D::IT_EDGE_VERTEX:
            updateNonTilting_EDGE_VERTEX();
            break;
          case IntSec2D::IT_VERTEX_EDGE:
            updateNonTilting_VERTEX_EDGE();
            break;
          case IntSec2D::IT_VERTEX_VERTEX:
            updateNonTilting_VERTEX_VERTEX();
            break;
          default:
            break;
        }
    }
    else if (mIntsec.size()==2) {

        auto& oe1 = mIntsec[0];
        auto& oe2 = mIntsec[1];
        switch (oe1.mType) {
          case IntSec2D::IT_EDGE_VERTEX:
            switch (oe2.mType) {
              case IntSec2D::IT_EDGE_VERTEX:
                updateNonTilting_EDGE_VERTEX_EDGE_VERTEX();
                break;
              case IntSec2D::IT_VERTEX_EDGE:
                updateNonTilting_EDGE_VERTEX_VERTEX_EDGE();
                break;
              case IntSec2D::IT_VERTEX_VERTEX:
                updateNonTilting_EDGE_VERTEX_VERTEX_VERTEX();
                break;
              default:
                break;
            }
            break;
          case IntSec2D::IT_VERTEX_EDGE:
            switch (oe2.mType) {
              case IntSec2D::IT_EDGE_VERTEX:
                updateNonTilting_VERTEX_EDGE_EDGE_VERTEX();
                break;
              case IntSec2D::IT_VERTEX_EDGE:
                updateNonTilting_VERTEX_EDGE_VERTEX_EDGE();
                break;
              case IntSec2D::IT_VERTEX_VERTEX:
                updateNonTilting_VERTEX_EDGE_VERTEX_VERTEX();
                break;
              default:
                break;
            }

            break;
          case IntSec2D::IT_VERTEX_VERTEX:
            switch (oe2.mType) {
              case IntSec2D::IT_EDGE_VERTEX:
                updateNonTilting_VERTEX_VERTEX_EDGE_VERTEX();
                break;
              case IntSec2D::IT_VERTEX_EDGE:
                updateNonTilting_VERTEX_VERTEX_VERTEX_EDGE();
                break;
              case IntSec2D::IT_VERTEX_VERTEX:
                updateNonTilting_VERTEX_VERTEX_VERTEX_VERTEX();
                break;
              default:
                break;
            }
            break;
          default:
            break;
        }
    }
}


void ContactUpdater_FACE_FACE::updateNonTilting_VERTEX_VERTEX()
{
    auto& oe = mIntsec[0];
    auto vit1    = mVertices1 [oe.mIndexA];
    auto vit2    = mVertices2 [oe.mIndexB];

    mInfo.mType1 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit1  = vit1;

    mInfo.mType2 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit2  = vit2;
}


void ContactUpdater_FACE_FACE::updateNonTilting_VERTEX_EDGE()
{
    auto& oe = mIntsec[0];
    auto vit1    = mVertices1 [oe.mIndexA];
    auto heit2   = mHalfEdges2[edgeIndexB(oe)];

    mInfo.mType1 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit1  = vit1;
    mInfo.mType2 = ContactPairInfo::FT_EDGE;
    mInfo.mEit2  = (*heit2)->edge();
}


void ContactUpdater_FACE_FACE::updateNonTilting_EDGE_VERTEX()
{
    auto& oe = mIntsec[0];
    auto heit1   = mHalfEdges1[edgeIndexA(oe)];
    auto vit2    = mVertices2 [oe.mIndexB];

    mInfo.mType1 = ContactPairInfo::FT_EDGE;
    mInfo.mEit1  = (*heit1)->edge();
    mInfo.mType2 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit2  = vit2;
}


void ContactUpdater_FACE_FACE::updateNonTilting_EDGE_VERTEX_EDGE_VERTEX()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto heit11  = mHalfEdges1[edgeIndexA(oe1)];
    auto heit21  = mHalfEdges1[edgeIndexA(oe2)];

    HalfEdgeIt heit2;

    if ((*heit11)->edge()!=(*heit21)->edge()) {

        log(ERROR, __FILE__, __LINE__, "Logic Error");

    }
    else if (!findHalfEdgeBody2(oe1.mIndexB, oe2.mIndexB, heit2)) {

        log(ERROR, __FILE__, __LINE__, "Logic Error");

    }
    else {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*heit11)->edge();

        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = (*heit2)->edge();

    }
}


void ContactUpdater_FACE_FACE::updateNonTilting_EDGE_VERTEX_VERTEX_EDGE()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto heit11  = mHalfEdges1[edgeIndexA(oe1)];
    auto vit12   = mVertices2 [oe1.mIndexB];
    auto vit21   = mVertices1 [oe2.mIndexA];
    auto heit22  = mHalfEdges2[edgeIndexB(oe2)];

    if (!isVertexIncidentToEdge(vit21, heit11)) {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
    if (!isVertexIncidentToEdge(vit12, heit22)) {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }

    mInfo.mType1 = ContactPairInfo::FT_EDGE;
    mInfo.mEit1  = (*heit11)->edge();

    mInfo.mType2 = ContactPairInfo::FT_EDGE;
    mInfo.mEit2  = (*heit22)->edge();
}


void ContactUpdater_FACE_FACE::updateNonTilting_EDGE_VERTEX_VERTEX_VERTEX()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto heit11  = mHalfEdges1[edgeIndexA(oe1)];
    auto vit21   = mVertices1 [oe2.mIndexA];

    HalfEdgeIt heit2;

    if (!isVertexIncidentToEdge(vit21, heit11)) {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
    else if (!findHalfEdgeBody2(oe1.mIndexB, oe2.mIndexB, heit2)) {

        log(ERROR, __FILE__, __LINE__, "Logic Error");

    }

    mInfo.mType1 = ContactPairInfo::FT_EDGE;
    mInfo.mEit1  = (*heit11)->edge();

    mInfo.mType2 = ContactPairInfo::FT_EDGE;
    mInfo.mEit2  = (*heit2)->edge();
}


void ContactUpdater_FACE_FACE::updateNonTilting_VERTEX_EDGE_EDGE_VERTEX()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto vit11   = mVertices1 [oe1.mIndexA];
    auto heit12  = mHalfEdges2[edgeIndexB(oe1)];
    auto heit21  = mHalfEdges1[edgeIndexA(oe2)];
    auto vit22   = mVertices2 [oe2.mIndexB];

    if (!isVertexIncidentToEdge(vit11, heit21)) {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
    if (!isVertexIncidentToEdge(vit22, heit12)) {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }

    mInfo.mType1 = ContactPairInfo::FT_EDGE;
    mInfo.mEit1  = (*heit21)->edge();

    mInfo.mType2 = ContactPairInfo::FT_EDGE;
    mInfo.mEit2  = (*heit12)->edge();
}


void ContactUpdater_FACE_FACE::updateNonTilting_VERTEX_EDGE_VERTEX_EDGE()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto heit12  = mHalfEdges2[edgeIndexB(oe1)];
    auto heit22  = mHalfEdges2[edgeIndexB(oe2)];

    HalfEdgeIt heit1;

    if ((*heit12)->edge()!=(*heit22)->edge()) {

        log(ERROR, __FILE__, __LINE__, "Logic Error");

    }
    else if (!findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        log(ERROR, __FILE__, __LINE__, "Logic Error");

    }
    else {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*heit1)->edge();

        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = (*heit12)->edge();

    }
}


void ContactUpdater_FACE_FACE::updateNonTilting_VERTEX_EDGE_VERTEX_VERTEX()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto heit12  = mHalfEdges2[edgeIndexB(oe1)];
    auto vit22   = mVertices2 [oe2.mIndexB];

    HalfEdgeIt heit1;

    if (!isVertexIncidentToEdge(vit22, heit12)) {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
    else if (!findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        log(ERROR, __FILE__, __LINE__, "Logic Error");

    }

    mInfo.mType1 = ContactPairInfo::FT_EDGE;
    mInfo.mEit1  = (*heit1)->edge();

    mInfo.mType2 = ContactPairInfo::FT_EDGE;
    mInfo.mEit2  = (*heit12)->edge();
}


void ContactUpdater_FACE_FACE::updateNonTilting_VERTEX_VERTEX_EDGE_VERTEX()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto vit11   = mVertices1 [oe1.mIndexA];
    auto heit21  = mHalfEdges1[edgeIndexA(oe2)];

    HalfEdgeIt heit2;

    if (!isVertexIncidentToEdge(vit11, heit21)) {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
    else if (!findHalfEdgeBody2(oe1.mIndexB, oe2.mIndexB, heit2)) {

        log(ERROR, __FILE__, __LINE__, "Logic Error");

    }

    mInfo.mType1 = ContactPairInfo::FT_EDGE;
    mInfo.mEit1  = (*heit21)->edge();

    mInfo.mType2 = ContactPairInfo::FT_EDGE;
    mInfo.mEit2  = (*heit2)->edge();
}


void ContactUpdater_FACE_FACE::updateNonTilting_VERTEX_VERTEX_VERTEX_EDGE()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto vit12   = mVertices2 [oe1.mIndexB];
    auto heit22  = mHalfEdges2[edgeIndexB(oe2)];
    HalfEdgeIt heit1;

    if (!isVertexIncidentToEdge(vit12, heit22)) {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
    else if (!findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        log(ERROR, __FILE__, __LINE__, "Logic Error");

    }

    mInfo.mType1 = ContactPairInfo::FT_EDGE;
    mInfo.mEit1  = (*heit1)->edge();

    mInfo.mType2 = ContactPairInfo::FT_EDGE;
    mInfo.mEit2  = (*heit22)->edge();
}


void ContactUpdater_FACE_FACE::updateNonTilting_VERTEX_VERTEX_VERTEX_VERTEX()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    HalfEdgeIt heit1, heit2;
    if (!findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        log(ERROR, __FILE__, __LINE__, "Logic Error");

    }
    else if (!findHalfEdgeBody2(oe1.mIndexB, oe2.mIndexB, heit2)) {

        log(ERROR, __FILE__, __LINE__, "Logic Error");

    }

    mInfo.mType1 = ContactPairInfo::FT_EDGE;
    mInfo.mEit1  = (*heit1)->edge();

    mInfo.mType2 = ContactPairInfo::FT_EDGE;
    mInfo.mEit2  = (*heit2)->edge();

}


Vec2 ContactUpdater_FACE_FACE::findTilt()
{
    auto fit1 = mInfo.mFit1;
    auto fit2 = mInfo.mFit2;
    auto n1   = (*fit1)->nGCS(mRotMat1);
    auto n2   = (*fit2)->nGCS(mRotMat2);
    Vec2 d(n1.x() + n2.x(), n1.y() + n2.y());
    d.normalize();
    d.scale(-1.0);
    return d;

}


void ContactUpdater_FACE_FACE::rotateWorld(const Vec3& zDir)
{
    Mat3x3 rMat = rotMatAlignZDirToZAxis(zDir);

    mRotMat1 = rMat * mBody1.Qmat();
    mCom1    = rMat * mBody1.CoM();
    mRotMat2 = rMat * mBody2.Qmat();
    mCom2    = rMat * mBody2.CoM();

}


Mat3x3 ContactUpdater_FACE_FACE::rotMatAlignZDirToZAxis(const Vec3& zDir)
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


void ContactUpdater_FACE_FACE::checkForTilt(bool& tilted, bool& abort)
{

    auto fit1 = mInfo.mFit1;

    auto fit2 = mInfo.mFit2;

    auto n1   = (*fit1)->nGCS(mBody1.Qmat());

    auto n2   = (*fit2)->nGCS(mBody2.Qmat());

    if (n1.dot(n2) > -1.0 * mEpsilonAngle) {

        // Too much tilt between two normals.
        abort = true;
    }
    else {

        abort = false;
    }

    auto cr   = n1.cross(n2);

    if (cr.squaredNorm2() < mEpsilonAngle) {

        tilted = false;
    }
    else {

        tilted = true;
    }
}


bool ContactUpdater_FACE_FACE::findIntersection2D()
{
    // Prepare InputElems
    auto fit1 = mInfo.mFit1;
    auto fit2 = mInfo.mFit2;

    vector< IntSec2D::InputElem> inputElem1;
    vector< IntSec2D::InputElem> inputElem2;

    mVertices1.clear();
    mVertices2.clear();
    mHalfEdges1.clear();
    mHalfEdges2.clear();
    mIntsec.clear();

    /* 
     *         FACE 1 (CCW)         FACE 2(CW)
     *
     *          i+1     i             i     i+1
     *           *<-----* src     src *----->*
     *          /        \           /        \
     *         /    N/    \         /    ^     \          
     *         \    /     /         \   /N     /
     *          \  v     /           \________/
     *           \______/
     *        Y
     *        ^
     *        |  /
     *        | /
     *        |/
     *    ----+----------->X
     *       /.
     *      / .
     *     v  .
     *     Z = Common Normal
     */

    mLastIndexA = ((*fit1)->halfEdges().size()) - 1 ;
    long index = 0;
    for (auto he : (*fit1)->halfEdges()) {
        auto vit = (*he)->src();
        mHalfEdges1.push_back(he);
        mVertices1.push_back(vit);
        auto p   = (*vit)->pGCS(mRotMat1, mCom1);
	Vec2 p2d(p.x(), p.y());
	IntSec2D::InputElem ie(p2d, index);
        inputElem1.push_back(ie);

        index++;
    }


    // We have to reverse the order of the surrounding half edges
    // since the face of Body 2 is in CW in 2D plane.
    for (long i = 0; i < (*fit2)->halfEdges().size();i++) {
	IntSec2D::InputElem ie(Vec2(0.0, 0.0), 0);
        inputElem2.push_back(ie);
    }
    auto rIndex = ((*fit2)->halfEdges().size()) - 1 ;
    mLastIndexB = rIndex;
    index = 0;
    for (auto he : (*fit2)->halfEdges()) {
        auto vit = (*he)->src();
        mHalfEdges2.push_back(he);
        mVertices2.push_back(vit);
        auto p = (*vit)->pGCS(mRotMat2, mCom2);
        Vec2 p2d(p.x(), p.y());
        auto& ie  = inputElem2[rIndex];
        ie.mP     = p2d;
        ie.mIndex = index;
        index++;
        rIndex--;
    }



    IntSec2D IntsecFinder2D(
                      mEpsilonZero, mEpsilonZero, mEpsilonAngle, mLogStream);
    IntsecFinder2D.findIntersection(inputElem1, inputElem2, mIntsec);
   
    return mIntsec.size() > 0;

}


void ContactUpdater_FACE_FACE::findExtremeFeature2D(
    const Vec2& dir,
    long&       index1,
    long&       index2
) {

    double maxDist  = dir.dot(mIntsec[0].mP);
    long   maxIndex = 0;
    const  long numOutElems = mIntsec.size();

    // Find the points in the leaning direction.
    for (long i = 1; i < numOutElems; i++) {

        double curDist = dir.dot(mIntsec[i].mP);
        if (curDist > maxDist) {
            maxDist  = curDist;
            maxIndex = i;
        }
    }

    if (mIntsec.size() ==1) {
        index1 = maxIndex;
        index2 = -1;
        return;
    }

    // Check the neighbors.
   
    long prevIndex = (maxIndex + numOutElems- 1) % numOutElems;
    long nextIndex = (maxIndex + 1) % numOutElems;

    Vec2 vPrev     = mIntsec[maxIndex].mP  - mIntsec[prevIndex].mP;
    Vec2 vNext     = mIntsec[nextIndex].mP - mIntsec[maxIndex].mP;

    if (vPrev.squaredNorm2() <= mEpsilonZero) {
        index1 = prevIndex;
        index2 = maxIndex;
        return;
    }

    if (vNext.squaredNorm2() <= mEpsilonZero) {
        index1 = maxIndex;
        index2 = nextIndex;
        return;
    }

    vPrev.normalize();
    vNext.normalize();

    if (fabs(dir.dot(vPrev)) <= mEpsilonAngle) {

        index1 = prevIndex;
        index2 = maxIndex;
        return;
    }
    if (fabs(dir.dot(vNext)) <= mEpsilonAngle) {

        index1 = maxIndex;
        index2 = nextIndex;
        return;
    }
    index1 = maxIndex;
    index2 = -1;
}


void ContactUpdater_FACE_FACE::dispatchAndUpdate(
    const long index1,
    const long index2
) {
    if (index2==-1) {

        auto& oe = mIntsec[index1];
        switch(oe.mType) {


          case IntSec2D::IT_EDGE_EDGE:
            processIntsec_EDGE_EDGE(index1);
            break;

          case IntSec2D::IT_EDGE_VERTEX:
            processIntsec_EDGE_VERTEX(index1);
            break;

          case IntSec2D::IT_VERTEX_EDGE:
            processIntsec_VERTEX_EDGE(index1);
            break;

          case IntSec2D::IT_VERTEX_INTERIOR:
            processIntsec_VERTEX_INTERIOR(index1);
            break;

          case IntSec2D::IT_INTERIOR_VERTEX:
            processIntsec_INTERIOR_VERTEX(index1);
            break;

          case IntSec2D::IT_VERTEX_VERTEX:
            processIntsec_VERTEX_VERTEX(index1);
            break;

          default:
            break;
        }
    }
    else {

        auto& oe1 = mIntsec[index1];
        auto& oe2 = mIntsec[index2];
        switch(oe1.mType) {

          case IntSec2D::IT_EDGE_EDGE:
            switch(oe2.mType) {
              case IntSec2D::IT_EDGE_EDGE:
                processIntsec_EDGE_EDGE__EDGE_EDGE(index1, index2);
                break;
              case IntSec2D::IT_EDGE_VERTEX:
                processIntsec_EDGE_EDGE__EDGE_VERTEX(index1, index2);
                break;
              case IntSec2D::IT_VERTEX_EDGE:
                processIntsec_EDGE_EDGE__VERTEX_EDGE(index1, index2);
                break;
              case IntSec2D::IT_VERTEX_INTERIOR:
                processIntsec_EDGE_EDGE__VERTEX_INTERIOR(index1, index2);
                break;
              case IntSec2D::IT_INTERIOR_VERTEX:
                processIntsec_EDGE_EDGE__INTERIOR_VERTEX(index1, index2);
                break;
              case IntSec2D::IT_VERTEX_VERTEX:
                processIntsec_EDGE_EDGE__VERTEX_VERTEX(index1, index2);
                break;
              default:
                break;
            }
            break;

          case IntSec2D::IT_EDGE_VERTEX:
            switch(oe2.mType) {
              case IntSec2D::IT_EDGE_EDGE:
                processIntsec_EDGE_EDGE__EDGE_VERTEX(index2, index1);
                break;
              case IntSec2D::IT_EDGE_VERTEX:
                processIntsec_EDGE_VERTEX__EDGE_VERTEX(index1, index2);
                break;
              case IntSec2D::IT_VERTEX_EDGE:
                processIntsec_EDGE_VERTEX__VERTEX_EDGE(index1, index2);
                break;
              case IntSec2D::IT_VERTEX_INTERIOR:
                processIntsec_EDGE_VERTEX__VERTEX_INTERIOR(index1, index2);
                break;
              case IntSec2D::IT_INTERIOR_VERTEX:
                processIntsec_EDGE_VERTEX__INTERIOR_VERTEX(index1, index2);
                break;
              case IntSec2D::IT_VERTEX_VERTEX:
                processIntsec_EDGE_VERTEX__VERTEX_VERTEX(index1, index2);
                break;
              default:
                break;
            }
            break;

          case IntSec2D::IT_VERTEX_EDGE:
            switch(oe2.mType) {
              case IntSec2D::IT_EDGE_EDGE:
                processIntsec_EDGE_EDGE__VERTEX_EDGE(index2, index1);
                break;
              case IntSec2D::IT_EDGE_VERTEX:
                processIntsec_EDGE_VERTEX__VERTEX_EDGE(index2, index1);
                break;
              case IntSec2D::IT_VERTEX_EDGE:
                processIntsec_VERTEX_EDGE__VERTEX_EDGE(index1, index2);
                break;
              case IntSec2D::IT_VERTEX_INTERIOR:
                processIntsec_VERTEX_EDGE__VERTEX_INTERIOR(index1, index2);
                break;
              case IntSec2D::IT_INTERIOR_VERTEX:
                processIntsec_VERTEX_EDGE__INTERIOR_VERTEX(index1, index2);
                break;
              case IntSec2D::IT_VERTEX_VERTEX:
                processIntsec_VERTEX_EDGE__VERTEX_VERTEX(index1, index2);
                break;
              default:
                break;
            }
            break;

          case IntSec2D::IT_VERTEX_INTERIOR:
            switch(oe2.mType) {
              case IntSec2D::IT_EDGE_EDGE:
                processIntsec_EDGE_EDGE__VERTEX_INTERIOR(index2, index1);
                break;
              case IntSec2D::IT_EDGE_VERTEX:
                processIntsec_EDGE_VERTEX__VERTEX_INTERIOR(index2, index1);
                break;
              case IntSec2D::IT_VERTEX_EDGE:
                processIntsec_VERTEX_EDGE__VERTEX_INTERIOR(index2, index1);
                break;
              case IntSec2D::IT_VERTEX_INTERIOR:
                processIntsec_VERTEX_INTERIOR__VERTEX_INTERIOR(index1, index2);
                break;
              case IntSec2D::IT_INTERIOR_VERTEX:
                processIntsec_VERTEX_INTERIOR__INTERIOR_VERTEX(index1, index2);
                break;
              case IntSec2D::IT_VERTEX_VERTEX:
                processIntsec_VERTEX_VERTEX__VERTEX_INTERIOR(index2, index1);
                break;
              default:
                break;
            }
            break;

          case IntSec2D::IT_INTERIOR_VERTEX:
            switch(oe2.mType) {
              case IntSec2D::IT_EDGE_EDGE:
                processIntsec_EDGE_EDGE__INTERIOR_VERTEX(index2, index1);
                break;
              case IntSec2D::IT_EDGE_VERTEX:
                processIntsec_EDGE_VERTEX__INTERIOR_VERTEX(index2, index1);
                break;
              case IntSec2D::IT_VERTEX_EDGE:
                processIntsec_VERTEX_EDGE__INTERIOR_VERTEX(index2, index1);
                break;
              case IntSec2D::IT_VERTEX_INTERIOR:
                processIntsec_VERTEX_INTERIOR__INTERIOR_VERTEX(index2, index1);
                break;
              case IntSec2D::IT_INTERIOR_VERTEX:
                processIntsec_INTERIOR_VERTEX__INTERIOR_VERTEX(index1, index2);
                break;
              case IntSec2D::IT_VERTEX_VERTEX:
                processIntsec_VERTEX_VERTEX__INTERIOR_VERTEX(index2, index1);
                break;
              default:
                break;
            }
            break;

          case IntSec2D::IT_VERTEX_VERTEX:
            switch(oe2.mType) {
              case IntSec2D::IT_EDGE_EDGE:
                processIntsec_EDGE_EDGE__VERTEX_VERTEX(index2, index1);
                break;
              case IntSec2D::IT_EDGE_VERTEX:
                processIntsec_EDGE_VERTEX__VERTEX_VERTEX(index2, index1);
                break;
              case IntSec2D::IT_VERTEX_EDGE:
                processIntsec_VERTEX_EDGE__VERTEX_VERTEX(index2, index1);
                break;
              case IntSec2D::IT_VERTEX_INTERIOR:
                processIntsec_VERTEX_VERTEX__VERTEX_INTERIOR(index1, index2);
                break;
              case IntSec2D::IT_INTERIOR_VERTEX:
                processIntsec_VERTEX_VERTEX__INTERIOR_VERTEX(index1, index2);
                break;
              case IntSec2D::IT_VERTEX_VERTEX:
                processIntsec_VERTEX_VERTEX__VERTEX_VERTEX(index1, index2);
                break;
              default:
                break;
            }
            break;

          default:
            break;
        }
    }
}


/*
 *          - If Pint = IT_EDGE_EDGE (Eint1,Eint2) => Eint1-Eint2 (EDGE-EDGE)
 *
 *                               \
 *                                *                ^
 *                                |Eint1          /
 *                                |              /
 *                        o-------+-----o      D12
 *                       /.Eint2..|      \
 *                       .........|
 *                      ..........*
 *                      ........./
 *
 */
void ContactUpdater_FACE_FACE::processIntsec_EDGE_EDGE(const long index)
{
    auto& oe     = mIntsec[index];
    auto heit1   = mHalfEdges1[edgeIndexA(oe)];
    auto heit2   = mHalfEdges2[edgeIndexB(oe)];  
    mInfo.mType1 = ContactPairInfo::FT_EDGE;
    mInfo.mEit1  = (*heit1)->edge();
    mInfo.mType2 = ContactPairInfo::FT_EDGE;
    mInfo.mEit2  = (*heit2)->edge();
}


/*
 *          - If Pint = IT_EDGE_VERTEX(Eint1, Vint2)
 *            There are many possible configurations of edges and vertices
 *            as follows. However, they all end up in the same update, 
 *            which is EDGE-VERTEX.
 *
 *              o  *       o  *        *  o       *o          *o  
 *               \ |      ..\ |    ....| /        ||TP     ...||
 *                \|      ...\|    ....|/         ||       ...||
 *                 oTP    ....o    ....o          |o       ...|o
 *                /|      .../|    .../|         /|        ../|
 *               / |      ../ |    ../ |        / |        ./ |
 *              o  *       o  *     o  *       o  *        o  *
 *
 *
 *              o  *       o  *      o  *        *o         *o
 *              .\ |        \ |      .\ |        ||TP    ...||
 *              ..\|         \|      ..\|        ||      ...||
 *              ...o          o      ...o        |o      ...|o
 *              ...|\         ||TP   ...|        ||TP    ...||
 *              ...| \        ||     ...|        ||      ...||
 *                 *  o       *o        *o       *o         *o
 *
 *                                    ==> Eint1-Vint2  (EDGE-VERTEX)
 */
void ContactUpdater_FACE_FACE::processIntsec_EDGE_VERTEX(const long index)
{

    auto& oe     = mIntsec[index];
    auto heit1   = mHalfEdges1[edgeIndexA(oe)];
    auto vit2    = mVertices2 [oe.mIndexB];  
    mInfo.mType1 = ContactPairInfo::FT_EDGE;
    mInfo.mEit1  = (*heit1)->edge();
    mInfo.mType2 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit2  = vit2;
}


void ContactUpdater_FACE_FACE::processIntsec_VERTEX_EDGE(const long index)
{

    auto& oe     = mIntsec[index];
    auto vit1    = mVertices1 [oe.mIndexA];
    auto heit2   = mHalfEdges2[edgeIndexB(oe)];
    mInfo.mType1 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit1  = vit1;
    mInfo.mType2 = ContactPairInfo::FT_EDGE;
    mInfo.mEit2  = (*heit2)->edge();
}


/*
 *         - If Pint = VERTEX-VERTEX (Vint1, Vint2)
 *            There are many possible configurations of edges and vertices
 *            as follows. However, they all end up in the same update, 
 *            which is VERTEX-VERTEX.
 *                                    ==> Vint1-Vint2  (VERTEX-VERTEX)
 *
 *              *          *          *         *          *o
 *               \          \          \         \         ..\
 *            *-_ \TP    o-_ \      o-_ \         \TP      ...\
 *              _-*o     .._-*o     .._-*o    *o==*o       *--*o
 *            o-  /      *-  /      o-  /         /           /
 *               /          /          /         /           /
 *              o          o          *         o           o
 *
 *              *o         *          *        *     *o
 *             ..\          \          \        \    .\\
 *             ...\          \          \        \   ..\\
 *             o--*o     *o==*oTP    o--*o    o--*o  ...*o
 *                /          /       .../     .../   ..//
 *               /          /        ../      ../    .//
 *              *          o          o*       *o    *o
 */
void ContactUpdater_FACE_FACE::processIntsec_VERTEX_VERTEX(const long index)
{
    auto& oe     = mIntsec[index];
    auto vit1    = mVertices1 [oe.mIndexA];
    auto vit2    = mVertices2 [oe.mIndexB];
    mInfo.mType1 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit1  = vit1;
    mInfo.mType2 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit2  = vit2;
}


/*          - If Pint = IT_VERTEX_INTERIOR (Vint1, FACE2)
 *                                             => Vint1-FACE2 (VERTEX-FACE)
 *
 *                              Vint1            ^
 *                       *------*               /
 *                      Eint1...|              /
 *                       .......|           D12
 *                       .......|Eint2      
 *                              *
 */
void ContactUpdater_FACE_FACE::processIntsec_VERTEX_INTERIOR(const long index)
{
    auto& oe     = mIntsec[index];
    auto vit1    = mVertices1 [oe.mIndexA];
    mInfo.mType1 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit1  = vit1;
}


void ContactUpdater_FACE_FACE::processIntsec_INTERIOR_VERTEX(const long index)
{
    auto& oe     = mIntsec[index];
    auto vit2    = mVertices2 [oe.mIndexB];
    mInfo.mType2 = ContactPairInfo::FT_VERTEX;
    mInfo.mVit2  = vit2;
}


/*
 *          - If Pint1 = IT_EDGE_EDGE (Eint11, Eint12)
 *            - If Pint2 = IT_EDGE_EDGE (Eint21, Eint22)
 *
 *              - Eint11==Eint21 on FACE1 => Eint11-FACE2 (EDGE-FACE)
 *
 *                D12 -->                      D12 -->
 *
 *                   o  *                       *  o
 *                ....\ |Eint11==Eint21         | /Eint12
 *               Eint12\|                       |/
 *                ......+                       +
 *                ......|\                     /|
 *                ......| \                   /.|
 *                ......|  o                 o..|
 *                ......|  |                 |..|Eint11==Eint21
 *                ......|  |                 |..|
 *                ......|  o                 o..|
 *                ......| /                   \.|
 *                ......|/                     \|
 *                ......+                       +
 *               Eint22/|                       |\
 *                ..../ |                       | \Eint22
 *                   o  *                       *  o
 *
 *              - Eint12==Eint22 on FACE2 => FACE1-Eint12 (FACE-EDGE)
 *
 *                D12 -->                      D12 -->
 *
 *                   *  o                       o  *
 *                ....\ |Eint12==Eint22         | /Eint11
 *               Eint11\|                       |/
 *                ......+                       +
 *                ......|\                     /|
 *                ......| \                   /.|
 *                ......|  *                 *..|
 *                ......|  |                 |..|Eint12==Eint22
 *                ......|  |                 |..|
 *                ......|  *                 *..|
 *                ......| /                   \.|
 *                ......|/                     \|
 *                ......+                       +
 *               Eint21/|                       |\
 *                ..../ |                       | \Eint21
 *                   *  o                       o  *
 */
void ContactUpdater_FACE_FACE::processIntsec_EDGE_EDGE__EDGE_EDGE(
    const long index1,
    const long index2
) {
    auto& oe1    = mIntsec[index1];
    auto& oe2    = mIntsec[index2];
    auto heit11  = mHalfEdges1[edgeIndexA(oe1)];
    auto heit21  = mHalfEdges1[edgeIndexA(oe2)];
    auto eit11   = (*heit11)->edge();
    auto eit21   = (*heit21)->edge();

    if (eit11==eit21) {
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = eit11;
    }
    else {
        auto heit12  = mHalfEdges2[edgeIndexB(oe1)];  
        auto heit22  = mHalfEdges2[edgeIndexB(oe2)];  
        auto eit12   = (*heit12)->edge();
        auto eit22   = (*heit22)->edge();
        if (eit12!=eit22) {
            log(ERROR, __FILE__, __LINE__, "Logic error");
        }
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = eit12;
    }
}


/*
 *          - If Pint1 = IT_EDGE_EDGE (Eint11, Eint12)
 *            - If Pint2 = IT_EDGE_VERTEX (Eint21, Vint22)
 *
 *              - Vint22 incident to Eint12 => FACE1-Eint12 (FACE-EDGE)
 *
 *                D12 -->                      D12 -->
 *
 *                   .* o                      o *
 *                   ..\|                      |/Eint11
 *                   ...+                      +
 *                   ...|\Eint11              /|
 *                Eint12| *                  *.|
 *                   ...| |                  |.|Eint12
 *                   ...| *                  *.|
 *                   ...|/                    \|
 *                   o--oVint22          Vint22o--o
 *                     /Eint21                  \Eint21
 *                    *                          *
 *
 *              - Eint11==Eint21 on FACE1 => Eint11-FACE2 (EDGE-FACE)
 *
 *                D12 -->                      D12 -->
 *
 *                   o *                       * o
 *                  ..\|                       |/Eint12
 *                  ...+                       +
 *                  ...|\Eint12               /|
 *               Eint11| \                   /.|Eint11==Eint21
 *             ==Eint21|  o                 o..|
 *                  ...|  |                 |..|
 *                  ...|  o                 o..|
 *                  ...| /                   \.|
 *                  ...|/                     \|
 *                  o--oVint22           Vint22o--o
 *                     |                       |
 *                     *                       *
 */
void ContactUpdater_FACE_FACE::processIntsec_EDGE_EDGE__EDGE_VERTEX(
    const long index1,
    const long index2
) {

    auto& oe1    = mIntsec[index1];
    auto& oe2    = mIntsec[index2];
    auto heit11  = mHalfEdges1[edgeIndexA(oe1)];
    auto heit12  = mHalfEdges2[edgeIndexB(oe1)];
    auto heit21  = mHalfEdges1[edgeIndexA(oe2)];
    auto vit22   = mVertices2 [oe2.mIndexB];

    auto eit11   = (*heit11)->edge();
    auto eit21   = (*heit21)->edge();

    if (eit11==eit21) {
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = eit11;
    }
    else {
        auto eit12   = (*heit12)->edge();
        if (!isVertexIncidentToEdge(vit22, heit12)) {
            log(ERROR, __FILE__, __LINE__, "Logic error");
        }
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = eit12;
    }
}



void ContactUpdater_FACE_FACE::processIntsec_EDGE_EDGE__VERTEX_EDGE(
    const long index1,
    const long index2
) {
    auto& oe1    = mIntsec[index1];
    auto& oe2    = mIntsec[index2];
    auto heit11  = mHalfEdges1[edgeIndexA(oe1)];
    auto heit12  = mHalfEdges2[edgeIndexB(oe1)];
    auto vit21   = mVertices1 [oe2.mIndexA];
    auto heit22  = mHalfEdges2[edgeIndexB(oe2)];
    auto eit12   = (*heit12)->edge();
    auto eit22   = (*heit22)->edge();

    if (eit12==eit22) {
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = eit12;
    }
    else {
        auto eit11   = (*heit11)->edge();
        if (!isVertexIncidentToEdge(vit21, heit11)) {
            log(ERROR, __FILE__, __LINE__, "Logic error");
        }
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = eit11;
    }
}


/*
 *          - If Pint1 = IT_EDGE_EDGE (Eint11, Eint12)
 *            - If Pint2 = IT_VERTEX__VERTEX (Vint21, Vint22)
 *
 *              - Vint21 is incident to Eint11 => Eint11-FACE2 (EDGE-FACE)
 *
 *                  D12 -->               D12 -->
 *
 *                   o  *                     *  o
 *              Eint12\ |Eint11         Eint11| /Eint12
 *                   ..\|                     |/
 *                   ...+                     +
 *                   ...|\                   /|
 *                   ...| \                 /.|
 *                   ...|  o               o..|
 *                   ...|  |               |..|
 *                   ...|  |               |..|
 *                   ...|  o               o..|
 *                   ...| /                 \.|
 *                   ...|/                   \|
 *                   *--*oVint21==     Vint21o*--*
 *                     /  Vint22     ==Vint22  \
 *                    o                         o
 *
 *
 *              - Vint22 is incident to Eint12 => FACE1-Eint12 (FACE-EDGE)
 *
 *                  D12 -->               D12 -->
 *
 *                   *  o                     o  *
 *              Eint11\ |Eint12         Eint12| /Eint11
 *                   ..\|                     |/
 *                   ...+                     +
 *                   ...|\                   /|
 *                   ...| \                 /.|
 *                   ...|  *               *..|
 *                   ...|  |               |..|
 *                   ...|  |               |..|
 *                   ...|  *               *..|
 *                   ...| /                 \.|
 *                   ...|/                   \|
 *                   o--*oVint21==     Vint21*o--o
 *                     /  Vint22      ==Vint22 \
 *                    *                         *
 *
 */
void ContactUpdater_FACE_FACE::processIntsec_EDGE_EDGE__VERTEX_VERTEX(
    const long index1,
    const long index2
) {
    auto& oe1    = mIntsec[index1];
    auto& oe2    = mIntsec[index2];
    auto heit11  = mHalfEdges1[edgeIndexA(oe1)];
    auto heit12  = mHalfEdges2[edgeIndexB(oe1)];
    auto vit21   = mVertices1 [oe2.mIndexA];
    auto vit22   = mVertices2 [oe2.mIndexB];

    auto eit11   = (*heit11)->edge();
    auto eit12   = (*heit12)->edge();

    if (isVertexIncidentToEdge(vit21, heit11)) {
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = eit11;
    }
    else{ 
        if (!isVertexIncidentToEdge(vit22, heit12)) {
            log(ERROR, __FILE__, __LINE__, "Logic error");
        }
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = eit12;
    }
}

/*
 *          - If Pint1 = IT_EDGE_EDGE (Eint11, Eint12)
 *            - If Pint2 = IT_VERTEX_INTERIOR (Vint21, *)
 *                                          => Eint11-FACE2 (EDGE-FACE)
 *              (Invariant: Vint21 is incident to Eint11.)
 *
 *                  D12 -->
 *
 *                 |..\
 *                 |...*Vint21
 *                 o...|
 *                  \..|
 *                   \.|Eint11
 *                    \|
 *                     +
 *                     |\
 *                     | \Eint12
 *                     *  o
 */
void ContactUpdater_FACE_FACE::processIntsec_EDGE_EDGE__VERTEX_INTERIOR(
    const long index1,
    const long index2
) {
    auto& oe1    = mIntsec[index1];
    auto& oe2    = mIntsec[index2];
    auto heit11  = mHalfEdges1[edgeIndexA(oe1)];
    auto vit21   = mVertices1 [oe2.mIndexA];

    auto eit11   = (*heit11)->edge();

    if (!isVertexIncidentToEdge(vit21, heit11)) {
        log(ERROR, __FILE__, __LINE__, "Logic error");
    }

    mInfo.mType1 = ContactPairInfo::FT_EDGE;
    mInfo.mEit1  = eit11;
}


void ContactUpdater_FACE_FACE::processIntsec_EDGE_EDGE__INTERIOR_VERTEX(
    const long index1,
    const long index2
) {
    auto& oe1    = mIntsec[index1];
    auto& oe2    = mIntsec[index2];
    auto heit12  = mHalfEdges2[edgeIndexB(oe1)];
    auto vit22   = mVertices2 [oe2.mIndexB];

    auto eit12   = (*heit12)->edge();

    if (!isVertexIncidentToEdge(vit22, heit12)) {
        log(ERROR, __FILE__, __LINE__, "Logic error");
    }

    mInfo.mType2 = ContactPairInfo::FT_EDGE;
    mInfo.mEit2  = eit12;
}


/*
 *          - If Pint1 = IT_EDGE_VERTEX (Eint11, Vint12)
 *
 *            - If Pint2 = IT_EDGE_VERTEX (Eint21, Vint22) 
 *
 *              - If (Vint12, Vint22) is on FACE2, 
 *
 *                - If Eint11 == Eint21  => Eint11-(Vint12,Vint22) (EDGE-EDGE)
 *                                  
 *                    D12 -->               D12 -->
 *
 *                       \                   /
 *                      o *               o *
 *                     ..\|                \|
 *                     ...o                 o
 *                     ...|                 |
 *                     ...o                 o
 *                     ../|                /|
 *                      o *               o *
 *                       /                   \
 *
 *                - If NOT ((Vint12,Vint22) //{eps} Eint11 ) &&
 *                     NOT ((Vint12,Vint22) //{eps} Eint21 )
 *
 *                                    => FACE1-(Vint12,Vint22) (FACE-EDGE)
 *
 *                    D12 -->               D12 <--
 *          
 *                     *  o                 *     o
 *                     .\ |                  \   /
 *                     ..\|                   \ /
 *                     ...oVint12              oVint12
 *                     ...|\Eint11             |\Eint11
 *                     ...| *                  |.*
 *                     ...| |                  |.|
 *                     ...| *                  |.*
 *                     ...|/Eint21             |/Eint21
 *                     ...oVint22              oVint22
 *                     ../|                   / \
 *                     ./ |                  /   \
 *                     *  o                 *     o
 *
 *                - If ((Vint12,Vint22) //{eps} Eint11 )
 *
 *                                    => Eint11-(Vint12,Vint22) (EDGE-EDGE)
 *
 *                    D12 -->               D12 -->
 *
 *                        \                 \
 *                    ..o  *                 *  o
 *                    ...\ |Eint11     Eint11| /
 *                  Vint12o|                 |oVint12
 *                    ....||                 ||
 *                    ....|*                 *|
 *                    ....||                 ||
 *                  Vint22o|                 |oVint22
 *                    .../ |Eint21     Eint21| \
 *                    ..o  *                 *  o
 *                        /                 /
 *
 *                - If ((Vint12,Vint22) //{eps} Eint12 )
 *
 *                                    => Eint12-(Vint12,Vint22) (EDGE-EDGE)
 *
 *                    D12 -->               D12 -->
 *
 *                        \                 \
 *                    ..o  *                 *  o
 *                    ...\ |Eint11     Eint11| /
 *                  Vint12o|                 |oVint12
 *                    ....||                 ||
 *                    ....|*                 *|
 *                    ....||                 ||
 *                  Vint22o|                 |oVint22
 *                    .../ |Eint21     Eint21| \
 *                    ..o  *                 *  o
 *                        /                 /
 *
 *              - If Eint11 == Eint21 && (Vint12, Vint22) is not on FACE2
 *
 *                                               => Eint11-FACE2 (EDGE-FACE)
 *                    D12 -->                     D12 --> 
 *                      \                           \
 *                    o  *                           *  o
 *                    .\ |Eint11==Eint21             | /
 *                    ..\|                           |/
 *                    ...oVint12                     o
 *                    ...|\                         /|
 *                    ...| o                       o.|
 *                    ...| |                       |.|
 *                    ...| o                       o.|
 *                    ...|/                         \|
 *                    ...oVint22                     o
 *                    ../|                           |\
 *                    ./ |                           | \
 *                    o  *                           *  o
 *                      /                           /
 */
void ContactUpdater_FACE_FACE::processIntsec_EDGE_VERTEX__EDGE_VERTEX(
    const long index1,
    const long index2
) {
    auto& oe1    = mIntsec[index1];
    auto& oe2    = mIntsec[index2];
    auto  heit11 = mHalfEdges1[edgeIndexA(oe1)];
    auto  heit21 = mHalfEdges1[edgeIndexA(oe2)];
    auto  eit11  = (*heit11)->edge();
    auto  eit21  = (*heit21)->edge();

    HalfEdgeIt heit2;
    if (findHalfEdgeBody2(oe1.mIndexB, oe2.mIndexB, heit2)) {
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = (*heit2)->edge();

        if (eit11 == eit21) {
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = eit11;

        }
        else if (areParallel((*heit2)->edge(), 2, eit11, 1)) {
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = eit11;

        }
        else if (areParallel((*heit2)->edge(), 2, eit21, 1)) {
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = eit21;
        }
    }
    else {
        if (eit11 != eit21) {
            log(ERROR, __FILE__, __LINE__, "Logic Error");
        }
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = eit11;

    }
}


/*
 *          - If Pint1 = IT_EDGE_VERTEX (Eint11, Vint12)
 *
 *            - If Pint2 = IT_VERTEX_EDGE (Vint21, Eint22)
 *
 *              - If Vint21 is incident to Eint11 &&
 *                   Vint12 is incident to Eint22
 *                (Invariant:  Eint11 //{eps} Eint22)
 *
 *                                               => Eint11-Eint22 (EDGE-EDGE)
 *
 *                  D12 -->             D12 -->
 *
 *                   \                   
 *                    *                   *  o
 *                ..\ |Eint11       Eint11| /
 *                ...\|                   |/
 *                ....oVint12             oVint12
 *                ....|                   ||
 *                ....|                   ||
 *              Eint11|                   ||
 *                ....|                   ||
 *                ....|                   ||
 *                ... *Vint21             *Vint21
 *                .../|                  /|
 *                ../ |Eint22           / |Eint22
 *                 *  o                *  o
 *
 *
 *              - If Vint21 is incident to Eint11
 *
 *                                               => Eint11-FACE2 (EDGE-FACE)
 *                  D12 -->             D12 -->
 *
 *                    *                   *
 *                    |                   |
 *                 o--oVint12       Vint12o--o
 *                ....|\                 /|
 *                ....| \               /.|
 *              Eint11|  o             o..|Eint11
 *                ....| /Eint22   Eint22\.|
 *                ....|/                 \|
 *                ... *Vint21       Vint21*
 *                .../|                   |\ 
 *                  o |                   | o
 *                    *                   *
 *
 *
 *              - If Vint12 is incident to Eint22
 *
 *                                               => FACE1-Eint22 (FACE-EDGE)
 *
 *                  D12 -->             D12 -->
 *
 *                    o                   o
 *                    |                   |
 *                 *--*Vint21       Vint21*--*
 *                ....|\                 /|
 *                ....| \               /.|
 *              Eint22|  *             *..|Eint22
 *                ....| /Eint11   Eint11\.|
 *                ....|/                 \|
 *                ... oVint12       Vint12o
 *                .../|                   |\ 
 *                  * |                   | *
 *                    o                   o
 */
void ContactUpdater_FACE_FACE::processIntsec_EDGE_VERTEX__VERTEX_EDGE(
    const long index1,
    const long index2
) {

    auto& oe1   = mIntsec[index1];
    auto& oe2   = mIntsec[index2];

    auto heit11 = mHalfEdges1[edgeIndexA(oe1)];
    auto vit12  = mVertices2 [oe1.mIndexB];
    auto vit21  = mVertices1 [oe2.mIndexA];
    auto heit22 = mHalfEdges2[edgeIndexB(oe2)];

    auto eit11  = (*heit11)->edge();
    auto eit22  = (*heit22)->edge();

    if (isVertexIncidentToEdge(vit21, heit11)) {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = eit11;

        if (isVertexIncidentToEdge(vit12, heit22)) {

            if (!areParallel(eit11, 1, eit22, 2)) {
                log(ERROR, __FILE__, __LINE__, "Logic Error");
            }

            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = eit22;
        }
    }
    else {
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = eit22;

    }
}


/*
 *          - If Pint1 = IT_EDGE_VERTEX (Eint11, Vint12)
 *
 *            - If Pint2 = IT_VERTEX_INTERIOR (Vint21, *)
 *
 *                                               => Eint11-FACE2 (EDGE-FACE)
 *                  D12 -->
 *
 *                  o  *
 *                 ..\ |Eint11
 *                 ...\|
 *                 ....oVint12
 *                 ....|\
 *                 ....| \
 *               Vint21*  \
 *                 .../    o
 *                   *
 */
void ContactUpdater_FACE_FACE::processIntsec_EDGE_VERTEX__VERTEX_INTERIOR(
    const long index1,
    const long index2
) {

    auto& oe1   = mIntsec[index1];
    auto heit11 = mHalfEdges1[edgeIndexA(oe1)];
    auto eit11  = (*heit11)->edge();

    mInfo.mType1 = ContactPairInfo::FT_EDGE;
    mInfo.mEit1  = eit11;

}


/*
 *          - If Pint1 = IT_EDGE_VERTEX (Eint11, Vint12)
 *
 *            - If Pint2 = IT_INTERIOR_VERTEX (*, Vint22)
 *              (Invariant: Edge (Vint12,Vint22) is on FACE2)
 *                                        => FACE1-(Vint12,Vint22) (FACE-EDGE)
 *                   D12 -->
 *
 *                   o     *             
 *                    \   /Eint11
 *                     \ /
 *                      oVint12
 *                     /||
 *                    /.||
 *                   /..||
 *                  /...oVint22
 *                 *  ./
 *                    /
 *                   o
 */
void ContactUpdater_FACE_FACE::processIntsec_EDGE_VERTEX__INTERIOR_VERTEX(
    const long index1,
    const long index2
) {

    auto& oe1  = mIntsec[index1];
    auto& oe2  = mIntsec[index2];

    HalfEdgeIt heit2;
    if (findHalfEdgeBody2(oe1.mIndexB, oe2.mIndexB, heit2)) {

        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = (*heit2)->edge();

    }
    else {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }

}


/*
 *          - If Pint1 = IT_EDGE_VERTEX (Eint11, Vint12)
 *
 *            - If Pint2 = IT_VERTEX_VERTEX (Vint21, Vint22)
 *
 *              - If Edge (Vint12,Vint22) is on FACE2 &&
 *                        Vint21 is incident to Eint11
 *                                      => Eint11-(Vint12,Vint22) (EDGE-EDGE)
 *
 *                 D12 -->                  D12 -->
 *
 *                     *  o                 o  *
 *               Eint11| /                  .\ |Eint11
 *                     |/                   ..\|
 *                     |oVint12             ..o|Vint12
 *                     ||                   ..||
 *                     ||                   ..||
 *                     ||                   ..||
 *                     *oVint21==Vint22     ..o*Vint21==Vint22
 *                    /  \                  ./|
 *                   *    o                 * o
 *
 *              - If Edge (Vint12,Vint22) is on FACE2
 *                                        => FACE1-(Vint12,Vint22) (FACE-EDGE)
 *
 *                   D12 <--               D12 -->
 *
 *                  *     o                  *
 *             Eint11\   /                    \Eint11
 *                    \ /                      \
 *               Vint12o                    o---oVint12
 *                    ||\                   ....|\
 *                    ||.\                  ....| \
 *                    ||..*                 ....|  *
 *                    ||..|                 ....|  |
 *                    ||..|                 ....|  |
 *                    ||..*                 ....|  *
 *                    ||./                  ....| /
 *                    ||/                   ....|/
 *                     o*Vint21==Vint22     ...o*Vint21==Vint22
 *                    / \                   ../|
 *                   /   \                  ./ |
 *                  *     o                 o  *
 *
 *              - If Vint21 is incident to Eint11
 *                                               => Eint11-FACE2 (EDGE-FACE)
 *
 *                   D12 -->            D12 -->
 *
 *                  o  *                       *  o
 *                  .\ |Eint11                 | /
 *                  ..\|                       |/
 *                  ...oVint12                 oVint12
 *                  ...|\                     /|
 *                  ...| \                   /.|
 *                  ...|  o                 o..|
 *                  ...|  |                 |..|Eint11
 *                  ...|  |                 |..|
 *                  ...|  o                 o..|
 *                  ...| /                   \.|
 *                  ...|/                     \|
 *                  o--*oVint21==Vint22        *oVint21==Vint22
 *                    /                       /
 *                   /                       /
 *                  *                       *
 */
void ContactUpdater_FACE_FACE::processIntsec_EDGE_VERTEX__VERTEX_VERTEX(
    const long index1,
    const long index2
) {
    auto& oe1   = mIntsec[index1];
    auto& oe2   = mIntsec[index2];
    auto heit11 = mHalfEdges1[edgeIndexA(oe1)];
    auto vit21  = mVertices1 [oe2.mIndexA];

    auto eit11  = (*heit11)->edge();

    HalfEdgeIt heit2;
    if (findHalfEdgeBody2(oe1.mIndexB, oe2.mIndexB, heit2)) {

        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = (*heit2)->edge();

        if (isVertexIncidentToEdge(vit21, heit11)) {

            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = eit11;
        }
    }
    else {
        if (!isVertexIncidentToEdge(vit21, heit11)) {
            log(ERROR, __FILE__, __LINE__, "Logic Error");
        }

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = eit11;
        
    }

}


void ContactUpdater_FACE_FACE::processIntsec_VERTEX_EDGE__VERTEX_EDGE(
    const long index1,
    const long index2
) {

    auto& oe1    = mIntsec[index1];
    auto& oe2    = mIntsec[index2];
    auto heit12  = mHalfEdges2[edgeIndexB(oe1)];
    auto heit22  = mHalfEdges2[edgeIndexB(oe2)];

    auto eit12   = (*heit12)->edge();
    auto eit22   = (*heit22)->edge();

    HalfEdgeIt heit1;
    if (findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*heit1)->edge();

        if (eit12 == eit22) {

            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = eit12;

        }
        else if (areParallel((*heit1)->edge(), 1, eit12, 2)) {

            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = eit12;

        }
        else if (areParallel((*heit1)->edge(), 1, eit22, 2)) {

            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = eit22;

        }
    }
    else {

        if (eit12 != eit22) {

            log(ERROR, __FILE__, __LINE__, "Logic Error");
        }
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = eit12;

    }
}


void ContactUpdater_FACE_FACE::processIntsec_VERTEX_EDGE__VERTEX_VERTEX(
    const long index1,
    const long index2
) {
    auto& oe1   = mIntsec[index1];
    auto& oe2   = mIntsec[index2];
    auto heit12 = mHalfEdges2[edgeIndexB(oe1)];
    auto vit22  = mVertices2 [oe2.mIndexB];

    auto eit12  = (*heit12)->edge();

    HalfEdgeIt heit1;
    if (findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*heit1)->edge();

        if (isVertexIncidentToEdge(vit22, heit12)) {

            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = eit12;
        }
    }
    else {
        if (!isVertexIncidentToEdge(vit22, heit12)) {
            log(ERROR, __FILE__, __LINE__, "Logic Error");
        }

        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = eit12;
        
    }
}


void ContactUpdater_FACE_FACE::processIntsec_VERTEX_EDGE__VERTEX_INTERIOR(
    const long index1,
    const long index2
) {
    auto& oe1  = mIntsec[index1];
    auto& oe2  = mIntsec[index2];

    HalfEdgeIt heit1;
    if (findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*heit1)->edge();

    }
    else {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }

}


void ContactUpdater_FACE_FACE::processIntsec_VERTEX_EDGE__INTERIOR_VERTEX(
    const long index1,
    const long index2
) {
    auto& oe1   = mIntsec[index1];
    auto heit12 = mHalfEdges2[edgeIndexB(oe1)];
    auto eit12  = (*heit12)->edge();

    mInfo.mType2 = ContactPairInfo::FT_EDGE;
    mInfo.mEit2  = eit12;

}


/*
 *          - If Pint1 = IT_VERTEX_VERTEX (Vint11, Vint12)
 *
 *            - If Pint2 = IT_VERTEX_VERTEX (Vint21, Vint22)
 *
 *                           ==> (Vint11,Vint21)-(Vint12,Vint22) (EDGE-EDGE)
 *
 *                D12 -->     D12 -->         D12 -->
 *
 *                *    o         *              o
 *                 \  /           \              \
 *                  \/             \              \
 *                  *o          o--o*         *---*o
 *                  ||          ...||          ...||
 *                  ||          ...||          ...||
 *                  ||          ...||          ...||
 *                  *o          o--o*         o---o*
 *                  /\             /              /
 *                 /  \           /              /
 *                *    o         *              *
 */
void ContactUpdater_FACE_FACE::processIntsec_VERTEX_VERTEX__VERTEX_VERTEX(
    const long index1,
    const long index2
) {
    auto& oe1  = mIntsec[index1];
    auto& oe2  = mIntsec[index2];

    HalfEdgeIt heit1;
    if (findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*heit1)->edge();

        HalfEdgeIt heit2;
        if (findHalfEdgeBody2(oe1.mIndexB, oe2.mIndexB, heit2)) {

            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = (*heit2)->edge();

        }
        else {
            log(ERROR, __FILE__, __LINE__, "Logic Error");
        }
    }
    else {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
}


/*
 *          - If Pint1 = IT_VERTEX_VERTEX (Vint11, Vint12)
 *            - If Pint2 = IT_VERTEX_INTERIOR (Vint21, *)
 *                (Invariant: Edge (Vint11,Vint21) is on FACE1)
 *                                     ==> (Vint11,Vint21)-FACE2  (EDGE-FACE)
 *
 *                  D12 -->
 *
 *                           o
 *                 ...\     /
 *                 ....*Vint21
 *                 ....|  /
 *                 ....| /
 *                 ....|/
 *                 o--*oVint11==Vint12
 *                    /
 *                   *
 */
void ContactUpdater_FACE_FACE::processIntsec_VERTEX_VERTEX__VERTEX_INTERIOR(
    const long index1,
    const long index2
) {
    auto& oe1  = mIntsec[index1];
    auto& oe2  = mIntsec[index2];
    HalfEdgeIt heit1;
    if (findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*heit1)->edge();

    }
    else {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
}


void ContactUpdater_FACE_FACE::processIntsec_VERTEX_VERTEX__INTERIOR_VERTEX(
    const long index1,
    const long index2
) {
    auto& oe1  = mIntsec[index1];
    auto& oe2  = mIntsec[index2];

    HalfEdgeIt heit2;
    if (findHalfEdgeBody2(oe1.mIndexB, oe2.mIndexB, heit2)) {

        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = (*heit2)->edge();

    }
    else {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
}


/*
 *          - If Pint1 = IT_VERTEX_INTERIOR (Vint11, *)
 *
 *            - If Pint2 = IT_VERTEX_INTERIOR (Vint21, *)
 *              (Invariant: Edge(Vint11,Vint21) is on FACE1)
 *                                     ==> (Vint11,Vint21)-FACE2  (EDGE-FACE)
 *               D12 -->
 *
 *                *
 *                .\
 *                ..\
 *                ...*Vint11
 *                ...|
 *                ...|
 *                ...|
 *                ...*Vint21
 *                ../
 *                ./
 *                *
 */
void ContactUpdater_FACE_FACE::processIntsec_VERTEX_INTERIOR__VERTEX_INTERIOR(
    const long index1,
    const long index2
) {
    auto& oe1  = mIntsec[index1];
    auto& oe2  = mIntsec[index2];

    HalfEdgeIt heit1;
    if (findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*heit1)->edge();

    }
    else {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
}


void ContactUpdater_FACE_FACE::processIntsec_INTERIOR_VERTEX__INTERIOR_VERTEX(
    const long index1,
    const long index2
) {
    auto& oe1  = mIntsec[index1];
    auto& oe2  = mIntsec[index2];

    HalfEdgeIt heit2;
    if (findHalfEdgeBody2(oe1.mIndexB, oe2.mIndexB, heit2)) {

        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = (*heit2)->edge();

    }
    else {
        log(ERROR, __FILE__, __LINE__, "Logic Error");
    }
}


void ContactUpdater_FACE_FACE::processIntsec_VERTEX_INTERIOR__INTERIOR_VERTEX(
    const long index1,
    const long index2
) {
    log(ERROR, __FILE__, __LINE__, "Logic Error");
}


}// namespace Makena


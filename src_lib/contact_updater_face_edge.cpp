#include "contact_updater_face_edge.hpp"

/**
 * @file contact_updater_face_edge.cpp
 *
 * @brief
 */
namespace Makena {


long ContactUpdater_FACE_EDGE::edgeIndex(IntSec2D::OutputElem& e)
{
    if (e.mIndexA==0&&e.mIndexAaux==mLastIndex) {
        return e.mIndexAaux;
    }
    else {
        return e.mIndexA;  
    }
}


bool ContactUpdater_FACE_EDGE::isVertexIncidentToEdge(VertexIt v, HalfEdgeIt h)
{
    return (*h)->src()==v ||(*h)->dst()==v;
}


bool ContactUpdater_FACE_EDGE::findHalfEdgeBody1(
    const long  index1,
    const long  index2,
    HalfEdgeIt& heOut
) {
    auto numEdges = mHalfEdges.size();
    auto heit1 = mHalfEdges[(index1+numEdges-1)%numEdges];
    auto heit2 = mHalfEdges[index1];
    auto vit1  = mVertices[index1];
    auto vit2  = mVertices[index2];

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


ContactUpdater_FACE_EDGE::ContactUpdater_FACE_EDGE(
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
     mEpsilonAngle(epsilonAngle),
     mFaceOnBody1(mInfo.mType1==ContactPairInfo::FT_FACE)
     {;}


ContactUpdater_FACE_EDGE::~ContactUpdater_FACE_EDGE(){;}


bool ContactUpdater_FACE_EDGE:: update()
{
    orientWorld();
    bool intersectionFound = findIntersection2D();

    if (intersectionFound) {

        if (mTiltedTowardVertex1||mTiltedTowardVertex2) {
            auto index = findExtremeFeature2D();
            dispatchAndUpdate(index);
        }
        else {
            checkAndUpdateEdgeCases();
        }

        return false;
    }

    else {
        return true;

    }
}


void ContactUpdater_FACE_EDGE::checkAndUpdateEdgeCases()
{
    if (mIntsec.size()==1) {
        auto& oe = mIntsec[0];
        switch (oe.mType) {
          case IntSec2D::IT_VERTEX_VERTEX:
            updateNonTilting_VERTEX_VERTEX();
            break;
          case IntSec2D::IT_VERTEX_EDGE:
            updateNonTilting_VERTEX_EDGE();
            break;
          case IntSec2D::IT_EDGE_VERTEX:
            updateNonTilting_EDGE_VERTEX();
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


void ContactUpdater_FACE_EDGE::updateNonTilting_VERTEX_VERTEX()
{
    auto& oe     = mIntsec[0];
    auto  vit1   = mVertices[oe.mIndexA];
    auto  vit2   = (oe.mIndexB==1)?mVit1:mVit2;
    mInfo.mType1 = ContactPairInfo::FT_VERTEX;
    mInfo.mType2 = ContactPairInfo::FT_VERTEX;

    if (mFaceOnBody1) {
        mInfo.mVit1  = vit1;
        mInfo.mVit2  = vit2;
    }
    else {
        mInfo.mVit1  = vit2;
        mInfo.mVit2  = vit1;
    }
}


void ContactUpdater_FACE_EDGE::updateNonTilting_VERTEX_EDGE()
{
    auto& oe  = mIntsec[0];
    auto vit1 = mVertices[oe.mIndexA];

    if (mFaceOnBody1) {
        mInfo.mType1 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit1  = vit1;
    }
    else {
        mInfo.mType2 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit2  = vit1;
    }
}


void ContactUpdater_FACE_EDGE::updateNonTilting_EDGE_VERTEX()
{
    auto& oe    = mIntsec[0];

    auto  heit1 = mHalfEdges[edgeIndex(oe)];
    auto  vit2  = (oe.mIndexB==1)?mVit1:mVit2;

    if (mFaceOnBody1) {
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*heit1)->edge();
        mInfo.mType2 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit2  = vit2;
    }
    else {
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = (*heit1)->edge();
        mInfo.mType1 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit1  = vit2;
    }
}


void ContactUpdater_FACE_EDGE::updateNonTilting_EDGE_VERTEX_EDGE_VERTEX()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto heit11  = mHalfEdges[edgeIndex(oe1)];
    auto heit21  = mHalfEdges[edgeIndex(oe2)];

    if ((*heit11)->edge()==(*heit21)->edge()) {

        if (mFaceOnBody1) {   
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = (*heit11)->edge();
        }
        else {
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = (*heit11)->edge();
        }
    }
}


void ContactUpdater_FACE_EDGE::updateNonTilting_EDGE_VERTEX_VERTEX_EDGE()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto heit11 = mHalfEdges[edgeIndex(oe1)];
    auto vit21  = mVertices [oe2.mIndexA];

    if (isVertexIncidentToEdge(vit21, heit11)) {

        if (mFaceOnBody1) {   
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = (*heit11)->edge();
        }
        else {
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = (*heit11)->edge();
        }
    }
}


void ContactUpdater_FACE_EDGE::updateNonTilting_EDGE_VERTEX_VERTEX_VERTEX()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto heit11 = mHalfEdges[edgeIndex(oe1)];
    auto vit21  = mVertices [oe2.mIndexA];

    if (isVertexIncidentToEdge(vit21, heit11)) {

        if (mFaceOnBody1) {   
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = (*heit11)->edge();
        }
        else {
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = (*heit11)->edge();
        }
    }
}


void ContactUpdater_FACE_EDGE::updateNonTilting_VERTEX_EDGE_EDGE_VERTEX()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto vit11  = mVertices [oe1.mIndexA];
    auto heit21 = mHalfEdges[edgeIndex(oe2)];

    if (isVertexIncidentToEdge(vit11, heit21)) {

        if (mFaceOnBody1) {   
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = (*heit21)->edge();
        }
        else {
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = (*heit21)->edge();
        }
    }
}


void ContactUpdater_FACE_EDGE::updateNonTilting_VERTEX_EDGE_VERTEX_EDGE()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    HalfEdgeIt heit1;
    if (findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        if (mFaceOnBody1) {   
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = (*heit1)->edge();
        }
        else {
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = (*heit1)->edge();
        }
    }
}


void ContactUpdater_FACE_EDGE::updateNonTilting_VERTEX_EDGE_VERTEX_VERTEX()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    HalfEdgeIt heit1;
    if (findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        if (mFaceOnBody1) {   
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = (*heit1)->edge();
        }
        else {
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = (*heit1)->edge();
        }
    }
}


void ContactUpdater_FACE_EDGE::updateNonTilting_VERTEX_VERTEX_EDGE_VERTEX()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    auto vit11  = mVertices [oe1.mIndexA];
    auto heit21 = mHalfEdges[edgeIndex(oe2)];

    if (isVertexIncidentToEdge(vit11, heit21)) {

        if (mFaceOnBody1) {   
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = (*heit21)->edge();
        }
        else {
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = (*heit21)->edge();
        }
    }
}


void ContactUpdater_FACE_EDGE::updateNonTilting_VERTEX_VERTEX_VERTEX_EDGE()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    HalfEdgeIt heit1;
    if (findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        if (mFaceOnBody1) {   
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = (*heit1)->edge();
        }
        else {
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = (*heit1)->edge();
        }
    }
}


void ContactUpdater_FACE_EDGE::updateNonTilting_VERTEX_VERTEX_VERTEX_VERTEX()
{
    auto& oe1 = mIntsec[0];
    auto& oe2 = mIntsec[1];

    HalfEdgeIt heit1;
    if (findHalfEdgeBody1(oe1.mIndexA, oe2.mIndexA, heit1)) {

        if (mFaceOnBody1) {   
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = (*heit1)->edge();
        }
        else {
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = (*heit1)->edge();
        }
    }
}



void ContactUpdater_FACE_EDGE::orientWorld()
{
    checkForTilt();

    rotateWorld(mInfo.mContactNormal1To2);
}


void ContactUpdater_FACE_EDGE::rotateWorld(const Vec3& zDir)
{
    Mat3x3 rMat = rotMatAlignZDirToZAxis(zDir);

    mRotMat1 = rMat * mBody1.Qmat();
    mCom1    = rMat * mBody1.CoM();
    mRotMat2 = rMat * mBody2.Qmat();
    mCom2    = rMat * mBody2.CoM();

}


Mat3x3 ContactUpdater_FACE_EDGE::rotMatAlignZDirToZAxis(const Vec3& zDir)
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


void ContactUpdater_FACE_EDGE::checkForTilt() {

    Vec3 n, p1, p2;

    if (mFaceOnBody1) {
        auto fit1 = mInfo.mFit1;
        n         = (*fit1)->nGCS(mBody1.Qmat());
        auto eit2 = mInfo.mEit2;
        auto heit = (*eit2)->he1();
        auto vit1 = (*heit)->src();
        auto vit2 = (*heit)->dst();
        p1        = (*vit1)->pGCS(mBody2.Qmat(), mBody2.CoM());
        p2        = (*vit2)->pGCS(mBody2.Qmat(), mBody2.CoM());
    }
    else {
        auto fit2 = mInfo.mFit2;
        n         = (*fit2)->nGCS(mBody2.Qmat());
        auto eit1 = mInfo.mEit1;
        auto heit = (*eit1)->he1();
        auto vit1 = (*heit)->src();
        auto vit2 = (*heit)->dst();
        p1        = (*vit1)->pGCS(mBody1.Qmat(), mBody1.CoM());
        p2        = (*vit2)->pGCS(mBody1.Qmat(), mBody1.CoM());
    }

    auto v12  = p2 - p1;

    v12.normalize();
    double dot = v12.dot(n);
    if (dot <= -2.0 * mEpsilonAngle) {
        mTiltedTowardVertex1 = false;
        mTiltedTowardVertex2 = true;
    }
    else if (dot >=  2.0 * mEpsilonAngle) {
        mTiltedTowardVertex1 = true;
        mTiltedTowardVertex2 = false;
    }
    else {
        mTiltedTowardVertex1 = false;
        mTiltedTowardVertex2 = false;
    }
}


bool ContactUpdater_FACE_EDGE::findIntersection2D()
{
    // Prepare InputElems
    FaceIt fit;
    EdgeIt eit;

    if (mFaceOnBody1) {
        fit = mInfo.mFit1;
        eit = mInfo.mEit2;
    }
    else {
        fit = mInfo.mFit2;
        eit = mInfo.mEit1;
    }

    vector< IntSec2D::InputElem> inputElem1;
    vector< IntSec2D::InputElem> inputElem2;

    mVertices.clear();
    mHalfEdges.clear();
    mIntsec.clear();

    /* 
     *         FACE  (CCW)
     *
     *          i+1     i
     *           *<-----* src
     *          /        \
     *         /    N/    \
     *         \    /     /
     *          \  v     /
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
    mLastIndex = ((*fit)->halfEdges().size()) - 1 ;
    long index = 0;
    for (auto he : (*fit)->halfEdges()) {
        auto vit = (*he)->src();
        mHalfEdges.push_back(he);
        mVertices.push_back(vit);
        Vec3 p = mFaceOnBody1?
                     ((*vit)->pGCS(mRotMat1, mCom1)):
                     ((*vit)->pGCS(mRotMat2, mCom2));
        Vec2 p2d(p.x(), p.y());
        IntSec2D::InputElem ie(p2d, index);
        inputElem1.push_back(ie);
        index++;
    }

    auto heit = (*eit)->he1();
    mVit1     = (*heit)->src();
    mVit2     = (*heit)->dst();
    Vec3 p1   = mFaceOnBody1?
                     ((*mVit1)->pGCS(mRotMat2, mCom2)):
                     ((*mVit1)->pGCS(mRotMat1, mCom1));

    Vec3 p2   = mFaceOnBody1?
                     ((*mVit2)->pGCS(mRotMat2, mCom2)):
                     ((*mVit2)->pGCS(mRotMat1, mCom1));

    Vec2 p1_2d(p1.x(), p1.y());
    Vec2 p2_2d(p2.x(), p2.y());
    IntSec2D::InputElem ie1(p1_2d, 1);
    inputElem2.push_back(ie1);
    IntSec2D::InputElem ie2(p2_2d, 2);
    inputElem2.push_back(ie2);

    mTiltingDir    = mTiltedTowardVertex1 ? (p1_2d - p2_2d) : (p2_2d - p1_2d);
    mTiltingVertex = mTiltedTowardVertex1 ? mVit1 : mVit2;

    IntSec2D IntsecFinder2D(
                      mEpsilonZero, mEpsilonZero, mEpsilonAngle, mLogStream);
    
    auto res=IntsecFinder2D.findIntersection(inputElem1, inputElem2, mIntsec);
    return res;

}


long ContactUpdater_FACE_EDGE::findExtremeFeature2D()
{
    double maxDist  = mTiltingDir.dot(mIntsec[0].mP);
    long   maxIndex = 0;

    const  long numOutElems = mIntsec.size();

    // Find the points in the leaning direction.
    for (long i = 1; i < numOutElems; i++) {
        double curDist = mTiltingDir.dot(mIntsec[i].mP);
        if (curDist > maxDist) {
            maxDist  = curDist;
            maxIndex = i;
        }
    }
    return maxIndex;
}


void ContactUpdater_FACE_EDGE::dispatchAndUpdate(const long index)
{

    auto& oe = mIntsec[index];
    switch(oe.mType) {

      case IntSec2D::IT_EDGE_EDGE:
        processIntsec_EDGE_EDGE(index);
        break;

      case IntSec2D::IT_EDGE_VERTEX:
        processIntsec_EDGE_VERTEX(index);
        break;

      case IntSec2D::IT_VERTEX_EDGE:
        processIntsec_VERTEX_EDGE(index);
        break;

      case IntSec2D::IT_VERTEX_INTERIOR:
        processIntsec_VERTEX_INTERIOR(index);
        break;

      case IntSec2D::IT_INTERIOR_VERTEX:
        processIntsec_INTERIOR_VERTEX(index);
        break;

      case IntSec2D::IT_VERTEX_VERTEX:
        processIntsec_VERTEX_VERTEX(index);
        break;

      default:
        break;
    }
}


/*
 *          - If Pint = IT_EDGE_EDGE (Eint1,Eint2) => Eint1-EDGE2 (EDGE-EDGE)
 *
 *                         ......\
 *                        ........*
 *                        ........|Eint1
 *                        ........|
 *                        o-------+-----o   D12 ==>
 *                        .EDGE2..|
 *                       .........|
 *                      ..........*
 *                      ........./
 */
void ContactUpdater_FACE_EDGE::processIntsec_EDGE_EDGE(const long index)
{
    auto& oe  = mIntsec[index];
    auto heit = mHalfEdges[edgeIndex(oe)];
    if (mFaceOnBody1) {
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*heit)->edge();
    }
    else {
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = (*heit)->edge();
    }
}


/*
 *          - If Pint = IT_EDGE_VERTEX(Eint1, Vint2)
 *            There are many possible configurations of edges and vertices
 *            as follows. However, they all end up in the same update,
 *            which is EDGE-VERTEX.
 *
 *                                    ==> Eint1-Vint2  (EDGE-VERTEX)
 *
 *
 *                              *          \............/
 *                              |           \..Eint1.../
 *                       o------oVint2    o--*-----o--*   D12 ==>
 *                        EDGE2 |                Vint2
 *                              |Eint1
 *                              |
 *                              *
 */
void ContactUpdater_FACE_EDGE::processIntsec_EDGE_VERTEX(const long index)
{

    auto& oe  = mIntsec[index];
    auto heit = mHalfEdges[edgeIndex(oe)];
    if (mFaceOnBody1) {

        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1  = (*heit)->edge();
        mInfo.mType2 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit2  = mTiltingVertex;

    }
    else {

        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2  = (*heit)->edge();
        mInfo.mType1 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit1  = mTiltingVertex;

    }
}


/*
 *         - If Pint = IT_VERTEX_EDGE(Vint1, Eint2)
 *
 *                                    ==> Vint1-EDGE2  (VERTEX-EDGE)
 *
 *                          *
 *                       ....\
 *                       .....\Vint1
 *                       o-----*----o   D12 ==>
 *                       ...../ Eint2
 *                       ..../
 *                          *
 */
void ContactUpdater_FACE_EDGE::processIntsec_VERTEX_EDGE(const long index)
{

    auto& oe = mIntsec[index];
    auto vit = mVertices [oe.mIndexA];
    if (mFaceOnBody1) {
        mInfo.mType1 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit1  = vit;
    }
    else {
        mInfo.mType2 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit2  = vit;
    }
}


/*
 *         - If Pint = VERTEX-VERTEX (Vint1, Vint2)
 *
 *                                    ==> Vint1-Vint2  (VERTEX-VERTEX) *
 *
 *                 *
 *                  \                      \.........../
 *                   \                      \........./
 *          o--------o*Vint1==Vint2      o---*------o*Vint1==Vint2
 *            EDGE2  /
 *                  /                     D12 ==>
 *                 *
 */
void ContactUpdater_FACE_EDGE::processIntsec_VERTEX_VERTEX(const long index)
{
    auto& oe = mIntsec[index];
    auto vit1 = mVertices [oe.mIndexA];
    auto vit2  = (oe.mIndexB==1)?mVit1:mVit2;
    if (mFaceOnBody1) {
        mInfo.mType1 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit1  = vit1;
        mInfo.mType2 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit2  = vit2;
    }
    else {
        mInfo.mType2 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit2  = vit1;
        mInfo.mType1 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit1  = vit2;
    }
}


/*
 *          - If Pint = IT_VERTEX_INTERIOR (Vint1, EDGE2)
 *                                             => Vint1-FACE2 (VERTEX-EDGE)
 *
 *                       .....\
 *                       ......\ Vint1
 *                       o------*-------o   D12 ==>
 *                      EDGE2../
 *                       ...../
 *                       ..../
 */
void ContactUpdater_FACE_EDGE::processIntsec_VERTEX_INTERIOR(const long index)
{
    auto& oe = mIntsec[index];
    auto vit = mVertices[oe.mIndexA];
    if (mFaceOnBody1) {
        mInfo.mType1 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit1  = vit;
    }
    else {
        mInfo.mType2 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit2  = vit;
    }
}


/*
 *          - If Pint = IT_INTERIOR_VERTEX (FACE1, Vint2)
 *                                             => FACE1-Vint2 (FACE-VERTEX)
 *
 *                     ..............|
 *                     .......Vint2..|
 *                     ..o------o....|  D12 ==>
 *                     ...EDGE2......|
 *                     ..............|
 *                     ..............|
 */
void ContactUpdater_FACE_EDGE::processIntsec_INTERIOR_VERTEX(const long index)
{
    if (mFaceOnBody1) {
        mInfo.mType2 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit2  = mTiltingVertex;
    }
    else {
        mInfo.mType1 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit1  = mTiltingVertex;
    }
}



}// namespace Makena


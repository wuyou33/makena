#include "contact_updater_face_vertex.hpp"

/**
 * @file contact_updater_face_vertex.cpp
 *
 * @brief
 */
namespace Makena {


long ContactUpdater_FACE_VERTEX::edgeIndex(IntSec2D::OutputElem& e)
{
    if (e.mIndexA==0&&e.mIndexAaux==mLastIndex) {
        return e.mIndexAaux;
    }
    else {
        return e.mIndexA;  
    }
}


bool ContactUpdater_FACE_VERTEX::findHalfEdgeBody1(
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


ContactUpdater_FACE_VERTEX::ContactUpdater_FACE_VERTEX(
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


ContactUpdater_FACE_VERTEX::~ContactUpdater_FACE_VERTEX(){;}


bool ContactUpdater_FACE_VERTEX::update()
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


void ContactUpdater_FACE_VERTEX::rotateWorld(const Vec3& zDir)
{
    Mat3x3 rMat = rotMatAlignZDirToZAxis(zDir);

    mRotMat1 = rMat * mBody1.Qmat();
    mCom1    = rMat * mBody1.CoM();
    mRotMat2 = rMat * mBody2.Qmat();
    mCom2    = rMat * mBody2.CoM();

}


Mat3x3 ContactUpdater_FACE_VERTEX::rotMatAlignZDirToZAxis(const Vec3& zDir)
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


bool ContactUpdater_FACE_VERTEX::findIntersection2D()
{
    // Prepare InputElems
    FaceIt   fit;
    VertexIt vit;

    if (mFaceOnBody1) {
        fit = mInfo.mFit1;
        vit = mInfo.mVit2;
    }
    else {
        fit = mInfo.mFit2;
        vit = mInfo.mVit1;
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

    Vec3 p = mFaceOnBody1?
                     ((*vit)->pGCS(mRotMat2, mCom2)):
                     ((*vit)->pGCS(mRotMat1, mCom1));

    Vec2 p_2d(p.x(), p.y());
    IntSec2D::InputElem ie1(p_2d, 0);
    inputElem2.push_back(ie1);

    IntSec2D IntsecFinder2D(
                      mEpsilonZero, mEpsilonZero, mEpsilonAngle, mLogStream);

    return IntsecFinder2D.findIntersection(inputElem1, inputElem2, mIntsec);

}


void ContactUpdater_FACE_VERTEX::dispatchAndUpdate()
{
    auto& oe = mIntsec[0];
    switch(oe.mType) {

      case IntSec2D::IT_EDGE_VERTEX:
        processIntsec_EDGE_VERTEX();
        break;

      case IntSec2D::IT_VERTEX_VERTEX:
        processIntsec_VERTEX_VERTEX();
        break;

      default:
        break;
    }
}


/*
 *          - If Pint = IT_EDGE_VERTEX (Eint1,Vint2)
 *                                       => Eint1-Vint2 (EDGE-VERTEX)
 *
 *                         ......\
 *                        ........*
 *                        ........|Eint1
 *                        ........|
 *                        ........oVit2
 *                       .........|
 *                       .........|
 *                      ..........*
 *                      ........./
 */
void ContactUpdater_FACE_VERTEX::processIntsec_EDGE_VERTEX()
{

    auto& oe  = mIntsec[0];
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
 *          - If Pint = IT_VERTEX_VERTEX (Vint1, Vint2)
 *                                             => Vint1-Vint2 (VERTEX-VERTEX)
 *
 *                       .....\
 *                       ......\ Vint1==Vint2
 *                       ......*o
 *                      ......./
 *                       ...../
 *                       ..../
 */
void ContactUpdater_FACE_VERTEX::processIntsec_VERTEX_VERTEX()
{
    auto& oe = mIntsec[0];
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


}// namespace Makena


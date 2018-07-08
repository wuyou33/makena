#include "bd_boundary_simplex_finder.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file bd_boundary_simplex_finder.cpp
 *
 * @brief this finds the minimal boundary triangle (or in rare cases edge or
 *        vertex)  of the A-B convex manifold that intersects the ray from
 *        the origin along the exit direction.
 *        It assumes a simplex of the A-B convex manifold that encloses the
 *        origin has been found by GJKOriginFounder. It is usually a 
 *        tetrahedron but it can be a triangle, an edge, or a vertex.
 *       
 *        First, we find the initial triangle as the starting point.
 *        If we can't find a triangle, this means the current edge or the 
 *        vertex is already the minimal boundary feature that encloses the ray.
 *        We return this the boundary simplex.
 *        Otherwise, we try to find the new vertex toward the normal of the 
 *        triangle. Of two possible directions of the normal, we choose the
 *        one along the exit direction.
 *        If we can't find a vertex, then the current triangle is already
 *        the minimal boundary simplex. We return this as the boundary.
 *        Otherwise, the new vertex together with each of the 3 edges of
 *        the triangle introduces 3 new triangles. We check which one of those
 *        intersects the ray. We pick that one, and discard the rest.
 *        We repeat the process above until we fail to find the new vertex
 *        toward the normal.
 * 
 *        Implementation note: To simplify the tests and pivoting we rotate
 *        the two rigit bodies and the A-B convex manifold such that the 
 *        exit direction coincides with the z-axis. This reduces the 
 *        intersection tests to (x,y) 2 dimensional tests against the origin.
 */
namespace Makena {

using namespace std;


BDBoundarySimplexFinder::BDBoundarySimplexFinder(
    ConvexRigidBody&         body1,
    ConvexRigidBody&         body2,
    vector<BDVertex>&        Q,
    const long&              maxNumPivots,
    const long&              maxNumCycles,
    const double&            epsilonZero,
    bool                     useGeomConfigTemp,
    const double&            scaling,
    std::ostream&            logStream
)
    :Loggable(logStream),
     mBody1(body1),
     mBody2(body2),
     mQ(Q),
     mMaxNumPivots(maxNumPivots),
     mMaxNumCycles(maxNumCycles),
     mEpsilonZero(epsilonZero),
     mUseGeomConfigTemp(useGeomConfigTemp),
     mScaling(scaling)
     {;}


BDBoundarySimplexFinder::~BDBoundarySimplexFinder(){;}


bool BDBoundarySimplexFinder::find(const Vec3& Zdir)
{
    mZdir = Zdir;
    rotateWorld();

    auto found = findInitial2Simplex();
    if (!found) {
        return true;
    }
    double maxDepth = averageDepth();
    long numCycles = 0;
    for (long numIter = 0; numIter < mMaxNumPivots; numIter++) {

        auto dir = findSupportDirectionOfTriangle();

        VertexIt vit1, vit2;

        bool found=findExtremePointsAlongSupportDirection(dir, vit1, vit2);

        if (!found) {
            return true;
        }

        selectNewTriangle(vit1, vit2);

        auto curDepth =averageDepth();
        if (curDepth <= maxDepth) {
            if (numCycles++ >= mMaxNumCycles) {
                log(WARNING, __FILE__, __LINE__, 
                     "Maximum number of cycles exceeded [%ld]", mMaxNumCycles);
                logQ(WARNING);
                return false;
            }
        }
    }

    log(WARNING, __FILE__, __LINE__, 
                "Maximum number of iteration exceeded [%ld]", mMaxNumPivots);
    logQ(WARNING);

    return false;
}


double BDBoundarySimplexFinder::averageDepth()
{
    const auto& p1 = mQ[0].p();
    const auto& p2 = mQ[1].p();
    const auto& p3 = mQ[2].p();
    return (p1.z() + p2.z() + p3.z()) / 3.0;
}


/** @brief this is basically the normal of the triangle oriented toward
 *         the departing direction.
 *         If the triangle is degenerate, take the departing direction.
 */
Vec3 BDBoundarySimplexFinder::findSupportDirectionOfTriangle()
{
    const auto& p1 = mQ[0].p();
    const auto& p2 = mQ[1].p();
    const auto& p3 = mQ[2].p();
    const Vec3 v12 = p2 - p1;
    const Vec3 v13 = p3 - p1;
    const Vec3 v12_cr_v13 = v12.cross(v13);

    if (v12_cr_v13.squaredNorm2()< mEpsilonZero) {
        return Vec3(0.0, 0.0, 1.0);
    }
    else {
        if (v12_cr_v13.dot(Vec3(0.0, 0.0, 1.0)) < 0.0) {
            return v12_cr_v13 * -1.0;
        }
        else {
            return v12_cr_v13;
        }
    }    
}


void BDBoundarySimplexFinder::rotateWorld()
{
    // Find the new coordinate axes with zDir for new z-axis.
    Vec3 r3 = mZdir;
    r3.normalize();
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
    const Vec3 r2 = r3.cross(r1);

    // Find the rotation matrix to rotate a point in GCS into the new CS.
    Mat3x3 rMat(r1, r2, r3);
    rMat.transposeInPlace();
    if (mUseGeomConfigTemp) {
        mRotMat1 = rMat * mBody1.QmatTemp();
        mCoM1    = rMat * mBody1.CoMTemp();
        mRotMat2 = rMat * mBody2.QmatTemp();
        mCoM2    = rMat * mBody2.CoMTemp();

        for (auto& q : mQ) {
            q.updateP( mBody1.ConvexHull(),  mRotMat1,  mCoM1,
                       mBody2.ConvexHull(),  mRotMat2,  mCoM2 );
        }
    }
    else {
        mRotMat1 = rMat * mBody1.Qmat();
        mCoM1    = rMat * mBody1.CoM();
        mRotMat2 = rMat * mBody2.Qmat();
        mCoM2    = rMat * mBody2.CoM();

        for (auto& q : mQ) {
            q.updateP( mBody1.ConvexHull(), mScaling, mRotMat1,  mCoM1,
                       mBody2.ConvexHull(), mScaling, mRotMat2,  mCoM2 );
        }
    }
}


bool BDBoundarySimplexFinder::isOriginInsideTriangleXY(
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& p3
) {

    const Vec3 v12 = p2 - p1;
    const Vec2 v12perp(v12.y(), -1.0*v12.x());// Rotate v12 CW.
    const Vec3 v23 = p3 - p2;
    const Vec2 v23perp(v23.y(), -1.0*v23.x());// Rotate v23 CW.
    const Vec3 v31 = p1 - p3;
    const Vec2 v31perp(v31.y(), -1.0*v31.x());// Rotate v31 CW.
    const Vec2 v1t(-1.0*p1.x(), -1.0*p1.y());
    const Vec2 v2t(-1.0*p2.x(), -1.0*p2.y());
    const Vec2 v3t(-1.0*p3.x(), -1.0*p3.y());

    bool isInside12, isInside23, isInside31;
    if (v12.x() * v31.y() - v31.x() * v12.y() <= 0.0) {
        // (1,2,3) is in CCW looking at it from positive z.
        isInside12 = v12perp.dot(v1t) <= 0.0;
        isInside23 = v23perp.dot(v2t) <= 0.0;
        isInside31 = v31perp.dot(v3t) <= 0.0;

    }
    else {
        // (1,2,3) is in CW looking at it from positive z.
        isInside12 = v12perp.dot(v1t) >= 0.0;
        isInside23 = v23perp.dot(v2t) >= 0.0;
        isInside31 = v31perp.dot(v3t) >= 0.0;

    }

    return isInside12 && isInside23 && isInside31;
}


bool BDBoundarySimplexFinder::isOriginInsideTriangleXYrelaxed(
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& p3
) {

    const Vec3 v12 = p2 - p1;
    const Vec2 v12perp(v12.y(), -1.0*v12.x());// Rotate v12 CW.
    const Vec3 v23 = p3 - p2;
    const Vec2 v23perp(v23.y(), -1.0*v23.x());// Rotate v23 CW.
    const Vec3 v31 = p1 - p3;
    const Vec2 v31perp(v31.y(), -1.0*v31.x());// Rotate v31 CW.
    const Vec2 v1t(-1.0*p1.x(), -1.0*p1.y());
    const Vec2 v2t(-1.0*p2.x(), -1.0*p2.y());
    const Vec2 v3t(-1.0*p3.x(), -1.0*p3.y());

    bool isInside12, isInside23, isInside31;
    if (v12.x() * v31.y() - v31.x() * v12.y() <= 0.0) {
        // (1,2,3) is in CCW looking at it from positive z.
        isInside12 = v12perp.dot(v1t) <= mEpsilonZero;
        isInside23 = v23perp.dot(v2t) <= mEpsilonZero;
        isInside31 = v31perp.dot(v3t) <= mEpsilonZero;

    }
    else {
        // (1,2,3) is in CW looking at it from positive z.
        isInside12 = v12perp.dot(v1t) >= -1.0 * mEpsilonZero;
        isInside23 = v23perp.dot(v2t) >= -1.0 * mEpsilonZero;
        isInside31 = v31perp.dot(v3t) >= -1.0 * mEpsilonZero;

    }

    return isInside12 && isInside23 && isInside31;
}


bool BDBoundarySimplexFinder::findInitial2Simplex()
{
    // If the current simplex is 3-simplex, then find the intersecting
    // 2-simplex.
    if (mQ.size()==4) {
        return findInitial2SimplexFrom3Simplex();
    }
    else if  (mQ.size()==2) {
        return findInitial2SimplexFrom1Simplex();
    }
    else if  (mQ.size()==1) {
        return findInitial2SimplexFrom0Simplex();
    }
    return true;
}


bool BDBoundarySimplexFinder::findInitial2SimplexFrom3Simplex()
{
    vector<BDVertex> Qnew;

    const auto& p1 = mQ[0].p();
    const auto& p2 = mQ[1].p();
    const auto& p3 = mQ[2].p();
    const auto& p4 = mQ[3].p();

    const Vec3 v12 = p2 - p1;
    const Vec3 v13 = p3 - p1;
    const Vec3 v14 = p4 - p1;
    const Vec3 v23 = p3 - p2;
    const Vec3 v24 = p4 - p2;

    const bool t132visible((v13.x()*v12.y() - v13.y()*v12.x()) >= 0.0);
    const bool t124visible((v12.x()*v14.y() - v12.y()*v14.x()) >= 0.0);
    const bool t143visible((v14.x()*v13.y() - v14.y()*v13.x()) >= 0.0);
    const bool t234visible((v23.x()*v24.y() - v23.y()*v24.x()) >= 0.0);

    if (t132visible) {
        if(isOriginInsideTriangleXY(p1, p3, p2)) {
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[2]);
            Qnew.push_back(mQ[1]);
            mQ = std::move(Qnew);           
            return true;
        }       
    }
    if (t124visible) {
        if(isOriginInsideTriangleXY(p1, p2, p4)) {
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[1]);
            Qnew.push_back(mQ[3]);
            mQ = std::move(Qnew);           
            return true;
        }       
    }
    if (t143visible) {
        if(isOriginInsideTriangleXY(p1, p4, p3)) {
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[3]);
            Qnew.push_back(mQ[2]);
            mQ = std::move(Qnew);           
            return true;
        }       
    }
    if (t234visible) {
        if(isOriginInsideTriangleXY(p2, p3, p4)) {
            Qnew.push_back(mQ[1]);
            Qnew.push_back(mQ[2]);
            Qnew.push_back(mQ[3]);
            mQ = std::move(Qnew);           
            return true;
        }       
    }


    log(WARNING, __FILE__, __LINE__, 
        "No visible triangles. Relaxing the condition.");
    logQ(WARNING);
    bool isInside132relaxed = isOriginInsideTriangleXYrelaxed(p1, p3, p2);
    bool isInside124relaxed = isOriginInsideTriangleXYrelaxed(p1, p2, p4);
    bool isInside143relaxed = isOriginInsideTriangleXYrelaxed(p1, p4, p3);
    bool isInside234relaxed = isOriginInsideTriangleXYrelaxed(p2, p3, p4);
    Vec3 n132 = v13.cross(v12); n132.normalize();
    Vec3 n124 = v12.cross(v14); n124.normalize();
    Vec3 n143 = v14.cross(v13); n143.normalize();
    Vec3 n234 = v23.cross(v24); n234.normalize();
    double ori132 = isInside132relaxed?n132.z():-1.0;
    double ori124 = isInside124relaxed?n124.z():-1.0;
    double ori143 = isInside143relaxed?n143.z():-1.0;
    double ori234 = isInside234relaxed?n234.z():-1.0;
    if (ori132 >= ori124 && ori132 >= ori143 && ori132 >= ori234) {
        Qnew.push_back(mQ[0]);
        Qnew.push_back(mQ[2]);
        Qnew.push_back(mQ[1]);
        mQ = std::move(Qnew);           
        return true;
    }
    if (ori124 >= ori132 && ori124 >= ori143 && ori124 >= ori234) {
        Qnew.push_back(mQ[0]);
        Qnew.push_back(mQ[1]);
        Qnew.push_back(mQ[3]);
        mQ = std::move(Qnew);           
        return true;
    }
    if (ori143 >= ori132 && ori143 >= ori124 && ori143 >= ori234) {
        Qnew.push_back(mQ[0]);
        Qnew.push_back(mQ[3]);
        Qnew.push_back(mQ[2]);
        mQ = std::move(Qnew);           
        return true;
    }
    if (ori234 >= ori132 && ori234 >= ori124 && ori234 >= ori143) {
        Qnew.push_back(mQ[1]);
        Qnew.push_back(mQ[2]);
        Qnew.push_back(mQ[3]);
        mQ = std::move(Qnew);           
        return true;
    }

    // Error
    log(ERROR, __FILE__, __LINE__, 
        "No visible triangles found.");
    logQ(ERROR);

    return false;
}


bool BDBoundarySimplexFinder::findInitial2SimplexFrom1Simplex()
{
    // Current 1-simplex intersects the origin.
    // Pick an arbitrary vertex and make a triangle.
    const auto& p1  = mQ[0].p();
    const auto& p2  = mQ[1].p();
    auto&       ch1 = mBody1.ConvexHull();
    auto&       ch2 = mBody2.ConvexHull();

    const Vec3 v12 = p2 - p1;

    VertexIt q3v1, q3v2;
    double   maxSq = 0.0;
    for (auto v1  = ch1.vertices().first;
              v1 != ch1.vertices().second;
              v1++                          ) {
        
        for (auto v2  = ch2.vertices().first;
                  v2 != ch2.vertices().second;
                  v2++                        ) {
            if ( v1 != mQ[0].v1() && v2 != mQ[0].v2() && 
                 v1 != mQ[1].v1() && v2 != mQ[1].v2()    ) {

                const Vec3 p3cand = mUseGeomConfigTemp?
                                      ((*v1)->pGCS(mRotMat1, mCoM1) - 
                                       (*v2)->pGCS(mRotMat2, mCoM2)   ):
                                      ((*v1)->pGCS(mScaling, mRotMat1, mCoM1)- 
                                       (*v2)->pGCS(mScaling, mRotMat2, mCoM2));
        
                const Vec3   cr = v12.cross(p3cand - p1);
                const double sq = cr.squaredNorm2();
                if (maxSq==0.0) {
                    maxSq = sq;
                    q3v1 = v1;
                    q3v2 = v2;
                }
                else if (sq > maxSq) {
                    maxSq = sq;
                    q3v1 = v1;
                    q3v2 = v2;
                }
            }
        }        
    }

    BDVertex q3(q3v1, q3v2);
    if (mUseGeomConfigTemp) {
        q3.updateP( mBody1.ConvexHull(), mRotMat1, mCoM1,
                    mBody2.ConvexHull(), mRotMat2, mCoM2 );
    }
    else {
        q3.updateP( mBody1.ConvexHull(), mScaling, mRotMat1, mCoM1,
                    mBody2.ConvexHull(), mScaling, mRotMat2, mCoM2 );
    }
    mQ.push_back(q3);

    return true;
}


bool BDBoundarySimplexFinder::findInitial2SimplexFrom0Simplex()
{
    // Current 0-simplex is the origin.
    // Pick an arbitrary vertex and make an edge first, and
    // then triangle.
    const auto& p1  = mQ[0].p();
    auto&       ch1 = mBody1.ConvexHull();
    auto&       ch2 = mBody2.ConvexHull();

    VertexIt q2v1, q2v2;
    double   maxSq = 0.0;
    for (auto v1  = ch1.vertices().first;
              v1 != ch1.vertices().second;
              v1++                          ) {
        
        for (auto v2  = ch2.vertices().first;
                  v2 != ch2.vertices().second;
                  v2++                        ) {
            if ( v1 != mQ[0].v1() && v2 != mQ[0].v2() ) {

                const Vec3 p2cand = mUseGeomConfigTemp?
                                      ((*v1)->pGCS(mRotMat1, mCoM1) - 
                                       (*v2)->pGCS(mRotMat2, mCoM2)   ):
                                      ((*v1)->pGCS(mScaling, mRotMat1, mCoM1)- 
                                       (*v2)->pGCS(mScaling, mRotMat2, mCoM2));
                const Vec3  v12cand = p2cand - p1;
                const double sq = v12cand.squaredNorm2();
                if (sq > maxSq) {
                    maxSq = sq;
                    q2v1 = v1;
                    q2v2 = v2;
                }
            }
        }        
    }

    BDVertex q2(q2v1, q2v2);
    if (mUseGeomConfigTemp) {
        q2.updateP( mBody1.ConvexHull(), mRotMat1, mCoM1,
                    mBody2.ConvexHull(), mRotMat2, mCoM2 );
    }
    else {
        q2.updateP( mBody1.ConvexHull(), mScaling, mRotMat1, mCoM1,
                    mBody2.ConvexHull(), mScaling, mRotMat2, mCoM2 );
    }
    mQ.push_back(q2);
    findInitial2SimplexFrom1Simplex();

    return true;
}


bool BDBoundarySimplexFinder::findExtremePointsAlongSupportDirection(
    const Vec3& dir,
    VertexIt&   vit1,
    VertexIt&   vit2
) {
    // Find the extremal vertex in the current Q.
    double prj;

    for (long i = 0; i < mQ.size(); i++) {
        if (i == 0) {
            vit1 = mQ[i].v1();
            prj  = dir.dot((*vit1)->pGCS(mRotMat1));
        }
        else {
            VertexIt cit = mQ[i].v1();
            double curPrj  = dir.dot((*cit)->pGCS(mRotMat1));
            if (prj < curPrj) {
                vit1 = cit;
                prj  = curPrj;
            }
        }
    }

    // Find the extremal vertex in manifold 1.
    bool updated = true;
    while (updated) {
        updated = false;
        VertexIt newVit1;
        for (auto heit : (*vit1)->halfEdges()) {
            if ((*heit)->src()==vit1) {
                VertexIt adjit  = (*heit)->dst();
                double   curPrj = dir.dot((*adjit)->pGCS(mRotMat1));
                if (prj < curPrj) {
                    newVit1   = adjit;
                    prj       = curPrj;
                    updated   = true;
                }
            }
        }
        if (updated) {
            vit1 = newVit1;
        }
    }


    for (long i = 0; i < mQ.size(); i++) {
        if (i == 0) {
            vit2 = mQ[i].v2();
            prj  = dir.dot((*vit2)->pGCS(mRotMat2));
        }
        else {
            VertexIt cit = mQ[i].v2();
            double curPrj  = dir.dot((*cit)->pGCS(mRotMat2));
            if (prj > curPrj) {
                vit2 = cit;
                prj  = curPrj;
            }
        }
    }

    // Find the extremal vertex in manifold 2.
    updated = true;
    while (updated) {
         updated = false;
        VertexIt newVit2;
        for (auto heit : (*vit2)->halfEdges()) {
            if ((*heit)->src()==vit2) {
                VertexIt adjit  = (*heit)->dst();
                double   curPrj = dir.dot((*adjit)->pGCS(mRotMat2));
                if (prj > curPrj) {
                    newVit2   = adjit;
                    prj       = curPrj;
                    updated   = true;
                }
            }
        }
        if (updated) {
            vit2 = newVit2;
        }
    }

    for (auto& q : mQ) {
        if ((q.v1()==vit1)&&(q.v2()==vit2)){
            return false;
        }
    }

    const Vec3 v1t = mUseGeomConfigTemp?
                         (( (*vit1)->pGCS(mRotMat1, mCoM1) - 
                            (*vit2)->pGCS(mRotMat2, mCoM2)   ) -
                          mQ[0].p()):
                         (( (*vit1)->pGCS(mScaling, mRotMat1, mCoM1) - 
                            (*vit2)->pGCS(mScaling, mRotMat2, mCoM2)   ) -
                          mQ[0].p());
 
    if (dir.dot(v1t) > 0.0) {
        return true;
    }
    else {
        return false;
    }
}


void BDBoundarySimplexFinder::selectNewTriangle(
    VertexIt vit1, 
    VertexIt vit2 
) {
    BDVertex vNew(vit1, vit2);
    auto& ch1 = mBody1.ConvexHull();
    auto& ch2 = mBody2.ConvexHull();
    if (mUseGeomConfigTemp) {
        vNew.updateP( ch1, mRotMat1, mCoM1,
                      ch2, mRotMat2, mCoM2 );
    }
    else {
        vNew.updateP( ch1, mScaling, mRotMat1, mCoM1,
                      ch2, mScaling, mRotMat2, mCoM2 );
    }

    // Test 3 new triangles.    

    const auto& p1   = mQ[0].p();
    const auto& p2   = mQ[1].p();
    const auto& p3   = mQ[2].p();
    const auto& pNew = vNew.p();

    vector<BDVertex> Qnew;
    bool inside12 = isOriginInsideTriangleXY(p1, p2, pNew);
    bool inside23 = isOriginInsideTriangleXY(p2, p3, pNew);
    bool inside31 = isOriginInsideTriangleXY(p3, p1, pNew);

    if (!inside12 && !inside23 && !inside31) {

        log(WARNING, __FILE__, __LINE__, "No intersecting triangle");
        logQ(WARNING);
        logPoint(WARNING, pNew);

        inside12 = isOriginInsideTriangleXYrelaxed(p1, p2, pNew);
        inside23 = isOriginInsideTriangleXYrelaxed(p2, p3, pNew);
        inside31 = isOriginInsideTriangleXYrelaxed(p3, p1, pNew);

        if (!inside12 && !inside23 && !inside31) {

            log(WARNING, __FILE__, __LINE__, 
                         "No intersecting triangle after relaxation");
            logQ(WARNING);
            logPoint(WARNING, pNew);

            // Picking two points that have greater z.
            if (p1.z() < p2.z() && p1.z() < p3.z()) {
                inside23 = true;
            }
            else if (p2.z() < p1.z() && p2.z() < p3.z()) {
                inside31 = true;
            }
            else {
                inside12 = true;
            }
        }
    }

    double avg12, avg23, avg31;
    // Arbitrary low number not to be chosen.
    double low = (p1.z() + p2.z() + p3.z()) * 2.0 / 3.0 - 1000.0;

    avg12 = inside12?(p1.z() + p2.z()):low;
    avg23 = inside23?(p2.z() + p3.z()):low;
    avg31 = inside31?(p3.z() + p1.z()):low;

    if (avg12 > avg23 && avg12 > avg31) {
        Qnew.push_back(mQ[0]);
        Qnew.push_back(mQ[1]);
        Qnew.push_back(vNew);
        mQ = std::move(Qnew);
        return;
    }
    else if (avg23 > avg12 && avg23 > avg31) {
        Qnew.push_back(mQ[1]);
        Qnew.push_back(mQ[2]);
        Qnew.push_back(vNew);
        mQ = std::move(Qnew);
        return;
    }
    else{
        Qnew.push_back(mQ[2]);
        Qnew.push_back(mQ[0]);
        Qnew.push_back(vNew);
        mQ = std::move(Qnew);
        return;
    }
}


double BDBoundarySimplexFinder::squaredDistanceToOriginXY(
    const Vec3& p1_3d, 
    const Vec3& p2_3d
) {
    const Vec2 p1(p1_3d.x(), p1_3d.y());
    const Vec2 p2(p2_3d.x(), p2_3d.y());
    const Vec2 v12 = p2 - p1;
    if (v12.squaredNorm2() <= mEpsilonZero) {
        return p1.squaredNorm2();
    }
    double s = -1.0 * p1.dot(v12) / v12.squaredNorm2();

    s = std::min(std::max(s, 0.0), 1.0);
    const Vec2 r = p1 - v12 * s;
    return r.squaredNorm2();
}


#ifdef UNIT_TESTS


void BDBoundarySimplexFinder::findSequencerInit() {
    mST_stepNo    = 0;
    mST_stepType  = ST_INIT;
}

void BDBoundarySimplexFinder::findSetParam(const Vec3& Zdir){ mZdir = Zdir; }

void BDBoundarySimplexFinder::findSequencerOneStep() {

    switch (mST_stepType) {

      case ST_INIT:
        
        rotateWorld();
        mST_stepType = ST_AFTER_rotateWorld;
        break;

      case ST_AFTER_rotateWorld:

        mST_found = findInitial2Simplex();
        if (!mST_found) {
            mST_stepType = ST_RETURN_TYPE1;
        }
        else {
            mST_stepType = ST_BEFORE_loop;
        }
        break;

      case ST_BEFORE_loop:

        mST_numIter = 0;
        if (mST_numIter < mMaxNumPivots) {
            mST_stepType = ST_BEGIN_loop;
        }
        else {
            mST_stepType = ST_AFTER_loop;
        }
        break;

      case ST_BEGIN_loop:

        mST_dir      = findSupportDirectionOfTriangle();
        mST_stepType = ST_AFTER_findSupportDirectionOfTriangle;
        break;

      case ST_AFTER_findSupportDirectionOfTriangle:

        mST_found = findExtremePointsAlongSupportDirection(
                                                  mST_dir, mST_vit1, mST_vit2);

        if (!mST_found) {
            mST_stepType = ST_RETURN_TYPE2;
        }
        else {
            mST_stepType = ST_AFTER_findExtremePointsAlongSupportDirection;
        }
        break;

      case ST_AFTER_findExtremePointsAlongSupportDirection:

        selectNewTriangle(mST_vit1, mST_vit2);
        mST_stepType =  ST_AFTER_selectNewTriangle;
        break;

      case ST_AFTER_selectNewTriangle:
        mST_numIter++;
        if (mST_numIter < mMaxNumPivots) {
            mST_stepType = ST_BEGIN_loop;
        }
        else {
            mST_stepType = ST_AFTER_loop;
        }
        break;

      case ST_AFTER_loop:
        mST_stepType = ST_RETURN_TYPE3;
        break;

      default:
        break;
    }
}


void BDBoundarySimplexFinder::renderObj1Surface(
    Vec3&         color,
    float         alpha,
    vector<Vec3>& vertices,
    vector<Vec3>& colors,
    vector<float>&alphas,
    vector<Vec3>& normals
) {
    auto& mani = mBody1.ConvexHull();
    mani.makeOpenGLVerticesColorsNormalsForTriangles(
        color,
        alpha,
        vertices,
        colors,
        alphas,
        normals,
        false,
        mScaling,
        mUseGeomConfigTemp?mBody1.QmatTemp():mBody1.Qmat(),
        mUseGeomConfigTemp?mBody1.CoMTemp():mBody1.CoM() );
}


void BDBoundarySimplexFinder::renderObj2Surface(
    Vec3&         color,
    float         alpha,
    vector<Vec3>& vertices,
    vector<Vec3>& colors,
    vector<float>&alphas,
    vector<Vec3>& normals
) {
    auto& mani = mBody2.ConvexHull();
    mani.makeOpenGLVerticesColorsNormalsForTriangles(
        color,
        alpha,
        vertices,
        colors,
        alphas,
        normals,
        false,
        mScaling,
        mUseGeomConfigTemp?mBody2.QmatTemp():mBody2.Qmat(),
        mUseGeomConfigTemp?mBody2.CoMTemp():mBody2.CoM() );
}


void BDBoundarySimplexFinder::renderObj1Lines(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
) {
    auto& mani = mBody1.ConvexHull();
    vector<HalfEdgeIt> halfEdges;
    for (auto eit = mani.edges().first; eit != mani.edges().second; eit++) {
        halfEdges.push_back((*eit)->he1());
    }
    mani.makeOpenGLVerticesColorsForLines(
        halfEdges,
        color,
        vertices,
        colors,
        mScaling,
        mUseGeomConfigTemp?mBody1.QmatTemp():mBody1.Qmat(),
        mUseGeomConfigTemp?mBody1.CoMTemp():mBody1.CoM() );
}


void BDBoundarySimplexFinder::renderObj2Lines(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
) {
    auto& mani = mBody2.ConvexHull();
    vector<HalfEdgeIt> halfEdges;
    for (auto eit = mani.edges().first; eit != mani.edges().second; eit++) {
        halfEdges.push_back((*eit)->he1());
    }
    mani.makeOpenGLVerticesColorsForLines(
        halfEdges,
        color,
        vertices,
        colors,
        mScaling,
        mUseGeomConfigTemp?mBody2.QmatTemp():mBody2.Qmat(),
        mUseGeomConfigTemp?mBody2.CoMTemp():mBody2.CoM() );
}


void BDBoundarySimplexFinder::renderBDManifoldSurface(
    Vec3&         color,
    float         alpha,
    vector<Vec3>& vertices,
    vector<Vec3>& colors,
    vector<float>&alphas,
    vector<Vec3>& normals
) {
    auto& mani1 = mBody1.ConvexHull();
    auto& mani2 = mBody2.ConvexHull();
    vector<Vec3> diffPoints;
    auto Vpair1 = mani1.vertices();
    auto Vpair2 = mani2.vertices();
    for (auto vit1 = Vpair1.first; vit1 != Vpair1.second; vit1++) {

        auto p1 = mUseGeomConfigTemp?
                  ((*vit1)->pGCS(mBody1.QmatTemp(), mBody1.CoMTemp())):
                  ((*vit1)->pGCS(mScaling, mBody1.Qmat(), mBody1.CoM()));
        for (auto vit2 = Vpair2.first; vit2 != Vpair2.second; vit2++) {
            auto p2 = mUseGeomConfigTemp?
                  ((*vit2)->pGCS(mBody2.QmatTemp(), mBody2.CoMTemp())):
                  ((*vit2)->pGCS(mScaling, mBody2.Qmat(), mBody2.CoM()));
            diffPoints.push_back(p1 - p2);
        }
    }

    Manifold       maniDiff;
    enum predicate pred;
    maniDiff.findConvexHull(diffPoints, pred);

    maniDiff.makeOpenGLVerticesColorsNormalsForTriangles(
        color,
        alpha,
        vertices,
        colors,
        alphas,
        normals,
        false );

}


void BDBoundarySimplexFinder::renderBDManifoldLines(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
) {
    auto& mani1 = mBody1.ConvexHull();
    auto& mani2 = mBody2.ConvexHull();

    vector<Vec3> diffPoints;
    auto Vpair1 = mani1.vertices();
    auto Vpair2 = mani2.vertices();
    for (auto vit1 = Vpair1.first; vit1 != Vpair1.second; vit1++) {
        auto p1 = mUseGeomConfigTemp?
                  ((*vit1)->pGCS(mBody1.QmatTemp(), mBody1.CoMTemp())):
                  ((*vit1)->pGCS(mScaling, mBody1.Qmat(), mBody1.CoM()));
        for (auto vit2 = Vpair2.first; vit2 != Vpair2.second; vit2++) {
            auto p2 = mUseGeomConfigTemp?
                  ((*vit2)->pGCS(mBody2.QmatTemp(), mBody2.CoMTemp())):
                  ((*vit2)->pGCS(mScaling, mBody2.Qmat(), mBody2.CoM()));
            diffPoints.push_back(p1 - p2);
            
            const Vec3 p2p1 = p1 - p2;
        }
    }


    Manifold       maniDiff;
    enum predicate pred;
    maniDiff.findConvexHull(diffPoints, pred);

    vector<HalfEdgeIt> halfEdges;
    for (auto eit  = maniDiff.edges().first; 
              eit != maniDiff.edges().second; eit++) {
        halfEdges.push_back((*eit)->he1());
    }

    maniDiff.makeOpenGLVerticesColorsForLines(
        halfEdges,
        color,
        vertices,
        colors
    );

}


void BDBoundarySimplexFinder::renderActiveFeature1Points(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
) {
    vector<Vec3> points;
    for (auto& q : mQ) {
        auto vit = q.v1();
        auto p1 = mUseGeomConfigTemp?
                      ((*vit)->pGCS(mBody1.QmatTemp(), mBody1.CoMTemp())):
                      ((*vit)->pGCS(mScaling, mBody1.Qmat(), mBody1.CoM()));
        bool skip = false;
        for (auto& p2 : points) {
            if (p1==p2) {
                skip = true;
            }
        }
        if (!skip) {
            points.push_back(p1);
        }
    }

    Manifold::makeOpenGLVerticesColorsForPoints(
        points,
        color,
        vertices,
        colors
    );
}


void BDBoundarySimplexFinder::renderActiveFeature2Points(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
) {
    vector<Vec3> points;
    for (auto& q : mQ) {
        auto vit = q.v2();
        auto p1 = mUseGeomConfigTemp?
                      ((*vit)->pGCS(mBody2.QmatTemp(), mBody2.CoMTemp())):
                      ((*vit)->pGCS(mScaling, mBody2.Qmat(), mBody2.CoM()));
        bool skip = false;
        for (auto& p2 : points) {
            if (p1==p2) {
                skip = true;
            }
        }
        if (!skip) {
            points.push_back(p1);
        }
    }

    Manifold::makeOpenGLVerticesColorsForPoints(
        points,
        color,
        vertices,
        colors
    );
}


void BDBoundarySimplexFinder::renderActiveFeatureBDManifoldPoints(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
){
    vector<Vec3> points;
    for (auto& q : mQ) {
        auto vit1 = q.v1();
        auto vit2 = q.v2();
        auto p1 = mUseGeomConfigTemp?
                      ((*vit1)->pGCS(mBody1.QmatTemp(), mBody1.CoMTemp())):
                      ((*vit1)->pGCS(mScaling, mBody1.Qmat(), mBody1.CoM()));

        auto p2 = mUseGeomConfigTemp?
                      ((*vit2)->pGCS(mBody2.QmatTemp(), mBody2.CoMTemp())):
                      ((*vit2)->pGCS(mScaling, mBody2.Qmat(), mBody2.CoM()));
        points.push_back(p1-p2);
    }

    Manifold::makeOpenGLVerticesColorsForPoints(
        points,
        color,
        vertices,
        colors
    );
}


void BDBoundarySimplexFinder::renderDepartingDirLine(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
) {
    vector<Vec3> points;
    Vec3 n = mZdir;
    n.normalize();
    points.push_back(Vec3(0.0, 0.0, 0.0));
    points.push_back(n);
    Manifold::makeOpenGLVerticesColorsForPoints(
                                              points, color, vertices, colors);
}


#endif

}// namespace Makena

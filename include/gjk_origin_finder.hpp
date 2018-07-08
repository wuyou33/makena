#ifndef _MAKENA_GJK_ORIGIN_FINDER_HPP_
#define _MAKENA_GJK_ORIGIN_FINDER_HPP_

#include <memory>
#include <array>
#include <iostream>
#include <string>
#include <list>
#include <map>
#include <vector>
#include <exception>
#include <stdexcept>
#include <cmath>
#include <cstdarg>
#include <iomanip> 

#include "convex_rigid_body.hpp"
#include "binary_dilation.hpp"
#include "voronoi_3simplex_graph.hpp"
#include "loggable.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file gjk_origin_finder.hpp
 *
 * @brief GJK collision detector on a pair of convex polytopes.

 */
namespace Makena {

using namespace std;


class GJKOriginFinder : public Loggable {

  public:

    /** @brief constructor. This is a lightweight process with just two 
     *         references.
     *
     *  @param body1        (in): Rigid body 1
     *
     *  @param body2        (in): Rigid body 2
     *
     *  @param Q  (inout): set of active vertices in the convex binary dilation
     *                     manifold that represents the current simplex.
     *                     Upon exit, it contains the 3-simlex that contains
     *                     the origin iside, or {0,1,2)-simplex that contains
     *                     the closest point to the origin in the convex 
     *                     binary dilation manifold.
     *
     *  @param maxNumIter   (in): Max number of iterations allowed before
     *                            giving up. This is a safety mechanism.
     *
     *  @param maxNumCycles (in): Max number of occurences where the distance
     *                            to the closest point goes up.
     *                            This is to absorb numerical fluctuation.
     *
     *  @param epsilonZero  (in): numerical tolerance to consider the
     *                            distance value zero.
     */
    inline GJKOriginFinder(
        ConvexRigidBody&         body1,
        ConvexRigidBody&         body2,
        vector<BDVertex>&        Q,
        const long&              maxNumIter         = 100,
        const long&              maxNumCycles       = 8,
        const double&            epsilonZero        = EPSILON_SQUARED,
        bool                     useGeomConfigTemp = true,
        const double&            scaling            = 1.0,
        std::ostream&            logStream          = std::cerr
    );

    /** @brief destructor. Nothing to be done.
     */
    inline ~GJKOriginFinder();

    /** @brief main function to try to explore toward orign in the 
     *         convex binary dilation manifold.
     *
     *         If the final simplex enclosing the origin is 3-simplex, then
     *         the order of the vertices in Q are such that 
     *         the normal of triangle (Q[0], Q[1], Q[2]) is facing Q[3],
     *         i.e., the 3-simplex morphs to the following tetrahedron.
     *
     *                ((0,0,0), (1,0.0), (0,1,0), (0,0,1)).
     *                 Origin    X-unit   Y-unit   Z-unit
     *
     *  @param coldStart (in): set to true if this is a cold start from 
     *                     without starting simplex set in Q.
     *
     *  @return  true  if the origin is contained in the convex binary dilation
     *                 manifold.
     *           false otherwise
     */
    bool findOrigin(const bool coldStart);

#ifndef UNIT_TESTS
  private:
#else
  public:
#endif

    inline void findInitialSimplex();


    void findClosestPointToOrigin(
        Vec3&           closestP,
        enum predicate& pred
    );


    void findClosestPointToOrigin3Simplex(
        Vec3&           closestP,
        enum predicate& pred
    );

    void findClosestPointToOrigin2Simplex(
        Vec3&           closestP,
        enum predicate& pred
    );


    void findClosestPointToOrigin1Simplex(
        Vec3&           closestP,
        enum predicate& pred
    );


    void findClosestPointToOrigin0Simplex(
        Vec3&           closestP,
        enum predicate& pred
    );


    inline void projectPointOnPlane(
        const Vec3& p1,
        const Vec3& p2,
        const Vec3& p3,
        const Vec3& pTest,
        Vec3&       proj
    );


    /** @brief classify a point against an edge
     *
     *  @param p1    (in): point 1 of the edge
     *
     *  @param p2    (in): point 2 of the edge
     *
     *  @param pTest (in): point to be classified
     *
     *  @param proj  (out): closest point on the edge to pTest
     *
     *  @return  VORONOI_INSIDE_TRIANGLE- pTest is in the triangle inclusive
     *           VORONOI_OVER_EDGE_1_2  - pTest is over (1,2) inclusive
     *           VORONOI_OVER_EDGE_2_3  - pTest is over (2,3) inclusive
     *           VORONOI_OVER_EDGE_3_1  - pTest is over (3,1) inclusive
     *           VORONOI_OVER_VERTEX_1  - pTest is over p1
     *           VORONOI_OVER_VERTEX_2  - pTest is over p2
     *           VORONOI_OVER_VERTEX_3  - pTest is over p3
     *           NONE                   - Error. Ill-conditioned triangle.
     */
    inline enum predicate testPointAgainstEdge(
        const Vec3& p1,
        const Vec3& p2,
        const Vec3& pTest,
        Vec3&       proj
    );// UNIT TESTED


    inline Vec3 projectPointOnLine(
        const Vec3& p1,
        const Vec3& p2,
        const Vec3& pTest
    );// UNIT TESTED


    /** @brief classify a point against a triangle
     *
     *  @param p1    (in): point 1 of triangle
     *
     *  @param p2    (in): point 2 of triangle
     *
     *  @param p3    (in): point 3 of triangle
     *
     *  @param pTest (in): point to be classified
     *
     *  @return  VORONOI_INSIDE_TRIANGLE- pTest is in the triangle inclusive
     *           VORONOI_OVER_EDGE_1_2  - pTest is over (1,2) inclusive
     *           VORONOI_OVER_EDGE_2_3  - pTest is over (2,3) inclusive
     *           VORONOI_OVER_EDGE_3_1  - pTest is over (3,1) inclusive
     *           VORONOI_OVER_VERTEX_1  - pTest is over p1
     *           VORONOI_OVER_VERTEX_2  - pTest is over p2
     *           VORONOI_OVER_VERTEX_3  - pTest is over p3
     *           NONE                   - Error. Ill-conditioned triangle.
     */
    inline enum predicate testPointAgainstTriangle(
        const Vec3& p1,
        const Vec3& p2,
        const Vec3& p3,
        const Vec3& pTest
    );// UNIT TESTED

    
    void removeUnsupportingVertices(const enum predicate pred);


    bool findExtremePointsAlongSupportDirection(
        const Vec3& closestP,
        VertexIt&   vit1,
        VertexIt&   vit2
    );


    void updateSimplexTowardOrigin(
        const enum predicate     pred,
        VertexIt&                vit1,
        VertexIt&                vit2
    );


    inline bool isInWrongOrderForTetrahedron();// UNIT_TESTED

    inline void logPoints(
        enum LogLevel  lvl,
        const Vec3&    p1, 
        const Vec3&    p2,
        const Vec3&    p3,
        const Vec3&    p4
    );


    inline void logPoints(
        enum LogLevel  lvl,
        const Vec3&    p1, 
        const Vec3&    p2,
        const Vec3&    p3
    );

    inline void logPoints(
        enum LogLevel  lvl,
        const Vec3&    p1, 
        const Vec3&    p2
    );

    inline void logPoint(
        enum LogLevel  lvl,
        const Vec3&    p1 
    );

    inline void logPredicate(
        enum LogLevel        lvl,
        const enum predicate pred
    );

    inline void logQ(enum LogLevel lvl);


    ConvexRigidBody&         mBody1;
    ConvexRigidBody&         mBody2;
    vector<BDVertex>&        mQ;
    const long               mMaxNumIter;
    const long               mMaxNumCycles;
    const double             mEpsilonZero;
    Voronoi3SimplexGraph     mGraph;
    const bool               mUseGeomConfigTemp;
    const double             mScaling;


#ifdef UNIT_TESTS

    enum stepType {
        ST_INIT,
        ST_BEFORE_loop,
        ST_BEGIN_loop,
        ST_AFTER_findClosestPointToOrigin,
        ST_AFTER_removeUnsupportingVertices,
        ST_AFTER_findExtremePointsAlongSupportDirection,
        ST_AFTER_updateSimplexTowardOrigin,
        ST_AFTER_loop,
        ST_RETURN_TYPE1,
        ST_RETURN_TYPE2,
        ST_RETURN_TYPE3,
        ST_RETURN_TYPE4,
        ST_RETURN_TYPE5,
        ST_TERMINATED
    }              mST_stepType;
    long           mST_stepNo;
    long           mST_coldStart;
    long           mST_numIter;
    long           mST_numCycles;
    enum predicate mST_pred;
    Vec3           mST_closestP;
    VertexIt       mST_vit1, mST_vit2;
    bool           mST_found;
    double         mST_sqDist;
    double         mST_prevSqDist;
    Vec3           mST_prevSqP;
     vector<pair<long,long> > 
                   mST_prevQ;

    void findOriginSequencerInit(const bool coldStart);
    void findOriginSequencerOneStep();


    // Rendering Helpers.

    void renderObj1Surface(
        Vec3&         color,
        float         alpha,
        vector<Vec3>& vertices,
        vector<Vec3>& colors,
        vector<float>&alphas,
        vector<Vec3>& normals
    );

    void renderObj2Surface(
        Vec3&         color,
        float         alpha,
        vector<Vec3>& vertices,
        vector<Vec3>& colors,
        vector<float>&alphas,
        vector<Vec3>& normals
    );


    void renderObj1Lines(
        Vec3&         color,
        vector<Vec3>& vertices,
        vector<Vec3>& colors
    );


    void renderObj2Lines(
        Vec3&         color,
        vector<Vec3>& vertices,
        vector<Vec3>& colors
    );


    void renderBDManifoldSurface(
        Vec3&         color,
        float         alpha,
        vector<Vec3>& vertices,
        vector<Vec3>& colors,
        vector<float>&alphas,
        vector<Vec3>& normals
    );

    void renderBDManifoldLines(
        Vec3&         color,
        vector<Vec3>& vertices,
        vector<Vec3>& colors
    );

    void renderActiveFeature1Points(
        Vec3&         color,
        vector<Vec3>& vertices,
        vector<Vec3>& colors
    );

    void renderActiveFeature2Points(
        Vec3&         color,
        vector<Vec3>& vertices,
        vector<Vec3>& colors
    );

    void renderActiveFeatureBDManifoldPoints(
        Vec3&         color,
        vector<Vec3>& vertices,
        vector<Vec3>& colors
    );

    void renderSupportVectorLine(
        Vec3&               color,
        vector<Vec3>&       vertices,
        vector<Vec3>&       colors
    );

    void renderClosestPoint(
        Vec3&               color,
        vector<Vec3>&       vertices,
        vector<Vec3>&       colors
    );

#endif

};


GJKOriginFinder::GJKOriginFinder(

    ConvexRigidBody&         body1,
    ConvexRigidBody&         body2,
    vector<BDVertex>&        Q,
    const long&              maxNumIter,
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
     mMaxNumIter(maxNumIter),
     mMaxNumCycles(maxNumCycles),
     mEpsilonZero(epsilonZero),
     mUseGeomConfigTemp(useGeomConfigTemp),
     mScaling(scaling)
     {;}
     

GJKOriginFinder::~GJKOriginFinder(){;}


void GJKOriginFinder::logQ(enum LogLevel lvl)
{
    if (mLogLevel>0 && mLogLevel>=lvl) {
        mLogStream << "Q size: " << mQ.size() << "\n";

        for (int i = 0; i < mQ.size() ; i++) {

            auto vit1 = mQ[i].v1();
            auto vit2 = mQ[i].v2();
            auto id1  = (*vit1)->id();
            auto id2  = (*vit2)->id();
            auto p1 = mUseGeomConfigTemp?
                      ((*vit1)->pGCS(mBody1.QmatTemp(), mBody1.CoMTemp())):
                      ((*vit1)->pGCS(mScaling, mBody1.Qmat(), mBody1.CoM()));

            auto p2 = mUseGeomConfigTemp?
                      ((*vit2)->pGCS(mBody2.QmatTemp(), mBody2.CoMTemp())):
                      ((*vit2)->pGCS(mScaling, mBody2.Qmat(), mBody2.CoM()));

            mLogStream << "Enrty [" << i << "] "
                       << "V1[" << id1 << ": " << p1 << " "
                       << "V2[" << id2 << ": " << p2 << " "
                       << "P: " <<  mQ[i].p() << "\n";
        }
        mLogStream << "\n";
    }
}


void GJKOriginFinder::logPoints(
    enum LogLevel lvl,
    const Vec3&    p1, 
    const Vec3&    p2,
    const Vec3&    p3,
    const Vec3&    p4
) {
    if (mLogLevel>0 && mLogLevel>=lvl) {
        mLogStream << p1 << "\n";
        mLogStream << p2 << "\n";
        mLogStream << p3 << "\n";
        mLogStream << p4 << "\n";
    }
}


void GJKOriginFinder::logPoints(
    enum LogLevel lvl,
    const Vec3&   p1, 
    const Vec3&   p2,
    const Vec3&   p3
) {
    if (mLogLevel>0 && mLogLevel>=lvl) {
        mLogStream << p1 << "\n";
        mLogStream << p2 << "\n";
        mLogStream << p3 << "\n";
    }
}


void GJKOriginFinder::logPoints(
    enum LogLevel lvl,
    const Vec3&   p1, 
    const Vec3&   p2
) {
    if (mLogLevel>0 && mLogLevel>=lvl) {
        mLogStream << p1 << "\n";
        mLogStream << p2 << "\n";
    }
}


void GJKOriginFinder::logPoint(
    enum LogLevel lvl,
    const Vec3&   p1
) {
    if (mLogLevel>0 && mLogLevel>=lvl) {
        mLogStream << p1 << "\n";
    }
}


void GJKOriginFinder::logPredicate(
    enum LogLevel        lvl,
    const enum predicate pred
) {

    if (mLogLevel>0 && mLogLevel>=lvl) {

        switch (pred) {

          case VORONOI_INSIDE_TETRAHEDRON:
            mLogStream << "VORONOI_INSIDE_TETRAHEDRON\n";
            break;

          case VORONOI_OVER_TRIANGLE_1_3_2:
            mLogStream << "VORONOI_OVER_TRIANGLE_1_3_2\n";
            break;

          case VORONOI_OVER_TRIANGLE_1_2_4:
            mLogStream << "VORONOI_OVER_TRIANGLE_1_2_4\n";
            break;

          case VORONOI_OVER_TRIANGLE_1_4_3:
            mLogStream << "VORONOI_OVER_TRIANGLE_1_4_3\n";
            break;

          case VORONOI_OVER_TRIANGLE_2_3_4:
            mLogStream << "VORONOI_OVER_TRIANGLE_2_3_4\n";
            break;

          case VORONOI_OVER_EDGE_1_2:
            mLogStream << "VORONOI_OVER_EDGE_1_2\n";
            break;

          case VORONOI_OVER_EDGE_2_3:
            mLogStream << "VORONOI_OVER_EDGE_2_3\n";
            break;

          case VORONOI_OVER_EDGE_3_1:
            mLogStream << "VORONOI_OVER_EDGE_3_1\n";
            break;

          case VORONOI_OVER_EDGE_1_4:
            mLogStream << "VORONOI_OVER_EDGE_1_4\n";
            break;

          case VORONOI_OVER_EDGE_2_4:
            mLogStream << "VORONOI_OVER_EDGE_2_4\n";
            break;

          case VORONOI_OVER_EDGE_3_4:
            mLogStream << "VORONOI_OVER_EDGE_3_4\n";
            break;

          case VORONOI_OVER_VERTEX_1:
            mLogStream << "VORONOI_OVER_VERTEX_1\n";
            break;

          case VORONOI_OVER_VERTEX_2:
            mLogStream << "VORONOI_OVER_VERTEX_2\n";
            break;

          case VORONOI_OVER_VERTEX_3:
            mLogStream << "VORONOI_OVER_VERTEX_3\n";
            break;

          case VORONOI_OVER_VERTEX_4:
            mLogStream << "VORONOI_OVER_VERTEX_4\n";
            break;

          case VORONOI_INSIDE_TRIANGLE:
            mLogStream << "VORONOI_INSIDE_TRIANGLE\n";
            break;

          case BETWEEN_1_AND_2:
            mLogStream << "BETWEEN_1_AND_2\n";
            break;

          case ON_POINT1:
            mLogStream << "ON_POINT1\n";
            break;

          case ON_POINT2:
            mLogStream << "ON_POINT2\n";
            break;

          case NONE:
            mLogStream << "NONE\n";
            break;

          default:
            mLogStream << "<Unknown>: " << pred << "\n";

        }
    }
}


void GJKOriginFinder::findInitialSimplex()
{

    mQ.clear();

    auto& ch1 = mBody1.ConvexHull();
    auto& ch2 = mBody2.ConvexHull();

    auto vit1 = ch1.vertices().first;
    auto vit2 = ch2.vertices().first;

    BDVertex mv(vit1, vit2);
    if (mUseGeomConfigTemp) {
        mv.updateP(ch1, mBody1.QmatTemp(), mBody1.CoMTemp(),
                   ch2, mBody2.QmatTemp(), mBody2.CoMTemp());
    }
    else {
        mv.updateP(ch1, mScaling, mBody1.Qmat(), mBody1.CoM(),
                   ch2, mScaling, mBody2.Qmat(), mBody2.CoM());
    }
    mQ.push_back(std::move(mv));

}


void GJKOriginFinder::projectPointOnPlane(
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& p3,
    const Vec3& pTest,
    Vec3&       proj
) {

    const Vec3 v12        = p2    - p1;
    const Vec3 v13        = p3    - p1;
    const Vec3 v1Test     = pTest - p1;
    const Vec3 v12_cr_v13 = v12.cross(v13);

    const Vec3 v12perp    = v12_cr_v13.cross(v12);
    const Vec3 v13perp    = v12_cr_v13.cross(v13);

    const double v12_dot_v13perp = v12.dot(v13perp);
    const double v13_dot_v12perp = v13.dot(v12perp);

    // Checking for degeneracy
    if ( fabs(v12_dot_v13perp) < EPSILON_SQUARED &&
         fabs(v13_dot_v12perp) < EPSILON_SQUARED    ) {

        // v12 and v13 are parallel to each other.

        log(WARNING, __FILE__, __LINE__, "Degenerate triangle");
        logPoints(ERROR, p1, p2, p3);
        logQ(ERROR);

        if (v12.dot(v13) >= 0.0) {
            if (v12.squaredNorm2() > v13.squaredNorm2()) {
                proj = projectPointOnLine(p1, p2, pTest);
            }
            else {
                proj = projectPointOnLine(p1, p3, pTest);
            }        
        }
        else {
            proj= projectPointOnLine(p2, p3, pTest);
        }
    }
    else {

        const double s = v1Test.dot(v13perp) / v12_dot_v13perp;
        const double t = v1Test.dot(v12perp) / v13_dot_v12perp;

        proj =  p1 + (v12 * s) + (v13 * t);
    }
}


enum predicate GJKOriginFinder::testPointAgainstEdge(
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& pTest,
    Vec3&       proj
) {
    const Vec3  v12     = p2   - p1;
    const Vec3  v1t     = pTest - p1;
    const double sqDist = v12.squaredNorm2();

    if ( fabs(sqDist) < EPSILON_SQUARED ) {

        log(WARNING, __FILE__, __LINE__, "Degenerate edge (p1,p2)");
        logPoints(WARNING, p1, p2);
        logQ(ERROR);

        const double sqDist1 = v1t.squaredNorm2();

        if (sqDist1 < EPSILON_SQUARED) {
            proj = p1;
            return BETWEEN_1_AND_2;
        }
        else {
            const Vec3 v2t = pTest - p2;
            const double sqDist2 = v2t.squaredNorm2();
            if (sqDist1 < sqDist2) {
                proj = p1;
                return ON_POINT1;
            }
            else {
                proj = p2;
                return ON_POINT2;
            }
        }
    }

    const double t = v1t.dot(v12) / sqDist;

    if (t > 1.0 ) {
        proj = p2;
        return ON_POINT2;
    }
    else if (t < 0.0) {
        proj = p1;
        return ON_POINT1;
    }
    else {

        proj = p1 + (v12 * t);
        return BETWEEN_1_AND_2;
    }
}


Vec3 GJKOriginFinder::projectPointOnLine(
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& pTest
) {
    const Vec3  v12     = p2   - p1;
    const Vec3  v1t     = pTest - p1;
    const double sqDist = v12.squaredNorm2();

    if ( fabs(sqDist) < EPSILON_SQUARED ) {

        log(WARNING, __FILE__, __LINE__, "Degenerate edge (p1,p2)");
        logPoints(WARNING, p1, p2);
        logQ(ERROR);

        const double sqDist1 = v1t.squaredNorm2();

        if (sqDist1 < EPSILON_SQUARED) {
            return p1;
        }
        else {
            const Vec3 v2t = pTest - p2;
            const double sqDist2 = v2t.squaredNorm2();
            if (sqDist1 < sqDist2) {
                return p1;
            }
            else {
                return p2;
            }
        }
    }

    const double t = v1t.dot(v12) / sqDist;

    return p1 + (v12 * t);
}


enum predicate GJKOriginFinder::testPointAgainstTriangle(
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& p3,
    const Vec3& pTest
) {

    const Vec3 v12 = p2 - p1;
    const Vec3 v13 = p3 - p1;
    const Vec3 v23 = p3 - p2;
    const Vec3 v31 = p1 - p3;
    const Vec3 v1t = pTest - p1;
    const Vec3 v2t = pTest - p2;
    const Vec3 v3t = pTest - p3;

    const Vec3 n123       = v12.cross(v13);
    const Vec3 v12_cr_v1t = v12.cross(v1t);
    const Vec3 v23_cr_v2t = v23.cross(v2t);
    const Vec3 v31_cr_v3t = v31.cross(v3t);

    bool isInside12 = v12_cr_v1t.dot(n123) >= 0.0;
    bool isInside23 = v23_cr_v2t.dot(n123) >= 0.0;
    bool isInside31 = v31_cr_v3t.dot(n123) >= 0.0;

    if (isInside12 && isInside23 && isInside31) {

        return VORONOI_INSIDE_TRIANGLE;
    }

    // Cycle in CCW: 1, (1,2), 2, (2,3), 3, (3,1)
    //
    //                  1
    //                 / \
    //                /   \
    //               2-----3
    //
    enum _testType {
        UNKNOWN,
        CW,
        CCW
    };

    enum _testType test_1_12 = UNKNOWN;
    enum _testType test_12_2 = UNKNOWN;
    enum _testType test_2_23 = UNKNOWN;
    enum _testType test_23_3 = UNKNOWN;
    enum _testType test_3_31 = UNKNOWN;
    enum _testType test_31_1 = UNKNOWN;

    if (isInside12) {
        test_1_12 = CW;
        test_12_2 = CCW;
    }
    else {
        test_1_12 = (v12.dot(v1t)>0.0)?CCW:CW;
        test_12_2 = (v12.dot(v2t)>0.0)?CCW:CW;
    }

    if (isInside23) {
        test_2_23 = CW;
        test_23_3 = CCW;
    }
    else {
        test_2_23 = (v23.dot(v2t)>0.0)?CCW:CW;
        test_23_3 = (v23.dot(v3t)>0.0)?CCW:CW;
    }

    if (isInside31) {
        test_3_31 = CW;
        test_31_1 = CCW;
    }
    else {
        test_3_31 = (v31.dot(v3t)>0.0)?CCW:CW;
        test_31_1 = (v31.dot(v1t)>0.0)?CCW:CW;
    }

    if (test_31_1 == CCW && test_1_12 == CW && !(isInside31 && isInside12)) {
        return VORONOI_OVER_VERTEX_1;
    }
    else if (test_1_12 == CCW && test_12_2 == CW && !isInside12) {
        return VORONOI_OVER_EDGE_1_2;
    }
    else if (test_12_2 == CCW && test_2_23 == CW && 
             !(isInside12 && isInside23)            ) {
        return VORONOI_OVER_VERTEX_2;
    }
    else if (test_2_23 == CCW && test_23_3 == CW && !isInside23) {
        return VORONOI_OVER_EDGE_2_3;
    }
    else if (test_23_3 == CCW && test_3_31 == CW && 
             !(isInside23 && isInside31)            ) {
        return VORONOI_OVER_VERTEX_3;
    }
    else if (test_3_31 == CCW && test_31_1 == CW && !isInside31) {
        return VORONOI_OVER_EDGE_3_1;
    }

    log(ERROR, __FILE__, __LINE__, "Test inconsistent.");
    logPoints(ERROR, p1, p2, p3, pTest);
    logQ(ERROR);

    return NONE;
}


bool GJKOriginFinder::isInWrongOrderForTetrahedron()
{
    const Vec3& p1  = mQ[0].p();
    const Vec3& p2  = mQ[1].p();
    const Vec3& p3  = mQ[2].p();
    const Vec3& p4  = mQ[3].p();
    const Vec3  v12 = p2 - p1;
    const Vec3  v13 = p3 - p1;
    const Vec3  v14 = p4 - p1;
    const Vec3  v12_cr_v13 = v12.cross(v13);

    return v12_cr_v13.dot(v14) < 0.0;
}

}// namespace Makena


#endif /*_MAKENA_GJK_ORIGIN_FINDER_HPP_*/

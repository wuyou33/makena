#ifndef _MAKENA_BD_BOUNDARY_SIMPLEX_FINDER_HPP_
#define _MAKENA_BD_BOUNDARY_SIMPLEX_FINDER_HPP_

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
#include "loggable.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file bd_boundary_simplex_finder.hpp
 *
 * @brief after finiding a {0,1,2,3}-simplex of the convex binary dilation
 *        manifold that encloses the origin by GJKOriginFinder, this performs
 *        pivoting along the departing direction
 *        until it finds the boundary minimal {0,1,2}-simplex that intersects
 *        the ray from the origin along the departing direction.
 */
namespace Makena {

using namespace std;


class BDBoundarySimplexFinder : public Loggable {

  public:

    /** @brief constructor. This is a lightweight process with just two 
     *         references.
     *
     *  @param body1        (in): Rigid body 1
     *
     *  @param body2        (in): Rigid body 2
     *
     *  @param Q  (in/out):set of active vertices in the convex binary dilation
     *                     manifold that represents the minimal boundary
     *                     simplex that intersects the ray from the origin
     *                     along the departing direction.
     *
     *  @param maxNumPivots (in): Max number of pivoting performed before
     *                            giving up. This is a safety mechanism.
     *
     *  @param epsilonZero  (in): numerical tolerance to consider the
     *                            distance value zero.
     *  @param logStream    (in): the output stream for logging.
     */
    BDBoundarySimplexFinder(
        ConvexRigidBody&         body1,
        ConvexRigidBody&         body2,
        vector<BDVertex>&        Q,
        const long&              maxNumPivots      = 100,
        const long&              maxNumCycles      = 8,
        const double&            epsilonZero       = EPSILON_SQUARED,
        bool                     useGeomConfigTemp = true,
        const double&            scaling           = 1.0,
        std::ostream&            logStream         = std::cerr
    );


    /** @brief destructor. Nothing to be done.
     */
    ~BDBoundarySimplexFinder();


    /** @brief main function that finds the minimal boundary simlex of 
     *         the convex binary dilation manifold that intersects the ray
     *         from the origin along the departing direction
     *
     *  @param Zdir         (in): the departing direction of the binary 
     *                            from the origin. Usually, it is the negative
     *                            of the departing direction of body 1 relative
     *                            to body 2. (i.e.,) Contact Normal of 1 to 2,
     *                            of relative velocity of body 1 relative to 2.
     */
    bool find(const Vec3& Zdir);


#ifndef UNIT_TESTS
  private:
#else
  public:
#endif

    double averageDepth();

    void rotateWorld();

    bool findInitial2Simplex();

    bool findInitial2SimplexFrom3Simplex();

    bool findInitial2SimplexFrom1Simplex();

    bool findInitial2SimplexFrom0Simplex();

    bool isOriginInsideTriangleXY(
        const Vec3& p1,
        const Vec3& p2,
        const Vec3& p3
    );

    bool isOriginInsideTriangleXYrelaxed(
        const Vec3& p1,
        const Vec3& p2,
        const Vec3& p3
    );

    Vec3 findSupportDirectionOfTriangle();

    bool findExtremePointsAlongSupportDirection(
        const Vec3& dir,
        VertexIt&   vit1,
        VertexIt&   vit2
    );

    void selectNewTriangle(VertexIt vit1, VertexIt vit2);

    double squaredDistanceToOriginXY(const Vec3& p1_3d, const Vec3& p2_3d);

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
    Vec3                     mZdir;
    const long               mMaxNumPivots;
    const long               mMaxNumCycles;
    const double             mEpsilonZero;
    bool                     mUseGeomConfigTemp;
    const double             mScaling;
    Mat3x3                   mRotMat1;
    Vec3                     mCoM1;
    Mat3x3                   mRotMat2;
    Vec3                     mCoM2;

#ifdef UNIT_TESTS

    enum stepType {
        ST_INIT,
        ST_AFTER_rotateWorld,
        ST_BEFORE_loop,
        ST_BEGIN_loop,
        ST_AFTER_findSupportDirectionOfTriangle,
        ST_AFTER_findExtremePointsAlongSupportDirection,
        ST_AFTER_selectNewTriangle,
        ST_AFTER_loop,
        ST_RETURN_TYPE1,
        ST_RETURN_TYPE2,
        ST_RETURN_TYPE3,
        ST_TERMINATED
    }              mST_stepType;
    long           mST_stepNo;
    bool           mST_found;
    Vec3           mST_dir;
    long           mST_numIter;
    VertexIt       mST_vit1, mST_vit2;

    void findSequencerInit();
    void findSetParam(const Vec3& Zdir);
    void findSequencerOneStep();

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

    void renderDepartingDirLine(
        Vec3&         color,
        vector<Vec3>& vertices,
        vector<Vec3>& colors
    );

#endif

};


void BDBoundarySimplexFinder::logQ(enum LogLevel lvl)
{
    if (mLogLevel>0 && mLogLevel>=lvl) {
        mLogStream << "Q size: " << mQ.size() << "\n";

        for (int i = 0; i < mQ.size() ; i++) {

            auto vit1 = mQ[i].v1();
            auto vit2 = mQ[i].v2();
            auto id1  = (*vit1)->id();
            auto id2  = (*vit2)->id();
            auto p1 = (*vit1)->pGCS(mRotMat1, mCoM1);
            auto p2 = (*vit2)->pGCS(mRotMat2, mCoM2);

            mLogStream << "Enrty [" << i << "] "
                       << "V1[" << id1 << ": " << p1 << " "
                       << "V2[" << id2 << ": " << p2 << " "
                       << "P: " <<  mQ[i].p() << "\n";
        }
        mLogStream << "\n";
    }
}


void BDBoundarySimplexFinder::logPoints(
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


void BDBoundarySimplexFinder::logPoints(
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


void BDBoundarySimplexFinder::logPoints(
    enum LogLevel lvl,
    const Vec3&   p1, 
    const Vec3&   p2
) {
    if (mLogLevel>0 && mLogLevel>=lvl) {
        mLogStream << p1 << "\n";
        mLogStream << p2 << "\n";
    }
}


void BDBoundarySimplexFinder::logPoint(
    enum LogLevel lvl,
    const Vec3&   p1
) {
    if (mLogLevel>0 && mLogLevel>=lvl) {
        mLogStream << p1 << "\n";
    }
}


void BDBoundarySimplexFinder::logPredicate(
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

}// namespace Makena

#endif /*_MAKENA_BD_BOUNDARY_SIMPLEX_FINDER_HPP_*/

#include "gjk_origin_finder.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file gjk_origin_finder.cpp
 *
 * @brief performs GJK collision detection algorithm to try to find 
 *        the origin in the difference configuration space.
 */
namespace Makena {

using namespace std;


bool GJKOriginFinder::findOrigin(const bool coldStart)
{
    if (coldStart) {

        findInitialSimplex();
    }

    double prevSqDist;
    Vec3   prevSqP;
    long   numCycles = 0;
    vector<pair<long,long> > prevQ;

    for (long numIter = 0; numIter < mMaxNumIter; numIter++) {
        enum predicate pred;
        Vec3           closestP;

        findClosestPointToOrigin(closestP, pred);
        removeUnsupportingVertices(pred);

        double sqDist = closestP.squaredNorm2();

        if (sqDist <= mEpsilonZero) {
#ifdef UNIT_TESTS 
            mST_closestP = closestP;

#endif
            return true;
        }

        if (numIter > 0) {
            if (prevSqDist <= sqDist) {
                if (numCycles > 0 ) {

                    if (mQ.size() == prevQ.size()) {
                        bool matched = true;
                        for (long i = 0; i < mQ.size(); i++) {
                            auto& q = mQ[i];
                            if ( (*(q.v1()))->id() != prevQ[i].first  || 
                                 (*(q.v2()))->id() != prevQ[i].second   ) {
                                matched = false;
                                break;
                            }
                        }
                        if (matched) {
                            logQ(WARNING);
                            log(WARNING, __FILE__, __LINE__, 
                                "Exact cycle detected. Aborting." );
                            return false;
                        }
                    }
                }
                if (numCycles >= mMaxNumCycles) {
                    logQ(WARNING);
                    log(WARNING, __FILE__, __LINE__, 
                        "Exceeding maxNumCycles. Aborting." );
                    return false;

                }

                prevQ.clear();
                for (auto& q : mQ) {
                    prevQ.push_back(make_pair<long,long>((*(q.v1()))->id(),
                                                         (*(q.v2()))->id()));
                }

                numCycles++;
            }
        }

        prevSqDist = sqDist;
        prevSqP    = closestP; 

        VertexIt vit1, vit2;
        bool found=findExtremePointsAlongSupportDirection(closestP,vit1,vit2);
                                                        
        if (!found) {
#ifdef UNIT_TESTS 
            mST_closestP = closestP;
#endif
            return false;
        }

        if (numIter + 1 == mMaxNumIter) {
#ifdef UNIT_TESTS 
            mST_closestP = closestP;
#endif
            return false;
        }

        updateSimplexTowardOrigin(pred, vit1, vit2);


    }

    log(WARNING, __FILE__, __LINE__, "Num iteration exceeded. ");
    logQ(WARNING);
    return false;

}


void GJKOriginFinder::findClosestPointToOrigin(
    Vec3&                    closestP,
    enum predicate&          pred
) {
    switch(mQ.size()) {

      case 4:
        findClosestPointToOrigin3Simplex(closestP, pred);
        break;

      case 3:
        findClosestPointToOrigin2Simplex(closestP, pred);
        break;

      case 2:
        findClosestPointToOrigin1Simplex(closestP, pred);
        break;

      case 1:
        findClosestPointToOrigin0Simplex(closestP, pred);
        break;

      default:
        log(ERROR, __FILE__, __LINE__, "Q size wrong [%ld]", mQ.size());
        logQ(ERROR);

    }
}


void GJKOriginFinder::findClosestPointToOrigin3Simplex(
    Vec3&           closestP,
    enum predicate& pred
) {

    const Vec3 zero(0.0, 0.0, 0.0);

    /*
     * It is assumed that the points are ordered such that
     * (p1p2 cross p1p3) dot p1p4 > 0, i.e., p4 is on the side of the normal
     * of the plane 1->2->3 in the right hand coordinate system.
     *
     *              * p3
     *             /^\
     *            / | \
     *           /  |  \
     *          /   |   \
     *      p4 * . .|. . * p1
     *          \   |   /
     *           \  |  /
     *            \ | /
     *             \|/
     *              * p2
     *
     */
    const auto& p1 = mQ[0].p();
    const auto& p2 = mQ[1].p();
    const auto& p3 = mQ[2].p();
    const auto& p4 = mQ[3].p();

    pred = mGraph.find(p1, p2, p3,p4, zero);

    switch (pred) {

      case VORONOI_INSIDE_TETRAHEDRON:
        closestP = zero;
        break;

      case VORONOI_OVER_TRIANGLE_1_3_2:
        projectPointOnPlane(p1, p3, p2, zero, closestP);
        break;

      case VORONOI_OVER_TRIANGLE_1_2_4:
        projectPointOnPlane(p1, p2, p4, zero, closestP);
        break;

      case VORONOI_OVER_TRIANGLE_1_4_3:
        projectPointOnPlane(p1, p4, p3, zero, closestP);
        break;

      case VORONOI_OVER_TRIANGLE_2_3_4:
        projectPointOnPlane(p2, p3, p4, zero, closestP);
        break;

      case VORONOI_OVER_EDGE_1_2:
        closestP = projectPointOnLine(p1, p2, zero);
        break;

      case VORONOI_OVER_EDGE_2_3:
        closestP = projectPointOnLine(p2, p3, zero);
        break;

      case VORONOI_OVER_EDGE_3_1:
        closestP = projectPointOnLine(p3, p1, zero);
        break;

      case VORONOI_OVER_EDGE_1_4:
        closestP = projectPointOnLine(p1, p4, zero);
        break;

      case VORONOI_OVER_EDGE_2_4:
        closestP = projectPointOnLine(p2, p4, zero);
        break;

      case VORONOI_OVER_EDGE_3_4:
        closestP = projectPointOnLine(p3, p4, zero);
        break;

      case VORONOI_OVER_VERTEX_1:
        closestP = p1;
        break;

      case VORONOI_OVER_VERTEX_2:
        closestP = p2;
        break;

      case VORONOI_OVER_VERTEX_3:
        closestP = p3;
        break;

      case VORONOI_OVER_VERTEX_4:
        closestP = p4;
        break;

      case NONE:
      default:
        closestP = zero;
        log(ERROR, __FILE__, __LINE__, "Unknown predicate");
        logPredicate(ERROR, pred);
        logQ(ERROR);
    }
}


void GJKOriginFinder::findClosestPointToOrigin2Simplex(
    Vec3&           closestP,
    enum predicate& pred
) {

    const Vec3 zero(0.0, 0.0, 0.0);

    const auto& p1 = mQ[0].p();
    const auto& p2 = mQ[1].p();
    const auto& p3 = mQ[2].p();

    pred = testPointAgainstTriangle(p1, p2, p3, zero);

    switch (pred) {

      case VORONOI_INSIDE_TRIANGLE:
        projectPointOnPlane(p1, p2, p3, zero, closestP);
        break;

      case VORONOI_OVER_EDGE_1_2:
        closestP = projectPointOnLine(p1, p2, zero);
        break;

      case VORONOI_OVER_EDGE_2_3:
        closestP = projectPointOnLine(p2, p3, zero);
        break;

      case VORONOI_OVER_EDGE_3_1:
        closestP = projectPointOnLine(p3, p1, zero);
        break;

      case VORONOI_OVER_VERTEX_1:
        closestP = p1;
        break;

      case VORONOI_OVER_VERTEX_2:
        closestP = p2;
        break;

      case VORONOI_OVER_VERTEX_3:
        closestP = p3;
        break;

      case NONE:
      default:
        closestP = zero;
    }

}


void GJKOriginFinder::findClosestPointToOrigin1Simplex(
    Vec3&                    closestP,
    enum predicate&          pred
) {
    const Vec3 zero(0.0, 0.0, 0.0);

    const auto& p1 = mQ[0].p();
    const auto& p2 = mQ[1].p();

    pred = testPointAgainstEdge(p1, p2, zero, closestP);

}


void GJKOriginFinder::findClosestPointToOrigin0Simplex(
    Vec3&                    closestP,
    enum predicate&          pred
) {

    const auto& p1 = mQ[0].p();

    pred = ON_POINT1;

    closestP = p1;
}


void GJKOriginFinder::removeUnsupportingVertices(const enum predicate pred)
{

    vector<BDVertex> Qnew;

    switch(mQ.size()) {

      case 4:

        switch (pred) {

          case VORONOI_INSIDE_TETRAHEDRON:
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[1]);
            Qnew.push_back(mQ[2]);
            Qnew.push_back(mQ[3]);
            break;

          case VORONOI_OVER_TRIANGLE_1_3_2:
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[1]);
            Qnew.push_back(mQ[2]);
            break;

          case VORONOI_OVER_TRIANGLE_1_2_4:
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[1]);
            Qnew.push_back(mQ[3]);
            break;

          case VORONOI_OVER_TRIANGLE_1_4_3:
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[2]);
            Qnew.push_back(mQ[3]);
            break;

          case VORONOI_OVER_TRIANGLE_2_3_4:
            Qnew.push_back(mQ[1]);
            Qnew.push_back(mQ[2]);
            Qnew.push_back(mQ[3]);
            break;

          case VORONOI_OVER_EDGE_1_2:
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[1]);
            break;

          case VORONOI_OVER_EDGE_2_3:
            Qnew.push_back(mQ[1]);
            Qnew.push_back(mQ[2]);
            break;

          case VORONOI_OVER_EDGE_3_1:
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[2]);
            break;

          case VORONOI_OVER_EDGE_1_4:
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[3]);
            break;

          case VORONOI_OVER_EDGE_2_4:
            Qnew.push_back(mQ[1]);
            Qnew.push_back(mQ[3]);
            break;

          case VORONOI_OVER_EDGE_3_4:
            Qnew.push_back(mQ[2]);
            Qnew.push_back(mQ[3]);
            break;

          case VORONOI_OVER_VERTEX_1:
            Qnew.push_back(mQ[0]);
            break;

          case VORONOI_OVER_VERTEX_2:
            Qnew.push_back(mQ[1]);
            break;

          case VORONOI_OVER_VERTEX_3:
            Qnew.push_back(mQ[2]);
            break;

          case VORONOI_OVER_VERTEX_4:
            Qnew.push_back(mQ[3]);
            break;

          default:
            break;
        }
        break;

      case 3:

        switch (pred) {

          case VORONOI_INSIDE_TRIANGLE:
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[1]);
            Qnew.push_back(mQ[2]);
            break;

          case VORONOI_OVER_EDGE_1_2:
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[1]);
            break;

          case VORONOI_OVER_EDGE_2_3:
            Qnew.push_back(mQ[1]);
            Qnew.push_back(mQ[2]);
            break;

          case VORONOI_OVER_EDGE_3_1:
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[2]);
            break;

          case VORONOI_OVER_VERTEX_1:
            Qnew.push_back(mQ[0]);
            break;

          case VORONOI_OVER_VERTEX_2:
            Qnew.push_back(mQ[1]);
            break;

          case VORONOI_OVER_VERTEX_3:
            Qnew.push_back(mQ[2]);
            break;

          default:
            break;

        }
        break;

      case 2:

        switch (pred) {

          case BETWEEN_1_AND_2:
            Qnew.push_back(mQ[0]);
            Qnew.push_back(mQ[1]);
            break;

          case ON_POINT1:
            Qnew.push_back(mQ[0]);
            break;

          case ON_POINT2:
            Qnew.push_back(mQ[1]);
            break;

          default:
            break;

        }
        break;

      case 1:
        Qnew.push_back(mQ[0]);
        break;

      default:
        break;

    }

    mQ = std::move(Qnew);
}


bool GJKOriginFinder::findExtremePointsAlongSupportDirection(
    const Vec3& closestP,
    VertexIt&   vit1,
    VertexIt&   vit2
) {

    // Support direction is in the direction from body1 to body2.
    Vec3 dir = closestP;
    dir.normalize();
    dir.scale(-1.0);

    // Find the extremal vertex in the current Q.
    double prj;

    for (long i = 0; i < mQ.size(); i++) {
        if (i == 0) {
            vit1 = mQ[i].v1();
            prj  = mUseGeomConfigTemp?
                       (dir.dot((*vit1)->pGCS(mBody1.QmatTemp()))):
                       (dir.dot((*vit1)->pGCS(mBody1.Qmat())));
        }
        else {
            VertexIt cit = mQ[i].v1();
            double curPrj= mUseGeomConfigTemp?
                               (dir.dot((*cit)->pGCS(mBody1.QmatTemp()))):
                               (dir.dot((*cit)->pGCS(mBody1.Qmat())));
            if (prj < curPrj) {
                vit1 = cit;
                prj  = curPrj;
            }
        }
    }

    // Find the extremal vertex in manifold 1.
    bool updated   = true;
    while (updated) {
        updated = false;
        VertexIt newVit1;
        for (auto heit : (*vit1)->halfEdges()) {
            if ((*heit)->src()==vit1) {
                VertexIt adjit  = (*heit)->dst();
                double   curPrj = mUseGeomConfigTemp?
                                (dir.dot((*adjit)->pGCS(mBody1.QmatTemp()))):
                                (dir.dot((*adjit)->pGCS(mBody1.Qmat())));
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
            prj  = mUseGeomConfigTemp?
                        (dir.dot((*vit2)->pGCS(mBody2.QmatTemp()))):
                        (dir.dot((*vit2)->pGCS(mBody2.Qmat())));
        }
        else {
            VertexIt cit = mQ[i].v2();
            double curPrj= mUseGeomConfigTemp?
                               (dir.dot((*cit)->pGCS(mBody2.QmatTemp()))):
                               (dir.dot((*cit)->pGCS(mBody2.Qmat())));
            if (prj > curPrj) {
                vit2 = cit;
                prj  = curPrj;
            }
        }
    }

    // Find the extremal vertex in manifold 2.
    updated   = true;
    while (updated) {
        updated = false;
        VertexIt newVit2;
        for (auto heit : (*vit2)->halfEdges()) {
            if ((*heit)->src()==vit2) {
                VertexIt adjit  = (*heit)->dst();
                double   curPrj = mUseGeomConfigTemp?
                                 (dir.dot((*adjit)->pGCS(mBody2.QmatTemp()))):
                                 (dir.dot((*adjit)->pGCS(mBody2.Qmat())));
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

    return true;

}


void GJKOriginFinder::updateSimplexTowardOrigin(
    const enum predicate pred,
    VertexIt&            vit1,
    VertexIt&            vit2
) {

    BDVertex newV(vit1, vit2);
    auto& ch1 = mBody1.ConvexHull();
    auto& ch2 = mBody2.ConvexHull();
    if (mUseGeomConfigTemp) {
        newV.updateP( ch1, mScaling, mBody1.QmatTemp(), mBody1.CoMTemp(),
                      ch2, mScaling, mBody2.QmatTemp(), mBody2.CoMTemp() );
    }
    else {
        newV.updateP( ch1, mScaling, mBody1.Qmat(), mBody1.CoM(),
                      ch2, mScaling, mBody2.Qmat(), mBody2.CoM() );
    }
    mQ.push_back(newV);

    if (mQ.size()==4) {
        if (isInWrongOrderForTetrahedron()) {
            std::swap(mQ[1], mQ[2]);
        }
    }
}

#ifdef UNIT_TESTS


void GJKOriginFinder::findOriginSequencerInit(const bool coldStart) {
    mST_stepNo    = 0;
    mST_stepType  = ST_INIT;
    mST_coldStart = coldStart;
}


void GJKOriginFinder::findOriginSequencerOneStep() {

    switch (mST_stepType) {

      case ST_INIT:

        if (mST_coldStart) {
            findInitialSimplex();
        }
        mST_stepType = ST_BEFORE_loop;
        break;

      case ST_BEFORE_loop:
        mST_numIter   = 0;
        mST_numCycles = 0;
        mST_stepType  = ST_BEGIN_loop;
        break;

      case ST_BEGIN_loop:
        findClosestPointToOrigin(mST_closestP, mST_pred);
        mST_stepType = ST_AFTER_findClosestPointToOrigin;
        break;

      case ST_AFTER_findClosestPointToOrigin:
        removeUnsupportingVertices(mST_pred);            
        mST_sqDist = mST_closestP.squaredNorm2();

        if (mST_numIter > 0) {
            if (mST_prevSqDist < mST_sqDist) {

               if (mST_numCycles > 0) {
                    if (mQ.size() == mST_prevQ.size()) {
                        bool matched = true;
                        for (long i = 0; i < mQ.size(); i++) {
                            auto& q = mQ[i];
                            if ( (*(q.v1()))->id() != mST_prevQ[i].first  || 
                                 (*(q.v2()))->id() != mST_prevQ[i].second   ) {
                                matched = false;
                                break;
                            }
                        }
                        if (matched) {
                            mST_stepType = ST_RETURN_TYPE4;
                        }
                    }
                }
                if (mST_numCycles >= mMaxNumCycles) {
                    mST_stepType = ST_RETURN_TYPE5;
                }
                else {
                    mST_stepType = ST_AFTER_removeUnsupportingVertices;
                }

                mST_prevQ.clear();
                for (auto& q : mQ) {
                    mST_prevQ.push_back(make_pair<long,long>((*(q.v1()))->id(),
                                                             (*(q.v2()))->id()));
                }
                mST_numCycles++;
            }
            else {    
                mST_stepType = ST_AFTER_removeUnsupportingVertices;
            }
        }
        else {
            mST_stepType = ST_AFTER_removeUnsupportingVertices;
        }
        mST_prevSqDist = mST_sqDist;
        mST_prevSqP    = mST_closestP; 
        break;

      case ST_AFTER_removeUnsupportingVertices:
        mST_found = findExtremePointsAlongSupportDirection(
                                             mST_closestP, mST_vit1, mST_vit2);
        if (!mST_found) {
            mST_stepType = ST_RETURN_TYPE2;
        }
        else {
            mST_stepType = ST_AFTER_findExtremePointsAlongSupportDirection;
        }
        break;

      case ST_AFTER_findExtremePointsAlongSupportDirection:
        updateSimplexTowardOrigin(mST_pred, mST_vit1, mST_vit2);
        mST_stepType =  ST_AFTER_updateSimplexTowardOrigin;
        break;

      case ST_AFTER_updateSimplexTowardOrigin:
        mST_numIter++;
        mST_stepType = ST_BEGIN_loop;
        break;

      case ST_AFTER_loop:
        mST_stepType = ST_RETURN_TYPE3;
        break;

      default:
        break;
    }
}


void GJKOriginFinder::renderObj1Surface(
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
        mBody1.Qmat(),
        mBody1.CoM());
}



void GJKOriginFinder::renderObj2Surface(
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
        mBody2.Qmat(),
        mBody2.CoM());
}


void GJKOriginFinder::renderObj1Lines(
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
        mBody1.Qmat(),
        mBody1.CoM()
    );
}


void GJKOriginFinder::renderObj2Lines(
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
        mBody2.Qmat(),
        mBody2.CoM()
    );
}


void GJKOriginFinder::renderBDManifoldSurface(
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
        auto p1 =  mUseGeomConfigTemp?
                       ((*vit1)->pGCS(mBody1.QmatTemp(), mBody1.CoMTemp())):
                       ((*vit1)->pGCS(mBody1.Qmat(), mBody1.CoM()));
        for (auto vit2 = Vpair2.first; vit2 != Vpair2.second; vit2++) {
            auto p2 = mUseGeomConfigTemp?
                       ((*vit2)->pGCS(mBody2.QmatTemp(), mBody2.CoMTemp())):
                       ((*vit2)->pGCS(mBody2.Qmat(), mBody2.CoM()));
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


void GJKOriginFinder::renderBDManifoldLines(
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
                      ((*vit1)->pGCS(mBody1.Qmat(), mBody1.CoM()));
        for (auto vit2 = Vpair2.first; vit2 != Vpair2.second; vit2++) {
            auto p2 = mUseGeomConfigTemp?
                       ((*vit2)->pGCS(mBody2.QmatTemp(), mBody2.CoMTemp())):
                       ((*vit2)->pGCS(mBody2.Qmat(), mBody2.CoM()));
            diffPoints.push_back(p1 - p2);
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


void GJKOriginFinder::renderActiveFeature1Points(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
) {
    vector<Vec3> points;
    for (auto& q : mQ) {
        auto vit = q.v1();
        auto p1 = mUseGeomConfigTemp?
                      ((*vit)->pGCS(mBody1.QmatTemp(), mBody1.CoMTemp())):
                      ((*vit)->pGCS(mBody1.Qmat(), mBody1.CoM()));
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


void GJKOriginFinder::renderActiveFeature2Points(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
) {
    vector<Vec3> points;
    for (auto& q : mQ) {
        auto vit = q.v2();
        auto p1 = mUseGeomConfigTemp?
                      ((*vit)->pGCS(mBody2.QmatTemp(), mBody2.CoMTemp())):
                      ((*vit)->pGCS(mBody2.Qmat(), mBody2.CoM()));
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


void GJKOriginFinder::renderActiveFeatureBDManifoldPoints(
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
                      ((*vit1)->pGCS(mBody1.Qmat(), mBody1.CoM()));
        auto p2 = mUseGeomConfigTemp?
                      ((*vit2)->pGCS(mBody2.QmatTemp(), mBody2.CoMTemp())):
                      ((*vit2)->pGCS(mBody2.Qmat(), mBody2.CoM()));
        points.push_back(p1-p2);
    }

    Manifold::makeOpenGLVerticesColorsForPoints(
        points,
        color,
        vertices,
        colors
    );
}


void GJKOriginFinder::renderSupportVectorLine(
    Vec3&         color,
    vector<Vec3>& vertices,
    vector<Vec3>& colors
) {
    vector<Vec3> points;
    Vec3 n = mST_closestP;
    n.normalize();
    n.scale(-1.0);
    points.push_back(Vec3(0.0, 0.0, 0.0));
    points.push_back(n);

    Manifold::makeOpenGLVerticesColorsForPoints(
                                              points, color, vertices, colors);
}

void GJKOriginFinder::renderClosestPoint(
    Vec3&               color,
    vector<Vec3>&       vertices,
    vector<Vec3>&       colors
) {
    vertices.push_back(mST_closestP);
    colors.push_back(color);
}

#endif

}// namespace Makena

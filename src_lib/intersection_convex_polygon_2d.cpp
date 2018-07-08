#include <iomanip>
#include "intersection_convex_polygon_2d.hpp"

/**
 * @file intersection_convex_polygon_2d.cpp
 *
 */

namespace Makena {

using namespace std;


bool IntersectionFinderConvexPolygon2D::findIntersection_point_point(
    vector<InputElem*>& A,
    vector<InputElem*>& B,
    vector<OutputElem>& intsec
) {
    const Vec2 vAB = A[0]->mP - B[0]->mP;

    if (vAB.squaredNorm2() <= mEpsilonSquared) {

        OutputElem e(
            A[0]->mP,
            A[0]->mIndex,
            -1,
            B[0]->mIndex,
            -1,
            IT_VERTEX_VERTEX
        );
        intsec.push_back(e);

        return true;
    }
    else {
        return false;
    }
}


bool IntersectionFinderConvexPolygon2D::findIntersection_point_edge(
    vector<InputElem*>& A,
    vector<InputElem*>& B,
    vector<OutputElem>& intsec
) {
    if (!boxTest(B[0]->mP, B[1]->mP, A[0]->mP)) {
        return false;
    }

    // Point and an edge. Middle point of the point and the
    // perpendicular projection point onto the edge. 
    double s;
    Vec2 pProjA = findPerpendicularProjectionOfPointOntoEdge(
                                              A[0]->mP, B[0]->mP, B[1]->mP, s);

    Vec2 vDiff = A[0]->mP - pProjA;

    if (vDiff.squaredNorm2() <= mEpsilonSquared) {

        // A[0] is on the edge (B[0] , B[1])

        if (s <= mEpsilonSquared) {

            // A[0] is coincident to B[0]
            intsec.push_back(OutputElem(
                B[0]->mP,
                A[0]->mIndex,
                -1,
                B[0]->mIndex,
                -1,
                IT_VERTEX_VERTEX
            ));
        }
        else if (1.0 - mEpsilonSquared <= s) {

            // A[0] is coincident to B[1]
            intsec.push_back(OutputElem(
                B[1]->mP,
                A[0]->mIndex,
                -1,
                B[1]->mIndex,
                -1,
                IT_VERTEX_VERTEX
            ));
        }
        else {

            // A[0] is interior of (B[0], B[1])
            intsec.push_back(OutputElem(
                pProjA,
                A[0]->mIndex,
                -1,
                B[0]->mIndex,
                B[1]->mIndex,
                IT_VERTEX_EDGE
            ));
        }
        return true;
    }
    else {
        return false;
    }

}


bool IntersectionFinderConvexPolygon2D::findIntersection_edge_point(
    vector<InputElem*>& A,
    vector<InputElem*>& B,
    vector<OutputElem>& intsec
) {
    if (!boxTest(A[0]->mP, A[1]->mP, B[0]->mP)) {
        return false;
    }

    // Point and an edge. Middle point of the point and the
    // perpendicular projection point onto the edge. 
    double s;
    Vec2 pProjB = findPerpendicularProjectionOfPointOntoEdge(
                                              B[0]->mP, A[0]->mP, A[1]->mP, s);

    Vec2 vDiff = B[0]->mP - pProjB;

    if (vDiff.squaredNorm2() <= mEpsilonSquared) {

        // B[0] is on the edge (A[0] , A[1])

        if (s <= mEpsilonSquared) {

            // B[0] is coincident to A[0]
            intsec.push_back(OutputElem(
                A[0]->mP,
                A[0]->mIndex,
                -1,
                B[0]->mIndex,
                -1,
                IT_VERTEX_VERTEX
            ));
        }
        else if (1.0 - mEpsilonSquared <= s) {

            // B[0] is coincident to A[1]
            intsec.push_back(OutputElem(
                A[1]->mP,
                A[1]->mIndex,
                -1,
                B[0]->mIndex,
                -1,
                IT_VERTEX_VERTEX
            ));
        }
        else {

            // B[0] is interior of (A[0], A[1])
            intsec.push_back(OutputElem(
                A[1]->mP,
                A[0]->mIndex,
                A[1]->mIndex,
                B[0]->mIndex,
                -1,
                IT_EDGE_VERTEX
            ));
        }
        return true;
    }
    else {
        return false;
    }

}


bool IntersectionFinderConvexPolygon2D::findIntersection_point_polygon(
    vector<InputElem*>& A,
    vector<InputElem*>& B,
    vector<OutputElem>& intsec
) {

    long lengthB = B.size();

    const Vec2& pA = A[0]->mP;
    bool  properlyInside = true;

    // Perform orientation test to see if pA is on the port (left) side 
    // of all the edges of B.
    for (size_t i = 0 ; i < lengthB; i++) {

        const size_t indexB1 = i;
        const size_t indexB2 = (i+1)%lengthB;

        const Vec2& pB1 = B[indexB1]->mP;
        const Vec2& pB2 = B[indexB2]->mP;
        Vec2   v12     = pB2 - pB1;
        Vec2   v1p     = pA  - pB1;
        Vec2   v12perp = v12.perp();
        double s       = v12perp.dot(v1p);

        if (s <= -0.5 * mEpsilonSquared) {
            // The factor 0.5 above is to avoid failure to 
            // capture edge incidence below.

            // A[0] is outside of (B[i],B[i+1]). Early out.
            return false;

        }
        else if  (s <= 2.0 * mEpsilonSquared) {
            // The factor 2.0 above is to make sure to run the edge incidence 
            // test if it is sufficiently close.

            // A[0] is pretty close to (B[i], B[i+1]).
            // Better perform edge incidence tests below.
            properlyInside = false;

        }

    }
    

    if (properlyInside) {

        // A[0] is properly inside of B with margin mEpsilonSquared
        intsec.push_back(OutputElem(
            pA,
            A[0]->mIndex,
            -1,
            -1,
            -1,
            IT_VERTEX_INTERIOR
        ));
        return true;

    }


    // Check for incidence of pA on any of the edges on B.
    for (size_t i = 0 ; i < lengthB; i++) {

        const size_t indexB1 = i;
        const size_t indexB2 = (i+1)%lengthB;

        const Vec2& pB1 = B[indexB1]->mP;
        const Vec2& pB2 = B[indexB2]->mP;

        if (!boxTest(pB1, pB2, pA)) {
            continue;
        }

        double s;

        auto  pProjA = findPerpendicularProjectionOfPointOntoEdge(
                                                             pA, pB1, pB2, s);

        if ( (pProjA - pA).squaredNorm2() <= mEpsilonSquared) {

            // A[0] is on the edge (B[i] , B[i+1])

            if (s <= mEpsilonSquared) {

                // A[0] is coincident to B[i]
                intsec.push_back(OutputElem(
                    pB1, 
                    A[0      ]->mIndex,
                    -1,
                    B[indexB1]->mIndex, 
                    -1,
                    IT_VERTEX_VERTEX
                ));

            }
            else if (1.0 - mEpsilonSquared <= s) {

                // A[0] is coincident to B[i+1]
                intsec.push_back(OutputElem(
                    pB2,
                    A[0      ]->mIndex,
                    -1,
                    B[indexB2]->mIndex,
                    -1,
                    IT_VERTEX_VERTEX
                ));

            }
            else {
                // A[0] is interior of (B[i], B[i+1])
                intsec.push_back(OutputElem(
                    pProjA,
                    A[0      ]->mIndex,
                    -1,
                    B[indexB1]->mIndex,
                    B[indexB2]->mIndex,
                    IT_VERTEX_EDGE
                ));

            }
            return true;
        }
    }

    // A[0] is inside of B
    intsec.push_back(OutputElem(
        pA,
        A[0]->mIndex,
        -1,
        -1,
        -1,
        IT_VERTEX_INTERIOR
    ));
                             
    return true;
}


bool IntersectionFinderConvexPolygon2D::findIntersection_polygon_point(
    vector<InputElem*>& A,
    vector<InputElem*>& B,
    vector<OutputElem>& intsec
) {

    long lengthA = A.size();

    const Vec2& pB = B[0]->mP;
    bool  properlyInside = true;

    // Perform orientation test to see if pB is on the port (left) side 
    // of all the edges of A.
    for (size_t i = 0 ; i < lengthA; i++) {
        const size_t indexA1 = i;
        const size_t indexA2 = (i+1)%lengthA;

        const Vec2& pA1 = A[indexA1]->mP;
        const Vec2& pA2 = A[indexA2]->mP;

        Vec2   v12     = pA2 - pA1;
        Vec2   v1p     = pB  - pA1;
        Vec2   v12perp = v12.perp();
        double s       = v12perp.dot(v1p);

        if (s <= -0.5 * mEpsilonSquared) {
            // The factor 0.5 above is to avoid failure to 
            // capture edge incide below.
            // B[0] is outside of (A[i], A[i+1]). Early out.
            return false;

        }
        else if  (s <= 2.0 * mEpsilonSquared) {
            // The factor 2.0 above is to make sure to run the edge incidence 
            // test if it is sufficiently close.

            // B[0] is pretty close to (A[i], A[i+1]).
            // Better perform edge incidence tests below.
            properlyInside = false;

        }

    }
    

    if (properlyInside) {
        // B[0] is properly inside of A with margin mEpsilonSquared
        intsec.push_back(OutputElem(
            pB,
            -1,
            -1,
            B[0]->mIndex,
            -1,
            IT_INTERIOR_VERTEX
        ));
        return true;

    }

    // Check for incidence of pB on any of the edges on A.
    for (size_t i = 0 ; i < lengthA; i++) {

        const size_t indexA1 = i;
        const size_t indexA2 = (i+1)%lengthA;

        const Vec2& pA1 = A[indexA1]->mP;
        const Vec2& pA2 = A[indexA2]->mP;

        if (!boxTest(pA1, pA2, pB)) {
            continue;
        }

        double s;
        auto  pProjB = findPerpendicularProjectionOfPointOntoEdge(
                                                             pB, pA1, pA2, s);

        if ( (pProjB - pB).squaredNorm2() <= mEpsilonSquared) {

            // B[0] is on the edge (A[i] , A[i+1])

            if (s <= mEpsilonSquared) {

                // B[0] is coincident to A[i]
                intsec.push_back(OutputElem(
                    pA1, 
                    A[indexA1]->mIndex,
                    -1,
                    B[0      ]->mIndex, 
                    -1,
                    IT_VERTEX_VERTEX
                ));

            }
            else if (1.0 - mEpsilonSquared <= s) {

                // B[0] is coincident to A[i+1]
                intsec.push_back(OutputElem(
                    pA2,
                    A[indexA2]->mIndex,
                    -1,
                    B[0      ]->mIndex,
                    -1,
                    IT_VERTEX_VERTEX
                ));

            }
            else {

                // B[0] is interior of (A[i], A[i+1])
                intsec.push_back(OutputElem(
                    pProjB,
                    A[indexA1]->mIndex,
                    A[indexA2]->mIndex,
                    B[0      ]->mIndex,
                    -1,
                    IT_EDGE_VERTEX
                ));

            }
            return true;
        }
    }

    // B[0] is inside of A
    intsec.push_back(OutputElem(
        pB,
        -1,
        -1,
        B[0]->mIndex,
        -1,
        IT_INTERIOR_VERTEX
    ));
                             
    return true;
}


bool IntersectionFinderConvexPolygon2D::findIntersection_edge_edge(
    vector<InputElem*>& A,
    vector<InputElem*>& B,
    vector<OutputElem>& intsec
) {
    if (!boxTest(A[0]->mP, A[1]->mP, B[0]->mP, B[1]->mP)) {
        return false;
    }

    CHNode *topA, *bottomA, *topB, *bottomB;
    long topIndexA, bottomIndexA, topIndexB, bottomIndexB;
    constructCycle(A, &topA, &bottomA, true,  topIndexA, bottomIndexA);
    constructCycle(B, &topB, &bottomB, false, topIndexB, bottomIndexB);

    auto found  = sweepALBL(topA, topB, A.size(), B.size());
    if (found) {
        exploreInterior_edge_edge(topA, bottomA, intsec);
    }
    // Check if either top-bottom or bottom-top  are touching.
    else if (topA->isCoincident(bottomB)) {
        intsec.push_back(OutputElem(
            topA->mP,
            topA   ->mIndex,
            -1,
            bottomB->mIndex,
            -1,
            IT_VERTEX_VERTEX
        ));
        found = true;
    }
    else if (topB->isCoincident(bottomA)) {
        intsec.push_back(OutputElem(
            topB->mP,
            bottomA->mIndex,
            -1,
            topB   ->mIndex,
            -1,
            IT_VERTEX_VERTEX
        ));
        found = true;
    }

    destructCycle(topA);
    destructCycle(topB);
    return found;

}


void IntersectionFinderConvexPolygon2D::exploreInterior_edge_edge(
    CHNode*             topA,
    CHNode*             bottomA,
    vector<OutputElem>& intsec
) {
    if (topA->mPeer != nullptr) {
        enum VertexType t;
        if (topA->mType == CHNode::NT_SPLIT_POINT) {
            t = IT_EDGE_VERTEX;
        }
        else if (topA->mPeer->mType == CHNode::NT_SPLIT_POINT) {
            t = IT_VERTEX_EDGE;
        }
        else {
            t = IT_VERTEX_VERTEX;
        }
        intsec.push_back(OutputElem(
            topA->mP,
            topA->mIndex,
            topA->mIndexAux,
            topA->mPeer->mIndex,
            topA->mPeer->mIndexAux,
            t
        ));
    }
    
    for (auto* curA = topA->mNext->mDst;
         curA      != bottomA; 
         curA       = curA->mNext->mDst  ) {
        enum VertexType t;

        if (curA->mType == CHNode::NT_INTERSECTION) {
            t = IT_EDGE_EDGE;
        }
        else if (curA->mType == CHNode::NT_SPLIT_POINT) {
            t = IT_EDGE_VERTEX;
        }
        else if (curA->mPeer->mType == CHNode::NT_SPLIT_POINT) {
            t = IT_VERTEX_EDGE;
        }
        else {
            t = IT_VERTEX_VERTEX;
        }

        intsec.push_back(OutputElem(
            curA->mP,
            curA->mIndex,
            curA->mIndexAux,
            curA->mPeer->mIndex,
            curA->mPeer->mIndexAux,
            t
        ));
    }

    if (bottomA->mPeer != nullptr) {
        enum VertexType t;
        if (bottomA->mType == CHNode::NT_SPLIT_POINT) {
            t = IT_EDGE_VERTEX;
        }
        else if (bottomA->mPeer->mType == CHNode::NT_SPLIT_POINT) {
            t = IT_VERTEX_EDGE;
        }
        else {
            t = IT_VERTEX_VERTEX;
        }

        intsec.push_back(OutputElem(
            bottomA->mP,
            bottomA->mIndex,
            bottomA->mIndexAux,
            bottomA->mPeer->mIndex,
            bottomA->mPeer->mIndexAux,
            t
        ));
    }
}


bool IntersectionFinderConvexPolygon2D::findIntersection_edge_polygon(
    vector<InputElem*>& A,
    vector<InputElem*>& B,
    vector<OutputElem>& intsec
) {
    // First, check if both points of the edge are inside the polygon
    vector<InputElem*> P0;
    P0.push_back(A[0]);
    vector<InputElem*> P1;
    P1.push_back(A[1]);
    vector<OutputElem> intsec_pp0;
    vector<OutputElem> intsec_pp1;

    auto res_pp0 = findIntersection_point_polygon(P0, B, intsec_pp0);
    auto res_pp1 = findIntersection_point_polygon(P1, B, intsec_pp1);
               
    if (res_pp0 && res_pp1) {
        // The edge is inside the polygon.
        intsec.push_back(intsec_pp0[0]);
        intsec.push_back(intsec_pp1[0]);
        return true;
    }


    // One or both incident vertices are outside of the polygon

    CHNode *topA, *bottomA, *topB, *bottomB;
    long topIndexA, bottomIndexA, topIndexB, bottomIndexB;
    constructCycle(A, &topA, &bottomA, true,  topIndexA, bottomIndexA);

    constructCycle(B, &topB, &bottomB, false, topIndexB, bottomIndexB);

    auto resALBL = sweepALBL(topA, topB, A.size(), B.size());
    log(INFO, __FILE__, __LINE__, "After sweepALBL");
    logChains(INFO, __FILE__, __LINE__, topA, topB);

    auto resALBR = sweepALBR(topA, topB, A.size(), B.size());
    log(INFO, __FILE__, __LINE__, "After sweepALBR");
    logChains(INFO, __FILE__, __LINE__, topA, topB);

    if (topIndexA == 0 && res_pp0){
        intsec.push_back(intsec_pp0[0]);
    }
    else if (topIndexA == 1 && res_pp1){
        intsec.push_back(intsec_pp1[0]);
    }

    if (resALBL || resALBR) {
        exploreInterior_edge_polygon(topA, bottomA, intsec);
    }
    destructCycle(topA);
    destructCycle(topB);

    if (bottomIndexA == 0 && res_pp0){
        intsec.push_back(intsec_pp0[0]);
    }
    else if (bottomIndexA == 1 && res_pp1){
        intsec.push_back(intsec_pp1[0]);
    }

    return resALBL || resALBR || res_pp0 || res_pp1;
}


void IntersectionFinderConvexPolygon2D::exploreInterior_edge_polygon(
    CHNode*             topA,
    CHNode*             bottomA,
    vector<OutputElem>& intsec
) {

    for (auto* curA = topA->mNext->mDst;
         curA      != bottomA; 
         curA       = curA->mNext->mDst  ) {

        enum VertexType t;

        if (curA->mType == CHNode::NT_INTERSECTION) {
            t = IT_EDGE_EDGE;
        }
        else if (curA->mType == CHNode::NT_SPLIT_POINT) {
            t = IT_EDGE_VERTEX;
        }
        else {
            t = IT_VERTEX_VERTEX;
        }

        intsec.push_back(OutputElem(
            curA->mP,
            curA->mIndex,
            curA->mIndexAux,
            curA->mPeer->mIndex,
            curA->mPeer->mIndexAux,
            t
        ));
    }
}


bool IntersectionFinderConvexPolygon2D::findIntersection_polygon_edge(
    vector<InputElem*>& A,
    vector<InputElem*>& B,
    vector<OutputElem>& intsec
) {

    // First, check if both points of the edge are inside the polygon
    vector<InputElem*> P0;
    P0.push_back(B[0]);
    vector<InputElem*> P1;
    P1.push_back(B[1]);
    vector<OutputElem> intsec_pp0;
    vector<OutputElem> intsec_pp1;

    auto res_pp0 = findIntersection_polygon_point(A, P0, intsec_pp0);
    auto res_pp1 = findIntersection_polygon_point(A, P1, intsec_pp1);
               
    if (res_pp0 && res_pp1) {
        // The edge is inside the polygon.
        intsec.push_back(intsec_pp0[0]);
        intsec.push_back(intsec_pp1[0]);
        return true;
    }

    // One or both incident vertices are outside of the polygon

    CHNode *topA, *bottomA, *topB, *bottomB;
    long topIndexA, bottomIndexA, topIndexB, bottomIndexB;
    constructCycle(A, &topA, &bottomA, true,  topIndexA, bottomIndexA);
    constructCycle(B, &topB, &bottomB, false, topIndexB, bottomIndexB);

    auto resALBL = sweepALBL(topA, topB, A.size(), B.size());
    log(INFO, __FILE__, __LINE__, "After sweepALBL");
    logChains(INFO, __FILE__, __LINE__, topA, topB);

    auto resARBL = sweepARBL(topA, topB, A.size(), B.size());
    log(INFO, __FILE__, __LINE__, "After sweepARBL");
    logChains(INFO, __FILE__, __LINE__, topA, topB);

    if (topIndexB == 0 && res_pp0){
        intsec.push_back(intsec_pp0[0]);
    }
    else if (topIndexB == 1 && res_pp1){
        intsec.push_back(intsec_pp1[0]);
    }

    if (resALBL || resARBL) {
        exploreInterior_polygon_edge(topB, bottomB, intsec);
    }
    destructCycle(topA);
    destructCycle(topB);

    if (bottomIndexB == 0 && res_pp0){
        intsec.push_back(intsec_pp0[0]);
    }
    else if (bottomIndexB == 1 && res_pp1){
        intsec.push_back(intsec_pp1[0]);
    }

    return resALBL || resARBL || res_pp0 || res_pp1;
}


void IntersectionFinderConvexPolygon2D::exploreInterior_polygon_edge(
    CHNode*             topB,
    CHNode*             bottomB,
    vector<OutputElem>& intsec
) {
    for (auto* curB = topB->mNext->mDst;
         curB      != bottomB; 
         curB       = curB->mNext->mDst  ) {

        enum VertexType t;

        if (curB->mType == CHNode::NT_INTERSECTION) {
            t = IT_EDGE_EDGE;
        }
        else if (curB->mType == CHNode::NT_SPLIT_POINT) {
            t = IT_VERTEX_EDGE;
        }
        else {
            t = IT_VERTEX_VERTEX;
        }

        intsec.push_back(OutputElem(
            curB->mP,
            curB->mPeer->mIndex,
            curB->mPeer->mIndexAux,
            curB->mIndex,
            curB->mIndexAux,
            t
        ));
    }
}


/** @brief finds the intersection of two convex polygons
 *
 *  @param Ain (in): Points of convex hull A arranged in CCW ordering.
 *  @param Bin (in): Points of convex hull B arranged in CCW ordering.
 *
 *  @param intsec (out): Points of the vertices of the intersecion.
 *
 *  @return true if twho convex hulls intersect.
 */
bool IntersectionFinderConvexPolygon2D::findIntersection(
    vector<InputElem>&  A,
    vector<InputElem>&  B,
    vector<OutputElem>& intsec
) {
    intsec.clear();
    // Remove consecutive coincident vertices and colinear edges from the input
    removeDegeneracy(A, mAcleaned);
    removeDegeneracy(B, mBcleaned);

    if (mAcleaned.size()==0 || mBcleaned.size()==0) {
        // This might occur if either A or B is ill-conditioned.
        return false;
    }

    // If one side is a point and the other is either a point or an edge,
    // then they may not be exactly coincident or exactly aligned.
    // In this case, we take the middle point of two features as the 
    // intersection.
    if (mAcleaned.size() == 1) {
        if (mBcleaned.size() == 1) {

            return findIntersection_point_point(mAcleaned, mBcleaned, intsec);

        }
        else if (mBcleaned.size() == 2) {
            return findIntersection_point_edge(mAcleaned, mBcleaned, intsec);

        }
        else {

            return findIntersection_point_polygon(mAcleaned, mBcleaned,intsec);

        }
    }
    else if (mAcleaned.size() == 2) {
        if (mBcleaned.size() == 1) {

            return findIntersection_edge_point(mAcleaned, mBcleaned, intsec);

        }
        else if (mBcleaned.size() == 2) {
            return findIntersection_edge_edge(mAcleaned, mBcleaned, intsec);

        }
        else {

            return findIntersection_edge_polygon(mAcleaned, mBcleaned, intsec);

        }
    }
    else {
        if (mBcleaned.size() == 1) {
            return findIntersection_polygon_point(mAcleaned, mBcleaned,intsec);

        }
        else if (mBcleaned.size() == 2) {
            return findIntersection_polygon_edge(mAcleaned, mBcleaned, intsec);
                                                  
        }
        else {

           return findIntersection_polygon_polygon(mAcleaned,mBcleaned,intsec);
        }
    }

    return false;
}


void IntersectionFinderConvexPolygon2D::constructInnerOutput(
    vector<OutputElem>& intsec,
    vector<InputElem*>& points,
    const bool          innerIsA
) {
    for (long i = 0; i < points.size(); i++) {
        auto* P = points[i];
        if (innerIsA) {
            intsec.push_back(OutputElem(
                P->mP,
                P->mIndex,
                -1,
                -1,
                -1,
                IT_VERTEX_INTERIOR));
        }
        else {
            intsec.push_back(OutputElem(
                P->mP,
                -1,
                -1,
                P->mIndex,
                -1,
                IT_INTERIOR_VERTEX
            ));
        }
    }
}
                        

/** @brief finds the perpendicular projection of a point on to an edge.
 *
 *  @param    edge (in): two points on the edge
 *
 *  @param    pts  (in): points of the edge.
 *
 *  @param    intsec (out): vertices of intersection
 *
 *  @return  true if non empty intersection is found
 */



/** @brief remove consecutive coincident vertices and colinear edges.
 *
 *  @param raw     (in):  Input points
 *
 *  @param cleaned (out): Output points
 */
void IntersectionFinderConvexPolygon2D::removeDegeneracy(
    vector<InputElem>&  raw,
    vector<InputElem*>& cleaned
) {
    cleaned.clear();
    if (raw.size() == 0) {
        return;
    }
    // Remove consecutive coincident points first.
    vector<InputElem*> coincidenceRemoved;
    InputElem& Pfirst = raw[0];
    Vec2 Pprev  = raw[0].mP;
    coincidenceRemoved.push_back(&Pfirst);
    for (size_t i = 1; i < raw.size(); i++) {
        InputElem& P = raw[i];
        const Vec2  V(P.mP - Pprev);
        if (i < raw.size()-1) {
            if (V.squaredNorm2() >= mEpsilonSquared) {
                coincidenceRemoved.push_back(&P);
                Pprev = P.mP;
            }
        }
        else {
            const Vec2 V2 = P.mP - Pfirst.mP;
            if ( V.squaredNorm2()  >= mEpsilonSquared &&
                 V2.squaredNorm2() >= mEpsilonSquared    ) {
                coincidenceRemoved.push_back(&P);
                Pprev = P.mP;
            }
        }
    }

    if (coincidenceRemoved.size()<4)  {
        cleaned = std::move(coincidenceRemoved);
    }
    else {
        // Remove consecutive colinear edges second
        InputElem* Cfirst  = coincidenceRemoved[0];
        InputElem* Csecond = coincidenceRemoved[1];
        Vec2       Cprev   = Csecond->mP;
        Vec2       Vprev   = Csecond->mP - Cfirst->mP; 
        auto       len     = coincidenceRemoved.size();

        for (size_t i = 2; i < len+2; i++) {

            InputElem* C = coincidenceRemoved[i%len];
            const Vec2 V(C->mP - Cprev);
            const Vec2 Vperp = V.perp();
            if (fabs(Vperp.dot(Vprev)) >= mEpsilonAngle) {
                cleaned.push_back(C);
                Cprev = C->mP;
                Vprev = V;
            }
        }
    }
}


void IntersectionFinderConvexPolygon2D::removeDegeneracy(
    vector<OutputElem>& raw,
    vector<OutputElem>& cleaned
) {
    cleaned.clear();
    if (raw.size() == 0) {
        return;
    }
    // Remove consecutive coincident points first.
    vector<OutputElem> coincidenceRemoved;
    OutputElem& Pfirst = raw[0];
    Vec2 Pprev  = raw[0].mP;
    coincidenceRemoved.push_back(Pfirst);
    for (size_t i = 1; i < raw.size(); i++) {
        const OutputElem& P = raw[i];
        const Vec2  V(P.mP - Pprev);
        if (i < raw.size()-1) {
            if (V.squaredNorm2() >= mEpsilonSquared) {
                coincidenceRemoved.push_back(P);
                Pprev = P.mP;
            }
        }
        else {
            const Vec2 V2 = P.mP - Pfirst.mP;
            if ( V.squaredNorm2()  >= mEpsilonSquared &&
                 V2.squaredNorm2() >= mEpsilonSquared    ) {
                coincidenceRemoved.push_back(P);
                Pprev = P.mP;
            }
        }
    }

    if (coincidenceRemoved.size()<4)  {
        cleaned = std::move(coincidenceRemoved);
    }
    else {
        // Remove consecutive colinear edges second
        OutputElem& Cfirst  = coincidenceRemoved[0];
        OutputElem& Csecond = coincidenceRemoved[1];
        Vec2        Cprev   = Csecond.mP;
        Vec2        Vprev   = Csecond.mP - Cfirst.mP; 
        auto        len     = coincidenceRemoved.size();

        for (size_t i = 2; i < len+2; i++) {

            OutputElem& C = coincidenceRemoved[i%len];
            const Vec2 V(C.mP - Cprev);
            const Vec2 Vperp = V.perp();
            if (fabs(Vperp.dot(Vprev)) >= mEpsilonAngle) {
                cleaned.push_back(C);
                Cprev = C.mP;
                Vprev = V;
            }
        }
    }
}


/** @brief construct the cyclic chain of alternating CHNode and CHEdge
 *         for the given convext polygon.
 *
 *  @param pts     (in):  points on the polygon in CCW ordering.
 *
 *  @param topV    (out): vertex of the highest y-coordinate.
 *
 *  @param bottomV (out): vertex of the lowest y-coordinate.
 */
void IntersectionFinderConvexPolygon2D::constructCycle(
    vector<InputElem*>& pts,
    CHNode**            topV,
    CHNode**            bottomV,
    bool                thisIsA,
    long&               topIndex,
    long&               bottomIndex
) {
    vector<CHNode*> nodes;
    double          topY;
    double          topX;
    size_t          topYi;
    double          bottomY;
    double          bottomX;
    size_t          bottomYi;
    for (size_t i = 0; i < pts.size() ; i++) {
        CHNode* n = new CHNode(
                            pts[i]->mP, 
                            thisIsA,
                            pts[i]->mIndex,
                            -1,
                            CHNode::NT_ORIGINAL_VERTEX,
                            mEpsilonSquared
                        );
        nodes.push_back(n);
        if (i == 0) {
            topY     = n->mP.y();
            topX     = n->mP.x();
            bottomY  = n->mP.y();
            bottomX  = n->mP.x();
            topYi    = 0;
            bottomYi = 0;
        }
        else {
            if ((topY + mEpsilonLinear < n->mP.y()) || 
                (fabs(topY - n->mP.y()) <= mEpsilonLinear && topX > (n->mP.x() + mEpsilonLinear)) ) {
                topY  = n->mP.y();
                topX  = n->mP.x();
                topYi = i;
            }
            if (bottomY > (n->mP.y() + mEpsilonLinear) ||
                (fabs(bottomY - n->mP.y()) <= mEpsilonLinear && (bottomX + mEpsilonLinear) < n->mP.x())) {
                bottomY  = n->mP.y();
                bottomX  = n->mP.x();
                bottomYi = i;
            }
        }

    }

    for (size_t i = 0; i < pts.size() ; i++) {

        CHEdge* e   = new CHEdge(mEpsilonAngle);
        CHNode* src = nodes[i];
        CHNode* dst = nodes[(i+1)%pts.size()];
        e->mSrc     = src;
        src->mNext  = e;
        e->mDst     = dst;
        dst->mPrev  = e;
    }

    (*topV)             = nodes[topYi];
    (*topV)->mTop       = true;
    (*bottomV)          = nodes[bottomYi];
    (*bottomV)->mBottom = true;
    topIndex            = topYi;
    bottomIndex         = bottomYi;
}


/** @brief destruct the cyclic chain of CHNode and CHEdge
 *
 *  @param n  (in):  a CHNode on the chain
 */
void IntersectionFinderConvexPolygon2D::destructCycle(CHNode* n)
{
    n->mPrev->mDst = nullptr;
    n->mPrev       = nullptr;
    while (n != nullptr) {
        auto* nn = n->mNext->mDst;
        delete n->mNext;
        delete n;
        n = nn;
    }
}


void IntersectionFinderConvexPolygon2D::logChains(
    enum LogLevel lvl,
    const char*   _file,
    const int     _line,
    CHNode*       na,
    CHNode*       nb
) const {
    if (mLogLevel>0 && mLogLevel>=lvl) {
        auto* cur = na;
        mLogStream << "RES A:\n";
        cur->debugPrint(mLogStream);
        cur = cur->mNext->mDst;
        for (; cur != na; cur= cur->mNext->mDst) {
            cur->debugPrint(mLogStream);
        }
        cur = nb;
        mLogStream << "RES B:\n";
        cur->debugPrint(mLogStream);
        cur = cur->mNext->mDst;
        for (; cur != nb; cur= cur->mNext->mDst) {
            cur->debugPrint(mLogStream);
        }
    }
}


void IntersectionFinderConvexPolygon2D::logDumpInput(
    enum LogLevel lvl,
    const char*   _file,
    const int     _line
) const {

    if (mLogLevel>0 && mLogLevel>=lvl) {
        mLogStream << "DUMPING 2D POINTS A\n";
        for (auto* e : mAcleaned) {
            e->mP.decDump(mLogStream);
            mLogStream << "\n";
        }

        mLogStream << "DUMPING 2D POINTS B\n";
        for (auto* e : mBcleaned) {
            e->mP.decDump(mLogStream);
            mLogStream << "\n";
        }

    }
}


bool IntersectionFinderConvexPolygon2D::findIntersection_polygon_polygon(
    vector<InputElem*>& A,
    vector<InputElem*>& B,
    vector<OutputElem>& intsec
) {
    CHNode *topA, *bottomA, *topB, *bottomB;
    long topIndexA, bottomIndexA, topIndexB, bottomIndexB;
    constructCycle(A, &topA, &bottomA, true,  topIndexA, bottomIndexA);

    constructCycle(B, &topB, &bottomB, false, topIndexB, bottomIndexB);

    auto resALBL = sweepALBL(topA, topB, A.size(), B.size());
    log(INFO, __FILE__, __LINE__, "After sweepALBL");
    logChains(INFO, __FILE__, __LINE__, topA, topB);

    auto resALBR = sweepALBR(topA, topB, A.size(), B.size());
    log(INFO, __FILE__, __LINE__, "After sweepALBR");
    logChains(INFO, __FILE__, __LINE__, topA, topB);

    auto resARBL = sweepARBL(topA, topB, A.size(), B.size());
    log(INFO, __FILE__, __LINE__, "After sweepARBL");
    logChains(INFO, __FILE__, __LINE__, topA, topB);

    auto resARBR = sweepARBR(bottomA, bottomB, A.size(), B.size());
    log(INFO, __FILE__, __LINE__, "After sweepARBR");
    logChains(INFO, __FILE__, __LINE__, topA, topB);

    auto res = resALBL || resALBR || resARBL || resARBR;

    if (res) {
        vector<OutputElem> intsecRaw;
        exploreInterior_polygon_polygon(topA, intsecRaw, A.size(), B.size());
        destructCycle(topA);
        destructCycle(topB);
        removeDegeneracy(intsecRaw, intsec);
        return true;
    }
    else {
        destructCycle(topA);
        destructCycle(topB);

        auto* topElemA    = A[topIndexA];
        auto* bottomElemA = A[bottomIndexA];
        auto* topElemB    = B[topIndexB];
        auto* bottomElemB = B[bottomIndexB];

        vector<InputElem*> vTopElemA;
        vTopElemA.push_back(topElemA);
        vector<InputElem*> vBottomElemA;
        vBottomElemA.push_back(bottomElemA);
        vector<InputElem*> vTopElemB;
        vTopElemB.push_back(topElemB);
        vector<InputElem*> vBottomElemB;
        vBottomElemB.push_back(bottomElemB);

        vector<OutputElem> intsecTopA, intsecBottomA;
        vector<OutputElem> intsecTopB, intsecBottomB; 
        auto resTopA    =
                findIntersection_point_polygon(vTopElemA,    B, intsecTopA);
        auto resBottomA =
                findIntersection_point_polygon(vBottomElemA, B, intsecBottomA);
        auto resTopB    =
                findIntersection_polygon_point(A, vTopElemB,    intsecTopB);
        auto resBottomB =
                findIntersection_polygon_point(A, vBottomElemB, intsecBottomB);

        if (resTopA && resBottomA && !resTopB && !resBottomB) {
            constructInnerOutput(intsec, A, true);
            return true;
        }
        else if (!resTopA && !resBottomA && resTopB && resBottomB) {
            constructInnerOutput(intsec, B, false);
            return true;
        }
        else if (resTopA && !resBottomA && !resTopB && resBottomB) {
            if ( (intsecTopA[0].mP - intsecBottomB[0].mP).squaredNorm2()
                                                        <= mEpsilonSquared ) {
                intsec.push_back(intsecTopA[0]);
            }
            else {
                intsec.push_back(intsecTopA[0]);
                intsec.push_back(intsecBottomB[0]);
            }
            return true;
        }
        else if (!resTopA && resBottomA && resTopB && !resBottomB) {
            if ( (intsecBottomA[0].mP - intsecTopB[0].mP).squaredNorm2()
                                                        <= mEpsilonSquared ) {
                intsec.push_back(intsecBottomA[0]);
            }
            else {
                intsec.push_back(intsecBottomA[0]);
                intsec.push_back(intsecTopB[0]);
            }
            return true;
        }
        else {
            return false;
        }
    }

    return false;
}


/** @brief explores the intertwined chains of A and B along the interior side
 *         and find the vertices of the intersection.
 *
 *  @param topA (in): initial point of the chain on A to explore.
 *
 *  @param intsec (out): vertices of intersection
 */
void IntersectionFinderConvexPolygon2D::exploreInterior_polygon_polygon(
    CHNode*             topA,
    vector<OutputElem>& intsec,
    long                sizeA,
    long                sizeB
) {

    // Find the initial intersection point of A and B along the chain.

    CHNode*   curA   = topA;
    CHNode* startV = nullptr;

    if (curA->mNext->mSrcInterior || curA->mPeer != nullptr) {

        startV = topA;
    }
    else {

        curA = curA->mNext->mDst;
        while(curA != topA) {

            if (curA->mNext->mSrcInterior || curA->mPeer != nullptr) {
                startV = curA;
                break;
            }
            curA = curA->mNext->mDst;       
        }

    }
    if (startV==nullptr) {
        return;// Safety exit. Must not reach here.
    }

    auto* curV                = startV;    
    bool  startIsTopA         = (startV==topA);
    bool  touchingPoint       = false;
    bool  interiorFound       = false;
    bool  finished            = false;
    bool  intsecFound         = false;
    bool  needToVisitBackward = true;

    intsec.push_back(constructOutputElem(curV));
    curV = visitVertexAndAdvanceForward( curV,
                                         finished,
                                         touchingPoint,
                                         interiorFound,
                                         intsecFound,
                                         needToVisitBackward );

    long maxIter = (sizeA + sizeB)*2.0 + 2.0;
    while (!finished && curV != startV && curV->mPeer != startV && maxIter >=0) {

        intsec.push_back(constructOutputElem(curV));
        curV = visitVertexAndAdvanceForward( curV, 
                                             finished, 
                                             touchingPoint, 
                                             interiorFound, 
                                             intsecFound,
                                             needToVisitBackward );
        --maxIter;
    }

    // Scan backward if necessary. This means if the starting
    // valid vertex was topA in the forward search, and it has not
    // found any intersection of proper non-zero area.
    if ( startIsTopA           && 
         !touchingPoint        &&
         !interiorFound        && 
         needToVisitBackward   && 
         curV != startV        && 
         curV != startV->mPeer    ) {

        for ( curV        =  topA->mPrev->mSrc;
              curV->mPeer != nullptr && curV != topA;
              curV        =  curV->mPrev->mSrc       ) {
            intsec.insert(intsec.begin(), constructOutputElem(curV));
        }
    }
}



/** @brief sweeps the left hull of A and the left hull of B from top to bottom
 *         and mark intersections.
 *
 *  @param topA (in): top vertex of the left hull of A
 *
 *  @param topB (in): top vertex of the left hull of B
 *
 *  @return true if a crossing or touching point is found
 */
bool IntersectionFinderConvexPolygon2D::sweepALBL(
    CHNode*  topA, 
    CHNode*  topB,
    long     sizeA,
    long     sizeB
) {
    // Go down until both edges vertically overlap.
//cerr << "top a: "; (topA)->debugPrint(std::cerr);
//cerr << "top b: "; (topB)->debugPrint(std::cerr);

    auto* curB = topB->mNext->mDst;
    while (curB->isHigherThan(topA) && !(curB->mBottom)) {
        curB = curB->mNext->mDst;
    }
    if (curB->isHigherThan(topA)) {
        return false;
    }
    curB = curB->mPrev->mSrc;

    auto* curA = topA->mNext->mDst;
    while (curA->isHigherThan(topB) && !(curA->mBottom)) {
        curA = curA->mNext->mDst;
    }
    if (curA->isHigherThan(topB)) {
        return false;
    }
    curA = curA->mPrev->mSrc;

    bool found = false;
    long maxIter  = (sizeA + sizeB)*2 + 2;
    while (maxIter-- >= 0) {
//cerr << "a: "; curA->debugPrint(std::cerr);
//cerr << "b: "; curB->debugPrint(std::cerr);

        auto mayIntersect = horizontalTest(curA->mNext, curB->mNext);
        if (curA->mBottom || curB->mBottom) {
            if (mayIntersect) {
                if (curA->isCoincident(curB)) {
                    handleCoincidence(&curA, &curB, true, true, true);
                    found = true;
                }
                else if (curA->mBottom || !curB->mBottom) {
                    if (curA->isCoincident(curB->mNext->mDst)) {
                        handleCoincidence(
                                &curA, &(curB->mNext->mDst), true, true, true);
                        found = true;
                    }
                    else if (curA->isOnEdge(curB->mNext)) {
                        handleVertexOnEdge(
                                       &curA, &curB, false, true, true, true);
                        found = true;
                    }
                }
                else if (!curA->mBottom || curB->mBottom) {
                    if (curA->mNext->mDst->isCoincident(curB)) {
                        handleCoincidence(
                                &(curA->mNext->mDst), &curB, true, true, true);
                        found = true;
                    }
                    else if (curB->isOnEdge(curA->mNext)) {
                        handleVertexOnEdge(
                                       &curA, &curB, true,  true, true, true);
                        found = true;
                    }
                }
            }
            break;
        }

        auto* nextA = curA->mNext->mDst;
        auto* nextB = curB->mNext->mDst;

        if (!mayIntersect) {
            updateCurACurBOE(&curA, &curB, true, true, true);
        }
        else if (curA->isCoincident(curB)) {
            handleCoincidence(&curA, &curB, true, true, true);
            found = true;
        }
        else if (curA->isCoincident(nextB)) {
            curB = nextB;
            found = true;
        }
        else if (curB->isCoincident(nextA)) {
            curA = nextA;
            found = true;
        }
        else if (curA->isOnEdge(curB->mNext)) {
            handleVertexOnEdge(&curA, &curB, false, true, true, true);
            found = true;
        }
        else if (curB->isOnEdge(curA->mNext)) {
            handleVertexOnEdge(&curA, &curB, true,  true, true, true);
            found = true;
        }
        else if (nextB->isCoincident(nextA)) {
            curA = nextA;
            curB = nextB;
            found = true;
        }
        else if (nextA->isOnEdge(curB->mNext)) {
            curA = nextA;
            handleVertexOnEdge(&curA, &curB, false, true, true, true);
            found = true;
        }
        else if (nextB->isOnEdge(curA->mNext)) {
            curB = nextB;
            handleVertexOnEdge(&curA, &curB, true,  true, true, true);
            found = true;
        }
        else if (curA->mNext->doesIntersect(curB->mNext)) {
            handleCrossing(&curA, &curB, true, true, true);
            found = true;
        }
        else {
            updateCurACurBOE(&curA, &curB, true, true, true);
        }
    }

    if (maxIter==0) {
        log(ERROR, __FILE__, __LINE__, "Infinite loop aborting");
        logDumpInput(ERROR, __FILE__, __LINE__);
    }
 

    return found;
}


/** @brief sweeps the right hull of A and the right hull of B from bottom to
 *         top and mark intersections.
 *
 *  @param bottomA (in): bottom vertex of the right hull of A
 *
 *  @param bottomB (in): bottom vertex of the right hull of B
 *
 *  @return true if a crossing or touching point is found
 */
bool IntersectionFinderConvexPolygon2D::sweepARBR(
    CHNode*  bottomA, 
    CHNode*  bottomB,
    long     sizeA,
    long     sizeB
) {
    auto* curB = bottomB->mNext->mDst;

    while (!curB->isHigherThan(bottomA) && !(curB->mTop)) {

        curB = curB->mNext->mDst;
    }
    if (!curB->isHigherThan(bottomA)) {
        return false;
    }
    curB = curB->mPrev->mSrc;

    auto* curA = bottomA->mNext->mDst;

    while (!curA->isHigherThan(bottomB) && !(curA->mTop)) {
        curA = curA->mNext->mDst;

    }
    if (!curA->isHigherThan(bottomB)) {
        return false;
    }
    curA = curA->mPrev->mSrc;

    bool found = false;

    long maxIter  = (sizeA + sizeB)*2 + 2;

    while (maxIter-- >= 0) {

        auto mayIntersect = horizontalTest(curA->mNext, curB->mNext);

        if (curA->mTop || curB->mTop) {
            if (mayIntersect) {
                // Check the touching point.
                if (curA->isCoincident(curB)) {
                    handleCoincidence(&curA, &curB, false, false, false);
                    found = true;
                }
                else if (curA->mTop || !curB->mTop) {
                    if (curA->isCoincident(curB->mNext->mDst)) {
                        handleCoincidence(
                             &curA, &(curB->mNext->mDst), false, false, false);
                        found = true;
                    }
                    else if (curA->isOnEdge(curB->mNext)) {
                        handleVertexOnEdge(
                                         &curA, &curB,false,false,false,false);
                        found = true;
                    }
                }
                else if (!curA->mTop || curB->mTop) {
                    if (curA->mNext->mDst->isCoincident(curB)) {
                        handleCoincidence(
                             &(curA->mNext->mDst), &curB, false, false, false);
                        found = true;
                    }
                    else if (curB->isOnEdge(curA->mNext)) {
                        handleVertexOnEdge(
                                        &curA, &curB,true,false,false,false);
                        found = true;
                    }
                }
            }
            break;
        }

        auto* nextA = curA->mNext->mDst;
        auto* nextB = curB->mNext->mDst;

        if (!mayIntersect) {
            updateCurACurBOE(&curA, &curB, false, false, false);
        }
        else if (curA->isCoincident(curB)) {
            handleCoincidence(&curA, &curB, false, false, false);
            found = true;
        }
        else if (curA->isCoincident(nextB)) {
            curB = nextB;
            found = true;
        }
        else if (curB->isCoincident(nextA)) {
            curA = nextA;
            found = true;
        }
        else if (curA->isOnEdge(curB->mNext)) {
            handleVertexOnEdge(&curA, &curB, false, false, false, false);
            found = true;
        }
        else if (curB->isOnEdge(curA->mNext)) {
            handleVertexOnEdge(&curA, &curB, true,  false, false, false);
            found = true;
        }
        else if (nextB->isCoincident(nextA)) {
            curA = nextA;
            curB = nextB;
            found = true;
        }
        else if (nextA->isOnEdge(curB->mNext)) {
            curA = nextA;
            handleVertexOnEdge(&curA, &curB, false, false, false, false);
            found = true;
        }
        else if (nextB->isOnEdge(curA->mNext)) {
            curB = nextB;
            handleVertexOnEdge(&curA, &curB, true,  false, false, false);
            found = true;
        }
        else if (curA->mNext->doesIntersect(curB->mNext)) {
            handleCrossing(&curA, &curB, false, false, false);
            found = true;
        }
        else {
            updateCurACurBOE(&curA, &curB, false, false, false);
        }
    }

    if (maxIter==0) {
        log(ERROR, __FILE__, __LINE__, "Infinite loop aborting");
        logDumpInput(ERROR, __FILE__, __LINE__);
    }

    return found;
}


/** @brief sweeps the left hull of A and the right hull of B from top to bottom
 *         and mark intersections.
 *
 *  @param topA (in): top vertex of the left hull of A
 *
 *  @param topB (in): top vertex of the right hull of B
 *
 *  @return true if a crossing or touching point is found
 */
bool IntersectionFinderConvexPolygon2D::sweepALBR(
    CHNode*  topA, 
    CHNode*  topB,
    long     sizeA,
    long     sizeB
) {
    auto* curB = topB->mPrev->mSrc;
    while (curB->isHigherThan(topA) && !(curB->mBottom)) {
        curB = curB->mPrev->mSrc;
    }
    if (curB->isHigherThan(topA)) {
        return false;
    }
    curB = curB->mNext->mDst;

    auto* curA = topA->mNext->mDst;
    while (curA->isHigherThan(topB) && !(curA->mBottom)) {
        curA = curA->mNext->mDst;
    }
    if (curA->isHigherThan(topB)) {
        return false;
    }
    curA = curA->mPrev->mSrc;

    bool found = false;
    long maxIter  = (sizeA + sizeB)*2 + 2;
    while (maxIter-- >= 0) {
        auto mayIntersect = horizontalTest(curA->mNext, curB->mPrev);

        if (curA->mBottom || curB->mBottom) {
            if (mayIntersect) {
                if (curA->isCoincident(curB)) {
                    handleCoincidence(&curA, &curB, true, false, true);
                    found = true;
                }
                else if (curA->mBottom || !curB->mBottom) {
                    if (curA->isCoincident(curB->mPrev->mSrc)) {
                        handleCoincidence(
                               &curA, &(curB->mPrev->mSrc), true, false, true);
                        found = true;
                    }
                    else if (curA->isOnEdge(curB->mPrev)) {
                        handleVertexOnEdge(
                                      &curA, &curB, false, true, false, true);
                        found = true;
                    }
                }
                else if (!curA->mBottom || curB->mBottom) {
                    if (curA->mNext->mDst->isCoincident(curB)) {
                        handleCoincidence(
                               &(curA->mNext->mDst), &curB, true, false, true);
                        found = true;
                    }
                    else if (curB->isOnEdge(curA->mNext)) {
                        handleVertexOnEdge(
                                      &curA, &curB, true,  true, false, true);
                        found = true;
                    }
                }
            }
            break;
        }

        auto* nextA = curA->mNext->mDst;
        auto* prevB = curB->mPrev->mSrc;

        if (!mayIntersect) {
            updateCurACurBOE(&curA, &curB, true, false, true);
        }
        else if (curA->isCoincident(curB)) {
            handleCoincidence(&curA, &curB, true, false, true);
            found = true;
        }
        else if (curA->isCoincident(prevB)) {
            curB = prevB;
            found = true;
        }
        else if (curB->isCoincident(nextA)) {
            curA = nextA;
            found = true;
        }
        else if (curA->isOnEdge(curB->mPrev)) {
            handleVertexOnEdge(&curA, &curB, false, true, false, true);
            found = true;
        }
        else if (curB->isOnEdge(curA->mNext)) {
            handleVertexOnEdge(&curA, &curB, true,  true, false, true);
            found = true;
        }
        else if (prevB->isCoincident(nextA)) {
            curA = nextA;
            curB = prevB;
            found = true;
        }
        else if (nextA->isOnEdge(curB->mPrev)) {
            curA = nextA;
            handleVertexOnEdge(&curA, &curB, false, true, false, true);
            found = true;
        }
        else if (prevB->isOnEdge(curA->mNext)) {
            curB = prevB;
            handleVertexOnEdge(&curA, &curB, true,  true, false, true);
            found = true;
        }
        else if (curA->mNext->doesIntersect(curB->mPrev)) {
            handleCrossing(&curA, &curB, true, false, true);
            found = true;
        }
        else {
            updateCurACurBOE(&curA, &curB, true, false, true);
        }
    }

    if (maxIter==0) {
        log(ERROR, __FILE__, __LINE__, "Infinite loop aborting");
        logDumpInput(ERROR, __FILE__, __LINE__);
    }

    return found;
}


/** @brief sweeps the right hull of A and the left hull of B from top to bottom
 *         and mark intersections.
 *
 *  @param topA (in): top vertex of the right hull of A
 *
 *  @param topB (in): top vertex of the left hull of B
 *
 *  @return true if a crossing or touching point is found
 */
bool IntersectionFinderConvexPolygon2D::sweepARBL(
    CHNode*  topA, 
    CHNode*  topB,
    long     sizeA,
    long     sizeB
) {

    auto* curB = topB->mNext->mDst;
    while (curB->isHigherThan(topA) && !(curB->mBottom)) {
        curB = curB->mNext->mDst;
    }
    if (curB->isHigherThan(topA)) {
        return false;
    }
    curB = curB->mPrev->mSrc;

    auto* curA = topA->mPrev->mSrc;
    while (curA->isHigherThan(topB) && !(curA->mBottom)) {
        curA = curA->mPrev->mSrc;
    }
    if (curA->isHigherThan(topB)) {
        return false;
    }
    curA = curA->mNext->mDst;

    bool found = false;
    long maxIter  = (sizeA + sizeB)*2 + 2;
    while (maxIter-- >= 0) {
        auto mayIntersect = horizontalTest(curA->mPrev, curB->mNext);
        if (curA->mBottom || curB->mBottom) {
            if (mayIntersect) {
                if (curA->isCoincident(curB)) {
                    handleCoincidence(&curA, &curB, false, true, true);
                    found = true;
                }
                else if (curA->mBottom || !curB->mBottom) {
                    if (curA->isCoincident(curB->mNext->mDst)) {
                        handleCoincidence(
                              &curA, &(curB->mNext->mDst), false, true, true);
                        found = true;
                    }
                    else if (curA->isOnEdge(curB->mNext)) {
                        handleVertexOnEdge(
                                      &curA, &curB, false, false, true, true);
                        found = true;
                    }
                }
                else if (!curA->mBottom || curB->mBottom) {

                    if (curA->mPrev->mSrc->isCoincident(curB)) {
                        handleCoincidence(
                               &(curA->mPrev->mSrc), &curB, false, true, true);
                        found = true;
                    }
                    else if (curB->isOnEdge(curA->mPrev)) {
                        handleVertexOnEdge(
                                      &curA, &curB, true,  false, true, true);
                        found = true;
                    }
                }
            }
            break;
        }


        auto* prevA = curA->mPrev->mSrc;
        auto* nextB = curB->mNext->mDst;

        if (!mayIntersect) {
            updateCurACurBOE(&curA, &curB, false, true, true);
        }
        else if (curA->isCoincident(curB)) {
            handleCoincidence(&curA, &curB, false, true, true);
            found = true;
        }
        else if (curA->isCoincident(nextB)) {
            curB = nextB;
            found = true;
        }
        else if (curB->isCoincident(prevA)) {
            curA = prevA;
            found = true;
        }
        else if (curA->isOnEdge(curB->mNext)) {
            handleVertexOnEdge(&curA, &curB, false, false, true, true);
            found = true;
        }
        else if (curB->isOnEdge(curA->mPrev)) {
            handleVertexOnEdge(&curA, &curB, true,  false, true, true);
            found = true;
        }
        else if (nextB->isCoincident(prevA)) {
            curA = prevA;
            curB = nextB;
            found = true;
        }
        else if (prevA->isOnEdge(curB->mNext)) {
            curA = prevA;
            handleVertexOnEdge(&curA, &curB, false, false, true, true);
            found = true;
        }
        else if (nextB->isOnEdge(curA->mPrev)) {
            curB = nextB;
            handleVertexOnEdge(&curA, &curB, true,  false, true, true);
            found = true;
        }
        else if (curA->mPrev->doesIntersect(curB->mNext)) {
            handleCrossing(&curA, &curB, false, true, true);
            found = true;
        }
        else {
            updateCurACurBOE(&curA, &curB, false, true, true);
        }
    }


    if (maxIter==0) {
        log(ERROR, __FILE__, __LINE__, "Infinite loop aborting");
        logDumpInput(ERROR, __FILE__, __LINE__);
    }

    return found;
}


/*
 * Special care must be taken for the following 16 cases if curA and curB are
 * found to be coincident but either of them already has a peer.
 * In this case, previous assumption marked in the pair will have to be 
 * invalidated.
 *
 * sweepALBL & sweepARBR
 * =====================
 *
 *        CV_10          CV_7          CV_8          CV_9
 *       *     o       *     o       o     *       o     *
 *        \   /...      \   /...      \   /...      \   /...
 *         v v....       v v....       v v....       v v....
 *        +---+...      +---+...      +---+...      +---+...
 *        |   |...      |   |...      |   |...      |   |...
 *        +---+...      +---+...      +---+...      +---+...
 *         / \....       / \....       / \....       / \....
 *        v   v...      v   v...      v   v...      v   v...
 *       *     o       o     *       *     o       o     *
 *
 *                                  (CURRENTLY NOT CONSIDERED)
 *        CV_1_3        CV_1_2       CV_4_1         CV_4_2
 *          *o            *o         *     o       o     *
 *          ||...         ||...       \   /...      \   /...
 *          vv....        vv....       v v....       v v....
 *        +---+...      +---+...      +---+...      +---+...
 *        |   |...      |   |...      |   |...      |   |...
 *        +---+...      +---+...      +---+...      +---+...
 *         / \....       / \....        ||....        ||....
 *        v   v...      v   v...        vv...         vv...
 *       *     o       o     *          *o            o*
 *
 *
 *  curA is coincident to curB    
 *
 *            curB has a peer        curA has a peer
 *           +----------------+    +----------------+  
 *           |                |    |                |
 *           |    *    o      |    |      *    o    |
 *           |    |    |      |    |      |    |    |
 *           |    |    |      |    |      |    |    |
 *           |    VpeerV      |    |      VpeerV    |
 *           |    *====o curB |    | curA *====o    |
 *           |    |     \     |    |     /     |    |
 *           |    |      \    |    |    /      |    |
 *           |    V       V   |    |   V       V    |
 *           |    *curA   o   |    |  *    curBo    |
 *           |                |    |                |
 *           +----------------+    +----------------+
 *
 * sweepALBR
 * =========
 *
 *         CV_7          CV_8
 *       o     *       *.....o
 *        ^   /         \...^
 *         \ v           v./
 *        +---+         +---+
 *        |   |         |   |
 *        +---+         +---+
 *         /.^           ^ \
 *        v...\         /   v
 *       *.....o       o     *
 *
 *
 *  curA is coincident to curB    
 *
 *            curB has a peer        curA has a peer
 *           +----------------+    +----------------+  
 *           |                |    |                |
 *           |    *    o      |    |      *    o    |
 *           |    |    ^      |    |      |    ^    |
 *           |    |    |      |    |      |    |    |
 *           |    Vpeer|      |    |      Vpeer|    |
 *           |    *====o curB |    | curA *====o    |
 *           |    |    ^      |    |     /     ^    |
 *           |    |     \     |    |    /      |    |
 *           |    V      \    |    |   V       |    |
 *           |    *curA   o   |    |  *    curBo    |
 *           |                |    |                |
 *           +----------------+    +----------------+
 *
 * sweepARBL
 * =========
 *
 *         CV_7          CV_8
 *       o.....*       *     o
 *        \...^         ^   /
 *         v./           \ v
 *        +---+         +---+
 *        |   |         |   |
 *        +---+         +---+
 *         ^ \           /.^
 *        /   v         v...\
 *       *     o       o.....*
 *
 *
 *            curB has a peer        curA has a peer
 *           +----------------+    +----------------+  
 *           |                |    |                |
 *           |    *    o      |    |      *    o    |
 *           |    ^    |      |    |      ^    |    |
 *           |    |    |      |    |      |    |    |
 *           |    |peerV      |    |      |peerV    |
 *           |    *====o curB |    | curA *====o    |
 *           |    ^     \     |    |     ^     |    |
 *           |    |      \    |    |    /      |    |
 *           |    |       V   |    |   /       V    |
 *           |    *curA   o   |    |  *    curBo    |
 *           |                |    |                |
 *           +----------------+    +----------------+
 */
void IntersectionFinderConvexPolygon2D::handleCoincidenceAhasPeer(
    enum _classifier cls,
    CHNode**         curA, 
    CHNode**         curB, 
    bool             leftHullA, 
    bool             leftHullB, 
    bool             sweepingDownward
) {

    if ( (!leftHullA && !leftHullB && !sweepingDownward) ||
         ( leftHullA &&  leftHullB &&  sweepingDownward)    ) {

        switch(cls) {

          case TYPE_CV_7:
          case TYPE_CV_9:
          case TYPE_CV_1_2:
            (*curA)->mNext->mSrcInterior = true;
            (*curA)->mPeer->mNext->mSrcInterior = false;
            break;

          case TYPE_CV_8:
          case TYPE_CV_10:
          case TYPE_CV_1_3:
            (*curA)->mNext->mSrcInterior = false;
            (*curA)->mPeer->mNext->mSrcInterior = true;
            break;

          case TYPE_CV_4_1:
          case TYPE_CV_4_2:
            // Currently not considered.
            //(*curA) = splitEdge((*curA)->mNext, *curB, CHNode::NT_SPLIT_POINT);
            break;

          default:
            break;
        }
    }
    if (leftHullA && !leftHullB && sweepingDownward) {

        switch(cls) {
          case TYPE_CV_7:
            (*curA)->mNext->mSrcInterior = true;
            (*curA)->mPrev->mDstInterior = false;
            (*curA)->mPeer->mNext->mSrcInterior = false;
            (*curA)->mPeer->mPrev->mDstInterior = true;
            break;

          case TYPE_CV_8:
            (*curA)->mNext->mSrcInterior = false;
            (*curA)->mPrev->mDstInterior = true;
            (*curA)->mPeer->mNext->mSrcInterior = true;
            (*curA)->mPeer->mPrev->mDstInterior = false;
            break;

          default:
            break;
        }
    }
    else if (!leftHullA && leftHullB && sweepingDownward) {

        switch(cls) {
          case TYPE_CV_7:
            (*curA)->mNext->mSrcInterior = true;
            (*curA)->mPrev->mDstInterior = false;
            (*curA)->mPeer->mNext->mSrcInterior = false;
            (*curA)->mPeer->mPrev->mDstInterior = true;
            break;

          case TYPE_CV_8:
            (*curA)->mNext->mSrcInterior = false;
            (*curA)->mPrev->mDstInterior = true;
            (*curA)->mPeer->mNext->mSrcInterior = true;
            (*curA)->mPeer->mPrev->mDstInterior = false;
            break;

          default:
            break;
        }
    }
    updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);

}


void IntersectionFinderConvexPolygon2D::handleCoincidenceBhasPeer(
    enum _classifier cls,
    CHNode**         curA,
    CHNode**         curB,
    bool             leftHullA, 
    bool             leftHullB, 
    bool             sweepingDownward
) {

    if ( (!leftHullA && !leftHullB && !sweepingDownward) ||
         ( leftHullA &&  leftHullB &&  sweepingDownward)    ) {
        switch(cls) {

          case TYPE_CV_7:
          case TYPE_CV_9:
          case TYPE_CV_1_2:
            (*curB)->mNext->mSrcInterior = false;
            (*curB)->mPeer->mNext->mSrcInterior = true;
            break;

          case TYPE_CV_8:
          case TYPE_CV_10:
          case TYPE_CV_1_3:
            (*curB)->mNext->mSrcInterior = true;
            (*curB)->mPeer->mNext->mSrcInterior = false;
            break;

          case TYPE_CV_4_1:
          case TYPE_CV_4_2:
            // Currently not considered.
            //(*curA) = splitEdge((*curA)->mNext, *curB, CHNode::NT_SPLIT_POINT);
            break;

          default:
            break;
        }
    }
    if (leftHullA && !leftHullB && sweepingDownward) {
        switch(cls) {
          case TYPE_CV_7:
            (*curB)->mNext->mSrcInterior = false;
            (*curB)->mPrev->mDstInterior = true;
            (*curB)->mPeer->mNext->mSrcInterior = true;
            (*curB)->mPeer->mPrev->mDstInterior = false;
            break;

          case TYPE_CV_8:
            (*curB)->mNext->mSrcInterior = true;
            (*curB)->mPrev->mDstInterior = false;
            (*curB)->mPeer->mNext->mSrcInterior = false;
            (*curB)->mPeer->mPrev->mDstInterior = true;
            break;

          default:
            break;
        }
    }
    else if (!leftHullA && leftHullB && sweepingDownward) {
        switch(cls) {
          case TYPE_CV_7:
            (*curB)->mNext->mSrcInterior = false;
            (*curB)->mPrev->mDstInterior = true;
            (*curB)->mPeer->mNext->mSrcInterior = true;
            (*curB)->mPeer->mPrev->mDstInterior = false;
            break;

          case TYPE_CV_8:
            (*curB)->mNext->mSrcInterior = true;
            (*curB)->mPrev->mDstInterior = false;
            (*curB)->mPeer->mNext->mSrcInterior = false;
            (*curB)->mPeer->mPrev->mDstInterior = true;
            break;

          default:
            break;
        }
    }

    updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
}


/** @brief handle a case where the current vertices on A and B coincide.
 *         It marks the intersection information along the chains, and
 *         advance the current vertices.
 *         This is a dispatcher function.
 *
 *  @param curA (in/out):         current vertex of the hull of A
 *
 *  @param curB (in/out):         current vertex of the hull of B
 *
 *  @param leftHullA (in):        true if the hull of A is the left hull.
 *
 *  @param leftHullB (in):        true if the hull of B is the left hull
 *
 *  @param sweepingDownward (in): true if the sweeping direction is downward
 */
void IntersectionFinderConvexPolygon2D::handleCoincidence(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {

    auto cls = classifyCoincidentNodes(*curA, *curB);
    if ((*curA)->mPeer != nullptr) {
        handleCoincidenceAhasPeer(
                    cls, curA, curB, leftHullA, leftHullB, sweepingDownward);
    }
    else if ((*curB)->mPeer != nullptr) {
        handleCoincidenceBhasPeer(
                    cls, curA, curB, leftHullA, leftHullB, sweepingDownward);
    }
    else {
        switch (cls) {
          case TYPE_CV_1_1:
            handleCV_1_1(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          case TYPE_CV_1_2:
            handleCV_1_2(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          case TYPE_CV_1_3:
            handleCV_1_3(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          case TYPE_CV_2_1:
            handleCV_2_1(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          case TYPE_CV_2_2:
            handleCV_2_2(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          case TYPE_CV_3_1:
            handleCV_3_1(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          case TYPE_CV_4_1:
            handleCV_4_1(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          case TYPE_CV_4_2:
            handleCV_4_2(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          case TYPE_CV_5:
            handleCV_5(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          case TYPE_CV_7:
            handleCV_7(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          case TYPE_CV_8:
            handleCV_8(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          case TYPE_CV_9:
            handleCV_9(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          case TYPE_CV_10:
            handleCV_10(curA, curB, leftHullA, leftHullB, sweepingDownward);
            break;
          default:
            break;
        }
    }
}


/** @brief handle a case where the current vertices on A or B is on the
 *         edge of the other.
 *         It marks the intersection information along the chains, and
 *         advance the current vertices.
 *         This is a dispatcher function.
 *
 *  @param curA (in/out): current vertex of the hull of A
 *
 *  @param curB (in/out): current vertex of the hull of B
 *
 *  @param leftHullA (in): true if the hull of A is the left hull.
 *
 *  @param leftHullB (in): true if the hull of B is the left hull
 *
 *  @param sweepingDownward (in): true if the sweeping direction is downward
 */
void IntersectionFinderConvexPolygon2D::handleVertexOnEdge(
    CHNode** curA,
    CHNode** curB,
    bool     edgeOnA,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    enum _classifier cls;
    if (edgeOnA) {
        if ((*curB)->mPeer != nullptr) {
            updateCurACurBOE(
                           curA, curB, leftHullA, leftHullB, sweepingDownward);
            return;
        }
        else if ( ( leftHullA &&  sweepingDownward) ||
             (!leftHullA && !sweepingDownward)   ) {
            cls = classifyVertexOnEdge(*curB, (*curA)->mNext);
        }
        else {
            cls = classifyVertexOnEdge(*curB, (*curA)->mPrev);
        }
    }
    else {

        if ((*curA)->mPeer != nullptr) {
            updateCurACurBOE(
                           curA, curB, leftHullA, leftHullB, sweepingDownward);
            return;
        }
        else if ( ( leftHullB &&  sweepingDownward) ||
             (!leftHullB && !sweepingDownward)   ) {
            cls = classifyVertexOnEdge(*curA, (*curB)->mNext);
        }
        else {
            cls = classifyVertexOnEdge(*curA, (*curB)->mPrev);
        }
    }

    switch (cls) {

      case TYPE_OE_1_1:
        handleOE_1_1(
                  curA, curB, edgeOnA, leftHullA, leftHullB, sweepingDownward);
        break;

      case TYPE_OE_1_2:
        handleOE_1_2(
                  curA, curB, edgeOnA, leftHullA, leftHullB, sweepingDownward);
        break;

      case TYPE_OE_2_1:
        handleOE_2_1(
                  curA, curB, edgeOnA, leftHullA, leftHullB, sweepingDownward);
        break;

      case TYPE_OE_3_1:
        handleOE_3_1(
                  curA, curB, edgeOnA, leftHullA, leftHullB, sweepingDownward);
        break;

      case TYPE_OE_3_2:
        handleOE_3_2(
                  curA, curB, edgeOnA, leftHullA, leftHullB, sweepingDownward);
        break;

      case TYPE_OE_4_1:
        handleOE_4_1(
                  curA, curB, edgeOnA, leftHullA, leftHullB, sweepingDownward);
        break;

      case TYPE_OE_5_1:
        handleOE_5_1(
                  curA, curB, edgeOnA, leftHullA, leftHullB, sweepingDownward);
        break;

      case TYPE_OE_5_2:
        handleOE_5_2(
                  curA, curB, edgeOnA, leftHullA, leftHullB, sweepingDownward);
        break;

      case TYPE_OE_6_1:
        handleOE_6_1(
                  curA, curB, edgeOnA, leftHullA, leftHullB, sweepingDownward);
        break;

      case TYPE_OE_6_2:
      default:
        handleOE_6_2(
                  curA, curB, edgeOnA, leftHullA, leftHullB, sweepingDownward);
        break;

    }
}            


/** @brief handle a case where the edges incident to the current vertices of
 *         A and B cross each other.
 *         It marks the intersection information along the chains, and
 *         advance the current vertices.
 *         This is a dispatcher function.
 *
 *  @param curA (in/out): current vertex of the hull of A
 *
 *  @param curB (in/out): current vertex of the hull of B
 *
 *  @param leftHullA (in): true if the hull of A is the left hull.
 *
 *  @param leftHullB (in): true if the hull of B is the left hull
 *
 *  @param sweepingDownward (in): true if the sweeping direction is downward
 */
void IntersectionFinderConvexPolygon2D::handleCrossing(
    CHNode** curA,
    CHNode** curB,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    CHEdge* ea;
    CHEdge* eb;

    if ( ( leftHullA &&  sweepingDownward) ||
         (!leftHullA && !sweepingDownward)   ) {
        ea = (*curA)->mNext;
    }
    else {
        ea = (*curA)->mPrev;
    }

    if ( ( leftHullB &&  sweepingDownward) ||
         (!leftHullB && !sweepingDownward)   ) {
        eb = (*curB)->mNext;
    }
    else {
        eb = (*curB)->mPrev;
    }

    const Vec2& p1 = ea->mSrc->mP;
    Vec2        v1 = ea->mDst->mP - ea->mSrc->mP;
    const Vec2& p2 = eb->mSrc->mP;
    Vec2        v2 = eb->mDst->mP - eb->mSrc->mP;

    auto cls    = classifyCrossing(ea, eb);
    auto intSec = findIntersectionOfTwoLines(p1, v1, p2, v2);
    auto* newA  = splitEdge(
                      ea,
                      intSec,
                      true,
                      CHNode::NT_INTERSECTION
                  );

    auto* newB  = splitEdge(eb, newA, CHNode::NT_INTERSECTION);

    switch (cls) {

      case TYPE_IS_3:
        newB->mNext->mSrcInterior = true;
        newA->mPrev->mDstInterior = true;
        break;

      case TYPE_IS_4:
      default:
        newA->mNext->mSrcInterior = true;
        newB->mPrev->mDstInterior = true;
        break;
    }

    (*curA) = newA;
    (*curB) = newB;
    updateCurACurBCV(curA, curB, leftHullA, leftHullB, sweepingDownward);

}


/** @brief classifies the arrangement of the edges around the current nodes
 *         that are coincident.
 *
 *  @param  ca (in): current node on A
 *
 *  @param  ca (in): current node on B
 *
 *  @returns a classifier
 */
enum _classifier IntersectionFinderConvexPolygon2D::classifyCoincidentNodes(
    CHNode* ca,
    CHNode* cb
) {

    auto* pa = ca->mPrev;
    auto* pb = cb->mPrev;
    auto* na = ca->mNext;
    auto* nb = cb->mNext;

    auto col_papb = pa->isColinear(pb);
    auto col_panb = pa->isColinear(nb);
    auto col_pbna = pb->isColinear(na);
    auto col_nanb = na->isColinear(nb);
    auto dir_papb = pa->isInTheSameDirection(pb);
    auto dir_panb = pa->isInOppositeDirection(nb);
    auto dir_pbna = pb->isInOppositeDirection(na);
    auto dir_nanb = na->isInTheSameDirection(nb);

    if (col_papb && dir_papb) {
        if (col_nanb && dir_nanb) {
            return TYPE_CV_1_1;
        }    
        else if (areInCCW(pb, false, nb, true, na, true)) {
            return TYPE_CV_1_2;
        }
        else {
            return TYPE_CV_1_3;
        }
    }
    else if (col_panb && dir_panb) {
        if (col_pbna && dir_pbna) {
            return TYPE_CV_2_1;
        }
        else {
            return TYPE_CV_2_2;
        }
    }
    else if (col_pbna && dir_pbna) {
        return TYPE_CV_3_1;
    }
    else if (col_nanb && dir_nanb) {
        if (areInCCW(na, true, pb, false, pa, false)) {
            return TYPE_CV_4_1;
        }
        else {
            return TYPE_CV_4_2;
        }
    }
    else if (areInCCW(pa, false, nb, true,  pb, false) && 
             areInCCW(nb, true,  pb, false, na, true ) &&
             areInCCW(pb, false, na, true,  pa, false)    ) {
        return TYPE_CV_5;
    }
    else if (areInCCW(pa, false, nb, true,  na, true ) && 
             areInCCW(nb, true,  na, true,  pb, false) && 
             areInCCW(na, true,  pb, false, pa, false)    ) {
        return TYPE_CV_7;
    }
    else if (areInCCW(na, true,  nb, true,  pa, false) && 
             areInCCW(nb, true,  pa, false, pb, false) && 
             areInCCW(pa, false, pb, false, na, true )    ) {
        return TYPE_CV_8;
    }
    else if (areInCCW(nb, true,  na, true,  pa, false) && 
             areInCCW(na, true,  pa, false, pb, false) && 
             areInCCW(pa, false, pb, false, nb, true )    ) {
        return TYPE_CV_9;
    }
    else{
        return TYPE_CV_10;
    }
}


/** @brief classifies the arrangement of the incident edges of the current 
 *         node that are on the edge of the other hull.
 *
 *  @param  n1 (in): current node on A or B that are on e2.
 *
 *  @param  e2 (in): edge on A or B. n1 is on it.
 *
 *  @returns a classifier
 */
enum _classifier IntersectionFinderConvexPolygon2D::classifyVertexOnEdge(
    CHNode* n1,
    CHEdge* e2
) {
    auto* ep1 = n1->mPrev;
    auto* en1 = n1->mNext;

    auto col_ep1e2 = ep1->isColinear(e2);
    auto col_en1e2 = en1->isColinear(e2);
    bool dir_ep1e2 = ep1->isInTheSameDirection(e2);
    bool dir_en1e2 = en1->isInTheSameDirection(e2);

    if (col_ep1e2) {
        if (col_en1e2) {
            if (dir_ep1e2) { // dir_en1e1==true
                return TYPE_OE_3_1;
            }
            else {
                return TYPE_OE_1_1;
            }
        }
        else{ 
            if (dir_ep1e2) {
                return TYPE_OE_3_2;
            }
            else {
                return TYPE_OE_1_2;
            }
        }
    }
    else {
        if (col_en1e2) {
            if (dir_en1e2) {
                return TYPE_OE_4_1;
            }
            else {
                return TYPE_OE_2_1;
            }
        }
        else{ 
            if (areInCCW(en1,  true, e2,  true, ep1, false)) {
                return TYPE_OE_5_1;
            }
            else if (areInCCW(en1,  true, e2, false, ep1, false)) {
                return TYPE_OE_6_1;
            }
            else {
                auto ve2Perp = (e2->mDst->mP  - e2->mSrc->mP ).perp();
                auto ven1    = (en1->mDst->mP - en1->mSrc->mP);
                if (ve2Perp.dot(ven1) > 0.0) {
                    return TYPE_OE_5_2;
                } 
                else {
                    return TYPE_OE_6_2;
                }
            }
        }
    }
}


/** @brief classifies the arrangement of the two crossing edges
 *
 *  @param  ea (in): edge on A
 *
 *  @param  eb (in): edge on B
 *
 *  @returns a classifier
 */
enum _classifier IntersectionFinderConvexPolygon2D::classifyCrossing(CHEdge* ea, CHEdge* eb)
{
    if (areInCCW(ea, false, eb, false, ea,  true)) {
        return TYPE_IS_3;
    }
    else {
        return TYPE_IS_4;
    }
}


/** @brief tests the circular arrangement of three incident edges.
 *
 *  @param e1       (in):  edge 1
 *
 *  @param dir1     (in):  the src of edge1 is incident.
 *
 *  @param e2       (in):  edge 2
 *
 *  @param dir2     (in):  the src of edge2 is incident.
 *
 *  @param e3       (in):  edge 3
 *
 *  @param dir3     (in):  the src of edge3 is incident.
 *
 *  @return  true if e1 , e2, e3 are ordered in CCW.
 */
bool IntersectionFinderConvexPolygon2D::areInCCW(
    CHEdge* e1, 
    bool    dir1,
    CHEdge* e2, 
    bool    dir2,
    CHEdge* e3,
    bool    dir3
) {
    Vec2 v1 = e1->mDst->mP - e1->mSrc->mP;
    Vec2 v2 = e2->mDst->mP - e2->mSrc->mP;
    Vec2 v3 = e3->mDst->mP - e3->mSrc->mP;

    if (!dir1) v1.scale(-1.0);
    if (!dir2) v2.scale(-1.0);
    if (!dir3) v3.scale(-1.0);

    Vec2 v1p = v1.perp();
    Vec2 v2p = v2.perp();

    // If any two are parallel, then consider it true.
    // This is needed for polygon-edge and edge-polygon cases.
    if ((fabs(v1p.dot(v2)) <= mEpsilonAngle && v1.dot(v2)>0.0) ||
        (fabs(v2p.dot(v3)) <= mEpsilonAngle && v2.dot(v3)>0.0) ||
        (fabs(v1p.dot(v3)) <= mEpsilonAngle && v1.dot(v3)>0.0)   ) {
        return true;
    }

    long q3, q2;

    if (v1.dot(v3) >= 0.0) {
        if (v1p.dot(v3) >= 0.0) { q3 = 1; } // 1st quadrant
        else                    { q3 = 4; } // 4th quadrant
    }
    else {
        if (v1p.dot(v3) >= 0.0) { q3 = 2; } // 2nd quadrant
        else                    { q3 = 3; } // 3nd quadrant
    }

    if (v1.dot(v2) >= 0.0) {
        if (v1p.dot(v2) >= 0.0) { q2 = 1; } // 1st quadrant
        else                    { q2 = 4; } // 4th quadrant
    }
    else {
        if (v1p.dot(v2) >= 0.0) { q2 = 2; } // 2nd quadrant
        else                    { q2 = 3; } // 3nd quadrant
    }

    if (q3 == 1) {
        return v1.dot(v2p) < 0.0 && v3.dot(v2p) > 0.0;
    }
    else if (q3 == 2) {
        return (q2 == 1) || (q2 == 2 && v3.dot(v2p) > 0.0 );
    }
    else if (q3 == 3) {
        return (q2 == 1) || (q2 == 2) || (q2 == 3 && v3.dot(v2p) > 0.0 );   
    }
    else {
        return (q2 == 1) || (q2 == 2) || (q2 == 3) || 
                                               (q2 == 4 && v3.dot(v2p) > 0.0 );
    }
}


/** @brief returns the intersection point of two lines.
 *
 *  @param p1 (in): a point on line 1.
 *
 *  @param v1 (in): direction vector of line 1.
 *
 *  @param p2 (in): a point on line 2.
 *
 *  @param v2 (in): direction vector of line 2.
 *
 *  @return the intersection point
 */
Vec2 IntersectionFinderConvexPolygon2D::findIntersectionOfTwoLines(
    const Vec2& p1,
    const Vec2& v1,
    const Vec2& p2,
    const Vec2& v2
) {
    Vec2 v12 = p2 - p1;
    Vec2 v2perp = v2.perp();
    double s = v12.dot(v2perp) / v1.dot(v2perp);
    return p1 + v1 * s;
}

            
/** @brief splits the CHEdge into two CHEdges and one CHNode in the middle.
 *          
 *  @param e (in): edge to be split. This must be deleted after the call.
 *
 *  @param p (in): the coordinates of the middle point.
 *
 *  @return CHNode newly created.
 */
CHNode* IntersectionFinderConvexPolygon2D::splitEdge(
    CHEdge*                 e,
    const Vec2&             p,
    const bool              thisIsA,
    const enum CHNode::type t
)
{
    auto* src     = e->mSrc;
    auto* dst     = e->mDst;
    long index1 = std::min(src->mIndex, dst->mIndex);
    long index2 = std::max(src->mIndex, dst->mIndex);

    if (src->mIndexAux!=-1) {
        index1 = std::min(index1, src->mIndexAux);
        index2 = std::max(index2, src->mIndexAux);
    }
    if (dst->mIndexAux!=-1) {
        index1 = std::min(index1, dst->mIndexAux);
        index2 = std::max(index2, dst->mIndexAux);
    }

    CHNode* newN  = new CHNode(p, thisIsA, index1, index2, t, mEpsilonSquared);
    CHEdge* newEs = new CHEdge(mEpsilonAngle); 
    CHEdge* newEd = new CHEdge(mEpsilonAngle); 

    src->mNext    = newEs;
    newEs->mSrc   = src;
    newEs->mDst   = newN;
    newN->mPrev   = newEs;
    newN->mNext   = newEd;
    newEd->mSrc   = newN;
    newEd->mDst   = dst;
    dst->mPrev    = newEd;

    if (e->mSrcInterior) {
        newEs->mSrcInterior = true;
    }
    if (e->mDstInterior) {
        newEd->mDstInterior = true;
    }

    delete(e);
    return newN;
}


/** @brief splits the CHEdge into two CHEdges and one CHNode in the middle.
 *          
 *  @param e (in): edge to be split. This must be deleted after the call.
 *
 *  @param p (in): CHNode on the other hull that are coincident to the new node
 *                 to be created. It will be the peer of the new node.
 *
 *  @return CHNode newly created.
 */
CHNode* IntersectionFinderConvexPolygon2D::splitEdge(
    CHEdge*                 e,
    CHNode*                 n,
    const enum CHNode::type t
)
{
    auto* newN = splitEdge(e, n->mP, !(n->mThisIsA), t);
    newN->mPeer = n;
    n->mPeer = newN;

    return newN;
}


void IntersectionFinderConvexPolygon2D::handleCV_1_1(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;

    if (leftHullA && leftHullB) {
        if (sweepingDownward) {
            handleColinearForwardAForwardB(curA, curB); 
        }
        else {
            handleColinearBackwardABackwardB(curA, curB);
        }
    }
    else if (!leftHullA && !leftHullB) {
        if (sweepingDownward) {
            handleColinearBackwardABackwardB(curA, curB);
        }
        else {
            handleColinearForwardAForwardB(curA, curB);
        }
    }
    else {
        updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
    }
}


void IntersectionFinderConvexPolygon2D::handleCV_1_2(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;

    if (leftHullA && leftHullB) {
        if (sweepingDownward) {
            (*curA)->mNext->mSrcInterior = true;
            updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
        }
        else {
            handleColinearBackwardABackwardB(curA, curB);
        }
    }
    else if (!leftHullA && !leftHullB) {
        if (sweepingDownward) {
            handleColinearBackwardABackwardB(curA, curB);
        }
        else {
            (*curA)->mNext->mSrcInterior = true;
            updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
        }
    }
    else {
        if ((*curA)->mTop||(*curA)->mBottom) {
            (*curA)->mNext->mSrcInterior = true;
        }
        updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
    }
}


void IntersectionFinderConvexPolygon2D::handleCV_1_3(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;

    if (leftHullA && leftHullB) {
        if (sweepingDownward) {
            (*curB)->mNext->mSrcInterior = true;
            updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
        }
        else {
            handleColinearBackwardABackwardB(curA, curB);
        }
    }
    else if (!leftHullA && !leftHullB) {
        if (sweepingDownward) {
            handleColinearBackwardABackwardB(curA, curB);
        }
        else {
            (*curB)->mNext->mSrcInterior = true;
            updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
        }
    }
    else {
        if ((*curB)->mTop||(*curB)->mBottom) {
            (*curB)->mNext->mSrcInterior = true;
        }
        updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
    }

}


void IntersectionFinderConvexPolygon2D::handleCV_2_1(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;

    if (leftHullA && !leftHullB) { 
        if (sweepingDownward) {
            handleColinearForwardABackwardB(curA, curB);
        }
        else {
            handleColinearBackwardAForwardB(curA, curB);
        }
    }
    else if (!leftHullA && leftHullB) {
        if (sweepingDownward) {
            handleColinearBackwardAForwardB(curA, curB);
        }
        else {
            handleColinearForwardABackwardB(curA, curB);
        }
    }
    else {
        updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
    }

}


void IntersectionFinderConvexPolygon2D::handleCV_2_2(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;

    if (leftHullA && !leftHullB) { 
        if (sweepingDownward) {
            updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
        }
        else {
            handleColinearBackwardAForwardB(curA, curB);
        }
    }
    else if (!leftHullA && leftHullB) { 
        if (sweepingDownward) {
            handleColinearBackwardAForwardB(curA, curB);
        }
        else {
            updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
        }     
    }
    else {
        updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
    }
}


void IntersectionFinderConvexPolygon2D::handleCV_3_1(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;

    if (!leftHullA && leftHullB) { 
        if (sweepingDownward) {
            updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
        }
        else {
            handleColinearForwardABackwardB(curA, curB);
        }
    }
    else if ( leftHullA && !leftHullB) { 
        if (sweepingDownward) {
            handleColinearForwardABackwardB(curA, curB);
        }
        else {
            updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
        }     
    }
    else {
        updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
    }
}


void IntersectionFinderConvexPolygon2D::handleCV_4_1(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mPrev->mDstInterior = true;

    if (leftHullA && leftHullB) { 
        if (sweepingDownward) {
            handleColinearForwardAForwardB(curA, curB);
        }
        else {
            updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
        }
    }
    else if (!leftHullA && !leftHullB) { 
        if (sweepingDownward) {
            updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
        }
        else {
            handleColinearForwardAForwardB(curA, curB);
        }
    }
    else {
        updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
    }

}


void IntersectionFinderConvexPolygon2D::handleCV_4_2(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curA)->mPrev->mDstInterior = true;

    if (leftHullA && leftHullB) { 
        if (sweepingDownward) {
            handleColinearForwardAForwardB(curA, curB);
        }
        else {
            updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
        }
    }
    else if (!leftHullA && !leftHullB) { 
        if (sweepingDownward) {
            updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
        }
        else {
            handleColinearForwardAForwardB(curA, curB);
        }
    }
    else {
        updateCurACurBCV(curA,curB,leftHullA,leftHullB,sweepingDownward);
    }

}


void IntersectionFinderConvexPolygon2D::handleCV_5(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curA)->mIsTouchingPoint = true;
    (*curB)->mIsTouchingPoint = true;

    updateCurACurBCV(curA, curB, leftHullA, leftHullB, sweepingDownward);

}


void IntersectionFinderConvexPolygon2D::handleCV_7(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curA)->mNext->mSrcInterior = true;
    (*curB)->mPrev->mDstInterior = true;
    updateCurACurBCV(curA, curB, leftHullA, leftHullB, sweepingDownward);
}


void IntersectionFinderConvexPolygon2D::handleCV_8(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mNext->mSrcInterior = true;
    (*curA)->mPrev->mDstInterior = true;
    updateCurACurBCV(curA, curB, leftHullA, leftHullB, sweepingDownward);
}


void IntersectionFinderConvexPolygon2D::handleCV_9(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curA)->mNext->mSrcInterior = true;
    (*curA)->mPrev->mDstInterior = true;
    updateCurACurBCV(curA, curB, leftHullA, leftHullB, sweepingDownward);
}


void IntersectionFinderConvexPolygon2D::handleCV_10(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {
    (*curA)->mPeer = (*curB);
    (*curB)->mPeer = (*curA);
    (*curA)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mType = CHNode::NT_COINCIDENT_VERTEX;
    (*curB)->mNext->mSrcInterior = true;
    (*curB)->mPrev->mDstInterior = true;
    updateCurACurBCV(curA, curB, leftHullA, leftHullB, sweepingDownward);
}


/** @brief updates the curA and curB if they are coincident, and split
 *         the longer edge if necessary.
 *         A's next and B's next are in the sweeping direction.
 *
 *  @param curA  (in/out): current and updated vertex on A
 *
 *  @param curB  (in/out): current and updated vertex on B
 */
void IntersectionFinderConvexPolygon2D::handleColinearForwardAForwardB(
    CHNode** curA,
    CHNode** curB
) {
    auto*  nextA = (*curA)->mNext->mDst;
    auto*  nextB = (*curB)->mNext->mDst;
    Vec2   va    = nextA->mP - (*curA)->mP;
    Vec2   vb    = nextB->mP - (*curB)->mP;
    double na    = va.squaredNorm2();
    double nb    = vb.squaredNorm2();
    if (fabs(na - nb) <= mEpsilonSquared) {
        // cA+-------->+nA
        // cB+-------->+nB
        (*curA) = nextA;
        (*curB) = nextB;
    }
    else if (na > nb) {
        //split curA->mNext at nextB
        //
        //            split point
        // cA+---------*---->+nA
        // cB+-------->+nB
        (*curA) = splitEdge((*curA)->mNext, nextB, CHNode::NT_SPLIT_POINT);
        (*curB) = nextB;
    }
    else {
        //split curB->mNext at nextA
        //
        // cA+-------->+nA
        // cB+---------*---->+nB
        //            split point
        (*curA) = nextA;
        (*curB) = splitEdge((*curB)->mNext, nextA, CHNode::NT_SPLIT_POINT);
    }
}


/** @brief updates the curA and curB if they are coincident, and split
 *         the longer edge if necessary.
 *         A's prev and B's prev are in the sweeping direction.
 *
 *  @param curA  (in/out): current and updated vertex on A
 *
 *  @param curB  (in/out): current and updated vertex on B
 */
void IntersectionFinderConvexPolygon2D::handleColinearBackwardABackwardB(
    CHNode** curA,
    CHNode** curB
) {
    auto* prevB = (*curB)->mPrev->mSrc;
    auto* prevA = (*curA)->mPrev->mSrc;
    Vec2 va = prevA->mP - (*curA)->mP;
    Vec2 vb = prevB->mP - (*curB)->mP;
    double na = va.squaredNorm2();
    double nb = vb.squaredNorm2();
    if (fabs(na - nb) <= mEpsilonSquared) {
        // cA+<--------+pA
        // cB+<--------+pB
        (*curA) = prevA;
        (*curB) = prevB;
    }
    else if (na > nb) {
        //split curA->mPrev at prevB
        //
        //            split point
        // cA+<--------*-----+pA
        // cB+<--------+pB
        (*curA) = splitEdge((*curA)->mPrev, prevB, CHNode::NT_SPLIT_POINT);
        (*curB) = prevB;
    }
    else {
        //split curB->mPrev at prevA
        //
        // cA+<--------+pA
        // cB+<--------*-----+pB
        //            split point
        (*curA) = prevA;
        (*curB) = splitEdge((*curB)->mPrev, prevA, CHNode::NT_SPLIT_POINT);
    }
}


/** @brief updates the curA and curB if they are coincident, and split
 *         the longer edge if necessary.
 *         A's next and B's prev are in the sweeping direction.
 *
 *  @param curA  (in/out): current and updated vertex on A
 *
 *  @param curB  (in/out): current and updated vertex on B
 */
void IntersectionFinderConvexPolygon2D::handleColinearForwardABackwardB(
    CHNode** curA,
    CHNode** curB
) {
    auto* nextA = (*curA)->mNext->mDst;
    auto* prevB = (*curB)->mPrev->mSrc;
    Vec2 va = nextA->mP - (*curA)->mP;
    Vec2 vb = prevB->mP - (*curB)->mP;
    double na = va.squaredNorm2();
    double nb = vb.squaredNorm2();
    if (fabs(na - nb) <= mEpsilonSquared) {
        // cA+-------->+nA
        // cB+<--------+pB
        (*curA) = nextA;
        (*curB) = prevB;
    }
    else if (na > nb) {
        //split curA->mNext at prevB
        //
        //            split point
        // cA+---------*---->+nA
        // cB+<--------+pB
        (*curA) = splitEdge((*curA)->mNext, prevB, CHNode::NT_SPLIT_POINT);
        (*curB) = prevB;
    }
    else {
        //split curB->mPrev at nextA
        //
        // cA+-------->+nA
        // cB+<--------*-----+pB
        //            split point
        (*curA) = nextA;
        (*curB) = splitEdge((*curB)->mPrev, nextA, CHNode::NT_SPLIT_POINT);
    }
}


/** @brief updates the curA and curB if they are coincident, and split
 *         the longer edge if necessary.
 *         A's prev and B's next are in the sweeping direction.
 *
 *  @param curA  (in/out): current and updated vertex on A
 *
 *  @param curB  (in/out): current and updated vertex on B
 */
void IntersectionFinderConvexPolygon2D::handleColinearBackwardAForwardB(
    CHNode** curA,
    CHNode** curB
) {
    auto* prevA = (*curA)->mPrev->mSrc;
    auto* nextB = (*curB)->mNext->mDst;
    Vec2 va = prevA->mP - (*curA)->mP;
    Vec2 vb = nextB->mP - (*curB)->mP;
    double na = va.squaredNorm2();
    double nb = vb.squaredNorm2();
    if (fabs(na - nb) <= mEpsilonSquared) {
        // cA+<--------+pA
        // cB+-------->+nB
        (*curA) = prevA;
        (*curB) = nextB;
    }
    else if (na > nb) {
        //split curA->mPrev at nextB
        //
        //            split point
        // cA+<--------*-----+pA
        // cB+-------->+nB
        (*curA) = splitEdge((*curA)->mPrev, nextB, CHNode::NT_SPLIT_POINT);
        (*curB) = nextB;
    }
    else {
        //split curB->mNext at prevA
        //
        // cA+<--------+nA
        // cB+---------*---->+pB
        //            split point
        (*curA) = prevA;
        (*curB) = splitEdge((*curB)->mNext, prevA, CHNode::NT_SPLIT_POINT);
    }
}


void IntersectionFinderConvexPolygon2D::handleOE_1_1(
    CHNode** curA,
    CHNode** curB,
    bool     edgeOnA,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    if ( (leftHullA && !leftHullB) || (!leftHullA && leftHullB) ||
         (*curA)->mTop || (*curA)->mBottom  ||
         (*curB)->mTop || (*curB)->mBottom                           ) {

        handleOECommon(curA,curB,edgeOnA,leftHullA,leftHullB,sweepingDownward);
        updateCurACurBOE(curA, curB, leftHullA, leftHullB, sweepingDownward);
    }
}


void IntersectionFinderConvexPolygon2D::handleOE_1_2(
    CHNode** curA,
    CHNode** curB,
    bool     edgeOnA,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    if ( (leftHullA && !leftHullB) || (!leftHullA && leftHullB) ||
         (*curA)->mTop || (*curA)->mBottom  ||
         (*curB)->mTop || (*curB)->mBottom                           ) {

        handleOECommon(curA,curB,edgeOnA,leftHullA,leftHullB,sweepingDownward);
        updateCurACurBOE(curA, curB, leftHullA, leftHullB, sweepingDownward);
    }
}


void IntersectionFinderConvexPolygon2D::handleOE_2_1(
    CHNode** curA,
    CHNode** curB,
    bool     edgeOnA,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    if ( (leftHullA && !leftHullB) || (!leftHullA && leftHullB) ||
         (*curA)->mTop || (*curA)->mBottom  ||
         (*curB)->mTop || (*curB)->mBottom                           ) {

        handleOECommon(curA,curB,edgeOnA,leftHullA,leftHullB,sweepingDownward);
        updateCurACurBOE(curA, curB, leftHullA, leftHullB, sweepingDownward);
    }
}


void IntersectionFinderConvexPolygon2D::handleOE_3_1(
    CHNode** curA,
    CHNode** curB,
    bool     edgeOnA,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    if ( (leftHullA && leftHullB) || (!leftHullA && !leftHullB) ||
         (*curA)->mTop || (*curA)->mBottom  ||
         (*curB)->mTop || (*curB)->mBottom                           ) {

        handleOECommon(curA,curB,edgeOnA,leftHullA,leftHullB,sweepingDownward);
        updateCurACurBOE(curA, curB, leftHullA, leftHullB, sweepingDownward);
    }
}


void IntersectionFinderConvexPolygon2D::handleOE_3_2(
    CHNode** curA,
    CHNode** curB,
    bool     edgeOnA,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    if ( (leftHullA && leftHullB) || (!leftHullA && !leftHullB) ||
         (*curA)->mTop || (*curA)->mBottom  ||
         (*curB)->mTop || (*curB)->mBottom                           ) {
        handleOECommon(curA,curB,edgeOnA,leftHullA,leftHullB,sweepingDownward);
        if (edgeOnA) {
            (*curB)->mNext->mSrcInterior = true;
        }
        else {
            (*curA)->mNext->mSrcInterior = true;
        }
        updateCurACurBOE(curA, curB, leftHullA, leftHullB, sweepingDownward);
    }
}


void IntersectionFinderConvexPolygon2D::handleOE_4_1(
    CHNode** curA,
    CHNode** curB,
    bool     edgeOnA,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    if ( (leftHullA && leftHullB) || (!leftHullA && !leftHullB) ||
         (*curA)->mTop || (*curA)->mBottom  ||
         (*curB)->mTop || (*curB)->mBottom                           ) {

        if (edgeOnA) {
            (*curB)->mPrev->mDstInterior = true;
        }
        else {
            (*curA)->mPrev->mDstInterior = true;
        }
        handleOECommon(curA,curB,edgeOnA,leftHullA,leftHullB,sweepingDownward);
        updateCurACurBOE(curA, curB, leftHullA, leftHullB, sweepingDownward);
    }
}


void IntersectionFinderConvexPolygon2D::handleOE_5_1(
    CHNode** curA,
    CHNode** curB,
    bool     edgeOnA,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    handleOECommon(curA,curB,edgeOnA,leftHullA,leftHullB,sweepingDownward);
    if (edgeOnA) {
        (*curA)->mNext->mSrcInterior = true;
        (*curB)->mPrev->mDstInterior = true;
    }
    else {
        (*curB)->mNext->mSrcInterior = true;
        (*curA)->mPrev->mDstInterior = true;
    }
    updateCurACurBOE(curA, curB, leftHullA, leftHullB, sweepingDownward);
}


void IntersectionFinderConvexPolygon2D::handleOE_5_2(
    CHNode** curA,
    CHNode** curB,
    bool     edgeOnA,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    if ( (leftHullA && leftHullB) || (!leftHullA && !leftHullB) ||
         (*curA)->mTop || (*curA)->mBottom  ||
         (*curB)->mTop || (*curB)->mBottom                           ) {
        handleOECommon(curA,curB,edgeOnA,leftHullA,leftHullB,sweepingDownward);
        if (edgeOnA) {
            (*curB)->mNext->mSrcInterior = true;
            (*curB)->mPrev->mDstInterior = true;
        }
        else {
            (*curA)->mNext->mSrcInterior = true;
            (*curA)->mPrev->mDstInterior = true;
        }
        updateCurACurBOE(curA, curB, leftHullA, leftHullB, sweepingDownward);
    }
}


void IntersectionFinderConvexPolygon2D::handleOE_6_1(
    CHNode** curA,
    CHNode** curB,
    bool     edgeOnA,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    handleOECommon(curA,curB,edgeOnA,leftHullA,leftHullB,sweepingDownward);
    if (edgeOnA) {
        (*curB)->mNext->mSrcInterior = true;
        (*curA)->mPrev->mDstInterior = true;
    }
    else {
        (*curA)->mNext->mSrcInterior = true;
        (*curB)->mPrev->mDstInterior = true;
    }
    updateCurACurBOE(curA, curB, leftHullA, leftHullB, sweepingDownward);
}


void IntersectionFinderConvexPolygon2D::handleOE_6_2(
    CHNode** curA,
    CHNode** curB,
    bool     edgeOnA,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    if ( (leftHullA && !leftHullB) || (!leftHullA && leftHullB) ||
         (*curA)->mTop || (*curA)->mBottom  ||
         (*curB)->mTop || (*curB)->mBottom                           ) {
        handleOECommon(curA,curB,edgeOnA,leftHullA,leftHullB,sweepingDownward);
        (*curA)->mIsTouchingPoint = true;
        (*curB)->mIsTouchingPoint = true;
        updateCurACurBOE(curA, curB, leftHullA, leftHullB, sweepingDownward);
    }
}


/** @brief split the edge at the point where the other vertec coincides, and
 *         advance the current vertex on that edge.
 *
 *        if edgeOnA == true
 *
 *                \
 *                 \curB
 *         +--------*---------+
 *       curA      split    nextA/prevA
 *                 point
 *                (new curA)
 *
 *         -------------------->
 *           sweeping direction
 *
 *  @param curA (in/out): current and updated vertex of the hull of A
 *
 *  @param curB (in/out): current and updated vertex of the hull of B
 *
 *  @param edgeOnA (in): true if the edge is on hull A
 *
 *  @param leftHullA (in): true if the hull of A is the left hull.
 *
 *  @param leftHullB (in): true if the hull of B is the left hull
 *
 *  @param sweepingDownward (in): true if the sweeping direction is downward
 */
void IntersectionFinderConvexPolygon2D::handleOECommon(
    CHNode** curA,
    CHNode** curB,
    bool     edgeOnA,
    bool     leftHullA,
    bool     leftHullB, 
    bool     sweepingDownward
) {
    CHEdge* e;
    if (edgeOnA) {
        if ( ( leftHullA &&  sweepingDownward) ||
             (!leftHullA && !sweepingDownward)   ) {
             e = (*curA)->mNext;
        }
        else {
             e = (*curA)->mPrev;
        }
        (*curA) = splitEdge(e, *curB, CHNode::NT_SPLIT_POINT);
    }
    else {
       if ( ( leftHullB &&  sweepingDownward) ||
            (!leftHullB && !sweepingDownward)   ) {
             e = (*curB)->mNext;
        }
        else {
             e = (*curB)->mPrev;
        }
        (*curB) = splitEdge(e, *curA, CHNode::NT_SPLIT_POINT);
    }
}


}// namespace Makena


#ifndef _MAKENA_INTERSECTION_CONVEX_POLYGON_2D_HPP_
#define _MAKENA_INTERSECTION_CONVEX_POLYGON_2D_HPP_

#include<vector>
#include<iostream>
#include "primitives.hpp"
#include "loggable.hpp"

/**
 * @file intersection_convex_polygon_2d.hpp
 *
 * @brief Finds the intersection of two convex polygons in 2D by 4 iteration 
 *        of sweeping.
 */

namespace Makena {

using namespace std;

/** @brief find the intersection of two convex polygons.
 *
 *   This is a sweep algorithm. It first splits a convex polygon
 *   into two half hulls at the two extremal points along the y-axis.
 *   Let A and B two convex polygons.  The vertices of
 *   each convex polygon are arranged in CCW.
 *   A and B are further split into AL, AR, BL, and BR.
 *   The vertices in AL and BL go down from the top Y to the bottom Y in CCW 
 *   ordering,  and AR and BR go up from the bottom Y to the top Y.
 *
 *         Top Y
 *           __                            __
 *          /  \                  /          \
 *         /    \                /            \
 *        |      |   ==>        |              |
 *        |      |              |              |
 *         \    /                \            /
 *          \__/                  \__        /  
 *        Bottom Y                 AL        AR
 * 
 *   There are 4 sweeps: AL:BL, AL:BR, AR:BL, and AR:BR.
 *   In AL:BL and AR:BR, there can be multiple intersections points.
 *   In AL:BR and AR:BL, there can be at most two intersection points due
 *   to their convexity.
 *
 *              AL
 *             /
 *          --*--BL
 *         | /
 *         |/
 *         *
 *        /|
 *       | |
 *       | |
 *       | |
 *        \|
 *         *
 *         |\
 *         | \
 *         |  \
 *          \  \
 *            \ \
 *              \\               
 *                *
 *                 \\
 *                  \ \___
 *                   \ 
 *
 *   Each sweep detects the intersections and classify the vicinity into
 *   exterior/interior.
 *
 *   If there is no intersection, then either two polygons are disjoint, 
 *   or one is included in the other.
 *   This is detected by picking a representative vertex from a convex polygon
 *   and test its orientation around the directed edges along the other convex
 *   polygon.
 *
 *   If there is an intersection point, then explore the hulls along the
 *   interior side until it comes back to the original intersection point.
 *
 *
 *   ---------------------
 *   Type of intersections
 *   ---------------------
 *
 *   1. Coincident vertices  
 *
 *           o a1/b1
 *          / \
 *         /   \
 *        o     o
 *       a2     b2
 *
 *   2. One vertex on the interior of an edge.
 *  
 *              o b1
 *             /
 *        o---o---o
 *        a1  b2  a2
 *
 *   3. An edge intersects the other at their interior points.
 *
 *       a1     b1
 *        o     o
 *         \   /
 *          \ /
 *           * <- intersection
 *          / \
 *         /   \
 *        o     o
 *       b2     a2
 *
 *
 *   ------------------
 *   Ordering of vertex
 *   ------------------
 *
 *         Z   
 *         ^  o   o   
 *         |  p1  p2
 *         |
 *      ---+-----> Y
 *         |
 * 
 *    p1 == p2 <==>  p1.z==p2.z && p1.y == p2.y
 *
 *    p1 < p2  <==>  p1.z < p2.z || (p1.z==p2.z && p1.y > p2.y)
 *
 *  ----------------------------------------
 *  Classification onhf coincident vertices
 *  ----------------------------------------
 *  
 *  Legend
 *   O : touching point (intersection)
 *   = : touching edge (intersection)
 *   # : interior of the intersection
 * 
 *
 *  1.pa & pb are colinear
 *
 *    1.1. na & nb are colinear     TYPE_CV_1_1
 *         
 *                   na=nb
 *             #### /
 *             ### /
 *        papb ---o
 *
 *    1.2. pab -> nb -> na          TYPE_CV_1_2
 *
 *             ## na nb
 *             ## | /
 *             ## |/
 *        papb ---o
 *
 *    1.3. pab -> na -> nb          TYPE_CV_1_3
 *
 *             ## nb na
 *             ## | /
 *             ## |/
 *        papb ---o
 *
 *  2. pa & nb are colinear
 *    
 *     2.1. pb & na are colinear    TYPE_CV_2_1
 *          (they must be all colinear)
 *        panb ===o=== pbna
 *
 *     2.2. panb -> pb -> na        TYPE_CV_2_2
 *
 *                    na
 *                   /
 *                  /
 *         panb ===o
 *                  \
 *                   \
 *                    pb
 *
 *  3. pb & na are colinear
 *    
 *     3.1. pbna -> pa -> nb       TYPE_CV_3_1
 *
 *                    nb
 *                   /
 *                  /
 *         pbna ===o
 *                  \
 *                   \
 *                    pa
 *
 *
 *   4. na & nb are colinear
 *
 *    4.1. nab -> pb -> pa         TYPE_CV_4_1
 *
 *        nanb ---o
 *             ## |\
 *             ## | \
 *             ## pb pa
 *
 *    4.2. nab -> pa -> pb         TYPE_CV_4_2
 *
 *       nanb  ---o
 *             ## |\
 *             ## | \
 *             ## pa pb
 *
 *   5. na -> pa -> nb -> pb       TYPE_CV_5
 *
 *        pa   na
 *         \   /
 *          \ /
 *           O  Touching point
 *          / \
 *         /   \
 *        nb   pb
 *
 *   7.  pa -> nb -> na -> pb      TYPE_CV_7
 *
 *        pa   pb
 *         \   / ##
 *          \ / ###
 *           o ####     Interior: pb->c->na
 *          / \ ###
 *         /   \ ##
 *        nb   na
 *
 *   8. na -> nb -> pa -> pb       TYPE_CV_8
 *
 *        pb   pa
 *         \   / ##
 *          \ / ###
 *           o ####     Interior: pa->c->nb
 *          / \ ###
 *         /   \ ##
 *        na   nb
 *
 *   9. nb -> na -> pa -> pb       TYPE_CV_9
 *
 *          pb  pa
 *           | / ##
 *           |/ ###
 *           o ####     Interior: pa->c->na
 *           |\ ###
 *           | \ ##
 *          nb  na
 *
 *   10. na -> nb -> pb -> pa      TYPE_CV_10
 *
 *          pa  pb
 *           | / ##
 *           |/ ###
 *           o ####     Interior: pa->c->na
 *           |\ ###
 *           | \ ##
 *          na  nb
 *
 *  ---------------------------------------------------------------------
 *  Classification when a vertex cb is on the interior of an edge (pa,na)
 *  ---------------------------------------------------------------------
 *   
 *  1. na=pb
 *
 *  1.1. nb=pa                     TYPE_OE_1_1
 *
 *               cb
 *      panb ====+====na/pb
 *
 *  1.2. pa->nb->na=pb
 *
 *            cb
 *     pa ----+====na/pb           TYPE_OE_1_2
 *            |
 *            |
 *           nb
 *
 *  2. nb=pa                       TYPE_OE_2_1
 *               cb
 *      panb ====+---- na
 *               |
 *               |
 *               pb
 *
 *  3. pa=pb
 *
 *  3.1.  papb->nanb               TYPE_OE_3_1
 *
 *          #########
 *          #########
 *     papb ----+----nanb
 *
 *  3.2.  papb->na->nb             TYPE_OE_3_2
 *
 *          ###### nb
 *          ##### /
 *          #### /
 *     papb ----+----na
 *
 *  4. na=nb                       TYPE_OE_4_1
 *
 *           pb
 *            \ #####
 *             \ ####
 *       pa ----+----nanb
 *
 *
 *  5. na->pb->pa
 *
 *  5.1. na->pb->pa->nb            TYPE_OE_5_1
 *
 *               pb
 *              /##
 *             /###
 *     pa ----+---- na
 *            |
 *            |
 *            nb
 *
 *  5.2. na->nb->pb->pa            TYPE_OE_5_2
 *
 *        pb ### nb
 *          \###/
 *           \#/
 *     pa ----+---- na
 *
 *  6.1. pa->pb->na->nb            TYPE_OE_6_1
 *
 *        ### nb
 *        ### |
 *        ### |
 *     pa ----+---- na
 *           /
 *          /
 *        pb
 *
 *
 *  6.2. pa->nb->pb->na            TYPE_OE_6_2
 *
 *        Touching point
 *     pa ----O---- na
 *           /|
 *          / |
 *        nb  pb
 *
 *
 *  ---------------------------------------------------
 *  Classification when edge (pa,na) intersects (pb,nb)
 *  ---------------------------------------------------
 *
 *  3. pa -> pb -> na -> nb        TYPE_IS_3
 *
 *     pa ######
 *       \ #####
 *        \ ####
 *   pb----+----nb
 *          \
 *           \
 *            na
 *
 *  4. pa -> nb -> na -> pb        TYPE_IS_4
 *
 *     pa
 *       \
 *        \
 *   nb----+----pb
 *          \ ###
 *           \ ##
 *            na
 */


/** @brief describes the type of arrangement of a pair of two edges, each
 *         from the two given convex hulls.
 *         The details are given in intersection_convex_polygon_2d.hpp.
 */
enum _classifier {

    TYPE_CV_1_1,
    TYPE_CV_1_2,
    TYPE_CV_1_3,
    TYPE_CV_2_1,
    TYPE_CV_2_2,
    TYPE_CV_3_1,
    TYPE_CV_4_1,
    TYPE_CV_4_2,
    TYPE_CV_5,
    TYPE_CV_7,
    TYPE_CV_8,
    TYPE_CV_9,
    TYPE_CV_10,
    TYPE_OE_1_1,
    TYPE_OE_1_2,
    TYPE_OE_2_1,
    TYPE_OE_3_1,
    TYPE_OE_3_2,
    TYPE_OE_4_1,
    TYPE_OE_5_1,
    TYPE_OE_5_2,
    TYPE_OE_6_1,
    TYPE_OE_6_2,
    TYPE_IS_3,
    TYPE_IS_4

};

class CHEdge;

/** @brief represents the vertices of the given convex hulls, and also
 *         the intersection points that are added during the sweeps.
 */
class CHNode {

  public:

    enum type {

        NT_NONE,


        NT_ORIGINAL_VERTEX,   // Represents an original point given in one of 
                              // two polygons.

        NT_COINCIDENT_VERTEX, // An original vertex but it coincides with a
                              // vertex from the other polygon.

        NT_SPLIT_POINT,       // This resulted from a split of an edge, which
                              // is found to coincide with a vertex of the
                              // other polygon.
        NT_INTERSECTION,      // This resulted from a proper intersection of 
                              // two edges, each from one of two polygons.

        NT_END

    };

    CHEdge*   mPrev;            // Previous edge (CW)

    CHEdge*   mNext;            // Next edge (CCW)

    CHNode*   mPeer;            // Peer edge if this is an intersection of two 
                                // hulls.

    CHNode*   mExtraCoincidentLinkUp;
    CHNode*   mExtraCoincidentLinkDown;

    enum type mType;            // Type of this node. See above.


    bool      mIsTouchingPoint; // True if this is the only touching point of
                                // two convex hulls.

    bool      mTop;             // True if this is the vertex of the hull 
                                // that has the highest Y coordinate.

    bool      mBottom;          // True if this is the vertex of the hull
                                // that has the lowest Y coordinate.

    Vec2      mP;               // The coordinates of the vertex.

    bool      mThisIsA;         // True if this represents a point of Ain.

    long      mIndex;           // Index into into Ain or Bin given to 
                                // IntersectionFinderConvexPolygon2D::
                                // findIntersection().

    long      mIndexAux;        // Index into into Ain or Bin to identify
                                // the split edge with mIndex.


    const double mEpsilonSquared;
    const double mEpsilonLinear;
    CHNode(
        const Vec2& p,
        bool thisIsA,
        long i,
        long j,
        enum type t,
        const double epsilonSquared = EPSILON_SQUARED,
        const double epsilonLinear  = EPSILON_LINEAR
    ):
        mPrev(nullptr),
        mNext(nullptr),
        mPeer(nullptr),
        mExtraCoincidentLinkUp(nullptr),
        mExtraCoincidentLinkDown(nullptr),
        mType(t),
        mIsTouchingPoint(false),
        mTop(false),
        mBottom(false),
        mP(p),
        mThisIsA(thisIsA),
        mIndex(i),
        mIndexAux(j),
        mEpsilonSquared(epsilonSquared),
        mEpsilonLinear(epsilonLinear){;}

    ~CHNode(){;}

    bool isCoincident(CHNode* p) {
        auto v = this->mP - p->mP;
        return fabs(v.x())<= mEpsilonLinear && fabs(v.y())<= mEpsilonLinear;
    }

    bool isHigherThan(CHNode* p) {
        return  this->mP.y() > (p->mP.y() + mEpsilonLinear)        || 
                ( fabs(this->mP.y()- p->mP.y())<= mEpsilonLinear && 
                  (mEpsilonLinear + this->mP.x()) < p->mP.x()    ); 
    }

    inline bool isOnEdge(CHEdge* e);


    inline ostream&debugPrint(ostream&s);

};


/** @brief represents the edgesof the given convex hulls, and also
 *         the edges of the intersection that are made during the sweeps.
 */
class CHEdge {

  public:
    CHNode*       mSrc;          // Incident vertex CW side.

    bool          mSrcInterior;  // True if the part in the vicinity of 
                                 // the src node is part of the intersection.

    CHNode*       mDst;          // Incident vertex CCW side.

    bool          mDstInterior;  // True if the part in the vicinity of 
                                 // the dst node is part of the intersection.

    const double  mEpsilonAngle;

    CHEdge(const double epsilonAngle = EPSILON_ANGLE):
        mSrc(nullptr),
        mSrcInterior(false),
        mDst(nullptr),
        mDstInterior(false),
        mEpsilonAngle(epsilonAngle){;}

    ~CHEdge(){;}

    /** @brief test the orientation of this edge against another.
     *         True if two edges are facing the different direction in the
     *         sense that the dot product of those are negative.
     */
    bool isInOppositeDirection(CHEdge* e) {
        Vec2 v12 = mDst->mP    - mSrc->mP;
        Vec2 v34 = e->mDst->mP - e->mSrc->mP;
        return v12.dot(v34) < 0.0;
    }

    /** @brief test the orientation of this edge against another.
     *         True if two edges are facing the same direction in the
     *         sense that the dot product of those are positive.
     */
    bool isInTheSameDirection(CHEdge* e) {
        Vec2 v12 = mDst->mP    - mSrc->mP;
        Vec2 v34 = e->mDst->mP - e->mSrc->mP;
        return v12.dot(v34) > 0.0;
    }


    /** @brief test the orientation of this edge against another.
     *         True if two edges are in the same direction (parallel).
     */
    bool isColinear(CHEdge* e) {
        Vec2 v12     = mDst->mP - mSrc->mP;
        Vec2 v12perp = v12.perp();
        Vec2 v34     = e->mDst->mP - e->mSrc->mP;
        return fabs(v12perp.dot(v34)) < mEpsilonAngle;
    }


    /** @brief test the location of the given point on yz-plane relative to
     *         this edge. True if the point is on the left side if you stand
     *         at the src and look at dst.
     *
     *               _->p (CCW)
     *             _/
     *           _/
     *       src -----> dst
     */
    bool isCCWAroundSrcTo(const Vec2& p) {
        Vec2 v12     = mDst->mP - mSrc->mP;
        Vec2 v12perp = v12.perp();
        Vec2 v1p     = p - mSrc->mP;
        return v12perp.dot(v1p) > 0.0;
    }


    /** @brief test if this edge intersects with another in the interior.
     */
    bool doesIntersect(CHEdge* e) {
        const Vec2& v1 = e->mSrc->mP;
        const Vec2& v2 = e->mDst->mP;
        auto b1 = this->isCCWAroundSrcTo(v1);
        auto b2 = this->isCCWAroundSrcTo(v2);
        if ( (b1 && b2) || (!b1 && !b2) ) {
            return false;
        }
        const Vec2& v3 = this->mSrc->mP;
        const Vec2& v4 = this->mDst->mP;
        auto b3 = e->isCCWAroundSrcTo(v3);
        auto b4 = e->isCCWAroundSrcTo(v4);
        if ( (b3 && b4) || (!b3 && !b4) ) {
            return false;
        }
        return true;
    }
};


inline bool CHNode::isOnEdge(CHEdge* e) {
    const Vec2& p1 = e->mSrc->mP;
    const Vec2& p2 = e->mDst->mP;
    Vec2 v12 = p2 - p1;
    Vec2 v1p = mP - p1;
    Vec2 v12perp = v12.perp();
    if (fabs(v12perp.dot(v1p)) < mEpsilonSquared) {
        double s = v12.dot(v1p) / v12.squaredNorm2();
        if (-1.0 * mEpsilonSquared <= s && s <= 1.0 + mEpsilonSquared) {
            return true;
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}


inline ostream&CHNode::debugPrint(ostream&s)
 {
    s << "(" << mP.x() << "," << mP.y() << ") ";
    switch (mType) {
      case NT_NONE:              cerr << "[NONE] " ;              break;
      case NT_ORIGINAL_VERTEX:   cerr << "[ORIGINAL_VERTEX] " ;   break;
      case NT_COINCIDENT_VERTEX: cerr << "[COINCIDENT_VERTEX] " ; break;
      case NT_SPLIT_POINT:       cerr << "[SPLIT_POINT] " ;       break;
      case NT_INTERSECTION:      cerr << "[INTERSECTION] " ;      break;
      case NT_END:               cerr << "[END] " ;               break;
    };

    if (mIsTouchingPoint){ s << "TP ";}
    if (mPeer!=nullptr) { s << "HAS PEER ";}
    if (mNext!= nullptr) {
        if (mNext->mSrcInterior) { s << "NEXT INTERIOR "; }
    }
    if (mPrev!= nullptr) {
        if (mPrev->mDstInterior) { s << "PREV INTERIOR"; }
    }
    s << "\n";
    return s;
}


class IntersectionFinderConvexPolygon2D : public Loggable {

  public:

    enum VertexType {
        IT_NONE,
        IT_EDGE_VERTEX,
        IT_VERTEX_VERTEX,
        IT_VERTEX_EDGE,
        IT_EDGE_EDGE,
        IT_VERTEX_INTERIOR,
        IT_INTERIOR_VERTEX,
        IT_END
    };


    class InputElem {

      public:
        InputElem(const Vec2& p, long i):mP(p), mIndex(i){;}

        Vec2            mP;
        long            mIndex;

    };


    class OutputElem {

      public:

        inline OutputElem(
            const Vec2& p,
            long ia1,
            long ia2,
            long ib1,
            long ib2,
            enum VertexType t
        ):
            mP(p),
            mIndexA(ia1),
            mIndexAaux(ia2),
            mIndexB(ib1),
            mIndexBaux(ib2),
            mType(t){;}

        Vec2            mP;
        long            mIndexA;
        long            mIndexAaux;
        long            mIndexB;
        long            mIndexBaux;
        enum VertexType mType;
    };


  public:

    inline IntersectionFinderConvexPolygon2D(
        const double epsilonLinear,
        const double epsilonSquared,
        const double epsilonAngle,
        std::ostream& os = std::cerr
    );

    inline IntersectionFinderConvexPolygon2D(std::ostream& os = std::cerr);

    inline ~IntersectionFinderConvexPolygon2D();

    /** @brief finds the intersection of two convex polygons on yz-plane.
     *
     *  @param inA (in): Points of convex hull A arranged in CCW ordering.
     *  @param inB (in): Points of convex hull B arranged in CCW ordering.
     *
     *  @param out (out): Points of the vertices of the intersecion.
     *
     *  @return true if twho convex hulls intersect.
     */
    bool findIntersection(
        vector<InputElem>&  inA,
        vector<InputElem>&  inB,
        vector<OutputElem>& out
    );


#ifdef UNIT_TESTS
   public:
#else
   private:
#endif

    bool findIntersection_point_point(
        vector<InputElem*>& A,
        vector<InputElem*>& B,
        vector<OutputElem>& intsec
    );

    bool findIntersection_point_edge(
        vector<InputElem*>& A,
        vector<InputElem*>& B,
        vector<OutputElem>& intsec
    );

    bool findIntersection_point_polygon(
        vector<InputElem*>& A,
        vector<InputElem*>& B,
        vector<OutputElem>& intsec
    );

    bool findIntersection_edge_point(
        vector<InputElem*>& A,
        vector<InputElem*>& B,
        vector<OutputElem>& intsec
    );

    bool findIntersection_edge_edge(
        vector<InputElem*>& A,
        vector<InputElem*>& B,
        vector<OutputElem>& intsec
    );

    bool findIntersection_edge_polygon(
        vector<InputElem*>& A,
        vector<InputElem*>& B,
        vector<OutputElem>& intsec
    );

    bool findIntersection_polygon_point(
        vector<InputElem*>& A,
        vector<InputElem*>& B,
        vector<OutputElem>& intsec
    );

    bool findIntersection_polygon_edge(
        vector<InputElem*>& A,
        vector<InputElem*>& B,
        vector<OutputElem>& intsec
    );

    bool findIntersection_polygon_polygon(
        vector<InputElem*>& A,
        vector<InputElem*>& B,
        vector<OutputElem>& intsec
    );

    inline Vec2 findPerpendicularProjectionOfPointOntoEdge(
        const Vec2& PTest,
        const Vec2& p1,
        const Vec2& p2,
        double&     s
    );

    bool findIntersectionOfEdgeAndEdge(
        const Vec2&   a1,
        const Vec2&   a2,
        const Vec2&   b1,
        const Vec2&   b2,
        vector<Vec2>& intsec
    );
    
    bool findIntersectionOfPointAndPolygon(
        const Vec2&         P,
        const vector<Vec2>& pts
    );

    void removeDegeneracy(
        vector<InputElem>&  raw,
        vector<InputElem*>& cleaned
    );

    void removeDegeneracy(
        vector<OutputElem>& raw,
        vector<OutputElem>& cleaned
    );

    void constructCycle(
        vector<InputElem*>& pts,
        CHNode**            topV,
        CHNode**            bottomV,
        bool                thisIsA,
        long&               topIndex,
        long&               bottomIndex
    );

    void destructCycle(CHNode* n);

    void exploreInterior_polygon_polygon(
        CHNode*             topA,
        vector<OutputElem>& intsec,
        long                sizeA,
        long                sizeB
    );

    void exploreInterior_edge_polygon(
        CHNode*             topA,
        CHNode*             bottomA,
        vector<OutputElem>& intsec
    );

    void exploreInterior_polygon_edge(
        CHNode*             topB,
        CHNode*             bottomB,
        vector<OutputElem>& intsec
    );

    void exploreInterior_edge_edge(
        CHNode*             topA,
        CHNode*             bottomA,
        vector<OutputElem>& intsec
    );

    inline CHNode* visitVertexAndAdvanceForward(
        CHNode*     curV,
        bool&       finished, 
        bool&       touchingPoint,
        bool&       interior,
        bool&       intsecFound,
        bool&       needToVisitBackward
    );

    inline OutputElem constructOutputElem(CHNode* curV);

    void constructInnerOutput(
        vector<OutputElem>& intsec,
        vector<InputElem*>& points,
	const bool          innerIsA
    );

    bool sweepALBL(CHNode* topA, CHNode* topB, long sizeA, long sizeB);

    bool sweepARBR(CHNode* topA, CHNode* topB, long sizeA, long sizeB);

    bool sweepALBR(CHNode* topA, CHNode* topB, long sizeA, long sizeB);

    bool sweepARBL(CHNode* topA, CHNode* topB, long sizeA, long sizeB);

    void handleCoincidence(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCoincidenceAhasPeer(
        enum _classifier cls,
        CHNode**         curA,
        CHNode**         curB,
        bool             leftHullA,
        bool             leftHullB,
        bool             sweepingDownward
    );

    void handleCoincidenceBhasPeer(
        enum _classifier cls,
        CHNode**         curA,
        CHNode**         curB,
        bool             leftHullA,
        bool             leftHullB,
        bool             sweepingDownward
    );


    void handleVertexOnEdge(
        CHNode** curA,
        CHNode** curB,
        bool     edgeOnA,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCrossing(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    enum _classifier classifyCoincidentNodes(CHNode* ca, CHNode* cb);

    enum _classifier classifyVertexOnEdge   (CHNode* n1, CHEdge* e2);

    enum _classifier classifyCrossing       (CHEdge* ea, CHEdge* eb);

    bool areInCCW(
        CHEdge* e1,
        bool    dir1,
        CHEdge* e2,
        bool    dir2,
        CHEdge* e3,
        bool    dir3
    );

    Vec2 findIntersectionOfTwoLines(
        const Vec2& p1,
        const Vec2& v1,
        const Vec2& p2,
        const Vec2& v2
    );

    CHNode* splitEdge(
        CHEdge*                 e,
        const Vec2&             p,
        const bool              thisIsA,
        const enum CHNode::type t
    );

    CHNode* splitEdge(
        CHEdge*                 e,
        CHNode*                 n,
        const enum CHNode::type t
    );



    void handleCV_1_1(
        CHNode** curA,
        CHNode** curB, 
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCV_1_2(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCV_1_3(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCV_2_1(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCV_2_2(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCV_3_1(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCV_4_1(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCV_4_2(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCV_5(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCV_7(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCV_8(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCV_9(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleCV_10(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleColinearForwardAForwardB  (CHNode** curA, CHNode** curB);

    void handleColinearBackwardABackwardB(CHNode** curA, CHNode** curB);

    void handleColinearForwardABackwardB (CHNode** curA, CHNode** curB);

    void handleColinearBackwardAForwardB (CHNode** curA, CHNode** curB);

    void handleOE_1_1(
        CHNode** curA,
        CHNode** curB,
        bool     edgeOnA, 
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleOE_1_2(
        CHNode** curA,
        CHNode** curB,
        bool     edgeOnA,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleOE_2_1(
        CHNode** curA,
        CHNode** curB,
        bool     edgeOnA,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleOE_3_1(
        CHNode** curA,
        CHNode** curB,
        bool     edgeOnA,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleOE_3_2(
        CHNode** curA,
        CHNode** curB,
        bool     edgeOnA,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleOE_4_1(
        CHNode** curA,
        CHNode** curB,
        bool     edgeOnA,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleOE_5_1(
        CHNode** curA,
        CHNode** curB,
        bool     edgeOnA,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleOE_5_2(
        CHNode** curA,
        CHNode** curB,
        bool     edgeOnA,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleOE_6_1(
        CHNode** curA,
        CHNode** curB,
        bool     edgeOnA,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleOE_6_2(
        CHNode** curA,
        CHNode** curB,
        bool     edgeOnA,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    void handleOECommon(
        CHNode** curA,
        CHNode** curB,
        bool     edgeOnA,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    inline void updateCurACurBCV(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    inline void updateCurACurBOE(
        CHNode** curA,
        CHNode** curB,
        bool     leftHullA,
        bool     leftHullB,
        bool     sweepingDownward
    );

    inline CHNode* newVertexCandidate(
        CHNode* cur,
        bool    leftHull,
        bool    sweepingDownward
    );

    inline double clamp(double v, double lower, double upper);


    inline bool horizontalTest(CHEdge* ea, CHEdge* eb);

    inline bool boxTest(
        const Vec2& a1,
        const Vec2& a2,
        const Vec2& b1,
        const Vec2& b2
    );

    inline bool boxTest(
        const Vec2& a1,
        const Vec2& a2,
        const Vec2& b
    );

    void logChains(
        enum LogLevel lvl,
        const char*   _file,
        const int     _line,
        CHNode*       na,
        CHNode*       nb
    ) const;

    void logDumpInput(
        enum LogLevel lvl,
        const char*   _file,
        const int     _line
    ) const;

    vector<InputElem*> mAcleaned;
    vector<InputElem*> mBcleaned;

    const double  mEpsilonLinear;
    const double  mEpsilonSquared;
    const double  mEpsilonAngle;

};


IntersectionFinderConvexPolygon2D::IntersectionFinderConvexPolygon2D(
    const double epsilonLinear,
    const double epsilonSquared,
    const double epsilonAngle,
    std::ostream& os
)
    :Loggable(os),
     mEpsilonLinear(epsilonLinear),
     mEpsilonSquared(epsilonSquared),
     mEpsilonAngle(epsilonAngle){;}


IntersectionFinderConvexPolygon2D::IntersectionFinderConvexPolygon2D(
    std::ostream& os
)
    :Loggable(os),
     mEpsilonLinear(EPSILON_LINEAR),
     mEpsilonSquared(EPSILON_SQUARED),
     mEpsilonAngle(EPSILON_ANGLE){;}


IntersectionFinderConvexPolygon2D::~IntersectionFinderConvexPolygon2D(){;}


double IntersectionFinderConvexPolygon2D::clamp(
    double v,
    double lower,
    double upper
) {
    v = std::max(v, lower);
    v = std::min(v, upper);
    return v;
}


bool IntersectionFinderConvexPolygon2D::horizontalTest(CHEdge* ea, CHEdge* eb)
{
    double aL = ea->mSrc->mP.x();
    double aR = ea->mDst->mP.x();
    double bL = eb->mSrc->mP.x();
    double bR = eb->mDst->mP.x();
    if (aL>aR){
        std::swap(aL, aR);
    }
    if (bL>bR){
        std::swap(bL, bR);
    }

    if ( (aR + mEpsilonLinear < bL) || (bR + mEpsilonLinear < aL) ) {
        return false;
    }
    else {
        return true;
    }
}


bool IntersectionFinderConvexPolygon2D::boxTest(
    const Vec2& a1,
    const Vec2& a2,
    const Vec2& b1,
    const Vec2& b2
) {

    double aUpper = a1.y();
    double aLower = a2.y();
    double bUpper = b1.y();
    double bLower = b2.y();

    if (aLower>aUpper) {

        std::swap(aUpper, aLower);

    }
    if (bLower>bUpper) {

        std::swap(bUpper, bLower);

    }

    double aRight = a1.x();
    double aLeft  = a2.x();
    double bRight = b1.x();
    double bLeft  = b2.x();

    if (aLeft>aRight) {

        std::swap(aLeft, aRight);

    }
    if (bLeft>bRight) {

        std::swap(bLeft, bRight);

    }

    if ( (aUpper + mEpsilonLinear < bLower) || 
         (bUpper + mEpsilonLinear < aLower) ||
         (aRight + mEpsilonLinear < bLeft)  || 
         (bRight + mEpsilonLinear < aLeft)     ) {

        return false;

    }
    else {

        return true;

    }
}


bool IntersectionFinderConvexPolygon2D::boxTest(
    const Vec2& a1,
    const Vec2& a2,
    const Vec2& b
) {

    double aUpper = a1.y();
    double aLower = a2.y();
    double bY     = b.y();

    if (aLower>aUpper) {

        std::swap(aUpper, aLower);

    }

    double aRight = a1.x();
    double aLeft  = a2.x();
    double bX     = b.x();

    if (aLeft>aRight) {

        std::swap(aLeft, aRight);

    }

    if ( (aUpper + mEpsilonLinear < bY     ) || 
         (bY     + mEpsilonLinear < aLower ) ||
         (aRight + mEpsilonLinear < bX     ) || 
         (bX     + mEpsilonLinear < aLeft  )    ) {

        return false;

    }
    else {

        return true;

    }
}



/** @brief finds the perpendicular projection of a point on to an edge.
 *
 *  @param    PTest  (in): point to be tested.
 *
 *  @param    p1     (in): a terminal point on the edge.
 *
 *  @param    p2     (in): the other terminal point on the edge.
 *
 *  @param    s      (out):normalized directed length of the projected point
 *                         such that |pProj-p1| = s|p2 -p1| clamped to
 *                         [0.0, 1.0].
 *
 *  @return   the projected point on the edge.
 */
Vec2
IntersectionFinderConvexPolygon2D::findPerpendicularProjectionOfPointOntoEdge(
    const Vec2& PTest,
    const Vec2& p1,
    const Vec2& p2,
    double&     s
) {

    Vec2 v12 = p2 - p1;
    Vec2 v1p = PTest  - p1;
    double v12_sqNorm2 = v12.squaredNorm2();
    if (v12_sqNorm2 < mEpsilonSquared) {
        return p1;
    }
    else {
        s = clamp(v1p.dot(v12) / v12_sqNorm2, 0.0, 1.0);
        return p1 + (v12 * s);
    }

}

IntersectionFinderConvexPolygon2D::OutputElem 
IntersectionFinderConvexPolygon2D::constructOutputElem(CHNode* curV)
{
    long            indexA    = -1;
    long            indexAaux = -1;
    long            indexB    = -1;
    long            indexBaux = -1;
    enum VertexType type      = IT_NONE;
    if (curV->mThisIsA) {
        indexA    = curV->mIndex;
        indexAaux = curV->mIndexAux;
        if (curV->mPeer != nullptr) {
            indexB    = curV->mPeer->mIndex;
            indexBaux = curV->mPeer->mIndexAux;
        }
    }
    else {
        indexB    = curV->mIndex;
        indexBaux = curV->mIndexAux;
        if (curV->mPeer != nullptr) {
            indexA    = curV->mPeer->mIndex;
            indexAaux = curV->mPeer->mIndexAux;
        }
    }
    bool peerVmTypeIs_SPLIT_POINT = false;
    if (curV->mPeer != nullptr) {
        if (curV->mPeer->mType == CHNode::NT_SPLIT_POINT ) {
            peerVmTypeIs_SPLIT_POINT = true;
        }
    }
    if ( curV->mType == CHNode::NT_SPLIT_POINT ) {
        if (curV->mThisIsA) {
            type = IT_EDGE_VERTEX;
        }
        else {
            type = IT_VERTEX_EDGE;
        }
    }                 
    else if (peerVmTypeIs_SPLIT_POINT) {
        if (curV->mThisIsA) {
            type = IT_VERTEX_EDGE;
        }
        else {
            type = IT_EDGE_VERTEX;
        }
    }                 
    else if ( curV->mType == CHNode::NT_COINCIDENT_VERTEX) {
        type = IT_VERTEX_VERTEX;
    }
    else if ( curV->mType == CHNode::NT_INTERSECTION) {
        type = IT_EDGE_EDGE;
    }
    else {
        if (curV->mThisIsA) {
            type = IT_VERTEX_INTERIOR;
        }
        else {
            type = IT_INTERIOR_VERTEX;
        }
    }
    return OutputElem(
               curV->mP,
               indexA,
               indexAaux,
               indexB,
               indexBaux,
               type
           );
}


/** @brief visits the current vertex, finds the exploring direction and
 *         advance the current pointer curV.
 *
 *  @param curV     (in):  current vertex to be visited.
 *
 *  @param finished (out): indicates if the exploration is finished.
 *                         This is set if a touching point has been found
 *                         or no furuter cue for the next visit found.
 *
 *  @param the next vertex node to be visited.
 */
inline CHNode* IntersectionFinderConvexPolygon2D::visitVertexAndAdvanceForward(
    CHNode*     curV,
    bool&       finished,
    bool&       touchingPoint,
    bool&       interior,
    bool&       intsecFound,
    bool&       needToVisitBackward
) {

    auto* nextV = curV->mNext->mDst;

    if (curV->mPeer != nullptr) {
        // It has a peer. This means this is one of
        //  - intersection       (peer is intersection)
        //  - split point        (peer is a vertex on edge)
        //  - a vertex on edge   (peer is a split point)
        //  - coincident vertex  (peer is coincident vertex)

        intsecFound = true;

        if (curV->mIsTouchingPoint) {
            // There is only one touching point. Finish.
            touchingPoint  = true;
            finished       = true;
            return nullptr;

        }
        else if (nextV->mPeer != nullptr) {
            // Next vertex also has a peer. This means the forward incident
            // edge pair are colinear. Move foward.
            return nextV;

        }
        else {
            // Next vertex also has no peer.
            // Either curV's forward incident edge or peer's is
            // inside the intersection polygon.

            if (curV->mNext->mSrcInterior) {
                // curV's forward edge goes inside the other polygon
                interior            = true;
                needToVisitBackward = false;
                return nextV;

            }
            else if (curV->mPeer->mNext->mSrcInterior) {
                // The peer's forward edge goes inside the other polygon
                interior            = true;
                needToVisitBackward = false;
                return curV->mPeer->mNext->mDst;

            }
            else {
                // Neither curV or peer's forward edge goes inside.
                // They are going apart.
                finished = true;
                return nullptr;
            }
        }
    }
    else {
        // curV has no peer. This means curV is inside the other polygon.
        return nextV;
    }

}



/** @brief update the current vertices of A and B based on the current
 *         configuration and the sweeping direction.
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
inline void IntersectionFinderConvexPolygon2D::updateCurACurBCV(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {

    auto* candA = newVertexCandidate(*curA, leftHullA, sweepingDownward);
    auto* candB = newVertexCandidate(*curB, leftHullB, sweepingDownward);
//    cerr << "candA: (" << candA->mP.x() << "," << candA->mP.y() << ")\n";
//    cerr << "candB: (" << candB->mP.x() << "," << candB->mP.y() << ")\n";
    if (sweepingDownward) {
        if (!(*curA)->mBottom && !(*curB)->mBottom ) {
            if (candA->isHigherThan(candB)) {
                (*curA) = candA;
            }
            else {
                (*curB) = candB;
            }
        }
    }
    else {
        if (!(*curA)->mTop && !(*curB)->mTop) {
            if (candA->isHigherThan(candB)) {
                (*curB) = candB;
            }
            else {
                (*curA) = candA;            
            }
        }
    }
}


/** @brief update the current vertices of A and B based on the current
 *         configuration and the sweeping direction.
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
inline void IntersectionFinderConvexPolygon2D::updateCurACurBOE(
    CHNode** curA, 
    CHNode** curB, 
    bool     leftHullA, 
    bool     leftHullB, 
    bool     sweepingDownward
) {

    auto* candA = newVertexCandidate(*curA, leftHullA, sweepingDownward);
    auto* candB = newVertexCandidate(*curB, leftHullB, sweepingDownward);

    if (candA->isCoincident(candB)) {
        (*curA) = candA;
        (*curB) = candB;
    }
    if (sweepingDownward) {
        if (!(*curA)->mBottom && !(*curB)->mBottom ) {
            if (candA->isHigherThan(candB)) {
                (*curA) = candA;
            }
            else {
                (*curB) = candB;
            }
        }
    }
    else {
        if (!(*curA)->mTop && !(*curB)->mTop) {
            if (candA->isHigherThan(candB)) {
                (*curB) = candB;
            }
            else {
                (*curA) = candA;            
            }
        }
    }
}


/** @brief returns the vertex which will be visited next during a sweep.
 *         It depends on the side of the hull and the sweeping direction
 *
 *  @param cur (in): current vertex
 *
 *  @param leftHull (in): true if the hull is the left hull.
 *
 *  @param sweepingDownward (in): true if the sweeping direction is downward
 */
inline CHNode* IntersectionFinderConvexPolygon2D::newVertexCandidate(
                           CHNode* cur, bool leftHull, bool sweepingDownward) {
    if ( ( sweepingDownward && cur->mBottom) ||
         (!sweepingDownward && cur->mTop   )   ) {
        return cur;
    }
    if ( (leftHull && sweepingDownward) || (!leftHull && !sweepingDownward)) {
        
        return cur->mNext->mDst;
    }
    else {
        return cur->mPrev->mSrc;
    }
}


}// namespace Makena


#endif/*_MAKENA_INTERSECTION_CONVEX_POLYGON_2D_HPP_*/


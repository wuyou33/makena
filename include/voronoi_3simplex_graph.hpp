#ifndef _MAKENA_VORONOI_3SIMPLEX_GRAPH_HPP_
#define _MAKENA_VORONOI_3SIMPLEX_GRAPH_HPP_

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

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file voronoi_3simplex_graph.hpp
 *
 * @brief
 * Data structure and algorithm to perform
 * a point to 3-simplex classification test into one of Voronoi regions.
 */
namespace Makena {

using namespace std;

class Voronoi3SimplexGraph;
class Voronoi3SimplexEdge;


/** @class Voronoi3SimplexNode
 *
 *  @brief The node class for Voronoi3SimplexGraph.
 *         It represents a feature (triangle, edge, or vertex).
 */
class Voronoi3SimplexNode {

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    enum type {
        VERTEX,
        EDGE,
        FACE
    };

    inline Voronoi3SimplexNode(enum type t, long vertices);
    inline ~Voronoi3SimplexNode();

    inline void setNeighbors(
        Voronoi3SimplexEdge* e1,
        Voronoi3SimplexEdge* e2,
        Voronoi3SimplexEdge* e3,
        Voronoi3SimplexEdge* e4);

    inline void reset();

    inline enum predicate pred() const;


    inline void                 DFSreset();
    inline bool                 DFShasVisitedAll() const;
    inline Voronoi3SimplexEdge* DFSnextNeighbor();
    inline void                 DFSadvance();
    inline bool                 DFSisVisited() const;

    enum type             mType;

    /** @ brief represents the incident vertices in three digits.
     *
     *    1:  vertex 1
     *    2:  vertex 2
     *    3:  vertex 3
     *    4:  vertex 4
     *   12:  edge (1,2)
     *   23:  edge (2,3)
     *   31:  edge (3,1)
     *   14:  edge (1,4)
     *   24:  edge (2,4)
     *   34:  edge (3,4)
     *  132:  triangle(1,3,2)
     *  124:  triangle(1,2,4)
     *  143:  triangle(1,4,3)
     *  234:  triangle(2,3,4)
     */
    long                  mVertexSet;  

    /** @brief neighbor edges. */
    Voronoi3SimplexEdge*  (mNeighbors[4]);

    char                  mInward;
    /** @brief true if this feature is not visible from the test point.
     *         If this node is a face, then set to true if the plane is not
     *         facing the point.
     *         If this node is an edge, then set to true if both of the 
     *         incident faces are invisible.
     *         If this node is a vertex, then set to tru if all the incident 
     *         faces are invisible.
     *         If this is set to true, then this node is excluded from the
     *         BFS exploration.
     */
    bool                  mInvisible;     

    /** @brief used for BFS exploration */
    bool                  mQueued;

    /** @brief used for DFS exploration */
    bool                  mVisited;
    char                  mVisiting;

friend class Voronoi3SimplexGraph;
friend class Voronoi3SimplexEdge;
friend class GJKOriginFinder;

};


/** @class Voronoi3SimplexEdge
 *
 *  @brief The edge class for Voronoi3SimplexGraph.
 *         It represents a test at the boundary of two voronoi regions.
 */
class Voronoi3SimplexEdge {

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    enum type {
        UNKNOWN,
        FACING_NODE1,
        FACING_NODE2
    };

    inline Voronoi3SimplexEdge(Voronoi3SimplexGraph* g);
    inline ~Voronoi3SimplexEdge();
    inline void setNeighbors(Voronoi3SimplexNode* n1, Voronoi3SimplexNode* n2);

    inline void reset();
    inline Voronoi3SimplexNode* adjacentNode(Voronoi3SimplexNode* N) const;
    inline void performTest();
    inline enum type performTest_F132_E31();
    inline enum type performTest_F132_E23();
    inline enum type performTest_F132_E12();
    inline enum type performTest_F124_E12();
    inline enum type performTest_F124_E24();
    inline enum type performTest_F124_E14();
    inline enum type performTest_F143_E14();
    inline enum type performTest_F143_E34();
    inline enum type performTest_F143_E31();
    inline enum type performTest_F234_E23();
    inline enum type performTest_F234_E34();
    inline enum type performTest_F234_E24();
    inline enum type performTest_E12_V1();
    inline enum type performTest_E12_V2();
    inline enum type performTest_E23_V2();
    inline enum type performTest_E23_V3();
    inline enum type performTest_E31_V3();
    inline enum type performTest_E31_V1();
    inline enum type performTest_E14_V1();
    inline enum type performTest_E14_V4();
    inline enum type performTest_E24_V2();
    inline enum type performTest_E24_V4();
    inline enum type performTest_E34_V3();
    inline enum type performTest_E34_V4();

    /** @brief boundary test result */
    enum type     mType;

    Voronoi3SimplexGraph* mGraph;

    /** @brief incident face if this edge is (face, edge)
     *         incident edge if this edge is (edge, vertex)
     */
    Voronoi3SimplexNode*  mNode1;

    /** @brief incident edge if this edge is (face, edge)
     *         incident vertex if this edge is (edge, vertex)
     */
    Voronoi3SimplexNode*  mNode2; // Edge/Vertex

friend class Voronoi3SimplexGraph;
friend class Voronoi3SimplexNode;
friend class GJKOriginFinder;

};


/** @class Voronoi3SimplexGraph
 *  @brief peforms a classification test of a point into one of Voronoi
 *         regions of a given tetrahedron (3-simplex).
 *         The user of this class is expected to keep an object as
 *         construction of it is somewhat expensive for construction internal
 *         planar graph structure.
 *
 *  @remark
 *         Nodes represent features of tetrahedron.
 *         Edges indicate the boundary to be tested.
 *         We find the orientation of the edges based on the
 *         boundary test.
 *         There should be no cycle in the graph.
 *         There should be only one sink node, and it indicates
 *         the Voronoi feature of this tetrahedron for the given point.
 *
 *         +-----------------V2-----------------+
 *         |                 |                  |
 *         |                 |                  |
 *         |                 |                  |
 *         |       +--------E24---------+       |
 *         |       |         |          |       |
 *         | +---F234        |         F124---+ |
 *         | |      |        |         |      | |
 *         | |      |    +---V4---+    |      | |
 *         | |      |    |   |    |    |      | |
 *         E23---+  +--E34   |    E14--+  +---E12
 *           |   |     | |   |    | |     |   |
 *           |   |     | +--F143--+ |     |   |
 *           |   |     |     |      |     |   |
 *           |   +----V3     |      V1----+   |
 *           |         |     |      |         |
 *           |         +----E31-----+         |
 *           |               |                |
 *           |               |                |
 *           |               |                |
 *           +--------------F132--------------+
 *
 */
class Voronoi3SimplexGraph {

  public:
    inline Voronoi3SimplexGraph();
    inline ~Voronoi3SimplexGraph();

    inline enum predicate find(
        const Vec3& p1,
        const Vec3& p2,
        const Vec3& p3,
        const Vec3& p4,
        const Vec3& pTest
    );

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    inline void reset(
        const Vec3& p1,
        const Vec3& p2,
        const Vec3& p3,
        const Vec3& p4,
        const Vec3& pTest
    );


    /** @brief */
    inline void FindNormalsAndTestVisibility();


    /** @brief starting from visible feature, find the voronoi region
     *         by BFS on the graph oriented according to the boundary tests
     *         along the edges. Justification for thos BFS to work is that
     *         there is at most one sink, which is the boronoi region for the
     *         test point. From an arbitrary visible feature (node), if a BFS
     *         exploration did not reach the sink, then it means there would be
     *         another sink, which contradicts the assumption.
     */
    inline enum predicate findVoronoiRegionBFS();

    /** @brief if BFS has failed to detect a Voronoi region, there must be 
     *         a cycle in the graph. This function tries to detect a cycle
     *         by DFS and returns the appropriate Voronoi region.
     */
    inline enum predicate detectCycle();

    inline enum predicate findMajority(vector<Voronoi3SimplexNode*>& S);

    Vec3                mP1;
    Vec3                mP2;
    Vec3                mP3;
    Vec3                mP4;
    Vec3                mF132;
    Vec3                mF124;
    Vec3                mF143;
    Vec3                mF234;
    Vec3                mPT;

    Voronoi3SimplexNode V1;
    Voronoi3SimplexNode V2;
    Voronoi3SimplexNode V3;
    Voronoi3SimplexNode V4;
    Voronoi3SimplexNode E12;
    Voronoi3SimplexNode E23;
    Voronoi3SimplexNode E31;
    Voronoi3SimplexNode E14;
    Voronoi3SimplexNode E24;
    Voronoi3SimplexNode E34;
    Voronoi3SimplexNode F132;
    Voronoi3SimplexNode F124;
    Voronoi3SimplexNode F143;
    Voronoi3SimplexNode F234;

    Voronoi3SimplexEdge F132_E12;
    Voronoi3SimplexEdge F132_E23;
    Voronoi3SimplexEdge F132_E31;
    Voronoi3SimplexEdge F124_E12;
    Voronoi3SimplexEdge F124_E24;
    Voronoi3SimplexEdge F124_E14;
    Voronoi3SimplexEdge F143_E31;
    Voronoi3SimplexEdge F143_E14;
    Voronoi3SimplexEdge F143_E34;
    Voronoi3SimplexEdge F234_E23;
    Voronoi3SimplexEdge F234_E34;
    Voronoi3SimplexEdge F234_E24;
    Voronoi3SimplexEdge E12_V1;
    Voronoi3SimplexEdge E12_V2;
    Voronoi3SimplexEdge E23_V2;
    Voronoi3SimplexEdge E23_V3;
    Voronoi3SimplexEdge E31_V3;
    Voronoi3SimplexEdge E31_V1;
    Voronoi3SimplexEdge E14_V1;
    Voronoi3SimplexEdge E14_V4;
    Voronoi3SimplexEdge E24_V2;
    Voronoi3SimplexEdge E24_V4;
    Voronoi3SimplexEdge E34_V3;
    Voronoi3SimplexEdge E34_V4;

friend class Voronoi3SimplexNode;
friend class Voronoi3SimplexEdge;
friend class GJKOriginFinder;

};


Voronoi3SimplexNode::Voronoi3SimplexNode(enum type t, long vertices):
    mType      (t),
    mVertexSet (vertices),
    mInward    (0),
    mInvisible (false),
    mQueued    (false),
    mVisited   (false),
    mVisiting  (0)     {;}
    

Voronoi3SimplexNode::~Voronoi3SimplexNode(){;}


void Voronoi3SimplexNode::setNeighbors(
    Voronoi3SimplexEdge* e1,
    Voronoi3SimplexEdge* e2,
    Voronoi3SimplexEdge* e3,
    Voronoi3SimplexEdge* e4
) {
    mNeighbors[0] = e1;
    mNeighbors[1] = e2;
    mNeighbors[2] = e3;
    mNeighbors[3] = e4;
}


void Voronoi3SimplexNode::reset() {
    mInward    = 0;
    mInvisible = false;
    mQueued    = false;
    mVisited   = false;
    mVisiting  = 0;
}


enum predicate Voronoi3SimplexNode::pred() const
{
    switch (mVertexSet) {

      case 1:
        return VORONOI_OVER_VERTEX_1;
        break;

      case 2:
        return VORONOI_OVER_VERTEX_2;
        break;

      case 3:
        return VORONOI_OVER_VERTEX_3;
        break;

      case 4:
        return VORONOI_OVER_VERTEX_4;
        break;

      case 12:
        return VORONOI_OVER_EDGE_1_2;
        break;

      case 14:
        return VORONOI_OVER_EDGE_1_4;
        break;

      case 23:
        return VORONOI_OVER_EDGE_2_3;
        break;

      case 24:
        return VORONOI_OVER_EDGE_2_4;
        break;

      case 31:
        return VORONOI_OVER_EDGE_3_1;
        break;

      case 34:
        return VORONOI_OVER_EDGE_3_4;
        break;

      case 124:
        return VORONOI_OVER_TRIANGLE_1_2_4;
        break;

      case 132:
        return VORONOI_OVER_TRIANGLE_1_3_2;
        break;

      case 143:
        return VORONOI_OVER_TRIANGLE_1_4_3;
        break;

      case 234:
        return VORONOI_OVER_TRIANGLE_2_3_4;
        break;

      default:
        return NONE;

    }
    return NONE;
}   


void Voronoi3SimplexNode::DFSreset()
{
   mVisited  = true;
   mVisiting = 0;
}


bool Voronoi3SimplexNode::DFShasVisitedAll() const
{
    return (mType == EDGE && mVisiting == 4) ||
           (mType != EDGE && mVisiting == 3) ;
}


Voronoi3SimplexEdge* Voronoi3SimplexNode::DFSnextNeighbor()
{
    return mNeighbors[size_t(mVisiting)];
}


void Voronoi3SimplexNode::DFSadvance()
{
    mVisiting++;
}


bool Voronoi3SimplexNode::DFSisVisited() const
{
    return mVisited;
}


Voronoi3SimplexEdge::Voronoi3SimplexEdge(Voronoi3SimplexGraph* g):
    mType  (UNKNOWN),
    mGraph (g),
    mNode1 (nullptr),
    mNode2 (nullptr) {;}


Voronoi3SimplexEdge::~Voronoi3SimplexEdge(){;}


void Voronoi3SimplexEdge::setNeighbors(
    Voronoi3SimplexNode* n1, 
    Voronoi3SimplexNode* n2
) {
    mNode1 = n1;
    mNode2 = n2;
}


void Voronoi3SimplexEdge::reset() {
    mType = UNKNOWN;
}


Voronoi3SimplexNode* Voronoi3SimplexEdge::adjacentNode(Voronoi3SimplexNode* N)
const
{
    return  (N == mNode1) ? mNode2 : mNode1;
}


void Voronoi3SimplexEdge::performTest()
{
    // Perform one of 24 tests.
    if (mNode1->mInvisible) {
        mType = FACING_NODE2;
        mNode2->mInward++;
        return;
    }
    else if (mNode2->mInvisible) {
        mType = FACING_NODE1;
        mNode1->mInward++;
        return;
    }

    auto v2  = mNode2->mVertexSet;

    switch (mNode1->mVertexSet) {

      case 132:

        if (v2==31) {
            mType = performTest_F132_E31();
        }
        else if (v2==23) {
            mType = performTest_F132_E23();
        }
        else {//v2 == 12
            mType = performTest_F132_E12();
        }
        break;

      case 124:

        if (v2==12) {
            mType = performTest_F124_E12();
        }
        else if (v2==24) {
            mType = performTest_F124_E24();
        }
        else {//v2 == 14
            mType = performTest_F124_E14();
        }
        break;

      case 143:

        if (v2==14) {
            mType = performTest_F143_E14();
        }
        else if (v2==34) {
            mType = performTest_F143_E34();
        }
        else {//v2 == 31
            mType = performTest_F143_E31();
        }
        break;

      case 234:
        if (v2==23) {
            mType = performTest_F234_E23();
        }
        else if (v2==34) {
            mType = performTest_F234_E34();
        }
        else {//v2 == 24
            mType = performTest_F234_E24();
        }
        break;

      case 12:

        if (v2==1) {
            mType = performTest_E12_V1();
        }
        else {//v2 == 2
            mType = performTest_E12_V2();
        }
        break;

      case 23:

        if (v2==2) {
            mType = performTest_E23_V2();
        }
        else {//v2 == 3
            mType = performTest_E23_V3();
        }
        break;

      case 31:

        if (v2==3) {
            mType = performTest_E31_V3();
        }
        else {//v2 == 1
            mType = performTest_E31_V1();
        }
        break;

      case 14:

        if (v2==1) {
            mType = performTest_E14_V1();
        }
        else {//v2 == 4
            mType = performTest_E14_V4();
        }
        break;

      case 24:

        if (v2==2) {
            mType = performTest_E24_V2();
        }
        else {//v2 == 4
            mType = performTest_E24_V4();
        }
        break;

      case 34:

        if (v2==3) {
            mType = performTest_E34_V3();
        }
        else {//v2 == 4
            mType = performTest_E34_V4();
        }
        break;

      default:
        break;      
    }    

    
    if (mType == FACING_NODE1) {
        mNode1->mInward++;
    }
    else if (mType == FACING_NODE2) {
        mNode2->mInward++;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_F132_E31()
{
    // The testing plane is along edge (3,1) perpendicular to face (1,3,2)
    // Base vertex is 1.

    Vec3   v1t  = mGraph->mPT - mGraph->mP1;
    Vec3   v13  = mGraph->mP3 - mGraph->mP1;
    Vec3   vRef = v13.cross(mGraph->mF132);
    double dRef = vRef.dot(v1t);

    if (dRef > 0.0) {
        return FACING_NODE2;
    }
    else {
        return FACING_NODE1;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_F132_E23()
{
    // The testing plane is along edge (2,3) perpendicular to face (1,3,2)
    // Base vertex is 2.

    Vec3   v2t  = mGraph->mPT - mGraph->mP2;
    Vec3   v23  = mGraph->mP3 - mGraph->mP2;
    Vec3   vRef = mGraph->mF132.cross(v23);
    double dRef = vRef.dot(v2t);

    if (dRef > 0.0) {
        return FACING_NODE2;
    }
    else {
        return FACING_NODE1;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_F132_E12()
{
    // The testing plane is along edge (1,2) perpendicular to face (1,3,2)
    // Base vertex is 1.
    Vec3   v1t  = mGraph->mPT - mGraph->mP1;
    Vec3   v12  = mGraph->mP2 - mGraph->mP1;
    Vec3   vRef = mGraph->mF132.cross(v12);
    double dRef = vRef.dot(v1t);
    if (dRef > 0.0) {
        return FACING_NODE2;
    }
    else {
        return FACING_NODE1;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_F124_E12()
{
    // The testing plane is along edge (1,2) perpendicular to face (1,2,4)
    // Base vertex is 1.
    Vec3   v1t  = mGraph->mPT - mGraph->mP1;
    Vec3   v12  = mGraph->mP2 - mGraph->mP1;
    Vec3   vRef = v12.cross(mGraph->mF124);
    double dRef = vRef.dot(v1t);
    if (dRef > 0.0) {
        return FACING_NODE2;
    }
    else {
        return FACING_NODE1;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_F124_E24()
{
    // The testing plane is along edge (2,4) perpendicular to face (1,2,4)
    // Base vertex is 2.
    Vec3   v2t  = mGraph->mPT - mGraph->mP2;
    Vec3   v24  = mGraph->mP4 - mGraph->mP2;
    Vec3   vRef = v24.cross(mGraph->mF124);
    double dRef = vRef.dot(v2t);
    if (dRef > 0.0) {
        return FACING_NODE2;
    }
    else {
        return FACING_NODE1;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_F124_E14()
{
    // The testing plane is along edge (1,4) perpendicular to face (1,2,4)
    // Base vertex is 1.
    Vec3   v1t  = mGraph->mPT - mGraph->mP1;
    Vec3   v14  = mGraph->mP4 - mGraph->mP1;
    Vec3   vRef = mGraph->mF124.cross(v14);
    double dRef = vRef.dot(v1t);
    if (dRef > 0.0) {
        return FACING_NODE2;
    }
    else {
        return FACING_NODE1;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_F143_E14()
{
    // The testing plane is along edge (1,4) perpendicular to face (1,4,3)
    // Base vertex is 1.
    Vec3   v1t  = mGraph->mPT - mGraph->mP1;
    Vec3   v14  = mGraph->mP4 - mGraph->mP1;
    Vec3   vRef = v14.cross(mGraph->mF143);
    double dRef = vRef.dot(v1t);
    if (dRef > 0.0) {
        return FACING_NODE2;
    }
    else {
        return FACING_NODE1;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_F143_E34()
{
    // The testing plane is along edge (3,4) perpendicular to face (1,4,3)
    // Base vertex is 3.
    Vec3   v3t  = mGraph->mPT - mGraph->mP3;
    Vec3   v34  = mGraph->mP4 - mGraph->mP3;
    Vec3   vRef = mGraph->mF143.cross(v34);
    double dRef = vRef.dot(v3t);
    if (dRef > 0.0) {
        return FACING_NODE2;
    }
    else {
        return FACING_NODE1;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_F143_E31()
{
    // The testing plane is along edge (3,1) perpendicular to face (1,4,3)
    // Base vertex is 1.
    Vec3   v1t  = mGraph->mPT - mGraph->mP1;
    Vec3   v13  = mGraph->mP3 - mGraph->mP1;
    Vec3   vRef = mGraph->mF143.cross(v13);
    double dRef = vRef.dot(v1t);
    if (dRef > 0.0) {
        return FACING_NODE2;
    }
    else {
        return FACING_NODE1;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_F234_E23()
{
    // The testing plane is along edge (2,3) perpendicular to face (2,3,4)
    // Base vertex is 2.
    Vec3   v2t  = mGraph->mPT - mGraph->mP2;
    Vec3   v23  = mGraph->mP3 - mGraph->mP2;
    Vec3   vRef = v23.cross(mGraph->mF234);
    double dRef = vRef.dot(v2t);
    if (dRef > 0.0) {
        return FACING_NODE2;
    }
    else {
        return FACING_NODE1;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_F234_E34()
{
    // The testing plane is along edge (3,4) perpendicular to face (2,3,4)
    // Base vertex is 3.
    Vec3   v3t  = mGraph->mPT - mGraph->mP3;
    Vec3   v34  = mGraph->mP4 - mGraph->mP3;
    Vec3   vRef = v34.cross(mGraph->mF234);
    double dRef = vRef.dot(v3t);
    if (dRef > 0.0) {
        return FACING_NODE2;
    }
    else {
        return FACING_NODE1;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_F234_E24()
{
    // The testing plane is along edge (2,4) perpendicular to face (2,3,4)
    // Base vertex is 2.
    Vec3   v2t  = mGraph->mPT - mGraph->mP2;
    Vec3   v24  = mGraph->mP4 - mGraph->mP2;
    Vec3   vRef = mGraph->mF234.cross(v24);
    double dRef = vRef.dot(v2t);
    if (dRef > 0.0) {
        return FACING_NODE2;
    }
    else {
        return FACING_NODE1;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_E12_V1()
{
    // The testing plane perpendicularly intersects vertex 1.
    // Base vertex is 1.
    Vec3   v1t  = mGraph->mPT - mGraph->mP1;
    Vec3   v12  = mGraph->mP2 - mGraph->mP1;
    double dRef = v12.dot(v1t);
    if (dRef > 0.0) {
        return FACING_NODE1;
    }
    else {
        return FACING_NODE2;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_E12_V2()
{
    // The testing plane perpendicularly intersects vertex 2.
    // Base vertex is 2.
    Vec3   v2t  = mGraph->mPT - mGraph->mP2;
    Vec3   v21  = mGraph->mP1 - mGraph->mP2;
    double dRef = v21.dot(v2t);
    if (dRef > 0.0) {
        return FACING_NODE1;
    }
    else {
        return FACING_NODE2;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_E23_V2()
{
    // The testing plane perpendicularly intersects vertex 2.
    // Base vertex is 2.
    Vec3   v2t  = mGraph->mPT - mGraph->mP2;
    Vec3   v23  = mGraph->mP3 - mGraph->mP2;
    double dRef = v23.dot(v2t);
    if (dRef > 0.0) {
        return FACING_NODE1;
    }
    else {
        return FACING_NODE2;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_E23_V3()
{
    // The testing plane perpendicularly intersects vertex 3.
    // Base vertex is 3.
    Vec3   v3t  = mGraph->mPT - mGraph->mP3;
    Vec3   v32  = mGraph->mP2 - mGraph->mP3;
    double dRef = v32.dot(v3t);
    if (dRef > 0.0) {
        return FACING_NODE1;
    }
    else {
        return FACING_NODE2;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_E31_V3()
{
    // The testing plane perpendicularly intersects vertex 3.
    // Base vertex is 3.
    Vec3   v3t  = mGraph->mPT - mGraph->mP3;
    Vec3   v31  = mGraph->mP1 - mGraph->mP3;
    double dRef = v31.dot(v3t);
    if (dRef > 0.0) {
        return FACING_NODE1;
    }
    else {
        return FACING_NODE2;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_E31_V1()
{
    // The testing plane perpendicularly intersects vertex 1.
    // Base vertex is 1.
    Vec3   v1t  = mGraph->mPT - mGraph->mP1;
    Vec3   v13  = mGraph->mP3 - mGraph->mP1;
    double dRef = v13.dot(v1t);
    if (dRef > 0.0) {
        return FACING_NODE1;
    }
    else {
        return FACING_NODE2;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_E14_V1()
{
    // The testing plane perpendicularly intersects vertex 1.
    // Base vertex is 1.
    Vec3   v1t  = mGraph->mPT - mGraph->mP1;
    Vec3   v14  = mGraph->mP4 - mGraph->mP1;
    double dRef = v14.dot(v1t);
    if (dRef > 0.0) {
        return FACING_NODE1;
    }
    else {
        return FACING_NODE2;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_E14_V4()
{
    // The testing plane perpendicularly intersects vertex 4.
    // Base vertex is 4.
    Vec3   v4t  = mGraph->mPT - mGraph->mP4;
    Vec3   v41  = mGraph->mP1 - mGraph->mP4;
    double dRef = v41.dot(v4t);
    if (dRef > 0.0) {
        return FACING_NODE1;
    }
    else {
        return FACING_NODE2;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_E24_V2()
{
    // The testing plane perpendicularly intersects vertex 2.
    // Base vertex is 2.
    Vec3   v2t  = mGraph->mPT - mGraph->mP2;
    Vec3   v24  = mGraph->mP4 - mGraph->mP2;
    double dRef = v24.dot(v2t);
    if (dRef > 0.0) {
        return FACING_NODE1;
    }
    else {
        return FACING_NODE2;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_E24_V4()
{
    // The testing plane perpendicularly intersects vertex 4.
    // Base vertex is 4.
    Vec3   v4t  = mGraph->mPT - mGraph->mP4;
    Vec3   v42  = mGraph->mP2 - mGraph->mP4;
    double dRef = v42.dot(v4t);
    if (dRef > 0.0) {
        return FACING_NODE1;
    }
    else {
        return FACING_NODE2;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_E34_V3()
{
    // The testing plane perpendicularly intersects vertex 3.
    // Base vertex is 3.
    Vec3   v3t  = mGraph->mPT - mGraph->mP3;
    Vec3   v34  = mGraph->mP4 - mGraph->mP3;
    double dRef = v34.dot(v3t);
    if (dRef > 0.0) {
        return FACING_NODE1;
    }
    else {
        return FACING_NODE2;
    }
}


enum Voronoi3SimplexEdge::type Voronoi3SimplexEdge::performTest_E34_V4()
{
    // The testing plane perpendicularly intersects vertex 4.
    // Base vertex is 4.
    Vec3   v4t  = mGraph->mPT - mGraph->mP4;
    Vec3   v43  = mGraph->mP3 - mGraph->mP4;
    double dRef = v43.dot(v4t);
    if (dRef > 0.0) {
        return FACING_NODE1;
    }
    else {
        return FACING_NODE2;
    }
}


Voronoi3SimplexGraph::Voronoi3SimplexGraph():
    V1  (Voronoi3SimplexNode::VERTEX,   1),
    V2  (Voronoi3SimplexNode::VERTEX,   2),
    V3  (Voronoi3SimplexNode::VERTEX,   3),
    V4  (Voronoi3SimplexNode::VERTEX,   4),
    E12 (Voronoi3SimplexNode::EDGE,    12),
    E23 (Voronoi3SimplexNode::EDGE,    23),
    E31 (Voronoi3SimplexNode::EDGE,    31),
    E14 (Voronoi3SimplexNode::EDGE,    14),
    E24 (Voronoi3SimplexNode::EDGE,    24),
    E34 (Voronoi3SimplexNode::EDGE,    34),
    F132(Voronoi3SimplexNode::FACE,   132),
    F124(Voronoi3SimplexNode::FACE,   124),
    F143(Voronoi3SimplexNode::FACE,   143),
    F234(Voronoi3SimplexNode::FACE,   234),
    F132_E12(this),
    F132_E23(this),
    F132_E31(this),
    F124_E12(this),
    F124_E24(this),
    F124_E14(this),
    F143_E31(this),
    F143_E14(this),
    F143_E34(this),
    F234_E23(this),
    F234_E34(this),
    F234_E24(this),
    E12_V1(this),
    E12_V2(this),
    E23_V2(this),
    E23_V3(this),
    E31_V3(this),
    E31_V1(this),
    E14_V1(this),
    E14_V4(this),
    E24_V2(this),
    E24_V4(this),
    E34_V3(this),
    E34_V4(this)
{

    V1.  setNeighbors(&E12_V1,   &E31_V1,   &E14_V1,   nullptr);
    V2.  setNeighbors(&E12_V2,   &E23_V2,   &E24_V2,   nullptr);
    V3.  setNeighbors(&E23_V3,   &E31_V3,   &E34_V3,   nullptr);
    V4.  setNeighbors(&E14_V4,   &E24_V4,   &E34_V4,   nullptr);

    E12. setNeighbors(&F132_E12, &F124_E12, &E12_V1,   &E12_V2);
    E23. setNeighbors(&F132_E23, &F234_E23, &E23_V2,   &E23_V3);
    E31. setNeighbors(&F132_E31, &F143_E31, &E31_V3,   &E31_V1);
    E14. setNeighbors(&F124_E14, &F143_E14, &E14_V1,   &E14_V4);
    E24. setNeighbors(&F124_E24, &F234_E24, &E24_V2,   &E24_V4);
    E34. setNeighbors(&F143_E34, &F234_E34, &E34_V3,   &E34_V4);

    F132.setNeighbors(&F132_E12, &F132_E23, &F132_E31, nullptr);
    F124.setNeighbors(&F124_E12, &F124_E24, &F124_E14, nullptr);
    F143.setNeighbors(&F143_E31, &F143_E14, &F143_E34, nullptr);
    F234.setNeighbors(&F234_E23, &F234_E34, &F234_E24, nullptr);

    F132_E12.setNeighbors(&F132, &E12);
    F132_E23.setNeighbors(&F132, &E23);
    F132_E31.setNeighbors(&F132, &E31);
    F124_E12.setNeighbors(&F124, &E12);
    F124_E24.setNeighbors(&F124, &E24);
    F124_E14.setNeighbors(&F124, &E14);
    F143_E31.setNeighbors(&F143, &E31);
    F143_E14.setNeighbors(&F143, &E14);
    F143_E34.setNeighbors(&F143, &E34);
    F234_E23.setNeighbors(&F234, &E23);
    F234_E34.setNeighbors(&F234, &E34);
    F234_E24.setNeighbors(&F234, &E24);
    E12_V1  .setNeighbors(&E12,  &V1);
    E12_V2  .setNeighbors(&E12,  &V2);
    E23_V2  .setNeighbors(&E23,  &V2);
    E23_V3  .setNeighbors(&E23,  &V3);
    E31_V3  .setNeighbors(&E31,  &V3);
    E31_V1  .setNeighbors(&E31,  &V1);
    E14_V1  .setNeighbors(&E14,  &V1);
    E14_V4  .setNeighbors(&E14,  &V4);
    E24_V2  .setNeighbors(&E24,  &V2);
    E24_V4  .setNeighbors(&E24,  &V4);
    E34_V3  .setNeighbors(&E34,  &V3);
    E34_V4  .setNeighbors(&E34,  &V4);

}


Voronoi3SimplexGraph::~Voronoi3SimplexGraph(){;}


enum predicate Voronoi3SimplexGraph::find(
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& p3,
    const Vec3& p4,
    const Vec3& pTest
) {

    reset(p1, p2, p3, p4, pTest);

    FindNormalsAndTestVisibility();

    return findVoronoiRegionBFS();

}


void Voronoi3SimplexGraph::reset(
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& p3,
    const Vec3& p4,
    const Vec3& pTest
) {

    mP1   = p1;
    mP2   = p2;
    mP3   = p3;
    mP4   = p4;
    mPT   = pTest;

    V1.reset();
    V2.reset();
    V3.reset();
    V4.reset();
    E12.reset();
    E23.reset();
    E31.reset();
    E14.reset();
    E24.reset();
    E34.reset();
    F132.reset();
    F124.reset();
    F143.reset();
    F234.reset();

    F132_E12.reset();
    F132_E23.reset();
    F132_E31.reset();
    F124_E12.reset();
    F124_E24.reset();
    F124_E14.reset();
    F143_E31.reset();
    F143_E14.reset();
    F143_E34.reset();
    F234_E23.reset();
    F234_E34.reset();
    F234_E24.reset();
    E12_V1.reset();
    E12_V2.reset();
    E23_V2.reset();
    E23_V3.reset();
    E31_V3.reset();
    E31_V1.reset();
    E14_V1.reset();
    E14_V4.reset();
    E24_V2.reset();
    E24_V4.reset();
    E34_V3.reset();
    E34_V4.reset();
}


void Voronoi3SimplexGraph::FindNormalsAndTestVisibility()
{

    // F132
    const Vec3 v12    = mP2 - mP1;
    const Vec3 v13    = mP3 - mP1;
    const Vec3 v1Test = mPT - mP1;

    mF132 = v13.cross(v12);
    if (mF132.dot(v1Test)<=0.0) {
        F132.mInvisible = true;
    }

    // F124
    const Vec3 v14        = mP4 - mP1;
    mF124 = v12.cross(v14);
    if (mF124.dot(v1Test)<=0.0) {
        F124.mInvisible = true;
    }

    // F143
    mF143 = v14.cross(v13);
    if (mF143.dot(v1Test)<=0.0) {
        F143.mInvisible = true;
    }

    // F234
    const Vec3 v23        = mP3 - mP2;
    const Vec3 v24        = mP4 - mP2;
    const Vec3 v2Test     = mPT - mP2;
    mF234 = v23.cross(v24);
    if (mF234.dot(v2Test)<=0.0) {
        F234.mInvisible = true;
    }

    // E12
    if (F132.mInvisible&&F124.mInvisible) {
        E12.mInvisible=true;
    }

    // E23
    if (F132.mInvisible&&F234.mInvisible) {
        E23.mInvisible=true;
    }

    // E31
    if (F132.mInvisible&&F143.mInvisible) {
        E31.mInvisible=true;
    }

    // E14
    if (F124.mInvisible&&F143.mInvisible) {
        E14.mInvisible=true;
    }

    // E24
    if (F124.mInvisible&&F234.mInvisible) {
        E24.mInvisible=true;
    }

    // E34
    if (F143.mInvisible&&F234.mInvisible) {
        E34.mInvisible=true;
    }

    // V1
    if (F132.mInvisible&&F143.mInvisible&&F124.mInvisible) {
        V1.mInvisible=true;
    }

    // V2
    if (F132.mInvisible&&F124.mInvisible&&F234.mInvisible) {
        V2.mInvisible=true;
    }

    // V3
    if (F132.mInvisible&&F143.mInvisible&&F234.mInvisible) {
        V3.mInvisible=true;
    }

    // V4
    if (F124.mInvisible&&F143.mInvisible&&F234.mInvisible) {
        V4.mInvisible=true;
    }
}


enum predicate Voronoi3SimplexGraph::findVoronoiRegionBFS()
{

    list<Voronoi3SimplexNode*> Q;

//cerr << "F132.mInvisible" << F132.mInvisible << "\n";
//cerr << "F124.mInvisible" << F124.mInvisible << "\n";
//cerr << "F143.mInvisible" << F143.mInvisible << "\n";
//cerr << "F234.mInvisible" << F234.mInvisible << "\n";

    if (!F132.mInvisible) {
        Q.push_back(&F132);
        F132.mQueued = true;
    }
    else if (!F124.mInvisible) {
        Q.push_back(&F124);
        F124.mQueued = true;
    }
    else if (!F143.mInvisible) {
        Q.push_back(&F143);
        F143.mQueued = true;
    }
    else if (!F234.mInvisible) {
        Q.push_back(&F234);
        F234.mQueued = true;
    }
    else {
        return VORONOI_INSIDE_TETRAHEDRON;
    }

    while (!Q.empty()) {
//cerr << "loop head\n";

        Voronoi3SimplexNode* N = Q.front();
        Q.pop_front();
//cerr << "N: " << N->mVertexSet << "\n";
        if ( (N->mType == Voronoi3SimplexNode::EDGE && N->mInward==4)||
             (N->mType != Voronoi3SimplexNode::EDGE && N->mInward==3)  ) {
            return N->pred();
        }

        for (auto* e : N->mNeighbors) {

            if (e!=nullptr) {

                if (e->mType == Voronoi3SimplexEdge::UNKNOWN) {
                    e->performTest();
                }

//cerr << "e->mNode1: " <<e->mNode1->mVertexSet << "\n";
//cerr << "e->mNode2: " << e->mNode2->mVertexSet << "\n";
//cerr << "e->mType: " << ((e->mType == Voronoi3SimplexEdge::FACING_NODE1)?
//                        " NODE1":" NODE2") << "\n";

                if ( ( e->mType == Voronoi3SimplexEdge::FACING_NODE2 &&
                       e->mNode1==N                             ) ||
                     ( e->mType == Voronoi3SimplexEdge::FACING_NODE1 &&
                       e->mNode2==N                             )    ) {
                    Voronoi3SimplexNode* A = e->adjacentNode(N);
                    if (!A->mQueued && !A->mInvisible) {
                        Q.push_back(A);
                        A->mQueued = true;
                    }
                }
            }
        }
        if ( (N->mType == Voronoi3SimplexNode::EDGE && N->mInward==4)||
             (N->mType != Voronoi3SimplexNode::EDGE && N->mInward==3)  ) {
            return N->pred();
        }
//cerr << "loop end\n";
    }

//cerr << "failed to detect\n";

    return detectCycle();
}


enum predicate Voronoi3SimplexGraph::detectCycle()
{
    vector<Voronoi3SimplexNode*> S;
    S.reserve(14);

    if (!F132.mInvisible) {
        S.push_back(&F132);
        F132.DFSreset();
    }
    else if (!F124.mInvisible) {
        S.push_back(&F124);
        F124.DFSreset();
    }
    else if (!F143.mInvisible) {
        S.push_back(&F143);
        F143.DFSreset();
    }
    else if (!F234.mInvisible) {
        S.push_back(&F234);
        F234.DFSreset();
    }

    while(!S.empty()) {

        Voronoi3SimplexNode* N = S.back();

        if (N->DFShasVisitedAll()) {
            S.pop_back();
        }
        else {

            auto* e = N->DFSnextNeighbor();

            N->DFSadvance();

            if ((e->mType==Voronoi3SimplexEdge::FACING_NODE2 && e->mNode1==N)||
                (e->mType==Voronoi3SimplexEdge::FACING_NODE1 && e->mNode2==N)){
                Voronoi3SimplexNode* A = e->adjacentNode(N);
                if (A->DFSisVisited()){
                    // Cycle detected
                    break;
                }
                else {
                    A->DFSreset();
                    S.push_back(A);
                }
            }
        }
    }

    return findMajority(S);
}


enum predicate Voronoi3SimplexGraph::findMajority(
    vector<Voronoi3SimplexNode*>& S
) {
    char v1Cnt = 0;
    char v2Cnt = 0;
    char v3Cnt = 0;
    char v4Cnt = 0;

    for (auto* n : S) {
        switch (n->mVertexSet) {

          case 1:
            v1Cnt++;
            break;

          case 2:
            v2Cnt++;
            break;

          case 3:
            v3Cnt++;
            break;

          case 4:
            v4Cnt++;
            break;

          case 12:
            v1Cnt++;
            v2Cnt++;
            break;

          case 23:
            v2Cnt++;
            v3Cnt++;
            break;

          case 31:
            v3Cnt++;
            v1Cnt++;
            break;

          case 14:
            v1Cnt++;
            v4Cnt++;
            break;

          case 24:
            v2Cnt++;
            v4Cnt++;
            break;

          case 34:
            v3Cnt++;
            v4Cnt++;
            break;

          case 132:
            v1Cnt++;
            v3Cnt++;
            v2Cnt++;
            break;

          case 124:
            v1Cnt++;
            v2Cnt++;
            v4Cnt++;
            break;

          case 143:
            v1Cnt++;
            v4Cnt++;
            v3Cnt++;
            break;

          case 234:
            v2Cnt++;
            v3Cnt++;
            v4Cnt++;
            break;

          default:
            break;
        }
    }

    char maxVal = v1Cnt;
    maxVal = std::max(v2Cnt,maxVal);
    maxVal = std::max(v3Cnt,maxVal);
    maxVal = std::max(v4Cnt,maxVal);

    if (v1Cnt==maxVal && v2Cnt <maxVal && v3Cnt <maxVal && v4Cnt <maxVal)
        return VORONOI_OVER_VERTEX_1;

    if (v1Cnt <maxVal && v2Cnt==maxVal && v3Cnt <maxVal && v4Cnt <maxVal)
        return VORONOI_OVER_VERTEX_2;

    if (v1Cnt <maxVal && v2Cnt <maxVal && v3Cnt==maxVal && v4Cnt <maxVal)
        return VORONOI_OVER_VERTEX_3;

    if (v1Cnt <maxVal && v2Cnt <maxVal && v3Cnt <maxVal && v4Cnt==maxVal)
        return VORONOI_OVER_VERTEX_4;

    if (v1Cnt==maxVal && v2Cnt==maxVal && v3Cnt <maxVal && v4Cnt <maxVal)
        return VORONOI_OVER_EDGE_1_2;

    if (v1Cnt <maxVal && v2Cnt==maxVal && v3Cnt==maxVal && v4Cnt <maxVal)
        return VORONOI_OVER_EDGE_2_3;

    if (v1Cnt==maxVal && v2Cnt <maxVal && v3Cnt==maxVal && v4Cnt <maxVal)
        return VORONOI_OVER_EDGE_3_1;

    if (v1Cnt==maxVal && v2Cnt <maxVal && v3Cnt <maxVal && v4Cnt==maxVal)
        return VORONOI_OVER_EDGE_1_4;

    if (v1Cnt <maxVal && v2Cnt==maxVal && v3Cnt <maxVal && v4Cnt==maxVal)
        return VORONOI_OVER_EDGE_2_4;

    if (v1Cnt <maxVal && v2Cnt <maxVal && v3Cnt==maxVal && v4Cnt==maxVal)
        return VORONOI_OVER_EDGE_3_4;

    if (v1Cnt==maxVal && v2Cnt==maxVal && v3Cnt==maxVal && v4Cnt <maxVal)
        return VORONOI_OVER_TRIANGLE_1_3_2;

    if (v1Cnt==maxVal && v2Cnt==maxVal && v3Cnt <maxVal && v4Cnt==maxVal)
        return VORONOI_OVER_TRIANGLE_1_2_4;

    if (v1Cnt==maxVal && v2Cnt <maxVal && v3Cnt==maxVal && v4Cnt==maxVal)
        return VORONOI_OVER_TRIANGLE_1_4_3;

    if (v1Cnt <maxVal && v2Cnt==maxVal && v3Cnt==maxVal && v4Cnt==maxVal)
        return VORONOI_OVER_TRIANGLE_2_3_4;

    return NONE;
}


}// namespace Makena

#endif /*_MAKENA_VORONOI_3SIMPLEX_GRAPH_HPP_*/

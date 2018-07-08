#ifndef _MAKENA_CONTACT_UPDATER_FACE_EDGE_HPP_
#define _MAKENA_CONTACT_UPDATER_FACE_EDGE_HPP_

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
#include "contact_pair_info.hpp"
#include "loggable.hpp"
#include "intersection_convex_polygon_2d.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file contact_updater_face_edge.hpp
 *
 * @brief assuming the current contact feature pair are face-edge or edge-face,
 *        check if it needs to be updated to another feature pair
 *        based on the current face normal directions and intersection of
 *        the two faces.
 *       
 *  Legend:
 *
 *          //{eps} : Parallel within angular tolerance eps
 *
 *         _|_{eps} : Perpendicular within angular tolerance eps
 *
 *          Nx : Normal vector. Usually a face normal.
 *
 *          Dx : Direction vector. Usually half edge direction from src to dir.
 *
 *          ·  : dot product
 *
 *          x  : cross product
 *
 *         |V| : Vector norm 2
 *
 *
 *  FACE1-FACE2
 *      
 *      Let Nf1 be the normal of FACE1
 *      Let Df2 be the unit direction vector of EDGE2
 *
 *      - If Nf1 _|_{2.0eps} Df2 => FACE1-EDGE2 (NO CHANGE)
 *
 *      - Else
 *
 *        Let D12 be the unit tilting direction vector of EDGE2, i.e.,
 *            D12 =   Df2 if Nf1 · Df2 < 0.0
 *            D12 = - Df2 otherwise
 *
 *        Let INT be the 2D convex intersection line segment of FACE1 and EDGE2
 *        on the XY-plane perpendicular to Nf1.
 *        Find the extremal vertex Pint of INT toward D12. 
 *        
 *
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
 *
 *          - If Pint = IT_VERTEX_INTERIOR (Vint1, EDGE2)
 *                                             => Vint1-FACE2 (VERTEX-EDGE)
 *                 
 *                       .....\
 *                       ......\ Vint1
 *                       o------*-------o   D12 ==>
 *                      EDGE2../
 *                       ...../
 *                       ..../
 *
 *          - If Pint = IT_INTERIOR_VERTEX (FACE1, Vint2)
 *                                             => FACE1-Vint2 (FACE-VERTEX)
 *
 *                     ..............|
 *                     .......Vint2..|
 *                     ..o------o....|  D12 ==>
 *                     ...EDGE2......|
 *                     ..............|
 *                     ..............|
 *
 *
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
 *
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
 *
 *
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
 *
 *
 */

namespace Makena {

using namespace std;

using IntSec2D = IntersectionFinderConvexPolygon2D;

class ContactUpdater_FACE_EDGE : public Loggable {

  public:

    ContactUpdater_FACE_EDGE(
        ConvexRigidBody&         body1,
        ConvexRigidBody&         body2,
        ContactPairInfo&         info,
        const double&            epsilonZero,
        const double&            epsilonAngle,
        std::ostream&            logStream
    );

    ~ContactUpdater_FACE_EDGE();

    /** @brief main function
     *
     *  @return true  : Need to run GJK to find the contact pair.
     *          false : The latest (updated) feature pair is valid.
     */
    bool update();

    void orientWorld();

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    void checkForTilt();

    Mat3x3 rotMatAlignZDirToZAxis(const Vec3& zDir);

    void rotateWorld(const Vec3& zDir);

    bool findIntersection2D();

    long findExtremeFeature2D();

    void dispatchAndUpdate(const long index);

    long edgeIndex(IntSec2D::OutputElem& e);

    bool isVertexIncidentToEdge(VertexIt v, HalfEdgeIt h);

    bool findHalfEdgeBody1(
                  const long  index1, const long  index2, HalfEdgeIt& heOut);


    void processIntsec_EDGE_EDGE       (const long index);
    void processIntsec_EDGE_VERTEX     (const long index);
    void processIntsec_VERTEX_EDGE     (const long index);
    void processIntsec_VERTEX_VERTEX   (const long index);
    void processIntsec_VERTEX_INTERIOR (const long index);
    void processIntsec_INTERIOR_VERTEX (const long index);

    void checkAndUpdateEdgeCases();
    void updateNonTilting_VERTEX_VERTEX();
    void updateNonTilting_VERTEX_EDGE();
    void updateNonTilting_EDGE_VERTEX();
    void updateNonTilting_EDGE_VERTEX_EDGE_VERTEX();
    void updateNonTilting_EDGE_VERTEX_VERTEX_EDGE();
    void updateNonTilting_EDGE_VERTEX_VERTEX_VERTEX();
    void updateNonTilting_VERTEX_EDGE_EDGE_VERTEX();
    void updateNonTilting_VERTEX_EDGE_VERTEX_EDGE();
    void updateNonTilting_VERTEX_EDGE_VERTEX_VERTEX();
    void updateNonTilting_VERTEX_VERTEX_EDGE_VERTEX();
    void updateNonTilting_VERTEX_VERTEX_VERTEX_EDGE();
    void updateNonTilting_VERTEX_VERTEX_VERTEX_VERTEX();

    ConvexRigidBody&             mBody1;
    ConvexRigidBody&             mBody2;
    ContactPairInfo&             mInfo;
    const double                 mEpsilonZero;
    const double                 mEpsilonAngle;
    bool                         mFaceOnBody1;
    bool                         mTiltedTowardVertex1;
    bool                         mTiltedTowardVertex2;
    Vec2                         mTiltingDir;
    VertexIt                     mTiltingVertex;
    VertexIt                     mVit1;
    VertexIt                     mVit2;

    Mat3x3                       mRotMat1;
    Mat3x3                       mRotMat2;
    Vec3                         mCom1;
    Vec3                         mCom2;

    vector<VertexIt>             mVertices;
    vector<HalfEdgeIt>           mHalfEdges;
    long                         mLastIndex;
    vector<IntSec2D::OutputElem> mIntsec;


};


}// namespace Makena

#endif /*_MAKENA_GJK_CONTACT_UPDATER_FACE_EDGE_HPP_*/

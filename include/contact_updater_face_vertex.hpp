#ifndef _MAKENA_CONTACT_UPDATER_FACE_VERTEX_HPP_
#define _MAKENA_CONTACT_UPDATER_FACE_VERTEX_HPP_

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
 * @file contact_updater_face_vertex.hpp
 *
 * @brief assuming the current contact feature pair are face-vertex or 
 *        vertex-face, check if it needs to be updated to another feature pair
 *        based on the current face normal directions and intersection of
 *        the face and the vertex.
 *       
 *        Let INT be the 2D convex intersection line segment of FACE1 and 
 *        VERTEX2 on the XY-plane perpendicular to the face normal.
 *       
 *
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
 *
 *          - If Pint = IT_VERTEX_VERTEX (Vint1, Vint2)
 *                                             => Vint1-Vint2 (VERTEX-VERTEX)
 *                 
 *                       .....\
 *                       ......\ Vint1==Vint2
 *                       ......*o
 *                      ......./
 *                       ...../
 *                       ..../
 *
 *          - If Pint = IT_INTERIOR_VERTEX (FACE1, Vint2)
 *                                             => FACE1-Vint2 (FACE-VERTEX)
 *
 *                     ..............|
 *                     .......Vint2..|
 *                     .........o....|
 *                     ..............|
 *                     ..............|
 *                     ..............|
 */
namespace Makena {

using namespace std;

using IntSec2D = IntersectionFinderConvexPolygon2D;

class ContactUpdater_FACE_VERTEX : public Loggable {

  public:

    ContactUpdater_FACE_VERTEX(
        ConvexRigidBody&         body1,
        ConvexRigidBody&         body2,
        ContactPairInfo&         info,
        const double&            epsilonZero,
        const double&            epsilonAngle,
        std::ostream&            logStream
    );

    ~ContactUpdater_FACE_VERTEX();

    /** @brief main function
     *
     *  @return true  : Need to run GJK to find the contact pair.
     *          false : The latest (updated) feature pair is valid.
     */
    bool update();

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    Mat3x3 rotMatAlignZDirToZAxis(const Vec3& zDir);

    void rotateWorld(const Vec3& zDir);

    bool findIntersection2D();

    long findExtremeFeature2D();

    void dispatchAndUpdate();

    long edgeIndex(IntSec2D::OutputElem& e);

    bool findHalfEdgeBody1(
                  const long  index1, const long  index2, HalfEdgeIt& heOut);

    void processIntsec_EDGE_VERTEX();
    void processIntsec_VERTEX_VERTEX();

    ConvexRigidBody&             mBody1;
    ConvexRigidBody&             mBody2;
    ContactPairInfo&             mInfo;
    const double                 mEpsilonZero;
    const double                 mEpsilonAngle;
    bool                         mFaceOnBody1;

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

#endif /*_MAKENA_GJK_CONTACT_UPDATER_FACE_VERTEX_HPP_*/

#ifndef _MAKENA_CONTACT_UPDATER_EDGE_EDGE_HPP_
#define _MAKENA_CONTACT_UPDATER_EDGE_EDGE_HPP_

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
 * @file contact_updater_edge_edge.hpp
 *
 * @brief assuming the current contact feature pair are face-edge or edge-face,
 *        check if it needs to be updated to another feature pair
 *        based on the current face normal directions and intersection of
 *        the two faces.
 *       
 *        If the old edges are parallel and now it's not, need contact update.
 *
 *  EDGE1-EDGE2
 *      
 *        - there is one intersectino feature.
 *
 *          - If Pint = IT_EDGE_EDGE (Eint1,Eint2) => Eint1-Eint2 (EDGE-EDGE)
 *
 *                                *
 *                                |EDGE1
 *                                |
 *                        o-------+-----o
 *                         EDGE2  |
 *                                |
 *                                *
 *
 *          - If Pint = IT_EDGE_VERTEX (Eint1,Vint2)
 *                                              => Eint1-Vint2 (EDGE-VERTEX)
 *                             Vint2
 *                        *-----o-----*
 *                       EDGE1   \
 *                                \EDGE2
 *                                 \
 *                                  o
 *
 *
 *          - If Pint = IT_VERTEX_EDGE (Vint1, Eint2)
 *                                              => Vint1-Eint2 (VERTEX-EDGE)
 *                             Vint1
 *                        o-----*-----o
 *                       EDGE2   \
 *                                \EDGE1
 *                                 \
 *                                  *
 *
 *
 *          - If Pint = IT_VERTEX_VERTEX (Vint1, Vint2)
 *                                              => Vint1-Vint2 (VERTEX-VERTEX)
 *                                   Vint1=Vint2
 *                        o-----------o*
 *                       EDGE2          \
 *                                       \EDGE1
 *                                        \
 *                                         *
 *
 *        - If there are two intersectino features.
 *
 *          - If Pint1 = IT_EDGE_VERTEX (Eint11,Vint12)
 *            - If Pint1 = IT_EDGE_VERTEX (Eint21,Vint22)
 *
 *                          Vint12    Vint22
 *                        *-----o========o-----*
 *                       EDGE1     EDGE2
 *
 *
 *            - If Pint1 = IT_VERTEX_EDGE (Vint21,Eint22)
 *
 *                          Vint12    Vint21
 *                        *-----o========*-----o
 *                       EDGE1        EDGE2
 *
 *            - If Pint1 = IT_VERTEX_VERTEX (Vint21,Vint22)
 *
 *                          Vint12    Vint21=Vint22
 *                        *-----o========*o
 *                       EDGE1     EDGE2
 *
 *
 *          - If Pint1 = IT_VERTEX_EDGE (Vint11,Eint12)
 *            - If Pint1 = IT_EDGE_VERTEX (Eint21,Vint22)
 *
 *                          Vint11    Vint22
 *                        o-----*========o-----*
 *                       EDGE2     EDGE1
 *
 *            - If Pint1 = IT_VERTEX_EDGE (Vint21,Eint22)
 *
 *                          Vint11    Vint21
 *                        o-----*========*-----o
 *                       EDGE2     EDGE1
 *
 *            - If Pint1 = IT_VERTEX_VERTEX (Vint21,Vint22)
 *
 *                          Vint11    Vint21=Vint22
 *                        o-----*========*o
 *                       EDGE2     EDGE1
 *
 *          - If Pint1 = IT_VERTEX_VERTEX (Vint11,Vint12)
 *            - If Pint1 = IT_EDGE_VERTEX (Eint21,Vint22)
 *
 *                   Vint11=Vint12    Vint22
 *                        o*============o-----*
 *                           EDGE2         EDGE1
 *
 *
 *            - If Pint1 = IT_VERTEX_EDGE (Vint21,Eint22)
 *
 *                    Vint11=VInt12    Vint21
 *                        o*==============*-----o
 *                          EDGE1            EDGE2
 *
 *            - If Pint1 = IT_VERTEX_VERTEX (Vint21,Vint22)
 *
 *                   Vint11=Vint12  Vint21=Vint22
 *                        o*=============*o
 *                           EDGE1/EDGE2
 *
 *
 *
 *
 */

namespace Makena {

using namespace std;

using IntSec2D = IntersectionFinderConvexPolygon2D;

class ContactUpdater_EDGE_EDGE : public Loggable {

  public:

    ContactUpdater_EDGE_EDGE(
        ConvexRigidBody&         body1,
        ConvexRigidBody&         body2,
        ContactPairInfo&         info,
        const double&            epsilonZero,
        const double&            epsilonAngle,
        std::ostream&            logStream
    );

    ~ContactUpdater_EDGE_EDGE();

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

    bool update_CROSS();
    bool update_PARALLEL();

    Mat3x3 rotMatAlignZDirToZAxis(const Vec3& zDir);

    void rotateWorld(const Vec3& zDir);
    void rotateWorld(const Vec3& zDir, const Vec3& xDir);

    bool findIntersection2D();

    void dispatchAndUpdate();

    void processIntsec_EDGE_VERTEX();
    void processIntsec_VERTEX_EDGE();
    void processIntsec_VERTEX_VERTEX();

    ConvexRigidBody&             mBody1;
    ConvexRigidBody&             mBody2;
    ContactPairInfo&             mInfo;
    const double                 mEpsilonZero;
    const double                 mEpsilonAngle;

    Mat3x3                       mRotMat1;
    Mat3x3                       mRotMat2;
    Vec3                         mCom1;
    Vec3                         mCom2;

    vector<IntSec2D::OutputElem> mIntsec;

    VertexIt                     mVit11;
    VertexIt                     mVit12;
    VertexIt                     mVit21;
    VertexIt                     mVit22;
};


}// namespace Makena

#endif /*_MAKENA_GJK_CONTACT_UPDATER_EDGE_EDGE_HPP_*/

#ifndef _MAKENA_CONTACT_UPDATER_FURTHER_CHECKER_HPP_
#define _MAKENA_CONTACT_UPDATER_FURTHER_CHECKER_HPP_

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
 * @file contact_updater_further_checker.hpp
 *
 * @brief it checks if the neighbor features of the current active features are
 *        penetrating.
 *        If the detected penetration is simple and shallow, it updates the
 *        current active features accordingly.
 *        If the detected penetration is deep or complex, it signals a need for
 *        a proper collision detector.
 */
namespace Makena {

class ContactUpdaterFurtherChecker : public Loggable {

  public:

    ContactUpdaterFurtherChecker(
        ConvexRigidBody&         body1,
        ConvexRigidBody&         body2,
        ContactPairInfo&         info,
        const double&            epsilonZero,
        const double&            epsilonAngle,
        std::ostream&            logStream
    );

    ~ContactUpdaterFurtherChecker();

    /** @brief checks for penetration and updates the current active feature
     *         pair if necessary.
     *
     *  @return true : need to run a proper collision detector as the detected
     *                 penetration is too deep or complex.
     *
     *  @return false: no need to run a collision detector.
     */
    bool checkFeatures();

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    enum WEDGE_CLASSIFIER {
        WC_NONE,

        WC_ON_APEX,
        WC_ON_EDGE1,
        WC_ON_EDGE2,
        WC_INSIDE,
        WC_OUTSIDE,
        WC_INTSEC_ON_APEX,
        WC_INTSEC_ON_EDGE1,
        WC_INTSEC_ON_EDGE2,
        WC_INTSEC_ON_EDGE12,

        WC_PARSE_RESULT_BEGIN,

        WC_COMPLETELY_INSIDE,
        WC_EDGE1_INSIDE_VERTICES,
        WC_EDGE1_TOUCHING,
        WC_EDGE2_INSIDE_VERTICES,
        WC_EDGE2_TOUCHING,
        WC_APEX_ONE_EDGE_CUT_ACROSS,
        WC_APEX_ONE_EDGE_TOUCHING,
        WC_APEX_INSIDE_VERTICES,
        WC_APEX_COINCIDENT_VERTEX,

        WC_END
    };

    bool checkFeatures_FACE_EDGE();
    bool checkFeatures_EDGE_FACE();
    bool checkFeatures_FACE_VERTEX();
    bool checkFeatures_VERTEX_FACE();
    bool checkFeatures_EDGE_EDGE();
    bool checkFeatures_EDGE_EDGE_CROSSING();
    bool checkFeatures_EDGE_EDGE_PARALLEL();
    bool checkFeatures_EDGE_VERTEX();
    bool checkFeatures_VERTEX_EDGE();
    bool checkFeatures_VERTEX_VERTEX();

    bool isCCW(
        const Vec3& v1,
        const Vec3& v2,
        const Vec3& v3,
        const Vec3& axis
    );

    bool isRayInsideCone(
        const Vec3& dir,
        VertexIt    vitApex,
        bool        coneOnBody1
    );

    void rotateWorld(
        const Vec3& yDir,
        const Vec3& zDir,
        Mat3x3&     rotMat1,
        Vec3&       com1,
        Mat3x3&     rotMat2,
        Vec3&       com2
    );

    void rotateWorld(
        const Vec3& zDir,
        Mat3x3&     rotMat1,
        Vec3&       com1,
        Mat3x3&     rotMat2,
        Vec3&       com2
    );

    Mat3x3 rotMatAlignZDirToZAxis(const Vec3& zDir);

    bool findIntersection(
        const Vec2& p11,
        const Vec2& p12,
        const Vec2& p21,
        const Vec2& p22,
        Vec2&       intsec
    );

    bool classifyPointsAgainstWedge(
        vector<Vec2>&                  fan,
        const Vec2&                    apex,
        const Vec2&                    side1,
        const Vec2&                    side1perp,
        const Vec2&                    side2,
        const Vec2&                    side2perp,
        vector<enum WEDGE_CLASSIFIER>& predMain,
        vector<enum WEDGE_CLASSIFIER>& predAux
    );

    enum WEDGE_CLASSIFIER parseClassifiersPointsAgainstWedge(
        vector<enum WEDGE_CLASSIFIER>& predMain,
        vector<enum WEDGE_CLASSIFIER>& predAux
    );

    bool updateFeatures_EDGE_VERTEX(
        const enum WEDGE_CLASSIFIER    parseResult,
        vector<enum WEDGE_CLASSIFIER>& predMain,
        vector<enum WEDGE_CLASSIFIER>& predAux,
        vector<HalfEdgeIt>&            hes,
        const Vec2&                    apex,
        vector<Vec2>&                  fan,
        const Vec2&                    side1perp,
        const Vec2&                    side2perp,
        FaceIt                         fit1,
        FaceIt                         fit2,
        const bool                     wedgeOnBody1
    );

    void findProjectedAreaAndDistance(
        const Vec2& p11,
        const Vec2& p12,
        const Vec2& p21,
        const Vec2& p22,
        double&     area,
        double&     distance,
        double&     alignment
    );

    void findAppropriateFacePair(
        long                                                  indexBody1,
        long                                                  indexBody2,
        vector<IntersectionFinderConvexPolygon2D::InputElem>& fan1_2D,
        vector<IntersectionFinderConvexPolygon2D::InputElem>& fan2_2D,
        vector<HalfEdgeIt>&                                   body1_hes,
        vector<HalfEdgeIt>&                                   body2_hes
    );

    bool guessPenetrationDirectionFromIntersection(
        vector<IntersectionFinderConvexPolygon2D::OutputElem>& intsec,
        Vec2&                                                  dir
    );

    void guessPenetrationDirectionFromSpread(
        vector<IntersectionFinderConvexPolygon2D::InputElem>&  fan1_2D,
        vector<IntersectionFinderConvexPolygon2D::InputElem>&  fan2_2D,
        vector<IntersectionFinderConvexPolygon2D::OutputElem>& intsec_2D,
        Vec2&                                                  dir
    );

    ConvexRigidBody&             mBody1;
    ConvexRigidBody&             mBody2;
    ContactPairInfo&             mInfo;
    const double                 mEpsilonZero;
    const double                 mEpsilonAngle;
};


}// namespace Makena

#endif /*_MAKENA_CONTACT_UPDATER_FURTHER_CHECKER_HPP_*/

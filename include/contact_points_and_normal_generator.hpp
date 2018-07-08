#ifndef _MAKENA_CONTACT_POINTS_AND_NORMAL_GENERATOR_HPP_
#define _MAKENA_CONTACT_POINTS_AND_NORMAL_GENERATOR_HPP_

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
#include "intersection_convex_polygon_2d.hpp"
#include "loggable.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file contact_points_and_normal_generator.hpp
 *
 * @brief Generates a setof point pairs in LCS from the active feature pair 
 *        and finds a contact direction. They are used to generate unilateral 
 *        constraints.
 */
namespace Makena {

using namespace std;

class ContactPointsAndNormalGenerator : public Loggable {

  public:

    /** @brief constructor.
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
     *  @param epsilonZero  (in): numerical tolerance to consider the
     *                            distance value zero.
     *
     *  @param epsilonAngle (in): numerical tolerance to consider the
     *                            angle zero.
     *
     *  @param useGeomConfigTemp (in): true if you want to use the 
     *                            preliminary uncommitted geometric 
     *                            configuration
     *
     *  @param maxNumPointsPerContact (in): maximum number of point pairs
     *                            per contact pair. 
     *  @param logStream    (in): the output stream for logging.
     */
    ContactPointsAndNormalGenerator(
        ConvexRigidBody&  body1,
        ConvexRigidBody&  body2,
        ContactPairInfo&  info,
        vector<BDVertex>& Q,
        const double&     epsilonZero             = EPSILON_SQUARED,
        const double&     epsilonAngle            = EPSILON_SQUARED,
        bool              useGeomConfigTemp       = true,
        const long        maxNumPointsPerContact  = 6,
        std::ostream&     logStream               = std::cerr
    );


    /** @brief destructor. Nothing to be done.
     */
    ~ContactPointsAndNormalGenerator();

    /** @brief finds the active feature pair on Body 1 and Body2 from 
     *         the minimal boundary simplex of the binary dilation.
     */
    void findActiveFeaturesFromBinaryDilation();

    /** @brief checks and updates the active feature pair based on tne
     *         geometric configuration of the incident features.
     */
    void adjustFeatures();

    /** @brief generates the contact points pairs and the contact normal
     *         from the acitive feature pair 
     */
    void generateContatctPointsAndNormalFromActiveFeatures();


#ifndef UNIT_TESTS
  private:
#else
  public:
#endif

    void process_FACE_FACE();
    void process_FACE_FACE_tryOneDirection(const Vec3& normalDir);
    void process_FACE_FACE_projection();

    void process_FACE_EDGE();
    void process_FACE_EDGE_tryOneDirection(const Vec3& normalDir);
    void process_FACE_EDGE_projection();

    void process_FACE_VERTEX();

    void process_EDGE_FACE();
    void process_EDGE_FACE_tryOneDirection(const Vec3& normalDir);
    void process_EDGE_FACE_projection();

    void process_EDGE_EDGE();
    void process_EDGE_EDGE_PARALLEL(const Vec3& axisDir);
    Vec3 interpolate(const Vec3&  p1, const Vec3&  p2, const double z);

    void process_EDGE_EDGE_CROSSING(const Vec3& zDir);

    void process_EDGE_VERTEX();
    void process_VERTEX_FACE();
    void process_VERTEX_EDGE();
    void process_VERTEX_VERTEX();

    Vec3 findRelativeVelocityOfBody1RelativeToBody2();

    void findFeatureOfRigidBody(
        Manifold&                          ch,
        vector<VertexIt>&                  vits,
        FaceIt&                            fit,
        EdgeIt&                            eit,
        VertexIt&                          vit,
        enum ContactPairInfo::FeatureType& type
    );

    void findVerticesOfRigidBody(
        vector<VertexIt>& vits1,
        vector<VertexIt>& vits2
    );

    void adjustFeatures_FACE_EDGE();
    void adjustFeatures_EDGE_FACE();
    void adjustFeatures_FACE_VERTEX();
    void adjustFeatures_VERTEX_FACE();
    void adjustFeatures_EDGE_EDGE();
    void adjustFeatures_EDGE_EDGE_CROSSING();
    void adjustFeatures_EDGE_EDGE_PARALLEL();
    void adjustFeatures_EDGE_VERTEX();
    void adjustFeatures_VERTEX_EDGE();
    void adjustFeatures_VERTEX_VERTEX();

    bool isCCW(const Vec3& v1,const Vec3& v2,const Vec3& v3,const Vec3& axis);

    void findContactNormal();

    void rotateWorld(
        Mat3x3&      rotMat1,
        Vec3&        com1,
        Mat3x3&      rotMat2,
        Vec3&        com2
    );

    void rotateWorld(
        const Vec3& zDir,
        Mat3x3&      rotMat1,
        Vec3&        com1,
        Mat3x3&      rotMat2,
        Vec3&        com2
    );

    Vec3 findClosestPointOnEdge(
        const Vec3& p1,
        const Vec3& p2,
        const Vec3& b
    );

    void findClosestPointsBetweenTwoEdges(
        const Vec3& a1,
        const Vec3& a2,
        const Vec3& b1,
        const Vec3& b2,
        Vec3&       pA,
        Vec3&       pB
    );

    Vec3 findProjectedPointOnPlane(
        const Vec3& pBase,
        const Vec3& N,
        const Vec3& pTest
    );

    void from2DXYto3DPoints(
        vector<IntersectionFinderConvexPolygon2D::OutputElem>&
                        points2D,
        Mat3x3&         rotMat1,
        Vec3&           com1,
        Mat3x3&         rotMat2,
        Vec3&           com2
    );

    void orthoProj(
        vector<IntersectionFinderConvexPolygon2D::OutputElem>&
                      points2D,
        const Vec3&   N,
        const Vec3&   pBase,
        vector<Vec3>& points3D
    );

    vector<bool> reduceContactPoints(
        vector<IntersectionFinderConvexPolygon2D::OutputElem>& points2D
    );
       
    ConvexRigidBody&         mBody1;
    ConvexRigidBody&         mBody2;
    ContactPairInfo&         mInfo;
    vector<BDVertex>&        mQ;
    const double             mEpsilonZero;
    const double             mEpsilonAngle;
    const bool               mUseGeomConfigTemp;
    const long               mMaxNumPointsPerContact;
#ifdef UNIT_TESTS

    // Needed for test_visualizer_contact_points_and_normal_generator.cpp,
    // which keeps an object of this for a long time.
    void updateGeomConfigInternally();

    Mat3x3                   mQmat1;
    Mat3x3                   mQmat2;
    Vec3                     mCoM1;
    Vec3                     mCoM2;
#else
    const Mat3x3             mQmat1;
    const Mat3x3             mQmat2;
    const Vec3               mCoM1;
    const Vec3               mCoM2;
#endif

#ifdef UNIT_TESTS

    static string FeatureTypeString(ContactPairInfo::FeatureType t);

    // Rendering Helpers.

    void renderObj1ActiveFace(
        Vec3&         color,
        float         alpha,
        vector<Vec3>& vertices,
        vector<Vec3>& colors,
        vector<float>&alphas,
        vector<Vec3>& normals
    );

    void renderObj2ActiveFace(
        Vec3&         color,
        float         alpha,
        vector<Vec3>& vertices,
        vector<Vec3>& colors,
        vector<float>&alphas,
        vector<Vec3>& normals
    );


    void renderObj1ActiveEdge(
        Vec3&         color,
        vector<Vec3>& vertices,
        vector<Vec3>& colors
    );


    void renderObj2ActiveEdge(
        Vec3&         color,
        vector<Vec3>& vertices,
        vector<Vec3>& colors
    );


    void renderObj1ActiveVertex(
        Vec3&         color,
        vector<Vec3>& vertices,
        vector<Vec3>& colors
    );

    void renderObj2ActiveVertex(
        Vec3&         color,
        vector<Vec3>& vertices,
        vector<Vec3>& colors
    );


    void renderContactPoints(
        Vec3&         color1,
        Vec3&         color2,
        vector<Vec3>& vertices,
        vector<Vec3>& colors
    );


#endif

};


}// namespace Makena


#endif /*_MAKENA_CONTACT_POINTS_AND_NORMAL_GENERATOR_HPP_*/

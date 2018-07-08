#ifndef _MAKENA_CONTACT_DISCOVERER_HPP_
#define _MAKENA_CONTACT_DISCOVERER_HPP_

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <exception>

#include "primitives.hpp"
#include "manifold.hpp"
#include "convex_rigid_body.hpp"
#include "contact_pair_info.hpp"
#include "broad_phase_aabb_collision_detector.hpp"
#include "obb_obb_test.hpp"
#include "gjk_origin_finder.hpp"
#include "bd_boundary_simplex_finder.hpp"
#include "contact_points_and_normal_generator.hpp"
#include "contact_updater.hpp"
#include "jacobian_constraint.hpp"
#include "constraint_manager.hpp"
#include "intersection_finder.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file contact_discoverer.hpp
 *
 * @brief
 *
 */
namespace Makena {

/** @brief NEEDS REWRITING
 *         find the penetration orientation and contact feature pair
 *         of two colliding rigid bodies.
 *         It is assumed that two rigid bodies are known to be
 *         in contact or penetrating into each other.
 *         It aims at bringing the two bodies in a contact state
 *         at a pair of features at the next step, if the
 *         penetration is not deep.
 *           
 *         1. Find a set of separation axes by scaling OBBs down.
 *            There are two purposes: 
 *            - To assist finding the penetratino manifold effectively.
 *            - Rough indicator of contact direction.
 *
 *         2. Find the penetration manifold PM
 *            (intersection of two rigid bodies)
 *
 *         3. Run the distribution analysis on PM.
 *            - Relative velocity at each vertex.
 *            - Weighted normal at each face (weight is the area of the face)
 *            - Relative normalized displacement in GCS relative to CoM.
 *              Normalization is over the distrubution of the entire polytope.
 *
 *         3. Classify the PM (and spatial relation between two rigid bodies)
 *            into one of the following.
 *
 *            A). Touching at a point
 *                Condition: PM's dimension is 0.
 *
 *            B). Touching at an edge.
 *                Condition: PM's dimension is 1.
 *           
 *            C). Touching in a polygon.
 *                Condition: PM's dimension is 2.
 *
 *            D1). Polytope 1 is inscribed or properly inside Polytope 2.
 *             
 *            D2). Polytope 2 is inscribed or properly inside Polytope 1.
 *
 *            E). Proper Intersection.
 *
 *
 *            # From D) onwards, PM's dimension is 3.
 *
 *            D). Neither CoM is in PM. (Shallow penetration)
 *
 *            E1). CoM of Polytope 1 is in PM. (Deep penetration)
 *               Condition: Point inclusion tests of CoM1 into PM is positive.
 *
 *              E1a). Polytope 1 is properly included in the interior of 
 *                    Polytope 2.
 *               Condition: The prediates of all the vertices of PM are 
 *                          INTERIOR_X.
 *
 *              E1b). Polytope 2 intersects or inscribed toPolytope 2.
 *               Condition: Complement of E1a.
 *
 *            E2). CoM of Polytope 2 is in PM. (Deep penetration)
 *               Condition: Point inclusion tests of CoM2 into PM is positive.
 *
 *              E2a). Polytope 2 is properly included in the interior of 
 *                    Polytope 1.
 *               Condition: The prediates of all the vertices of PM are 
 *                          X_INTERIOR.
 *
 *              E2b). Polytope 2 intersects or inscribed toPolytope 1.
 *               Condition: Complement of E2a.
 *
 *            F). Both CoMs are in PM. (Deep penetration)
 *               Condition: Point inclusion tests of CoM1 and CoM2 into PM is 
 *                          positive.
 *                
 *              F1). Polytope 1 and 2 intersect with each other, or one is 
 *                   inscribed to the other.
 *                  
 *               Condition: Complement of (F2 || F3).
 *
 *              F2). Polytope 1 is properly included in the interior of 
 *                   Polytope 2.
 *               Condition: The prediates of all the vertices of PM are 
 *                          INTERIOR_X.
 *
 *              F3). Polytope 2 is properly included in the interior of 
 *                   Polytope 1.
 *               Condition: The prediates of all the vertices of PM are 
 *                          X_INTERIOR.
 *
 *         4. If the classification is one of, D, E1a, E2a, and E1.
 *            then perform relative velocity analysis on PM.
 *            - Find the relative velocity of polytope 1 to polytope 2 each 
 *              vertex in PM.
 *
 *              relvel_12 = (vcom1 - a1_i X ω1) - (vcom2 - a2_i X ω2)
 *              * vcom1 and vcom2 are the linear velocities at CoMs 
 *                respectively.
 *              * ω1 and ω2 are the angular velocities of polytope 1 and 2 
 *                respectively.
 *              * a1_i is the point on polytope 1 in LCS that corresponds to 
 *                vertex i.
 *              * a2_i is the point on polytope 2 in LCS that corresponds to 
 *                vertex i.
 *              
 *            - Perform PCA.
 *              If PCA is 0 or 1 dimensional, the velocity is uniform.
 *              If PCA is 2 dimensional, the rotation axis is found to be
 *              the 3rd eigen vector. 
 *              Let dirPCA be a unit vector in the 3rd eigen vector.
 *              If the rotation axis is aligned with consider it 1 dimensional.
 *              IF the resultant dimention <= 1, then check if the dot product
 *              of the mean and the separation axis is positive.
 *             
 *
 *         5. Dispatch process based on the classification.
 *
 *         5-A. It can be one of the following contact feature pairs based on 
 *             the vertex attribute.
 *             The penetration direction can be determined from the feature 
 *             normals.
 * 
 *              Type                Penetratio direction
 *              -----------------------------------------
 *              IF_VERTEX_VERTEX    NT_VERTEX_VERTEX_AVG
 *              IF_VERTEX_EDGE      NT_EDGE_NORMAL_2
 *              IF_VERTEX_FACE      NT_FACE_NORMAL_2
 *              IF_EDGE_VERTEX      NT_EDGE_NORMAL_1
 *              IF_FACE_VERTEX      NT_FACE_NORMAL_1
 *              IF_EDGE_EDGE        NT_EDGE_CROSS_EDGE,
 *
 *         5-B. It can be one of the following contact feature pairs based on 
 *             the edge attribute
 *              Type                Penetratio direction
 *              -----------------------------------------
 *              IF_EDGE_EDGE        NT_EDGE_EDGE_AVG
 *              IF_EDGE_FACE        NT_FACE_NORMAL_2
 *              IF_FACE_EDGE        NT_FACE_NORMAL_1
 *              IF_FACE_FACE        NT_FACE_NORMAL_1
 *
 *         5-C. It is the following contact feature pair based on 
 *             the face attribute
 *              Type                Penetratio direction
 *              -----------------------------------------
 *              IF_FACE_FACE        NT_FACE_NORMAL_1
 *
 *         5-D). Shallow penetration.
 *
 *            5-D1). If the relative velocity is uniform, and the mean of 
 *                   relative velocity is within the the separating axes cone,
 *                   then use the mean as the collision direction.
 *                   Otherwise, we try to find the collision direction from
 *                   the shape of PM.
 *                   For each face of PM, we check if the normal is within the
 *                   separation axes cone.
 *                   Then we take the weighted average of normals for each of 
 *                   two polytopes, where the weight correspoinds to the area 
 *                   of the face.
 *                   Then, use those two averages as the collision derection.
 *
 *            5-D2). If one polytope is inscribed to, or properly in the 
 *                   interior of the other polytope, then we find the relative
 *                   velocity at CoM of interior polytope, and use it as the 
 *                   penetration direction. To find the feature on the other 
 *                   polytope, we perform a variant of GJK pivoting to find
 *                   the intersecting feature.
 *                  
 *
 *
 *         2. If PM is a point, that corresponds to
 *            a pair of vertices on the rigid bodies.
 *            make a unilateral constraint with it.
 *
 *            If PM is an edge, that corresponds to
 *            a pair of edges on the rigid bodies.
 *            make a unilateral constraint with it. *
 *
 *            If PM is a convex polygon, the corresponds
 *            to a pair of faces on the rigid bodies.
 *            make a unilateral constraint with it.
 *
 *            If PM is a convex polytope, then perform the following tests
 *            - Find the relative velocity of each vertex in PM
 *              and perform PCA analysis.
 *              If PCA is 0 or 1 dimensional, the velocity is uniform.
 *              If PCA is 2 dimensional, the rotation axis is found to be
 *              the 3rd eigen vector. If the rotation axis is aligned with
 *              the mean relative velocity consider it 1 dimensional.
 *              Otherwise, consider it 2 dimensional.
 *              (If PCA is 3 dimensional, something weird is happening.)
 *            - Tests if CoMs are included in PM.
 *
 *         3a. If the dimension of the relative velocities is 0- or
 *            1-dimensional, then use the mean relative velocity as the
 *            penetration direction.
 *            And we can find the contact feature pair using the penetration
 *            direction on PM.
 *
 *
 *         3b. If the dimension of the relative velocity is 1- or
 *             2-dimensional, then:
 *                Go to 4a if neither CoM is included in PM.
 *                Go to 4b if CoM of Polytope 1 is included in PM.
 *                Go to 4c if CoM of Polytope 2 is included in PM.
 *                Go to 4d if both CoMs are included in PM.
 *
 *         4a. Conceptually remove the degenerate vertices, edges, and the
 *            faces whose normals are not facing the opposing CoM from PM.
 *            After this clean up, it can be assumed that PM is decomposed
 *            into three:
 *            1. Surface of body 1, which is a compact bounded 2-manifold. (S1)
 *            2. Boundary of body 1 and 2, which is a compact 2-manifold
 *               partially bounded. Topologically is it is a square [0,1]x[0,1]
 *               with its top and bottom sides identified by the relation
 *               (x,0)=(x,1).
 *            3. Surface of body 2, which is a compact bounded 2-manifold. (S2)
 *
 *            In other words, we can assume the following.
 *            S1 and S2 are connected and facing the CoM of body 2 and body 1
 *            respectively.
 *            Find the weighted average of face normals N1 and N2 of S1 and S2
 *            respectively. The weight is the area of the face.
 *            We can assume such average normals of body 1 and 2 are opposing
 *            against each other, i.e., the dot product is negative.
 *            We use those normals to find the contact feature pair.
 *
 *            The rationale behind this is that we know the following.
 *            - The relative velocity varies on S1 and S2.
 *            - The penetration is not deep.
 *            These imply two objects are colliding on a wide surface and PM
 *            has a pancake like shape. It is quite likely to have a dominant
 *            face on either S1 and S2, which can be detected using the average
 *            normals.
 *            If there is no dominant face, then it is likely that we can
 *            detect a dominant edge of an extreme vertex on each of S1 and S2
 *            using such average normals.
 *
 *         4b. It is likely that Polytope 1 is smaller and mobile while
 *            Polytope 2 is larger and stationary. This is a deep and wild
 *            collision. We try to solve the penetration rather than focusing
 *            on finding a proper contact pair for the next step.
 *            For each face of PM that corresponds to a face of Polytope 2,
 *            we find a score whose formular is given as follows:
 *
 *            score = area / translation
 *
 *            Area is the area of the face and translation is the distance
 *            from a point of the face to the extreme point on Polytope 1
 *            along the negative direction of the normal of the face, which
 *            can be considered to be the distance required to solve the
 *            penetration.
 *            We use the normal as the penetration direction, and the contact
 *            pair is the feature incident to the extreme point on Polytope 1
 *            and the face on Polytope 2.
 *
 *         4c. Similar to 4b but reverse the role of Polytope 1 and 2.
 *
 *         4d. The wildest collision that rarely occurs.
 *            This is likely to be a fast and wild highspeed collision between
 *            two objects. We focus on solving the penetration rather than
 *            on finding a proper contact pair for the next step.
 *
 *            We run the same face scoring as 4b on all the faces of PM and
 *            the contact pair is found in an similar way.
 */
class ContactDiscoverer :public Loggable {

  public:

    /** @brief constructor
     *
     *  @param logStream(in): Log output
     */
    ContactDiscoverer(
        ConvexRigidBody& body1,
        ConvexRigidBody& body2,
        const double     epsilonZero,
        const double     epsilonZeroPCA,
        const double     epsilonAngle,
        const double     obbScalingStep,
        const long       obbScalingNumExtra,
        const bool       useTempGeomConfig,
        std::ostream&    logStream
    );

    /** @brief descructor. */
    ~ContactDiscoverer();

    ContactPairInfo discover();

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    void     findSeparationAxesFromOBBscaling();
    void     processPM_0D();
    void     processPM_1D();
    void     processPM_2D();
    void     isOnePolytopeInteriorOfTheOther(
        bool& polytope1Inside,
        bool& polytope2Inside
    );
    VertexIt findExtremePointOnPMInTheDirection(const Vec3& dir);
    Mat3x3   rotMatAlignZDirToZAxis(const Vec3& zDir);
    bool     doesFaceIncludeOriginXY(
        FaceIt        fit,
        const Mat3x3& Q,
        const Vec3&   com
    );
    void     processPM_Poly1Inside();
    FaceIt   findIntersectingFaceOnBody2FromPointInDirection(
        const Vec3& pBase,
        const Vec3& vDir
    );
    void     processPM_Poly2Inside();
    FaceIt   findIntersectingFaceOnBody1FromPointInDirection(
        const Vec3& pBase,
        const Vec3& vDir
    );
    void     processPM_3D();
    void     findRelativeVelocitiesOf1To2InPenetrationManifold(
        Vec3& mean,
        long& dimension,
        Vec3& primaryEigenVector
    );
    bool     isWithinSepAxesCone1To2(const Vec3& n);
    bool     isWithinSepAxesCone2To1(const Vec3& n);

    Vec3     findDirectionFromWeightedAverageSide1();
    Vec3     findDirectionFromWeightedAverageSide2();

    VertexIt findExtremeVertexSide1(const Vec3& dir);
    VertexIt findExtremeVertexSide2(const Vec3& dir);

    void     adjustFeatureOnPMSide1(
        const Vec3& dir,
        VertexIt    vitBest,
        EdgeIt&     eitBest,
        FaceIt&     fitBest
    );

    void     adjustFeatureOnPMSide2(
        const Vec3& dir,
        VertexIt    vitBest,
        EdgeIt&     eitBest,
        FaceIt&     fitBest
    );

    void adjustFeatureAndSetContactInfoSide1(
        VertexIt    vit,
        const Vec3& dir
    );

    void adjustFeatureAndSetContactInfoSide2(
        VertexIt    vit,
        const Vec3& dir
    );

    void updateContactNormal();

    ConvexRigidBody&         mBody1;
    ConvexRigidBody&         mBody2;
    ContactPairInfo          mInfo;
    const bool               mUseTempGeomConfig;
    const Mat3x3             mQmat1;
    const Vec3               mCoM1;
    const Mat3x3             mQmat2;
    const Vec3               mCoM2;
    const Vec3               mVlin1;
    const Vec3               mVang1;
    const Vec3               mVlin2;
    const Vec3               mVang2;
    const double             mEpsilonZero;
    const double             mEpsilonZeroPCA;
    const double             mEpsilonAngle;
    const double             mObbScalingStep;
    const long               mObbScalingNumExtra;
    IntersectionFinder       mISFinder;
    vector<Vec3>             mSeparationAxes1To2;

};



}// namespace Makena


#endif /*_MAKENA_CONTACT_DISCOVERER_HPP_*/

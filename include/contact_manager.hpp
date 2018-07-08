#ifndef _MAKENA_CONTACT_MANAGER_HPP_
#define _MAKENA_CONTACT_MANAGER_HPP_

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
#include "contact_discoverer.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file contact_manager.hpp
 *
 * @brief top-level class to manage the contacts among convex rigid bodies.
 *        The object of this class is meant to be a perpetual singleton.
 *
 *        Roughly speaking, it takes the geometric configuration of the
 *        rigid bodies, and produces unilateral and bilateral boxed Jacobian 
 *        constraints.
 *        The unilateral constraints are used for penetration avoidance, and
 *        the boxed bilateral constraints are used to approximate static or
 *        dynamic friction.
 *
 *        [Current Geometric Config at time step t]
 *
 *               ||
 *              \||/
 *
 *        Update Contact Points and the Normal and generate Jacobian
 *        constraints for contacts and frictions.
 *
 *               ||
 *              \||/
 *
 *
 *        [Temporarily Updated Geometric Config for time step t + Δt]
 *
 *               ||
 *              \||/
 *
 *        Discover new Contact Points and the Normal and generate Jacobian
 *        constraints for new contacts.
 *
 *               ||
 *              \||/
 *
 *        [Officially Updated Geometric Config for time step t + Δt]
 *
 *               ||
 *              \||/
 *
 *        If the contact constraints are all inactive, remove the contacts.
 *        If the contact constraints are active, 
 *        estimate the contact pressure to estimate the maximum friction force
 *        for the next step.
 *
 *
 *        *  It detects collisions based on the temporary geometric
 *           configuration for the current time step.
 *
 *        *  It generates information required to generate contact constraints
 *           and friction constraints. It stores them into ContactPairInfo in
 *           JacobianConstraint.
 *          
 *        *  It keeps track of the contacts and their constraints in 
 *           ContactPairInfo.
 *
 *        *  It removes ContactPairInfo if is is no longer needed.
 *
 *        [Contact Discovery]
 *        The collisions are discovered in the following 6 steps.
 *
 *        1. Broad-phase AABB collision pairs culling.
 *        2. OBB separation axis test.
 *        3. GJK test.
 *        4. Find intersection of two penetrating polytopes and find a pair of
 *           active features.
 *        5. Find pairs of contact points and the contact normal
 *        6. Generate unilateral constraints
 *
 *        At step 3, GJK test finds a simplex of the binary dilation that 
 *        encloses the origin if there is a collision.
 *        The origin in the binary dilation indicates the intersection of 
 *        two convex hulls.
 *
 *        In order to generate appropriate unilateral constraints to avoid
 *        collision, first, we woud like to find an appropriate pair of
 *        colliding features of the two colliding rigid bodies from which we
 *        will generate contact points pair and contact normal.
 *        Here 'feature' means any of vertex, edge or a face of the convex
 *        manifold.
 *        An appropriate feature pair sits on the boundary of the intersection
 *        of two polytopes that depends on the exit direction
 *        (negative of penetration  direction).
 *       
 *        We try to use se the negative of the relative velocity of two 
 *        objects at the feature pair as the exit direction.
 *        We first check if the relative velocity is consistent within the
 *        intersection. We calcualte the covariance of the relative velocity
 *        in the penetration intersection, and then perform PCA to check
 *        for consistency.
 *        If the relative velocity varies a lot, then we can't use it.
 *        We use the shape of the intersection and the separation axis
 *        found on the scaled (shrunk) polytopes.
 *
 *        The contact features found at Step 4, called active feature pair is
 *        used at Step 5. We first check the neighbor features for promotion.
 *        For example if the active feature pair is Face-Edge and one of the
 *        incident face of the edge is parallel to the face, then it is 
 *        promoted to Face-Face.
 *        After the update the contact normal is calculated and the active
 *        features are projected onto a plane perpendicular to the contact 
 *        normal, then we find the 2D intersecion (convex polygon) of the 
 *        features.
 *        The vertices of the 2D intersecsion are transformed back to the LCS
 *        of each rigid bodies and they become the contact points pairs.
 *        All the necessary information is stored in a ContactPairInfo object.
 *       
 *        Step 6 generates Jacobian unilateral constraints in the velocity
 *        space, which is used to generate PSD Matrix for MLCP.
 *        To limit the number of constraints per contact (usually 6), we
 *        reduce the number of pairs to at most maxNumConstsPerContact.
 *        The pairs to be removed are selected such that those points make
 *        the most obtuse angles in the intersection.
 *
 *        [Contact Update/Tracking]
 *        To update the existing contact pair, we first try to udpate the 
 *        active feature pair we first try to demote the pair based on the
 *        updated feature normals. Then we try to promote it.
 *        If the change in the configuration is significant enough, we run 
 *        GJK to rediscover the active features using the current normal 
 *        direction as the penetration direction.
 *        If GJK does not find penetration, the contact is removed.
 *        The contact point pairs and the contact normal are found in the same
 *        way as Step 6 above.
 *
 *        [Additional Contact Removal]
 *        After MLCP, if all the unilateral constraints are non-active, i.e., 
 *        all the lambdas are zero, ContactPairInfo is removed.
 *       
 *        [Clean up]
 *        At the end of the simulation step, the geometric configuration will
 *        have been updated and committed, and then the newly discovered
 *        contact pairs are moved to the active pairs list to be tracked
 *        from the next step onwards.
 *        The lambdas of the unilateral constraints are summed up to estimate
 *        the contact pressure along the normal direction, it is used to
 *        approximate the maximum friction force at the next step.
 *
 *        [Estimate Friction]
 *        The approximation formula is given as follows.
 *
 *          |λ_max(t+Δt)| <= μ Σ|λ_n(t)|
 *        
 *        In some engines, the friction is usually handled using a friction
 *        cone which is complementary to the normal pressure force on the
 *        same time step using additional Jacobian constraints.
 *        However, the author thinks it is impossible to handle such
 *        constraints keeping the matrix PSD (and numerically stable),
 *        which is required for the iterative Gauss-Seidel iterative MLCP 
 *        solver. (Makena does not use pivotting solver such as Lemke for
 *        inefficiency, numerical instability, and unpredictability of time)
 *        Hence, Makena delays the handling of friction by one time step
 *        as above. Since the contact surfaces that are affected by
 *        some friction are normally resting orin uniform linear motion, 
 *        the delay by one step will not be a significant issue.
 *
 *        Reasoning:
 *               Pressure:       f_n = Σ[ (1/Δt)Jλ_n ]
 *               Friction force: |f_d| <= μ|f_n|
 *                                f_d = (1/Δt)Jλ_d
 *
 *                               |λ_d| =   |f_d|/|(1/Δt)J|
 *                                     >=  μ|f_n|/|(1/Δt)J|
 *                                     =   μ|Σ[(1/Δt)Jλ_n]|/|(1/Δt)J|
 *                                     =   μΣ[|λ_n|]
 *
 */
namespace Makena {


class ContactManager {

  public:


    using CollisionMapType  = std::map<

        pair<ConvexRigidBody*const, ConvexRigidBody*const>,
        ContactPairInfo*

    >;

    using CollisionMapItType = CollisionMapType::iterator;

    /** @brief constructor
     *
     *  @param KCorr    (in): coefficient for RHS of the constraints.
     *                        0.9 seems to work good for most of the cases.
     *  @param velocityThreshold
     *                  (in): the threshold in the velocity perpendicular
     *                        to the contact normal used to determine if
     *                        the static or dynamic friction coefficient
     *                        is used.
     *  @param maxNumPointsPerContact
     *                  (in): The maximum numper of point pairs used to
     *                        generate unilateral constraints per contact.
     *                        Usually 4 (quadrilateral) or  (hexagon) will
     *                        be sufficient.
     *  @param epsilonZero
     *                  (in): Numerical tolerance to be considered zero.
     *
     *  @param epsilonAngle
     *                  (in): Angular tolerance to be considered zero radian.
     *  @param epsilonZeroGJK
     *                  (in): Numerical tolerance to be considered zero used
     *                        for contact rediscovery. It can be relaxed to
     *                        make the discovery easier.
     *  @param scalingGJK
     *                  (in): Scaling of the convex hull of the rigid body
     *                        used for contact rediscovery. It can be greater
     *                        than 1.0 to make the discovery easier.
     *
     *  @param maxNumIter
     *                  (in): Max number of iterations allowed before
     *                        giving up. This is a safety mechanism.
     *
     *  @param maxNumCycles 
     *                  (in): Max number of occurences where the distance
     *                        to the closest point goes up.
     *                        This is to absorb numerical fluctuation.
     *
     *  @param logStream(in): Log output
     */
    ContactManager(const double      KCorr,
                   const double      velocityThreshold,
                   const long        maxNumPointsPerContact,
                   const double      epsilonZero,
                   const double      epsilonAngle,
                   const double      epsilonZeroGJK,
                   const double      scalingGJK,
                   const long        maxNumIter,
                   const long        maxNumCycles,
                   std::ostream&     logStream
    );


    /** @brief descructor. */
    ~ContactManager();

    /** @brief registers a new ConvexRigidBody which can collide with others.*/
    void registerRigidBody(ConvexRigidBody* p);

    /** @brief unregistere a ConvexRigidBody when it is dying,
     *         or when it is no longer posible to collide with others.
     */
    void unregisterRigidBody(ConvexRigidBody* p);

    /** @brief initialize AABB before the first simulation step.
     *         It makes the first sorted list of objects along the 3 axes.
     */
    void initializeAABB();

    /** @brief initialize for one simulation step. After this call, 
     *         register/unregister of RigidBodies are  not allowed.
     */
    void lockForNextStep(const double deltaT);


    /** @brief terminate one simulation step. After this call, 
     *         register/unregister of RigidBodies are allowed again.
     */
    void unlockAfterStep(ConstraintManager& m);

    /** @brief updates the existing contacts and generates the constraints
     *         based on the current geometric configuration.
     */
    void updateActiveContacts(ConstraintManager& m);

    /** @brief discover the new contacts and generates contraints 
     *         based on the current (temporary) geometric configuration.
     *
     *  @return true  : if at least one new contact is discovered
     *          false : if no new contact is discovered
     */
    bool discoverNewContacts(bool useGeomConfigTemp = true);


    /** @brief assuming the two rigid bodies are now intersecting, 
     *         find the separation axis along which two bodies are
     *         separated with minimum translation.
     *         Such an axis is found by gradually shrinking the 
     *         orienting bounding boxes.
     *
     *  @param body1 (in): body 1
     *
     *  @param body2 (in): body 2
     *
     *  @param separationAxis (out): separating axis
     *
     *  @return true if such a separation axis is found
     */
    bool findSeparationAxis(
        ConvexRigidBody& body1,
        ConvexRigidBody& body2,
        Vec3&            separationAxis
    );


    /** @brief register the constraints of the existing contacts
     *         to the constraint manager.
     */
    void registerActiveConstraints(ConstraintManager& m);

    /** @brief register the constraints of the newly discovered contacts
     *         to the constraint manager.
     */
    void registerNewConstraints(ConstraintManager& m);


    /** @brief unregister and remove all the constraints
     */
    void cleanupJacobianConstraints(ConstraintManager& m);

    CollisionMapType& activePairs() { return mActivePairs; }

    CollisionMapType& newPairs() { return mNewPairs; }

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif
    void removeInactivePairs(ConstraintManager& m);


    void mergeNewPairs();


    void findContactPressure();


    bool discoverNewContactPair(
        bool             useGeomConfigTemp,
        ConvexRigidBody& body1,
        ConvexRigidBody& body2,
        ContactPairInfo& info
    );


    Vec3 findRelativeVelocityOfBody1RelativeToBody2(
        ConvexRigidBody& body1,
        ConvexRigidBody& body2,
        ContactPairInfo& info,
        bool             useGeomConfigTemp
    );


    void getPerpendicularVectors(const Vec3& n, Vec3& u ,Vec3& w);


    bool isRotationalFrictionRequired(ContactPairInfo& info);


    void generateFriction(
        ConvexRigidBody& body1,
        ConvexRigidBody& body2,
        ContactPairInfo& info
    );

    void findFeatureOfRigidBody(
        Manifold&                          ch,
        vector<VertexIt>&                  vits,
        FaceIt&                            fit,
        EdgeIt&                            eit,
        VertexIt&                          vit,
        enum ContactPairInfo::FeatureType& type
    );

    void decomposeQ(
        vector<BDVertex>& Q,
        vector<VertexIt>& vertices1,
        vector<VertexIt>& vertices2
    );

    void generateUnilateralConstraints(
        ConvexRigidBody& body1,
        ConvexRigidBody& body2,
        ContactPairInfo& info
    );

    /** @brief coefficient for RHS of the constraints.*/
    const double                    mKCorr;


    /** @brief the threshold in the velocity perpendicular to the contact 
     *        normal used to determine if the static or dynamic friction
     *        coefficient is used.
     */
    const double                    mVelocityThresholdSquared;


    /** @brief maximum numper of point pairs used to generate unilateral
     *         constraints per contact.
     */
    const long                      mMaxNumPointsPerContact;


    /** @brief numerical tolerance to be considered zero. */
    const double                    mEpsilonZero;


    /** @brief angular tolerance to be considered zero radian.*/
    const double                    mEpsilonAngle;


    /** @brief numerical tolerance to be considered zero used for contact 
     *         rediscovery.
     */
    const double                    mEpsilonZeroGJK;


    /** @brief scaling of the convex hull of the rigid body used for contact 
     *         rediscovery. It can be greater than 1.0 to make the discovery
     *         easier.
     */
    const double                    mScalingGJK;


    /** @brief max number of iterations allowed before giving up. This is a 
     *         safety mechanism.
     */              
    const long                      mMaxNumIter;


    /** @brief max number of occurences where the distance to the closest point
     *         goes up. This is to absorb numerical fluctuation.
     */
    const long                      mMaxNumCycles;


    /** @brief log output stream */
    std::ostream&                   mLogStream;

    /** @brief true if one simulation step has started. During that locked
     *         period, you are not allowed to register or unregister
     *         rigid bodies.
     */
    bool                            mLocked;


    /** @brief 1.0/deltaT is stored here to avoid recalculation.*/
    double                          mDeltaTinv;


    /** @brief KCorr / deltaT is stored here to avoid recalculation.*/
    double                          mKCorrOverDeltaT;


    /** @brief Broad-phase AABB collision detector. */
    BroadPhaseAABBCollisionDetector mAABBDetector;


    /** @brief the set of currently active contact pairs indexed by the pair
     *        of ConvexRigidBodies. The objects are owned here.
     */
    CollisionMapType                mActivePairs;


    /** @brief the set of newly discovered contact pairs indexed by the pair
     *        of ConvexRigidBodies. The objects are owned here.
     */
    CollisionMapType                mNewPairs;

};


}// namespace Makena

#endif /*_MAKENA_CONTACT_MANAGER_HPP_*/

#ifndef _MAKENA_CONTACT_UPDATER_HPP_
#define _MAKENA_CONTACT_UPDATER_HPP_

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
#include "contact_discoverer.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file contact_updater.hpp
 *
 * @brief update a contact pair based on the current geometric configuration.
 *        It tries to demote and then promote the feature pair, and then
 *        update the contact point pairs and the contact normal.
 *
 *        F-F  <---[promotion]-+-[demotion]---> V-V

 *        Strategy to update.
 *
 *        Decompose the direction into 2 types: Normal and Tangential.
 *        First, we check the two bodies in the tangential directions.
 *        We check if the current two features intersects in the XY plane
 *        perpendicular to the current normal.
 *        If not then use the currnent normal as the penetration direction,
 *        run GJK.
 *        Second, we check if the two normals of the current bodies are
 *        permitted for the current feature pairs.
 *        If not, update.
 *        After the feature update, check for penetration around the new
 *        features. If there is a penetration, then find the relative 
 *        penetration direction A->B, and run GJK with dilated convex hull
 *        with the direction.
 *
 *        Process flow
 *        ============
 *
 *        - Active feature pair from the previous time step
 *        - Current Geometry for the new time step.
 * 
 *        1. Update the contact normal based on the current active feature
 *           pair.
 *
 *        2. Decompose the geometric configurations of the two rigid bodys 
 *           into two: One along the contact normal and the one along the
 *           plane perpendicular to the contact normal.
 *
 *        3. Check the 2D overlap of two features.
 *           If there is no overlap of the features on the projected 2D plane,
 *           run GJK with the contact normal as a hint for the penetration 
 *           direction.
 *
 *        4. Check if the normals of the two rigid bodys permit the current
 *           feature pair. Update the feature pair and correct the normal
 *           if necessary.
 *
 *        5. Check if the incident features of the current feature pair are
 *           penetrating to each other. Update the feature pair and correct the
 *           normal if necessary.
 *           If the penetration is too complex, then run GJK with the current
 *           relative velocity around the current active feature as the
 *           penetration direction.
 *
 *        6. Find the 2D intersection of two features into the ordered pairs
 *           of contact points in LCS.
 *
 *        7. Generate Unilateral constraints from the contact points and the
 *           contact normal.
 *
 *        8. If the pressure along the contact normal has been calculated
 *           in the previous step, find the representative contact points pair
 *           (avg) and find the relative velocity. if the velocity is zero,
 *           apply the static friction coefficient and find the max tangential
 *           force, if not, then apply the dynamic friction coefficient and
 *           find the max tangential force.
 *           From the max force, calculate the max lambda for the bilateral
 *           contraints and generate boxed bilateral constraints.
 *
 *
 *   Update Normal
 *
 *   Previous contact generation is irregular => Rediscover.
 *
 *   F-F  => checkForUpdate() => Non-tilting edge case => 
 *                                   E-E, E-V, V-E, V-V
 *                               Tilting case => 
 *                                   F-E, F-V, E-F, E-E, E-V, V-F, V-E, V-V
 *                               Rediscover if there is no overlap along
 *                               face normal 1, or if the tilt is more than 
 *                               90 degrees.
 *
 *   F-E  => checkForUpdate() => Non-tilting edge case => 
 *                                   E-E, E-V, V-E, V-V
 *                               Tilting case => 
 *                                   F-V, E-E, E-V, V-E, V-V
 *                               Rediscover if there is no overlap along
 *                               face normal 1.
 *                                     
 *   F-V  => checkForUpdate() => Non-tilting edge case => 
 *                                   E-V, V-V
 *                               Tilting case => 
 *                                   E-V, V-V
 *                               Rediscover if there is no overlap along
 *                               face normal 1.
 *
 *                                     
 *   E-E  Cross => checkForUpdate() => E-V, V-E, V-V
 *                               Those get the default normal.
 *                               Rediscover if there is no overlap along
 *                               the cross normal.
 *
 *   E-E  Parallel => checkForUpdate() =>
 *                               Rediscover if there is no overlap.
 *
 *   E-V  => checkForUpdate() => V-V 
 *                               Those get the default normal.
 *                               Rediscover if the vertex is not on the edge.
 *   
 *   V-V =>  checkForUpdate() => none
 *                               Rediscover if the vertices don't coincide
 *                               with each other.
 *
 *
 *   Further Check:
 *   Basic Strategy:
 *   - Touching => Update
 *   - Penetration => Rediscovery
 *
 *   F-E  => If an E's incident face is touching => F-F
 *           If en E's incident face is deeply penetrating => Rediscover.
 *
 *   F-V  => If one incident edge is touching => F-E
 *           If two consecutive edges are touching => F-F
 *           If more than two consecutive edges are touching => Rediscover.
 *           If an incident edge is deeply penetrating => Rediscover.
 *
 *   E-E Crossing => If incident face is touching to an edge => E-F, F-E, F-F
 *           If incident face is deeply penetrating edge => Rediscover.
 *
 *   E-E Parallel => If two faces are touching => F-F
 *                   If any two faces are penetrating => Rediscover.
 *
 *   E-V, V-V  => Separation axes test fails => Rediscover.
 *
 *
 *   Adjust Features:
 *   Basic Strategy:
 *   - Features are clean (from GJK or cleaned update)
 *   - No Rediscovery
 *   - Update on touching or deep penetration.
 *
 *   F-E  => If an E's incident face is touching  or penetrating => F-F
 *
 *   F-V  => If two or more consecutive edges are touching => F-F 
 *           If one incident edge is touching => F-E
 *
 *   E-E Crossing => If incident face is touching or penetrating 
 *                   to an edge => E-F, F-E, F-F
 *
 *   E-E Parallel => If two faces are touching or penetrating => F-F
 *
 *   E-V, V-E, V-V => No check.
 *
 *   Process Features
 *
 *
 */

namespace Makena {

using namespace std;

class ContactUpdater : public Loggable {

  public:

    /** @brief constructor.
     *
     *  @param body1        (in): Rigid body 1
     *
     *  @param body2        (in): Rigid body 2
     *
     *  @param info     (in/out): the contact info to be updated
     *
     *  @param epsilonZero  (in): numerical tolerance to consider the
     *                            distance value zero.
     *
     *  @param epsilonAngle (in): numerical tolerance to consider the
     *                            angle zero.
     *
     *  @param logStream    (in): the output stream for logging.
     */
    ContactUpdater(
        ConvexRigidBody&  body1,
        ConvexRigidBody&  body2,
        ContactPairInfo&  info,
        const double&     epsilonZero,
        const double&     epsilonAngle,
        const double&     epsilonZeroGJK,
        const double&     scalingGJK,
        const double&     maxNumPointsPerContact,
        std::ostream&     logStream
    );

    /** @brief destructor. Nothing to be done.
     */
    ~ContactUpdater();


    /** @brief main function that generates a ContactPairInfo.
     *
     *  @return  true if the contact pair is no longer active and
     *           should be removed.
     */
    bool update();


#ifndef UNIT_TESTS
  private:
#else
  public:
#endif

    void updateContactNormal();

    void checkForUpdate(bool& needRediscovery, bool& featurePairUpdated);

    void checkFurther(bool& needRediscovery, bool& featurePairUpdated);

    bool findContactPointPairs();

    void makeInitialBDVertices();

    bool rediscover(bool preferExternalExitDir, const Vec3& exitDir);

    Vec3 findRelativeVelocityOfBody1RelativeToBody2();

    void decomposeQ(
        vector<BDVertex>& Q,
        vector<VertexIt>& vertices1,
        vector<VertexIt>& vertices2
    );


    ConvexRigidBody&         mBody1;
    ConvexRigidBody&         mBody2;
    ContactPairInfo&         mInfo;
    const double             mEpsilonZero;
    const double             mEpsilonAngle;
    const double             mEpsilonZeroGJK;
    const double             mScalingGJK;
    const double             mMaxNumPointsPerContact;
#ifdef UNIT_TESTS
    void showResult();
#endif
};


}// namespace Makena

#endif /*_MAKENA_CONTACT_UPDATER_HPP_*/

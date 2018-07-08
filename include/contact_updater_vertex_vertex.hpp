#ifndef _MAKENA_CONTACT_UPDATER_VERTEX_VERTEX_HPP_
#define _MAKENA_CONTACT_UPDATER_VERTEX_VERTEX_HPP_

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
 * @file contact_updater_vertex_vertex.hpp
 *
 * @brief
 *
 */

namespace Makena {

using namespace std;

class ContactUpdater_VERTEX_VERTEX : public Loggable {

  public:

    ContactUpdater_VERTEX_VERTEX(
        ConvexRigidBody&         body1,
        ConvexRigidBody&         body2,
        ContactPairInfo&         info,
        const double&            epsilonZero,
        const double&            epsilonAngle,
        std::ostream&            logStream
    );

    ~ContactUpdater_VERTEX_VERTEX();

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

    ConvexRigidBody&             mBody1;
    ConvexRigidBody&             mBody2;
    ContactPairInfo&             mInfo;
    const double                 mEpsilonZero;
    const double                 mEpsilonAngle;
};


}// namespace Makena

#endif /*_MAKENA_GJK_CONTACT_UPDATER_VERTEX_VERTEX_HPP_*/

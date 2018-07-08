#ifndef _MAKENA_OBB_OBB_TEST_HPP_
#define _MAKENA_OBB_OBB_TEST_HPP_

#include <memory>
#include <array>
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <exception>
#include <stdexcept>
#include <cmath>

#include "primitives.hpp"
#include "manifold.hpp"


#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file obb_obb_test.hpp
 *
 * @brief Performs separation axis test on a pair of OBBs.
 *
 * @reference
 *
 */
namespace Makena {

using namespace std;


/** @brief performs the 15 separating axis tests on two 
 *         OBBS situated and oriented in GCS.
 *
 *  @param  cob1gcs  (in): center of the bounding box 1 in GCS
 *
 *  @param  axes1gcs (in): The columns of the matrix represents the
 *                         3 orthonormal axes of box 1 oriented in GCS
 *
 *  @param  e1       (in): 3 lengths of the box 1
 *
 *  @param  cob2gcs  (in): center of the bounding box 2 in GCS
 *
 *  @param  axes1gcs (in): The columns of the matrix represents the
 *                         3 orthonormal axes of box 2 oriented in GCS
 *
 *  @param  e2       (in): 3 lengths of the box 2
 *
 *  @param  scaling  (in): Scaling to extents
 *
 *  @param  separationAxis 
 *                   (out) the separating axis if they do not intersect
 *
 *  @return true if the two boxes intersect
 */
bool doIntersect(
    const Vec3&   cob1gcs,
    const Mat3x3& axes1gcs,
    const Vec3&   e1,

    const Vec3&   cob2gcs,
    const Mat3x3& axes2gcs,
    const Vec3&   e2,

    const double  scaling,

    Vec3&         separationAxis
);


void findSeparationAxes(
    const Vec3&   cob1gcs,
    const Mat3x3& axes1gcs,
    const Vec3&   e1,

    const Vec3&   cob2gcs,
    const Mat3x3& axes2gcs,
    const Vec3&   e2,

    const double  scaling,

    vector<Vec3>& separationAxes,
    const double  epsilonZero
);


}// namespace Makena


#endif/*_MAKENA_OBB_OBB_TEST_HPP_*/

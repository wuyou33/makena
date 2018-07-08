#ifndef _MAKENA_BINARY_DILATION_HPP_
#define _MAKENA_BINARY_DILATION_HPP_

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

#include "primitives.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file binary_dilation.hpp
 *
 * @brief  Represents a binary dilation (also known as Minkowski sum) of A- B.
 *         Used by GJK collision detection algorithm.
 */
namespace Makena {

using namespace std;


class BDVertex
{
  public:

    inline BDVertex(VertexIt v1, VertexIt v2):m1(v1),m2(v2){;}
    inline ~BDVertex(){;}

    inline void updateP(
        Manifold&     mani1, 
        const Mat3x3& rotMat1,
        const Vec3&   com1,
        Manifold&     mani2, 
        const Mat3x3& rotMat2,
        const Vec3&   com2
    ) {
        mP = (*m1)->pGCS(rotMat1, com1) - (*m2)->pGCS(rotMat2, com2) ;
    }

    inline void updateP(
        Manifold&     mani1, 
        const double& mScaling1,
        const Mat3x3& rotMat1,
        const Vec3&   com1,
        Manifold&     mani2, 
        const double& mScaling2,
        const Mat3x3& rotMat2,
        const Vec3&   com2
    ) {
        mP = (*m1)->pGCS(mScaling1, rotMat1, com1) - 
             (*m2)->pGCS(mScaling2, rotMat2, com2) ;
    }
    inline const Vec3& p() const {return mP;}
    inline VertexIt& v1() { return m1; }
    inline VertexIt& v2() { return m2; }

#ifndef UNIT_TESTS
  private:
#else
  public:
#endif

    VertexIt m1;   /** @brief vertex on rigid body 1. */
    VertexIt m2;   /** @brief vertex on rigid body 2. */
    Vec3     mP;   /** @brief point in difference configuration space */

friend std::ostream& operator<<(std::ostream& os, const BDVertex& V);

};

inline std::ostream& operator<<(std::ostream& os, const BDVertex& V)
{
    os << "Body1: " << (*(*V.m1)) << "\n";
    os << "Body2: " << (*(*V.m2)) << "\n";
    os << "Point in DS: " << V.mP << "\n";
    return os;
}

} // namespace Makena

#endif /*_MAKENA_BINARY_DILATION_HPP_*/

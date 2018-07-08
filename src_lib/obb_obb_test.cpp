#include "obb_obb_test.hpp"

/**
 * @file obb_obb_test.cpp
 *
 * @brief intersection test between two orienting bounding boxes.
 *
 * @reference
 *
 */
namespace Makena {

using namespace std;

inline bool testOneAxis(
    Vec3&   d,     // normalized direction of test
    Vec3&   t,     // vector from the center of OBB1 to OBB2.
    Vec3&   x1,    // axis X of OBB1
    Vec3&   y1,    // axis Y of OBB1
    Vec3&   z1,    // axis Z of OBB1
    Vec3&   he1,   // half extents of OBB1
    Vec3&   x2,    // axis X of OBB2
    Vec3&   y2,    // axis Y of OBB2
    Vec3&   z2,    // axis Z of OBB2
    Vec3&   he2,   // half extents of OBB2
    double  scaling
) {
    double ra = fabs(d.dot(x1) * he1.x()) + 
                fabs(d.dot(y1) * he1.y()) + 
                fabs(d.dot(z1) * he1.z()) ;
    double rb = fabs(d.dot(x2) * he2.x()) + 
                fabs(d.dot(y2) * he2.y()) + 
                fabs(d.dot(z2) * he2.z()) ;
    return fabs(d.dot(t)) <= (ra + rb)*scaling;
}


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
) {

    Vec3 t   = cob2gcs - cob1gcs;
    Vec3 x1  = axes1gcs.col(1);
    Vec3 y1  = axes1gcs.col(2);
    Vec3 z1  = axes1gcs.col(3);
    Vec3 he1 = e1 * 0.5;

    Vec3 x2  = axes2gcs.col(1);
    Vec3 y2  = axes2gcs.col(2);
    Vec3 z2  = axes2gcs.col(3);
    Vec3 he2 = e2 * 0.5;

    bool res;

    // Axis X of OBB1 (X1)
    x1.normalize();
    res = testOneAxis(x1, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
    if (!res) {
        if ( x1.dot(t) < 0.0 ) {
            x1.scale(-1.0);
        }
        separationAxis = x1;
        return false;
    }

    // Axis Y of OBB1 (Y1)
    y1.normalize();
    res = testOneAxis(y1, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
    if (!res) {
        if ( y1.dot(t) < 0.0 ) {
            y1.scale(-1.0);
        }
        separationAxis = y1;
        return false;
    }


    // Axis Z of OBB1 (Z1)
    z1.normalize();
    res = testOneAxis(z1, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
    if (!res) {
        if ( z1.dot(t) < 0.0 ) {
            z1.scale(-1.0);
        }
        separationAxis = z1;
        return false;
    }


    // Axis X of OBB2 (X2)
    x2.normalize();
    res = testOneAxis(x2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
    if (!res) {
        if ( x2.dot(t) < 0.0 ) {
            x2.scale(-1.0);
        }
        separationAxis = x2;
        return false;
    }


    // Axis Y of OBB2 (Y2)
    y2.normalize();
    res = testOneAxis(y2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
    if (!res) {
        if ( y2.dot(t) < 0.0 ) {
            y2.scale(-1.0);
        }
        separationAxis = y2;
        return false;
    }


    // Axis Z of OBB3 (Z2)
    z2.normalize();
    res = testOneAxis(z2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
    if (!res) {
        if ( z2.dot(t) < 0.0 ) {
            z2.scale(-1.0);
        }
        separationAxis = z2;
        return false;
    }


    // X1 x X2
    auto x1_cr_x2 = x1.cross(x2);
    if (x1_cr_x2.squaredNorm2() > EPSILON_SQUARED) {
        x1_cr_x2.normalize();
        res = testOneAxis(
                     x1_cr_x2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( x1_cr_x2.dot(t) < 0.0 ) {
                x1_cr_x2.scale(-1.0);
            }
            separationAxis = x1_cr_x2;
            return false;
        }
    }


    // X1 x Y2
    auto x1_cr_y2 = x1.cross(y2);
    if (x1_cr_y2.squaredNorm2() > EPSILON_SQUARED) {
        x1_cr_y2.normalize();
        res = testOneAxis(
                     x1_cr_y2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( x1_cr_y2.dot(t) < 0.0 ) {
                x1_cr_y2.scale(-1.0);
            }
            separationAxis = x1_cr_y2;
            return false;
        }
    }


    // X1 x Z2
    auto x1_cr_z2 = x1.cross(z2);
    if (x1_cr_z2.squaredNorm2() > EPSILON_SQUARED) {
        x1_cr_z2.normalize();
        res = testOneAxis(
                     x1_cr_z2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( x1_cr_z2.dot(t) < 0.0 ) {
                x1_cr_z2.scale(-1.0);
            }
            separationAxis = x1_cr_z2;
            return false;
        }
    }


    // Y1 x X2
    auto y1_cr_x2 = y1.cross(x2);
    if (y1_cr_x2.squaredNorm2() > EPSILON_SQUARED) {
        y1_cr_x2.normalize();
        res = testOneAxis(
                     y1_cr_x2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( y1_cr_x2.dot(t) < 0.0 ) {
                y1_cr_x2.scale(-1.0);
            }
            separationAxis = y1_cr_x2;
            return false;
        }
    }


    // Y1 x Y2
    auto y1_cr_y2 = y1.cross(y2);
    if (y1_cr_y2.squaredNorm2() > EPSILON_SQUARED) {
        y1_cr_y2.normalize();
        res = testOneAxis(
                     y1_cr_y2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( y1_cr_y2.dot(t) < 0.0 ) {
                y1_cr_y2.scale(-1.0);
            }
            separationAxis = y1_cr_y2;
            return false;
        }
    }


    // Y1 x Z2
    auto y1_cr_z2 = y1.cross(z2);
    if (y1_cr_z2.squaredNorm2() > EPSILON_SQUARED) {
        y1_cr_z2.normalize();
        res = testOneAxis(
                     y1_cr_z2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( y1_cr_z2.dot(t) < 0.0 ) {
                y1_cr_z2.scale(-1.0);
            }
            separationAxis = y1_cr_z2;
            return false;
        }
    }


    // Z1 x X2
    auto z1_cr_x2 = z1.cross(x2);
    if (z1_cr_x2.squaredNorm2() > EPSILON_SQUARED) {
        z1_cr_x2.normalize();
        res = testOneAxis(
                     z1_cr_x2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            separationAxis = z1_cr_x2;
            return false;
        }
    }


    // Z1 x Y2
    auto z1_cr_y2 = z1.cross(y2);
    if (z1_cr_y2.squaredNorm2() > EPSILON_SQUARED) {
        z1_cr_y2.normalize();
        res = testOneAxis(
                     z1_cr_y2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( z1_cr_y2.dot(t) < 0.0 ) {
                z1_cr_y2.scale(-1.0);
            }
            separationAxis = z1_cr_y2;
            return false;
        }
    }


    // Z1 x Z2
    auto z1_cr_z2 = z1.cross(z2);
    if (z1_cr_z2.squaredNorm2() > EPSILON_SQUARED) {
        z1_cr_z2.normalize();
        res = testOneAxis(
                     z1_cr_z2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( z1_cr_z2.dot(t) < 0.0 ) {
                z1_cr_z2.scale(-1.0);
            }
            separationAxis = z1_cr_z2;
            return false;
        }
    }

    return true;
}


void findSeparationAxesCheckAndAdd(
    const Vec3&   v1,
    vector<Vec3>& separationAxes,
    const double  epsilonZero
) {

    bool found = false;

    for ( auto& v2 : separationAxes ) {

        auto diff = v2 - v1;
        if ( diff.squaredNorm2() <= epsilonZero ) {
            found = true;
            break;
        }
    }
    if (!found) {
        separationAxes.push_back(v1);
    }
}


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
) {

    Vec3 t   = cob2gcs - cob1gcs;
    Vec3 x1  = axes1gcs.col(1);
    Vec3 y1  = axes1gcs.col(2);
    Vec3 z1  = axes1gcs.col(3);
    Vec3 he1 = e1 * 0.5;

    Vec3 x2  = axes2gcs.col(1);
    Vec3 y2  = axes2gcs.col(2);
    Vec3 z2  = axes2gcs.col(3);
    Vec3 he2 = e2 * 0.5;

    bool res;

    // Axis X of OBB1 (X1)
    x1.normalize();
    res = testOneAxis(x1, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
    if (!res) {
        if ( x1.dot(t) < 0.0 ) {
            x1.scale(-1.0);
        }
        findSeparationAxesCheckAndAdd( x1, separationAxes, epsilonZero );
    }

    // Axis Y of OBB1 (Y1)
    y1.normalize();
    res = testOneAxis(y1, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
    if (!res) {
        if (y1.dot(t)<0.0) {
            y1.scale(-1.0);
        }
        findSeparationAxesCheckAndAdd( y1, separationAxes, epsilonZero );
    }


    // Axis Z of OBB1 (Z1)
    z1.normalize();
    res = testOneAxis(z1, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
    if (!res) {
        if ( z1.dot(t) < 0.0 ) {
            z1.scale(-1.0);
        }
        findSeparationAxesCheckAndAdd( z1, separationAxes, epsilonZero );
    }


    // Axis X of OBB2 (X2)
    x2.normalize();
    res = testOneAxis(x2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
    if (!res) {
        if ( x2.dot(t) < 0.0 ) {
             x2.scale(-1.0);
        }
        findSeparationAxesCheckAndAdd( x2, separationAxes, epsilonZero );
    }


    // Axis Y of OBB2 (Y2)
    y2.normalize();
    res = testOneAxis(y2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
    if (!res) {
        if ( y2.dot(t) < 0.0 ) {
            y2.scale(-1.0);
        }
        findSeparationAxesCheckAndAdd( y2, separationAxes, epsilonZero );
    }


    // Axis Z of OBB3 (Z2)
    z2.normalize();
    res = testOneAxis(z2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
    if (!res) {
        if ( z2.dot(t) < 0.0 ) {
            z2.scale(-1.0);
        }
        findSeparationAxesCheckAndAdd( z2, separationAxes, epsilonZero );
    }


    // X1 x X2
    auto x1_cr_x2 = x1.cross(x2);
    if (x1_cr_x2.squaredNorm2() > EPSILON_SQUARED) {
        x1_cr_x2.normalize();
        res = testOneAxis(
                     x1_cr_x2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( x1_cr_x2.dot(t) < 0.0 ) {
                x1_cr_x2.scale(-1.0);
            }
            findSeparationAxesCheckAndAdd(x1_cr_x2,separationAxes,epsilonZero);
        }
    }


    // X1 x Y2
    auto x1_cr_y2 = x1.cross(y2);
    if (x1_cr_y2.squaredNorm2() > EPSILON_SQUARED) {
        x1_cr_y2.normalize();
        res = testOneAxis(
                     x1_cr_y2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( x1_cr_y2.dot(t) < 0.0 ) {
                x1_cr_y2.scale(-1.0);
            }
            findSeparationAxesCheckAndAdd(x1_cr_y2,separationAxes,epsilonZero);
        }
    }


    // X1 x Z2
    auto x1_cr_z2 = x1.cross(z2);
    if (x1_cr_z2.squaredNorm2() > EPSILON_SQUARED) {
        x1_cr_z2.normalize();
        res = testOneAxis(
                     x1_cr_z2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( x1_cr_z2.dot(t) < 0.0 ) {
                x1_cr_z2.scale(-1.0);
            }
            findSeparationAxesCheckAndAdd(x1_cr_z2,separationAxes,epsilonZero);
        }
    }


    // Y1 x X2
    auto y1_cr_x2 = y1.cross(x2);
    if (y1_cr_x2.squaredNorm2() > EPSILON_SQUARED) {
        y1_cr_x2.normalize();
        res = testOneAxis(
                     y1_cr_x2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( y1_cr_x2.dot(t) < 0.0 ) {
                y1_cr_x2.scale(-1.0);
            }
            findSeparationAxesCheckAndAdd(y1_cr_x2,separationAxes,epsilonZero);
        }
    }


    // Y1 x Y2
    auto y1_cr_y2 = y1.cross(y2);
    if (y1_cr_y2.squaredNorm2() > EPSILON_SQUARED) {
        y1_cr_y2.normalize();
        res = testOneAxis(
                     y1_cr_y2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( y1_cr_y2.dot(t) < 0.0 ) {
                y1_cr_y2.scale(-1.0);
            }
            findSeparationAxesCheckAndAdd(y1_cr_y2,separationAxes,epsilonZero);
        }
    }


    // Y1 x Z2
    auto y1_cr_z2 = y1.cross(z2);
    if (y1_cr_z2.squaredNorm2() > EPSILON_SQUARED) {
        y1_cr_z2.normalize();
        res = testOneAxis(
                     y1_cr_z2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( y1_cr_z2.dot(t) < 0.0 ) {
                y1_cr_z2.scale(-1.0);
            }
            findSeparationAxesCheckAndAdd(y1_cr_z2,separationAxes,epsilonZero);
        }
    }


    // Z1 x X2
    auto z1_cr_x2 = z1.cross(x2);
    if (z1_cr_x2.squaredNorm2() > EPSILON_SQUARED) {
        z1_cr_x2.normalize();
        res = testOneAxis(
                     z1_cr_x2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( z1_cr_x2.dot(t) < 0.0 ) {
                z1_cr_x2.scale(-1.0);
            }
            findSeparationAxesCheckAndAdd(z1_cr_x2,separationAxes,epsilonZero);
        }
    }


    // Z1 x Y2
    auto z1_cr_y2 = z1.cross(y2);
    if (z1_cr_y2.squaredNorm2() > EPSILON_SQUARED) {
        z1_cr_y2.normalize();
        res = testOneAxis(
                     z1_cr_y2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( z1_cr_y2.dot(t) < 0.0 ) {
                z1_cr_y2.scale(-1.0);
            }
            findSeparationAxesCheckAndAdd(z1_cr_y2,separationAxes,epsilonZero);
        }
    }


    // Z1 x Z2
    auto z1_cr_z2 = z1.cross(z2);
    if (z1_cr_z2.squaredNorm2() > EPSILON_SQUARED) {
        z1_cr_z2.normalize();
        res = testOneAxis(
                     z1_cr_z2, t, x1, y1, z1, he1, x2, y2, z2, he2, scaling);
        if (!res) {
            if ( z1_cr_z2.dot(t) < 0.0 ) {
                z1_cr_z2.scale(-1.0);
            }
            findSeparationAxesCheckAndAdd(z1_cr_z2,separationAxes,epsilonZero);
        }
    }
}

}// namespace Makena


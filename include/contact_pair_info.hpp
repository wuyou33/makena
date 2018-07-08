#ifndef _MAKENA_CONTACT_PAIR_INFO_HPP_
#define _MAKENA_CONTACT_PAIR_INFO_HPP_

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <exception>

#include "primitives.hpp"
#include "manifold.hpp"
#include "binary_dilation.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file contact_pair_info.hpp
 *
 * @brief Represents a contact between two convex rigid bodies. It is used to
 *        track contacts over multiple iterations.  Makena tries tokeep track
 *        of acontactby a feature pair. A feature is one of vertex, edge, or 
 *        face. It alsoholds the current set of contact point pairs inLCS 
 *        foreach body, and thecontactdirection, which will be used to generate
 *        a setof unilateral constraints for the current time step, and to 
 *        generate a set of boxed bilateral constraints forthe next step if the
 *        current feature pair is still in effect.
 *
 *        NOTE: It owns JacobianConstraints kept in mUniConstraints &
 *              mBiConstraints. They are constructed by new, and destructed
 *              by ContactPairInfo::deleteConstraints() in 
 *              ContactManager. Their lifetime is within a simulation step.
 */
namespace Makena {


class ContactPairInfo {

  public:    

    enum FeatureType {
        FT_NONE,
        FT_VERTEX,
        FT_EDGE,
        FT_FACE
    };

    enum NormalType {
        NT_NONE,
        NT_FACE_NORMAL_1,
        NT_FACE_NORMAL_2,
        NT_EDGE_NORMAL_1,
        NT_EDGE_NORMAL_2,
        NT_EDGE_CROSS_EDGE,
        NT_EDGE_EDGE_AVG,
        NT_VERTEX_VERTEX_AVG,
        NT_RELATIVE_VELOCITY
    };

    ContactPairInfo():
        mType1(FT_NONE), 
        mType2(FT_NONE), 
        mContactNormalType(NT_NONE),
        mInactiveCount(0){;}

    ~ContactPairInfo() {
        for (auto* p : mUniConstraints) {
            delete p;
        }
        for (auto* p : mBiConstraints) {
            delete p;
        }
    }

    void deleteConstraints() {

        for (auto* p : mUniConstraints) {
            delete p;
        }
        for (auto* p : mBiConstraints) {
            delete p;
        }
        mUniConstraints.clear();
        mBiConstraints.clear();

    }

    enum FeatureType                     mType1;
    FaceIt                               mFit1;
    EdgeIt                               mEit1;
    VertexIt                             mVit1;

    enum FeatureType                     mType2;
    FaceIt                               mFit2;
    EdgeIt                               mEit2;
    VertexIt                             mVit2;

    Vec3                                 mContactNormal1To2;
    enum NormalType                      mContactNormalType;
    bool                                 mContactGenerationIrregular;
    vector<Vec3>                         mIntsect1;
    vector<Vec3>                         mIntsect2;
    vector<Vec3>                         mIntsect1reduced;
    vector<Vec3>                         mIntsect2reduced;

    std::vector<Makena::BDVertex>        mLastVertices;

    /* The force along the contact normal. Used to generate the 
     * friction constraints based on Hooke's law.
     * This is equal (but oppsite direction) to the repulsive force
     * induced by ConstraintEngine to avoid penetration.
     */
    double                               mSumLambdas;

    /* Conact constraints are owned here */        
    std::vector<JacobianConstraint*>     mUniConstraints;
    std::vector<JacobianConstraint*>     mBiConstraints;

    /* Number of consecutive steps in the immdeiate past
     * at which this pair is inactive.
     */
    long                                 mInactiveCount;

};



}// namespace Makena

#endif /*_MAKENA_CONTACT_PAIR_INFO_HPP_*/

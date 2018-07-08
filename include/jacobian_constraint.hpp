#ifndef _MAKENA_JACOBIAN_CONSTRAINT_HPP_
#define _MAKENA_JACOBIAN_CONSTRAINT_HPP_

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <exception>

#include "primitives.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file jacobian_constraint.hpp
 *
 * @brief It represents aJacobian constraint in the velocity space for
 *        constraint-based simulation. It is one of bilateral free, bilateral
 *        boxed, and unilateral constraint.
 */
namespace Makena {

class ConvexRigidBody;

class JacobianConstraint {

  public:

    enum type {
        UNILATERAL,
        BILATERAL_FREE,
        BILATERAL_BOX
    };

    enum engineState {
        UNREGISTERED,
        NOT_CHECKED_IN,
        CHECKED_IN,
        END
    };


    JacobianConstraint(
        enum type        t,
        ConvexRigidBody* body1,
        ConvexRigidBody* body2
    ):
        mType(t),
        mEngineState(UNREGISTERED),
        mBody1(body1),
        mBody2(body2),
        mMatIndex(-1),
        mBoxLimit(0.0),
        mLambda(0.0) {;}

    ~JacobianConstraint(){;}

    enum type type() { return mType; }

    void setLHSlin1 (const Vec3&  v){ 
        mLHSlin1    = v; 
    }

    void setLHSlin1 (const double& x, const double& y, const double& z) {
        mLHSlin1.setX(x);
        mLHSlin1.setY(y);
        mLHSlin1.setZ(z);
    }

    void setLHSlin2 (const Vec3&  v){
        mLHSlin2    = v; 
    }

    void setLHSlin2 (const double& x, const double& y, const double& z) {
        mLHSlin2.setX(x);
        mLHSlin2.setY(y);
        mLHSlin2.setZ(z);
    }

    void setLHSang1 (const Vec3&  v){ 
        mLHSang1    = v;
    }

    void setLHSang1 (const double& x, const double& y, const double& z) {
        mLHSang1.setX(x);
        mLHSang1.setY(y);
        mLHSang1.setZ(z);
    }

    void setLHSang2 (const Vec3&  v){
        mLHSang2  = v;
    }

    void setLHSang2 (const double& x, const double& y, const double& z) {
        mLHSang2.setX(x);
        mLHSang2.setY(y);
        mLHSang2.setZ(z);
    }

    void setRHS     (const double v){
        mRHS       = v;
    }

    void setBoxLimit(const double v){ mBoxLimit = v; }

    const Vec3& LHSlin1()  const { return mLHSlin1;  }
    const Vec3& LHSlin2()  const { return mLHSlin2;  }
    const Vec3& LHSang1()  const { return mLHSang1;  }
    const Vec3& LHSang2()  const { return mLHSang2;  }
    double      RHS ()     const { return mRHS;      }
    double      boxLimit() const { return mBoxLimit; }

    void setMatIndex(long i)     { mMatIndex     = i; }
    long matIndex()        const { return mMatIndex;  }

    void        setLambda(const double& v) {mLambda = v;}
    double      lambda() const {return mLambda;}

    ConvexRigidBody* body1() {return mBody1;}
    ConvexRigidBody* body2() {return mBody2;}

    void setEngineState(enum engineState s) { mEngineState =s; }
    enum engineState engineState() const { return mEngineState; }

     void setBackItConstraintManager(list<JacobianConstraint*>::iterator it)
    {
        mBackItConstraintManager = it;
    }

    list<JacobianConstraint*>::iterator backItConstraintManager()
    {
        return mBackItConstraintManager;
    }


    void setBackItBody1(list<JacobianConstraint*>::iterator it) {
        mBackItBody1 = it;
    }

    void setBackItBody2(list<JacobianConstraint*>::iterator it) {
        mBackItBody2 = it;
    }

    list<JacobianConstraint*>::iterator backItBody1() { return mBackItBody1; }
    list<JacobianConstraint*>::iterator backItBody2() { return mBackItBody2; }

  private:

    enum type              mType;  
    enum engineState       mEngineState;
    ConvexRigidBody*       mBody1;
    ConvexRigidBody*       mBody2;
    long                   mMatIndex;
    Vec3                   mLHSlin1;
    Vec3                   mLHSlin2;
    Vec3                   mLHSang1;
    Vec3                   mLHSang2;
    double                 mRHS;
    double                 mBoxLimit;
    double                 mLambda;

    list<JacobianConstraint*>::iterator mBackItConstraintManager;
    list<JacobianConstraint*>::iterator mBackItBody1;
    list<JacobianConstraint*>::iterator mBackItBody2;
};


}// namespace Makena

#endif /*_MAKENA_JACOBIAN_CONSTRAINT_HPP_*/

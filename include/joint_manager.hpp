#ifndef _MAKENA_JOINT_MANAGER_HPP_
#define _MAKENA_JOINT_MANAGER_HPP_

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <exception>

#include "primitives.hpp"
#include "convex_rigid_body.hpp"
#include "jacobian_constraint.hpp"
#include "constraint_manager.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file joint_manager.hpp
 *
 * @brief
 *  Manages the joints, and generates a set of Jacobian constraints at
 *  each simulation step.
 *
 *   - Piston joint (DoF 2)
 *   - Slider joint (DoF 1)
 *   - Fixed  joint (DoF 0)
 *   - Hinge 2/Universal joint (DoF 2)
 *   - Hinge 1 joint (DoF 1)
 *   - Ball joint (DoF 3)
 *  The linkages are implemented by bilateral free constraints.
 *  The joint limits are implemented by unilateral constraints.
 *  The joint frictions/motors areimplemented by bilateral boxed constraints.
 *  The joint velocity limiters areimplemented by unilateral constraints.

 */
namespace Makena {

class JointManager;

class Joint {

  public:

    Joint(ConvexRigidBody* body1, ConvexRigidBody* body2):
        mBody1(body1),
        mBody2(body2) {;}

    virtual ~Joint() {;
        for (auto* p : mConstraints) {
            delete p;
        }
    }

    virtual void update(const double deltaTinv){;}

    vector<JacobianConstraint*>& constraints() {
        return mConstraints;
    }

    inline bool isBody1Fixed() { 
        if (mBody1==nullptr) {
            return true;
        }
        else {
            return mBody1->isFixed();
        }
    }

    inline bool isBody2Fixed() { 
        if (mBody2==nullptr) {
            return true;
        }
        else {
            return mBody2->isFixed();
        }
    }

    inline JacobianConstraint* newUnilateral() {
        auto* p = new JacobianConstraint(
                            JacobianConstraint::UNILATERAL, mBody1, mBody2);
        clearConstraint(p);
        mConstraints.push_back(p);
        return p;
    }

    inline JacobianConstraint* newBilateralFree() {
        auto* p = new JacobianConstraint(
                          JacobianConstraint::BILATERAL_FREE, mBody1, mBody2);
        clearConstraint(p);
        mConstraints.push_back(p);
        return p;
    }

    inline JacobianConstraint* newBilateralBox() {
        auto* p = new JacobianConstraint(
                           JacobianConstraint::BILATERAL_BOX, mBody1, mBody2);
        clearConstraint(p);
        mConstraints.push_back(p);
        return p;
    }

    inline void clearConstraint(JacobianConstraint* c) {
        c->setLHSlin1(0.0, 0.0, 0.0);
        c->setLHSlin2(0.0, 0.0, 0.0);
        c->setLHSang1(0.0, 0.0, 0.0);
        c->setLHSang2(0.0, 0.0, 0.0);
        c->setRHS(0.0);
    }

    inline void updateConstraint(JacobianConstraint* c,
                                 const Vec3&         lin1,
                                 const Vec3&         ang1,
                                 const Vec3&         lin2,
                                 const Vec3&         ang2,
                                 const double&       rhs   )
    {
        if (!isBody1Fixed()) {
            c->setLHSlin1(lin1);
            c->setLHSang1(ang1);
        }
        if (!isBody2Fixed()) {
            c->setLHSlin2(lin2);
            c->setLHSang2(ang2);
        }
        c->setRHS(rhs);
    }

    inline void updateConstraint(JacobianConstraint* c,
                                 const Vec3&         lin1,
                                 const Vec3&         ang1,
                                 const Vec3&         lin2,
                                 const Vec3&         ang2 )
    {
        if (!isBody1Fixed()) {
            c->setLHSlin1(lin1);
            c->setLHSang1(ang1);
        }
        if (!isBody2Fixed()) {
            c->setLHSlin2(lin2);
            c->setLHSang2(ang2);
        }
    }

    inline void updateConstraintAng(JacobianConstraint* c,
                                    const Vec3&         ang1,
                                    const Vec3&         ang2,
                                    const double&       rhs   )
    {
        if (!isBody1Fixed()) {
            c->setLHSang1(ang1);
        }
        if (!isBody2Fixed()) {
            c->setLHSang2(ang2);
        }
        c->setRHS(rhs);
    }

    inline void updateConstraintAng(JacobianConstraint* c,
                                    const Vec3&         ang1,
                                    const Vec3&         ang2 )
    {
        if (!isBody1Fixed()) {
            c->setLHSang1(ang1);
        }
        if (!isBody2Fixed()) {
            c->setLHSang2(ang2);
        }
    }

    inline double asinSafe(const double& v) {
        return asin(std::max(-1.0, std::min(v, 1.0)));
    }

    inline double angleRad(const Vec3& base, const Vec3& v, const Vec3& axis) {
        Vec3   b_cr_v  = base.cross(v);
        double b_dot_v = base.dot(v);
        double theta   = asinSafe(b_cr_v.norm2());
        if (b_cr_v.dot(axis) >= 0.0) {
            if (b_dot_v >= 0.0) {
                return theta;
            }
            else {
                return M_PI - theta;
            }
        }
        else {
            if (b_dot_v >= 0.0) {
                return -1.0 * theta;
            }
            else {
                return theta - M_PI;
            }
        }
    }

  protected:
    ConvexRigidBody*            mBody1;
    ConvexRigidBody*            mBody2;
    vector<JacobianConstraint*> mConstraints;
    list<Joint*>::iterator      mBackIt;

friend class JointManager;

};


class JointManager {
  public:
    void registerJoint(Joint* j, ConstraintManager& m) {
        j->mBackIt = mJoints.insert(mJoints.end(),j);
        for (auto* p : j->constraints()) {
            m.registerConstraint(p);
        }
    }

    void unregisterJoint(Joint* j, ConstraintManager& m) {
        for (auto* p : j->constraints()) {
            m.unregisterConstraint(p);
        }
        mJoints.erase(j->mBackIt);
    }

    void update(const double deltaTinv) {;
        for (auto* j : mJoints) {
            j->update(deltaTinv);
        }
    }

  private:
    list<Joint*> mJoints;    
    
};


class PistonJoint : public Joint {


  public:

    /** @brief constructor
     *
     *  @param   body1        (in): ConvexRigidBody 1.
     *
     *  @param   body2        (in): ConvexRigidBody 2.
     *
     *  @param   a1          (in): Anchor point on Body1 in LCS. Or an anchor
     *                             point in GCS if ConvexRigidBody1 is null.
     *
     *  @param   a2          (in): Anchor point on Body2 in LCS. Or an anchor
     *                             point in GCS if ConvexRigidBody2 is null.
     *
     *  @param   s1          (in): unit vector for 
     *                             Sliding direction in Body1's LCS, or in GCS
     *                             if ConvexRigidBody1 is null.
     *  @param   u1          (in): unit direction vector perpendicular to s1
     *                             in Body1's LCS, or in GCS if 
     *                             ConvexRigidBody1 is null.
     *
     *  @param   s2          (in): unit vector for 
     *                             Sliding direction in Body2's LCS, or in GCS
     *                             if ConvexRigidBody2 is null.
     *  @param   u2          (in): unit direction vector perpendicular to s2
     *                             in Body2's LCS, or in GCS if 
     *                             ConvexRigidBody2 is null.
     *
     *  @param   useLimitsLengthwise
     *                       (in): Set to true if you want to use the limit
     *                             along the linear movement along s1.
     *
     *  @param   dMin        (in): The lowest point along s1 (=s2 in GCS).
     *
     *  @param   dMax        (in): The highest point along s1 (=s2 in GCS).
     *
     *  @param   useFrictionOrMotorLengthwise
     *                       (in): Set to true if you want to activate
     *                             friction or motor along s1.
     *
     *  @param   wDesiredLengthwise 
     *                       (in): the desired linear speed along s1 (=s2).
     *                             Must be zero if you want to implement
     *                             friction, rather than a motor.
     *
     *  @param   lambdaMaxLengthwise 
     *                       (in): The range 
     *                             [-lamdaMaxLengthwise, lambdaMaxLengthwise]
     *                             is used for the boxed constraints to 
     *                             implement friction or motor. 
     *                             This must be calculated by the application
     *                             based on the viscosity of the joint or the
     *                             maximum force of the linear motor.
     *
     *  @param   useVelocityLimitLengthwise
     *                       (in): Set to true if you want to use velocity 
     *                             limit. The relative angular velocity around
     *                             s1 is restricted in
     *                             [-wLimitLengthwise, wLimitLengthwise].
     *
     *  @param   wLimitLengthwise
     *                       (in): Maximum velocity around around s1.
     *                             Used if useVelocityLimit is set to true.
     *
     *  @param   useLimitsRotational
     *                       (in): Set to true if you want to use the angular 
     *                             limit.
     *
     *  @param   phiMin      (in): The lowest point along s1 (=s2 in GCS).
     *
     *  @param   phiMax      (in): The highest point along s1 (=s2 in GCS).
     *
     *  @param   useFrictionOrMotorRotational
     *                       (in): Set to true if you want to activate
     *                             friction or motor.
     *
     *  @param   wDesiredRotational
     *                       (in): the desired angular speed around s1 (=s2).
     *                             Must be zero if you want to implement
     *                             friction, rather than a motor.
     *
     *  @param   lambdaMaxRotational
     *                       (in): The range
     *                             [-lamdaMaxRotational, lambdaMaxRotational]
     *                             is used for the boxed constraints to 
     *                             implement friction or motor. 
     *                             This must be calculated by the application
     *                             based on the viscosity of the joint or the
     *                             maximum torque of the motor.
     *
     *  @param   useVelocityLimitRotational
     *                       (in): Set to true if you want to use velocity 
     *                             limit. The relative angular velocity around
     *                             s1 is restricted in 
     *                             [-wLimitRotational, wLimitRotational].
     *
     *  @param   wLimitRotational
     *                       (in): Maximum velocity around around s1.
     *                             Used if useVelocityLimit is set to true.
     *
     *  @param   kCoeff      (in): The coefficient to be multiplied to the RHS
     *                             of the linkage constraints to avoid over
     *                             compensation by numerical error.
     *                             Usually 0.8 or 0.9.
     */
    PistonJoint(

        ConvexRigidBody*         body1,
        ConvexRigidBody*         body2,
        const Vec3&   a1,
        const Vec3&   a2,
        const Vec3&   s1,
        const Vec3&   u1,
        const Vec3&   s2,
        const Vec3&   u2,
        const bool&   useLimitsLengthwise,
        const double& dMin,
        const double& dMax,
        const bool&   useFrictionOrMotorLengthwise,
        const double& wDesiredLengthwise,
        const double& lambdaMaxLengthwise,
        const bool&   useVelocityLimitLengthwise,
        const double& wLimitLengthwise,
        const bool&   useLimitsRotational,
        const double& phiMin,
        const double& phiMax,
        const bool&   useFrictionOrMotorRotational,
        const double& wDesiredRotational,
        const double& lambdaMaxRotational,
        const bool&   useVelocityLimitRotational,
        const double& wLimitRotational,
        const double& kCoeff

    ):

        Joint(body1, body2),
        ma1(a1),
        ma2(a2),
        ms1(s1),
        mu1(u1),
        ms2(s2),
        mu2(u2),
        mUseLimitsLengthwise(useLimitsLengthwise),
        mdMin(dMin),
        mdMax(dMax),
        mUseFrictionOrMotorLengthwise(useFrictionOrMotorLengthwise),
        mwDesiredLengthwise(wDesiredLengthwise),
        mLambdaMaxLengthwise(lambdaMaxLengthwise),
        mUseVelocityLimitLengthwise(useVelocityLimitLengthwise),
        mwLimitLengthwise(wLimitLengthwise),
        mUseLimitsRotational(useLimitsRotational),
        mPhiMin(phiMin),
        mPhiMax(phiMax),
        mUseFrictionOrMotorRotational(useFrictionOrMotorRotational),
        mwDesiredRotational(wDesiredRotational),
        mLambdaMaxRotational(lambdaMaxRotational),
        mUseVelocityLimitRotational(useVelocityLimitRotational),
        mwLimitRotational(wLimitRotational),
        mkCoeff(kCoeff)

    {

        ms1.normalize();
        ms2.normalize();
        mu1.normalize();
        mu2.normalize();

        mcLink1 = newBilateralFree();
        mcLink2 = newBilateralFree();
        mcLink3 = newBilateralFree();
        mcLink4 = newBilateralFree();

        if (mUseLimitsLengthwise) {

            mcLimit1 = newUnilateral();
            mcLimit2 = newUnilateral();

        }

        if (mUseLimitsRotational) {

            mcLimit3 = newUnilateral();
            mcLimit4 = newUnilateral();

        }

        if (mUseFrictionOrMotorLengthwise) {

            mcFric1 = newBilateralBox();
            mcFric1->setRHS(mwDesiredLengthwise);
            mcFric1->setBoxLimit(mLambdaMaxLengthwise);

        }

        if (mUseFrictionOrMotorRotational) {

            mcFric2 = newBilateralBox();
            mcFric2->setRHS(mwDesiredRotational);
            mcFric2->setBoxLimit(mLambdaMaxRotational);

        }

        if (mUseVelocityLimitLengthwise) {

            mcVelLimit1 = newUnilateral();
            mcVelLimit1->setRHS(-1.0 * mwLimitLengthwise);

            mcVelLimit2 = newUnilateral();
            mcVelLimit2->setRHS(-1.0 * mwLimitLengthwise);

        }

        if (mUseVelocityLimitRotational) {

            mcVelLimit3 = newUnilateral();
            mcVelLimit3->setRHS(-1.0 * mwLimitRotational);

            mcVelLimit4 = newUnilateral();
            mcVelLimit4->setRHS(-1.0 * mwLimitRotational);

        }

        if (isBody1Fixed()) {

            ma1_gcs = ma1;
            mr1_gcs = ma1;
            ms1_gcs = ms1;
            mu1_gcs = mu1;
            mRa1.zero();

        }

        if (isBody2Fixed()) {

            ma2_gcs = ma2;
            mr2_gcs = ma2;
            ms2_gcs = ms2;
            mu2_gcs = mu2;
            mRa2.zero();
        }
    }


    void update(const double deltaTinv) {

        double coeffNeg = -1.0 * mkCoeff * deltaTinv;

        if (!isBody1Fixed()) {

            ma1_gcs = mBody1->Qmat() * ma1;
            mr1_gcs = mBody1->CoM() + ma1_gcs;
            ms1_gcs = mBody1->Qmat() * ms1;
            mu1_gcs = mBody1->Qmat() * mu1;
            mRa1    = ma1_gcs.crossMat();
        }

        if (!isBody2Fixed()) {

            ma2_gcs = mBody2->Qmat() * ma2;
            mr2_gcs = mBody2->CoM() + ma2_gcs;
            ms2_gcs = mBody2->Qmat() * ms2;
            mu2_gcs = mBody2->Qmat() * mu2;
            mRa2    = ma2_gcs.crossMat();
        }

        Vec3 r21   = mr1_gcs - mr2_gcs;
        Vec3 Ra1u2 = mRa1 * mu2_gcs;
        Vec3 Ra2u2 = mRa2 * mu2_gcs;

        updateConstraint(mcLink1, mu2_gcs *  1.0,  Ra1u2   *  1.0,
                                  mu2_gcs * -1.0,  Ra2u2   * -1.0,
                                  mu2_gcs.dot(r21) * coeffNeg  );
        auto w2_gcs = ms2_gcs.cross(mu2_gcs);
        Vec3 Ra1w2 = mRa1 * w2_gcs;
        Vec3 Ra2w2 = mRa2 * w2_gcs;

        updateConstraint(mcLink2, w2_gcs *  1.0,  Ra1w2   *  1.0,
                                  w2_gcs * -1.0,  Ra2w2   * -1.0,
                                  w2_gcs.dot(r21) * coeffNeg  );

        Vec3 s2_cr_s1 = ms2_gcs.cross(ms1_gcs);

        double phi21u = asinSafe(mu2_gcs.dot(s2_cr_s1));
        double phi21w = asinSafe( w2_gcs.dot(s2_cr_s1));

        updateConstraintAng(mcLink3,mu2_gcs*1.0,mu2_gcs*-1.0,phi21u*coeffNeg);
        updateConstraintAng(mcLink4, w2_gcs*1.0, w2_gcs*-1.0,phi21w*coeffNeg);

        Vec3 Ra1s2 = mRa1 * ms2_gcs;
        Vec3 Ra2s2 = mRa2 * ms2_gcs;

        if (mUseLimitsLengthwise) {

            updateConstraint(mcLimit1, ms2_gcs,      Ra1s2,
                                       ms2_gcs*-1.0, Ra2s2*-1.0,
                                       (ms2_gcs.dot(r21) - mdMin) * coeffNeg); 


            updateConstraint(mcLimit2, ms2_gcs*-1.0, Ra1s2*-1.0, 
                                       ms2_gcs,      Ra2s2,
                                       (mdMax - ms2_gcs.dot(r21)) * coeffNeg);
        }

        if (mUseLimitsRotational) {

            double phi21 = angleRad(mu2_gcs, mu1_gcs, ms2_gcs);

            updateConstraintAng(mcLimit3,     ms2_gcs,        ms2_gcs * -1.0,
                                              (phi21 - mPhiMin) * coeffNeg );

            updateConstraintAng(mcLimit4,     ms2_gcs * -1.0, ms2_gcs,
                                              (mPhiMax - phi21) * coeffNeg );

        }

        if (mUseFrictionOrMotorLengthwise) {

            updateConstraint(mcFric1, ms2_gcs,      Ra1s2* 1.0, 
                                      ms2_gcs*-1.0, Ra2s2*-1.0  );
        }

        if (mUseFrictionOrMotorRotational) {

            updateConstraintAng(mcFric2, ms2_gcs, ms2_gcs * -1.0);

        }

        if (mUseVelocityLimitLengthwise) {

            updateConstraint(mcVelLimit1, ms2_gcs,        Ra1s2 *  1.0, 
                                          ms2_gcs * -1.0, Ra2s2 * -1.0   );

            updateConstraint(mcVelLimit2, ms2_gcs * -1.0, Ra1s2 * -1.0,
                                          ms2_gcs,        Ra2s2 *  1.0 );
        }

        if (mUseVelocityLimitRotational) {

            updateConstraintAng(mcVelLimit3, ms2_gcs,        ms2_gcs * -1.0 );
            updateConstraintAng(mcVelLimit4, ms2_gcs * -1.0, ms2_gcs);

        }

    }

    virtual ~PistonJoint() {;}


  private:

    const Vec3          ma1;
    const Vec3          ma2;
    Vec3                ms1;
    Vec3                mu1;
    Vec3                ms2;
    Vec3                mu2;

    Vec3                ma1_gcs;
    Vec3                ma2_gcs;
    Vec3                mr1_gcs;
    Vec3                mr2_gcs;
    Vec3                ms1_gcs;
    Vec3                ms2_gcs;
    Vec3                mu1_gcs;
    Vec3                mu2_gcs;
    Mat3x3              mRa1;
    Mat3x3              mRa2;

    const bool          mUseLimitsLengthwise;
    const double        mdMin;
    const double        mdMax;
    const bool          mUseFrictionOrMotorLengthwise;
    const double        mwDesiredLengthwise;
    const double        mLambdaMaxLengthwise;
    const bool          mUseVelocityLimitLengthwise;
    const double        mwLimitLengthwise;
    const bool          mUseLimitsRotational;
    const double        mPhiMin;
    const double        mPhiMax;
    const bool          mUseFrictionOrMotorRotational;
    const double        mwDesiredRotational;
    const double        mLambdaMaxRotational;
    const bool          mUseVelocityLimitRotational;
    const double        mwLimitRotational;
    const double        mkCoeff;

    JacobianConstraint* mcLink1;
    JacobianConstraint* mcLink2;
    JacobianConstraint* mcLink3;
    JacobianConstraint* mcLink4;
    JacobianConstraint* mcLimit1;
    JacobianConstraint* mcLimit2;
    JacobianConstraint* mcLimit3;
    JacobianConstraint* mcLimit4;
    JacobianConstraint* mcFric1;
    JacobianConstraint* mcFric2;
    JacobianConstraint* mcVelLimit1;
    JacobianConstraint* mcVelLimit2;
    JacobianConstraint* mcVelLimit3;
    JacobianConstraint* mcVelLimit4;

    friend class JointManager;

};


class SliderJoint : public Joint {

  public:
    /** @brief constructor
     *
     *  @param   body1        (in): ConvexRigidBody 1.
     *
     *  @param   body2        (in): ConvexRigidBody 2.
     *
     *  @param   a1          (in): Anchor point on Body1 in LCS.
     *
     *  @param   a2          (in): Anchor point on Body2 in LCS.
     *
     *  @param   q1          (in): Initial orientation of ConvexRigidBody1 in 
     *                             GCS. It must be (1,0,0,0) if 
     *                             ConvexRigidBody1 is null.
     *                            
     *  @param   q2          (in): Initial orientation of ConvexRigidBody2 in 
     *                             GCS. It must be (1,0,0,0) if 
     *                             ConvexRigidBody2 is null.
     *
     *  @param   s1          (in): unit vector for 
     *                             Sliding direction in Body1's LCS, or in GCS
     *                             if ConvexRigidBody1 is null.
     *  @param   u1          (in): unit direction vector perpendicular to s1
     *                             in Body1's LCS, or in GCS if 
     *                             ConvexRigidBody1 is null.
     *
     *  @param   s2          (in): unit vector for 
     *                             Sliding direction in Body2's LCS, or in GCS
     *                             if ConvexRigidBody2 is null.
     *  @param   u2          (in): unit direction vector perpendicular to s2
     *                             in Body2's LCS, or in GCS if 
     *                             ConvexRigidBody2 is null.
     *
     *  @param   dMin        (in): The lowest point along s1 (=s2 in GCS).
     *
     *  @param   dMax        (in): The highest point along s1 (=s2 in GCS).
     *
     *  @param   useFrictionOrMotor 
     *                       (in): Set to true if you want to activate
     *                             friction or motor.
     *
     *  @param   wDesired    (in): the desired angular speed around s1 (=s2).
     *                             Must be zero if you want to implement
     *                             friction, rather than a motor.
     *
     *  @param   lambdaMax   (in): The range [-lamdaMax, lambdaMax] is used
     *                             for the boxed constraints to implement
     *                             friction or motor. 
     *                             This must be calculated by the application
     *                             based on the viscosity of the joint or the
     *                             maximum torque of the motor.
     *
     *  @param   useVelocityLimit
     *                       (in): Set to true if you want to use velocity 
     *                             limit. The relative angular velocity around
     *                             s1 is restricted in [-wLimit, wLimit].
     *
     *  @param   wLimit      (in): Maximum velocity around around s1.
     *                             Used if useVelocityLimit is set to true.
     *
     *  @param   kCoeff      (in): The coefficient to be multiplied to the RHS
     *                             of the linkage constraints to avoid over
     *                             compensation by numerical error.
     *                             Usually 0.8 or 0.9.
     */
    SliderJoint(

        ConvexRigidBody*             body1,
        ConvexRigidBody*             body2,
        const Vec3&       a1,
        const Vec3&       a2,
        const Vec3&       s1,
        const Vec3&       u1,
        const Vec3&       s2,
        const Vec3&       u2,
        const bool&       useLimits,
        const double&     dMin,
        const double&     dMax,
        const bool&       useFrictionOrMotor,
        const double&     wDesired,
        const double&     lambdaMax,
        const bool&       useVelocityLimit,
        const double&     wLimit,
        const double&     kCoeff
    ):
        Joint(body1, body2),
        ma1(a1),
        ma2(a2),
        mq1(1.0, 0.0, 0.0, 0.0),
        mq2(1.0, 0.0, 0.0, 0.0),
        mq12cInit(1.0, 0.0, 0.0, 0.0),
        ms1(s1),
        mu1(u1),
        ms2(s2),
        mu2(u2),
        mUseLimits(useLimits),
        mdMin(dMin),
        mdMax(dMax),
        mUseFrictionOrMotor(useFrictionOrMotor),
        mwDesired(wDesired),
        mLambdaMax(lambdaMax),
        mUseVelocityLimit(useVelocityLimit),
        mwLimit(wLimit),       
        mkCoeff(kCoeff),
        mcLink1(nullptr),
        mcLink2(nullptr),
        mcLink3(nullptr),
        mcLink4(nullptr),
        mcLink5(nullptr),
        mcLimit1(nullptr),
        mcLimit2(nullptr),
        mcFric1(nullptr),
        mcVelLimit1(nullptr),
        mcVelLimit2(nullptr)
    {
        ms1.normalize();
        mu1.normalize();
        ms2.normalize();
        mu2.normalize();

        if (body1!=nullptr) {
            mq1 = body1->Q();
        }
        if (body2!=nullptr) {
            mq2 = body2->Q();
        }
        mq12cInit = mq1 * mq2.conjugate();
          
        mcLink1 = newBilateralFree();
        mcLink2 = newBilateralFree();
        mcLink3 = newBilateralFree();
        mcLink4 = newBilateralFree();
        mcLink5 = newBilateralFree();

        if (!isBody1Fixed()) {
            mcLink3->setLHSang1(  1.0,  0.0,  0.0);
            mcLink4->setLHSang1(  0.0,  1.0,  0.0);
            mcLink5->setLHSang1(  0.0,  0.0,  1.0);
        }

        if (!isBody2Fixed()) {
            mcLink3->setLHSang2( -1.0,  0.0,  0.0);
            mcLink4->setLHSang2(  0.0, -1.0,  0.0);
            mcLink5->setLHSang2(  0.0,  0.0, -1.0);
        }

        if (mUseLimits) {

            mcLimit1 = newUnilateral();
            mcLimit2 = newUnilateral();

        }

        if (mUseFrictionOrMotor) {

            mcFric1 = newBilateralBox();
            mcFric1->setRHS(mwDesired);
            mcFric1->setBoxLimit(mLambdaMax);
        }

        if (mUseVelocityLimit) {

            mcVelLimit1 = newUnilateral();
            mcVelLimit1->setRHS(-1.0 * mwLimit);
            mcVelLimit2 = newUnilateral();
            mcVelLimit2->setRHS(-1.0 * mwLimit);
        }

        if (isBody1Fixed()) {

            ma1_gcs = ma1;
            mr1_gcs = ma1;
            ms1_gcs = ms1;
            mu1_gcs = mu1;
            mRa1.zero();
        }

        if (isBody2Fixed()) {

            ma2_gcs = ma2;
            mr2_gcs = ma2;
            ms2_gcs = ms2;
            mu2_gcs = mu2;
            mRa2.zero();
        }
    }


    void update(const double deltaTinv) {

        double coeffNeg = -1.0 * mkCoeff * deltaTinv;

        if (!isBody1Fixed()) {

            ma1_gcs = mBody1->Qmat() * ma1;
            mr1_gcs = mBody1->CoM() + ma1_gcs;
            ms1_gcs = mBody1->Qmat() * ms1;
            mu1_gcs = mBody1->Qmat() * mu1;
            mRa1    = ma1_gcs.crossMat();
        }

        if (!isBody2Fixed()) {

            ma2_gcs = mBody2->Qmat() * ma2;
            mr2_gcs = mBody2->CoM() + ma2_gcs;
            ms2_gcs = mBody2->Qmat() * ms2;
            mu2_gcs = mBody2->Qmat() * mu2;
            mRa2    = ma2_gcs.crossMat();
        }

        Vec3 r21   = mr1_gcs - mr2_gcs;

        Vec3 Ra1u2 = mRa1 * mu2_gcs;
        Vec3 Ra2u2 = mRa2 * mu2_gcs;

        updateConstraint(mcLink1, mu2_gcs *  1.0,  Ra1u2   *  1.0,
                                  mu2_gcs * -1.0,  Ra2u2   * -1.0,
                                  mu2_gcs.dot(r21) * coeffNeg  );

        auto w2_gcs = ms2_gcs.cross(mu2_gcs);
        Vec3 Ra1w2 = mRa1 * w2_gcs;
        Vec3 Ra2w2 = mRa2 * w2_gcs;

        updateConstraint(mcLink2, w2_gcs *  1.0,  Ra1w2   *  1.0,
                                  w2_gcs * -1.0,  Ra2w2   * -1.0,
                                  w2_gcs.dot(r21) * coeffNeg  );

        if (!isBody1Fixed()) {
            mq1 = mBody1->Q();
        }
        if (!isBody2Fixed()) {
            mq2 = mBody2->Q();
        }

        Quaternion q12c   = mq1 * mq2.conjugate();
        Quaternion qDev12 = q12c * mq12cInit.conjugate();
        Vec3 RHSori(qDev12.x(), qDev12.y(), qDev12.z());
        if (qDev12.s() < 0.0) {
            RHSori.scale(-2.0 * coeffNeg);
        }
        else {
            RHSori.scale( 2.0 * coeffNeg);
        }

        mcLink3->setRHS( RHSori.x() );
        mcLink4->setRHS( RHSori.y() );
        mcLink5->setRHS( RHSori.z() );

        Vec3 Ra1s2 = mRa1 * ms2_gcs;
        Vec3 Ra2s2 = mRa2 * ms2_gcs;

        if (mUseLimits) {

            updateConstraint(mcLimit1, ms2_gcs,        Ra1s2 *  1.0,
                                       ms2_gcs * -1.0, Ra2s2 * -1.0,
                                       (ms2_gcs.dot(r21) - mdMin) * coeffNeg);

            updateConstraint(mcLimit2, ms2_gcs * -1.0, Ra1s2 * -1.0,
                                       ms2_gcs,        Ra2s2 *  1.0,
                                       (mdMax - ms2_gcs.dot(r21)) * coeffNeg);
        }

        if (mUseFrictionOrMotor) {

            updateConstraint(mcFric1, ms2_gcs,        Ra1s2 *  1.0, 
                                      ms2_gcs * -1.0, Ra2s2 * -1.0    );
        }

        if (mUseVelocityLimit) {

            updateConstraint(mcVelLimit1, ms2_gcs,        Ra1s2 *  1.0,
                                          ms2_gcs * -1.0, Ra2s2 * -1.0  );

            updateConstraint(mcVelLimit2, ms2_gcs * -1.0, Ra1s2 * -1.0,
                                          ms2_gcs,        Ra2s2 *  1.0 );
        }

    }


    virtual ~SliderJoint() {;}


  private:

    const Vec3          ma1;
    const Vec3          ma2;
    Quaternion          mq1;
    Quaternion          mq2;
    Quaternion          mq12cInit;
    Vec3                ms1;
    Vec3                mu1;
    Vec3                ms2;
    Vec3                mu2;

    Vec3                ma1_gcs;
    Vec3                ma2_gcs;
    Vec3                mr1_gcs;
    Vec3                mr2_gcs;
    Vec3                ms1_gcs;
    Vec3                ms2_gcs;
    Vec3                mu1_gcs;
    Vec3                mu2_gcs;
    Mat3x3              mRa1;
    Mat3x3              mRa2;

    const bool          mUseLimits;
    const double        mdMin;
    const double        mdMax;
    const bool          mUseFrictionOrMotor;
    const double        mwDesired;
    const double        mLambdaMax;
    const bool          mUseVelocityLimit;
    const double        mwLimit;
    const double        mkCoeff;

    JacobianConstraint* mcLink1;
    JacobianConstraint* mcLink2;
    JacobianConstraint* mcLink3;
    JacobianConstraint* mcLink4;
    JacobianConstraint* mcLink5;
    JacobianConstraint* mcLimit1;
    JacobianConstraint* mcLimit2;
    JacobianConstraint* mcFric1;
    JacobianConstraint* mcVelLimit1;
    JacobianConstraint* mcVelLimit2;

    friend class JointManager;
};


class FixedJoint  : public Joint {

  public:
    /** @brief constructor
     *
     *  @param   body1        (in): ConvexRigidBody 1. Can not be nullptr.
     *
     *  @param   body2        (in): ConvexRigidBody 2. Can not be nullptr.
     *
     *  @param   a1          (in): Anchor point on Body1 in LCS.
     *
     *  @param   a2          (in): Anchor point on Body2 in LCS.
     *
     *  @param   q1          (in): Initial orientation of ConvexRigidBody1 in
     *                             GCS.
     *
     *  @param   q2          (in): Initial orientation of ConvexRigidBody2 in
     *                             GCS.
     *
     *  @param   kCoeff      (in): The coefficient to be multiplied to the RHS
     *                             of the linkage constraints to avoid over
     *                             compensation by numerical error.
     *                             Usually 0.8 or 0.9.
     */
    FixedJoint(

        ConvexRigidBody*             body1,
        ConvexRigidBody*             body2,
        const Vec3&       a1,
        const Vec3&       a2,
        const Quaternion& q1,
        const Quaternion& q2,
        const double      kCoeff

    ):
        Joint(body1, body2),
        ma1(a1),
        ma2(a2),
        mq12cInit(q1*q2.conjugate()),
        mkCoeff(kCoeff)
    {

        mcLink1 = newBilateralFree();
        mcLink2 = newBilateralFree();
        mcLink3 = newBilateralFree();
        mcLink4 = newBilateralFree();
        mcLink5 = newBilateralFree();
        mcLink6 = newBilateralFree();

        mcLink1->setLHSlin1(         1.0,          0.0,          0.0);
        mcLink2->setLHSlin1(         0.0,          1.0,          0.0);
        mcLink3->setLHSlin1(         0.0,          0.0,          1.0);
        mcLink4->setLHSang1(         1.0,          0.0,          0.0);
        mcLink5->setLHSang1(         0.0,          1.0,          0.0);
        mcLink6->setLHSang1(         0.0,          0.0,          1.0);

        mcLink1->setLHSlin2(        -1.0,          0.0,          0.0);
        mcLink2->setLHSlin2(         0.0,         -1.0,          0.0);
        mcLink3->setLHSlin2(         0.0,          0.0,         -1.0);
        mcLink4->setLHSang2(        -1.0,          0.0,          0.0);
        mcLink5->setLHSang2(         0.0,         -1.0,          0.0);
        mcLink6->setLHSang2(         0.0,          0.0,         -1.0);

    }


    void update(const double deltaTinv) {

        double coeffNeg = -1.0 * mkCoeff * deltaTinv;

        Vec3 Ra1 = mBody1->Qmat() * ma1;
        Vec3 r1  = mBody1->CoM() + Ra1;

        mcLink1->setLHSang1(         0.0,  1.0*Ra1.z(), -1.0*Ra1.y());
        mcLink2->setLHSang1(-1.0*Ra1.z(),          0.0,  1.0*Ra1.x());
        mcLink3->setLHSang1( 1.0*Ra1.y(), -1.0*Ra1.x(),          0.0);

        Vec3 Ra2 = mBody2->Qmat() * ma2;
        Vec3 r2  = mBody2->CoM() + Ra2;

        mcLink1->setLHSang2(         0.0, -1.0*Ra2.z(),      Ra2.y());
        mcLink2->setLHSang2(     Ra2.z(),          0.0, -1.0*Ra2.x());
        mcLink3->setLHSang2(-1.0*Ra2.y(),      Ra2.x(),          0.0);

        Vec3 RHSlin(r1 - r2);
        RHSlin.scale(coeffNeg);

        mcLink1->setRHS    ( RHSlin.x() );
        mcLink2->setRHS    ( RHSlin.y() );
        mcLink3->setRHS    ( RHSlin.z() );

        Quaternion q12c   = mBody1->Q() * (mBody2->Q()).conjugate();
        Quaternion qDev12 = q12c * mq12cInit.conjugate();

        Vec3 RHSori(qDev12.x(), qDev12.y(), qDev12.z());
        if (qDev12.s() < 0.0) {
            RHSori.scale(-2.0 * coeffNeg);
        }
        else {
            RHSori.scale( 2.0 * coeffNeg);
        }

        mcLink4->setRHS    ( RHSori.x() );
        mcLink5->setRHS    ( RHSori.y() );
        mcLink6->setRHS    ( RHSori.z() );

    }


    virtual ~FixedJoint() {;}


  private:

    const Vec3          ma1;
    const Vec3          ma2;
    Quaternion          mq12cInit;
    const double        mkCoeff;

    JacobianConstraint* mcLink1;
    JacobianConstraint* mcLink2;
    JacobianConstraint* mcLink3;
    JacobianConstraint* mcLink4;
    JacobianConstraint* mcLink5;
    JacobianConstraint* mcLink6;

    friend class JointManager;
};



class Hinge2Joint : public Joint {

  public:

    /** @brief constructor
     *
     *  @param   body1       (in): ConvexRigidBody 1. Null if this side is 
     *                              fixed.
     *
     *  @param   body2       (in): ConvexRigidBody 2. Null if this side is 
     *                              fixed.
     *
     *  @param   a1          (in): Anchor point on Body1 in LCS, or the anchor
     *                             point in GCS if Body1 is null.
     *
     *  @param   a2          (in): Anchor point on Body2 in LCS, or the anchor
     *                             point in GCS if Body2 is null.
     *
     *  @param   s1          (in): unit Axis direction vector on Body1 in LCS, 
     *                             or in GCS if Body1 is null.
     *
     *  @param   s2          (in): unit Axis direction vector on Body2 in LCS, 
     *                             or in GCS if Body2 is null.
     *                             s2 x s1 makes the axis of rotation if
     *                             you see this joint as a universal joint
     *                             between two rotational shafts.
     *
     *  @param   phiHinge    (in): The angle that s1 and s2 makes around the
     *                             two rotational shafts. 0 < phiHinge < PI.
     *                             It is usually 0.5*PI, i.e., s1 and s2 meet
     *                             at the right angle.
     *
     *  @param   u1          (in): unit direction vector on Body1 in LCS,
     *                             or in GCS if Body1 is null.
     *                             This is comopared with (s2 x s1) to indicate
     *                             the relative orientation of body1 along s1.
     *
     *                             u1 and s2 x s1 are used to indicate the 
     *                             relative orientation of two bodies around
     *                             s1 such that asin ((u1 x (s2 x s1))·s1)
     *                             indicates the directed angle  ϕ_21,
     *                             which is compared with phiMin1 and phiMax1.
     *
     *  @param   u2          (in): unit direction vector on Body2 in LCS, 
     *                             or in GCS if Body2 is null.
     *                             This is comopared with -(s2 x s1) to 
     *                             indicate the relative orientation of body2
     *                             along s2.
     *
     *                             u2 and -(s2 x s1) are used to indicate the 
     *                             relative orientation of two bodies around
     *                             s1 such that asin ((u2 x -(s2 x s1))·s2)
     *                             indicates the directed angle  ϕ_12,
     *                             which is compared with phiMin2 and phiMax2.
     *
     *  @param   useLimits   (in): Set to true if you want to use the angular 
     *                             limit.
     *
     *  @param   phiMin1     (in): Minimum angle between u1 and (s2 x s1).
     *                            Used if useLimits is set to true.
     *
     *  @param   phiMax1     (in): Maximum angle between u1 and (s2 x s1).
     *                             Used if useLimits is set to true.
     *
     *  @param   phiMin2     (in): Minimum angle between u2 and -(s2 x s1).
     *                             Used if useLimits is set to true.
     *
     *  @param   phiMax2     (in): Maximum angle between u2 and -(s2 x s1).
     *                             Used if useLimits is set to true.
     *
     *  @param   useFriction (in): Set to true if you want to activate
     *                             friction
     *
     *  @param   lambdaMax   (in): The range [-lamdaMax, lambdaMax] is used
     *                             for the boxed constraints to implement
     *                             friction or motor. 
     *                             This must be calculated by the application
     *                             based on the viscosity of the joint or the
     *                             maximum torque of the motor.
     *
     *  @param   useVelocityLimit
     *                       (in): Set to true if you want to use velocity 
     *                             limit. The relative angular velocity around
     *                             s1 is restricted in [-wLimit, wLimit].
     *
     *  @param   wLimit      (in): Maximum velocity around around s1.
     *                             Used if useVelocityLimit is set to true.
     *
     *  @param   kCoeff      (in): The coefficient to be multiplied to the RHS
     *                             of the linkage constraints to avoid over
     *                             compensation by numerical error.
     *                             Usually 0.8 or 0.9.
     */
    Hinge2Joint(

        ConvexRigidBody*        body1,
        ConvexRigidBody*        body2,
        const Vec3&  a1,
        const Vec3&  a2,
        const Vec3&  s1,
        const Vec3&  s2,
        const double phiHinge,
        const Vec3&  u1,
        const Vec3&  u2,
        const bool   useLimits,
        const double phiMin1,
        const double phiMax1,
        const double phiMin2,
        const double phiMax2,
        const bool   useFriction,
        const double lambdaMax,
        const bool   useVelocityLimit,
        const double wLimit,
        const double kCoeff

    ):
        Joint(body1, body2),
        ma1(a1),
        ma2(a2),
        ms1(s1),
        ms2(s2),
        mPhiHinge(phiHinge),
        mu1(u1),
        mu2(u2),
        mUseLimits(useLimits),
        mPhiMin1(phiMin1),
        mPhiMax1(phiMax1),
        mPhiMin2(phiMin2),
        mPhiMax2(phiMax2),
        mUseFriction(useFriction),
        mLambdaMax(lambdaMax),
        mUseVelocityLimit(useVelocityLimit),
        mwLimit(wLimit),
        mkCoeff(kCoeff),
        mcLink1(nullptr),
        mcLink2(nullptr),
        mcLink3(nullptr),
        mcLink4(nullptr),
        mcLimit1(nullptr),
        mcLimit2(nullptr),
        mcLimit3(nullptr),
        mcLimit4(nullptr),
        mcFric1(nullptr),
        mcFric2(nullptr),
        mcVelLimit1(nullptr),
        mcVelLimit2(nullptr),
        mcVelLimit3(nullptr),
        mcVelLimit4(nullptr)
    {
        ms1.normalize();
        ms2.normalize();
        mu1.normalize();
        mu2.normalize();

        mcLink1 = newBilateralFree();
        mcLink2 = newBilateralFree();
        mcLink3 = newBilateralFree();
        mcLink4 = newBilateralFree();

        if (!isBody1Fixed()) {
            mcLink1->setLHSlin1(         1.0,          0.0,          0.0);
            mcLink2->setLHSlin1(         0.0,          1.0,          0.0);
            mcLink3->setLHSlin1(         0.0,          0.0,          1.0);
        }

        if (!isBody2Fixed()) {
            mcLink1->setLHSlin2(        -1.0,          0.0,          0.0);
            mcLink2->setLHSlin2(         0.0,         -1.0,          0.0);
            mcLink3->setLHSlin2(         0.0,          0.0,         -1.0);
        }

        if (mUseLimits) {

            mcLimit1 = newUnilateral();
            mcLimit2 = newUnilateral();
            mcLimit3 = newUnilateral();
            mcLimit4 = newUnilateral();

        }

        if (mUseFriction) {

            mcFric1 = newBilateralBox();
            mcFric1->setBoxLimit(mLambdaMax);

            mcFric2 = newBilateralBox();
            mcFric2->setBoxLimit(mLambdaMax);

        }

        if (mUseVelocityLimit) {

            mcVelLimit1 = newUnilateral();
            mcVelLimit2 = newUnilateral();
            mcVelLimit3 = newUnilateral();
            mcVelLimit4 = newUnilateral();

            mcVelLimit1->setRHS(-1.0 * mwLimit);
            mcVelLimit2->setRHS(-1.0 * mwLimit);
            mcVelLimit3->setRHS(-1.0 * mwLimit);
            mcVelLimit4->setRHS(-1.0 * mwLimit);
        }

        if (isBody1Fixed()) {

            ma1_gcs = ma1;
            mr1_gcs = ma1;
            ms1_gcs = ms1;
            mu1_gcs = mu1;
            mRa1.zero();
        }

        if (isBody2Fixed()) {

            ma2_gcs = ma2;
            mr2_gcs = ma2;
            ms2_gcs = ms2;
            mu2_gcs = mu2;
            mRa2.zero();

        }
    }
    

    void update(const double deltaTinv) {

        double coeffNeg = -1.0 * mkCoeff * deltaTinv;

        if (!isBody1Fixed()) {

            ma1_gcs = mBody1->Qmat() * ma1;
            mr1_gcs = mBody1->CoM() + ma1_gcs;
            ms1_gcs = mBody1->Qmat() * ms1;
            mu1_gcs = mBody1->Qmat() * mu1;
            mRa1    = ma1_gcs.crossMat();
        }

        if (!isBody2Fixed()) {

            ma2_gcs = mBody2->Qmat() * ma2;
            mr2_gcs = mBody2->CoM() + ma2_gcs;
            ms2_gcs = mBody2->Qmat() * ms2;
            mu2_gcs = mBody2->Qmat() * mu2;
            mRa2    = ma2_gcs.crossMat();
        }

        Vec3 RHSlin(mr1_gcs - mr2_gcs);
        RHSlin.scale(coeffNeg);

        if (!isBody1Fixed()) {
            auto Ra1n = mRa1;
            Ra1n.scale(-1.0);
            mcLink1->setLHSang1(Ra1n.cell(1,1),Ra1n.cell(1,2),Ra1n.cell(1,3));
            mcLink2->setLHSang1(Ra1n.cell(2,1),Ra1n.cell(2,2),Ra1n.cell(2,3));
            mcLink3->setLHSang1(Ra1n.cell(3,1),Ra1n.cell(3,2),Ra1n.cell(3,3));
        }

        if (!isBody2Fixed()) {
            mcLink1->setLHSang2(mRa2.cell(1,1),mRa2.cell(1,2),mRa2.cell(1,3));
            mcLink2->setLHSang2(mRa2.cell(2,1),mRa2.cell(2,2),mRa2.cell(2,3));
            mcLink3->setLHSang2(mRa2.cell(3,1),mRa2.cell(3,2),mRa2.cell(3,3));
        }

        mcLink1->setRHS    ( RHSlin.x() );
        mcLink2->setRHS    ( RHSlin.y() );
        mcLink3->setRHS    ( RHSlin.z() );

        Vec3   s2_cr_s1 = ms2_gcs.cross(ms1_gcs);
        double norm2    = s2_cr_s1.norm2();
        double theta    = asinSafe(norm2);
        s2_cr_s1.scale(1.0/norm2);
        double s1_dot_s2 = ms1_gcs.dot(ms2_gcs);
        double phi21link = (s1_dot_s2>=0.0)?theta:(M_PI - theta);

        updateConstraintAng(mcLink4, s2_cr_s1,  s2_cr_s1 * -1.0, 
                                   -1.0*(mPhiHinge - phi21link) * coeffNeg);

        if (mUseLimits) {

            double phi21 = angleRad(s2_cr_s1, mu1_gcs, ms1_gcs);
            
            updateConstraintAng(mcLimit1, ms1_gcs,         ms1_gcs * -1.0, 
                                          (phi21 - mPhiMin1) * coeffNeg) ;

            updateConstraintAng(mcLimit2, ms1_gcs * -1.0,  ms1_gcs, 
                                          (mPhiMax1 - phi21) * coeffNeg) ;

            double phi12 = angleRad(s2_cr_s1* -1.0, mu2_gcs, ms2_gcs);

            updateConstraintAng(mcLimit3, ms2_gcs * -1.0,  ms2_gcs,
                                          (phi12 - mPhiMin2) * coeffNeg) ;

            updateConstraintAng(mcLimit4, ms2_gcs,         ms2_gcs * -1.0, 
                                          (mPhiMax2 - phi12) * coeffNeg) ;
        }
        
        if (mUseFriction) {

            updateConstraintAng(mcFric1, ms2_gcs, ms2_gcs * -1.0);
            updateConstraintAng(mcFric2, ms1_gcs, ms1_gcs * -1.0);

        }

        if (mUseVelocityLimit) {

            updateConstraintAng(mcVelLimit1, ms1_gcs,        ms1_gcs * -1.0);
            updateConstraintAng(mcVelLimit2, ms1_gcs * -1.0, ms1_gcs       );
            updateConstraintAng(mcVelLimit3, ms2_gcs,        ms2_gcs * -1.0);
            updateConstraintAng(mcVelLimit4, ms2_gcs * -1.0, ms2_gcs       );

        }
    }

    virtual ~Hinge2Joint() {;}


  private:

    const Vec3          ma1;
    const Vec3          ma2;
    Vec3                ms1;
    Vec3                ms2;
    const double        mPhiHinge;
    Vec3                mu1;
    Vec3                mu2;

    Vec3                ma1_gcs;
    Vec3                ma2_gcs;
    Vec3                mr1_gcs;
    Vec3                mr2_gcs;
    Vec3                ms1_gcs;
    Vec3                ms2_gcs;
    Vec3                mu1_gcs;
    Vec3                mu2_gcs;
    Mat3x3              mRa1;
    Mat3x3              mRa2;

    const bool          mUseLimits;
    const double        mPhiMin1;
    const double        mPhiMax1;
    const double        mPhiMin2;
    const double        mPhiMax2;
    const bool          mUseFriction;
    const double        mLambdaMax;
    const double        mUseVelocityLimit;
    const double        mwLimit;
    const double        mkCoeff;

    JacobianConstraint* mcLink1;
    JacobianConstraint* mcLink2;
    JacobianConstraint* mcLink3;
    JacobianConstraint* mcLink4;
    JacobianConstraint* mcLimit1;
    JacobianConstraint* mcLimit2;
    JacobianConstraint* mcLimit3;
    JacobianConstraint* mcLimit4;
    JacobianConstraint* mcFric1;
    JacobianConstraint* mcFric2;
    JacobianConstraint* mcVelLimit1;
    JacobianConstraint* mcVelLimit2;
    JacobianConstraint* mcVelLimit3;
    JacobianConstraint* mcVelLimit4;

    friend class JointManager;
};


class Hinge1Joint : public Joint {

  public:


    /** @brief constructor
     *
     *  @param   body1       (in): ConvexRigidBody 1. Null if this side is 
     *                             fixed.
     *
     *  @param   body2       (in): ConvexRigidBody 2. Null if this side is 
     *                             fixed.
     *
     *  @param   a1          (in): Anchor point on Body1 in LCS, or the anchor
     *                             point in GCS if Body1 is null.
     *
     *  @param   a2          (in): Anchor point on Body2 in LCS, or the anchor
     *                             point in GCS if Body2 is null.
     *
     *  @param   s1          (in): unit Axis direction vector on Body1 in LCS, 
     *                             or in GCS if Body1 is null.
     *
     *  @param   s2          (in): unit Axis direction vector on Body2 in LCS, 
     *                             or in GCS if Body2 is null.
     *                             s1 is supposed to be aligned with s2 in GCS.
     *
     *  @param   u1          (in): unit direction vector on Body1 in LCS,
     *                             or in GCS if Body1 is null.
     *
     *  @param   u2          (in): unit direction vector on Body2 in LCS, 
     *                             or in GCS if Body2 is null.
     *                             u1 and u2 are used to indicate the relative
     *                             orientation of two bodies around s1 (=s2)
     *                             such that asin((u2 x u1)·s1) indicates the
     *                             directed angle ϕ_21, wchih is compared 
     *                             with phiMin and phiMax.
     *  
     *  @param   useLimits   (in): Set to true if you want to use the angular 
     *                             limit.
     *
     *  @param   phiMin      (in): Minimum angle between u1 and u2.
     *                             Used if useLimits is set to true.
     *
     *  @param   phiMax      (in): Maximum angle between u1 and u2.
     *                             Used if useLimits is set to true.
     *
     *  @param   useFrictionOrMotor 
     *                       (in): Set to true if you want to activate
     *                             friction or motor.
     *
     *  @param   wDesired    (in): the desired angular speed around s1 (=s2).
     *                             Must be zero if you want to implement
     *                             friction, rather than a motor.
     *
     *  @param   lambdaMax   (in): The range [-lamdaMax, lambdaMax] is used
     *                             for the boxed constraints to implement
     *                             friction or motor. 
     *                             This must be calculated by the application
     *                             based on the viscosity of the joint or the
     *                             maximum torque of the motor.
     *
     *  @param   useVelocityLimit
     *                       (in): Set to true if you want to use velocity 
     *                             limit. The relative angular velocity around
     *                             s1 is restricted in [-wLimit, wLimit].
     *
     *  @param   wLimit      (in): Maximum velocity around around s1.
     *                             Used if useVelocityLimit is set to true.
     *
     *  @param   kCoeff      (in): The coefficient to be multiplied to the RHS
     *                             of the linkage constraints to avoid over
     *                             compensation by numerical error.
     *                             Usually 0.8 or 0.9.
     */
    Hinge1Joint(

        ConvexRigidBody*        body1,
        ConvexRigidBody*        body2,
        const Vec3&  a1,
        const Vec3&  a2,
        const Vec3&  s1,
        const Vec3&  s2,
        const Vec3&  u1,
        const Vec3&  u2,
        const bool   useLimits,
        const double phiMin,
        const double phiMax,
        const bool   useFrictionOrMotor,
        const double wDesired,
        const double lambdaMax,
        const bool   useVelocityLimit,
        const double wLimit,
        const double kCoeff

    ): 
        Joint(body1, body2),
        ma1(a1),
        ma2(a2),
        ms1(s1),
        ms2(s2),
        mu1(u1),
        mu2(u2),
        mUseLimits(useLimits),
        mPhiMin(phiMin),
        mPhiMax(phiMax),
        mUseFrictionOrMotor(useFrictionOrMotor),
        mwDesired(wDesired),
        mLambdaMax(lambdaMax),
        mUseVelocityLimit(useVelocityLimit),
        mwLimit(wLimit),
        mkCoeff(kCoeff),
        mcLink1(nullptr),
        mcLink2(nullptr),
        mcLink3(nullptr),
        mcLink4(nullptr),
        mcLink5(nullptr),
        mcLimit1(nullptr),
        mcLimit2(nullptr),
        mcFric1(nullptr),
        mcVelLimit1(nullptr),
        mcVelLimit2(nullptr)
    {
        ms1.normalize();
        ms2.normalize();
        mu1.normalize();
        mu2.normalize();

        mcLink1 = newBilateralFree();
        mcLink2 = newBilateralFree();
        mcLink3 = newBilateralFree();
        mcLink4 = newBilateralFree();
        mcLink5 = newBilateralFree();

        if (!isBody1Fixed()) {
            mcLink1->setLHSlin1(         1.0,          0.0,          0.0);
            mcLink2->setLHSlin1(         0.0,          1.0,          0.0);
            mcLink3->setLHSlin1(         0.0,          0.0,          1.0);
        }

        if (!isBody2Fixed()) {
            mcLink1->setLHSlin2(        -1.0,          0.0,          0.0);
            mcLink2->setLHSlin2(         0.0,         -1.0,          0.0);
            mcLink3->setLHSlin2(         0.0,          0.0,         -1.0);
        }

        if (mUseLimits) {
            mcLimit1 = newUnilateral();
            mcLimit2 = newUnilateral();
        }

        if (mUseFrictionOrMotor) {
            mcFric1 = newBilateralBox();
            mcFric1->setRHS(mwDesired);
            mcFric1->setBoxLimit(mLambdaMax);
        }

        if (mUseVelocityLimit) {
            mcVelLimit1 = newUnilateral();
            mcVelLimit2 = newUnilateral();
            mcVelLimit1->setRHS(-1.0 * mwLimit);
            mcVelLimit2->setRHS(-1.0 * mwLimit);
        }

        if (isBody1Fixed()) {

            ma1_gcs = ma1;
            mr1_gcs = ma1;
            ms1_gcs = ms1;
            mu1_gcs = mu1;
            mRa1.zero();
        }

        if (isBody2Fixed()) {

            ma2_gcs = ma2;
            mr2_gcs = ma2;
            ms2_gcs = ms2;
            mu2_gcs = mu2;
            mRa2.zero();
        }
    }


    virtual ~Hinge1Joint() {;}


    void update(const double deltaTinv) {


        double coeffNeg = -1.0 * mkCoeff * deltaTinv;

        if (!isBody1Fixed()) {

            ma1_gcs = mBody1->Qmat() * ma1;
            mr1_gcs = mBody1->CoM() + ma1_gcs;
            ms1_gcs = mBody1->Qmat() * ms1;
            mu1_gcs = mBody1->Qmat() * mu1;
            mRa1    = ma1_gcs.crossMat();
        }

        if (!isBody2Fixed()) {

            ma2_gcs = mBody2->Qmat() * ma2;
            mr2_gcs = mBody2->CoM() + ma2_gcs;
            ms2_gcs = mBody2->Qmat() * ms2;
            mu2_gcs = mBody2->Qmat() * mu2;
            mRa2    = ma2_gcs.crossMat();
        }

        Vec3 RHSlin(mr1_gcs - mr2_gcs);
        RHSlin.scale(coeffNeg);

        if (!isBody1Fixed()) {
            auto Ra1n = mRa1;
            Ra1n.scale(-1.0);
            mcLink1->setLHSang1(Ra1n.cell(1,1),Ra1n.cell(1,2),Ra1n.cell(1,3));
            mcLink2->setLHSang1(Ra1n.cell(2,1),Ra1n.cell(2,2),Ra1n.cell(2,3));
            mcLink3->setLHSang1(Ra1n.cell(3,1),Ra1n.cell(3,2),Ra1n.cell(3,3));
        }

        if (!isBody2Fixed()) {
            mcLink1->setLHSang2(mRa2.cell(1,1),mRa2.cell(1,2),mRa2.cell(1,3));
            mcLink2->setLHSang2(mRa2.cell(2,1),mRa2.cell(2,2),mRa2.cell(2,3));
            mcLink3->setLHSang2(mRa2.cell(3,1),mRa2.cell(3,2),mRa2.cell(3,3));
        }

        mcLink1->setRHS( RHSlin.x() );
        mcLink2->setRHS( RHSlin.y() );
        mcLink3->setRHS( RHSlin.z() );

        Vec3 s2_cr_s1 = ms2_gcs.cross(ms1_gcs);
        Vec3 w2_gcs = ms2_gcs.cross(mu2_gcs);
        double phi21u = angleRad(ms2_gcs, ms1_gcs, mu2_gcs);
        double phi21w = angleRad(ms2_gcs, ms1_gcs,  w2_gcs);

        updateConstraintAng(mcLink4, mu2_gcs, mu2_gcs * -1.0, phi21u*coeffNeg);
        updateConstraintAng(mcLink5,  w2_gcs,  w2_gcs * -1.0, phi21w*coeffNeg);

        if (mUseLimits) {
            double phi21 = angleRad(mu2_gcs, mu1_gcs, ms2_gcs);

            updateConstraintAng(mcLimit1, ms2_gcs,        ms2_gcs * -1.0, 
                                          (phi21 - mPhiMin) * coeffNeg);

            updateConstraintAng(mcLimit2, ms2_gcs * -1.0, ms2_gcs,
                                          (mPhiMax - phi21) * coeffNeg);
        }

        if (mUseFrictionOrMotor) {
            updateConstraintAng(mcFric1, ms2_gcs, ms2_gcs * -1.0);
        }

        if (mUseVelocityLimit) {
            updateConstraintAng(mcVelLimit1, ms2_gcs,        ms2_gcs * -1.0);
            updateConstraintAng(mcVelLimit2, ms2_gcs * -1.0, ms2_gcs       );
        }
    }

  private:

    const Vec3          ma1;
    const Vec3          ma2;
    Vec3                ms1;
    Vec3                ms2;
    Vec3                mu1;
    Vec3                mu2;

    Vec3                ma1_gcs;
    Vec3                ma2_gcs;
    Vec3                mr1_gcs;
    Vec3                mr2_gcs;
    Vec3                ms1_gcs;
    Vec3                ms2_gcs;
    Vec3                mu1_gcs;
    Vec3                mu2_gcs;
    Mat3x3              mRa1;
    Mat3x3              mRa2;

    const bool          mUseLimits;
    const double        mPhiMin;
    const double        mPhiMax;
    const bool          mUseFrictionOrMotor;
    const double        mwDesired;
    const double        mLambdaMax;
    const double        mUseVelocityLimit;
    const double        mwLimit;
    const double        mkCoeff;

    JacobianConstraint* mcLink1;
    JacobianConstraint* mcLink2;
    JacobianConstraint* mcLink3;
    JacobianConstraint* mcLink4;
    JacobianConstraint* mcLink5;
    JacobianConstraint* mcLimit1;
    JacobianConstraint* mcLimit2;
    JacobianConstraint* mcFric1;
    JacobianConstraint* mcVelLimit1;
    JacobianConstraint* mcVelLimit2;

    friend class JointManager;
};


class BallJoint : public Joint {

  public:

    /** @brief constructor
     *
     *  @param   body1       (in): ConvexRigidBody 1. Null if this side is 
     *                             fixed.
     *
     *  @param   body2       (in): ConvexRigidBody 2. Null if this side is
     *                             fixed.
     *
     *  @param   a1          (in): Anchor point on Body1 in LCS, or the anchor
     *                             point in GCS if Body1 is null.
     *
     *  @param   a2          (in): Anchor point on Body2 in LCS, or the anchor
     *                             point in GCS if Body2 is null.
     *
     *  @param   useLimits   (in): Set to true if you want to use the angular 
     *                             limit.
     *
     *  @param   s1          (in): Reference direction vector on Body1 in LCS, 
     *                             or in GCS if Body1 is null.
     *
     *  @param   s2          (in): Reference direction vector on Body2 in LCS, 
     *                             or in GCS if Body2 is null.
     *
     *  @param   phi         (in): Minimum angle between s1 and s2.
     *                             Used if useLimits is set to true.
     *
     *  @param   useFriction (in): Set to true if you want to activate
     *                             isometric rotational friction.
     *
     *  @param   lambdaMax   (in): The range [-lamdaMax, lambdaMax] is used
     *                             for the boxed constraints to implement
     *                             friction. This must be calculated by the 
     *                             application based on the viscosity of the
     *                             joint.
     *
     *  @param   kCoeff      (in): The coefficient to be multiplied to the RHS
     *                             of the linkage constraints to avoid over
     *                             compensation by numerical error.
     *                             Usually 0.8 or 0.9.
     */
    BallJoint(

        ConvexRigidBody*        body1,
        ConvexRigidBody*        body2,
        const Vec3&  a1,
        const Vec3&  a2,
        const bool   useLimits,
        const Vec3&  s1,
        const Vec3&  s2,
        const double phi,
        const bool   useFriction,
        const double lambdaMax,
        const double kCoeff

    ): 
        Joint(body1, body2),
        ma1(a1), 
        ma2(a2), 
        mUseLimits(useLimits),
        ms1(s1),
        ms2(s2),
        mPhi(phi),
        mUseFriction(useFriction),
        mLambdaMax(lambdaMax),
        mkCoeff(kCoeff),
        mcLink1(nullptr),
        mcLink2(nullptr),
        mcLink3(nullptr),
        mcLimit1(nullptr),
        mcFric1(nullptr),
        mcFric2(nullptr),
        mcFric3(nullptr)
    {

        mcLink1 = newBilateralFree();
        mcLink2 = newBilateralFree();
        mcLink3 = newBilateralFree();

        if (!isBody1Fixed()) {
            mcLink1->setLHSlin1(         1.0,          0.0,          0.0);
            mcLink2->setLHSlin1(         0.0,          1.0,          0.0);
            mcLink3->setLHSlin1(         0.0,          0.0,          1.0);
        }

        if (!isBody2Fixed()) {
            mcLink1->setLHSlin2(        -1.0,          0.0,          0.0);
            mcLink2->setLHSlin2(         0.0,         -1.0,          0.0);
            mcLink3->setLHSlin2(         0.0,          0.0,         -1.0);
        }

        if (mUseLimits) {
            mcLimit1 = newUnilateral();
        }

        if (mUseFriction) {

            mcFric1 = newBilateralBox();
            mcFric2 = newBilateralBox();
            mcFric3 = newBilateralBox();

            if (!isBody1Fixed()) {
                mcFric1->setLHSang1(   1.0,   0.0,   0.0);
                mcFric2->setLHSang1(   0.0,   1.0,   0.0);
                mcFric3->setLHSang1(   0.0,   0.0,   1.0);
            }

            if (!isBody2Fixed()) {
                mcFric1->setLHSang2(  -1.0,   0.0,   0.0);
                mcFric2->setLHSang2(   0.0,  -1.0,   0.0);
                mcFric3->setLHSang2(   0.0,   0.0,  -1.0);
            }
            mcFric1->setBoxLimit(mLambdaMax);
            mcFric2->setBoxLimit(mLambdaMax);
            mcFric3->setBoxLimit(mLambdaMax);
        }

        if (isBody1Fixed()) {

            ma1_gcs = ma1;
            mr1_gcs = ma1;
            ms1_gcs = ms1;
            mRa1.zero();
        }

        if (isBody2Fixed()) {

            ma2_gcs = ma2;
            mr2_gcs = ma2;
            ms2_gcs = ms2;
            mRa2.zero();
        }

    }

    virtual ~BallJoint() {;}


    void update(const double deltaTinv) {

        double coeffNeg = -1.0 * mkCoeff * deltaTinv;

        if (!isBody1Fixed()) {

            ma1_gcs = mBody1->Qmat() * ma1;
            mr1_gcs = mBody1->CoM() + ma1_gcs;
            ms1_gcs = mBody1->Qmat() * ms1;
            mRa1    = ma1_gcs.crossMat();
        }

        if (!isBody2Fixed()) {

            ma2_gcs = mBody2->Qmat() * ma2;
            mr2_gcs = mBody2->CoM() + ma2_gcs;
            ms2_gcs = mBody2->Qmat() * ms2;
            mRa2    = ma2_gcs.crossMat();
        }


        Vec3 RHSlin(mr1_gcs - mr2_gcs);
        RHSlin.scale(coeffNeg);

        if (!isBody1Fixed()) {
            auto Ra1n = mRa1;
            Ra1n.scale(-1.0);
            mcLink1->setLHSang1(Ra1n.cell(1,1),Ra1n.cell(1,2),Ra1n.cell(1,3));
            mcLink2->setLHSang1(Ra1n.cell(2,1),Ra1n.cell(2,2),Ra1n.cell(2,3));
            mcLink3->setLHSang1(Ra1n.cell(3,1),Ra1n.cell(3,2),Ra1n.cell(3,3));
        }

        if (!isBody2Fixed()) {
            mcLink1->setLHSang2(mRa2.cell(1,1),mRa2.cell(1,2),mRa2.cell(1,3));
            mcLink2->setLHSang2(mRa2.cell(2,1),mRa2.cell(2,2),mRa2.cell(2,3));
            mcLink3->setLHSang2(mRa2.cell(3,1),mRa2.cell(3,2),mRa2.cell(3,3));
        }

        mcLink1->setRHS( RHSlin.x() );
        mcLink2->setRHS( RHSlin.y() );
        mcLink3->setRHS( RHSlin.z() );

        if (mUseLimits) {

            Vec3   s2_cr_s1  = ms2_gcs.cross(ms1_gcs);
            s2_cr_s1.normalize();

            double phi21 = angleRad(ms2_gcs, ms1_gcs, s2_cr_s1);
            updateConstraintAng(mcLimit1, s2_cr_s1,  s2_cr_s1 * -1.0, 
                                                (phi21 - mPhi)*coeffNeg);
        }

    }

  private:

    const Vec3          ma1;
    const Vec3          ma2;
    const bool          mUseLimits;
    Vec3                ms1;
    Vec3                ms2;
    const double        mPhi;
    const bool          mUseFriction;
    const double        mLambdaMax;
    const double        mkCoeff;

    Vec3                ma1_gcs;
    Vec3                ma2_gcs;
    Vec3                mr1_gcs;
    Vec3                mr2_gcs;
    Vec3                ms1_gcs;
    Vec3                ms2_gcs;
    Mat3x3              mRa1;
    Mat3x3              mRa2;

    JacobianConstraint* mcLink1;
    JacobianConstraint* mcLink2;
    JacobianConstraint* mcLink3;
    JacobianConstraint* mcLimit1;
    JacobianConstraint* mcFric1;
    JacobianConstraint* mcFric2;
    JacobianConstraint* mcFric3;

    friend class JointManager;
};


}// namespace Makena

#endif /* _MAKENA_JOINT_MANAGER_HPP_*/

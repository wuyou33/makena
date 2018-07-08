#ifndef _MAKENA_CONVEX_RIGID_BODY_HPP_
#define _MAKENA_CONVEX_RIGID_BODY_HPP_

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <exception>

#include "primitives.hpp"
#include "manifold.hpp"
#include "aabb.hpp"
#include "orienting_bounding_box.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file convex_rigid_body.hpp
 *
 * @brief It represents aconvex ridig body for physics simulation.
 *        It has the shape in LCSrepresented by Manifold, a geometric
 *        configuration (position, orientation,linear velocity, and angular 
 *        velocity), forces, torques, mass, inertia matrix, and friction 
 *        coefficients.
 *        It alsohas AABB and OBB internally up-to-date,depending on the 
 *        current position and orientation in GCS.
 */
namespace Makena {

/** @brief limitter to the maximum RPM allowed for rigid body (5RPM) */
static constexpr double MAX_ANGULAR_SPEED  = 1000000.0*2.0*M_PI;

class JacobianConstraint;

/** @class ConvexRigidBody
 *
 *  @brief the top level representation of a convex rigid body for simulation.
 *         It contains the following essential information.
 *
 *         - Shape (local coordinates) in a compact convex manifold.
 *           Planar graph and its dual of vertex, edge, and face.
 *           with local 3D geometric coordinates , 2D UV texture coordinates,
 *           and normals.
 *
 *         - Mass & inertia
 *
 *         - Geometric configuration at time step t.
 *               CoM, orientation, linear and angular velocities.
 *
 *         - External force and torque at time step t + Δt
 * 
 *         It also holds following auxiliary information for collision
 *         management and for the constraint engine.
 *
 *         - Tentative geometric configuration.
 *
 *         - Axis-aligned bounding box
 *
 *         - Oriented bounding box
 *
 *         - State information for constraint engine
 *
 *         - Some values for constraint engine, such as inverse of inertia,
 *           angular momentum ω x Iω, induced force and torque.
 */
class ConvexRigidBody {

  public:

    enum engineState {

        UNREGISTERED,  /** @brief not registered to ConstraintEngine  */

        FREE,          /** @brief not participating in any constraint */

        RESET,         /** @brief reset for the current simulation step
                        *         mUlin and mUang not yet set.
                        */

        MANAGED,       /** @brief internal values for the current 
                        *         simulation step for the constraints 
                        *         that use this body are all set.
                        */

        NEEDS_UPDATE,  /** @brief internal values for the current 
                        *         simulation step need to be updated
                        *         as some new constraints that use this
                        *         body have been added.
                        */
        END,
    };

    /** @brief constructor
     *
     *  @param id (in): Unique ID for this object.
     */
    inline ConvexRigidBody(const long id);


    /** @brief descructor */
    inline ~ConvexRigidBody();


    /** @brief sets the attributes of this object that are supposed to be
     *         invariant in the real world.
     *
     *  @param CH        (in): Convex hull (martialized)
     *
     *  @param mass      (in): Mass of this object
     *
     *  @param I         (in): Inertial matrix in LCS
     *
     *  @param MuDynamic (in): Dynamic friction coefficient
     *
     *  @param MuStatic  (in): Static friction coefficient
     *
     *  @param fixed     (in): True if this object is stationary and 
     *                         fixed to the simulation world
     */
    inline void setObjectAttributes(
        Manifold::Martialled& CH,
        const double&         mass,
        const Mat3x3&         I,
        const double&         MuDynamic,
        const double&         MuStatic,
        bool                  fixed
    );


    /** @brief sets the geometric configuration of this object.
     *         Called before the simulation starts for time 0.
     *
     *  @param CoM       (in): Center of Mass in GCS .
     *
     *  @param Qori      (in): Orientation of this object in GCS 
     *
     *  @param Vlin      (in): Linear velocity at time 0. Must be zero if 
     *                         this object is fixed.
     *
     *  @param Vang      (in): Angular velocity at time 0. Must be zero if 
     *                         this object is fixed.
     */
    inline void setGeomConfig(
        const Vec3&       CoM,
        const Quaternion& q,
        const Vec3&       Vlin,
        const Vec3&       Vang
    );


    /** @brief perform Euler integration and update the geometry config
     *         using the external force & torque as well as induced ones.
     *         The updated configurations are stored to:
     *             mVlin_uncommitted
     *             mVang_uncommitted
     *             mCoM_uncommitted
     *             mQ_uncommitted
     *         It also updates the information required for collision detection
     *         based on the updated values.
     */      
    inline void updateGeomConfig(const double& deltaT);


    /** @brief commits the geometric configuration stored in the temporary
     *         variables to the official ones.
     *         It also updates the inertia matrix based on the committed 
     *         orientation.
     */      
    inline void commitGeomConfig();


    inline long              id() const;
    inline bool              isFixed() const;
    inline const Mat3x3&     Iinv() const;
    inline double            mInv() const;

    inline void              resetExternalForceTorque();
    inline void              addExternalForce(const Vec3& F);
    inline void              addExternalTorque(const Vec3& T);

    inline Manifold&         ConvexHull();
    inline const Quaternion& Q() const;
    inline const Mat3x3      Qmat() const;
    inline const Vec3&       CoM()  const;
    inline const Vec3&       Vlin() const;
    inline const Vec3&       Vang() const;
    inline const Vec3&       Fext() const;
    inline const Vec3&       Text() const;
    inline const Vec3&       FintDt() const;
    inline const Vec3&       TintDt() const;

    inline double            muStatic()  const;
    inline double            muDynamic() const;


    /** @brief following functions are used for collision detection */
    inline AABBelem&         aabb();
    inline const Vec3&       OBBcenter()  const;
    inline const Mat3x3&     OBBaxes()    const;
    inline const Vec3&       OBBextents() const;
    inline const Mat3x3      QmatTemp()  const;
    inline const Vec3&       CoMTemp()   const;
    inline const Vec3&       VlinTemp()  const;
    inline const Vec3&       VangTemp()  const;


    /** @brief following function are used by ContraintEngine */
    
    inline enum engineState  engineState() const;
    inline void              setEngineState(enum engineState s);

    inline list<ConvexRigidBody*>::iterator
                             backItConstraintManager() const;

    inline void              setBackItConstraintManager(
                                         list<ConvexRigidBody*>::iterator it);


    inline void              addConstraintNotCheckedIn(JacobianConstraint* c);

    inline void              moveToCheckedIn(JacobianConstraint* c);

    inline void              moveToNotCheckedIn(JacobianConstraint* c);

    inline void              removeConstraint(JacobianConstraint* c);

    inline list<JacobianConstraint*>& 
                             constraintsNotCheckedIn();

    inline list<JacobianConstraint*>& 
                             constraintsCheckedIn();

    inline long              numConstraintsNotCheckedIn() const;

    inline long              numConstraintsCheckedIn()const;

    inline void              setU(const double& deltaT);

    inline const Vec3&       Ulin() const;
    inline const Vec3&       Uang() const;

    inline void              resetInducedForceTorque();
    inline void              addInducedForceTorque(
                       const Vec3& lin, const Vec3& ang, const double& lambda);

  private:

    /** @brief unique ID for this object. For CollisionContact pairs in 
     *         pair<long,long>, the first ID is always smaller than the second
     *         to avoid ambiguity.
     */
    const long mId;

    /** @brief true if this object is stationary and fixed in GCS */
    bool       mFixed;

    /** @brief reciprocal of mass */
    double     mMinv;

    /** @brief Inertia matrix in LCS*/
    Mat3x3     mI_lcs;

    /** @brief static friction coefficient for Coulomb model friction */
    double     mMuStatic;

    /** @brief dynamic friction coefficient for Coulomb model friction */
    double     mMuDynamic;

    /** @brief Axis-Aligned Bounding Box.  Used by the broad-phase
     *         collision culling.
     */
    AABBelem   mAABB;

    /** @brief The center of the OBB in LCS*/
    Vec3       mOBBcenter_lcs;

    /** @brief The 3 orthonormal vectors that represent the direction of
     *         the extents of the OBB.
     */
    Mat3x3     mOBBAxes_lcs;

    /** @brief The center of the OBB in GCS*/
    Vec3       mOBBcenter;

    /** @brief The 3 orthonormal extent vectors in GCS*/
    Mat3x3     mOBBAxes;

    /** @brief The extents of the OBB. The lengths of the edges of the OBB.*/
    Vec3       mOBBExtents;

    /** @brief The convex hull that defines the shape of this object */
    Manifold   mCH;

    /** @brief internally induced momentum in GCS due to rotation */
    Vec3       mwxIw;

    /** @brief Center of Mass in GCS */
    Vec3       mCoM;
    Vec3       mCoM_uncommitted;

    /** @brief orientation of this object in GCS in quaternion */
    Quaternion mQ;
    Quaternion mQ_uncommitted;

    /** @brief orientation of this object in GCS in rotation matrix */
    Mat3x3     mQMat;
    Mat3x3     mQMat_uncommitted;

    /** @brief linear velocity of this object */
    Vec3       mVlin;
    Vec3       mVlin_uncommitted;

    /** @brief angular velocity of this object */
    Vec3       mVang;
    Vec3       mVang_uncommitted;

    /** @brief external force applied to this object such as gravity */
    Vec3       mFext;

    /** @brief external force applied to this object */
    Vec3       mText;

    /** @brief Inertia in GCS */
    Mat3x3     mI;

    /** @brief Inverse of the inertia in GCS */
    Mat3x3     mIinv;

    /** @brief linear impulse induced by ConstraintManager */
    Vec3       mFintDt;

    /** @brief angular impulse induced by ConstraintManager */
    Vec3       mTintDt;

    /** @brief used by ConstraintManager. It holds the estimated linear
     *         velocity at the next step without induced force and torque.
     *         This is used to generate q vector for MLCP.
     */
    Vec3       mUlin;

    /** @brief used by ConstraintManager. It holds the estimated angular
     *         velocity at the next step without induced force and torque.
     *         This is used to generate q vector for MLCP.
     */
    Vec3       mUang;

    /** @brief the state of this object for ConstraintManager */
    enum engineState 
               mEngineState;

    /** @brief constraints that involve this object */
    list<JacobianConstraint*> mConstraintsNotCheckedIn;
    list<JacobianConstraint*> mConstraintsCheckedIn;

    /** @brief back iterator to a list in Engine */
    list<ConvexRigidBody*>::iterator mBackItConstraintManager;

};

} // namespace Makena

#include "jacobian_constraint.hpp"

namespace Makena {


ConvexRigidBody::ConvexRigidBody(const long id):
    mId(id),
    mAABB(this),
    mEngineState(UNREGISTERED)
{
    mAABB.setId(id);
}


ConvexRigidBody:: ~ConvexRigidBody() {;}


void ConvexRigidBody::setObjectAttributes(
    Manifold::Martialled& CH,
    const double&         mass,
    const Mat3x3&         I,
    const double&         MuDynamic,
    const double&         MuStatic,
    bool                  fixed
) {
    mCH.importData(CH);
    // Make OBB in LCS.
    Manifold obb;
    double   obbVolume;
    findOBB3D(mCH, obb, mOBBAxes_lcs, mOBBcenter_lcs, mOBBExtents, obbVolume);
    mFixed     = fixed;
    mMinv      = 1.0/mass;
    mI_lcs     = I;
    mMuStatic  = MuStatic;
    mMuDynamic = MuDynamic;
    mAABB.setFixed(fixed);
}


void ConvexRigidBody::setGeomConfig(
    const Vec3&       CoM,
    const Quaternion& q,
    const Vec3&       Vlin,
    const Vec3&       Vang
) {
    mCoM  = CoM;
    mQ    = q;
    mQ.normalize();
    mQMat = q.rotationMatrix();

    if (!mFixed) {
        mVlin = Vlin;
        mVang = Vang;
    }
    else {
        mVlin.zero();
        mVang.zero();
    }

    mFext.zero();
    mText.zero();
    mFintDt.zero();
    mTintDt.zero();
  
    Mat3x3 Rt = mQMat.transpose();
    mI    = mQMat * mI_lcs * Rt;
    mIinv = mI.inverse();
    mwxIw = mVang.cross(mI * mVang);

    // Set OBB in GCS
    mOBBAxes    = mQMat * mOBBAxes_lcs;
    mOBBcenter  = mCoM + mQ.rotationMatrix() * mOBBcenter_lcs;

    // Set AABB
    mAABB.updateFromOBB(mOBBcenter, mOBBAxes, mOBBExtents);

    mVlin_uncommitted = mVlin;
    mVang_uncommitted = mVang;
    mCoM_uncommitted  = mCoM;
    mQ_uncommitted    = mQ;
    mQMat_uncommitted = mQ_uncommitted.rotationMatrix();
}


void ConvexRigidBody::updateGeomConfig(const double& deltaT)
{
    if (mFixed) {
        return;
    }

    mVlin_uncommitted = mVlin + (mFext * deltaT  + mFintDt) * mMinv;
    mVang_uncommitted = mVang + mIinv * ((mText - mwxIw) * deltaT  + mTintDt);

    // Limit the angular velocity
    double revolutionSq = mVang_uncommitted.squaredNorm2();
    if (revolutionSq>(MAX_ANGULAR_SPEED*MAX_ANGULAR_SPEED)) {
        double revolution = sqrt(revolutionSq);
        mVang_uncommitted.scale(MAX_ANGULAR_SPEED/revolution);
    }

    mCoM_uncommitted  = mCoM + mVlin_uncommitted * deltaT;

    double Q4x3[12];
    mQ.matrix4x3(Q4x3);
    const double & wx = mVang_uncommitted.x();
    const double & wy = mVang_uncommitted.y();
    const double & wz = mVang_uncommitted.z();

    Quaternion Qdelta(
            Q4x3[ 0]*wx + Q4x3[ 1]*wy + Q4x3[ 2]*wz,
            Q4x3[ 3]*wx + Q4x3[ 4]*wy + Q4x3[ 5]*wz,
            Q4x3[ 6]*wx + Q4x3[ 7]*wy + Q4x3[ 8]*wz,
            Q4x3[ 9]*wx + Q4x3[10]*wy + Q4x3[11]*wz
    );

    Qdelta.scale(deltaT);
    mQ_uncommitted = mQ + Qdelta;
    mQ_uncommitted.normalize();
    mQMat_uncommitted = mQ_uncommitted.rotationMatrix();

    // Update OBB
    mOBBAxes    = mQMat_uncommitted * mOBBAxes_lcs;
    mOBBcenter  = mCoM_uncommitted + mQ_uncommitted.rotate(mOBBcenter_lcs);

    // Update AABB
    mAABB.updateFromOBB(mOBBcenter, mOBBAxes, mOBBExtents);


}


void ConvexRigidBody::commitGeomConfig()
{
    if (mFixed) {
        return;
    }
    mVlin = mVlin_uncommitted;
    mVang = mVang_uncommitted;
    mCoM  = mCoM_uncommitted;
    mQ    = mQ_uncommitted;
    mQMat = mQMat_uncommitted;
    Mat3x3 Rt = mQMat.transpose();
    mI    = mQMat * mI_lcs * Rt;
    mIinv = mI.inverse();
    mwxIw = mVang.cross(mI * mVang);
}


long          ConvexRigidBody::id() const { return mId; }
bool          ConvexRigidBody::isFixed() const {return mFixed; }
const Mat3x3& ConvexRigidBody::Iinv() const { return mIinv; }
double        ConvexRigidBody::mInv() const { return mMinv;  }
void          ConvexRigidBody::resetExternalForceTorque()
{
    mFext.zero(); mText.zero();
}
void          ConvexRigidBody::addExternalForce(const Vec3& F){ mFext += F; }
void          ConvexRigidBody::addExternalTorque(const Vec3& T){ mText += T; }
Manifold&     ConvexRigidBody::ConvexHull() { return mCH; }
const Quaternion&
              ConvexRigidBody::Q() const { return mQ; }
const Mat3x3  ConvexRigidBody::Qmat() const { return mQMat; }
const Vec3&   ConvexRigidBody::CoM()  const { return mCoM;  }
const Vec3&   ConvexRigidBody::Vlin() const { return mVlin; }
const Vec3&   ConvexRigidBody::Vang() const { return mVang; }
const Vec3&   ConvexRigidBody::Fext() const { return mFext; }
const Vec3&   ConvexRigidBody::Text() const { return mText; }
const Vec3&   ConvexRigidBody::FintDt() const { return mFintDt; }
const Vec3&   ConvexRigidBody::TintDt() const { return mTintDt; }
double        ConvexRigidBody::muStatic()  const { return mMuStatic;  }
double        ConvexRigidBody::muDynamic() const { return mMuDynamic; }
AABBelem&     ConvexRigidBody::aabb()             { return mAABB; };
const Vec3&   ConvexRigidBody::OBBcenter()  const { return mOBBcenter; }
const Mat3x3& ConvexRigidBody::OBBaxes()    const { return mOBBAxes; }
const Vec3&   ConvexRigidBody::OBBextents() const { return mOBBExtents; }
const Mat3x3  ConvexRigidBody::QmatTemp()  const { return mQMat_uncommitted; }
const Vec3&   ConvexRigidBody::CoMTemp()   const { return mCoM_uncommitted;  }
const Vec3&   ConvexRigidBody::VlinTemp()  const { return mVlin_uncommitted; }
const Vec3&   ConvexRigidBody::VangTemp()  const { return mVang_uncommitted; }


enum ConvexRigidBody::engineState ConvexRigidBody::engineState() const
{
    return mEngineState;
}


void ConvexRigidBody::setEngineState(enum engineState s) { mEngineState = s; }


list<ConvexRigidBody*>::iterator ConvexRigidBody::backItConstraintManager()
const
{ 
    return mBackItConstraintManager;
}


void ConvexRigidBody::setBackItConstraintManager(
    list<ConvexRigidBody*>::iterator it
) {
    mBackItConstraintManager = it;
}


void ConvexRigidBody::addConstraintNotCheckedIn(JacobianConstraint* c)
{
    auto it=mConstraintsNotCheckedIn.insert(mConstraintsNotCheckedIn.end(),c);
                                          
    if (c->body1()==this) {

        c->setBackItBody1(it);
    }
    else {

        c->setBackItBody2(it);
    }
}


void ConvexRigidBody::moveToCheckedIn(JacobianConstraint* c)
{
    if (c->body1()==this) {
    
        mConstraintsNotCheckedIn.erase(c->backItBody1());

        auto it = mConstraintsCheckedIn.insert(mConstraintsCheckedIn.end(), c);

        c->setBackItBody1(it);                                       

    }
    else {

        mConstraintsNotCheckedIn.erase(c->backItBody2());

        auto it = mConstraintsCheckedIn.insert(mConstraintsCheckedIn.end(), c);

        c->setBackItBody2(it);
    }
}                                                


void ConvexRigidBody::moveToNotCheckedIn(JacobianConstraint* c)
{
    if (c->body1()==this) {

        mConstraintsCheckedIn.erase(c->backItBody1());

        auto it = mConstraintsNotCheckedIn.insert(
                                            mConstraintsNotCheckedIn.end(), c);

        c->setBackItBody1(it);
    }
    else {
        mConstraintsCheckedIn.erase(c->backItBody2());

        auto it = mConstraintsNotCheckedIn.insert(
                                            mConstraintsNotCheckedIn.end(), c);

        c->setBackItBody2(it);
    }
}


void ConvexRigidBody::removeConstraint(JacobianConstraint* c)
{
    if (c->engineState()==JacobianConstraint::NOT_CHECKED_IN) {
        mConstraintsNotCheckedIn.erase(
                     (c->body1()==this)?(c->backItBody1()):(c->backItBody2()));
    }
    else {
        mConstraintsCheckedIn.erase(
                     (c->body1()==this)?(c->backItBody1()):(c->backItBody2()));
    }
}


list<JacobianConstraint*>& ConvexRigidBody::constraintsNotCheckedIn() {
    return mConstraintsNotCheckedIn;
}


list<JacobianConstraint*>& ConvexRigidBody::constraintsCheckedIn() {
    return mConstraintsCheckedIn;
}


long ConvexRigidBody::numConstraintsNotCheckedIn() const {
    return mConstraintsNotCheckedIn.size();
}


long ConvexRigidBody::numConstraintsCheckedIn()const{
    return mConstraintsCheckedIn.size();
}

inline void   ConvexRigidBody::setU(const double& deltaT) {
    mUlin = mVlin + ((mFext * deltaT ) * mMinv);
    mUang = mVang + (mIinv * ((mText - mwxIw) * deltaT));
}


const Vec3&   ConvexRigidBody::Ulin() const { return mUlin; }


const Vec3&   ConvexRigidBody::Uang() const { return mUang; }


void   ConvexRigidBody::resetInducedForceTorque()
{
    mFintDt.zero(); mTintDt.zero();
}


void   ConvexRigidBody::addInducedForceTorque(
                       const Vec3& lin, const Vec3& ang, const double& lambda)
{
    mFintDt += (lin * lambda);
    mTintDt += (ang * lambda);
}


}// namespace Makena

#endif /*_MAKENA_CONVEX_RIGID_BODY_HPP_*/

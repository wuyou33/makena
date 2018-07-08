#ifndef _MAKENA_AABB_HPP_
#define _MAKENA_AABB_HPP_

#include <exception>
#include <stdexcept>
#include <set>
#include <map>
#include "primitives.hpp"


#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file aapp.hpp
 *
 * @brief Represents an axis-alignedbounding box.
 */
namespace Makena {

using namespace std;

/* forward declaration */
class AABBelem;
class ConvexRigidBody;


/** @class AABBScanElem
 *  @brief an element for the scan array.
 */
class AABBScanElem {
  public:

    inline AABBScanElem(
        AABBelem&       ref, 
        double&         valRef,
        const bool      isEnter,
        bool&           disabled
    );

    inline ~AABBScanElem();

    inline bool operator< (const AABBScanElem& rhs) const;

    inline bool operator<=(const AABBScanElem& rhs) const;

    inline bool operator> (const AABBScanElem& rhs) const;

    inline bool operator>=(const AABBScanElem& rhs) const;

    AABBelem&                     mRef;
    double&                       mValRef;
    const bool                    mIsEnter;
    bool&                         mDisabled;
    list<AABBScanElem*>::iterator mScanLineBackIt;

};


/** @class AABBelem
 *
 *  @brief Represents one solid object for AABB collision detection.
 *         RedBlackTreeNode used for sweeps.
 */

class AABBelem {

  public:

    inline AABBelem(ConvexRigidBody* const body);

    inline void setCoords(const Vec3& low, const Vec3& high);

    inline void setFixed(bool fixed);
    inline bool isFixed()const;
    inline void disable();
    inline void enable();

    inline void updateFromOBB(
        const Vec3&   com_gcs,
        const Mat3x3& axes_gcs,
        const Vec3&   extents
    );

    inline std::set<AABBelem*>& exclusionList();

    inline ConvexRigidBody* rigidBody();

    inline void setId(long id);

    inline long id() const;

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    /** @brief pointer to the solid object. */
    ConvexRigidBody*               mRigidBody;

    long                           mId;

    /** @brief back it to the list in BroadPhaseAABBCollisionDetector */
    list<AABBelem*>::iterator      mBackIt;

    /** @brief set to true if this object has to be excluded from 
     *         collision detection.
     */
    bool                           mDisabled;

    /** @brief true if this object is a fixed object.
     */
    bool                           mFixed;

    /** @brief holds the IDs of other objects that can't be collided
     *         by design (e.g. the incident arm over a joint.)
     */
    std::set<AABBelem*>            mExclusionList;

    /** @brief coordinates of the lower side of AABB. */
    double                         mLowX;
    double                         mLowY;
    double                         mLowZ;

    /** @brief coordinates of the higher side of AABB. */
    double                         mHighX;
    double                         mHighY;
    double                         mHighZ;

    AABBScanElem                   mScanElemXE;
    AABBScanElem                   mScanElemXL;
    AABBScanElem                   mScanElemYE;
    AABBScanElem                   mScanElemYL;
    AABBScanElem                   mScanElemZE;
    AABBScanElem                   mScanElemZL;


friend class BroadPhaseAABBCollisionDetector;

};


AABBScanElem::AABBScanElem(
    AABBelem&  ref, 
    double&    valRef,
    const bool isEnter,
    bool&      disabled
):
    mRef(ref),
    mValRef(valRef),
    mIsEnter(isEnter),
    mDisabled(disabled){;}


AABBScanElem::~AABBScanElem(){;}


bool AABBScanElem::operator<(const AABBScanElem& rhs) const {
    return (mValRef < rhs.mValRef) ||
           (mValRef == rhs.mValRef && mIsEnter && !(rhs.mIsEnter));
}


bool AABBScanElem::operator<=(const AABBScanElem& rhs) const {
    return (mValRef < rhs.mValRef) ||
           (mValRef == rhs.mValRef && mIsEnter && !(rhs.mIsEnter))||
           (mValRef == rhs.mValRef && mIsEnter == rhs.mIsEnter);
}


bool AABBScanElem::operator>(const AABBScanElem& rhs) const {
    return (mValRef > rhs.mValRef) ||
           (mValRef == rhs.mValRef && !mIsEnter && rhs.mIsEnter);
}


bool AABBScanElem::operator>=(const AABBScanElem& rhs) const {
    return (mValRef > rhs.mValRef) ||
           (mValRef == rhs.mValRef && !mIsEnter && rhs.mIsEnter)||
           (mValRef == rhs.mValRef && mIsEnter == rhs.mIsEnter);
}


AABBelem::AABBelem(ConvexRigidBody* const body):
    mRigidBody (body),
    mId        (-1),
    mDisabled  (false),
    mFixed     (false),
    mLowX      (0.0),
    mLowY      (0.0),
    mLowZ      (0.0),
    mHighX     (0.0),
    mHighY     (0.0),
    mHighZ     (0.0),
    mScanElemXE( *this, mLowX,  true,  mDisabled ),
    mScanElemXL( *this, mHighX, false, mDisabled ),
    mScanElemYE( *this, mLowY,  true,  mDisabled ),
    mScanElemYL( *this, mHighY, false, mDisabled ),
    mScanElemZE( *this, mLowZ,  true,  mDisabled ),
    mScanElemZL( *this, mHighZ, false, mDisabled )  {;}


void AABBelem::setId(long id) { mId = id;   }

long AABBelem::id() const     { return mId; }

void AABBelem::setCoords(const Vec3& low, const Vec3& high)
{

    mLowX  = low.x();
    mLowY  = low.y();
    mLowZ  = low.z();
    mHighX = high.x();
    mHighY = high.y();
    mHighZ = high.z();

}

void AABBelem::setFixed(bool fixed) { mFixed = fixed;  }

bool AABBelem::isFixed() const { return mFixed; }

void AABBelem::disable() { mDisabled = true;  }

void AABBelem::enable()  { mDisabled = false; }

void AABBelem::updateFromOBB(
              const Vec3& com_gcs, const Mat3x3& axes_gcs, const Vec3& extents)
{
    double extX = fabs(axes_gcs.val(1,1))*extents.x() +
                  fabs(axes_gcs.val(1,2))*extents.y() +
                  fabs(axes_gcs.val(1,3))*extents.z() ;

    double extY = fabs(axes_gcs.val(2,1))*extents.x() +
                  fabs(axes_gcs.val(2,2))*extents.y() +
                  fabs(axes_gcs.val(2,3))*extents.z() ;

    double extZ = fabs(axes_gcs.val(3,1))*extents.x() +
                  fabs(axes_gcs.val(3,2))*extents.y() +
                  fabs(axes_gcs.val(3,3))*extents.z() ;

    mLowX  = com_gcs.x() - extX;
    mHighX = com_gcs.x() + extX;
    mLowY  = com_gcs.y() - extY;
    mHighY = com_gcs.y() + extY;
    mLowZ  = com_gcs.z() - extZ;
    mHighZ = com_gcs.z() + extZ;

}


std::set<AABBelem*>& AABBelem::exclusionList() { return mExclusionList; }


ConvexRigidBody* AABBelem::rigidBody() { return mRigidBody; }


} // namespace Makena

#endif /*_MAKENA_AABB_HPP_*/

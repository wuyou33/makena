#ifndef _MAKENA_BROAD_PHASE_AABB_COLLISION_DETECTOR_HPP_
#define _MAKENA_BROAD_PHASE_AABB_COLLISION_DETECTOR_HPP_

#include <exception>
#include <stdexcept>
#include <set>
#include <map>
#include "primitives.hpp"
#include "aabb.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file broad_phase_aabb_collision_detector.hpp
 *
 * @brief Provides broad-phase collision culling using AABBs.
 *        It usesenter-leave sweeping paradigm on X, Y, and Z axes.
 *        For each axis, the order or theobjectsis updated using bubble sort.
 */
namespace Makena {

using namespace std;

/** @class BroadPhaseAABBCollisionDetector
 *
 *  @brief 
 *
 *  @remark  complexity on speed is O(N), and  O(N^2) for memory.
 */
class BroadPhaseAABBCollisionDetector {

  public: 

    BroadPhaseAABBCollisionDetector()  {;}
    ~BroadPhaseAABBCollisionDetector() {;}

    /** @brief register a solid object to be tested
     *
     *  @param e  (in): the object to be tested.
     *
     *  @remark a registered object has mIndex set internally.
     */
    inline void registerElem(AABBelem* e);


    /** @brief unregister a solid object permanently.
     *
     *  @param e  (in): the object to be unregistered.
     *
     *  @remark This is an O(N) operation and should not be called frequently.
     *          To exclude this object from testing temporarily use
     *          AABBelem::disable()/enable() instead.
     */
    inline void unregisterElem(AABBelem* e);


    /** @brief updates an internal structure and generate colliding pairs.
     *         This is expected to be called in the rendering loop per frame.
     *
     *  @remark This is asymptotically O(N + C) operation where C is the 
     *          number of the colliding pairs.
     */
    inline void update();


    /** @brief returs a reference to the current colliding pairs.
     */
    inline std::list<pair<ConvexRigidBody*const, ConvexRigidBody*const> >& 
                                                   pairs() { return mPairs; }

    inline list<AABBelem*>& elems() { return mElems; }

  private:    

    list<AABBelem*>     mElems;
    list<AABBScanElem*> mScanLineX;
    list<AABBScanElem*> mScanLineY;
    list<AABBScanElem*> mScanLineZ;

    class TableElem {
      public:
        TableElem(std::list<pair<ConvexRigidBody*const, 
                                 ConvexRigidBody*const> >::iterator it):
                                                         mCnt(1),mPairIt(it){;}
        ~TableElem(){;}

        unsigned char                                     mCnt;
        std::list<pair<ConvexRigidBody*const,
                       ConvexRigidBody*const> >::iterator mPairIt;
    };


    std::map<std::pair<AABBelem*, AABBelem*>, TableElem>           mTable;
    std::list<pair<ConvexRigidBody*const, ConvexRigidBody*const> > mPairs;

    inline void countUpTable(AABBelem* i, AABBelem* j);

    inline void countDownTable(AABBelem* i, AABBelem* j);

    inline void scanSortAndUpdate(list<AABBScanElem*>& line) {
        for(auto it = line.begin(); it != line.end(); it++) {
            swapToTheLeft(line, it);
        }
    }

    inline void swapToTheLeft(
        list<AABBScanElem*>&          line,
        list<AABBScanElem*>::iterator it
    ) {
        while (it != line.begin()) {
            auto prevIt = it;
            prevIt--;
            if (*(*prevIt) > *(*it)) {
                auto* eR = *it;
                auto* eL = *prevIt;
                if (eR->mIsEnter && !(eL->mIsEnter)) {
                    countUpTable(   &(eL->mRef), &(eR->mRef) );
                }
                else if (!(eR->mIsEnter) && eL->mIsEnter) {
                    countDownTable( &(eL->mRef), &(eR->mRef) );
                }
                std::swap(*it, *prevIt);
                it--;
            }
            else {
                break;
            }
        }
    }


#ifdef UNIT_TESTS
  public:
    inline std::vector<pair<ConvexRigidBody*const, ConvexRigidBody*const> >
                findCollidingPairsNbyN();
    inline bool AABBcrudeTestPair(AABBelem* e1, AABBelem* e2);
#endif

};


inline void BroadPhaseAABBCollisionDetector::registerElem(AABBelem* e)
{
    e->mBackIt = mElems.insert(mElems.end(), e);

    e->mScanElemXE.mScanLineBackIt =
                        mScanLineX.insert(mScanLineX.end(), &(e->mScanElemXE));
    e->mScanElemXL.mScanLineBackIt =
                        mScanLineX.insert(mScanLineX.end(), &(e->mScanElemXL));
    e->mScanElemYE.mScanLineBackIt =
                         mScanLineY.insert(mScanLineY.end(),&(e->mScanElemYE));
    e->mScanElemYL.mScanLineBackIt =
                         mScanLineY.insert(mScanLineY.end(),&(e->mScanElemYL));
    e->mScanElemZE.mScanLineBackIt =
                         mScanLineZ.insert(mScanLineZ.end(),&(e->mScanElemZE));
    e->mScanElemZL.mScanLineBackIt =
                         mScanLineZ.insert(mScanLineZ.end(),&(e->mScanElemZL));
}


void BroadPhaseAABBCollisionDetector::unregisterElem(AABBelem* e)
{
    mElems.erase(e->mBackIt);
    mScanLineX.erase(e->mScanElemXE.mScanLineBackIt);
    mScanLineX.erase(e->mScanElemXL.mScanLineBackIt);
    mScanLineY.erase(e->mScanElemYE.mScanLineBackIt);
    mScanLineY.erase(e->mScanElemYL.mScanLineBackIt);
    mScanLineZ.erase(e->mScanElemZE.mScanLineBackIt);
    mScanLineZ.erase(e->mScanElemZL.mScanLineBackIt);
}


inline void BroadPhaseAABBCollisionDetector::update()
{
    scanSortAndUpdate(mScanLineX);
    scanSortAndUpdate(mScanLineY);
    scanSortAndUpdate(mScanLineZ);
}


inline void BroadPhaseAABBCollisionDetector::countUpTable(
    AABBelem* i,
    AABBelem* j
) {
    if (i->id() > j->id()){
        std::swap(i,j);
    }

    pair<AABBelem*, AABBelem*> key(i,j);
    auto it = mTable.find(key);
    if (it == mTable.end()) {
        TableElem e(mPairs.end());
        mTable.insert(make_pair(key,e));
    }
    else {
        auto& e = it->second;
        if (++(e.mCnt) == 3) {
            auto& ki = key.first;
            auto& kj = key.second;
            if (!(ki->mFixed) || !(kj->mFixed)) {
                if (!(ki->mDisabled) && !(kj->mDisabled)) {
                    if (ki->mExclusionList.find(kj)==ki->mExclusionList.end()){
                        e.mPairIt = mPairs.insert( mPairs.end(), 
                                  make_pair(ki->rigidBody(), kj->rigidBody()));
                    }
                }
            }
        }
    }
}


inline void BroadPhaseAABBCollisionDetector::countDownTable(
    AABBelem* i,
    AABBelem* j
) {
    if ((i->mFixed) && (j->mFixed)) {
        return;
    }

    if (i->id() > j->id()){
        std::swap(i,j);
    }

    pair<AABBelem*, AABBelem*> key(i,j);
    auto it = mTable.find(key);
    if (it == mTable.end()) {
        throw std::logic_error("countDowntable can't find element");
    }
    else {
        auto& e = it->second;
        --(e.mCnt);
        if (e.mCnt == 2) {
            if (e.mPairIt != mPairs.end()) {
                mPairs.erase(e.mPairIt);
                e.mPairIt = mPairs.end();
            }
        }
        else if (e.mCnt == 0) {
            mTable.erase(it);
        }
    }
}


#ifdef UNIT_TESTS

inline std::vector<pair<ConvexRigidBody*const, ConvexRigidBody*const> >
                      BroadPhaseAABBCollisionDetector::findCollidingPairsNbyN()
{
    std::vector<pair<ConvexRigidBody*const, ConvexRigidBody*const> > pairs;
    for (auto it1 = mElems.begin(); it1 != mElems.end(); it1++)  {
        if (!(*it1)->mDisabled) {
            auto it2 = it1;
            for(it2++; it2 != mElems.end(); it2++) {
                if (!(*it2)->mDisabled) {
                    if ( !((*it1)->mFixed) || !((*it2)->mFixed)) {
                        if (AABBcrudeTestPair(*it1, *it2)) {
                            if ((*it1)->mExclusionList.find(*it2)==
                                                (*it1)->mExclusionList.end()) {
                                if ((*it1)->id() < (*it2)->id()) {
                                
                                    pairs.push_back(make_pair(
                                        (*it1)->rigidBody(),
                                        (*it2)->rigidBody()
                                    ));
                                }
                                else {
                                    pairs.push_back(make_pair(
                                        (*it2)->rigidBody(),
                                        (*it1)->rigidBody()
                                    ));
                                }                                        
                            }
                        }
                    }
                }
            }
        }
    }
    return pairs;
}


inline bool BroadPhaseAABBCollisionDetector::AABBcrudeTestPair(
    AABBelem* e1,
    AABBelem* e2
) {
    if ((e2->mHighX < e1->mLowX)||(e1->mHighX < e2->mLowX)||
        (e2->mHighY < e1->mLowY)||(e1->mHighY < e2->mLowY)||
        (e2->mHighZ < e1->mLowZ)||(e1->mHighZ < e2->mLowZ)  ) {
        return false;
    }
    return true;
}

#endif


} // namespace Makena

#endif /*_MAKENA_BROAD_PHASE_AABB_COLLISION_DETECTOR_HPP_*/

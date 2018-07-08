#ifndef _MAKENA_CONSTRAINT_MANAGER_HPP_
#define _MAKENA_CONSTRAINT_MANAGER_HPP_

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <exception>

#include "primitives.hpp"
#include "jacobian_constraint.hpp"
#include "mlcp.hpp"
#include "broad_phase_aabb_collision_detector.hpp"
#include "convex_rigid_body.hpp"
#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file constraint_engine.hpp
 *
 * @brief Core ofthe engine. It takes the constraints, and calculates the 
 *        lambdas (lagrange multipliers) by MLCP,which are used to calculate 
 *        theinduced internal forcesand torques at each simulation step. 
 *        Then it updates the geometric configurationof eachconvex rigid body.
 */
namespace Makena {


/** @class ConstraintManager
 *
 *  @brief Generates internal forces and torques for each convex rigid body 
 *         based on the given constraints.
 *        
 *         It interacts with the following three types of objects (classes).
 *
 *         ConvexRigidBody :    This represents a rigid body in terms of the 
 *                              geometric configuration.
 *
 *         JacobianConstraint : This represents one bilateral or unilateral
 *                              constraint between two ConvexRigidBodies in 
 *                              the velocity space. Usually it takes multiple 
 *                              constraints to implement a joint or contact.
 *
 *         MLCP :               This is an implementation of Mixed Linear
 *                              Complementarity Problem solver.
 *
 *         Usually one simulation step requires multiple runs of MLCP solver.
 *         ConstraintManager tries to utilize the previous matrix M and vector
 *         q as much as possible for efficiency.
 *
 *  An Example Usage
 *
 *    // ...
 *
 *    ConstraintManager cm;
 *    ContactManager    ct;
 *    JointManager      jm;
 *
 *    for (ConvexRigidBody* p : allObjects) {
 *        cm.registerConvexRigidBody(p);
 *    }
 *
 *    for (ConvexRigidBody* p : collidableObjects) {
 *        ct.registerConvexRigidBody(p);
 *    }
 *
 *    for (Joint* p : permanentJoints) {
 *        jm.registerJoint(p, cm);
 *    }
 *
 *    ct.initializeAABB();
 *
 *    // ...
 *
 *    // In One simulation loop:
 *    void step(double deltaT)
 *    {
 *
 *        // You can register and unregister ConvexRigidBody here.
 *        // Also external forces and torques for this step must be set here.
 *
 *        cm.initializeStep(deltaT);
 *
 *        ct.lockForNextStep(deltaT);
 *        ct.udpateActiveContacts();
 *        ct.registerActiveConstraints(ce);
 *
 *        jm.update(1.0/deltaT);
 *
 *        cm.update();
 *
 *        ct.discoverNewContacts(true);
 *        ct.registerNewConstraints(ce);
 *
 *        cm.update();
 *        cm.commitUpdate();
 *
 *        cm.terminateStep();
 *        ct.unlockAfterStep(ce);
 *
 *        // Render Objects to the screen here.
 *
 *    }
 *
 *    for (Joint* p : permanentJoints) {
 *        jm.unregisterJoint(p, cm);
 *    }
 *
 *    for (ConvexRigidBody* p : allObjects) {
 *        cm.unregisterConvexRigidBody(p);
 *    }
 *
 *    for (ConvexRigidBody* p : collidableObjects) {
 *        ct.unregisterConvexRigidBody(p);
 *    }
 */
class ConstraintManager {

  public:

    /** @brief default constructor
     *         mMLCP_cfm(1.0e-08),
     *         mMLCP_maxIter(3),
     *         mMLCP_kgs(4),
     *         mMLCP_ksm(3),
     *         mMLCP_tolerance(1.0e-10),
     *         mMemMgr(200, 50){;}
     */
    inline ConstraintManager();


    /** @brief constructor
     *
     *  @param  MLCP_cfm        (in): Constraint Force Mixing. 
     *  @param  MLCP_maxIter    (in): Max iterations of outer loop
     *  @param  MLCP_kgs        (in): Max iterations of PGS inner loop
     *  @param  MLCP_ksm        (in): Max iterations of SM inner loop
     *  @param  MLCP_tolerance  (in): Threshold against the residual value
     *  @param  colsMinCapacity (in): Initial memory capacity in number of 
     *                                columns/rows of matrix
     *  @param  colsReserve     (in): Extra margin set aside in memory in
     *                                number of columns/rows of matrix.
     */
    inline ConstraintManager(
        double MLCP_cfm,
        long   MLCP_maxIter,
        long   MLCP_kgs,
        long   MLCP_ksm,
        double MLCP_tolerance,
        long   colsMinCapacity, 
        long   colsReserve
    );


    /** @brief destructor */
    inline ~ConstraintManager();

    /** @brief register a ConvexRigidBody. */
    inline void registerConvexRigidBody(ConvexRigidBody* p);

    /** @brief unregister ConvexRigidBody. */
    inline void unregisterConvexRigidBody(ConvexRigidBody* p);

    /** @brief register a Jacobian constraint */
    inline void registerConstraint(JacobianConstraint* c);

    /** @brief unregister a Jacobian constraint */
    inline void unregisterConstraint(JacobianConstraint* c);

    /** @brief resets the constraints and objects before one simulation step */
    void initializeStep(const double deltaT);

    /** @brief update the geometric configuration according to the external
     *         forces and torques applied to each ConvexRigidBody under given
     *         constraints
     */
    void update();

    /** @brief make the geometric configuration generated by the latest
     *        update() permanent.
     */
    void commitUpdate();

    /** @brief returns the duration of the current simulation step.
     *         Used by the subcomponents of ConstraintManager such as
     *         CollisionManager.
     */
    inline double deltaT() const;

    /** @brief terminates the current simulation step */
    void terminateStep();

    /** @brief helper utilities */
    inline bool isBody1Fixed(JacobianConstraint* c);
    inline bool isBody2Fixed(JacobianConstraint* c);

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif
    /** @class ColumnManager
     *
     *  @brief This is a helper subcomponent.
     *         It manages the indices for the rows & columns of the matrices 
     *         and vectors for MLCP. ConstraintManager calls MLCP multiple 
     *         times with slightly different M and q. During that time some
     *         rows/cols are added and some others are deleted.
     *         This component keeps track of what indices are active and what
     *         are not.
     *         The indices are of 1-origin.
     */
    class ColumnManager {

      public:

        inline ColumnManager();

        inline ~ColumnManager();

        /** @brief it allocates an available index number for the new 
         *         row/column
         */
        inline long allocateOne();

        /** @brief it receives an index from the user for reuse */
        inline void returnOne(long i);

        /** @brief reset the allocation. All the indices allocated are
         *         to be reclaimed.
         */
        inline void reset();

        /** @brief returns the highest index currently allocated */
        inline long highestIndex() const;

#ifdef UNIT_TESTS
      public:
#else
      private:
#endif
        long       mCurrentHighest;
        list<long> mDisusedIndices;

    };


    /** @class MemoryManager 
     *
     *  @brief Memory management subcomponent for ConstraintManager.
     *         It handles the memory required for MLCP to solve Mz = q.
     *         M is a symmetric PD matrix, and z & q are vectors.
     *         MLCP also requires work memory pointed by mWM.
     *         ConstraintManager calls MLCP multiple times with slightly 
     *         different M and q. It manages the memory allocation and the
     *         reuse of the contents of M and q to avoid recalculation of them.
     *         The memory for the matrix is allocated only for the lower
     *         diagnoal part in row major. This means we can use the same
     *         memory and its content for any upper submatrix, and we can
     *         copy the content 'as is' to the upper part of the bigger matrix
     *
     *  @dependencies It internally calls malloc() and free().
     */
    class MemoryManager {

      public:

        /** @brief constructor
         *
         *  @param colsMinCapacity (in): minimum number of columns/rows
         *                               allocated regardless of the requested
         *                               size.
         *  @param colsReserve     (in): extra number of columns/rows set aside
         *                               beyond the requested size to avoid
         *                               reallocation in the future.
         *                               For example if the current allocated
         *                               size is 50, the requested size is
         *                               64. In this case we have to reallocate
         *                               memory to accommodate the new size.
         *                               When we reallocate if colsReserve is
         *                               20, then 84 columns are allocated.
         */
        inline MemoryManager(long colsMinCapacity, long colsReserve);

        /** @brief descructor */
        inline ~MemoryManager();

        /** @brief sets the size of the constraints at the beginning of 
         *         a simulation step. 
         */
        inline void setSizeAndInit(long colsNewEnd);

        /** @brief updates the size of the constraints during a simulation step
         */
        inline void adjustSize(long colsNewEnd);

        /** @brief fills the specified column with zero in the internal 
         *         matrix mM
         *
         *  @param in (in): The index of the column
         *
         *
         *       mM (symmetric)
         *    +-----------------+
         *    |*                |
         *    |* *              |
         *    |* * *            |
         *    |* * * *          | <= i-1
         *    |0 0 0 0 0        | <= i
         *    |* * * * * *      |
         *    |* * * * * * *    |
         *    |* * * * * * * *  |
         *    |* * * * * * * * *|
         *    +-----------------+
         */
        inline void clearMcol(long i);

        /** @brief sets the type of constraint for the specified index.*/
        inline void setMLCPtype(long i, enum MLCP::columnType t);

        /** @brief returns the type of constraint for the specified index.*/
        inline enum MLCP::columnType MLCPtype(long i);


        inline void* M()      { return mM; }
        inline void* q()      { return mq; }
        inline void* z()      { return mz; }
        inline void* z_low()  { return mz_low;  }
        inline void* z_high() { return mz_high; }
        inline void* WM()     { return mWM; }
        inline vector<enum MLCP::columnType>&  MLCPtypes() {return mMLCPtypes;}

#ifdef UNIT_TESTS
      public:
#else
      private:
#endif

        inline void allocateMemory(
             long cols, void** M, void** q, 
             void** z, void** z_low, void** z_high, void** WM
        );

        inline void freeMemory(
            void** M, void** q, 
            void** z, void** z_low, void** z_high, void** WM
        );

        /** @brief pads the Martix M and the vectors q, z, z_low, and z_high 
         *         with zero  in the specified range.
         *
         *  @param indexLow (in): the first index of the range to be filled 
         *                        with zero.
         *
         *  @param indexHigh(in): the last index of the range to be filled with
         *                        zero plus 1. I.e. The first index of the 
         *                        next range.
         *
         *  @param mdst     (in): Matrix M
         *
         *  @param qdst     (in): RHS vector q
         *
         *  @param zdst     (in): Solution vector z
         *
         *  @param z_lowdst (in): Box range vector low side
         *
         *  @param z_highdst(in): Box range vector high side
         */
        inline void fillWithZero(
            long indexLow, long indexHigh, char* Mdst, char* qdst, 
            char* zdst, char* z_lowdst, char* z_highdst            );


        /** @brief copies the contents of the internal matrix and the vectors
         *         to the new locations specified by the parameters.
         *
         *    Range: [1, indexHigh]
         *
         *       mM      => Mdst
         *       mq      => qdst
         *       mz      => zdst
         *       mz_low  => z_lowdst
         *       mz_high => z_highdst
         */
        inline void copyInExisting(
             long indexHigh, void* Mdst, void* qdst, 
             void* zdst, void* z_lowdst, void* z_highdst );


        inline void calcMemorySizes(long n, long& Mat, long& Vec, long& WM);


        /** @brief [1, mColsCapacity] can be used with the current memory */
        long  mColsCapacity;

        /** @brief [1, mColsActualEnd] are actively used currently.*/
        long  mColsActualEnd;

        /** @brief extra columns allocated on top of the requested size as
         *         margin.
         */
        long  mColsReserve;

        /** @brief at least memory for mColsMinimumCapacity columns is
         *         allocated.*/
        long  mColsMinimumCapacity;

        /** @brief memory for M matrix for MLCP */
        void* mM;

        /** @brief memory for q vector for MLCP */
        void* mq;

        /** @brief memory for z vector for MLCP */
        void* mz;

        /** @brief memory for z_low vector for MLCP */
        void* mz_low;

        /** @brief memory for z_high vector for MLCP */
        void* mz_high;

        /** @brief memory required for MLCP solver */
        void* mWM;

        /** @brief vector to specify the type of column*/
        vector<enum MLCP::columnType>  mMLCPtypes;

    };


    /** @brief StateManager
     *         This is a helper subcomponent.
     *
     *         It manages the states of ConvexRigidBodies and 
     *         JacobianConstraints for 
     *         efficient recalculation of M and q for MLCP. A simulation step
     *         for the time step t + Δt may take multiple calculations of MLCP 
     *         in the form of Mz = q. Those calculations use the same geometric
     *         configuration for time step t, but with slightly different set
     *         of constraints.
     *
     *         Across those calculations the matrix M and the RHS vector q 
     *         share many elements, as they share many constraints.
     *         Realistically, the first calculation may have 100 constraints 
     *         and the second calculation has 10 new constraints on top of 
     *         those 100.
     *        
     *         To avoid recalculation of those shared, we keep the contents of
     *         M and q and we keep track of the state of the elements.
     *
     *         An element M(i,j) is generated from the jacobian constraints
     *         for i and j as well as the mass/inertia of the one or two 
     *         common ConvexRigidBodies.
     *         The mass/inertia for a ConvexRigidBody is expressed by 
     *         6x6 mass/inertia matrix.
     *        
     *         An element q(i) depends on the jacobian constraint for i as 
     *         well as the hypothetical velocity vector U(s) calculated only
     *         with external forces and torques applied to the 
     *         ConvexRigidBody s at the time t + Δt.
     *         We should calculate U(s) only for the ConvexRigidBodies involved
     *         in the constraints.
     *
     *         ConvexRigidBody is in one of the following states.
     *
     *            UNREGISTERED  :  Not registered to ConstraintManager yet or
     *                             it has been unregistered from it.
     *            FREE          :  Not involved in any constraints.
     *                             The movement is governed by the external
     *                             forces and torques only.
     *
     *            RESET         :  Involved in a constraint.
     *                             U(s) vector is not calculated.
     *                             None of the relevant constraints has been
     *                             processed (i.e. element M(i,*) and q(i) 
     *                             have not been calculated for constraint i
     *                             involving this ConvexRigidBody).
     *
     *            MANAGED       :  Involved in a constraint.
     *                             U(s) vector has been calculated
     *                             All of the relevant constraints have been
     *                             processed for M and q.
     *                            
     *            NEEDS_UPDATE  :  Involved in a constraint.
     *                             U(s) vector has been calculated.
     *                             Some of the relevant constraints have been
     *                             processed, however some new constraints
     *                             involving this ConvexRigidBody have been 
     *                             added since the last MLCP calculation.
     *         
     *         Also each ConvexRigidBody maintains two lists of participating 
     *         constraints
     *
     *             mConstraintsNotCheckedIn - list of constraints for which 
     *                             M(i,*) and q(i) have not been calculated.
     *                                        
     *             mConstraintsCheckedIn    - list of constraints for which
     *                             M(i,*) and q(i) have been calculated.
     *
     *         If it is in FREE state, there must be no element in both lists.
     *         If it is in MANAGED state, there must be no element in 
     *         mJacobianConstraintsNotCheckedIn and at least one element in 
     *         mJacobianConstraintsCheckedIn.
     *         If it is in RESET or NEEDS_UPDATE state, at least one element
     *         must be in mJacobianConstraintsNotCheckedIn.
     *
     *         It is not efficient to calculate M from the matrix
     *         multiplications Jt·Minv·J directly.
     *         Also, looping over columns and rows of M to generate its
     *         elements as follows is not effieicnet either, as  many of M(i,j)
     *         are zero.
     *        
     *             for (i : columns of M)
     *                 for (j>=i : rows of M)
     *                     M(i,j) = 0
     *                     for (s : common ConvexRigidBody of i and j)
     *                         M(i,j) += Jt(i)Minv(s)J(j)
     *
     *         Instead, we loop over the list of ConvexRigidBodies that are 
     *         either in  RESET or NEEDS_UPDATE state.
     *             for (s : ConvexRigidBodies in RESET or NEEDS_UPDATE)
     *                 for (i : mConstraintsNotCheckedIn(s))
     *                     for (j : mConstraintsNotCheckedIn(s) not explored by
     *                              i, including i, 
     *                              and mConstraintsCheckedIn(s) )
     *                             
     *                         M(i,j) += Jt(i)Minv(s)J(j)
     *
     *         The purpose of mConstraintsNotCheckedIn mConstraintsCheckedIn is
     *         to realize this process.
     *
     *         JacobianConstraint is in one of the following states.
     *
     *             UNREGISTERED   : Not registered to ConstraintManager yet, or
     *                              if has been unregistered from it.
     *
     *             NOT_CHECKED_IN : M(i,*) and q(i) have not been calculated 
     *
     *             CHECKED_IN     : M(i,*) and q(i) have been calculated 
     */
    class StateManager {

      public:

        inline void registerBody(ConvexRigidBody* p);
        inline void unregisterBody(ConvexRigidBody* p);
        inline void registerConstraint(JacobianConstraint* c);
        inline void linkConstraintToBody(
                                    JacobianConstraint* c, ConvexRigidBody* p);
        inline void unregisterConstraint(JacobianConstraint* c);
        inline void unlinkConstraintFromBodyNotCheckedIn(
                                    JacobianConstraint* c, ConvexRigidBody* p);
        inline void unlinkConstraintFromBodyCheckedIn(
                                    JacobianConstraint* c, ConvexRigidBody* p);
        inline void moveBodyFromFreeToReset(ConvexRigidBody* p);
        inline void moveBodyFromResetToFree(ConvexRigidBody* p);
        inline void moveBodyFromResetToManaged(ConvexRigidBody* p);
        inline void moveBodyFromResetToNeedsUpdate(ConvexRigidBody* p);
        inline void moveBodyFromManagedToFree(ConvexRigidBody* p);
        inline void moveBodyFromManagedToReset(ConvexRigidBody* p);
        inline void moveBodyFromManagedToNeedsUpdate(ConvexRigidBody* p);
        inline void moveBodyFromNeedsUpdateToReset(ConvexRigidBody* p);
        inline void moveBodyFromNeedsUpdateToFree(ConvexRigidBody* p);
        inline void moveBodyFromNeedsUpdateToManaged(ConvexRigidBody* p);
        inline void resetConstraints();
        inline void checkInConstraints();
        inline long numConstraintsNotCheckedIn() const;
        inline long numConstraintsCheckedIn() const;
        inline list<JacobianConstraint* >&  constraintsNotCheckedIn();
        inline list<JacobianConstraint* >&  constraintsCheckedIn();
        inline list<ConvexRigidBody* >& bodiesFree();
        inline list<ConvexRigidBody* >& bodiesReset();
        inline list<ConvexRigidBody* >& bodiesManaged();
        inline list<ConvexRigidBody* >& bodiesNeedsUpdate();

        inline bool isBody1Fixed(JacobianConstraint* c);
        inline bool isBody2Fixed(JacobianConstraint* c);


      private:

        list<JacobianConstraint* > mConstraintsNotCheckedIn;
        list<JacobianConstraint* > mConstraintsCheckedIn;
        list<ConvexRigidBody* >    mBodiesFree;
        list<ConvexRigidBody* >    mBodiesReset;
        list<ConvexRigidBody* >    mBodiesManaged;
        list<ConvexRigidBody* >    mBodiesNeedsUpdate;

        static const char*         EXCEPTION_MSG_01;

    };

    enum state {
        NONE,
        IN_SIMULATION_STEP,
        OUT_OF_SIMULATION_STEP,
        END
    };

    /** @brief state of the simulation */
    enum state    mState;

    /** @brief the length of the current simulation step.*/
    double        mDeltaT;

    /** @brief simulation parameters used by MLCP solver */
    double        mMLCP_cfm;       // Constraint Force Mixing. 
    long          mMLCP_maxIter;   // Max iterations of outer loop
    long          mMLCP_kgs;       // Max iterations of PGS inner loop
    long          mMLCP_ksm;       // Max iterations of SM inner loop
    double        mMLCP_tolerance; // Threshold against the residual value
                                   // to declare solution convergence.
    // Helper subcomponents
    MemoryManager mMemMgr;
    ColumnManager mColMgr;
    StateManager  mStateMgr;

    static const char* EXCEPTION_MSG_01;

#ifdef UNIT_TESTS
    void printDebug();
#endif
};


/*************************************************************/
/*                                                           */
/*  Implementation of ColumnManager inline functions BEGIN   */
/*                                                           */
/*************************************************************/


ConstraintManager::ColumnManager::ColumnManager():mCurrentHighest(0){;}


ConstraintManager::ColumnManager::~ColumnManager(){;}


long ConstraintManager::ColumnManager::allocateOne()
{
    if (mDisusedIndices.empty()) {

        return ++mCurrentHighest;
    }
    else {

        long allocated = mDisusedIndices.back();
        mDisusedIndices.pop_back();

        return allocated;
    }
}


void ConstraintManager::ColumnManager::returnOne(long i)
{
    if (i == mCurrentHighest) {

        mCurrentHighest--;
    }
    else{

        mDisusedIndices.push_back(i);
    }
}


void ConstraintManager::ColumnManager::reset()
{ 
    mCurrentHighest = 0;
    mDisusedIndices.clear();
}


long ConstraintManager::ColumnManager::highestIndex() const
{
    return mCurrentHighest;
}


/*************************************************************/
/*                                                           */
/*  Implementation of ColumnManager inline functions END     */
/*                                                           */
/*************************************************************/


/*************************************************************/
/*                                                           */
/*  Implementation of MemoryManager inline functions BEGIN   */
/*                                                           */
/*************************************************************/


ConstraintManager::MemoryManager::MemoryManager(
    long colsMinCapacity,
    long colsReserve
):
    mColsCapacity(0),
    mColsActualEnd(0),
    mColsReserve(colsReserve),
    mColsMinimumCapacity(colsMinCapacity),
    mM(nullptr),
    mq(nullptr),
    mz(nullptr),
    mz_low(nullptr),
    mz_high(nullptr),
    mWM(nullptr){;}


ConstraintManager::MemoryManager::~MemoryManager()
{
    freeMemory(&mM, &mq, &mz, &mz_low, &mz_high, &mWM);
}


void ConstraintManager::MemoryManager::setSizeAndInit(long colsNewEnd)
{

    mColsActualEnd = colsNewEnd;

    long colsNewCapacity = std::max( mColsActualEnd + mColsReserve, 
                                     mColsMinimumCapacity          );

    if (colsNewCapacity > mColsCapacity) {

        // Reallocate memory and fill with zero.
        mColsCapacity = colsNewCapacity;

        freeMemory(&mM, &mq, &mz, &mz_low, &mz_high, &mWM);

        allocateMemory(mColsCapacity, &mM, &mq, &mz, &mz_low, &mz_high, &mWM);

    }

    // Elements in [1, mColsCapacity] are filled with 0.
    fillWithZero(1, mColsCapacity+1, (char*)mM, (char*)mq, (char*)mz,
                     (char*)mz_low, (char*)mz_high                     );

    mMLCPtypes.resize(mColsCapacity, MLCP::NOT_USED);

}


void ConstraintManager::MemoryManager::adjustSize(long colsNewEnd)
{
    if (colsNewEnd <= mColsCapacity) {

        // The range [1, colsNewEnd] will be used.
        // Elements in [colsNewEnd+1, mColsCapacity] are filled with 0.
        // 
        //  1     colsNewEnd
        //  |     |
        //  |     |       mColsCapacity
        //  |     |       |
        //  *******00000000
        fillWithZero( colsNewEnd+1, mColsCapacity+1, 
               (char*)mM, (char*)mq, (char*)mz, (char*)mz_low, (char*)mz_high);
               
        mColsActualEnd = colsNewEnd;

    }
    else {// colsNewEnd > mColsCapacity >= mColsActualEnd) {
        long colsCapacityNew = std::max( colsNewEnd + mColsReserve, 
                                         mColsMinimumCapacity       );

        void *Mnew, *qnew, *znew, *z_lownew, *z_highnew, *WMnew;
        allocateMemory( colsCapacityNew, &Mnew, &qnew, 
                        &znew, &z_lownew, &z_highnew, &WMnew );

        // The range [1, mColsActualEnd] will be copied to the new memory.
        // Elements in [mColsActualEnd+1, colsCapacityNew] are filled with 0.
        // 
        //  1     mColsActualEnd
        //  |     |   colsNewEnd
        //  |     |   |   colsCapacityNew
        //  |     |   |   |
        //  *******00000000

        // The elements in [1, mColsActualEnd] are copied 
        // from the current M and q.
        // The elements in [mColsActualEnd+1, colsCapacityNew] 
        // are filled with 0.
        copyInExisting(mColsActualEnd, Mnew, qnew, znew, z_lownew, z_highnew);

        fillWithZero( mColsActualEnd+1, colsCapacityNew+1, 
               (char*)mM, (char*)mq, (char*)mz, (char*)mz_low, (char*)mz_high);

        mColsCapacity  = colsCapacityNew;
        mColsActualEnd = colsNewEnd;

        freeMemory(&mM, &mq, &mz, &mz_low, &mz_high, &mWM);

        mM      = Mnew;
        mq      = qnew;
        mz      = znew;
        mz_low  = z_lownew;
        mz_high = z_highnew;
        mWM     = WMnew;

    }

    mMLCPtypes.resize(mColsCapacity, MLCP::NOT_USED);

}

/** @brief fills the specified column with zero in the internal matrix mM
 *
 *  @param in (in): The index of the column
 *
 *
 *       mM (symmetric)
 *    +-----------------+
 *    |*                |
 *    |* *              |
 *    |* * *            |
 *    |* * * *          | <= i-1
 *    |0 0 0 0 0        | <= i
 *    |* * * * * *      |
 *    |* * * * * * *    |
 *    |* * * * * * * *  |
 *    |* * * * * * * * *|
 *    +-----------------+
 */
void ConstraintManager::MemoryManager::clearMcol(long i)
{

    // memset 0 to M.
    long bytesLow  = SymMat::requiredMemorySizeInBytes(i-1);
    long bytesHigh = SymMat::requiredMemorySizeInBytes(i);

    memset((char*)mM + bytesLow, (int)0, bytesHigh - bytesLow);


    bytesLow  = VarVec::requiredMemorySizeInBytes(i-1);
    bytesHigh = VarVec::requiredMemorySizeInBytes(i);

    memset((char*)mq      + bytesLow, (int)0, bytesHigh - bytesLow);
    memset((char*)mz      + bytesLow, (int)0, bytesHigh - bytesLow);
    memset((char*)mz_low  + bytesLow, (int)0, bytesHigh - bytesLow);
    memset((char*)mz_high + bytesLow, (int)0, bytesHigh - bytesLow);

}


void ConstraintManager::MemoryManager::allocateMemory(
    long cols, void** M, void** q, 
    void** z,  void** z_low, void** z_high, void** WM
) {

    long numMat,  numVec,  numWM;

    calcMemorySizes( cols,numMat, numVec, numWM );

    (*M) = malloc(numMat);
    if ((*M)==nullptr) {
        throw std::bad_alloc(); 
    }

    (*q) = malloc(numVec);
    if ((*q)==nullptr) {
        free(*M);
        throw std::bad_alloc(); 
    }

    (*z) = malloc(numVec);
    if ((*z)==nullptr) {
        free(*M);
        free(*q);
        throw std::bad_alloc(); 
    }

    (*z_low) = malloc(numVec);
    if ((*z_low)==nullptr) {
        free(*M);
        free(*q);
        free(*z);
        throw std::bad_alloc(); 
    }

    (*z_high) = malloc(numVec);
    if ((*z_high)==nullptr) {
        free(*M);
        free(*q);
        free(*z);
        free(*z_low);
        throw std::bad_alloc(); 
    }

    (*WM) = malloc(numWM);
    if ((*WM)==nullptr) {
        free(*M);
        free(*q);
        free(*z);
        free(*z_low);
        free(*z_high);
        throw std::bad_alloc(); 
    }
}


void ConstraintManager::MemoryManager::freeMemory(
    void** M,
    void** q,
    void** z,
    void** z_low,
    void** z_high,
    void** WM
) {
    if (*M!=nullptr) {
        free(*M);
        *M = nullptr;
    }

    if (*q!=nullptr) {
        free(*q);
        *q = nullptr;
    }

    if (*z!=nullptr) {
        free(*z);
        *z = nullptr;
    }

    if (*z_low!=nullptr) {
        free(*z_low);
        *z_low = nullptr;
    }

    if (*z_high!=nullptr) {
        free(*z_high);
        *z_high = nullptr;
    }

    if (*WM!=nullptr) {
        free(*WM);
        *WM = nullptr;
    }
}


void ConstraintManager::MemoryManager::setMLCPtype(
    long                  i,
    enum MLCP::columnType t
) {
    mMLCPtypes[i-1] = t;
}


enum MLCP::columnType ConstraintManager::MemoryManager::MLCPtype(long i) {
    return mMLCPtypes[i-1];
}


/** @brief pads the Martix M and the vectors q, z, z_low, and z_high with zero 
 *         in the specified range.
 *
 *  @param indexLow (in): the first index of the range to be filled with zero.
 *
 *  @param indexHigh(in): the last index of the range to be filled with zero
 *                        plus 1. I.e. The first index of the next range.
 *
 *  @param mdst     (in): Matrix M
 *
 *  @param qdst     (in): RHS vector q
 *
 *  @param zdst     (in): Solution vector z
 *
 *  @param z_lowdst (in): Box range vector low side
 *
 *  @param z_highdst(in): Box range vector high side
 */
void ConstraintManager::MemoryManager::fillWithZero(
    long  indexLow,
    long  indexHigh,
    char* Mdst,
    char* qdst,
    char* zdst,
    char* z_lowdst,
    char* z_highdst
) {

    long numMatLow,  numVecLow,  numWMLow;
    calcMemorySizes( indexLow-1, numMatLow, numVecLow, numWMLow );

    long numMatHigh, numVecHigh, numWMHigh;
    calcMemorySizes( indexHigh-1, numMatHigh, numVecHigh, numWMHigh );

    memset(Mdst      + numMatLow, (int)0, numMatHigh - numMatLow);
    memset(qdst      + numVecLow, (int)0, numVecHigh - numVecLow);
    memset(zdst      + numVecLow, (int)0, numVecHigh - numVecLow);
    memset(z_lowdst  + numVecLow, (int)0, numVecHigh - numVecLow);
    memset(z_highdst + numVecLow, (int)0, numVecHigh - numVecLow);

}


/** @brief copies the contents of the internal matrix and the vectors
 *         to the new locations specified by the parameters.
 *
 *    Range: [1, indexHigh]
 *
 *    mM      => Mdst
 *    mq      => qdst
 *    mz      => zdst
 *    mz_low  => z_lowdst
 *    mz_high => z_highdst
 */
void ConstraintManager::MemoryManager::copyInExisting(
    long  indexHigh,
    void* Mdst,
    void* qdst,
    void* zdst,
    void* z_lowdst,
    void* z_highdst
) {
    long numMat, numVec, numWM;
    calcMemorySizes( indexHigh, numMat, numVec, numWM );

    memcpy(Mdst,      mM,      numMat);
    memcpy(qdst,      mq,      numVec);
    memcpy(zdst,      mz,      numVec);
    memcpy(z_lowdst,  mz_low,  numVec);
    memcpy(z_highdst, mz_high, numVec);
}


void ConstraintManager::MemoryManager::calcMemorySizes(
    long  n,
    long& Mat,
    long& Vec,
    long& WM
) {
    Mat = SymMat::requiredMemorySizeInBytes(n);
    Vec = VarVec::requiredMemorySizeInBytes(n);
    WM  = MLCP::WorkMemory::requiredMemorySizeInBytes(n);
}


/*************************************************************/
/*                                                           */
/*  Implementation of MemoryManager inline functions END     */
/*                                                           */
/*************************************************************/


/*************************************************************/
/*                                                           */
/*  Implementation of StateManager inline functions BEGIN    */
/*                                                           */
/*************************************************************/


void ConstraintManager::StateManager::registerBody(ConvexRigidBody* p) {

    if(p->engineState()!=ConvexRigidBody::UNREGISTERED) {
        throw logic_error(EXCEPTION_MSG_01);
    }

    p->setEngineState(ConvexRigidBody::FREE);

    p->setBackItConstraintManager(mBodiesFree.insert(mBodiesFree.end(), p));
}


void ConstraintManager::StateManager::unregisterBody(ConvexRigidBody* p) {

    if(p->engineState()!=ConvexRigidBody::FREE) {
        throw logic_error(EXCEPTION_MSG_01);
    }

    p->setEngineState(ConvexRigidBody::UNREGISTERED);

    mBodiesFree.erase(p->backItConstraintManager());
}


void ConstraintManager::StateManager::registerConstraint(
    JacobianConstraint* c
) {

    if (c->engineState()!=JacobianConstraint::UNREGISTERED) {
        throw logic_error(EXCEPTION_MSG_01);
    }

    c->setBackItConstraintManager(mConstraintsNotCheckedIn.insert(
                                           mConstraintsNotCheckedIn.end(), c));

    if (!isBody1Fixed(c)) {

        linkConstraintToBody(c, c->body1());
    }

    if (!isBody2Fixed(c)) {

        linkConstraintToBody(c, c->body2());
    }

    c->setEngineState(JacobianConstraint::NOT_CHECKED_IN);
}


void ConstraintManager::StateManager::linkConstraintToBody(
    JacobianConstraint* c,
    ConvexRigidBody*    p
) {

    p->addConstraintNotCheckedIn(c);

    if (p->engineState()==ConvexRigidBody::FREE) {

        moveBodyFromFreeToReset(p);           
    }
    else if (p->engineState()==ConvexRigidBody::MANAGED) {

        moveBodyFromManagedToNeedsUpdate(p);
    }
}


void ConstraintManager::StateManager::unregisterConstraint(
    JacobianConstraint* c
) {

    auto state = c->engineState();

    if (state==JacobianConstraint::UNREGISTERED) {
        throw logic_error(EXCEPTION_MSG_01);
    }
    else if (state==JacobianConstraint::NOT_CHECKED_IN) {

        mConstraintsNotCheckedIn.erase(c->backItConstraintManager());

        if (!isBody1Fixed(c)) {

            unlinkConstraintFromBodyNotCheckedIn(c, c->body1());
        }
        if (!isBody2Fixed(c)) {

            unlinkConstraintFromBodyNotCheckedIn(c, c->body2());
        }

    }
    else if (state==JacobianConstraint::CHECKED_IN) {

        mConstraintsCheckedIn.erase(c->backItConstraintManager());

        if (!isBody1Fixed(c)) {

            unlinkConstraintFromBodyCheckedIn(c, c->body1());
        }
        if (!isBody2Fixed(c)) {

            unlinkConstraintFromBodyCheckedIn(c, c->body2());
        }

    }

    c->setEngineState(JacobianConstraint::UNREGISTERED);
}


void ConstraintManager::StateManager::unlinkConstraintFromBodyNotCheckedIn(
    JacobianConstraint* c,
    ConvexRigidBody*    p
) {

    p->removeConstraint(c);

    const auto state = p->engineState();

    if (p->numConstraintsNotCheckedIn()==0) {

        if (p->numConstraintsCheckedIn()==0) {

            // Prev: CheckedIn=0, NotCheckedIn=1 => RESET / NEEDS_UPDATE
            // Now:  CheckedIn=0, NotCheckedIn=0 => FREE
            if (state == ConvexRigidBody::NEEDS_UPDATE) {

                moveBodyFromNeedsUpdateToFree(p);
            }
            else if (state == ConvexRigidBody::RESET) {

                moveBodyFromResetToFree(p);
            }
            else {
                throw logic_error(EXCEPTION_MSG_01);
            }
        }
        else {
            // Prev: CheckedIn>0, NotCheckedIn=1 => NEEDS_UPDATE
            // Now:  CheckedIn>0, NotCheckedIn=0 => MANAGED
            if (state != ConvexRigidBody::NEEDS_UPDATE) {
                throw logic_error(EXCEPTION_MSG_01);
            }
            moveBodyFromNeedsUpdateToManaged(p);
        }
    }
}


void ConstraintManager::StateManager::unlinkConstraintFromBodyCheckedIn(
    JacobianConstraint* c,
    ConvexRigidBody*    p
) {

    p->removeConstraint(c);

    const auto state = p->engineState();

    if (p->numConstraintsCheckedIn()==0) {

        if (p->numConstraintsNotCheckedIn()==0) {
            // Prev: CheckedIn=1, NotCheckedIn=0 => MANAGED
            // Now:  CheckedIn=0, NotCheckedIn=0 => FREE
            if (state != ConvexRigidBody::MANAGED) {
                throw logic_error(EXCEPTION_MSG_01);
            }
            moveBodyFromManagedToFree(p);
        }
        else {
            if (state != ConvexRigidBody::NEEDS_UPDATE) {
                throw logic_error(EXCEPTION_MSG_01);
            }
            // Prev: CheckedIn=1, NotCheckedIn>0 => NEEDS_UPDATE
            // Now:  CheckedIn=0, NotCheckedIn>0 => NEEDS_UPDATE
            ;
        }
    }
}


void ConstraintManager::StateManager::moveBodyFromFreeToReset(
    ConvexRigidBody* p
) {

    p->setEngineState(ConvexRigidBody::RESET);

    mBodiesFree.erase(p->backItConstraintManager());

    p->setBackItConstraintManager(mBodiesReset.insert(mBodiesReset.end(), p));

}


void ConstraintManager::StateManager::moveBodyFromResetToFree(
    ConvexRigidBody* p
) {

    p->setEngineState(ConvexRigidBody::FREE);

    mBodiesReset.erase(p->backItConstraintManager());

    p->setBackItConstraintManager(mBodiesFree.insert(mBodiesFree.end(), p));
}                           


void ConstraintManager::StateManager::moveBodyFromResetToManaged(
    ConvexRigidBody* p
) {

    p->setEngineState(ConvexRigidBody::MANAGED);

    mBodiesReset.erase(p->backItConstraintManager());        

    p->setBackItConstraintManager(
                            mBodiesManaged.insert(mBodiesManaged.end(), p)  );
}


void ConstraintManager::StateManager::moveBodyFromResetToNeedsUpdate(
    ConvexRigidBody* p
) {

    p->setEngineState(ConvexRigidBody::NEEDS_UPDATE);

    mBodiesReset.erase(p->backItConstraintManager());        

    p->setBackItConstraintManager(
                    mBodiesNeedsUpdate.insert(mBodiesNeedsUpdate.end(), p)  );
}


void ConstraintManager::StateManager::moveBodyFromManagedToFree(
    ConvexRigidBody* p
) {

    p->setEngineState(ConvexRigidBody::FREE);

    mBodiesManaged.erase(p->backItConstraintManager());

    p->setBackItConstraintManager(mBodiesFree.insert(mBodiesFree.end(), p));
}                           


void ConstraintManager::StateManager::moveBodyFromManagedToReset(
    ConvexRigidBody* p
) {

    p->setEngineState(ConvexRigidBody::RESET);

    mBodiesManaged.erase(p->backItConstraintManager());

    p->setBackItConstraintManager(mBodiesReset.insert(mBodiesReset.end(), p));
}                           


void ConstraintManager::StateManager::moveBodyFromManagedToNeedsUpdate(
    ConvexRigidBody* p
) {

    p->setEngineState(ConvexRigidBody::NEEDS_UPDATE);

    mBodiesManaged.erase(p->backItConstraintManager());

    p->setBackItConstraintManager(
                      mBodiesNeedsUpdate.insert(mBodiesNeedsUpdate.end(), p));
}                           


void ConstraintManager::StateManager::moveBodyFromNeedsUpdateToReset(
    ConvexRigidBody* p
) {

    p->setEngineState(ConvexRigidBody::RESET);

    mBodiesNeedsUpdate.erase(p->backItConstraintManager());

    p->setBackItConstraintManager(mBodiesReset.insert(mBodiesReset.end(), p));
}                           


void ConstraintManager::StateManager::moveBodyFromNeedsUpdateToFree(
    ConvexRigidBody* p
) {

    p->setEngineState(ConvexRigidBody::FREE);

    mBodiesNeedsUpdate.erase(p->backItConstraintManager());

    p->setBackItConstraintManager( mBodiesFree.insert(mBodiesFree.end(), p) );
}                           


void ConstraintManager::StateManager::moveBodyFromNeedsUpdateToManaged(
    ConvexRigidBody* p
) {

    p->setEngineState(ConvexRigidBody::MANAGED);

    mBodiesNeedsUpdate.erase(p->backItConstraintManager());

    p->setBackItConstraintManager(
                          mBodiesManaged.insert(mBodiesManaged.end(), p)  );
}                           


void ConstraintManager::StateManager::resetConstraints() {

    vector<JacobianConstraint* > checkedInCopy(
                   mConstraintsCheckedIn.begin(), mConstraintsCheckedIn.end());
        
    for (auto* c : checkedInCopy) {

        mConstraintsCheckedIn.erase(c->backItConstraintManager());

        c->setBackItConstraintManager(mConstraintsNotCheckedIn.insert(
                                          mConstraintsNotCheckedIn.end(), c) );

        c->setEngineState(JacobianConstraint::NOT_CHECKED_IN);

        if (!isBody1Fixed(c)) {
 
           c->body1()->moveToNotCheckedIn(c);
        }

        if (!isBody2Fixed(c)) {

            c->body2()->moveToNotCheckedIn(c);
        }
    }

    vector<ConvexRigidBody*> mBodiesNeedsUpdateCopy(
                        mBodiesNeedsUpdate.begin(), mBodiesNeedsUpdate.end() );

    for (auto* p : mBodiesNeedsUpdateCopy) {

        moveBodyFromNeedsUpdateToReset(p);

    }

    vector<ConvexRigidBody*> mBodiesManagedCopy(
                           mBodiesManaged.begin(), mBodiesManaged.end() );
                                     
    for (auto* p : mBodiesManagedCopy) {

        moveBodyFromManagedToReset(p);

    }
}


void ConstraintManager::StateManager::checkInConstraints() {

    vector<JacobianConstraint* > notCheckedInCopy
          (mConstraintsNotCheckedIn.begin(), mConstraintsNotCheckedIn.end() );

    for (auto* c : notCheckedInCopy) {

        mConstraintsNotCheckedIn.erase(c->backItConstraintManager());

        c->setBackItConstraintManager(mConstraintsCheckedIn.insert(
                                             mConstraintsCheckedIn.end(), c) );

        c->setEngineState(JacobianConstraint::CHECKED_IN);

        if (!isBody1Fixed(c)) {

            c->body1()->moveToCheckedIn(c);
        }

        if (!isBody2Fixed(c)) {

            c->body2()->moveToCheckedIn(c);
        }
    }

    vector<ConvexRigidBody*> mBodiesResetCopy(
                                  mBodiesReset.begin(), mBodiesReset.end());
                                        
    for (auto* p : mBodiesResetCopy) {
        moveBodyFromResetToManaged(p);
    }

    vector<ConvexRigidBody*> mBodiesNeedsUpdateCopy(
                       mBodiesNeedsUpdate.begin(), mBodiesNeedsUpdate.end()  );

    for (auto* p : mBodiesNeedsUpdateCopy) {

        moveBodyFromNeedsUpdateToManaged(p);
    }
}


long ConstraintManager::StateManager::numConstraintsNotCheckedIn() const
{
    return mConstraintsNotCheckedIn.size();
}


long ConstraintManager::StateManager::numConstraintsCheckedIn() const
{
    return mConstraintsCheckedIn.size();
}


list<JacobianConstraint* >&
ConstraintManager::StateManager::constraintsNotCheckedIn()
{
    return mConstraintsNotCheckedIn;
}


list<JacobianConstraint* >& 
ConstraintManager::StateManager::constraintsCheckedIn()
{
    return mConstraintsCheckedIn;
}


list<ConvexRigidBody* >& ConstraintManager::StateManager::bodiesFree()
{
    return mBodiesFree; 
}


list<ConvexRigidBody* >& ConstraintManager::StateManager::bodiesReset()
{
    return mBodiesReset;
}


list<ConvexRigidBody* >& ConstraintManager::StateManager::bodiesManaged()
{
    return mBodiesManaged;
}


list<ConvexRigidBody* >& ConstraintManager::StateManager::bodiesNeedsUpdate()
{
    return mBodiesNeedsUpdate;
}


bool ConstraintManager::StateManager::isBody1Fixed(JacobianConstraint* c)
{
    if (c->body1()==nullptr) {
        return true;
    }
    return c->body1()->isFixed();
}


bool ConstraintManager::StateManager::isBody2Fixed(JacobianConstraint* c)
{
    if (c->body2()==nullptr) {
        return true;
    }
    return c->body2()->isFixed();
}


/*************************************************************/
/*                                                           */
/*  Implementation of StateManager inline functions END      */
/*                                                           */
/*************************************************************/


/*************************************************************/
/*                                                           */
/* Implementation of ConstraintManager inline functions BEGIN */
/*                                                           */
/*************************************************************/


ConstraintManager::ConstraintManager(
    double MLCP_cfm,
    long   MLCP_maxIter,
    long   MLCP_kgs,
    long   MLCP_ksm,
    double MLCP_tolerance,
    long   colsMinCapacity, 
    long   colsReserve
):
    mState(OUT_OF_SIMULATION_STEP),
    mMLCP_cfm(MLCP_cfm),
    mMLCP_maxIter(MLCP_maxIter),
    mMLCP_kgs(MLCP_kgs),
    mMLCP_ksm(MLCP_ksm),
    mMLCP_tolerance(MLCP_tolerance),
    mMemMgr(colsMinCapacity,colsReserve)
{;}


ConstraintManager::ConstraintManager():
    mState(OUT_OF_SIMULATION_STEP),
    mMLCP_cfm(1.0e-4),
    mMLCP_maxIter(10),
    mMLCP_kgs(4),
    mMLCP_ksm(3),
    mMLCP_tolerance(1.0e-4),
    mMemMgr(200, 50){;}


ConstraintManager::~ConstraintManager(){;}


void ConstraintManager::registerConvexRigidBody(ConvexRigidBody* p)
{
    if (mState!=OUT_OF_SIMULATION_STEP) {
        throw logic_error(EXCEPTION_MSG_01);
    }
    mStateMgr.registerBody(p);
}


void ConstraintManager::unregisterConvexRigidBody(ConvexRigidBody* p)
{
    if (mState!=OUT_OF_SIMULATION_STEP) {
        throw logic_error(EXCEPTION_MSG_01);
    }
    mStateMgr.unregisterBody(p);
}


void ConstraintManager::registerConstraint(JacobianConstraint* c)
{
    mStateMgr.registerConstraint(c);
}


void ConstraintManager::unregisterConstraint(JacobianConstraint* c)
{
    if (c->engineState()==JacobianConstraint::CHECKED_IN) {

        mMemMgr.setMLCPtype(c->matIndex(), MLCP::NOT_USED);
        mMemMgr.clearMcol(c->matIndex());
        mColMgr.returnOne(c->matIndex());           
    }

    mStateMgr.unregisterConstraint(c);
}

double ConstraintManager::deltaT() const { return mDeltaT; }


bool ConstraintManager::isBody1Fixed(JacobianConstraint* c)
{
    if (c->body1()==nullptr) {
        return true;
    }
    return c->body1()->isFixed();
}


bool ConstraintManager::isBody2Fixed(JacobianConstraint* c)
{
    if (c->body2()==nullptr) {
        return true;
    }
    return c->body2()->isFixed();
}


/*************************************************************/
/*                                                           */
/*  Implementation of ConstraintManager inline functions END  */
/*                                                           */
/*************************************************************/


}// namespace Makena

#endif /*_MAKENA_CONSTRAINT_MANAGER_HPP_*/

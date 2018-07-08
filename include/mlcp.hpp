#ifndef _MAKENA_MLCP_HPP_
#define _MAKENA_MLCP_HPP_

#include <array>
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <exception>
#include <stdexcept>
#include <cmath>

#include "variable_primitives.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file mlcp.hpp
 *
 * @brief Solves Mixed Linear Complementarity Problem or MLCP for
 *        a symmetric PD matrix with the algorithm presented in [MNS08].
 *
 *         |     |     | | |   | |   | |
 *         |  A  |  C  | |u|   |a|   |0|
 *         |     |     | | |   | |   | |
 *         |-----+-----|*|-| + |-| = |-|
 *         |   t |     | | |   | |   | |    
 *         |  C  |  B  | |v|   |b|   |w|  
 *         |     |     | | |   | |   | |, b >= 0, w >= 0, b·w = 0
 *
 *        In a simplified form, we express the above as Mz + q = w.
 *
 *        The algorithm is called projected Gauss-Seidel with subspace
 *        minimization proposed in MNS08.
 *        The partition of the matrix above is for illustration purpose only
 *        and both unilateral and bilateral constraints can be ordered 
 *        arbitrarily.
 *        The solutions for bilateral constraints can have limits (boxed 
 *        constraints.) The boxed constraints must be carefully used to keep
 *        the feasibility of the entire problem.
 *
 * @usage
 *
 *    SymMat M     (matrixSize, true, true); // PD matrix
 *    VarVec q     (matrixSize, true, true); // RHS vector
 *    VarVec z     (matrixSize, true, true); // The solution is to be filled.
 *    VarVec z_low (matrixSize, true, true); // The lower limits for bilateral
 *                                           // constraints
 *    VarVec z_high(matrixSize, true, true); // The lower limits for bilateral
 *                                           // constraints
 *
 *    vector<enum columnType>constraintTypes;// Fill this with matrixSize 
 *                                           // elements of an appropriate
 *                                           // type.
 *                                           // NOTE: this is std::vector and
 *                                           // it uses 0-based index.
 *        // NOT_USED : Specify this if this entry is not used.
 *        // BILATERAL_FREE : Bilateral free constraint
 *        // BILATERAL_BOXED_NOT_CLAMPED : Bilateral boxed constraint
 *        // UNILATERAL_NON_ACTIVE : Unilateral constraint
 *
 *    MLCP::WorkMemory(matrixSize,true,true); // Internal memory required for 
 *                                            // subspace minimization
 *
 *    long maxIter;     // Max number of iterations for the outer-most loop.
 *                      // Typically 5
 *    long kgs;         // Max number of iterations for one PGS.
 *                      // Typically 3
 *    long ksm;         // Max number of iterations for one subspace
 *                      // minimization. Typically 3
 *    double tolerance; // Numerical tolerance for solution convergence
 *                      // This is compared with the residual described in 
 *                      // [MNS08]. Generally within [0.001, 0.1].
 *   
 *
 *    MLCP::solve(M, q, z, z_low, z_high, matrixSize, constraintTypes, 
 *                                    maxIter, kgs, ksm, tolerance, wm );
 *    // Now the solution is in z.
 *    // Also the solution type is returned in constraintTypes for the 
 *    // unilateral and bilateral boxed constraints as follows.
 *    //    BILATERAL_BOXED_NOT_CLAMPED : The solution is between the limits.
 *    //    BILATERAL_BOXED_CLAMPED_LOW : The solution is at the lowest limit
 *    //    BILATERAL_BOXED_CLAMPED_HIGH: The solution is at the highest limit
 *    //    UNILATERAL_NON_ACTIVE : The solution is zero.
 *    //    UNILATERAL_ACTIVE : The solution is non-zero
 *
 * Please note that the subspace minimization is not invoked if the 
 * results from PGS indicates a solution for a bilateral boxed constraint
 * is clamped at either lower or higher limit.
 *
 * @reference 
 * [MNS08] An algorithm for the fast solution of symmetric linear 
 *         complementarity problems,
 *         J.L.Morales, J.Nocedal, and M.Smelyanskiy, 
 *         Numerische Mathematik, vol.111, no.2, pp.251-266, 2008. 
 */
namespace Makena {

using namespace std;


class MLCP {

  public:

    enum columnType {
        NOT_USED,
        BILATERAL_FREE,
        BILATERAL_BOXED_NOT_CLAMPED,
        BILATERAL_BOXED_CLAMPED_LOW,
        BILATERAL_BOXED_CLAMPED_HIGH,
        UNILATERAL_NON_ACTIVE,
        UNILATERAL_ACTIVE,
    };

    class WorkMemory {
      public:

        /** @brief constructor
         *
         *  @param rows         (in): number of rows/columns of the matrix.
         *
         *  @param allocateSelf (in): set to true if you want it to  allocate
         *                            internal memory by itself using malloc.
         *  @param padWithZero  (in): set to true if allocateSelf is true, and
         *                            if you want it to clear all the elements
         *                            with zero.
         */
        WorkMemory(
            const size_t rows,
            const bool   allocateSelf,
            const bool   padWithZero
        ):
            mMemAllocated(allocateSelf),
            mMcopy(rows, allocateSelf, padWithZero),
            mz_0  (rows, allocateSelf, padWithZero),
            mz_1  (rows, allocateSelf, padWithZero),
            mz_b  (rows, allocateSelf, padWithZero),
            mz_s  (rows, allocateSelf, padWithZero){;}

        /** @brief it returns the user the size of the memory required for this
         *         object in bytes.
         */
        size_t requiredMemorySizeInBytes() const {
            return  mMcopy.capacity() + 
                    mz_0.capacity() + 
                    mz_1.capacity() + 
                    mz_b.capacity() + 
                    mz_s.capacity() ;
        }
        static size_t requiredMemorySizeInBytes(const size_t rows) {
            return  SymMat::requiredMemorySizeInBytes(rows)+
                    VarVec::requiredMemorySizeInBytes(rows)*4;
        }


        /** @brief it receives from the user the allocated memory.
         *
         *  @param p           (in): pointer to the allocated memory
         *
         *  @param padWithZero (in): set to true if the user want it to clear
         *                           all the elements with zero.
         */
        void receiveMemory(const void* p, const bool padWithZero)
        {
            if (mMemAllocated) {
                throw std::logic_error("MLCP memory already allocated");
            }
            char* cp = (char*)p;
            mMcopy.receiveMemory(cp, padWithZero);
            cp += mMcopy.capacity();
            mz_0.receiveMemory(cp, padWithZero);
            cp += mz_0.capacity();
            mz_1.receiveMemory(cp, padWithZero);
            cp += mz_1.capacity();
            mz_b.receiveMemory(cp, padWithZero);
            cp += mz_b.capacity();
            mz_s.receiveMemory(cp, padWithZero);
        }

        const bool    mMemAllocated;

        /** @brief copy of the main square PD matrix used for chokesky 
         *         factorization. The result of the factorization is stored
         *         to this matrix. The contents is copied from mM before every
         *         cholesky factorization.
         */
        SymMat mMcopy;


        /** @brief the solution vector that stores the initial feasible 
         *         solution at the beginning of each iteration of the 
         *         subspace minimization.
         *         The active bits are made based on this vector.
         */
        VarVec mz_0;


        /** @brief the solution vector that stores the result of subspace SLE
         *         at each iteration of the subspace minimization.
         *         This may not be feasible.
         */
        VarVec mz_1;


        /** @brief the solution vector that stores the interpolated
         *         result of the first iteration of the subspace minimization.
         */
        VarVec mz_b;


        /** @brief the solution vector that stores the interpolated
         *         result of the latest iteration of the subspace minimization.
         */
        VarVec mz_s;

    };

    /** @brief the main function. Solves the given MLCP problem.
     *         the unilateral constraints can be placed anywhere in the 
     *         matrix. The type of constraints are specified by types vector,
     *         which is also used internally to specify the active elements
     *         for SLE for unilateral constraints.
     *
     *         It can perform kind of warm start with the previous values
     *         kept in z if M and q  do not change much from the previous call.
     *
     *         On the cold start with the unilateral constraints not aggregated
     *         this version is approximately 7% slower than the previous one.
     *         This 7% is expected to be absorbed in reconstruction of M and q.
     *
     *  @param M         (in):     Square PD matrix M.
     *
     *  @param q         (in):     RHS coefficient vector q
     *
     *  @param z         (in/out): Solution vector z.
     *
     *  @param z_low     (in):     The lowest possible value for z. Used if
     *                             it is for the boxed bilateral constraint.
     *
     *  @param z_high    (in):     The highest possible value for z. Used if
     *                             it is for the boxed bilateral constraint.
     *
     *  @param effectiveNumCols 
     *                   (in):     If the capacity of M, q, and z is larger
     *                             than the actual size of them, this specifies
     *                             the last 1-based index of them.
     *
     *  @param types     (in/out): the type of columns. This is std::vector
     *                             and hence it is in 0-based index.
     *
     *         NOT_USED:               the corresponding column/row is not used
     *                                        
     *         BILATERAL_FREE:         the corresponding column/row is a 
     *                                 bilateral constraint with free z
     *
     *         BILATERAL_BOXED_CLAMPED_LOW the corresponding column/row is a 
     *                                 bilateral boxed constraint. On exit
     *                                 this constraint is found to be clamped
     *                                 to the loweest allowed value
     *
     *         BILATERAL_BOXED_CLAMPED_HIGH the corresponding column/row is a 
     *                                 bilateral boxed constraint. On exit
     *                                 this constraint is found to be clamped
     *                                 to the highest allowed value
     *
     *         BILATERAL_BOXED_NOT_CLAMPED: 
     *                                 the corresponding column/row is a 
     *                                 bilateral boxed constraint. On exit
     *                                 this constraint is found to be between
     *                                 the specified boxed range.
     *
     *         UNILATERAL_NON_ACTIVE:  the corresponding column/row is a 
     *                                 unilateral constraint. On exit
     *                                 this constraint is found to be
     *                                 non-active.
     *
     *         UNILATERAL_ACTIVE:      the corresponding column/row is a 
     *                                 unilateral constraint. On exit
     *                                 this constraint is found to be active.
     *
     *  @param maxIter   (in): maximum number of iterations for the outer loop.
     *
     *  @param kgs       (in): maximum number of iterations for one PGS.
     *
     *  @param ksm       (in): maximum number of iterations for one 
     *                         subspace minimization.
     *  @param tolerance (in): convergence tolerance for residual.
     *
     *  @param wm        (in): working memory required for cholesky 
     *                         factorization of M and the subspace minimization
     */
    static long solve(
         const SymMat&            M,
         const VarVec&            q,
         VarVec&                  z,
         const VarVec&            z_low,
         const VarVec&            z_high,
         const long               effectiveNumCols,
         vector<enum columnType>& types,
         const long               maxIter,
         const long               kgs,
         const long               ksm,
         const double             tolerance,
         WorkMemory&              wm
    );


#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    /** @brief tolerance to be considered zero.
     *
     *         The minimum positive normalized number for a float is 1.18*e-38.
     *         For a double it is 2.225074e-308.
     */
    constexpr static const double EPSILON_LCP_ZERO = +1.0e-100;


    /** @brief projected Gauss-Seidel iterative solver for MLCP.
     */
    static void pgs(
        const SymMat&            M,
        const VarVec&            q,
        VarVec&                  z,
        const VarVec&            z_low,
        const VarVec&            z_high,
        const long               num,
        vector<enum columnType>& types,
        const long               kgs
    );


    /** @brief checks if any of the current solution for bilateral boxed 
     *         constraints is clamped to the lower or the higher limit, in
     *         which case subspace minimization cannot be invoked.
     */
    static bool isEligibleForSubspaceMinimization(
        vector<enum columnType>& types,
        const long               num
    );

    /** @brief improve the solution by solving SLE of the problem in a
     *         subspace and then interpolating the solution with the base
     *         solution found by PGS.
     */
    static bool minimizeSubspace(
        const SymMat&            M,
        const VarVec&            q,
        VarVec&                  z,
        const VarVec&            z_low,
        const VarVec&            z_high,
        const long               num,
        vector<enum columnType>& types,
        const long               kgs,
        WorkMemory&              wm
    );

    /** @brief construct a boolean array for the current solution stored in z_0
     *         in which true indicates the element is active.
     */
    static void findActiveBits_z_0(
        const long               num,
        vector<enum columnType>& types,
        const VarVec&            z_0,
        const VarVec&            z_low,
        const VarVec&            z_high
    );


    /** @brief finds the value alpha for the following interpolation
     *
     *         z = (1 - alpha) * z_0 + alpha * z_1
     *         
     *         where z_0 is the latest feasible solution and 
     *         z_1 is the solution found by SLE with Cholesky factorization
     *         which may not be feasible.
     *         We try to raise alpha as much as possible as long as z
     *         stays feasible.
     *
     *  @param  feasible (out): true z_1 itself is feasible
     *
     *  @return alpha value found
     */
    static double findAlpha_z_0_z_1(
        const long               num,
        vector<enum columnType>& types,
        bool&                    feasible,
        VarVec&                  z_0,
        VarVec&                  z_1
    );


    /** @brief performs the following interpolation.
     *         This is used to find z_b or z_s.
     *
     *         zOut = (1 - alpha) * z_0 + alpha * z_1
     *         
     *  @param  alpha (in):  alpha
     *
     *  @param  zOut  (out): interpolated solution vector
     */
    static void interpolate_z_0_z_1(
        const double             alpha,
        VarVec&                  zOut,
        const long               num,
        VarVec&                  z_0,
        VarVec&                  z_1,
        vector<enum columnType>& types,
        const VarVec&            z_low,
        const VarVec&            z_high

    );


    /** @brief calculates phi value of the given solution
     *         This is used to find z_b or z_s.
     *
     *                t        t
     *         phi = z  M z + q  z
     *   
     *         The smaller the values the better the solution
     *         
     *  @param  a (in): input solution vector
     *
     *  @param  zOut  (out): interpolated solution vector
     */
    static double phi(
        const SymMat&            M,
        const VarVec&            q,
        const long               num,
        vector<enum columnType>& types,
        const VarVec&            z
    );


    /** @brief calculates the error with the following formula:
     *
     *          Ra = max {Au + Cv + a}
     *          Rb = max { min {v_i ,  [C^tu + Bv + b]_i} }
     *          Rc = max { max {0,   - [C^tu + Bv + b]_i} }
     *
     *          Error = max{ ||Ra|| * coeff_a,
     *                       ||Rb|| * coeff_b,
     *                       ||Rc|| * coeff_c  }
     *
     * NOTE: [MNS08] states Rb uses min, instead of max, but it seems wrong.
     */
    static double residual(
        const SymMat&            M,
        const VarVec&            q,
        VarVec&                  z,
        const long               num,
        vector<enum columnType>& types,
        const double             coeff_a,
        const double             coeff_b,
        const double             coeff_c
    );


    /** @brief solves the simultaneous linear equations of the following form.
     *
     *          M'z' + q' = 0;
     *
     *              |A  |C'|       |a |        |u |
     *          M'= |---+--|  q' = |--|   z' = |--|
     *              |C't|B'|       |b'|        |v'|
     *
     *        where B', C', b', and v' are made by removing non-active 
     *        rows/columns from B, C, b, and v respectively.
     *        The non-active rows/columns are specified by the active bit array
     *        It uses Cholesky factorization: M' = L * L^t.
     *        It overwrites the lower triangular matrix to M'.
     *        Then it performs forward and backward substitution to get 
     *        solution into x.
     *
     *        (L·L^t) z' + q' = 0  
     *      
     *        Let L^t·z' = y, then
     *        L·y   + q'= 0  (forward substitution)
     *        and
     *        L^t·z' - y = 0  (backward substitution)
     * 
     *        It uses mMcopy. It first copies mM into mMcopy and then 
     *        factroize it in-place in mMcopy.
     *        Then the active elements are specified by mActiveBits.
     *        The solution is stored in mz_1.
     */
    static void solveSLECholeskyInto_z_1(
        SymMat&                  M,
        const VarVec&            q,
        const long               num,
        vector<enum columnType>& types,
        VarVec&                  z_1
    );


    /** @brief decomposes mMcopy into L*Lt. This algorithm uses only the lower 
     *         diagonal part.
     */
    static void decompose_cholesky_submatrix(
        const long               num,
        vector<enum columnType>& types,
        SymMat&                  M
    );


    /** @brief solves the simultaneous linear equations  L'y + q' = 0,
     *         which is a subproblem of  L'y + q' = 0,  by forward 
     *         substitution. L' is a lower diagonal matrix of mMcopy, and
     *         y is mz_1.
     */
    static void solve_lower_diagonal_submatrix(
        const VarVec&            q,
        const long               num,
        vector<enum columnType>& types,
        const SymMat&            M,
        VarVec&                  z_1
    );

    /** @brief solves the simultaneous linear equations  U'z' - y = 0, 
     *         which is a subproblem of  Uz' - y = 0, by backward substitution.
     *         U' is an upper diagonal matrix, which is a transpose of mMcopy,
     *         and z' and y are mz_1.
     *         The non-active elemens of mz_1 are set to zero.
     */
    static void solve_upper_diagonal_submatrix(
        const long               num,
        vector<enum columnType>& types,
        const SymMat&            M,
        VarVec&                  z_1
    );

};


class SLECholesky {

  public:

    /** @brief solve the following SLE with Chokesky factorization.
     *              Ax + b = 0
     *          where A is a symmetric PD matrix.
     *
     *  @param A      (in): PD matrix A
     *
     *  @param b      (in): Vector b
     *
     *  @param x      (in): Solution vector x
     *
     *  Decompose A into L·L^t, then the problem becomes:
     *
     *        (L·L^t) x + b = 0  
     *      
     *  Let L^t·x = y, then
     *
     *        L·y + b= 0.    (forward substitution)
     *  And
     *        L^t·x - y = 0. (backward substitution)
     */
    inline static void solve(SymMat& A, const VarVec& b, VarVec& x)
    {
        long dim = A.numRows();

        // Decomposition
        for (long k = 1 ; k <= dim; k++) {
            A.val(k,k) = sqrt(A.const_val(k,k));
            for (long m = k + 1 ; m <= dim ; m++) {
                A.val(m,k) = A.const_val(m,k) / A.const_val(k,k);
            }
            for (long j = k + 1; j <= dim; j++) {
                for (long m = j; m <= dim; m++) {
                    A.val(m,j)= A.const_val(m,j) -
                                A.const_val(m,k) * A.const_val(j,k);
                }
            }
        }

        // Forward substitution
        for (long i = 1; i <= dim; i++) {
            double sum = 0.0;
            for (long j = 1; j <= i-1; j++) {
                sum += (A.const_val(i,j) * x.const_val(j));
            }
            x.val(i) = ( 0.0 - b.const_val(i) - sum ) / A.const_val(i,i);
        }

        // Backward substitution
        for (long i = dim; i >= 1; i--) {
            double sum = 0.0;
            for (long j = dim; j >= i+1; j--) {
                sum += (A.const_val(j,i) * x.const_val(j));
            }
            x.val(i) = ( x.const_val(i) - sum ) / A.const_val(i,i);
        }
    }

};



}// namespace Makena


#endif/*_MAKENA_MLCP_HPP_*/

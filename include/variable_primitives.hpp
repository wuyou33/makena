#ifndef _MAKENA_VARIABLE_PRIMITIVES_HPP_
#define _MAKENA_VARIABLE_PRIMITIVES_HPP_

#include <array>
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <exception>
#include <stdexcept>
#include <cmath>
#include <memory>

#include "primitives.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file variable_primitives.hpp
 *
 * @brief implements a vector anda symmetrix square matrix of variable lengths.
 *        This isused bythe iterative MLCP solver.
 */

namespace Makena {

using namespace std;


/** @class SymMat
 *
 *  @brief square symmetric matrix.
 *         Memory is allocated only for the lower diagonal part, and hence
 *         any continuous subsequence of the memory from the beginning
 *         (of an appropriate length) corresponds to an upper square diagonal
 *         submatrix.
 *
 *  @remark Note on the internal array.
 *          The elements are stored in a linear array for the lower diagonal
 *          part as follows. The number indicates the index into the linear
 *          array.
 *        
 *          | 0                |
 *          |                  |
 *          | 1   2            |
 *          |                  | 
 *          | 3   4   5        |
 *          |                  | 
 *          | 6   7   8   9    |
 *          |                  | 
 *          |10  11  12  13  14|
 *
 *          In this arrangement, if you want to copy in a matrix to a larger
 *          matrix at the upper diagonal part, you can utilize memcpy() into
 *          the beginning of the linear array.
 *          This is implemented in copyInUpperLeftSubmatrixAndClearRest().
 *
 */
class SymMat {

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
    SymMat(const size_t rows, const bool allocateSelf, const bool padWithZero):
        mMemAllocated(allocateSelf),
        mRows(rows),
        mCapacity(sizeof(double)*(rows+1)*rows/2)
    {
        if (allocateSelf) {
            mArray = (double*)malloc(mCapacity);
            if (mArray ==nullptr) {
                throw std::bad_alloc();
            }
            if (padWithZero) {
                memset(mArray, (int)0, mCapacity);
            }
        }
    }


    /** @brief descructor */
    ~SymMat() {if (mMemAllocated) { free(mArray); } }


    /** @brief it returns the user the size of the memory required for this
     *         object in bytes.
     */
    size_t requiredMemorySizeInBytes() const { return mCapacity; }
    static size_t requiredMemorySizeInBytes(const size_t rows) {
        return sizeof(double)*(rows+1)*rows/2;
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
            throw std::logic_error("SymmetricMat memory already allocated");
        }
        mArray = (double*)p;
        if (padWithZero) {
            memset (mArray, (int)0, mCapacity);
        }
    }


    /** @brief returns the internal array */
    double* array() { return mArray; }


    /** @brief returns the pointer to the specified element.
     *
     *  @param i   (in): row    (1-origin)
     *
     *  @param j   (in): column (1-origin)
     *
     *  @return pointer to the element
     */
    double* ptr(const size_t i, const size_t j)
    {
        return &(mArray[linearizedIndex(i,j)]);
    }


    /** @brief returns the reference to the specified element.
     *
     *  @param i   (in): row    (1-origin)
     *
     *  @param j   (in): column (1-origin)
     *
     *  @return reference to the element
     */
    double& val(const size_t i, const size_t j)
    {
        return mArray[linearizedIndex(i,j)];
    }


    /** @brief returns the const reference to the specified element.
     *         Used if you want to just read the element.
     *  @param i   (in): row    (1-origin)
     *
     *  @param j   (in): column (1-origin)
     *
     *  @return const reference to the element
     */
    const double& const_val(const size_t i, const size_t j) const
    {
        return mArray[linearizedIndex(i,j)];
    }


    /** @brief copies a submatrix to the upper diagonal part.
     *
     *   Before                 Gap start: 4 Gap width: 2
     *   +---------------+      +---------+
     *   |*              |      |%        |
     *   |* *            | copy |% %      |
     *   |* * *          |  <== |% % %    |
     *   |* * * *        |      |% % % %  |
     *   |* * * * *      |      |% % % % %|
     *   |* * * * * *    |      +---------+
     *   |* * * * * * *  |
     *   |* * * * * * * *|
     *   +---------------+
     *
     *   After
     *   +---------------+
     *   |%              | +
     *   |% %            | | Part 1
     *   |% % %          | +
     *   |0 0 0 0        | + Gap 1
     *   |0 0 0 0 0      | +
     *   |% % % 0 0 %    | + Part 2 | Gap 2 | Part 3
     *   |% % % 0 0 % %  | +
     *   |0 0 0 0 0 0 0 0| + Gap 3
     *   +---------------+
     *
     *  @param UL   (in): Matrix to be copied.
     *
     */
    long posRowBytes(size_t row) { return row *(row-1)*sizeof(double)/ 2; }

    void copyInUpperLeftSubmatrixAndClearRest(
                                  SymMat& UL, long gapStart, long gapWidth)
    {
        if (gapStart > UL.mRows + 1) {
            gapStart = UL.mRows + 1;
        }
        if (gapStart <  0) {
            gapStart = 0;
        }
        if (gapWidth + UL.mRows > mRows) {
            gapWidth = mRows - UL.mRows;
        }
        if (gapWidth < 0) {
            gapWidth = 0;
        }

        char*  pBaseSrcBytes = (char*)UL.mArray;
        char*  pBaseDstBytes = (char*)mArray;
        const size_t elemSizeBytes = sizeof(double);

        // Part 1
        long const part1SizeBytes = posRowBytes(gapStart);

        memcpy( pBaseDstBytes, pBaseSrcBytes, part1SizeBytes);

        // Gap 1
        const long gap1StartBytes  = part1SizeBytes;

        const long gap1LengthBytes = posRowBytes(gapStart+gapWidth) -
                                     part1SizeBytes;


        memset( pBaseDstBytes + gap1StartBytes, (int)0, gap1LengthBytes);

        // Loop for Part 2, Gap 2, and Part 3
        const long part2LengthBytes = (gapStart-1) * elemSizeBytes;

        const long gap2LengthBytes  = gapWidth * elemSizeBytes;

        for (long i = gapStart; i <= UL.mRows; i++) {

            // Part 2 
            const long part2StartOffsetSrcBytes = posRowBytes(i);
            const long part2StartOffsetDstBytes = posRowBytes(i+gapWidth);

            memcpy( pBaseDstBytes + part2StartOffsetDstBytes, 
                    pBaseSrcBytes + part2StartOffsetSrcBytes, 
                    part2LengthBytes               );

            // Gap 2
            const long gap2StartBytes = part2StartOffsetDstBytes +
                                        part2LengthBytes;

            memset( pBaseDstBytes+gap2StartBytes, (int)0, gap2LengthBytes);

            // Part 3
            const long part3StartOffsetSrcBytes = part2StartOffsetSrcBytes + 
                                                  part2LengthBytes;

            const long part3StartOffsetDstBytes = part2StartOffsetDstBytes + 
                                                  part2LengthBytes +
                                                  gap2LengthBytes;

            const long part3LengthBytes   = (i + 1 - gapStart) * elemSizeBytes;

            memcpy( pBaseDstBytes + part3StartOffsetDstBytes, 
                    pBaseSrcBytes + part3StartOffsetSrcBytes, 
                    part3LengthBytes                          );
        }


        // Gap 3
        const long gap3StartBytes  = posRowBytes(UL.mRows + gapWidth + 1);
        const long gap3LengthBytes = posRowBytes(mRows+1) - gap3StartBytes;

        memset( pBaseDstBytes+gap3StartBytes, (int)0, gap3LengthBytes);
    }


    /** @brief returns the number of rows/columns of the matrix.*/
    size_t numCols() const { return mRows; }


    /** @brief returns the number of rows/columns of the matrix.*/
    size_t numRows() const { return mRows; }


    /** @brief copy the contents of the RHS into this. */
    SymMat& operator = (const SymMat& rhs)
    {
        if (mCapacity < rhs.mCapacity) {
            throw std::range_error("SymMat operator=");
        }
        memcpy(mArray, rhs.mArray, rhs.mCapacity);
        return *this;
    }

    size_t capacity() const { return mCapacity; }

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    /** @brief convert the two dimensional (i,j) to the index to mArray.*/
    size_t linearizedIndex(const size_t i, const size_t j) const
    {
        if (i<j) {
            // Make lower diagonal
            return j*(j-1)/2 + i - 1;
        }
        return i*(i-1)/2 + j - 1;
    }


    double*       mArray;
    const bool    mMemAllocated;
    const size_t  mRows;
    const size_t  mCapacity;

};


/** @class VarMat
 *
 *  @brief matrix of variable dimensions
 *         Memory is allocated either row-major or column major
 *
 *  @remark Note on the internal array.
 *
 *          Row-major
 *          | 0   1   2   3   4|
 *          |                  |
 *          | 5   6   7   8   9|
 *          |                  | 
 *          |10  11  12  13  14|
 *
 *          Column-major
 *          | 0   3   6   9  12|
 *          |                  |
 *          | 1   4   7  10  13|
 *          |                  | 
 *          | 2   5   8  11  14|
 */
class VarMat {

  public:


    /** @brief constructor
     *
     *  @param rows         (in): number of rows of the matrix.
     *
     *  @param cols         (in): number of columns of the matrix.
     *
     *  @param allocateSelf (in): set to true if you want it to  allocate
     *                            internal memory by itself using malloc.
     *  @param padWithZero  (in): set to true if allocateSelf is true, and if
     *                            you want it to clear all the elements with 0.
     *  @param rowMajor     (in): set to true if you want row-major.
     *                            set to false if you want column-major.
     */
    VarMat(
        const size_t rows,
        const size_t cols,
        const bool   allocateSelf,
        const bool   padWithZero,
        const bool   rowMajor
    ):
        mMemAllocated(allocateSelf),
        mRows(rows),
        mCols(cols),
        mRowMajor(rowMajor),
        mCapacity(sizeof(double)*rows*cols)
    {
        if (allocateSelf) {
            mArray = (double*)malloc(mCapacity);
            if (mArray ==nullptr) {
                throw std::bad_alloc();
            }
            if (padWithZero) {
                memset (mArray, (int)0, mCapacity);
            }
        }
    }

    /** @brief descructor */
    ~VarMat() { if (mMemAllocated) { free(mArray); } }


    /** @brief it returns the user the size of the memory required for this
     *         object in bytes.
     */
    size_t requiredMemorySizeInBytes() const { return mCapacity; }
    static size_t requiredMemorySizeInBytes(const size_t r, const size_t c) {
        return sizeof(double) * r * c;
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
        mArray = (double*)p;
        if (padWithZero) {
            memset (mArray, (int)0, mCapacity);
        }
    }


    /** @brief returns the internal array */
    double* array(){ return mArray; }


    /** @brief returns the pointer to the specified element.
     *
     *  @param i   (in): row    (1-origin)
     *
     *  @param j   (in): column (1-origin)
     *
     *  @return pointer to the element
     */
    double* ptr(const size_t i, const size_t j)
    {
        return &(mArray[linearizedIndex(i,j)]);
    }


    /** @brief returns the reference to the specified element.
     *
     *  @param i   (in): row    (1-origin)
     *
     *  @param j   (in): column (1-origin)
     *
     *  @return reference to the element
     */
    double& val(const size_t i, const size_t j)
    {
        return mArray[linearizedIndex(i,j)];
    }


    /** @brief returns the const reference to the specified element.
     *         Used if you want to just read the element.
     *  @param i   (in): row    (1-origin)
     *
     *  @param j   (in): column (1-origin)
     *
     *  @return const reference to the element
     */
    const double& const_val(const size_t i, const size_t j) const
    {
        return mArray[linearizedIndex(i,j)];
    }


    /** @brief returns the number of columns of the matrix.*/
    size_t numCols() const { return mCols; }


    /** @brief returns the number of rows of the matrix.*/
    size_t numRows() const { return mRows; }

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    /** @brief convert the two dimensional (i,j) to the index to mArray.*/
    size_t linearizedIndex(const size_t i, const size_t j) const
    {
        if (mRowMajor) {
            return (i-1)*numCols() + j - 1;
        }
        else {
            return (j-1)*numRows() + i - 1;
        }
    }

    double*       mArray;
    const bool    mMemAllocated;
    const size_t  mRows;
    const size_t  mCols;
    const bool    mRowMajor;
    const size_t  mCapacity;

};


/** @class VarVec
 *
 *  @brief vector of specified length
 */
class VarVec {

  public:


    /** @brief constructor
     *
     *  @param rows         (in): number of elements
     *
     *  @param allocateSelf (in): set to true if you want it to  allocate
     *                            internal memory by itself using malloc.
     *  @param padWithZero  (in): set to true if allocateSelf is true, and if
     *                            you want it to clear all the elements with 0.
     */
    VarVec(const size_t rows, const bool allocateSelf, const bool padWithZero):
        mMemAllocated(allocateSelf),
        mRows(rows),
        mCapacity(sizeof(double)*rows)
    {
        if (allocateSelf) {
            mArray = (double*)malloc(mCapacity);
            if (mArray ==nullptr) {
                throw std::bad_alloc();
            }
            if (padWithZero) {
                memset (mArray, (int)0, mCapacity);
            }
        }
    }


    /** @brief descructor */
    ~VarVec() {if (mMemAllocated) { free(mArray); } }


    /** @brief it returns the user the size of the memory required for this
     *         object in bytes.
     */
    size_t requiredMemorySizeInBytes() const { return mCapacity; }
    static size_t requiredMemorySizeInBytes(const size_t rows) {
        return sizeof(double) * rows;
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
        mArray = (double*)p;
        if (padWithZero) {
            memset (mArray, (int)0, mCapacity);
        }
    }


    /** @brief returns the internal array */
    double* array(){ return mArray; }


    /** @brief fills all the elements with 0.*/
    void zero(){ memset(mArray, (int)0, mCapacity); }


    /** @brief returns the pointer to the specified element.
     *
     *  @param i   (in): element index (1-origin)
     *
     *  @return pointer to the element
     */
    double* ptr(const size_t i) { return &(mArray[i-1]); }


    /** @brief returns the reference to the specified element.
     *
     *  @param i   (in): element index (1-origin)
     *
     *  @return reference to the element
     */
    double& val(const size_t i){ return mArray[i-1]; }


    /** @brief returns the const reference to the specified element.
     *         Used if you want to just read the element.
     *
     *  @param i   (in): element index (1-origin)
     *
     *  @return const reference to the element
     */
    const double& const_val(const size_t i) const { return mArray[i-1]; }


    /** @brief returns the number of elements of the vector.*/
    size_t numRows() const { return mRows; }


    /** @brief copy the contents of the RHS into this. */
    VarVec& operator = (const VarVec& rhs)
    {
        if (mCapacity < rhs.mCapacity) {
            throw std::range_error("VarVec operator=");
        }
        memcpy(mArray, rhs.mArray, rhs.mCapacity);
        return *this;
    }


    /** @brief returns the infinite norm of this vector.*/
    double normInf() const
    {
        double val = 0.0;
        for (int i = 0; i < mRows; i++) {
            val = std::max(val, fabs(mArray[i]));
        }
        return val;
    }

    size_t capacity() const {return mCapacity;}

    void copyInUpperSubvectorAndClearRest(
                                  VarVec& u, long gapStart, long gapWidth)
    {

        if (gapStart > u.mRows + 1) {
            gapStart = u.mRows + 1;
        }
        if (gapStart <  0) {
            gapStart = 0;
        }
        if (gapWidth + u.mRows > mRows) {
            gapWidth = mRows - u.mRows;
        }
        if (gapWidth < 0) {
            gapWidth = 0;
        }

        char*  pBaseSrcBytes = (char*)u.mArray;
        char*  pBaseDstBytes = (char*)mArray;
        const size_t eSize = sizeof(double);

        // Part 1
        long const part1LengthBytes = (gapStart-1) * eSize;
        memcpy( pBaseDstBytes, pBaseSrcBytes, part1LengthBytes);

        // Gap 1
        const long gap1StartBytes  = part1LengthBytes;
        const long gap1LengthBytes = gapWidth*eSize;
        memset( pBaseDstBytes + gap1StartBytes, (int)0, gap1LengthBytes);

        // Part 2
        const long part2StartSrcBytes = gap1StartBytes;
        const long part2StartDstBytes = gap1StartBytes + gap1LengthBytes;
        const long part2LengthBytes   = (u.mRows + 1 - gapStart)*eSize;
        memcpy( pBaseDstBytes + part2StartDstBytes, 
                pBaseSrcBytes + part2StartSrcBytes, 
                part2LengthBytes               );

        // Gap 2
        const long gap2StartBytes  = part2StartDstBytes + part2LengthBytes;
        const long gap2LengthBytes = mRows * eSize - gap2StartBytes;
        memset( pBaseDstBytes + gap2StartBytes, (int)0, gap2LengthBytes);

    }


#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    double*       mArray;
    const bool    mMemAllocated;
    const size_t  mRows;
    const size_t  mCapacity;

};


}// namespace Makena


#endif/*_MAKENA_VARIABLE_PRIMITIVES_HPP_*/




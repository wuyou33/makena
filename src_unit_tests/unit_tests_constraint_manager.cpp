#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include "constraint_manager.hpp"

namespace Makena {


class ConstraintManagerTest : public ::testing::Test {

  protected:  

    ConstraintManagerTest(){;}
    virtual ~ConstraintManagerTest(){;}
    virtual void SetUp() {;};
    virtual void TearDown() {;};

};


/**  @brief ConstraintManager::ColumnManager
 */
TEST_F(ConstraintManagerTest, Test01) {

    ConstraintManager::ColumnManager cm;

    for (long i = 1; i <= 1000; i++) {

        auto res_01 = cm.allocateOne();
        EXPECT_EQ(res_01, i);
        EXPECT_EQ(cm.mDisusedIndices.size(), 0);

    }

    for (long i = 1000; i >= 1; i--) {

        EXPECT_EQ(cm.highestIndex(), i);
        cm.returnOne(i);
        EXPECT_EQ(cm.highestIndex(), i-1);
        EXPECT_EQ(cm.mDisusedIndices.size(), 0);

    }

    for (long i = 1; i <= 1000; i++) {

        auto res_01 = cm.allocateOne();
        EXPECT_EQ(res_01, i);
        EXPECT_EQ(cm.mDisusedIndices.size(), 0);

    }

    for (long i = 999; i >= 1; i--) {

        EXPECT_EQ(cm.highestIndex(), 1000);
        cm.returnOne(i);
        EXPECT_EQ(cm.highestIndex(), 1000);
        EXPECT_EQ(cm.mDisusedIndices.size(), 1000-i);
    }

    for (long i = 1; i <= 500; i++) {

        EXPECT_EQ(cm.highestIndex(), 1000);
        auto res_01 = cm.allocateOne();
        EXPECT_EQ(res_01, i);
        EXPECT_EQ(cm.highestIndex(), 1000);
        EXPECT_EQ(cm.mDisusedIndices.size(), 999-i);
    }

    cm.reset();

    EXPECT_EQ(cm.highestIndex(), 0);
    EXPECT_EQ(cm.mDisusedIndices.size(), 0);

}


/**  @brief ConstraintManager::MemoryManager
 */
TEST_F(ConstraintManagerTest, Test02) {

    ConstraintManager::MemoryManager mm(1000, 100);

    long Mat1, Vec1, WM1;
    for (long i = 0; i < 1000; i++) {

        mm.calcMemorySizes(i, Mat1, Vec1, WM1);
        EXPECT_EQ(Mat1, i*(i+1)/2 * 8);
        EXPECT_EQ(Vec1, i * 8);
        EXPECT_EQ(WM1, i*(i+1)/2 * 8 + i * 8 * 4);
    }

    long cols1 = 100;
    mm.calcMemorySizes(cols1, Mat1, Vec1, WM1);
    mm.allocateMemory(
            cols1, &mm.mM, &mm.mq, &mm.mz, &mm.mz_low, &mm.mz_high, &mm.mWM);
    memset(mm.mM,     (int)0xef, Mat1);
    memset(mm.mq,     (int)0xef, Vec1);
    memset(mm.mz,     (int)0xef, Vec1);
    memset(mm.mz_low, (int)0xef, Vec1);
    memset(mm.mz_high,(int)0xef, Vec1);
    memset(mm.mWM,    (int)0xef, WM1);

    mm.fillWithZero(25, 76, 
             (char*)mm.mM, (char*)mm.mq, 
             (char*)mm.mz, (char*)mm.mz_low, (char*)mm.mz_high
    );


    double* dM      = (double*)mm.mM;
    double* dq      = (double*)mm.mq;
    double* dz      = (double*)mm.mz;
    double* dz_low  = (double*)mm.mz_low;
    double* dz_high = (double*)mm.mz_high;

    long cols1Low = 25;
    long cols1High = 76;
    for (long i = 1; i <= cols1; i++) {

        for (long iM = i*(i-1)/2; iM < i*(i-1)/2 + i; iM++) {
            if (cols1Low <= i && i < cols1High) {
                EXPECT_EQ(dM[iM], 0.0);
            }
            else {
                EXPECT_NE(dM[iM], 0.0);
            }
        }        

        if (cols1Low <= i && i < cols1High) {
            EXPECT_EQ(dq[i-1], 0.0);
            EXPECT_EQ(dz[i-1], 0.0);
            EXPECT_EQ(dz_low[i-1], 0.0);
            EXPECT_EQ(dz_high[i-1], 0.0);
        }
        else {
            EXPECT_NE(dq[i-1], 0.0);
            EXPECT_NE(dz[i-1], 0.0);
            EXPECT_NE(dz_low[i-1], 0.0);
            EXPECT_NE(dz_high[i-1], 0.0);
        }
    }

    long cols2 = 50;
    void *mM, *mq, *mz, *mz_low, *mz_high, *mWM;
    long Mat2, Vec2, WM2;
    mm.calcMemorySizes(cols2, Mat2, Vec2, WM2);
    mm.allocateMemory(cols2, &mM, &mq, &mz, &mz_low, &mz_high, &mWM);
            
    memset(mM,     (int)0xcd, Mat2);
    memset(mq,     (int)0xcd, Vec2);
    memset(mz,     (int)0xcd, Vec2);
    memset(mz_low, (int)0xcd, Vec2);
    memset(mz_high,(int)0xcd, Vec2);
    memset(mWM,    (int)0xcd, WM2);

    double* dM2      = (double*)mM;
    double* dq2      = (double*)mq;
    double* dz2      = (double*)mz;
    double* dz_low2  = (double*)mz_low;
    double* dz_high2 = (double*)mz_high;
    mm.copyInExisting(cols1Low-1, mM, mq, mz, mz_low, mz_high);

    long long vali1 = 0xefefefefefefefef;
    long long vali2 = 0xcdcdcdcdcdcdcdcd; 
    double*  vald1p = (double*)&vali1;
    double*  vald2p = (double*)&vali2;
    for (long i = 1; i <= cols2; i++) {

        for (long iM = i*(i-1)/2; iM < i*(i-1)/2 + i; iM++) {
            if (i <= cols1Low-1) {
                EXPECT_EQ(dM2[iM], *vald1p);
            }
            else {
                EXPECT_EQ(dM2[iM], *vald2p);
            }
        }        

        if (i <= cols1Low-1) {
            EXPECT_EQ(dq2[i-1],      *vald1p);
            EXPECT_EQ(dz2[i-1],      *vald1p);
            EXPECT_EQ(dz_low2[i-1],  *vald1p);
            EXPECT_EQ(dz_high2[i-1], *vald1p);
        }
        else {
            EXPECT_EQ(dq2[i-1],      *vald2p);
            EXPECT_EQ(dz2[i-1],      *vald2p);
            EXPECT_EQ(dz_low2[i-1],  *vald2p);
            EXPECT_EQ(dz_high2[i-1], *vald2p);
        }
    }

    long clearedCol = 60;
    mm.clearMcol(clearedCol);

    for (long i = 1; i <= cols1; i++) {

        for (long iM = i*(i-1)/2; iM < i*(i-1)/2 + i; iM++) {
            if (cols1Low <= i && i < cols1High) {
                EXPECT_EQ(dM[iM], 0.0);
            }
            else if (i==clearedCol) {
                EXPECT_EQ(dM[iM], 0.0);
            }
            else {
                EXPECT_NE(dM[iM], 0.0);
            }
        }        
    }

    mm.freeMemory(&mm.mM, &mm.mq, &mm.mz, &mm.mz_low, &mm.mz_high, &mm.mWM);

}



//
}// Makena

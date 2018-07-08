#include <time.h>
#include "gtest/gtest.h"
#include "mlcp.hpp"

namespace Makena {

 
class MLCPTests : public ::testing::Test {
 
  protected:

    MLCPTests(){;}
    virtual ~MLCPTests(){;}

    virtual void SetUp() {;};
    virtual void TearDown() {;};

};
 

TEST_F(MLCPTests, Test1) {

    MLCP::WorkMemory wm(5, true, true);
    SymMat M(5, true, true);
    VarVec q(5, true, true);
    VarVec z(5, true, true);
    VarVec z_low(5, true, true);
    VarVec z_high(5, true, true);
    vector<enum MLCP::columnType> types;
    types.push_back(MLCP::BILATERAL_FREE);
    types.push_back(MLCP::BILATERAL_FREE);
    types.push_back(MLCP::BILATERAL_FREE);
    types.push_back(MLCP::UNILATERAL_NON_ACTIVE);
    types.push_back(MLCP::UNILATERAL_NON_ACTIVE);

    SymMat a01(3, true, true);
    M.val(1,1) = 2.0;
    M.val(1,2) =-1.0;
    M.val(1,3) = 0.0;
    M.val(2,2) = 2.0;
    M.val(2,3) =-1.0;
    M.val(3,3) = 2.0;
    M.val(4,4) = 2.0;
    M.val(4,5) =-1.0;
    M.val(5,5) = 2.0;

    EXPECT_EQ(M.const_val(1,1),  2.0);
    EXPECT_EQ(M.const_val(1,2), -1.0);
    EXPECT_EQ(M.const_val(1,3),  0.0);
    EXPECT_EQ(M.const_val(1,4),  0.0);
    EXPECT_EQ(M.const_val(1,5),  0.0);
    EXPECT_EQ(M.const_val(2,1), -1.0);
    EXPECT_EQ(M.const_val(2,2),  2.0);
    EXPECT_EQ(M.const_val(2,3), -1.0);
    EXPECT_EQ(M.const_val(2,4),  0.0);
    EXPECT_EQ(M.const_val(2,5),  0.0);
    EXPECT_EQ(M.const_val(3,1),  0.0);
    EXPECT_EQ(M.const_val(3,2), -1.0);
    EXPECT_EQ(M.const_val(3,3),  2.0);
    EXPECT_EQ(M.const_val(3,4),  0.0);
    EXPECT_EQ(M.const_val(3,5),  0.0);
    EXPECT_EQ(M.const_val(4,5), -1.0);
    EXPECT_EQ(M.const_val(4,4),  2.0);
    EXPECT_EQ(M.const_val(4,5), -1.0);
    EXPECT_EQ(M.const_val(5,5),  2.0);

    q.val(1) = -1.5;
    q.val(2) =  3.5;
    q.val(3) =  2.5;
    q.val(4) = -0.5;
    q.val(5) = -3.7;

    EXPECT_EQ(q.const_val(1), -1.5);
    EXPECT_EQ(q.const_val(2),  3.5);
    EXPECT_EQ(q.const_val(3),  2.5);
    EXPECT_EQ(q.const_val(4), -0.5);
    EXPECT_EQ(q.const_val(5), -3.7);

    auto iter = MLCP::solve(
               M, q, z, z_low, z_high, 5, types, 5, 5, 3, 0.0000000000001, wm);

#ifdef PRINT_PROBLEM
    double r01 = M.const_val(1,1) * z.const_val(1) + 
                 M.const_val(1,2) * z.const_val(2) + 
                 M.const_val(1,3) * z.const_val(3) +
                 M.const_val(1,4) * z.const_val(4) +
                 M.const_val(1,5) * z.const_val(5) +

                 m01.mq.const_val(1);

    cerr << "R01: " << r01 << "\n";

    double r02 = M.const_val(2,1) * z.const_val(1) + 
                 M.const_val(2,2) * z.const_val(2) + 
                 M.const_val(2,3) * z.const_val(3) +
                 M.const_val(2,4) * z.const_val(4) +
                 M.const_val(2,5) * z.const_val(5) +
                 q.const_val(2);

    cerr << "R02: " << r02 << "\n";

    double r03 = M.const_val(3,1) * z.const_val(1) + 
                 M.const_val(3,2) * z.const_val(2) + 
                 M.const_val(3,3) * z.const_val(3) +
                 M.const_val(3,4) * z.const_val(4) +
                 M.const_val(3,5) * z.const_val(5) +
                 q.const_val(3);
    cerr << "R03: " << r03 << "\n";

    double r04 = M.const_val(4,1) * z.const_val(1) + 
                 M.const_val(4,2) * z.const_val(2) + 
                 M.const_val(4,3) * z.const_val(3) +
                 M.const_val(4,4) * z.const_val(4) +
                 M.const_val(4,5) * z.const_val(5) +
                 q.const_val(4);
    cerr << "R04: " << r04 << "\n";

    double r05 = M.const_val(5,1) * z.const_val(1) + 
                 M.const_val(5,2) * z.const_val(2) + 
                 M.const_val(5,3) * z.const_val(3) +
                 M.const_val(5,4) * z.const_val(4) +
                 M.const_val(5,5) * z.const_val(5) +
                 q.const_val(5);
    cerr << "R05: " << r05 << "\n";

    cerr << "Z[1]: " << z.const_val(1) << "\n";
    cerr << "Z[2]: " << z.const_val(2) << "\n";
    cerr << "Z[3]: " << z.const_val(3) << "\n";
    cerr << "Z[4]: " << z.const_val(4) << "\n";
    cerr << "Z[5]: " << z.const_val(5) << "\n";
#endif
}


static double frand() {
    double f = ((double)rand())/RAND_MAX;
    return -1.0 + f * 2.0;
}

//#define PRINT_PROBLEM 1

// Test2 deleted.

TEST_F(MLCPTests, Test3) {

    int numTrials  = 100;

    int numA = 100;

    SymMat A(numA, true, true);
    VarVec b(numA, true, true);
    VarVec x(numA, true, true);
    SymMat Acopy(numA, true, true);

    double sum_time = 0.0;
    for (int r = 0; r < numTrials; r++) {
        double* p1 = A.array();
        double* p2 = Acopy.array();
        for (int i = 0; i < numA*(numA+1)/2; i++) {
            p2[i] = p1[i] = frand();
        }

        for (int i = 1; i <= numA; i++) {
            A.val(i,i) = (((double)numA) + A.const_val(i,i));
            Acopy.val(i,i) = (((double)numA) + Acopy.const_val(i,i));
        }

        for (int i = 1; i <= numA; i++) {
            b.val(i) = 5.0 * frand();
        }

        clock_t tStart = clock();
        SLECholesky::solve(A, b, x);
        clock_t tEnd = clock();
        sum_time += (tEnd - tStart);
        double sum = 0.0;
        for (int i = 1; i <= numA; i++) {
            double row = 0.0;
            for (int j = 1; j <= numA; j++) {
                row += (Acopy.val(i,j) * x.val(j));
            }
            row += b.val(i);
            sum += fabs(row);
        }
        EXPECT_LE(sum, 0.00000000001);
    }
}


TEST_F(MLCPTests, Test4) {

    int numTrials  = 200;

    int numAB = 100;

    MLCP::WorkMemory wm(numAB, true, true);
    SymMat M(numAB, true, true);
    VarVec q(numAB, true, true);
    VarVec z(numAB, true, true);
    VarVec z_low(numAB, true, true);
    VarVec z_high(numAB, true, true);

    vector<enum MLCP::columnType> types;
    for (size_t i = 0; i < numAB; i++) {types.push_back(MLCP::NOT_USED);}

    double sum_time = 0.0;
    for (int r = 0; r < numTrials; r++) {

        double* a01 = M.array();

        for (int i = 0; i < numAB*(numAB+1)/2; i++) {
            a01[i] = frand();
        }

        for (int i = 1; i <= numAB; i++) {
            M.val(i,i) = (((double)numAB) + M.const_val(i,i));
        }

        for (int i = 1; i <= numAB; i++) {
            q.val(i) = 5.0 * frand();
        }

        for (int i = 1; i <= numAB; i++) {

            auto j = rand()%10;
            if (j==0||j==1||j==2||j==3) {
                types[i-1] = MLCP::UNILATERAL_NON_ACTIVE;
            }
            else if (j==4||j==5||j==6||j==7) {
                types[i-1] = MLCP::BILATERAL_FREE;
            }
            else if (j==8) { 
                types[i-1] = MLCP::BILATERAL_BOXED_NOT_CLAMPED;
                z_low.val(i)  =  -0.1;
                z_high.val(i) =   0.1;
            }
            else if (j==9) {
                types[i-1] = MLCP::NOT_USED;
            }
            else { cerr << "ERROR\n"; }

//             if (i<=50) {
//                 types[i] = MLCP::BILATERAL_FREE;
//             }
//             else {
//                 types[i] = MLCP::UNILATERAL_NON_ACTIVE;
//             }
        }

#ifdef PRINT_PROBLEM
        cerr << setprecision(2);
        for (int i = 1; i <= numAB; i++) {

            for (int j = 1; j <= numAB; j++) {
                cerr << M.val(i,j) << "\t";
            }
            cerr << "|\t" << q.val(i) << "\n";
        }
#endif
        z.zero();

        clock_t tStart = clock();
        auto iter = MLCP::solve(
           M, q, z, z_low, z_high, numAB, types, 5, 5, 3, 0.0000000000001, wm);
        clock_t tEnd = clock();
#ifdef PRINT_PROBLEM
        cerr << "Solved in " << iter << "\n";
        cerr << "Time: " << (tEnd - tStart) << "[micro sec]\n";
#endif
        sum_time += (tEnd - tStart);

//        cerr << "Bilateral solutions\n";

        double sum = 0.0;
        for (int i = 1; i <= numAB; i++) {
            auto& ti = types[i-1];
            if (ti!=MLCP::NOT_USED) {
                double row = 0.0;
                for (int j = 1; j <= numAB; j++) {
                    auto& tj = types[j-1];
                    if (tj!=MLCP::NOT_USED) {
                        row += (M.val(i,j) * z.val(j));
                    }
                }
                row += q.val(i);

#ifdef PRINT_PROBLEM
                cerr << "z[" << i << "]:";
                switch(types[i-1]) {

                  case MLCP::NOT_USED:
                    cerr << "NOT_USED";
                    break;
                  case MLCP::BILATERAL_FREE:
                    cerr << "BILATERAL_FREE";
                    break;

                  case MLCP::BILATERAL_BOXED_NOT_CLAMPED:
                    cerr << "BILATERAL_BOXED_NOT_CLAMPED";
                    break;

                  case MLCP::BILATERAL_BOXED_CLAMPED_LOW:
                    cerr << "BILATERAL_BOXED_CLAMPED_LOW";
                    break;

                  case MLCP::BILATERAL_BOXED_CLAMPED_HIGH:
                    cerr << "BILATERAL_BOXED_CLAMPED_HIGH";
                    break;

                  case MLCP::UNILATERAL_NON_ACTIVE:
                    cerr << "UNILATERAL_NON_ACTIVE";
                    break;

                  case MLCP::UNILATERAL_ACTIVE:
                  default:
                    cerr << "UNILATERAL_ACTIVE";
                    ;
                }
                cerr << "  " << z.val(i) << "\n";

                cerr << "w[" << i << "]:" << row << "\n";
#endif

                if( ti==MLCP::BILATERAL_FREE               ||
                    ti==MLCP::BILATERAL_BOXED_NOT_CLAMPED  ||
                    ti==MLCP::BILATERAL_BOXED_CLAMPED_LOW  ||
                    ti==MLCP::BILATERAL_BOXED_CLAMPED_HIGH    ) {

                    sum += fabs(row);
                }
                else {
                    EXPECT_GE(row , -0.00001); 
                    EXPECT_GE(z.val(i) , -0.00001);
                    sum += fabs(row * z.val(i));
                }
            }
        }

        EXPECT_LE(sum, 0.0000001);
#ifdef PRINT_PROBLEM
        cerr << "\n";
        cerr << "Error:" << sum << "\n";
#endif
    }
//#ifdef PRINT_PROBLEM
    cerr << "Avg time: " << sum_time / numTrials << "[micro sec]\n";
//#endif
}


} // namespace Makena

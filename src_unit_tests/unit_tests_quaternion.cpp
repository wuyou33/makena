#include "gtest/gtest.h"
#include "quaternion.hpp"

namespace Makena {

class QuaternionTests : public ::testing::Test {
 
  protected:

    QuaternionTests(){;}
    virtual ~QuaternionTests(){;}

    virtual void SetUp() {;};
    virtual void TearDown() {;};

    inline double mS(Quaternion& q){ return q.mS; }
    inline Vec3& mV(Quaternion& q){ return q.mV; }

};


TEST_F(QuaternionTests, Test1) {
    Quaternion q1;

    EXPECT_EQ(mS(q1), 0.0);
    EXPECT_EQ(mV(q1).x(), 0.0);
    EXPECT_EQ(mV(q1).y(), 0.0);
    EXPECT_EQ(mV(q1).z(), 0.0);
}


TEST_F(QuaternionTests, Test2) {
    Quaternion q1(0.1, 0.2, 0.3, 0.4);

    EXPECT_EQ(mS(q1), 0.1);
    EXPECT_EQ(mV(q1).x(), 0.2);
    EXPECT_EQ(mV(q1).y(), 0.3);
    EXPECT_EQ(mV(q1).z(), 0.4);
}


TEST_F(QuaternionTests, Test3) {
    Vec3 v1(0.2, 0.3, 0.4);
    Quaternion q1(0.1, v1);

    EXPECT_EQ(mS(q1), 0.1);
    EXPECT_EQ(mV(q1).x(), 0.2);
    EXPECT_EQ(mV(q1).y(), 0.3);
    EXPECT_EQ(mV(q1).z(), 0.4);
}


TEST_F(QuaternionTests, Test4) {
    Vec3 v1(1.0, 0.0, 0.0);
    Quaternion q1(v1, 0.5*M_PI);

//cerr << "mS:" << mS(q1) << "\n";
//cerr << "mV:" << mV(q1)[1] << "," << mV(q1)[2] << "," << mV(q1)[3] << "\n";

    EXPECT_EQ(fabs(mS(q1) - 1.0/sqrt(2.0))< EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).x() - 1.0/sqrt(2.0))< EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).y() - 0.0)< EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).z() - 0.0)< EPSILON_LINEAR, true);

    Vec3 v2(1.0, 1.0, 1.0);
    Quaternion q2(v2, 2.0*M_PI/3.0);

    EXPECT_EQ(fabs(mS(q2) - 0.5)< EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).x() - 0.5)< EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).y() - 0.5)< EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).z() - 0.5)< EPSILON_LINEAR, true);

}


TEST_F(QuaternionTests, Test5) {

    Mat3x3 m1( 1.0, 2.0, 3.0,
               4.0, 5.0, 6.0,
               7.0, 8.0, 9.0 );
    Quaternion q1(m1);
    // s = 0.5*sqrt(1+r11+r22+r33)
    // x = (r32 - r23)/4s
    // y = (r13 - r31)/4s
    // z = (r21 - r12)/4s

    EXPECT_EQ(fabs(mS(q1)     - 2.0 ) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).x() - 0.25) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).y() - (-0.5)) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).z() - 0.25)< EPSILON_LINEAR, true);

//  cerr << "mS:" << mS(q1) << "\n";
//  cerr << "mV:" << mV(q1)[1] << "," << mV(q1)[2] << "," << mV(q1)[3] << "\n";

    Mat3x3 m2( 1.0, 2.0, 3.0,
               4.0,-5.0, 6.0,
               7.0, 8.0,-9.0 );
    Quaternion q2(m2);

    // s = (r32-r23)/4x            = 0.25
    // x = 0.5*sqrt(1+r11-r22-r33) = 2.0
    // y = (r12+r21)/4x            = 0.75
    // z = (r13+r31)/4x            = 1.25

    EXPECT_EQ(fabs(mS(q2)     - 0.25) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).x() - 2.0 ) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).y() - 0.75) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).z() - 1.25) < EPSILON_LINEAR, true);

//  cerr << "mS:" << mS(q2) << "\n";
//  cerr << "mV:" << mV(q2)[1] << "," << mV(q2)[2] << "," << mV(q2)[3] << "\n";

    Mat3x3 m3(-1.0, 2.0, 3.0,
               4.0, 5.0, 6.0,
               7.0, 8.0,-9.0 );
    Quaternion q3(m3);

    // s = (r13-r31)/4y             = 0.5
    // x = (r12+r21)/4y             = 0.75
    // y = 0.5*sqrt(1-r11+r22-r33)  = 2.0
    // z = (r23+r32)/4y             = 3.5

    EXPECT_EQ(fabs(mS(q3)     + 0.50) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q3).x() - 0.75) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q3).y() - 2.0 ) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q3).z() - 1.75) < EPSILON_LINEAR, true);

    // s = (r21-r12)/4z             = 0.25
    // x = (r13+r31)/4z             = 1.25
    // y = (r23+r32)/4z             = 1.75
    // z = 0.5*sqrt(1-r11-r22+r33)  = 2.0

    Mat3x3 m4(-1.0, 2.0, 3.0,
               4.0,-5.0, 6.0,
               7.0, 8.0,+9.0 );
    Quaternion q4(m4);

    EXPECT_EQ(fabs(mS(q4)     - 0.25) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q4).x() - 1.25) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q4).y() - 1.75) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q4).z() - 2.00) < EPSILON_LINEAR, true);
}


TEST_F(QuaternionTests, Test6) {

    Quaternion q1(0.1, 0.2, 0.3, 0.4);

    Quaternion q2(0.5, 0.6, 0.7, 0.8);

    q2 += q1;

    EXPECT_EQ(fabs(mS(q1)     - 0.1) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).x() - 0.2) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).y() - 0.3) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).z() - 0.4) < EPSILON_LINEAR, true);

    EXPECT_EQ(fabs(mS(q2)     - 0.6) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).x() - 0.8) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).y() - 1.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).z() - 1.2) < EPSILON_LINEAR, true);
    
}


TEST_F(QuaternionTests, Test7) {

    Quaternion q1(0.1, 0.2, 0.3, 0.4);

    Quaternion q2(0.5, 0.6, 0.7, 0.8);

    q2 -= q1;

    EXPECT_EQ(fabs(mS(q1)     - 0.1) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).x() - 0.2) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).y() - 0.3) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).z() - 0.4) < EPSILON_LINEAR, true);

    EXPECT_EQ(fabs(mS(q2)     - 0.4) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).x() - 0.4) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).y() - 0.4) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).z() - 0.4) < EPSILON_LINEAR, true);
    
}


TEST_F(QuaternionTests, Test8) {

    Quaternion q1(0.1, 0.2, 0.3, 0.4);

    Quaternion q2(0.5, 0.6, 0.7, 0.8);

    Quaternion q3(0.1, 0.2, 0.3, 0.4);
    Quaternion q4(0.1+1.1*EPSILON_LINEAR, 0.2, 0.3, 0.4);
    Quaternion q5(0.1+0.9*EPSILON_LINEAR, 0.2, 0.3, 0.4);
    Quaternion q6(0.1, 0.2+1.1*EPSILON_LINEAR, 0.3, 0.4);
    Quaternion q7(0.1, 0.2+0.9*EPSILON_LINEAR, 0.3, 0.4);
    Quaternion q8(0.1, 0.2, 0.3+1.1*EPSILON_LINEAR, 0.4);
    Quaternion q9(0.1, 0.2, 0.3+0.9*EPSILON_LINEAR, 0.4);
    Quaternion q10(0.1, 0.2, 0.3, 0.4+1.1*EPSILON_LINEAR);
    Quaternion q11(0.1, 0.2, 0.3, 0.4+0.9*EPSILON_LINEAR);


    EXPECT_EQ(q1 == q2, false);
    EXPECT_EQ(q1 != q2, true);
    EXPECT_EQ(q1 == q3, true);
    EXPECT_EQ(q1 != q3, false);
    EXPECT_EQ(q1 == q4, false);
    EXPECT_EQ(q1 != q4, true);
    EXPECT_EQ(q1 == q5, true);
    EXPECT_EQ(q1 != q5, false);
    EXPECT_EQ(q1 == q6, false);
    EXPECT_EQ(q1 != q6, true);
    EXPECT_EQ(q1 == q7, true);
    EXPECT_EQ(q1 != q7, false);
    EXPECT_EQ(q1 == q8, false);
    EXPECT_EQ(q1 != q8, true);
    EXPECT_EQ(q1 == q9, true);
    EXPECT_EQ(q1 != q9, false);
    EXPECT_EQ(q1 == q10, false);
    EXPECT_EQ(q1 != q10, true);
    EXPECT_EQ(q1 == q11, true);
    EXPECT_EQ(q1 != q11, false);

}


TEST_F(QuaternionTests, Test9) {

    Quaternion q1(0.1, 0.2, 0.3, 0.4);

    Quaternion q2(0.5, 0.6, 0.7, 0.8);

    Quaternion q3 = q2 + q1;

    EXPECT_EQ(fabs(mS(q1)     - 0.1) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).x() - 0.2) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).y() - 0.3) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).z() - 0.4) < EPSILON_LINEAR, true);

    EXPECT_EQ(fabs(mS(q2)     - 0.5) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).x() - 0.6) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).y() - 0.7) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).z() - 0.8) < EPSILON_LINEAR, true);

    EXPECT_EQ(fabs(mS(q3)     - 0.6) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q3).x() - 0.8) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q3).y() - 1.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q3).z() - 1.2) < EPSILON_LINEAR, true);
}


TEST_F(QuaternionTests, Test10) {

    Quaternion q1(0.1, 0.2, 0.3, 0.4);

    Quaternion q2(0.5, 0.6, 0.7, 0.8);

    Quaternion q3 = q1 - q2;

    EXPECT_EQ(fabs(mS(q1)     - 0.1) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).x() - 0.2) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).y() - 0.3) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).z() - 0.4) < EPSILON_LINEAR, true);

    EXPECT_EQ(fabs(mS(q2)     - 0.5) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).x() - 0.6) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).y() - 0.7) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).z() - 0.8) < EPSILON_LINEAR, true);

    EXPECT_EQ(fabs(mS(q3)     + 0.4) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q3).x() + 0.4) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q3).y() + 0.4) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q3).z() + 0.4) < EPSILON_LINEAR, true);
}


TEST_F(QuaternionTests, Test11) {

    Quaternion q1(0.1, 0.2, 0.3, 0.4);

    Quaternion q2(0.5, 0.6, 0.7, 0.8);

    Quaternion q3 = q1 * q2;

    Vec3 v1(0.2, 0.3, 0.4);
    Vec3 v2(0.6, 0.7, 0.8);

    EXPECT_EQ(fabs(mS(q1)     - 0.1) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).x() - 0.2) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).y() - 0.3) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).z() - 0.4) < EPSILON_LINEAR, true);

    EXPECT_EQ(fabs(mS(q2)     - 0.5) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).x() - 0.6) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).y() - 0.7) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).z() - 0.8) < EPSILON_LINEAR, true);

    EXPECT_EQ(fabs(mS(q3)     + 0.6)  < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q3).x() - 0.12) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q3).y() - 0.30) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q3).z() - 0.24) < EPSILON_LINEAR, true);

}


TEST_F(QuaternionTests, Test12) {

    Quaternion q1(0.1, 0.2, 0.3, 0.4);

    q1.scale(0.5);

    EXPECT_EQ(fabs(mS(q1)     - 0.05) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).x() - 0.1) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).y() - 0.15) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).z() - 0.2) < EPSILON_LINEAR, true);

}


TEST_F(QuaternionTests, Test13) {

    Quaternion q1(0.1, 0.2, 0.3, 0.4);

    Quaternion q2 = q1.conjugate();

    EXPECT_EQ(fabs(mS(q1)     - 0.1) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).x() - 0.2) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).y() - 0.3) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q1).z() - 0.4) < EPSILON_LINEAR, true);

    EXPECT_EQ(fabs(mS(q2)     - 0.1) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).x() + 0.2) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).y() + 0.3) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).z() + 0.4) < EPSILON_LINEAR, true);

}


TEST_F(QuaternionTests, Test14) {
/*
 *    +---------------------------------------------------+
 *    |s^2+x^2-y^2-z^2      2*(xy-sz)         2*(sy+xz)   |
 *    |   2(xy+sz)       s^2-x^2+y^2-z^2      2*(yz-sx)   |
 *    |   2(xz-sy)          2*(sx+yz)      s^2-x^2-y^2+z^2|
 *    +---------------------------------------------------+
 */

    Quaternion q1(1.0, 2.0, 3.0, 4.0);

    Mat3x3 m1 = q1.rotationMatrix();

    EXPECT_EQ(fabs(m1.cell(1,1)  + 20.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(m1.cell(1,2)  -  4.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(m1.cell(1,3)  - 22.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(m1.cell(2,1)  - 20.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(m1.cell(2,2)  + 10.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(m1.cell(2,3)  - 20.0) < EPSILON_LINEAR, true);

}


TEST_F(QuaternionTests, Test15) {

/*
 *    +------------------------+
 *    |-0.5*x   -0.5*y   -0.5*z|
 *    | 0.5*s    0.5*z   -0.5*y|
 *    |-0.5*z    0.5*s    0.5*x|
 *    | 0.5*y   -0.5*x    0.5*s|
 *    +------------------------+
 */

    Quaternion q1(1.0, 2.0, 3.0, 4.0);

    double v1[12];
    q1.matrix4x3(v1);

    EXPECT_EQ(v1[ 0], -1.0);
    EXPECT_EQ(v1[ 1], -1.5);
    EXPECT_EQ(v1[ 2], -2.0);
    EXPECT_EQ(v1[ 3],  0.5);
    EXPECT_EQ(v1[ 4],  2.0);
    EXPECT_EQ(v1[ 5], -1.5);
    EXPECT_EQ(v1[ 6], -2.0);
    EXPECT_EQ(v1[ 7],  0.5);
    EXPECT_EQ(v1[ 8],  1.0);
    EXPECT_EQ(v1[ 9],  1.5);
    EXPECT_EQ(v1[10], -1.0);   
    EXPECT_EQ(v1[11],  0.5);

}


TEST_F(QuaternionTests, Test16) {

/*
 *    +------------------------+
 *    |-0.5*x   -0.5*y   -0.5*z|
 *    | 0.5*s    0.5*z   -0.5*y|
 *    |-0.5*z    0.5*s    0.5*x|
 *    | 0.5*y   -0.5*x    0.5*s|
 *    +------------------------+
 */

    Quaternion q1(1.0, 2.0, 3.0, 4.0);
    Vec3 w(0.7, 0.5, 0.3);

    Quaternion q2 = q1.derivative(w);

    // q2 = 0.5*[0,w]*q1 = 0.5*[0.0, 0.7, 0.5, 0.3]*[1.0, 2.0, 3.0, 4.0]
    // s = 0.5 * (0.0 - (1.4 + 1.5 + 1.2)) = -2.05
    // x = 0.5 * (0.7 + 1.1) = 0.9
    // y = 0.5 * (0.5 - 2.2) = 0.85
    // z = 0.5 * (0.3 + 1.1) = 0.7

    //   i    j    k
    // [0.7, 0.5, 0.3]
    // [2.0, 3.0, 4.0]

    EXPECT_EQ(fabs(mS(q2)     + 2.05) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).x() - 0.9 ) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).y() + 0.85) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(q2).z() - 0.7 ) < EPSILON_LINEAR, true);

}


TEST_F(QuaternionTests, Test17) {

    Vec3 v1(1.0, 0.0, 0.0);
    Quaternion q1(v1, 0.5*M_PI);
    Vec3 p1(1.0, 1.0, 1.0);

    Vec3 p2 = q1.rotate(p1);

    EXPECT_EQ(fabs(p2.x() - 1.0 ) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(p2.y() + 1.0 ) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(p2.z() - 1.0 ) < EPSILON_LINEAR, true);

    Vec3 v2(0.0, 1.0, 0.0);
    Quaternion q2(v2, 0.5*M_PI);

    Vec3 p3 = q2.rotate(p1);

    EXPECT_EQ(fabs(p3.x() - 1.0 ) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(p3.y() - 1.0 ) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(p3.z() + 1.0 ) < EPSILON_LINEAR, true);

    Vec3 v3(0.0, 0.0, 1.0);
    Quaternion q3(v3, 0.5*M_PI);

    Vec3 p4 = q3.rotate(p1);

    EXPECT_EQ(fabs(p4.x() + 1.0 ) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(p4.y() - 1.0 ) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(p4.z() - 1.0 ) < EPSILON_LINEAR, true);
}


TEST_F(QuaternionTests, Test18) {

    Quaternion q1(1.0, 2.0, 3.0, 4.0);

    EXPECT_EQ(fabs(q1.s()  - 1.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(q1.i()  - 2.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(q1.j()  - 3.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(q1.k()  - 4.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(q1.x()  - 2.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(q1.y()  - 3.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(q1.z()  - 4.0) < EPSILON_LINEAR, true);

}


TEST_F(QuaternionTests, Test19) {

    Quaternion q1(1.0, 2.0, 3.0, 4.0);
    q1.normalize();
    //    1.0+4.0+9.0+16+0 = 1.0/sqrt(30.0)

    EXPECT_EQ(fabs(q1.s()  - 1.0/sqrt(30.0)) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(q1.i()  - 2.0/sqrt(30.0)) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(q1.j()  - 3.0/sqrt(30.0)) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(q1.k()  - 4.0/sqrt(30.0)) < EPSILON_LINEAR, true);

}

} // namespace Makena

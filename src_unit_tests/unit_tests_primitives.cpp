#include "gtest/gtest.h"
#include "quaternion.hpp"
#include "primitives.hpp"

namespace Makena {

 
class PrimitivesTests : public ::testing::Test {
 
  protected:

    PrimitivesTests(){;}
    virtual ~PrimitivesTests(){;}

    virtual void SetUp() {;};
    virtual void TearDown() {;};

    std::array<double,3>& mV(Vec3& v) {return v.mV;}
    std::array<double,9>& mV(Mat3x3& v) {return v.mV;}
};


static double rand100()
{
    return ((rand()%10000)+1)/100.0;
}


static Vec3 randVec3D100()
{
    return Vec3(rand100(), rand100(), rand100());
}


static Quaternion randomRotQ()
{
    Vec3 v(0.0, 0.0, 0.0);
    while (v == Vec3(0.0, 0.0, 0.0)) {
        v = Vec3( (rand()%100)/100.0 - 0.5, 
                  (rand()%100)/100.0 - 0.5, 
                  (rand()%100)/100.0 - 0.5 );
    }
    
    Quaternion q(v, (rand()%100)/50.0 * M_PI);
    q.normalize();
    return q;
}


TEST_F(PrimitivesTests, Test1) {
    Vec3 v1;
    EXPECT_EQ(mV(v1)[0],0);
    EXPECT_EQ(mV(v1)[1],0);
    EXPECT_EQ(mV(v1)[2],0);
}


TEST_F(PrimitivesTests, Test2) {
    double x = 0.123;
    double y = 0.456;
    double z = 0.789;

    Vec3 v1(x, y, z);
    EXPECT_EQ(mV(v1)[0],0.123);
    EXPECT_EQ(mV(v1)[1],0.456);
    EXPECT_EQ(mV(v1)[2],0.789);
}

TEST_F(PrimitivesTests, Test3) {

    std::array<double,3> a1{{0.111,0.222,0.333}};

    Vec3 v1(a1);
    EXPECT_EQ(mV(v1)[0],0.111);
    EXPECT_EQ(mV(v1)[1],0.222);
    EXPECT_EQ(mV(v1)[2],0.333);
}

TEST_F(PrimitivesTests, Test4) {

    double a1[3] = {0.111,0.222,0.333};

    Vec3 v1(a1);
    EXPECT_EQ(mV(v1)[0],0.111);
    EXPECT_EQ(mV(v1)[1],0.222);
    EXPECT_EQ(mV(v1)[2],0.333);
}


TEST_F(PrimitivesTests, Test5) {

    Vec3 v1(1.0, 2.0, 3.0);
    Vec3 v2(5.0, 7.0, 9.0);

    v2 += v1;

    EXPECT_EQ(mV(v2)[0], 6.0);
    EXPECT_EQ(mV(v2)[1], 9.0);
    EXPECT_EQ(mV(v2)[2], 12.0);
    EXPECT_EQ(mV(v1)[0], 1.0);
    EXPECT_EQ(mV(v1)[1], 2.0);
    EXPECT_EQ(mV(v1)[2], 3.0);


}


TEST_F(PrimitivesTests, Test6) {

    Vec3 v1(1.0, 2.0, 3.0);
    Vec3 v2(5.0, 7.0, 9.0);

    v2 -= v1;

    EXPECT_EQ(mV(v2)[0], 4.0);
    EXPECT_EQ(mV(v2)[1], 5.0);
    EXPECT_EQ(mV(v2)[2], 6.0);

    EXPECT_EQ(mV(v1)[0], 1.0);
    EXPECT_EQ(mV(v1)[1], 2.0);
    EXPECT_EQ(mV(v1)[2], 3.0);

}


TEST_F(PrimitivesTests, Test7) {

    Vec3 v1(1.0, 2.0, 3.0);
    Vec3 v2(1.01, 2.0, 3.0);
    Vec3 v3(1.00, 2.01, 3.0);
    Vec3 v4(1.00, 2.0, 3.01);
    Vec3 v5(1.0, 2.0, 3.0);
    Vec3 v6(1.0, 2.0, 3.0 + 0.9*EPSILON_LINEAR);
    Vec3 v7(1.0, 2.0, 3.0 + 1.01*EPSILON_LINEAR);
    Vec3 v8(1.0, 2.0 + 0.9*EPSILON_LINEAR, 3.0);
    Vec3 v9(1.0 + 0.9*EPSILON_LINEAR, 2.0, 3.0);


    EXPECT_EQ(v1==v2, false);
    EXPECT_EQ(v1!=v2, true);
    EXPECT_EQ(v1==v3, false);
    EXPECT_EQ(v1!=v3, true);
    EXPECT_EQ(v1==v4, false);
    EXPECT_EQ(v1!=v4, true);
    EXPECT_EQ(v1==v5, true);
    EXPECT_EQ(v1!=v5, false);
    EXPECT_EQ(v1==v6, true);
    EXPECT_EQ(v1!=v6, false);
    EXPECT_EQ(v1==v7, false);
    EXPECT_EQ(v1!=v7, true);
    EXPECT_EQ(v1==v8, true);
    EXPECT_EQ(v1!=v8, false);
    EXPECT_EQ(v1==v9, true);
    EXPECT_EQ(v1!=v9, false);

}


TEST_F(PrimitivesTests, Test8) {

    Vec3 v1(1.0, 2.0, 3.0);
    Vec3 v2(4.0, 5.0, 6.0);


    Vec3 v3 = v1 + v2;

    EXPECT_EQ(mV(v1)[0], 1.0);
    EXPECT_EQ(mV(v1)[1], 2.0);
    EXPECT_EQ(mV(v1)[2], 3.0);
    EXPECT_EQ(mV(v2)[0], 4.0);
    EXPECT_EQ(mV(v2)[1], 5.0);
    EXPECT_EQ(mV(v2)[2], 6.0);
    EXPECT_EQ(mV(v3)[0], 5.0);
    EXPECT_EQ(mV(v3)[1], 7.0);
    EXPECT_EQ(mV(v3)[2], 9.0);

}


TEST_F(PrimitivesTests, Test9) {

    Vec3 v1(1.0, 2.0, 3.0);
    Vec3 v2(4.0, 5.0, 6.0);


    Vec3 v3 = v1 - v2;

    EXPECT_EQ(mV(v1)[0], 1.0);
    EXPECT_EQ(mV(v1)[1], 2.0);
    EXPECT_EQ(mV(v1)[2], 3.0);
    EXPECT_EQ(mV(v2)[0], 4.0);
    EXPECT_EQ(mV(v2)[1], 5.0);
    EXPECT_EQ(mV(v2)[2], 6.0);
    EXPECT_EQ(mV(v3)[0], -3.0);
    EXPECT_EQ(mV(v3)[1], -3.0);
    EXPECT_EQ(mV(v3)[2], -3.0);

}

TEST_F(PrimitivesTests, Test10) {

    Vec3 v1(1.0, 2.0, 3.0);
    Vec3 v2 = v1 * 0.5;

    EXPECT_EQ(mV(v1)[0], 1.0);
    EXPECT_EQ(mV(v1)[1], 2.0);
    EXPECT_EQ(mV(v1)[2], 3.0);

    EXPECT_EQ(mV(v2)[0], 0.5);
    EXPECT_EQ(mV(v2)[1], 1.0);
    EXPECT_EQ(mV(v2)[2], 1.5);
}


TEST_F(PrimitivesTests, Test11) {

    Vec3 v1(1.0, 2.0, 3.0);

    EXPECT_EQ(v1[1], 1.0);
    EXPECT_EQ(v1[2], 2.0);
    EXPECT_EQ(v1[3], 3.0);

}


TEST_F(PrimitivesTests, Test12) {

    Vec3 v1(1.0, 2.0, 3.0);
    Vec3 v2(4.0, 5.0, 6.0);

    double d1 = v1.dot(v2);

    EXPECT_EQ(d1, 32);
}


TEST_F(PrimitivesTests, Test13) {

    Vec3 v1(2.0, 3.0, 4.0);
    Vec3 v2(5.0, 6.0, 7.0);

    Vec3 v3 = v1.cross(v2);

    EXPECT_EQ(v3[1], -3.0);
    EXPECT_EQ(v3[2],  6.0);
    EXPECT_EQ(v3[3], -3.0);
}


TEST_F(PrimitivesTests, Test14) {

    Vec3 v1(2.0, 3.0, 4.0);

    Vec3 v2 = v1.perp();

    EXPECT_EQ( fabs(v1.dot(v2)) < EPSILON_LINEAR, true);

    EXPECT_EQ(v2[1],  0.0);
    EXPECT_EQ(v2[2],  4.0);
    EXPECT_EQ(v2[3], -3.0);

}


TEST_F(PrimitivesTests, Test15) {

    Vec3 v1(2.0, 3.0, 4.0);

    Mat3x3 m1 = v1.crossMat();

    EXPECT_EQ(m1.cell(1,1),  0.0);
    EXPECT_EQ(m1.cell(1,2), -4.0);
    EXPECT_EQ(m1.cell(1,3),  3.0);
    EXPECT_EQ(m1.cell(2,1),  4.0);
    EXPECT_EQ(m1.cell(2,2),  0.0);
    EXPECT_EQ(m1.cell(2,3), -2.0);
    EXPECT_EQ(m1.cell(3,1), -3.0);
    EXPECT_EQ(m1.cell(3,2),  2.0);
    EXPECT_EQ(m1.cell(3,3),  0.0);
}


TEST_F(PrimitivesTests, Test16) {

    Vec3 v1(3.0, 4.0, 5.0);

    double d1 = v1.norm2();

    EXPECT_EQ(fabs(d1 - 7.071067811865475)<EPSILON_LINEAR, true);
}


TEST_F(PrimitivesTests, Test17) {

    Vec3 v1(3.0, 4.0, 5.0);

    double d1 = v1.squaredNorm2();

    EXPECT_EQ(fabs(d1 - 50.0)<EPSILON_LINEAR, true);
}


TEST_F(PrimitivesTests, Test18) {

    Vec3 v1(3.0, 4.0, 5.0);

    v1.scale(0.5);

    EXPECT_EQ(v1[1],  1.5);
    EXPECT_EQ(v1[2],  2.0);
    EXPECT_EQ(v1[3],  2.5);
}


TEST_F(PrimitivesTests, Test19) {

    Vec3 v1(3.0, 4.0, 5.0);

    v1.set(1.23, 2.34, 3.45);

    EXPECT_EQ(v1[1],  1.23);
    EXPECT_EQ(v1[2],  2.34);
    EXPECT_EQ(v1[3],  3.45);
}


TEST_F(PrimitivesTests, Test20) {

    Vec3 v1(3.0, 4.0, 5.0);
    double a[3] = {1.23, 2.34, 3.45};
    v1.set(a);

    EXPECT_EQ(v1[1],  1.23);
    EXPECT_EQ(v1[2],  2.34);
    EXPECT_EQ(v1[3],  3.45);
}


TEST_F(PrimitivesTests, Test21) {

    Vec3 v1(3.0, 4.0, 5.0);

    v1.setX(1.23);

    EXPECT_EQ(v1[1],  1.23);
    EXPECT_EQ(v1[2],  4.0);
    EXPECT_EQ(v1[3],  5.0);
}


TEST_F(PrimitivesTests, Test22) {

    Vec3 v1(3.0, 4.0, 5.0);

    v1.setY(1.23);

    EXPECT_EQ(v1[1],  3.0);
    EXPECT_EQ(v1[2],  1.23);
    EXPECT_EQ(v1[3],  5.0);
}


TEST_F(PrimitivesTests, Test23) {

    Vec3 v1(3.0, 4.0, 5.0);

    v1.setZ(1.23);

    EXPECT_EQ(v1[1],  3.0);
    EXPECT_EQ(v1[2],  4.0);
    EXPECT_EQ(v1[3],  1.23);
}


TEST_F(PrimitivesTests, Test24) {

    Vec3 v1(3.0, 4.0, 5.0);

    v1.zero();

    EXPECT_EQ(v1[1],  0.0);
    EXPECT_EQ(v1[2],  0.0);
    EXPECT_EQ(v1[3],  0.0);
}


TEST_F(PrimitivesTests, Test25) {

    Vec3 v1(3.0, 4.0, 5.0);

    v1.normalize();

    EXPECT_EQ(fabs(v1[1]-3.0/7.071067811865475)< EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(v1[2]-4.0/7.071067811865475)< EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(v1[3]-5.0/7.071067811865475)< EPSILON_LINEAR, true);

}


TEST_F(PrimitivesTests, Test26) {

    Vec3 v1(3.0, 4.0, 5.0);

    EXPECT_EQ(v1.x(), 3.0);
    EXPECT_EQ(v1.y(), 4.0);
    EXPECT_EQ(v1.z(), 5.0);


}


TEST_F(PrimitivesTests, Test27) {

    Mat3x3 m1;

    EXPECT_EQ(mV(m1)[0],0);
    EXPECT_EQ(mV(m1)[1],0);
    EXPECT_EQ(mV(m1)[2],0);
    EXPECT_EQ(mV(m1)[3],0);
    EXPECT_EQ(mV(m1)[4],0);
    EXPECT_EQ(mV(m1)[5],0);
    EXPECT_EQ(mV(m1)[6],0);
    EXPECT_EQ(mV(m1)[7],0);
    EXPECT_EQ(mV(m1)[8],0);

}


TEST_F(PrimitivesTests, Test28) {

    double v[9] ={1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0};
    Mat3x3 m1(v);
    EXPECT_EQ(mV(m1)[0],1.0);
    EXPECT_EQ(mV(m1)[1],2.0);
    EXPECT_EQ(mV(m1)[2],3.0);
    EXPECT_EQ(mV(m1)[3],4.0);
    EXPECT_EQ(mV(m1)[4],5.0);
    EXPECT_EQ(mV(m1)[5],6.0);
    EXPECT_EQ(mV(m1)[6],7.0);
    EXPECT_EQ(mV(m1)[7],8.0);
    EXPECT_EQ(mV(m1)[8],9.0);

}


TEST_F(PrimitivesTests, Test29) {

    Mat3x3 m1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    EXPECT_EQ(mV(m1)[0],1.0);
    EXPECT_EQ(mV(m1)[1],2.0);
    EXPECT_EQ(mV(m1)[2],3.0);
    EXPECT_EQ(mV(m1)[3],4.0);
    EXPECT_EQ(mV(m1)[4],5.0);
    EXPECT_EQ(mV(m1)[5],6.0);
    EXPECT_EQ(mV(m1)[6],7.0);
    EXPECT_EQ(mV(m1)[7],8.0);
    EXPECT_EQ(mV(m1)[8],9.0);

}


TEST_F(PrimitivesTests, Test30) {

    Vec3 v1(1.0, 2.0, 3.0);
    Vec3 v2(4.0, 5.0, 6.0);
    Vec3 v3(7.0, 8.0, 9.0);

    Mat3x3 m1(v1, v2, v3);
    EXPECT_EQ(mV(m1)[0],1.0);
    EXPECT_EQ(mV(m1)[1],4.0);
    EXPECT_EQ(mV(m1)[2],7.0);
    EXPECT_EQ(mV(m1)[3],2.0);
    EXPECT_EQ(mV(m1)[4],5.0);
    EXPECT_EQ(mV(m1)[5],8.0);
    EXPECT_EQ(mV(m1)[6],3.0);
    EXPECT_EQ(mV(m1)[7],6.0);
    EXPECT_EQ(mV(m1)[8],9.0);

}


TEST_F(PrimitivesTests, Test31) {


    Mat3x3 m1(1.0,   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3 m2(13.0, 21.0, 23.0, 27.0, 31.0, 33.0, 37.0, 43.0, 47.0);

    m1 += m2;
    EXPECT_EQ(mV(m1)[0],14.0);
    EXPECT_EQ(mV(m1)[1],23.0);
    EXPECT_EQ(mV(m1)[2],26.0);
    EXPECT_EQ(mV(m1)[3],31.0);
    EXPECT_EQ(mV(m1)[4],36.0);
    EXPECT_EQ(mV(m1)[5],39.0);
    EXPECT_EQ(mV(m1)[6],44.0);
    EXPECT_EQ(mV(m1)[7],51.0);
    EXPECT_EQ(mV(m1)[8],56.0);

}


TEST_F(PrimitivesTests, Test32) {


    Mat3x3 m1(1.0,   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3 m2(13.0, 21.0, 23.0, 27.0, 31.0, 33.0, 37.0, 43.0, 47.0);

    m2 -= m1;
    EXPECT_EQ(mV(m2)[0],12.0);
    EXPECT_EQ(mV(m2)[1],19.0);
    EXPECT_EQ(mV(m2)[2],20.0);
    EXPECT_EQ(mV(m2)[3],23.0);
    EXPECT_EQ(mV(m2)[4],26.0);
    EXPECT_EQ(mV(m2)[5],27.0);
    EXPECT_EQ(mV(m2)[6],30.0);
    EXPECT_EQ(mV(m2)[7],35.0);
    EXPECT_EQ(mV(m2)[8],38.0);

}


TEST_F(PrimitivesTests, Test33) {


    Mat3x3  m1(1.0,   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3  m2(13.0, 21.0, 23.0, 27.0, 31.0, 33.0, 37.0, 43.0, 47.0);
    Mat3x3  m3(1.0 + 1.1*EPSILON_LINEAR,   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3  m4(1.0 + 0.9*EPSILON_LINEAR,   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3  m5(1.0,  2.0 + 1.1*EPSILON_LINEAR, 3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3  m6(1.0,  2.0 + 0.9*EPSILON_LINEAR, 3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3  m7(1.0,  2.0,  3.0 + 1.1*EPSILON_LINEAR, 4.0,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3  m8(1.0,  2.0,  3.0 + 0.9*EPSILON_LINEAR, 4.0,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3  m9(1.0,  2.0,  3.0, 4.0 + 1.1*EPSILON_LINEAR,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3 m10(1.0,  2.0,  3.0, 4.0 + 0.9*EPSILON_LINEAR,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3 m11(1.0,  2.0,  3.0, 4.0,  5.0 + 1.1*EPSILON_LINEAR,  6.0,  7.0,  8.0,  9.0);
    Mat3x3 m12(1.0,  2.0,  3.0, 4.0,  5.0 + 0.9*EPSILON_LINEAR,  6.0,  7.0,  8.0,  9.0);
    Mat3x3 m13(1.0,  2.0,  3.0, 4.0,  5.0,  6.0 + 1.1*EPSILON_LINEAR,  7.0,  8.0,  9.0);
    Mat3x3 m14(1.0,  2.0,  3.0, 4.0,  5.0,  6.0 + 0.9*EPSILON_LINEAR,  7.0,  8.0,  9.0);
    Mat3x3 m15(1.0,  2.0,  3.0, 4.0,  5.0,  6.0,  7.0 + 1.1*EPSILON_LINEAR,  8.0,  9.0);
    Mat3x3 m16(1.0,  2.0,  3.0, 4.0,  5.0,  6.0,  7.0 + 0.9*EPSILON_LINEAR,  8.0,  9.0);
    Mat3x3 m17(1.0,  2.0,  3.0, 4.0,  5.0,  6.0,  7.0,  8.0 + 1.1*EPSILON_LINEAR,  9.0);
    Mat3x3 m18(1.0,  2.0,  3.0, 4.0,  5.0,  6.0,  7.0,  8.0 + 0.9*EPSILON_LINEAR,  9.0);
    Mat3x3 m19(1.0,  2.0,  3.0, 4.0,  5.0,  6.0,  7.0,  8.0,  9.0 + 1.1*EPSILON_LINEAR);
    Mat3x3 m20(1.0,  2.0,  3.0, 4.0,  5.0,  6.0,  7.0,  8.0,  9.0 + 0.9*EPSILON_LINEAR);
    Mat3x3 m21(1.0,   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0);

    EXPECT_EQ(m1 == m2, false);
    EXPECT_EQ(m1 != m2, true);
    EXPECT_EQ(m1 == m3, false);
    EXPECT_EQ(m1 != m3, true);
    EXPECT_EQ(m1 == m4, true);
    EXPECT_EQ(m1 != m4, false);
    EXPECT_EQ(m1 == m5, false);
    EXPECT_EQ(m1 != m5, true);
    EXPECT_EQ(m1 == m6, true);
    EXPECT_EQ(m1 != m6, false);
    EXPECT_EQ(m1 == m7, false);
    EXPECT_EQ(m1 != m7, true);
    EXPECT_EQ(m1 == m8, true);
    EXPECT_EQ(m1 != m8, false);
    EXPECT_EQ(m1 == m9, false);
    EXPECT_EQ(m1 != m9, true);
    EXPECT_EQ(m1 == m10, true);
    EXPECT_EQ(m1 != m10, false);
    EXPECT_EQ(m1 == m11, false);
    EXPECT_EQ(m1 != m11, true);
    EXPECT_EQ(m1 == m12, true);
    EXPECT_EQ(m1 != m12, false);
    EXPECT_EQ(m1 == m13, false);
    EXPECT_EQ(m1 != m13, true);
    EXPECT_EQ(m1 == m14, true);
    EXPECT_EQ(m1 != m14, false);
    EXPECT_EQ(m1 == m15, false);
    EXPECT_EQ(m1 != m15, true);
    EXPECT_EQ(m1 == m16, true);
    EXPECT_EQ(m1 != m16, false);
    EXPECT_EQ(m1 == m17, false);
    EXPECT_EQ(m1 != m17, true);
    EXPECT_EQ(m1 == m18, true);
    EXPECT_EQ(m1 != m18, false);
    EXPECT_EQ(m1 == m19, false);
    EXPECT_EQ(m1 != m19, true);
    EXPECT_EQ(m1 == m20, true);
    EXPECT_EQ(m1 != m20, false);
    EXPECT_EQ(m1 == m21, true);
    EXPECT_EQ(m1 != m21, false);
    EXPECT_EQ(m1 == m1, true);
    EXPECT_EQ(m1 != m1, false);

}


TEST_F(PrimitivesTests, Test34) {


    Mat3x3 m1(1.0,   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3 m2(13.0, 21.0, 23.0, 27.0, 31.0, 33.0, 37.0, 43.0, 47.0);

    Mat3x3 m3 = m2 + m1;

    EXPECT_EQ(mV(m1)[0],1.0);
    EXPECT_EQ(mV(m1)[1],2.0);
    EXPECT_EQ(mV(m1)[2],3.0);
    EXPECT_EQ(mV(m1)[3],4.0);
    EXPECT_EQ(mV(m1)[4],5.0);
    EXPECT_EQ(mV(m1)[5],6.0);
    EXPECT_EQ(mV(m1)[6],7.0);
    EXPECT_EQ(mV(m1)[7],8.0);
    EXPECT_EQ(mV(m1)[8],9.0);

    EXPECT_EQ(mV(m2)[0],13.0);
    EXPECT_EQ(mV(m2)[1],21.0);
    EXPECT_EQ(mV(m2)[2],23.0);
    EXPECT_EQ(mV(m2)[3],27.0);
    EXPECT_EQ(mV(m2)[4],31.0);
    EXPECT_EQ(mV(m2)[5],33.0);
    EXPECT_EQ(mV(m2)[6],37.0);
    EXPECT_EQ(mV(m2)[7],43.0);
    EXPECT_EQ(mV(m2)[8],47.0);

    EXPECT_EQ(mV(m3)[0],14.0);
    EXPECT_EQ(mV(m3)[1],23.0);
    EXPECT_EQ(mV(m3)[2],26.0);
    EXPECT_EQ(mV(m3)[3],31.0);
    EXPECT_EQ(mV(m3)[4],36.0);
    EXPECT_EQ(mV(m3)[5],39.0);
    EXPECT_EQ(mV(m3)[6],44.0);
    EXPECT_EQ(mV(m3)[7],51.0);
    EXPECT_EQ(mV(m3)[8],56.0);

}


TEST_F(PrimitivesTests, Test35) {


    Mat3x3 m1(1.0,   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0);
    Mat3x3 m2(13.0, 21.0, 23.0, 27.0, 31.0, 33.0, 37.0, 43.0, 47.0);

    Mat3x3 m3 = m2 - m1;

    EXPECT_EQ(mV(m1)[0],1.0);
    EXPECT_EQ(mV(m1)[1],2.0);
    EXPECT_EQ(mV(m1)[2],3.0);
    EXPECT_EQ(mV(m1)[3],4.0);
    EXPECT_EQ(mV(m1)[4],5.0);
    EXPECT_EQ(mV(m1)[5],6.0);
    EXPECT_EQ(mV(m1)[6],7.0);
    EXPECT_EQ(mV(m1)[7],8.0);
    EXPECT_EQ(mV(m1)[8],9.0);

    EXPECT_EQ(mV(m2)[0],13.0);
    EXPECT_EQ(mV(m2)[1],21.0);
    EXPECT_EQ(mV(m2)[2],23.0);
    EXPECT_EQ(mV(m2)[3],27.0);
    EXPECT_EQ(mV(m2)[4],31.0);
    EXPECT_EQ(mV(m2)[5],33.0);
    EXPECT_EQ(mV(m2)[6],37.0);
    EXPECT_EQ(mV(m2)[7],43.0);
    EXPECT_EQ(mV(m2)[8],47.0);

    EXPECT_EQ(mV(m3)[0],12.0);
    EXPECT_EQ(mV(m3)[1],19.0);
    EXPECT_EQ(mV(m3)[2],20.0);
    EXPECT_EQ(mV(m3)[3],23.0);
    EXPECT_EQ(mV(m3)[4],26.0);
    EXPECT_EQ(mV(m3)[5],27.0);
    EXPECT_EQ(mV(m3)[6],30.0);
    EXPECT_EQ(mV(m3)[7],35.0);
    EXPECT_EQ(mV(m3)[8],38.0);

}


TEST_F(PrimitivesTests, Test36) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );

    Mat3x3 m2 = m1.transpose();

    EXPECT_EQ(mV(m1)[0],1.0);
    EXPECT_EQ(mV(m1)[1],2.0);
    EXPECT_EQ(mV(m1)[2],3.0);
    EXPECT_EQ(mV(m1)[3],4.0);
    EXPECT_EQ(mV(m1)[4],5.0);
    EXPECT_EQ(mV(m1)[5],6.0);
    EXPECT_EQ(mV(m1)[6],7.0);
    EXPECT_EQ(mV(m1)[7],8.0);
    EXPECT_EQ(mV(m1)[8],9.0);

    EXPECT_EQ(mV(m2)[0],1.0);
    EXPECT_EQ(mV(m2)[1],4.0);
    EXPECT_EQ(mV(m2)[2],7.0);
    EXPECT_EQ(mV(m2)[3],2.0);
    EXPECT_EQ(mV(m2)[4],5.0);
    EXPECT_EQ(mV(m2)[5],8.0);
    EXPECT_EQ(mV(m2)[6],3.0);
    EXPECT_EQ(mV(m2)[7],6.0);
    EXPECT_EQ(mV(m2)[8],9.0);

}


TEST_F(PrimitivesTests, Test37) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );

    m1.transposeInPlace();

    EXPECT_EQ(mV(m1)[0],1.0);
    EXPECT_EQ(mV(m1)[1],4.0);
    EXPECT_EQ(mV(m1)[2],7.0);
    EXPECT_EQ(mV(m1)[3],2.0);
    EXPECT_EQ(mV(m1)[4],5.0);
    EXPECT_EQ(mV(m1)[5],8.0);
    EXPECT_EQ(mV(m1)[6],3.0);
    EXPECT_EQ(mV(m1)[7],6.0);
    EXPECT_EQ(mV(m1)[8],9.0);

}


TEST_F(PrimitivesTests, Test38) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );

    Mat3x3 m2(10.0,  11.0, 12.0,
              13.0,  14.0, 15.0, 
              16.0,  17.0, 18.0 );


    Mat3x3 m3 = m1 * m2;

    EXPECT_EQ(mV(m3)[0],84.0);
    EXPECT_EQ(mV(m3)[1],90.0);
    EXPECT_EQ(mV(m3)[2],96.0);
    EXPECT_EQ(mV(m3)[3],201.0);
    EXPECT_EQ(mV(m3)[4],216.0);
    EXPECT_EQ(mV(m3)[5],231.0);
    EXPECT_EQ(mV(m3)[6],318.0);
    EXPECT_EQ(mV(m3)[7],342.0);
    EXPECT_EQ(mV(m3)[8],366.0);

}


TEST_F(PrimitivesTests, Test39) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );

    Mat3x3 m2 = m1.pow2();

    EXPECT_EQ(mV(m2)[0],30.0);
    EXPECT_EQ(mV(m2)[1],36.0);
    EXPECT_EQ(mV(m2)[2],42.0);
    EXPECT_EQ(mV(m2)[3],66.0);
    EXPECT_EQ(mV(m2)[4],81.0);
    EXPECT_EQ(mV(m2)[5],96.0);
    EXPECT_EQ(mV(m2)[6],102.0);
    EXPECT_EQ(mV(m2)[7],126.0);
    EXPECT_EQ(mV(m2)[8],150.0);

}


TEST_F(PrimitivesTests, Test40) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );

    Vec3 v1(10.0, 11.0, 12.0);

    Vec3 v2 = m1 * v1;

    EXPECT_EQ(mV(v2)[0],68.0);
    EXPECT_EQ(mV(v2)[1],167.0);
    EXPECT_EQ(mV(v2)[2],266.0);

}


TEST_F(PrimitivesTests, Test41) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0, 10.0 );

    double d1  = m1.det();

    EXPECT_EQ(d1, -3);

}


TEST_F(PrimitivesTests, Test42) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0, 10.0 );

    Mat3x3 m2 = m1.inverse();

    EXPECT_EQ(fabs(mV(m2)[0] + 2.0/3.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(m2)[1] + 4.0/3.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(m2)[2] - 1.0)     < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(m2)[3] + 2.0/3.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(m2)[4] -11.0/3.0) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(m2)[5] + 2.0)     < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(m2)[6] - 1.0)     < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(m2)[7] + 2.0)     < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(mV(m2)[8] - 1.0)     < EPSILON_LINEAR, true);

    Mat3x3 m3(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );

    EXPECT_THROW(m3.inverse(), std::underflow_error);
}


TEST_F(PrimitivesTests, Test43) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );

    m1.scale(0.5);

    EXPECT_EQ(mV(m1)[0],0.5);
    EXPECT_EQ(mV(m1)[1],1.0);
    EXPECT_EQ(mV(m1)[2],1.5);
    EXPECT_EQ(mV(m1)[3],2.0);
    EXPECT_EQ(mV(m1)[4],2.5);
    EXPECT_EQ(mV(m1)[5],3.0);
    EXPECT_EQ(mV(m1)[6],3.5);
    EXPECT_EQ(mV(m1)[7],4.0);
    EXPECT_EQ(mV(m1)[8],4.5);

}


TEST_F(PrimitivesTests, Test44) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );


    EXPECT_EQ(m1.cell(1,1),1.0);
    EXPECT_EQ(m1.cell(1,2),2.0);
    EXPECT_EQ(m1.cell(1,3),3.0);
    EXPECT_EQ(m1.cell(2,1),4.0);
    EXPECT_EQ(m1.cell(2,2),5.0);
    EXPECT_EQ(m1.cell(2,3),6.0);
    EXPECT_EQ(m1.cell(3,1),7.0);
    EXPECT_EQ(m1.cell(3,2),8.0);
    EXPECT_EQ(m1.cell(3,3),9.0);
}


TEST_F(PrimitivesTests, Test45) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );

    double* dp1 = m1.get_array();
    EXPECT_EQ(dp1[0],1.0);
    EXPECT_EQ(dp1[1],2.0);
    EXPECT_EQ(dp1[2],3.0);
    EXPECT_EQ(dp1[3],4.0);
    EXPECT_EQ(dp1[4],5.0);
    EXPECT_EQ(dp1[5],6.0);
    EXPECT_EQ(dp1[6],7.0);
    EXPECT_EQ(dp1[7],8.0);
    EXPECT_EQ(dp1[8],9.0);
}


TEST_F(PrimitivesTests, Test46) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );

    double d1 = m1.trace();
    EXPECT_EQ(d1, 15.0);

}


TEST_F(PrimitivesTests, Test47) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );

    double d1 = m1.trace();
    EXPECT_EQ(d1, 15.0);

}


TEST_F(PrimitivesTests, Test48) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );

    Vec3 c1 = m1.col(1);
    EXPECT_EQ(c1[1], 1.0);
    EXPECT_EQ(c1[2], 4.0);
    EXPECT_EQ(c1[3], 7.0);

    Vec3 c2 = m1.col(2);
    EXPECT_EQ(c2[1], 2.0);
    EXPECT_EQ(c2[2], 5.0);
    EXPECT_EQ(c2[3], 8.0);

    Vec3 c3 = m1.col(3);
    EXPECT_EQ(c3[1], 3.0);
    EXPECT_EQ(c3[2], 6.0);
    EXPECT_EQ(c3[3], 9.0);
}


TEST_F(PrimitivesTests, Test49) {


    Mat3x3 m1(1.0,  2.0,  3.0,
              4.0,  5.0,  6.0, 
              7.0,  8.0,  9.0 );

    Vec3 r1 = m1.row(1);
    EXPECT_EQ(r1[1], 1.0);
    EXPECT_EQ(r1[2], 2.0);
    EXPECT_EQ(r1[3], 3.0);

    Vec3 r2 = m1.row(2);
    EXPECT_EQ(r2[1], 4.0);
    EXPECT_EQ(r2[2], 5.0);
    EXPECT_EQ(r2[3], 6.0);

    Vec3 r3 = m1.row(3);
    EXPECT_EQ(r3[1], 7.0);
    EXPECT_EQ(r3[2], 8.0);
    EXPECT_EQ(r3[3], 9.0);

}


TEST_F(PrimitivesTests, Test50) {

    Mat3x3 m1( 3.0,  2.0,   4.0,
               2.0,  0.0,   2.0,
               4.0,  2.0,   5.0 );
    Vec3   v1;
    Mat3x3 m2 = m1.EigenVectorsIfSymmetric(v1);

//    cerr << "M1:" << m1.cell(1,1) << "," << m1.cell(1,2) << "," << m1.cell(1,3) << "\n";
//    cerr << "M1:" << m1.cell(2,1) << "," << m1.cell(2,2) << "," << m1.cell(2,3) << "\n";
//    cerr << "M1:" << m1.cell(3,1) << "," << m1.cell(3,2) << "," << m1.cell(3,3) << "\n";
//    cerr << "\n";
//    cerr << "M2:" << m2.cell(1,1) << "," << m2.cell(1,2) << "," << m2.cell(1,3) << "\n";
//    cerr << "M2:" << m2.cell(2,1) << "," << m2.cell(2,2) << "," << m2.cell(2,3) << "\n";
//    cerr << "M2:" << m2.cell(3,1) << "," << m2.cell(3,2) << "," << m2.cell(3,3) << "\n";

//    cerr << "lambda1:" << v1[1] << "\n";
//    cerr << "lambda2:" << v1[2] << "\n";
//    cerr << "lambda3:" << v1[3] << "\n";

    EXPECT_EQ(fabs(v1[1]-9.0)<EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(v1[2]-0.0)<EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(v1[3]+1.0)<EPSILON_LINEAR, true);

    Vec3 ev1 = m2.col(1);
    Vec3 ev2 = m2.col(2);
    Vec3 ev3 = m2.col(3);

//    cerr << "EV1:" << ev1[1] << "," << ev1[2] << "," << ev1[3] << "\n";
//    cerr << "EV2:" << ev2[1] << "," << ev2[2] << "," << ev2[3] << "\n";
//    cerr << "EV3:" << ev3[1] << "," << ev3[2] << "," << ev3[3] << "\n";

    Vec3 Mev1 = m1 * ev1;
    Mev1.scale(v1[1]);
    Mev1.normalize();
    Vec3 Mev2 = m1 * ev2;
    Mev2.scale(v1[2]);
    Mev2.normalize();
    Vec3 Mev3 = m1 * ev3;
    Mev3.scale(v1[3]);
    Mev3.normalize();

//    cerr << "M*EV1:" << Mev1[1] << "," << Mev1[2] << "," << Mev1[3] << "\n";
//    cerr << "M*EV2:" << Mev2[1] << "," << Mev2[2] << "," << Mev2[3] << "\n";
//    cerr << "M*EV3:" << Mev3[1] << "," << Mev3[2] << "," << Mev3[3] << "\n";

    EXPECT_EQ(ev1==Mev1, true);
    if (fabs(v1[2]) < EPSILON_LINEAR) {
        Vec3 zero;
        EXPECT_EQ(zero==Mev2, true);
    }
    else {
        EXPECT_EQ(ev2==Mev2, true);
    }
    EXPECT_EQ(ev3==Mev3, true);

}


TEST_F(PrimitivesTests, Test51) {

    Mat3x3 m1( 1.0,  3.0,   5.0,
               3.0, 11.0,  13.0,
               5.0, 13.0,  17.0 );
    Vec3   v1;
    Mat3x3 m2 = m1.EigenVectorsIfSymmetric(v1);

//    cerr << "M1:" << m1.cell(1,1) << "," << m1.cell(1,2) << "," << m1.cell(1,3) << "\n";
//    cerr << "M1:" << m1.cell(2,1) << "," << m1.cell(2,2) << "," << m1.cell(2,3) << "\n";
//    cerr << "M1:" << m1.cell(3,1) << "," << m1.cell(3,2) << "," << m1.cell(3,3) << "\n";
//    cerr << "\n";
//    cerr << "M2:" << m2.cell(1,1) << "," << m2.cell(1,2) << "," << m2.cell(1,3) << "\n";
//    cerr << "M2:" << m2.cell(2,1) << "," << m2.cell(2,2) << "," << m2.cell(2,3) << "\n";
//    cerr << "M2:" << m2.cell(3,1) << "," << m2.cell(3,2) << "," << m2.cell(3,3) << "\n";

//    cerr << "lambda1:" << v1[1] << "\n";
//    cerr << "lambda2:" << v1[2] << "\n";
//    cerr << "lambda3:" << v1[3] << "\n";

    EXPECT_EQ(fabs(v1[1]-28.5552) <EPSILON_LINEAR*10000.0, true);
    EXPECT_EQ(fabs(v1[2]-1.08832) <EPSILON_LINEAR*10000.0, true);
    EXPECT_EQ(fabs(v1[3]+0.643556)<EPSILON_LINEAR*10000.0, true);

    Vec3 ev1 = m2.col(1);
    Vec3 ev2 = m2.col(2);
    Vec3 ev3 = m2.col(3);

//    cerr << "EV1:" << ev1[1] << "," << ev1[2] << "," << ev1[3] << "\n";
//    cerr << "EV2:" << ev2[1] << "," << ev2[2] << "," << ev2[3] << "\n";
//    cerr << "EV3:" << ev3[1] << "," << ev3[2] << "," << ev3[3] << "\n";

    Vec3 Mev1 = m1 * ev1;
    Mev1.scale(v1[1]);
    Mev1.normalize();
    Vec3 Mev2 = m1 * ev2;
    Mev2.scale(v1[2]);
    Mev2.normalize();
    Vec3 Mev3 = m1 * ev3;
    Mev3.scale(v1[3]);
    Mev3.normalize();

//    cerr << "M*EV1:" << Mev1[1] << "," << Mev1[2] << "," << Mev1[3] << "\n";
//    cerr << "M*EV2:" << Mev2[1] << "," << Mev2[2] << "," << Mev2[3] << "\n";
//    cerr << "M*EV3:" << Mev3[1] << "," << Mev3[2] << "," << Mev3[3] << "\n";

    EXPECT_EQ(ev1==Mev1, true);
    if (fabs(v1[2]) < EPSILON_LINEAR) {
        Vec3 zero;
        EXPECT_EQ(zero==Mev2, true);
    }
    else {
        EXPECT_EQ(ev2==Mev2, true);
    }
    EXPECT_EQ(ev3==Mev3, true);

}



TEST_F(PrimitivesTests, Test52) {

    Mat3x3 m1( 5.0,  4.0,   3.0,
               4.0,  5.0,   3.0,
               3.0,  3.0,   2.0 );
    Vec3   v1;
    Mat3x3 m2 = m1.EigenVectorsIfSymmetric(v1);

//    cerr << "M1:" << m1.cell(1,1) << "," << m1.cell(1,2) << "," << m1.cell(1,3) << "\n";
//    cerr << "M1:" << m1.cell(2,1) << "," << m1.cell(2,2) << "," << m1.cell(2,3) << "\n";
//    cerr << "M1:" << m1.cell(3,1) << "," << m1.cell(3,2) << "," << m1.cell(3,3) << "\n";
//    cerr << "\n";
//    cerr << "M2:" << m2.cell(1,1) << "," << m2.cell(1,2) << "," << m2.cell(1,3) << "\n";
//    cerr << "M2:" << m2.cell(2,1) << "," << m2.cell(2,2) << "," << m2.cell(2,3) << "\n";
//    cerr << "M2:" << m2.cell(3,1) << "," << m2.cell(3,2) << "," << m2.cell(3,3) << "\n";

//    cerr << "lambda1:" << v1[1] << "\n";
//    cerr << "lambda2:" << v1[2] << "\n";
//    cerr << "lambda3:" << v1[3] << "\n";

//    EXPECT_EQ(fabs(v1[1]-28.5552) <EPSILON_LINEAR*10000.0, true);
//    EXPECT_EQ(fabs(v1[2]-1.08832) <EPSILON_LINEAR*10000.0, true);
//    EXPECT_EQ(fabs(v1[3]+0.643556)<EPSILON_LINEAR*10000.0, true);

    Vec3 ev1 = m2.col(1);
    Vec3 ev2 = m2.col(2);
    Vec3 ev3 = m2.col(3);

//    cerr << "EV1:" << ev1[1] << "," << ev1[2] << "," << ev1[3] << "\n";
//    cerr << "EV2:" << ev2[1] << "," << ev2[2] << "," << ev2[3] << "\n";
//    cerr << "EV3:" << ev3[1] << "," << ev3[2] << "," << ev3[3] << "\n";

    Vec3 Mev1 = m1 * ev1;
    Mev1.scale(v1[1]);
    Mev1.normalize();
    Vec3 Mev2 = m1 * ev2;
    Mev2.scale(v1[2]);
    Mev2.normalize();
    Vec3 Mev3 = m1 * ev3;
    Mev3.scale(v1[3]);
    Mev3.normalize();

//    cerr << "M*EV1:" << Mev1[1] << "," << Mev1[2] << "," << Mev1[3] << "\n";
//    cerr << "M*EV2:" << Mev2[1] << "," << Mev2[2] << "," << Mev2[3] << "\n";
//    cerr << "M*EV3:" << Mev3[1] << "," << Mev3[2] << "," << Mev3[3] << "\n";

    EXPECT_EQ(ev1==Mev1, true);
    if (fabs(v1[2]) < EPSILON_LINEAR) {
        Vec3 zero;
        EXPECT_EQ(zero==Mev2, true);
    }
    else {
        EXPECT_EQ(ev2==Mev2, true);
    }

    if (fabs(v1[3]) < EPSILON_LINEAR) {
        Vec3 zero;
        EXPECT_EQ(zero==Mev3, true);
    }
    else {
        EXPECT_EQ(ev3==Mev3, true);
    }

}


/*
 * Following python script is used to find the expected values for 
 * findPrincipalComponents(points, spread)
 */

//import numpy as np
//
//points = [[7.0, 4.0, 3.0],
//          [4.0, 1.0, 8.0],
//          [6.0, 3.0, 5.0],
//          [8.0, 6.0, 1.0],
//          [8.0, 5.0, 7.0],
//          [7.0, 2.0, 9.0],
//          [5.0, 3.0, 3.0],
//          [9.0, 5.0, 8.0],
//          [7.0, 4.0, 5.0],
//          [8.0, 2.0, 2.0] ]
//
//np_points = np.array(points)
//tps  = np.reshape(np_points.T, (3,10))
//
//cov_xy = np.cov(tps[0], tps[1])
//cov_xz = np.cov(tps[0], tps[2])
//cov_yz = np.cov(tps[1], tps[2])
//
//cov = [ [ cov_xy[0][0], cov_xy[0][1], cov_xz[0][1] ],
//        [ cov_xy[1][0], cov_xy[1][1], cov_yz[0][1] ],
//        [ cov_xz[1][0], cov_yz[0][1], cov_yz[1][1] ] ]
//
//np_cov = np.array(cov)
//evals, evec = np.linalg.eigh(np_cov)
//
//print 'Covariance'
//print np_cov
//print 'Eigen values: ' + str(evals[2]) + ' ' + str(evals[1]) + ' ' +  str(evals[0])
//print 'Eigen vector 1: ' + str(evec[:,2])
//print 'Eigen vector 2: ' + str(evec[:,1])
//print 'Eigen vector 3: ' + str(evec[:,0])


TEST_F(PrimitivesTests, Test53) {

    Vec3   v1(7.0, 4.0, 3.0);
    Vec3   v2(4.0, 1.0, 8.0);
    Vec3   v3(6.0, 3.0, 5.0);
    Vec3   v4(8.0, 6.0, 1.0);
    Vec3   v5(8.0, 5.0, 7.0);
    Vec3   v6(7.0, 2.0, 9.0);
    Vec3   v7(5.0, 3.0, 3.0);
    Vec3   v8(9.0, 5.0, 8.0);
    Vec3   v9(7.0, 4.0, 5.0);
    Vec3  v10(8.0, 2.0, 2.0);

    vector<Vec3> points;
    points.push_back(v1);
    points.push_back(v2);
    points.push_back(v3);
    points.push_back(v4);
    points.push_back(v5);
    points.push_back(v6);
    points.push_back(v7);
    points.push_back(v8);
    points.push_back(v9);
    points.push_back(v10);

    Vec3 spread;
    Vec3 mean;
    Mat3x3 pca = findPrincipalComponents(points, spread, mean);

//    cerr << "lambda1:" << spread[1] << "\n";
//    cerr << "lambda2:" << spread[2] << "\n";
//    cerr << "lambda3:" << spread[3] << "\n";

    Vec3 ev1 = pca.col(1);
    Vec3 ev2 = pca.col(2);
    Vec3 ev3 = pca.col(3);

//    cerr << "EV1:" << ev1[1] << "," << ev1[2] << "," << ev1[3] << "\n";
//    cerr << "EV2:" << ev2[1] << "," << ev2[2] << "," << ev2[3] << "\n";
//    cerr << "EV3:" << ev3[1] << "," << ev3[2] << "," << ev3[3] << "\n";

// Covariance
//[[ 2.32222222  1.61111111 -0.43333333]
// [ 1.61111111  2.5        -1.27777778]
// [-0.43333333 -1.27777778  7.87777778]]
//Eigen values: 8.27394258041 3.6761292668 0.749928152795
//Eigen vector 1: [-0.1375708  -0.25045969  0.95830278]
//Eigen vector 2: [ 0.69903712  0.66088917  0.27307986]
//Eigen vector 3: [ 0.70172743 -0.70745703 -0.08416157]

    EXPECT_EQ(fabs(spread[1]-8.27394258041  ) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(spread[2]-3.6761292668   ) < EPSILON_LINEAR, true);
    EXPECT_EQ(fabs(spread[3]-0.749928152795 ) < EPSILON_LINEAR, true);

    EXPECT_EQ( fabs(ev1[1] - (-0.1375708)  ) < EPSILON_LINEAR, true);
    EXPECT_EQ( fabs(ev1[2] - (-0.25045969) ) < EPSILON_LINEAR, true);
    EXPECT_EQ( fabs(ev1[3] - 0.95830278    ) < EPSILON_LINEAR, true);

    EXPECT_EQ( fabs(ev2[1] - 0.69903712    ) < EPSILON_LINEAR, true);
    EXPECT_EQ( fabs(ev2[2] - 0.66088917    ) < EPSILON_LINEAR, true);
    EXPECT_EQ( fabs(ev2[3] - 0.27307986    ) < EPSILON_LINEAR, true);

    EXPECT_EQ( fabs(ev3[1] - (-0.70172743) ) < EPSILON_LINEAR, true);
    EXPECT_EQ( fabs(ev3[2] - ( 0.70745703) ) < EPSILON_LINEAR, true);
    EXPECT_EQ( fabs(ev3[3] - ( 0.08416157) ) < EPSILON_LINEAR, true);

}


TEST_F(PrimitivesTests, Test54) {

    for (long i = 0; i < 10; i++) {

        auto q01 = randomRotQ();

        vector<Vec3> points;
        for (long j = 0; j < 1000; j++) {
            Vec3 p(rand100(), rand100(), 0.0);
            p = q01.rotate(p);
            points.push_back(p);
        }

        Vec3 spread;
        Vec3 mean;
        Mat3x3 pca = findPrincipalComponents(points, spread, mean);

        cerr << "\n";
        cerr << "lambda1:" << spread[1] << "\n";
        cerr << "lambda2:" << spread[2] << "\n";
        cerr << "lambda3:" << spread[3] << "\n";

        Vec3 ev1 = pca.col(1);
        Vec3 ev2 = pca.col(2);
        Vec3 ev3 = pca.col(3);

        cerr << "EV1:" << ev1[1] << "," << ev1[2] << "," << ev1[3] << "\n";
        cerr << "EV2:" << ev2[1] << "," << ev2[2] << "," << ev2[3] << "\n";
        cerr << "EV3:" << ev3[1] << "," << ev3[2] << "," << ev3[3] << "\n";
        cerr << "\n";
    }
}


TEST_F(PrimitivesTests, Test55) {

    for (long i = 0; i < 10; i++) {

        auto q01 = randomRotQ();

        vector<Vec3> points;
        for (long j = 0; j < 1000; j++) {
            Vec3 p(rand100(), 0.0, 0.0);
            p = q01.rotate(p);
            points.push_back(p);
        }

        Vec3 spread;
        Vec3 mean;
        Mat3x3 pca = findPrincipalComponents(points, spread, mean);

        cerr << "\n";
        cerr << "lambda1:" << spread[1] << "\n";
        cerr << "lambda2:" << spread[2] << "\n";
        cerr << "lambda3:" << spread[3] << "\n";

        Vec3 ev1 = pca.col(1);
        Vec3 ev2 = pca.col(2);
        Vec3 ev3 = pca.col(3);

        cerr << "EV1:" << ev1[1] << "," << ev1[2] << "," << ev1[3] << "\n";
        cerr << "EV2:" << ev2[1] << "," << ev2[2] << "," << ev2[3] << "\n";
        cerr << "EV3:" << ev3[1] << "," << ev3[2] << "," << ev3[3] << "\n";
        cerr << "\n";
    }
}

} // namespace Makena

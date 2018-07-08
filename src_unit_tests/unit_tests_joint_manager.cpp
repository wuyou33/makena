#include "gtest/gtest.h"
#include "joint_manager.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace Makena {


class JointManagerTest : public ::testing::Test {

  protected:  

    JointManagerTest(){;}
    virtual ~JointManagerTest(){;}
    virtual void SetUp() {;};
    virtual void TearDown() {;};

};


static double rand100()
{
    return ((rand()%10000)+1)/100.0;
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


/**  @brief 
 */
TEST_F(JointManagerTest, Test01) {

    ConvexRigidBody body1(1), body2(2);

    Joint jt01(&body1, &body2);

    for (long i=0; i < 1000; i++) {

        auto val01 = rand100()/50.0 - 1.0;
        auto res01 = jt01.asinSafe(val01);
        EXPECT_EQ(res01, asin(val01));

        auto val02  = rand100()/5.0 + 1.0;
        auto res02 = jt01.asinSafe(val02);

        EXPECT_EQ(res02, asin(1.0));

        auto val03  = -1.0* (rand100()/5.0 + 1.0);
        auto res03 = jt01.asinSafe(val03);

        EXPECT_EQ(res03, asin(-1.0));


    }
}


/**  @brief 
 */
TEST_F(JointManagerTest, Test02) {

    ConvexRigidBody body1(1), body2(2);

    Joint jt01(&body1, &body2);

    for (long i=0; i < 1000; i++) {

        Vec3   x01(1.0, 0.0, 0.0);
        Vec3   y01(0.0, 1.0, 0.0);
        double rand01 = M_PI*2.0*(rand100()/100.0) - M_PI;
        Vec3   v01(0.0, cos(rand01) , sin(rand01));

        auto   q01 = randomRotQ();
        x01 = q01.rotate(x01);
        y01 = q01.rotate(y01);
        v01 = q01.rotate(v01);

        auto res01 = jt01.angleRad(y01, v01, x01);
        EXPECT_EQ(fabs(res01-rand01)<=EPSILON_SQUARED*10.0, true);

    }
}


/**  @brief compilation (instantiation) tests
 */
TEST_F(JointManagerTest, Test03) {

    ConvexRigidBody body1(1), body2(2);
    Vec3 vDummy;
    PistonJoint pj01(&body1, &body2,
        vDummy, vDummy, vDummy, vDummy, vDummy, vDummy,
        false,  0.0, 0.0, 
        false,  0.0, 0.0, 
        false,  0.0,
        false,  0.0, 0.0, 
        false,  0.0, 0.0, 
        false,  0.0,
        0.0);

}


}// Makena

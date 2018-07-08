#include <time.h>
#include "gtest/gtest.h"
#include "broad_phase_aabb_collision_detector.hpp"
#include "convex_rigid_body.hpp"

namespace Makena {

 
class BroadPhaseAABBCollisionDetectorTests : public ::testing::Test {
 
  protected:
    BroadPhaseAABBCollisionDetectorTests(){;}
    virtual ~BroadPhaseAABBCollisionDetectorTests(){;}    
    virtual void SetUp() {;};
    virtual void TearDown() {;};

};


static bool AABBcrudeTestPair(AABBelem* e1, AABBelem* e2) {

//    cerr << "E1 Low:  (" << e1->mLow.x()  << ", " << e1->mLow.y()  << ", " << e1->mLow.z()  << ")\n";
//    cerr << "E1 High: (" << e1->mHigh.x() << ", " << e1->mHigh.y() << ", " << e1->mHigh.z() << ")\n";
//    cerr << "E2 Low:  (" << e2->mLow.x()  << ", " << e2->mLow.y()  << ", " << e2->mLow.z()  << ")\n";
//    cerr << "E2 High: (" << e2->mHigh.x() << ", " << e2->mHigh.y() << ", " << e2->mHigh.z() << ")\n";

    if ( (e2->mHighX < e1->mLowX) || (e1->mHighX < e2->mLowX)   ) {
//        cerr << "abort on X\n";
        return false;
    }
    if ( (e2->mHighY < e1->mLowY) || (e1->mHighY < e2->mLowY)   ) {
//        cerr << "abort on Y\n";
        return false;
    }
    if ( (e2->mHighZ < e1->mLowZ) || (e1->mHighZ < e2->mLowZ)   ) {
//        cerr << "abort on Z\n";
        return false;
    }
//    cerr << "survived\n";
    return true;

}

static double dRandPositive(){ return double(rand()%10000) / 10000.0; }

static void generateCoordinates(AABBelem* e) {

    Vec3 p1( dRandPositive()*1000.0,
             dRandPositive()*1000.0,
             dRandPositive()*1000.0 );

    Vec3 r1( dRandPositive()*50.0,
             dRandPositive()*50.0,
             dRandPositive()*50.0 );
    e->setCoords(p1, p1 + r1);
}


static void updateCoordinates(ConvexRigidBody* const p) {
    AABBelem* e = &(p->aabb());
    Vec3 offset( dRandPositive()*1.0 - 2.5,
                 dRandPositive()*1.0 - 2.5,
                 dRandPositive()*1.0 - 2.5 );
    Vec3 Low (e->mLowX,  e->mLowY,  e->mLowZ);
    Vec3 High(e->mHighX, e->mHighY, e->mHighZ);
    e->setCoords(Low + offset, High + offset);
}



static bool isEquivalentTo(
    std::vector<pair<ConvexRigidBody*const, ConvexRigidBody*const> >& v1,
    std::vector<pair<ConvexRigidBody*const, ConvexRigidBody*const> >& v2
) {
    if (v1.size()!=v2.size()) {
        return false;
    }
    else {
        std::map<pair<ConvexRigidBody*const, ConvexRigidBody*const>,long> m;
        for (auto& p : v1) {
            if (m.find(p)==m.end()) {
                m[p] = 1;
            }
            else {
                return false;
            }
        }
        for (auto& p : v2) {
            if (m.find(p)!=m.end()) {
                m[p] = m[p] + 1;
            }
            else {
                return false;
            }
        }
        for (auto& kv : m) {
            if (kv.second !=2) {
                return false;                
            }
        }
        return true;
    }
}


TEST_F(BroadPhaseAABBCollisionDetectorTests, Test1) {

    BroadPhaseAABBCollisionDetector detector;

    clock_t tStart, tEnd;
    double test2 = 0.0;
    double test3 = 0.0;
    double numCollisions = 0;
    double numIterations = 100.0;
    long   numObjects = 1000;

    vector<AABBelem*> elements;
    vector<ConvexRigidBody*> objs;

    for (size_t i = 0; i < numObjects ; i++) {
        ConvexRigidBody* po    = new ConvexRigidBody(i);
        generateCoordinates(&(po->aabb()));
        detector.registerElem(&(po->aabb()));
        elements.push_back(&(po->aabb()));
        objs.push_back(po);
    }    

    for (size_t r = 0; r < 100000 ; r++) {
        long i = rand()%numObjects;
        long j = rand()%numObjects;
        if (i!=j) {
            elements[i]->exclusionList().insert(elements[j]);
            elements[j]->exclusionList().insert(elements[i]);
        }
    }    

    detector.update();

    for (size_t k = 0.0; k < numIterations ; k+=1.0) {
        for (auto* obj : objs) {
            updateCoordinates(obj);
        }    
        tStart = clock();
        std::vector<pair<ConvexRigidBody*const, ConvexRigidBody*const> > 
                            collisionPairs = detector.findCollidingPairsNbyN();
        tEnd = clock();
        test2 += (tEnd - tStart);
        numCollisions += collisionPairs.size();

        tStart = clock();

        detector.update();

        tEnd = clock();
        test3 += (tEnd - tStart);
        std::list<pair<ConvexRigidBody*const, ConvexRigidBody*const> >& 
                                          collisionPairs3 = detector.pairs();
                             
        std::vector<pair<ConvexRigidBody*const, ConvexRigidBody*const> > 
               collisionPairs4(collisionPairs3.begin(), collisionPairs3.end());
                              
        EXPECT_EQ(isEquivalentTo(collisionPairs, collisionPairs4), true);

    }

    std::cerr << "Avg num collisions: : "<<numCollisions/numIterations << "\n";
    std::cerr << "Test NbyN: " << test2/(1000.0*numIterations) << "\n";
    std::cerr << "Test SweepAndPrune: " <<test3/(1000.0*numIterations) << "\n";

}


} // namespace Makena

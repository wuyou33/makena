#include <time.h>
#include "gtest/gtest.h"
#include "red_black_tree.hpp"

namespace Makena {

 
class RedBlackTreeTests : public ::testing::Test {
 
  protected:
    RedBlackTreeTests(){;}
    virtual ~RedBlackTreeTests(){;}

    virtual void SetUp() {;};
    virtual void TearDown() {;};
};



class RBNCompElem {

  public:
    double mVal;
    double mIsEnter;

    inline bool operator<(const RBNCompElem& rhs) const  {
        return (mVal < rhs.mVal) || 
                ((mVal == rhs.mVal) && mIsEnter && !(rhs.mIsEnter) );
    }

    inline bool operator>(const RBNCompElem& rhs) const  {
        return (mVal > rhs.mVal) || 
                ((mVal == rhs.mVal) && !mIsEnter && rhs.mIsEnter );
    }

    inline bool operator==(const RBNCompElem& rhs) const  {
        return mVal == rhs.mVal && mIsEnter == rhs.mIsEnter;
    }


};

static vector<RBNCompElem> valuesExplored;

class RedBlackTree01 : public RedBlackTree<RBNCompElem> {
    inline void visitLeaf(RedBlackTreeNode<RBNCompElem>* n) override {
        valuesExplored.push_back(n->key());
    }
    inline void visitMidOrder(RedBlackTreeNode<RBNCompElem>* n) override {
        valuesExplored.push_back(n->key()); 
    }
};
 
static const int LENGTH = 10;

TEST_F(RedBlackTreeTests, Test1) {

    clock_t tStart, tEnd;
    double test1 = 0.0;
    double test2 = 0.0;
    double test3 = 0.0;
    double test4 = 0.0;
    for (int i = 0; i < 100; i++) {

        valuesExplored.clear();
        RedBlackTree01 tree;
        vector<RBNCompElem> values;
        for (size_t i = 0; i < LENGTH; i++) {
            RBNCompElem e;
            e.mVal = double(rand()%100000) / 1000.0;
            values.push_back(e);
        }

        vector<RBNCompElem> valuesSorted(values);
        std::sort(valuesSorted.begin(), valuesSorted.end());

        tStart = clock();
        for (auto& v : values) {
            auto* n = new RedBlackTreeNode<RBNCompElem>();
            n->setKey(v);
            tree.insert(n);
        }
        tEnd = clock();
        test1 += (tEnd - tStart);

        tStart = clock();
        tree.explore();
        tEnd = clock();
        test2 += (tEnd - tStart);
        for (size_t i = 0; i < valuesSorted.size(); i++) {
            EXPECT_EQ(valuesSorted[i], valuesExplored[i]);
        }

    }

//    cerr << "Test 1: " << test1 << "\n";
//    cerr << "Test 2: " << test2 << "\n";

}


} // namespace Makena

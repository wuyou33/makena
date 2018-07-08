#ifndef _MAKENA_RED_BLACK_TREE_HPP_
#define _MAKENA_RED_BLACK_TREE_HPP_

#include <exception>
#include <stdexcept>

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file red_black_tree.hpp
 *
 * @brief 
 *
 * @reference https://en.wikipedia.org/wiki/Red–black_tree
 *
 */
namespace Makena {

using namespace std;

/** @class RedBlackTreeNode
 *
 *  @brief tree node representation. It is expected to be a base class.
 *
 *  @return points along the convex hull in counter-clockwise ordering.
 */
template <class T>
class RedBlackTreeNode {

  public:

    inline RedBlackTreeNode():
        mParent(nullptr),
        mLeftChild(nullptr),
        mRightChild(nullptr),
        mThisIsRed(false),
        mState(NONE){;}

    inline ~RedBlackTreeNode(){;}

    inline void reset(){
        mParent     = nullptr;
        mLeftChild  = nullptr;
        mRightChild = nullptr;
        mThisIsRed  = false;
        mState      = NONE;
    };

    inline T key() const { return mKey;}

    inline void setKey(const T& k) { mKey = k; }

    inline bool isRed() const { return mThisIsRed; }
    inline void setRed() { mThisIsRed = true; }
    inline void setBlack() { mThisIsRed = false; }

    enum explorationState {
        NONE,
        INIT,
        VISITING_LEFT,
        VISITING_RIGHT,
        DONE
    };

    inline enum explorationState state() const { return mState; }
    inline void setState(enum explorationState s) { mState = s; }

    inline RedBlackTreeNode<T>* parent() const { return mParent; }
    inline void setParent(RedBlackTreeNode<T>* n) { mParent = n; }

    inline RedBlackTreeNode<T>* leftChild() const { return mLeftChild; }
    inline void setLeftChild(RedBlackTreeNode<T>* n) { mLeftChild = n; }

    inline RedBlackTreeNode<T>* rightChild() const { return mRightChild; }
    inline void setRightChild(RedBlackTreeNode<T>* n) { mRightChild = n; }

    inline RedBlackTreeNode<T>* grandParent() const { return mParent->mParent; }

    inline RedBlackTreeNode<T>* uncle() const {
        auto* gp = grandParent();
        if (gp->mLeftChild == mParent) {
            return gp->mRightChild;
        }
        else {
            return gp->mLeftChild;
        }
    }

    inline void rotateLeft() {
        //   P            P
        //   |            |
        //   N            R
        //  / \          / \
        //     R    =>  N  RL
        //    /
        //   RL
        bool hasParent            = (mParent!=nullptr);
        bool isLeftChildOfParent;
        if (hasParent) {
            isLeftChildOfParent   = (mParent->mLeftChild==this);
        }
        auto rightLeftGrandChild  = mRightChild->mLeftChild;
        auto rightChildOrg        = mRightChild;
        auto parentOrg            = mParent;
        mRightChild               = rightLeftGrandChild;
        if (rightLeftGrandChild != nullptr) {
            rightLeftGrandChild->mParent = this;
        }
        mParent                      = rightChildOrg;
        rightChildOrg->mLeftChild = this;
        if (hasParent) {
            if (isLeftChildOfParent) {
                parentOrg->mLeftChild = rightChildOrg;
            }
            else {
                parentOrg->mRightChild = rightChildOrg;
            }
        }
        rightChildOrg->mParent = parentOrg;
    }

    inline void rotateRight() {
        //   P          P
        //   |          |
        //   N          L
        //  /          / \
        // L      =>  LR  N
        //  \
        //   LR
        bool hasParent            = (mParent!=nullptr);
        bool isLeftChildOfParent;
        if (hasParent) {
            isLeftChildOfParent   = (mParent->mLeftChild==this);
        }
        auto leftRightGrandChild = mLeftChild->mRightChild;
        auto leftChildOrg        = mLeftChild;
        auto parentOrg           = mParent;
        mLeftChild               = leftRightGrandChild;
        if (leftRightGrandChild != nullptr) {
            leftRightGrandChild->mParent = this;
        }
        mParent                      = leftChildOrg;

        leftChildOrg->mRightChild = this;
        if (hasParent) {
            if (isLeftChildOfParent) {
                parentOrg->mLeftChild = leftChildOrg;
            }
            else {
                parentOrg->mRightChild = leftChildOrg;
            }
        }
        leftChildOrg->mParent = parentOrg;
    }

  private:

    RedBlackTreeNode<T>*     mParent;
    RedBlackTreeNode<T>*     mLeftChild;
    RedBlackTreeNode<T>*     mRightChild;
    bool                  mThisIsRed;
    T                     mKey;

    enum explorationState mState;

#ifdef UNIT_TESTS
friend class RedBlackTreeTests;
#endif

};


/** @class RedBlackTree
 *
 *  @brief red black tree representation. It is expected to be a base class.
 *         It is designed for broad-phase AABB collision culling.
 *
 *  @remark remove operation not implemented.
 */
template <class T>
class RedBlackTree {

  public:

    RedBlackTree():mRoot(nullptr){;}

    inline void insert(RedBlackTreeNode<T>* n);
    inline void explore();

    inline virtual void visitLeaf(RedBlackTreeNode<T>* n){;}
    inline virtual void visitPreOrder(RedBlackTreeNode<T>* n){;}
    inline virtual void visitMidOrder(RedBlackTreeNode<T>* n){;}
    inline virtual void visitPostOrder(RedBlackTreeNode<T>* n){;}

  private:
    RedBlackTreeNode<T>* mRoot;
    inline void balance_after_insersion(RedBlackTreeNode<T>* n);
    inline bool balance_one_red_node(RedBlackTreeNode<T>* n);

#ifdef UNIT_TESTS
  public:
    long  mDepth;
friend class RedBlackTreeTests;
#endif
};

template <class T>
void RedBlackTree<T>::explore()
{
    if (mRoot==nullptr) {
        return;
    }

    mRoot->setState(RedBlackTreeNode<T>::INIT);

    auto cur = mRoot;

#ifdef UNIT_TESTS
    long curDepth = 1;

    mDepth = curDepth;
#endif
    while(cur != nullptr) {

        switch (cur->state()) {
          case RedBlackTreeNode<T>::INIT:

            if (cur->leftChild()==nullptr) {
                if (cur->rightChild()==nullptr) {
                    visitLeaf(cur);
                    cur = cur->parent();
#ifdef UNIT_TESTS
                    curDepth--;
#endif
                }
                else {
                    visitPreOrder(cur);
                    visitMidOrder(cur);
                    cur->setState(RedBlackTreeNode<T>::VISITING_RIGHT);
                    cur->rightChild()->setState(RedBlackTreeNode<T>::INIT);
                    cur = cur->rightChild();
#ifdef UNIT_TESTS
                    curDepth++;
                    mDepth = max(curDepth,mDepth);
#endif
                }
            }
            else {
                visitPreOrder(cur);
                cur->setState(RedBlackTreeNode<T>::VISITING_LEFT);
                cur->leftChild()->setState(RedBlackTreeNode<T>::INIT);
                cur = cur->leftChild();
#ifdef UNIT_TESTS
                curDepth++;
                mDepth = max(curDepth,mDepth);
#endif
            }
            break;
            
          case RedBlackTreeNode<T>::VISITING_LEFT:

            visitMidOrder(cur);

            if (cur->rightChild()!=nullptr) {
                cur->setState(RedBlackTreeNode<T>::VISITING_RIGHT);
                cur->rightChild()->setState(RedBlackTreeNode<T>::INIT);
                cur = cur->rightChild();
#ifdef UNIT_TESTS
                curDepth++;
                mDepth = max(curDepth,mDepth);
#endif
            }
            else {
                visitPostOrder(cur);
                cur = cur->parent();
#ifdef UNIT_TESTS
                curDepth--;
#endif
            }
            break;

          case RedBlackTreeNode<T>::VISITING_RIGHT:
          default:

            visitPostOrder(cur);
            cur = cur->parent();
#ifdef UNIT_TESTS
            curDepth--;
#endif
            break;

        }

    }
}

template <class T>
bool RedBlackTree<T>::balance_one_red_node(RedBlackTreeNode<T>* n)
{
    if (n==mRoot) {
        n->setBlack();
        return true;
    }
    if (n->grandParent()==nullptr) {
        return true;
    }
    else if (!(n->parent()->isRed())) {
        return true;
    }
    auto* uncle       = n->uncle();
    auto* parent      = n->parent();
    auto* grandParent = n->grandParent();

    bool uncleIsNullOrBlack = false;

    if (uncle==nullptr) {
        uncleIsNullOrBlack = true;
    }
    else if (uncle->isRed()) {
        uncleIsNullOrBlack = false;
    }

    if (!uncleIsNullOrBlack) {
        grandParent->setRed();
        uncle->setBlack();
        parent->setBlack();
        return false;
    }
    else {
        if (parent->leftChild()==n) {
            if (grandParent->leftChild() == parent) {
                //        g           p
                //       / \         / \
                //      p*  u* =>   x*  g*
                //     /                 \
                //    n*                  n
                grandParent->rotateRight();
                parent->setBlack();
                grandParent->setRed();
                if (grandParent==mRoot) {
                    mRoot = parent;
                }
            }
            else {
                //        g           g             n
                //       / \         / \           / \
                //      u   p*  =>  u   n*   =>   g*  p*
                //         /             \       /
                //        n*              p*    u
                parent->rotateRight();
                grandParent->rotateLeft();
                n->setBlack();
                grandParent->setRed();
                if (grandParent==mRoot) {
                    mRoot = n;
                }
            }
        }
        else {
            if (grandParent->leftChild() == parent) {
                //        g           g          n
                //       / \         / \        / \
                //      p*  u  =>   n*  u  =>  p*  g*
                //       \         /                \
                //        n*      p*                 u 
                parent->rotateLeft();
                grandParent->rotateRight();
                n->setBlack();
                grandParent->setRed();
                if (grandParent==mRoot) {
                    mRoot = n;
                }
            }
            else {
                //        g             p
                //       / \           / \
                //      u   p*   =>   g*  x*
                //           \       /
                //            x*    u
                grandParent->rotateLeft();
                parent->setBlack();
                grandParent->setRed();
                if (grandParent==mRoot) {
                    mRoot = parent;
                }
            }
        }
        return true;
    }
}

template <class T>
void RedBlackTree<T>::balance_after_insersion(RedBlackTreeNode<T>* n)
{
    auto* cur = n;
    while(true) {
        auto finished = balance_one_red_node(n);
        if (finished) {
            break;
        }
        cur = cur->grandParent();
    }
}

template <class T>
void RedBlackTree<T>::insert(RedBlackTreeNode<T>* n)
{
    n->reset();
    if (mRoot==nullptr) {
        mRoot = n;
        n->setParent(nullptr);
    }
    else {
        auto* cur = mRoot;
        while (true) {
            if (n->key() < cur->key()) {
                if (cur->leftChild()==nullptr) {
                    cur->setLeftChild(n);
                    n->setParent(cur);
                    break;
                }
                else {
                    cur = cur->leftChild();
                }
            }
            else {
                if (cur->rightChild()==nullptr) {
                    cur->setRightChild(n);
                    n->setParent(cur);
                    break;
                }
                else {
                    cur = cur->rightChild();
                }
            }
        }
        n->setRed();
        balance_after_insersion(n);
    }
}

} //namespace Makena

#endif /*_MAKENA_RED_BLACK_TREE_HPP_*/

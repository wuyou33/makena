#include "intersection_finder.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file intersection_finder.cpp
 *
 * @brief find the interseciton of two convex polytopes
 *        the result is a convex polytope, convex polygon, edge, or a point.
 */


namespace Makena {

using namespace std;


IntersectionFinder::IntersectionFinder(
    Manifold&     m1,
    const Mat3x3& Qmat1, 
    const Vec3&   CoM1,
    Manifold&     m2,
    const Mat3x3& Qmat2, 
    const Vec3&   CoM2,
    const double& epsilonZero,
    const double& epsilonZeroPCA,
    const double& epsilonAngle,
    std::ostream& logStream
):
        Loggable(logStream),
        mM1(m1),
        mM2(m2),
        mQmat1(Qmat1),
        mQmat2(Qmat2),
        mCoM1(CoM1),
        mCoM2(CoM2),
        mEpsilonZero(epsilonZero),
        mEpsilonZeroPCA(epsilonZeroPCA),
        mEpsilonAngle(epsilonAngle)
        {;}


IntersectionFinder::~IntersectionFinder(){;}



void IntersectionFinder::find(const Vec3 sepAxis)
{
    prepareVertices(sepAxis);

    prepareEdgesAndFaces();

    findIntersectionsEdgesOn1FacesOn2();

    findIntersectionsFacesOn1EdgesOn2();

    findRemainingInteriorVertices();

    constructVertexAttributes();

    findDimension();

    if (mDimension == -1) {        
        ;
    }
    else if (mDimension == 0) {        

        processZeroDimensionalIntersection();// Point

    }
    else if (mDimension == 1) {   // Edge : Find two extramal points.

        processOneDimensionalIntersection();

    }
    else if (mDimension == 2) {   // 2D Convex Polygon

        processTwoDimensionalIntersection();

    }
    else { // Convex Polytope

        processThreeDimensionalIntersection();

    }
}


void IntersectionFinder::find()
{
    prepareVertices();

    prepareEdgesAndFaces();

    findIntersectionsEdgesOn1FacesOn2();

    findIntersectionsFacesOn1EdgesOn2();

    findRemainingInteriorVertices();

    constructVertexAttributes();

    findDimension();

    if (mDimension == -1) {        
        ;
    }
    else if (mDimension == 0) {        

        processZeroDimensionalIntersection();// Point

    }
    else if (mDimension == 1) {   // Edge : Find two extramal points.

        processOneDimensionalIntersection();

    }
    else if (mDimension == 2) {   // 2D Convex Polygon

        processTwoDimensionalIntersection();

    }
    else { // Convex Polytope

        processThreeDimensionalIntersection();

    }
}


void IntersectionFinder::prepareVertices(const Vec3 sepAxis)
{
    auto vPair1 = mM1.vertices();
    double min1, max1;
    bool first = true;
    for (auto vit = vPair1.first; vit != vPair1.second; vit++) {
        auto dot = sepAxis.dot((*vit)->pGCS(mQmat1, mCoM1));
        (*vit)->IFsetDot(dot);
        (*vit)->IFreset();
        if (first) {
            min1 = dot;
            max1 = dot;
            first = false;
        }
        else {
            min1 = std::min(min1, dot);
            max1 = std::max(max1, dot);
        }
    }

    auto vPair2 = mM2.vertices();
    double min2, max2;
    first = true;
    for (auto vit = vPair2.first; vit != vPair2.second; vit++) {
        auto dot = sepAxis.dot((*vit)->pGCS(mQmat2, mCoM2));
        (*vit)->IFsetDot(dot);
        (*vit)->IFreset();
        if (first) {
            min2 = dot;
            max2 = dot;
            first = false;
        }
        else {
            min2 = std::min(min2, dot);
            max2 = std::max(max2, dot);
        }
    }

    if ( min1 <= (min2+mEpsilonZero) &&
         max1 <= (max2+mEpsilonZero) &&
         min2 <= (max1+mEpsilonZero)    ) {

        //       min1===============max1    
        //                min2=============max2

        for (auto vit = vPair1.first; vit != vPair1.second; vit++) {
            if ((*vit)->IFdot() + mEpsilonZero >= min2) {
                (*vit)->IFsetActive();
                mActiveVertices1.push_back(vit);
            }            
        }

        for (auto vit = vPair2.first; vit != vPair2.second; vit++) {
            if ((*vit)->IFdot() <= max1 + mEpsilonZero) {
                (*vit)->IFsetActive();
                mActiveVertices2.push_back(vit);
            }            
        }
    }
    else if ( min2 <= (min1+mEpsilonZero) &&
              max2 <= (max1+mEpsilonZero) &&
              min1 <= (max2+mEpsilonZero)    ) {

        //       min2===============max2
        //                min1=============max1

        for (auto vit = vPair1.first; vit != vPair1.second; vit++) {
            if ((*vit)->IFdot() <= max2 + mEpsilonZero) {
                (*vit)->IFsetActive();
                mActiveVertices1.push_back(vit);
            }            
        }

        for (auto vit = vPair2.first; vit != vPair2.second; vit++) {
            if ((*vit)->IFdot() + mEpsilonZero >= min1) {
                (*vit)->IFsetActive();
                mActiveVertices2.push_back(vit);
            }            
        }
    }
    else if (max2 < min1 || max1 < min2) {// No need to add mEpsilonZeo to
                                          // terms here as they have been
                                          // already taken care of above.
        //     min2======max2
        //                     min1======max1
        ;// No overlap between two convex polytopes.
    }
    else {
        // Could not find a proper overlap. Mark all the vertices.
        for (auto vit = vPair1.first; vit != vPair1.second; vit++) {
            (*vit)->IFsetActive();
            mActiveVertices1.push_back(vit);
        }
        for (auto vit = vPair2.first; vit != vPair2.second; vit++) {
            (*vit)->IFsetActive();
            mActiveVertices2.push_back(vit);
        }
    }
}


void IntersectionFinder::prepareVertices()
{
    auto vPair1 = mM1.vertices();
    for (auto vit = vPair1.first; vit != vPair1.second; vit++) {
        (*vit)->IFreset();
        (*vit)->IFsetActive();
        mActiveVertices1.push_back(vit);
    }

    auto vPair2 = mM2.vertices();
    for (auto vit = vPair2.first; vit != vPair2.second; vit++) {
        (*vit)->IFreset();
        (*vit)->IFsetActive();
        mActiveVertices2.push_back(vit);
    }
}


void IntersectionFinder::prepareEdgesAndFaces()
{
    auto ePair1 = mM1.edges();
    for (auto eit = ePair1.first; eit != ePair1.second; eit++) {
        auto heit = (*eit)->he1();
        auto src  = (*heit)->src();
        auto dst  = (*heit)->dst();
        if ((*src)->IFisActive()||(*dst)->IFisActive()) {
            mActiveEdges1.push_back(eit);
        }
    }

    auto ePair2 = mM2.edges();
    for (auto eit = ePair2.first; eit != ePair2.second; eit++) {
        auto heit = (*eit)->he1();
        auto src  = (*heit)->src();
        auto dst  = (*heit)->dst();
        if ((*src)->IFisActive()||(*dst)->IFisActive()) {
            mActiveEdges2.push_back(eit);
        }
    }

    auto fPair1 = mM1.faces();
    for (auto fit = fPair1.first; fit != fPair1.second; fit++) {
        for (auto heit : (*fit)->halfEdges()) {
            auto src = (*heit)->src();
            if ((*src)->IFisActive()) {
                mActiveFaces1.push_back(fit);
                break;
            }
        }
    }

    auto fPair2 = mM2.faces();
    for (auto fit = fPair2.first; fit != fPair2.second; fit++) {
        for (auto heit : (*fit)->halfEdges()) {
            auto src = (*heit)->src();
            if ((*src)->IFisActive()) {
                mActiveFaces2.push_back(fit);
                break;
            }
        }
    }
}


Mat3x3 IntersectionFinder::rotMatAlignZDirToZAxis(const Vec3& zDir)
{
    // Find the new coordinate axes with zDir for new z-axis.
    Vec3 r3 = zDir;
    double ax = fabs(r3.x());
    double ay = fabs(r3.y());
    double az = fabs(r3.z());

    Vec3 r1;
    if (ax <= ay && ay <= az) {
      Vec3 xUnit(1.0, 0.0, 0.0);
        r1 = r3.cross(xUnit);
    }
    else if (ay <= ax && ay <= az) {
        Vec3 yUnit(0.0, 1.0, 0.0);
        r1 = r3.cross(yUnit);
    }
    else {
        Vec3 zUnit(0.0, 0.0, 1.0);
        r1 = r3.cross(zUnit);
    }
    r1.normalize();
    // r3 x r1 = r2 => r1 x r2 = r3
    const Vec3 r2 = r3.cross(r1);

    // Find the rotation matrix to rotate a point in GCS into the new CS.
    Mat3x3 rMat(r1, r2, r3);
    rMat.transposeInPlace();
    return rMat;
}


long IntersectionFinder::edgeIndex(long i1, long i2, long numElems)
{
    if ( ( i1 == 0 && i2 == 1 ) || ( i1 == 1 && i2 == 0 ) ) {

        return 0;
    }
    else if ( i1 == 0 || i2 == 0 ) {

        return numElems - 1;
    }
    else {

        return std::min(i1, i2);
    }
}


void IntersectionFinder::rotateWorld(
    const Vec3& n,
    Mat3x3&     newQmat1,
    Vec3&       newCoM1,
    Mat3x3&     newQmat2,
    Vec3&       newCoM2,
    Mat3x3&     matAlignInv
) {
    auto matAlign = rotMatAlignZDirToZAxis(n);
    matAlignInv   = matAlign.transpose();
    newQmat1 = matAlign * mQmat1;
    newCoM1  = matAlign * mCoM1;
    newQmat2 = matAlign * mQmat2;
    newCoM2  = matAlign * mCoM2;
}


void IntersectionFinder::processIntersection2DMid(
    EdgeIt                       eit1,
    vector<IntSec2D::InputElem>& elems1,
    const Vec3&                  pGCS,
    FaceIt                       fit2,
    vector<IntSec2D::InputElem>& elems2,
    vector<VertexIt>&            vertices2,
    vector<EdgeIt>&              edges2
) {

    IntSec2D finder( mEpsilonZero, mEpsilonZero, mEpsilonAngle, mLogStream );
                     
    vector<IntSec2D::OutputElem> intsec;

    finder.findIntersection(elems1, elems2, intsec);

    for (auto& oe : intsec) {
        switch(oe.mType) {
          case IntSec2D::IT_VERTEX_EDGE:
            {
                auto eit2 = edges2[ edgeIndex(
                                   oe.mIndexB, oe.mIndexBaux, elems2.size()) ];
                auto keyEE = make_pair((*eit1)->id(), (*eit2)->id());
                if (mMapEE.find(keyEE) == mMapEE.end()) {
                    mMapEE[keyEE] = pGCS;
                }
            }
            break;

          case IntSec2D::IT_VERTEX_VERTEX:
            {
                auto vit2 = vertices2[oe.mIndexB];
                auto keyEV = make_pair((*eit1)->id(), (*vit2)->id());
                if (mMapEV.find(keyEV) == mMapEV.end()) {
                    mMapEV[keyEV] = pGCS;
                }
                (*vit2)->IFsetProcessed();
            }
            break;

          case IntSec2D::IT_VERTEX_INTERIOR:
            {
                auto keyEF = make_pair((*eit1)->id(), (*fit2)->id());
                if (mMapEF.find(keyEF) == mMapEF.end()) {
                    mMapEF[keyEF] = pGCS;
                }
            }
            break;
          default:
            break;
        }
    }
}


void IntersectionFinder::processIntersection2DMid(
    FaceIt                       fit1,
    vector<IntSec2D::InputElem>& elems1,
    vector<VertexIt>&            vertices1,
    vector<EdgeIt>&              edges1,
    EdgeIt                       eit2,
    vector<IntSec2D::InputElem>& elems2,
    const Vec3&                  pGCS
) {

    IntSec2D finder( mEpsilonZero, mEpsilonZero, mEpsilonAngle, mLogStream );
                     
    vector<IntSec2D::OutputElem> intsec;

    finder.findIntersection(elems1, elems2, intsec);

    for (auto& oe : intsec) {
        switch(oe.mType) {
          case IntSec2D::IT_EDGE_VERTEX:
            {
                auto eit1 = edges1[ edgeIndex(
                                   oe.mIndexA, oe.mIndexAaux, elems1.size()) ];
                auto keyEE = make_pair((*eit1)->id(), (*eit2)->id());
                if (mMapEE.find(keyEE) == mMapEE.end()) {
                    mMapEE[keyEE] = pGCS;
                }
            }
            break;

          case IntSec2D::IT_VERTEX_VERTEX:
            {
                auto vit1 = vertices1[oe.mIndexA];
                auto keyVE = make_pair((*vit1)->id(), (*eit2)->id());
                if (mMapVE.find(keyVE) == mMapVE.end()) {
                    mMapVE[keyVE] = pGCS;
                }
                (*vit1)->IFsetProcessed();
            }
            break;

          case IntSec2D::IT_INTERIOR_VERTEX:
            {
                auto keyFE = make_pair((*fit1)->id(), (*eit2)->id());
                if (mMapFE.find(keyFE) == mMapFE.end()) {
                    mMapFE[keyFE] = pGCS;
                }
            }
            break;
          default:
            break;
        }
    }
}


void IntersectionFinder::processIntersection2DEnd(
    FaceIt                       fit1,
    vector<IntSec2D::InputElem>& elems1,
    vector<VertexIt>&            vertices1,
    vector<EdgeIt>&              edges1,
    FaceIt                       fit2,
    vector<IntSec2D::InputElem>& elems2,
    vector<VertexIt>&            vertices2,
    vector<EdgeIt>&              edges2,
    const Vec3&                  pGCS
) {

    IntSec2D finder( mEpsilonZero, mEpsilonZero, mEpsilonAngle, mLogStream );
                     
    vector<IntSec2D::OutputElem> intsec;
    finder.findIntersection(elems1, elems2, intsec);

    for (auto& oe : intsec) {
        switch(oe.mType) {
          case IntSec2D::IT_EDGE_VERTEX:
            {
                auto eit1 = edges1[ edgeIndex(
                                   oe.mIndexA, oe.mIndexAaux, elems1.size()) ];
                auto vit2 = vertices2[oe.mIndexB];
                auto keyEV = make_pair((*eit1)->id(), (*vit2)->id());
                if (mMapEV.find(keyEV) == mMapEV.end()) {
                    mMapEV[keyEV] = pGCS;
                }
                (*vit2)->IFsetProcessed();
            }
            break;

          case IntSec2D::IT_VERTEX_EDGE:
            {
                auto vit1 = vertices1[oe.mIndexA];
                auto eit2 = edges2[ edgeIndex(
                                   oe.mIndexB, oe.mIndexBaux, elems2.size()) ];
                auto keyVE = make_pair((*vit1)->id(), (*eit2)->id());
                if (mMapVE.find(keyVE) == mMapVE.end()) {
                    mMapVE[keyVE] = pGCS;
                }
                (*vit1)->IFsetProcessed();
            }
            break;

          case IntSec2D::IT_VERTEX_VERTEX:
            {
                auto vit1 = vertices1[oe.mIndexA];
                auto vit2 = vertices2[oe.mIndexB];
                auto keyVV = make_pair((*vit1)->id(), (*vit2)->id());
                if (mMapVV.find(keyVV) == mMapVV.end()) {
                    mMapVV[keyVV] = pGCS;
                }
                (*vit1)->IFsetProcessed();
                (*vit2)->IFsetProcessed();
            }
            break;

          case IntSec2D::IT_VERTEX_INTERIOR:
            {
                auto vit1 = vertices1[oe.mIndexA];
                auto keyVF = make_pair((*vit1)->id(), (*fit2)->id());
                if (mMapVF.find(keyVF) == mMapVF.end()) {
                    mMapVF[keyVF] = pGCS;
                }
                (*vit1)->IFsetProcessed();
            }
            break;

          case IntSec2D::IT_INTERIOR_VERTEX:
            {
                auto vit2 = vertices2[oe.mIndexB];
                auto keyFV = make_pair((*fit1)->id(), (*vit2)->id());
                if (mMapFV.find(keyFV) == mMapFV.end()) {
                    mMapFV[keyFV] = pGCS;
                }
                (*vit2)->IFsetProcessed();
            }
            break;
          default:
            break;
        }
    }
}


Vec3 IntersectionFinder::interpolatePoint(
    IntSec2D::InputElem&  ie1,
    IntSec2D::InputElem&  ie2,
    IntSec2D::OutputElem& oe,
    const Vec3&           p1_3D,
    const Vec3&           p2_3D
) {                                 
    const auto& p1_2D   = ie1.mP;
    const auto& p2_2D   = ie2.mP;
    const auto& pOut_2D = oe.mP;

    const auto deltaX = p2_2D.x() - p1_2D.x();
    const auto deltaY = p2_2D.y() - p1_2D.y();

    double alpha = 0.0;
    if (fabs(deltaX) > fabs(deltaY)) {
        alpha = (pOut_2D.x() - p1_2D.x())/deltaX;
    }
    else if (fabs(deltaY) > mEpsilonZero) {
        alpha = (pOut_2D.y() - p1_2D.y())/deltaY;
    }
    return p2_3D * alpha + p1_3D * (1.0 - alpha);
}                                 


void IntersectionFinder::processIntersection2DCoplanar(
    FaceIt                       fit1,
    vector<IntSec2D::InputElem>& elems1,
    vector<VertexIt>&            vertices1,
    vector<EdgeIt>&              edges1,
    FaceIt                       fit2,
    vector<IntSec2D::InputElem>& elems2,
    vector<VertexIt>&            vertices2,
    vector<EdgeIt>&              edges2,
    bool                         edgeOn1,
    const Vec3&                  p0GCS,
    const Vec3&                  p1GCS
) {
    IntSec2D finder( mEpsilonZero, mEpsilonZero, mEpsilonAngle, mLogStream );
                     
    vector<IntSec2D::OutputElem> intsec;
    finder.findIntersection(elems1, elems2, intsec);

    for (auto& oe : intsec) {

        auto p = interpolatePoint( edgeOn1? elems1[0] : elems2[0],
                                   edgeOn1? elems1[1] : elems2[1],
                                   oe, p0GCS, p1GCS               );
        switch (oe.mType) {
          case IntSec2D::IT_EDGE_EDGE:
            {
                auto eit1 = edges1[ edgeIndex(
                                   oe.mIndexA, oe.mIndexAaux, elems1.size()) ];
                auto eit2 = edges2[ edgeIndex(
                                   oe.mIndexB, oe.mIndexBaux, elems2.size()) ];
                auto key = make_pair((*eit1)->id(), (*eit2)->id());
                if (mMapEE.find(key) == mMapEE.end()) {
                    mMapEE[key] = p;
                }
            }
            break;

          case IntSec2D::IT_EDGE_VERTEX:
            {
                auto eit1 = edges1[ edgeIndex(
                                   oe.mIndexA, oe.mIndexAaux, elems1.size()) ];
                auto vit2 = vertices2[oe.mIndexB];
                auto key = make_pair((*eit1)->id(), (*vit2)->id());
                if (mMapEV.find(key) == mMapEV.end()) {
                    mMapEV[key] = p;
                }
                (*vit2)->IFsetProcessed();
            }
            break;

          case IntSec2D::IT_VERTEX_EDGE:
            {
                auto vit1 = vertices1[oe.mIndexA];
                auto eit2 = edges2[ edgeIndex(
                                   oe.mIndexB, oe.mIndexBaux, elems2.size()) ];
                auto key = make_pair((*vit1)->id(), (*eit2)->id());
                if (mMapVE.find(key) == mMapVE.end()) {
                    mMapVE[key] = p;
                }
                (*vit1)->IFsetProcessed();
            }
            break;

          case IntSec2D::IT_VERTEX_VERTEX:
            {

                auto vit1 = vertices1[oe.mIndexA];
                auto vit2 = vertices2[oe.mIndexB];
                auto key = make_pair((*vit1)->id(), (*vit2)->id());
                if (mMapVV.find(key) == mMapVV.end()) {
                    mMapVV[key] = p;
                }
                (*vit1)->IFsetProcessed();
                (*vit2)->IFsetProcessed();
            }
            break;

          case IntSec2D::IT_VERTEX_INTERIOR:
            {
                auto vit1 = vertices1[oe.mIndexA];
                auto key = make_pair((*vit1)->id(), (*fit2)->id());
                if (mMapVF.find(key) == mMapVF.end()) {
                    mMapVF[key] = p;
                }
                (*vit1)->IFsetProcessed();
            }
            break;

          case IntSec2D::IT_INTERIOR_VERTEX:
            {
                auto vit2 = vertices2[oe.mIndexB];
                auto key = make_pair((*fit1)->id(), (*vit2)->id());
                if (mMapFV.find(key) == mMapFV.end()) {
                    mMapFV[key] = p;
                }
                (*vit2)->IFsetProcessed();
            }
            break;
          default:
            break;
        }
    }
}


void IntersectionFinder::findEdgePointsToBeTested(
    const Vec3& pBase,
    const Vec3& pSrc,
    const Vec3& pDst,
    bool&       processSrc,
    bool&       processDst,
    bool&       processMid,
    Vec3&       pMid
) {

    processSrc = false;
    processDst = false;
    processMid = false;

    if ( fabs(pBase.z() - pSrc.z()) <= mEpsilonZero ) {

        if ( fabs(pBase.z() - pDst.z()) <= mEpsilonZero ) {

            // The edge is coplanar to the face
            processSrc = true;
            processDst = true;
        }
        else {
            // The src vertex is coplanar to the face
            processSrc = true;
        }
    }
    else {
        if ( fabs(pBase.z() - pDst.z()) <= mEpsilonZero ) {

            // The dst vertex is coplanar to the face
            processDst = true;

        }
        else {

            // The edge is not coplanar to the face

            if ( ( pBase.z() <= pSrc.z() && pBase.z() >= pDst.z() ) ||
                 ( pBase.z() <= pDst.z() && pBase.z() >= pSrc.z() )    ) {

                // The edge virtically overlaps with the face.
                // Test the proper interior of the edge 
                // for intersection against face.
                double coeff = (pBase.z() - pSrc.z()) / (pDst.z() - pSrc.z()) ;
                double pX    = pSrc.x() + (pDst.x() - pSrc.x())*coeff;
                double pY    = pSrc.y() + (pDst.y() - pSrc.y())*coeff;

                Vec3 pLevelZ(pX, pY, pBase.z());

                if ((pSrc - pLevelZ).squaredNorm2()<=mEpsilonZero) {
                    // The src vertex is coplanar to the face
                    processSrc = true;
                }
                else if ((pDst-pLevelZ).squaredNorm2()<=mEpsilonZero) {
                    // The dst vertex is coplanar to the face
                    processDst = true;
                }
                else {
                    processMid = true;
                    pMid = pLevelZ;
                }
            }
        }
    }
}


void IntersectionFinder::construct2DInputElementsForFace(
    FaceIt                       fit,
    const Mat3x3&                Qmat,
    const Vec3&                  CoM,
    vector<IntSec2D::InputElem>& elems,
    vector<VertexIt>&            vertices,
    vector<EdgeIt>&              edges,
    Vec3&                        pBase
) {
    elems.clear();
    vertices.clear();
    edges.clear();
    
    for (auto heit : (*fit)->halfEdges()) {

        auto fSrc = (*heit)->src();
        pBase     = (*fSrc)->pGCS(Qmat, CoM);

        IntSec2D::InputElem ie( Vec2(pBase.x(), pBase.y()), elems.size() );

        elems.   push_back( ie );
        vertices.push_back( fSrc );
        edges.   push_back( (*heit)->edge() );
        
    }
}


void IntersectionFinder::findIntersectionsEdgesOn1FacesOn2()
{
    for (auto fit : mActiveFaces2) {

        Mat3x3 newQmat1, newQmat2;
        Vec3   newCoM1,  newCoM2;
        Mat3x3 matAlignInv;

        rotateWorld( (*fit)->nGCS(mQmat2), 
                     newQmat1, newCoM1,
                     newQmat2, newCoM2,
                     matAlignInv            );

        vector<IntSec2D::InputElem> elems2;
        vector<VertexIt>            vertices2;
        vector<EdgeIt>              edges2;
        Vec3                        pBase;

        construct2DInputElementsForFace( fit, newQmat2, newCoM2,
                                         elems2, vertices2, edges2, pBase );

        for (auto eit : mActiveEdges1) {
            vector<IntSec2D::InputElem> elems1;
            vector<VertexIt>            vertices1;
            vector<EdgeIt>              edges1;

            auto heit    = (*eit)->he1();
            auto src     = (*heit)->src();
            auto dst     = (*heit)->dst();
            auto pSrc    = (*src)->pGCS(newQmat1, newCoM1);
            auto pDst    = (*dst)->pGCS(newQmat1, newCoM1);
            bool processSrc, processDst, processMid;
            Vec3 pMid;
            findEdgePointsToBeTested(
                  pBase, pSrc, pDst, processSrc, processDst, processMid, pMid);

            if (processMid) {
                IntSec2D::InputElem ie1(Vec2(pMid.x(), pMid.y()), 0);
                elems1.push_back( ie1 );

                processIntersection2DMid( eit, elems1, matAlignInv * pMid,
                                          fit, elems2, vertices2, edges2  );
            }
            else if (processSrc) {
                if (processDst) {
                    IntSec2D::InputElem ie1(Vec2(pSrc.x(), pSrc.y()), 0);
                    elems1.   push_back( ie1 );
                    vertices1.push_back( src );
                    edges1.   push_back( eit );

                    IntSec2D::InputElem ie2(Vec2(pDst.x(), pDst.y()), 1);
                    elems1.   push_back( ie2 );
                    vertices1.push_back( dst );
                    edges1.   push_back( eit );

                    processIntersection2DCoplanar(
                               fit, elems1, vertices1, edges1,
                               fit, elems2, vertices2, edges2,
                               true,
                               matAlignInv  * pSrc, matAlignInv * pDst
                    ); 
                }
                else {
                    if ((*src)->IFisActive()) {
                        IntSec2D::InputElem ie1(Vec2(pSrc.x(), pSrc.y()), 0);
                        elems1.   push_back( ie1 );
                        vertices1.push_back( src );
                        edges1.   push_back( eit );

                        processIntersection2DEnd(
                               fit, elems1, vertices1, edges1,
                               fit, elems2, vertices2, edges2,
                               matAlignInv * pSrc
                        ); 
                    }
                }
            }
            else {
                if (processDst) {
                    if ((*dst)->IFisActive()) {
                        IntSec2D::InputElem ie1(Vec2(pDst.x(), pDst.y()), 0);
                        elems1.   push_back( ie1 );
                        vertices1.push_back( dst );
                        edges1.   push_back( eit );

                        processIntersection2DEnd(
                               fit, elems1, vertices1, edges1,
                               fit, elems2, vertices2, edges2,
                               matAlignInv * pDst
                        );
                    }
                }
            }
        }
    }
}


void IntersectionFinder::findIntersectionsFacesOn1EdgesOn2()
{
    for (auto fit : mActiveFaces1) {

        Mat3x3 newQmat1, newQmat2;
        Vec3   newCoM1,  newCoM2;
        Mat3x3 matAlignInv;

        rotateWorld( (*fit)->nGCS(mQmat1), 
                     newQmat1, newCoM1,
                     newQmat2, newCoM2,
                     matAlignInv            );

        vector<IntSec2D::InputElem> elems1;
        vector<VertexIt>            vertices1;
        vector<EdgeIt>              edges1;
        Vec3                        pBase;

        construct2DInputElementsForFace( fit, newQmat1, newCoM1,
                                         elems1, vertices1, edges1, pBase );

        for (auto eit : mActiveEdges2) {

            vector<IntSec2D::InputElem> elems2;
            vector<VertexIt>            vertices2;
            vector<EdgeIt>              edges2;

            auto heit    = (*eit)->he1();
            auto src     = (*heit)->src();
            auto dst     = (*heit)->dst();
            auto pSrc    = (*src)->pGCS(newQmat2, newCoM2);
            auto pDst    = (*dst)->pGCS(newQmat2, newCoM2);

            bool processSrc, processDst, processMid;
            Vec3 pMid;
            findEdgePointsToBeTested(
                  pBase, pSrc, pDst, processSrc, processDst, processMid, pMid);

            if (processMid) {

                IntSec2D::InputElem ie1(Vec2(pMid.x(), pMid.y()), 0);
                elems2.push_back( ie1 );

                processIntersection2DMid( fit, elems1, vertices1, edges1, 
                                          eit, elems2, matAlignInv * pMid );
            }
            else if (processSrc) {
                if (processDst) {
                    IntSec2D::InputElem ie1(Vec2(pSrc.x(), pSrc.y()), 0);
                    elems2.   push_back( ie1 );
                    vertices2.push_back( src );
                    edges2.   push_back( eit );

                    IntSec2D::InputElem ie2(Vec2(pDst.x(), pDst.y()), 1);
                    elems2.   push_back( ie2 );
                    vertices2.push_back( dst );
                    edges2.   push_back( eit );

                    processIntersection2DCoplanar(
                               fit, elems1, vertices1, edges1,
                               fit, elems2, vertices2, edges2,
                               false,
                               matAlignInv * pSrc , matAlignInv * pDst
                    ); 
                }
                else {
                    if ((*src)->IFisActive()) {
                        IntSec2D::InputElem ie1(Vec2(pSrc.x(), pSrc.y()), 0);
                        elems2.   push_back( ie1 );
                        vertices2.push_back( src );
                        edges2.   push_back( eit );

                        processIntersection2DEnd(
                               fit, elems1, vertices1, edges1,
                               fit, elems2, vertices2, edges2,
                               matAlignInv * pSrc
                        ); 
                    }
                }
            }
            else {
                if (processDst) {
                    if ((*dst)->IFisActive()) {
                        IntSec2D::InputElem ie1(Vec2(pDst.x(), pDst.y()), 0);
                        elems2.   push_back( ie1 );
                        vertices2.push_back( dst );
                        edges2.   push_back( eit );

                        processIntersection2DEnd(
                               fit, elems1, vertices1, edges1,
                               fit, elems2, vertices2, edges2,
                               matAlignInv * pDst
                        );
                    }
                }
            }
        }
    }
}


enum predicate IntersectionFinder::classifyVertexAgainstPlane(
    const Vec3& pTest,
    const Vec3& pPlane,
    const Vec3& nPlane
) {
    Vec3 v = pTest - pPlane;
    if (v.squaredNorm2() <= mEpsilonZero) {
        return IF_ON_PLANE;
    }
    else {
        v.normalize();
        auto dot = v.dot(nPlane);
        if (fabs(dot) <= mEpsilonAngle) {
            return IF_ON_PLANE;
        }
        else if (dot < 0.0) {
            return IF_BACK_OF_PLANE;
        }
        else {
            return IF_FRONT_OF_PLANE;
        }
    }
}


void IntersectionFinder::findRemainingInteriorVertices()
{
    auto numFaces2 = mActiveFaces2.size();

    for (auto fit : mActiveFaces2) {
        auto heit = *(((*fit)->halfEdges()).begin());
        auto src  = (*heit)->src();
        auto fp   = (*src)->pGCS(mQmat2, mCoM2);
        auto fn   = (*fit)->nGCS(mQmat2);
   
        for (auto vit : mActiveVertices1) {
            if (!(*vit)->IFisProcessed()) {
                auto p    = (*vit)->pGCS(mQmat1, mCoM1);
                auto pred = classifyVertexAgainstPlane(p, fp, fn);

                if (pred == IF_BACK_OF_PLANE || pred == IF_ON_PLANE) {
                    (*vit)->IFincrement();
                    if ((*vit)->IFcnt() == numFaces2) {
                        mMapVI[(*vit)->id()] = p;
                        (*vit)->IFsetProcessed();
                    }
                }

            }
        }
    }

    auto numFaces1 = mActiveFaces1.size();

    for (auto fit : mActiveFaces1) {
        auto heit = *(((*fit)->halfEdges()).begin());
        auto src  = (*heit)->src();
        auto fp   = (*src)->pGCS(mQmat1, mCoM1);
        auto fn   = (*fit)->nGCS(mQmat1);
   
        for (auto vit : mActiveVertices2) {
            if (!(*vit)->IFisProcessed()) {
                auto p    = (*vit)->pGCS(mQmat2, mCoM2);
                auto pred = classifyVertexAgainstPlane(p, fp, fn);

                if (pred == IF_BACK_OF_PLANE || pred == IF_ON_PLANE) {
                    (*vit)->IFincrement();
                    if ((*vit)->IFcnt() ==  numFaces1) {
                        mMapIV[(*vit)->id()] = p;
                        (*vit)->IFsetProcessed();
                    }
                }

            }
        }
    }
}


void IntersectionFinder::constructVertexAttributes()
{
    for (auto mit : mMapVV) {
        Attributes a;
        a.mP    = mit.second;
        a.mPred = IF_VERTEX_VERTEX;
        a.mVit1 = mM1.vertexIt(mit.first.first);
        a.mVit2 = mM2.vertexIt(mit.first.second);
        mVertexAttributes.push_back(a);
        mPoints.push_back(a.mP);
    }
    for (auto mit : mMapVE) {
        Attributes a;
        a.mP    = mit.second;
        a.mPred = IF_VERTEX_EDGE;
        a.mVit1 = mM1.vertexIt(mit.first.first);
        a.mEit2 = mM2.edgeIt  (mit.first.second);
        mVertexAttributes.push_back(a);
        mPoints.push_back(a.mP);
    }
    for (auto mit : mMapVF) {
        Attributes a;
        a.mP    = mit.second;
        a.mPred = IF_VERTEX_FACE;
        a.mVit1 = mM1.vertexIt(mit.first.first);
        a.mFit2 = mM2.faceIt(mit.first.second);
        mVertexAttributes.push_back(a);
        mPoints.push_back(a.mP);
    }
    for (auto mit : mMapEV) {
        Attributes a;
        a.mP    = mit.second;
        a.mPred = IF_EDGE_VERTEX;
        a.mEit1 = mM1.edgeIt  (mit.first.first);
        a.mVit2 = mM2.vertexIt(mit.first.second);
        mVertexAttributes.push_back(a);
        mPoints.push_back(a.mP);
    }
    for (auto mit : mMapEE) {
        Attributes a;
        a.mP    = mit.second;
        a.mPred = IF_EDGE_EDGE;
        a.mEit1 = mM1.edgeIt(mit.first.first);
        a.mEit2 = mM2.edgeIt(mit.first.second);
        mVertexAttributes.push_back(a);
        mPoints.push_back(a.mP);
    }
    for (auto mit : mMapEF) {
        Attributes a;
        a.mP    = mit.second;
        a.mPred = IF_EDGE_FACE;
        a.mEit1 = mM1.edgeIt(mit.first.first);
        a.mFit2 = mM2.faceIt(mit.first.second);
        mVertexAttributes.push_back(a);
        mPoints.push_back(a.mP);
    }
    for (auto mit : mMapFV) {
        Attributes a;
        a.mP    = mit.second;
        a.mPred = IF_FACE_VERTEX;
        a.mFit1 = mM1.faceIt  (mit.first.first);
        a.mVit2 = mM2.vertexIt(mit.first.second);
        mVertexAttributes.push_back(a);
        mPoints.push_back(a.mP);
    }
    for (auto mit : mMapFE) {
        Attributes a;
        a.mP    = mit.second;
        a.mPred = IF_FACE_EDGE;
        a.mFit1 = mM1.faceIt(mit.first.first);
        a.mEit2 = mM2.edgeIt(mit.first.second);
        mVertexAttributes.push_back(a);
        mPoints.push_back(a.mP);
    }
    for (auto mit : mMapVI) {
        Attributes a;
        a.mP    = mit.second;
        a.mPred = IF_VERTEX_INTERIOR;
        a.mVit1 = mM1.vertexIt(mit.first);
        mVertexAttributes.push_back(a);
        mPoints.push_back(a.mP);
    }
    for (auto mit : mMapIV) {
        Attributes a;
        a.mP    = mit.second;
        a.mPred = IF_INTERIOR_VERTEX;
        a.mVit2 = mM2.vertexIt(mit.first);
        mVertexAttributes.push_back(a);
        mPoints.push_back(a.mP);
    }

    mMapVV.clear();
    mMapVE.clear();
    mMapVF.clear();
    mMapEV.clear();
    mMapEE.clear();
    mMapEF.clear();
    mMapFV.clear();
    mMapFE.clear();
    mMapVI.clear();
    mMapIV.clear();
}


void IntersectionFinder::findDimension()
{
    mDimension = mPoints.size() - 1;
    if (mPoints.size() > 3) {
        mDimension = findDimensionByPCA();
    }
}


long IntersectionFinder::findDimensionByPCA()
{
    mEigenMatrix = findPrincipalComponents(mPoints, mSpread, mMean);
    if (fabs(mSpread.z()) <= mEpsilonZeroPCA) {

        if (fabs(mSpread.y()) <= mEpsilonZeroPCA) {

            if (fabs(mSpread.x()) <= mEpsilonZeroPCA) {
                return 0;
            }
            return 1;
        }
        return 2;
    }
    return 3;
}


void IntersectionFinder::processZeroDimensionalIntersection()
{
    if (mVertexAttributes.size() > 0) {
        mVertexAttributes0D1D2D.push_back(mVertexAttributes[0]);
    }
    mVertexAttributes.clear();
}


IntersectionFinder::Attributes& IntersectionFinder::vertexAttributes0D()
{
     return mVertexAttributes0D1D2D[0];
}


void IntersectionFinder::processOneDimensionalIntersection()
{
    if (mVertexAttributes.size() == 2) {

        mVertexAttributes0D1D2D.push_back(mVertexAttributes[0]);
        mVertexAttributes0D1D2D.push_back(mVertexAttributes[1]);
    }
    else {

        const auto axis = mEigenMatrix.col(1);
        long minInd = 0;
        auto minDot = axis.dot(mPoints[0]);
        long maxInd = 0;
        auto maxDot = minDot;
        for (long i = 1; i < mPoints.size(); i++) {
            const auto curP   = mPoints[i];
            const auto curDot = axis.dot(curP);
            if (curDot > maxDot) {
                maxDot = curDot;
                maxInd = i;
            }
            else if (curDot < minDot) {
                minDot = curDot;
                minInd = i;
            }
        }

        mVertexAttributes0D1D2D.push_back(mVertexAttributes[minInd]);
        mVertexAttributes0D1D2D.push_back(mVertexAttributes[maxInd]);

    }        

    constructEdgeAttributes1D2D();
    mVertexAttributes.clear();
}


IntersectionFinder::Attributes &IntersectionFinder::vertexAttributes1D(long index)
{
    return mVertexAttributes0D1D2D[std::max(std::min(index,1L),0L)];
}


IntersectionFinder::Attributes &IntersectionFinder::edgeAttributes1D()
{
    return mEdgeAttributes1D2D[0];
}


void IntersectionFinder::constructEdgeAttributes1D2D()
{
    for (long i1 = 0; i1 < mVertexAttributes0D1D2D.size(); i1++) {

        long i2 = (i1 + 1) % mVertexAttributes0D1D2D.size();

        auto& srcAttr = mVertexAttributes0D1D2D[i1];
        auto& dstAttr = mVertexAttributes0D1D2D[i2];

        auto vertexPreds1 = extractEdgePredicatesSide1(srcAttr, dstAttr);
        auto vertexPreds2 = extractEdgePredicatesSide2(srcAttr, dstAttr);

        EdgeIt         eit1,      eit2;
        FaceIt         fit1,      fit2;
        enum predicate edgePred1, edgePred2;

        processEdgePredicatesSide1(
                       srcAttr, dstAttr, vertexPreds1, eit1, fit1, edgePred1 );

        processEdgePredicatesSide2(
                       srcAttr, dstAttr, vertexPreds2, eit2, fit2, edgePred2 );

        Attributes attr;
        switch (edgePred1) {
          case IF_INTERIOR:
            switch (edgePred2) {
              case IF_INTERIOR:
                attr.mPred = IF_INTERIOR_INTERIOR;
                break;
              case IF_EDGE:
                attr.mPred = IF_INTERIOR_EDGE;
                attr.mEit2 = eit2;
                break;
              case IF_FACE:
                attr.mPred = IF_INTERIOR_FACE;
                attr.mFit2 = fit2;
                break;
              default:
                break;
            }

            break;
          case IF_EDGE:
            attr.mEit1 = eit1;
            switch (edgePred2) {
              case IF_INTERIOR:
                attr.mPred = IF_EDGE_INTERIOR;
                break;
              case IF_EDGE:
                attr.mPred = IF_EDGE_EDGE;
                attr.mEit2 = eit2;
                break;
              case IF_FACE:
                attr.mPred = IF_EDGE_FACE;
                attr.mFit2 = fit2;
                break;
              default:
                break;
            }
            break;
          case IF_FACE:
            attr.mFit1 = fit1;
            switch (edgePred2) {
              case IF_INTERIOR:
                attr.mPred = IF_FACE_INTERIOR;
                break;
              case IF_EDGE:
                attr.mPred = IF_FACE_EDGE;
                attr.mEit2 = eit2;
                break;
              case IF_FACE:
                attr.mPred = IF_FACE_FACE;
                attr.mFit2 = fit2;
                break;
              default:
                break;
            }
            break;
          default:
            break;
        }

        mEdgeAttributes1D2D.push_back(attr);

        if (mVertexAttributes0D1D2D.size()==2) {
            return;// If the size is 2, there is only one edge attrib needed.
        }        
    }
}


void IntersectionFinder::processEdgePredicatesSide1(
    Attributes&     srcAttr,
    Attributes&     dstAttr,
    enum predicate  predIn,
    EdgeIt&         eitOut,
    FaceIt&         fitOut,
    enum predicate& predOut
) {

    predOut = IF_INTERIOR;

    switch (predIn) {
      case IF_VERTEX_VERTEX:
        {
            auto eit = mM1.findEdge(srcAttr.mVit1, dstAttr.mVit1);
            if (eit != mM1.edges().second) {
                predOut = IF_EDGE;
                eitOut  = eit;
            }
            else {
                auto fit = mM1.findFace(srcAttr.mVit1, dstAttr.mVit1);
                if (fit != mM1.faces().second) {
                    predOut = IF_FACE;
                    fitOut  = fit;
                }
            }
        }
        break;

      case IF_VERTEX_EDGE: 
        {
            auto he = (*(dstAttr.mEit1))->he1();
            auto src  = (*he)->src();
            auto dst  = (*he)->dst();
            if (src==srcAttr.mVit1) {
                predOut = IF_EDGE;
                eitOut  = dstAttr.mEit1;
            }
            else if  (dst==srcAttr.mVit1) {
                predOut = IF_EDGE;
                eitOut  = dstAttr.mEit1;
            }
            else {
                auto fit = mM1.findFace(srcAttr.mVit1, src, dst);
                if (fit != mM1.faces().second) {
                    predOut = IF_FACE;
                    fitOut  = fit;
                }
            }
        }
        break;

      case IF_VERTEX_FACE:
        {
            for (auto& he : (*(dstAttr.mFit1))->halfEdges()) {
                auto src = (*he)->src();
                if (src == srcAttr.mVit1) {
                    predOut = IF_FACE;
                    fitOut  = dstAttr.mFit1;
                    break;
                }
            }       
        }
        break;

      case IF_EDGE_VERTEX:
        {
            auto he = (*(srcAttr.mEit1))->he1();
            auto src  = (*he)->src();
            auto dst  = (*he)->dst();
            if (src==dstAttr.mVit1) {
                predOut = IF_EDGE;
                eitOut  = srcAttr.mEit1;
            }
            else if  (dst==dstAttr.mVit1) {
                predOut = IF_EDGE;
                eitOut  = srcAttr.mEit1;
            }
            else {
                auto fit = mM1.findFace(dstAttr.mVit1, src, dst);
                if (fit != mM1.faces().second) {
                    predOut = IF_FACE;
                    fitOut  = fit;
                }
            }
        }
        break;

      case IF_FACE_VERTEX:
        {
            for (auto& he : (*(srcAttr.mFit1))->halfEdges()) {
                auto src = (*he)->src();
                if (src == dstAttr.mVit1) {
                    predOut = IF_FACE;
                    fitOut  = srcAttr.mFit1;
                    break;
                }
            }
        }       
        break;

      case IF_EDGE_EDGE:
        {
            if (srcAttr.mEit1 == dstAttr.mEit1) {
                predOut = IF_EDGE;
                eitOut  = srcAttr.mEit1;
            }
            else {
                auto f11 = (*((*(srcAttr.mEit1))->he1()))->face();
                auto f12 = (*((*(srcAttr.mEit1))->he2()))->face();
                auto f21 = (*((*(dstAttr.mEit1))->he1()))->face();
                auto f22 = (*((*(dstAttr.mEit1))->he2()))->face();
                if (f11 == f21 || f11 == f22) {
                    predOut = IF_FACE;
                    fitOut  = f11;
                }
                else if (f12 == f21 || f12 == f22) {
                    predOut = IF_FACE;
                    fitOut  = f12;
                }
            }
        }
        break;

      case IF_EDGE_FACE:
        {
            auto he1  = (*(srcAttr.mEit1))->he1();
            auto fit1 = (*he1)->face();
            auto he2  = (*(srcAttr.mEit1))->he2();
            auto fit2 = (*he2)->face();
            if (fit1 == dstAttr.mFit1) {
                predOut = IF_FACE;
                fitOut  = dstAttr.mFit1;
            }
            else if (fit2 == dstAttr.mFit1) {
                predOut = IF_FACE;
                fitOut  = dstAttr.mFit1;
            }
        }
        break;

      case IF_FACE_EDGE:
        {
            auto he1  = (*(dstAttr.mEit1))->he1();
            auto fit1 = (*he1)->face();
            auto he2  = (*(dstAttr.mEit1))->he2();
            auto fit2 = (*he2)->face();
            if (fit1 == srcAttr.mFit1) {
                predOut = IF_FACE;
                fitOut  = srcAttr.mFit1;
            }
            else if (fit2 == srcAttr.mFit1) {
                predOut = IF_FACE;
                fitOut  = srcAttr.mFit1;
            }
        }
        break;

      case IF_FACE_FACE:
        {
            if (srcAttr.mFit1 == dstAttr.mFit1) {
                predOut = IF_FACE;
                fitOut  = srcAttr.mFit1;
            }
        }
        break;
      case IF_EDGE_INTERIOR:
      case IF_FACE_INTERIOR:
      case IF_VERTEX_INTERIOR:
      case IF_INTERIOR_VERTEX:
      case IF_INTERIOR_EDGE:
      case IF_INTERIOR_FACE:
      case IF_INTERIOR_INTERIOR:
      default:
        break;
    }
}


void IntersectionFinder::processEdgePredicatesSide2(
    Attributes&     srcAttr,
    Attributes&     dstAttr,
    enum predicate  predIn,
    EdgeIt&         eitOut,
    FaceIt&         fitOut,
    enum predicate& predOut
) {

    predOut = IF_INTERIOR;

    switch (predIn) {
      case IF_VERTEX_VERTEX:
        {
            auto eit = mM2.findEdge(srcAttr.mVit2, dstAttr.mVit2);
            if (eit != mM2.edges().second) {
                predOut = IF_EDGE;
                eitOut  = eit;
            }
            else {
                auto fit = mM2.findFace(srcAttr.mVit2, dstAttr.mVit2);
                if (fit != mM2.faces().second) {
                    predOut = IF_FACE;
                    fitOut  = fit;
                }
            }
        }
        break;

      case IF_VERTEX_EDGE: 
        {
            auto he = (*(dstAttr.mEit2))->he1();
            auto src  = (*he)->src();
            auto dst  = (*he)->dst();
            if (src==srcAttr.mVit2) {
                predOut = IF_EDGE;
                eitOut  = dstAttr.mEit2;
            }
            else if  (dst==srcAttr.mVit2) {
                predOut = IF_EDGE;
                eitOut  = dstAttr.mEit2;
            }
            else {
                auto fit = mM2.findFace(srcAttr.mVit2, src, dst);
                if (fit != mM2.faces().second) {
                    predOut = IF_FACE;
                    fitOut  = fit;
                }
            }
        }
        break;

      case IF_VERTEX_FACE:
        {
            for (auto& he : (*(dstAttr.mFit2))->halfEdges()) {
                auto src = (*he)->src();
                if (src == srcAttr.mVit2) {
                    predOut = IF_FACE;
                    fitOut  = dstAttr.mFit2;
                    break;
                }
            }
        }       
        break;

      case IF_EDGE_VERTEX:
        {
            auto he = (*(srcAttr.mEit2))->he1();
            auto src  = (*he)->src();
            auto dst  = (*he)->dst();
            if (src==dstAttr.mVit2) {
                predOut = IF_EDGE;
                eitOut  = srcAttr.mEit2;
            }
            else if  (dst==dstAttr.mVit2) {
                predOut = IF_EDGE;
                eitOut  = srcAttr.mEit2;
            }
            else {
                auto fit = mM2.findFace(dstAttr.mVit2, src, dst);
                if (fit != mM2.faces().second) {
                    predOut = IF_FACE;
                    fitOut  = fit;
                }
            }
        }
        break;

      case IF_FACE_VERTEX:
        {
            for (auto& he : (*(srcAttr.mFit2))->halfEdges()) {
                auto src = (*he)->src();
                if (src == dstAttr.mVit2) {
                    predOut = IF_FACE;
                    fitOut  = srcAttr.mFit2;
                    break;
                }
            }       
        }
        break;

      case IF_EDGE_EDGE:
        {
            if (srcAttr.mEit2 == dstAttr.mEit2) {
                predOut = IF_EDGE;
                eitOut  = srcAttr.mEit2;
            }
            else {
                auto f11 = (*((*(srcAttr.mEit2))->he1()))->face();
                auto f12 = (*((*(srcAttr.mEit2))->he2()))->face();
                auto f21 = (*((*(dstAttr.mEit2))->he1()))->face();
                auto f22 = (*((*(dstAttr.mEit2))->he2()))->face();
                if (f11 == f21 || f11 == f22) {
                    predOut = IF_FACE;
                    fitOut  = f11;
                }
                else if (f12 == f21 || f12 == f22) {
                    predOut = IF_FACE;
                    fitOut  = f12;
                }
            }
        }
        break;

      case IF_EDGE_FACE:
        {
            auto he1  = (*(srcAttr.mEit2))->he1();
            auto fit1 = (*he1)->face();
            auto he2  = (*(srcAttr.mEit2))->he2();
            auto fit2 = (*he2)->face();
            if (fit1 == dstAttr.mFit2) {
                predOut = IF_FACE;
                fitOut  = dstAttr.mFit2;
            }
            else if (fit2 == dstAttr.mFit2) {
                predOut = IF_FACE;
                fitOut  = dstAttr.mFit2;
            }
        }
        break;

      case IF_FACE_EDGE:
        {
            auto he1  = (*(dstAttr.mEit2))->he1();
            auto fit1 = (*he1)->face();
            auto he2  = (*(dstAttr.mEit2))->he2();
            auto fit2 = (*he2)->face();
            if (fit1 == srcAttr.mFit2) {
                predOut = IF_FACE;
                fitOut  = srcAttr.mFit2;
            }
            else if (fit2 == srcAttr.mFit2) {
                predOut = IF_FACE;
                fitOut  = srcAttr.mFit2;
            }
        }
        break;

      case IF_FACE_FACE:
        {
            if (srcAttr.mFit2 == dstAttr.mFit2) {
                predOut = IF_FACE;
                fitOut  = srcAttr.mFit2;
            }
        }
        break;

      case IF_EDGE_INTERIOR:
      case IF_FACE_INTERIOR:
      case IF_VERTEX_INTERIOR:
      case IF_INTERIOR_VERTEX:
      case IF_INTERIOR_EDGE:
      case IF_INTERIOR_FACE:
      case IF_INTERIOR_INTERIOR:
      default:
        break;
    }
}


void IntersectionFinder::processTwoDimensionalIntersection()
{
    vector<Vec2> rotatedPoints;
    Mat3x3 rotMat = mEigenMatrix.transpose();
    for (auto& p : mPoints) {
        auto rP = rotMat * p;
        rotatedPoints.push_back(Vec2(rP.x(), rP.y()));
    }
    auto resultIndices = findConvexHull2D(rotatedPoints);

    constructVertexAttributes1D2D(resultIndices);

    constructEdgeAttributes1D2D();

    constructFaceAttributes2D();

    mVertexAttributes.clear();
}


void IntersectionFinder::constructVertexAttributes1D2D(vector<long>& indices)
{
    for (auto i : indices) {
        mVertexAttributes0D1D2D.push_back(mVertexAttributes[i]);
    }               
}


void IntersectionFinder::constructFaceAttributes2D()
{
    FaceIt fitOut1, fitOut2;
    auto res1 = findUnieuqFaceSide12D(fitOut1);
    if (res1) {
        res1 = checkFaceIncidenceSide12D(fitOut1);
    }

    auto res2 = findUnieuqFaceSide22D(fitOut2);
    if (res2) {
        res2 = checkFaceIncidenceSide22D(fitOut2);
    }

    auto& attr = mFaceAttributes2D;
    if (res1) {
        attr.mFit1 = fitOut1;
        if (res2) {
            attr.mFit2 = fitOut2;
            attr.mPred = IF_FACE_FACE;
        }
        else {
            attr.mPred = IF_FACE_INTERIOR;
        }
    }
    else {
        if (res2) {
            attr.mFit2 = fitOut2;
            attr.mPred = IF_INTERIOR_FACE;
        }
        else {
            attr.mPred = IF_INTERIOR_INTERIOR;
        }
    }
}


enum predicate IntersectionFinder::typeSide1(enum predicate p)
{
    switch (p) {

      case IF_VERTEX_VERTEX:
      case IF_VERTEX_EDGE:
      case IF_VERTEX_FACE:
      case IF_VERTEX_INTERIOR:
        return IF_VERTEX;

      case IF_EDGE_VERTEX:
      case IF_EDGE_EDGE:
      case IF_EDGE_FACE:
      case IF_EDGE_INTERIOR:
        return IF_EDGE;

      case IF_FACE_VERTEX:
      case IF_FACE_EDGE:
      case IF_FACE_FACE:
      case IF_FACE_INTERIOR:
        return IF_FACE;

      case IF_INTERIOR_VERTEX:
      case IF_INTERIOR_EDGE:
      case IF_INTERIOR_FACE:
      case IF_INTERIOR_INTERIOR:
      default:
        return IF_INTERIOR;
    }
}


enum predicate IntersectionFinder::typeSide2(enum predicate p)
{
    switch (p) {

      case IF_VERTEX_VERTEX:
      case IF_EDGE_VERTEX:
      case IF_FACE_VERTEX:
      case IF_INTERIOR_VERTEX:
        return IF_VERTEX;

      case IF_VERTEX_EDGE:
      case IF_EDGE_EDGE:
      case IF_FACE_EDGE:
      case IF_INTERIOR_EDGE:
        return IF_EDGE;

      case IF_VERTEX_FACE:
      case IF_EDGE_FACE:
      case IF_FACE_FACE:
      case IF_INTERIOR_FACE:
        return IF_FACE;

      case IF_VERTEX_INTERIOR:
      case IF_EDGE_INTERIOR:
      case IF_FACE_INTERIOR:
      case IF_INTERIOR_INTERIOR:
      default:
        return IF_INTERIOR;
    }
}


bool IntersectionFinder::checkFaceIncidenceSide12D(FaceIt fit)
{
    for (auto& vAttrib : mVertexAttributes0D1D2D) {

        auto t = typeSide1(vAttrib.mPred);

        if (t == IF_VERTEX) {
            if (!(*fit)->isIncident(vAttrib.mVit1)) {
                return false;
            }
        }
        if ( t == IF_EDGE) {
            if (!(*fit)->isIncident(vAttrib.mEit1)) {
                return false;
            }
        }
    }

    for (auto& eAttrib : mEdgeAttributes1D2D) {

        auto t = typeSide1(eAttrib.mPred);

        if ( t == IF_EDGE) {
            if (!(*fit)->isIncident(eAttrib.mEit1)) {
                return false;
            }
        }
    }
    return true;
}


bool IntersectionFinder::checkFaceIncidenceSide22D(FaceIt fit)
{
    for (auto& vAttrib : mVertexAttributes0D1D2D) {

        auto t = typeSide2(vAttrib.mPred);

        if (t == IF_VERTEX) {
            if (!(*fit)->isIncident(vAttrib.mVit2)) {
                return false;
            }
        }
        if ( t == IF_EDGE) {
            if (!(*fit)->isIncident(vAttrib.mEit2)) {
                return false;
            }
        }
    }

    for (auto& eAttrib : mEdgeAttributes1D2D) {

        auto t = typeSide2(eAttrib.mPred);

        if ( t == IF_EDGE) {
            if (!(*fit)->isIncident(eAttrib.mEit2)) {
                return false;
            }
        }
    }
    return true;
}


bool IntersectionFinder::findUnieuqFaceSide12D(FaceIt& fit)
{
    // Find the unique face
    long     faceCnt = 0;
    for (auto& vAttrib : mVertexAttributes0D1D2D) {

        auto t = typeSide1(vAttrib.mPred);

        if ( t == IF_INTERIOR) {

            return false;
        }
        if ( t == IF_FACE) {
            if (faceCnt == 0) {

                fit = vAttrib.mFit1;
                faceCnt++;
            }
            else if (fit != vAttrib.mFit1) {

                return false;
            }
        }
    }

    for (auto& eAttrib : mEdgeAttributes1D2D) {

        auto t = typeSide1(eAttrib.mPred);

        if ( t == IF_INTERIOR) {
            return false;
        }
        if ( t == IF_FACE) {
            if (faceCnt == 0) {

                fit = eAttrib.mFit1;
                faceCnt++;
            }
            else if (fit != eAttrib.mFit1) {
                return false;
            }
        }
    }

    if (faceCnt > 0) {
        return true;
    }

    long     vertexCnt = 0;
    VertexIt vit1, vit2, vit3;

    for (auto& vAttrib : mVertexAttributes0D1D2D) {

        auto t = typeSide1(vAttrib.mPred);
        if ( t == IF_VERTEX ) {
            if (vertexCnt==0) {
                vit1 = vAttrib.mVit1;
                vertexCnt++;
            }
            else if (vertexCnt==1 && vit1 != vAttrib.mVit1) {
                vit2 = vAttrib.mVit1;
                vertexCnt++;
            }
            else if (vertexCnt==2 && vit1 != vAttrib.mVit1 && 
                                     vit2 != vAttrib.mVit1     ) {
                vit3 = vAttrib.mVit1;
                vertexCnt++;
            }
        }

        if (vertexCnt==3) {
            break;
        }
    }

    if (vertexCnt<3) {

        for (auto& eAttrib : mEdgeAttributes1D2D) {

            auto t = typeSide1(eAttrib.mPred);

            if ( t == IF_EDGE) {

                auto heit1 = (*(eAttrib.mEit1))->he1();
                auto src   = (*heit1)->src();
                auto dst   = (*heit1)->dst();
                
                if (vertexCnt==0) {
                    vit1 = src;
                    vertexCnt++;
                }
                else if (vertexCnt==1 && vit1 != src) {
                    vit2 = src;
                    vertexCnt++;
                }
                else if (vertexCnt==2 && vit1 != src && vit2 != src) {
                    vit3 = src;
                    vertexCnt++;
                }

                if (vertexCnt==0) {
                    vit1 = dst;
                    vertexCnt++;
                }
                else if (vertexCnt==1 && vit1 != dst) {
                    vit2 = dst;
                    vertexCnt++;
                }
                else if (vertexCnt==2 && vit1 != dst && vit2 != dst) {
                    vit3 = dst;
                    vertexCnt++;
                }
            }

            if (vertexCnt==3) {
                break;
            }
        }
    }

    if (vertexCnt==3) {

        fit = mM1.findFace(vit1, vit2, vit3);
        if (fit != mM1.faces().second) {
            return true;
        }    

    }
    return false;
}


bool IntersectionFinder::findUnieuqFaceSide22D(FaceIt& fit)
{
    // Find the unique face
    long     faceCnt = 0;
    for (auto& vAttrib : mVertexAttributes0D1D2D) {

        auto t = typeSide2(vAttrib.mPred);

        if ( t == IF_INTERIOR) {
            return false;
        }
        if ( t == IF_FACE) {
            if (faceCnt == 0) {
                fit = vAttrib.mFit2;
                faceCnt++;
            }
            else if (fit != vAttrib.mFit2) {
                return false;
            }
        }
    }

    for (auto& eAttrib : mEdgeAttributes1D2D) {

        auto t = typeSide2(eAttrib.mPred);

        if ( t == IF_INTERIOR) {
            return false;
        }
        if ( t == IF_FACE) {
            if (faceCnt == 0) {
                fit = eAttrib.mFit2;
                faceCnt++;
            }
            else if (fit != eAttrib.mFit2) {
                return false;
            }
        }
    }

    if (faceCnt > 0) {
        return true;
    }

    long     vertexCnt = 0;
    VertexIt vit1, vit2, vit3;
    for (auto& vAttrib : mVertexAttributes0D1D2D) {

        auto t = typeSide2(vAttrib.mPred);

        if ( t == IF_VERTEX ) {
            if (vertexCnt==0) {
                vit1 = vAttrib.mVit2;
                vertexCnt++;
            }
            else if (vertexCnt==1 && vit1 != vAttrib.mVit2) {
                vit2 = vAttrib.mVit2;
                vertexCnt++;
            }
            else if (vertexCnt==2 && vit1 != vAttrib.mVit2 && 
                                     vit2 != vAttrib.mVit2     ) {
                vit3 = vAttrib.mVit2;
                vertexCnt++;
            }
        }

        if (vertexCnt==3) {
            break;
        }
    }

    if (vertexCnt<3) {

        for (auto& eAttrib : mEdgeAttributes1D2D) {

            auto t = typeSide2(eAttrib.mPred);

            if ( t == IF_EDGE) {

                auto heit1 = (*(eAttrib.mEit2))->he1();
                auto src   = (*heit1)->src();
                auto dst   = (*heit1)->dst();
                
                if (vertexCnt==0) {
                    vit1 = src;
                    vertexCnt++;
                }
                else if (vertexCnt==1 && vit1 != src) {
                    vit2 = src;
                    vertexCnt++;
                }
                else if (vertexCnt==2 && vit1 != src && vit2 != src) {
                    vit3 = src;
                    vertexCnt++;
                }

                if (vertexCnt==0) {
                    vit1 = dst;
                    vertexCnt++;
                }
                else if (vertexCnt==1 && vit1 != dst) {
                    vit2 = dst;
                    vertexCnt++;
                }
                else if (vertexCnt==2 && vit1 != dst && vit2 != dst) {
                    vit3 = dst;
                    vertexCnt++;
                }
            }

            if (vertexCnt==3) {
                break;
            }
        }
    }

    if (vertexCnt==3) {

        fit = mM2.findFace(vit1, vit2, vit3);
        if (fit != mM2.faces().second) {
            return true;
        }    

    }
    return false;
}


long IntersectionFinder::size2D() const {
    return mVertexAttributes0D1D2D.size();
}


IntersectionFinder::Attributes &IntersectionFinder::vertexAttributes2D(long index)
{
    return mVertexAttributes0D1D2D[index];
}


IntersectionFinder::Attributes &IntersectionFinder::edgeAttributes2D(long index)
{
    return mEdgeAttributes1D2D[index];
}


IntersectionFinder::Attributes &IntersectionFinder::faceAttributes2D()
{
    return mFaceAttributes2D;
}


IntersectionFinder::Attributes &IntersectionFinder::vertexAttributes3D(VertexIt vit)
{
    return mVertexAttributes3D[(*vit)->id()];
}


IntersectionFinder::Attributes &IntersectionFinder::edgeAttributes3D(EdgeIt eit)
{
    return mEdgeAttributes3D[(*eit)->id()];
}


IntersectionFinder::Attributes &IntersectionFinder::faceAttributes3D(FaceIt fit)
{
    return mFaceAttributes3D[(*fit)->id()];
}


long IntersectionFinder::dimension() const {
    return mDimension;
}


Manifold& IntersectionFinder::hull3D()
{ 
    return mResultMani;
}


void IntersectionFinder::processThreeDimensionalIntersection()
{
    enum predicate pred;
    mResultMani.findConvexHull(mPoints, pred, mEpsilonZero);
    constructVertexAttributes3D();
    constructEdgeAttributes3D();
    constructFaceAttributes3D();
}


void IntersectionFinder::constructVertexAttributes3D()
{
    auto vPair = mResultMani.vertices();
    for (auto vit = vPair.first; vit != vPair.second; vit++) {
        mVertexAttributes3D[(*vit)->id()] = mVertexAttributes[(*vit)->id()];
    }
    mVertexAttributes.clear();

}


enum predicate IntersectionFinder::extractEdgePredicatesSide1(
    Attributes& srcAttr,
    Attributes& dstAttr
) {
    enum predicate pred = NONE;

    switch (srcAttr.mPred) {

      case IF_VERTEX_VERTEX: case IF_VERTEX_EDGE:
      case IF_VERTEX_FACE:   case IF_VERTEX_INTERIOR:
     
        switch (dstAttr.mPred) {

          case IF_VERTEX_VERTEX: case IF_VERTEX_EDGE: 
          case IF_VERTEX_FACE:   case IF_VERTEX_INTERIOR:

            pred = IF_VERTEX_VERTEX;

            break;
          case IF_EDGE_VERTEX:   case IF_EDGE_EDGE:   case IF_EDGE_FACE:
            pred = IF_VERTEX_EDGE;
            break;

          case IF_FACE_VERTEX:   case IF_FACE_EDGE:
            pred = IF_VERTEX_FACE;
            break;

          case IF_INTERIOR_VERTEX: default:
            pred = IF_VERTEX_INTERIOR;
            break;
        }
        break;

      case IF_EDGE_VERTEX:   case IF_EDGE_EDGE:   case IF_EDGE_FACE:

        switch (dstAttr.mPred) {

          case IF_VERTEX_VERTEX: case IF_VERTEX_EDGE: 
          case IF_VERTEX_FACE:   case IF_VERTEX_INTERIOR:

            pred = IF_EDGE_VERTEX;
            break;

          case IF_EDGE_VERTEX:   case IF_EDGE_EDGE:   case IF_EDGE_FACE:

            pred = IF_EDGE_EDGE;
            break;
          case IF_FACE_VERTEX:   case IF_FACE_EDGE:

            pred = IF_EDGE_FACE;
            break;

          case IF_INTERIOR_VERTEX: default:
            pred = IF_EDGE_INTERIOR;
            break;
        }
        break;

      case IF_FACE_VERTEX:   case IF_FACE_EDGE:

        switch (dstAttr.mPred) {

          case IF_VERTEX_VERTEX: case IF_VERTEX_EDGE: 
          case IF_VERTEX_FACE:   case IF_VERTEX_INTERIOR:

            pred = IF_FACE_VERTEX;
            break;

          case IF_EDGE_VERTEX:   case IF_EDGE_EDGE:   case IF_EDGE_FACE:
            pred = IF_FACE_EDGE;
            break;

          case IF_FACE_VERTEX:   case IF_FACE_EDGE:
            pred = IF_FACE_FACE;
            break;

          case IF_INTERIOR_VERTEX: default:
            pred = IF_FACE_INTERIOR;
            break;
        }
        break;

      case IF_INTERIOR_VERTEX: default:

        switch (dstAttr.mPred) {

          case IF_VERTEX_VERTEX: case IF_VERTEX_EDGE: 
          case IF_VERTEX_FACE:   case IF_VERTEX_INTERIOR:
            pred = IF_INTERIOR_VERTEX;
            break;

          case IF_EDGE_VERTEX:   case IF_EDGE_EDGE:   case IF_EDGE_FACE:
            pred = IF_INTERIOR_EDGE;
            break;

          case IF_FACE_VERTEX:   case IF_FACE_EDGE:
            pred = IF_INTERIOR_FACE;
            break;

          case IF_INTERIOR_VERTEX: default:
            pred = IF_INTERIOR_INTERIOR;
            break;
        }
        break;
    }
    return pred;
}


enum predicate IntersectionFinder::extractEdgePredicatesSide2(
    Attributes& srcAttr,
    Attributes& dstAttr
) {

    enum predicate pred = NONE;

    switch (srcAttr.mPred) {

      case IF_VERTEX_VERTEX: case IF_EDGE_VERTEX:
      case IF_FACE_VERTEX:   case IF_INTERIOR_VERTEX:
     
        switch (dstAttr.mPred) {

          case IF_VERTEX_VERTEX: case IF_EDGE_VERTEX: 
          case IF_FACE_VERTEX:   case IF_INTERIOR_VERTEX:
            pred = IF_VERTEX_VERTEX;
            break;

          case IF_VERTEX_EDGE:   case IF_EDGE_EDGE:   case IF_FACE_EDGE:
            pred = IF_VERTEX_EDGE;
            break;

          case IF_VERTEX_FACE:   case IF_EDGE_FACE:
            pred = IF_VERTEX_FACE;
            break;

          case IF_VERTEX_INTERIOR: default:
            pred = IF_VERTEX_INTERIOR;
            break;
        }
        break;

      case IF_VERTEX_EDGE:   case IF_EDGE_EDGE:   case IF_FACE_EDGE:

        switch (dstAttr.mPred) {

          case IF_VERTEX_VERTEX: case IF_EDGE_VERTEX: 
          case IF_FACE_VERTEX:   case IF_INTERIOR_VERTEX:
            pred = IF_EDGE_VERTEX;
            break;

          case IF_VERTEX_EDGE:   case IF_EDGE_EDGE:   case IF_FACE_EDGE:
            pred = IF_EDGE_EDGE;
            break;

          case IF_VERTEX_FACE:   case IF_EDGE_FACE:
            pred = IF_EDGE_FACE;
            break;

          case IF_VERTEX_INTERIOR: default:
            pred = IF_EDGE_INTERIOR;
            break;
        }
        break;

      case IF_VERTEX_FACE:   case IF_EDGE_FACE:

        switch (dstAttr.mPred) {

          case IF_VERTEX_VERTEX: case IF_EDGE_VERTEX: 
          case IF_FACE_VERTEX:   case IF_INTERIOR_VERTEX:
            pred = IF_FACE_VERTEX;
            break;

          case IF_VERTEX_EDGE:   case IF_EDGE_EDGE:   case IF_FACE_EDGE:
            pred = IF_FACE_EDGE;
            break;

          case IF_VERTEX_FACE:   case IF_EDGE_FACE:
            pred = IF_FACE_FACE;
            break;

          case IF_VERTEX_INTERIOR: default:
            pred = IF_FACE_INTERIOR;
            break;
        }
        break;

      case IF_VERTEX_INTERIOR: default:

        switch (dstAttr.mPred) {

          case IF_VERTEX_VERTEX: case IF_EDGE_VERTEX: 
          case IF_FACE_VERTEX:   case IF_INTERIOR_VERTEX:
            pred = IF_INTERIOR_VERTEX;
            break;

          case IF_VERTEX_EDGE:   case IF_EDGE_EDGE:   case IF_FACE_EDGE:
            pred = IF_INTERIOR_EDGE;
            break;

          case IF_VERTEX_FACE:   case IF_EDGE_FACE:
            pred = IF_INTERIOR_FACE;
            break;

          case IF_VERTEX_INTERIOR: default:
            pred = IF_INTERIOR_INTERIOR;
            break;
        }
        break;
    }
    return pred;
}


void IntersectionFinder::constructEdgeAttributes3D()
{
    auto ePair = mResultMani.edges();

    for (auto eit = ePair.first; eit != ePair.second; eit++) {

        auto heit     = (*eit)->he1();
        auto src      = (*heit)->src();
        auto dst      = (*heit)->dst();
        auto& srcAttr = mVertexAttributes[(*src)->id()];
        auto& dstAttr = mVertexAttributes[(*dst)->id()];

        auto vertexPreds1 = extractEdgePredicatesSide1(srcAttr, dstAttr);
        auto vertexPreds2 = extractEdgePredicatesSide2(srcAttr, dstAttr);

        EdgeIt         eit1,      eit2;
        FaceIt         fit1,      fit2;
        enum predicate edgePred1, edgePred2;

        processEdgePredicatesSide1(
                    srcAttr, dstAttr, vertexPreds1, eit1, fit1, edgePred1 );

        processEdgePredicatesSide2(
                    srcAttr, dstAttr, vertexPreds2, eit2, fit2, edgePred2 );

        Attributes attr;
        switch (edgePred1) {
          case IF_INTERIOR:
            switch (edgePred2) {
              case IF_INTERIOR:
                attr.mPred = IF_INTERIOR_INTERIOR;
                break;
              case IF_EDGE:
                attr.mPred = IF_INTERIOR_EDGE;
                attr.mEit2 = eit2;
                break;
              case IF_FACE:
                attr.mPred = IF_INTERIOR_FACE;
                attr.mFit2 = fit2;
                break;
              default:
                break;
            }

            break;
          case IF_EDGE:
            attr.mEit1 = eit1;
            switch (edgePred2) {
              case IF_INTERIOR:
                attr.mPred = IF_EDGE_INTERIOR;
                break;
              case IF_EDGE:
                attr.mPred = IF_EDGE_EDGE;
                attr.mEit2 = eit2;
                break;
              case IF_FACE:
                attr.mPred = IF_EDGE_FACE;
                attr.mFit2 = fit2;
                break;
              default:
                break;
            }
            break;
          case IF_FACE:
            attr.mFit1 = fit1;
            switch (edgePred2) {
              case IF_INTERIOR:
                attr.mPred = IF_FACE_INTERIOR;
                break;
              case IF_EDGE:
                attr.mPred = IF_FACE_EDGE;
                attr.mEit2 = eit2;
                break;
              case IF_FACE:
                attr.mPred = IF_FACE_FACE;
                attr.mFit2 = fit2;
                break;
              default:
                break;
            }
            break;
          default:
            break;
        }
        mEdgeAttributes3D[(*eit)->id()] = attr;

    }
}


void IntersectionFinder::constructFaceAttributes3D()
{
    auto fPair = mResultMani.faces();
    for (auto fit = fPair.first; fit != fPair.second; fit++) {

        FaceIt fitOut1, fitOut2;
        auto res1 = findUnieuqFaceSide13D(fit, fitOut1);
        auto res2 = findUnieuqFaceSide23D(fit, fitOut2);

        Attributes attr;
        if (res1) {
            attr.mFit1 = fitOut1;
            if (res2) {
                attr.mFit2 = fitOut2;
                attr.mPred = IF_FACE_FACE;
            }
            else {
                attr.mPred = IF_FACE_INTERIOR;
            }
        }
        else {
            if (res2) {
                attr.mFit2 = fitOut2;
                attr.mPred = IF_INTERIOR_FACE;
            }
            else {
                attr.mPred = IF_INTERIOR_INTERIOR;
            }
        }

        mFaceAttributes3D[(*fit)->id()] = attr;

    }
}


bool IntersectionFinder::findUnieuqFaceSide13D(
    FaceIt  fitIn,
    FaceIt& fitOut
) {
    // Find the unique face
    long     faceCnt = 0;

    for (auto heit : (*fitIn)->halfEdges()) {
        auto  src     = (*heit)->src();
        auto& vAttrib = mVertexAttributes[(*src)->id()];
        auto t        = typeSide1(vAttrib.mPred);
        if ( t == IF_INTERIOR) {
            return false;
        }
        if ( t == IF_FACE) {
            if (faceCnt == 0) {
                fitOut = vAttrib.mFit1;
                faceCnt++;
            }
            else if (fitOut != vAttrib.mFit1) {
                return false;
            }
        }
    }

    for (auto heit : (*fitIn)->halfEdges()) {
        auto  eit     = (*heit)->edge();
        auto& eAttrib = mEdgeAttributes3D[(*eit)->id()];
        auto t        = typeSide1(eAttrib.mPred);
        if ( t == IF_INTERIOR) {
            return false;
        }
        if ( t == IF_FACE) {
            if (faceCnt == 0) {
                fitOut = eAttrib.mFit1;
                faceCnt++;
            }
            else if (fitOut != eAttrib.mFit1) {
                return false;
            }
        }
    }

    if (faceCnt > 0) {
        return true;
    }

    long     vertexCnt = 0;
    VertexIt vit1, vit2, vit3;
    for (auto heit : (*fitIn)->halfEdges()) {

        auto  src     = (*heit)->src();
        auto& vAttrib = mVertexAttributes[(*src)->id()];
        auto  t       = typeSide1(vAttrib.mPred);

        if ( t == IF_VERTEX ) {
            if (vertexCnt==0) {
                vit1 = vAttrib.mVit1;
                vertexCnt++;
            }
            else if (vertexCnt==1 && vit1 != vAttrib.mVit1) {
                vit2 = vAttrib.mVit1;
                vertexCnt++;
            }
            else if (vertexCnt==2 && vit1 != vAttrib.mVit1 && 
                                     vit2 != vAttrib.mVit1     ) {
                vit3 = vAttrib.mVit1;
                vertexCnt++;
            }
        }

        if (vertexCnt==3) {
            break;
        }
    }

    if (vertexCnt<3) {

        for (auto heit : (*fitIn)->halfEdges()) {
            auto  eit     = (*heit)->edge();
            auto& eAttrib = mEdgeAttributes3D[(*eit)->id()];
            auto t        = typeSide1(eAttrib.mPred);
            if ( t == IF_EDGE) {
                auto heit1 = (*(eAttrib.mEit1))->he1();
                auto src   = (*heit1)->src();
                auto dst   = (*heit1)->dst();
                
                if (vertexCnt==0) {
                    vit1 = src;
                    vertexCnt++;
                }
                else if (vertexCnt==1 && vit1 != src) {
                    vit2 = src;
                    vertexCnt++;
                }
                else if (vertexCnt==2 && vit1 != src && vit2 != src) {
                    vit3 = src;
                    vertexCnt++;
                }

                if (vertexCnt==0) {
                    vit1 = dst;
                    vertexCnt++;
                }
                else if (vertexCnt==1 && vit1 != dst) {
                    vit2 = dst;
                    vertexCnt++;
                }
                else if (vertexCnt==2 && vit1 != dst && vit2 != dst) {
                    vit3 = dst;
                    vertexCnt++;
                }
            }

            if (vertexCnt==3) {
                break;
            }
        }
    }

    if (vertexCnt==3) {

        fitOut = mM1.findFace(vit1, vit2, vit3);
        if (fitOut != mM1.faces().second) {
            return true;
        }    

    }
    return false;
}


bool IntersectionFinder::findUnieuqFaceSide23D(
    FaceIt  fitIn,
    FaceIt& fitOut
) {
    // Find the unique face
    long     faceCnt = 0;

    for (auto heit : (*fitIn)->halfEdges()) {
        auto  src     = (*heit)->src();
        auto& vAttrib = mVertexAttributes[(*src)->id()];
        auto t        = typeSide2(vAttrib.mPred);
        if ( t == IF_INTERIOR) {
            return false;
        }
        if ( t == IF_FACE) {
            if (faceCnt == 0) {
                fitOut = vAttrib.mFit2;
                faceCnt++;
            }
            else if (fitOut != vAttrib.mFit2) {
                return false;
            }
        }
    }

    for (auto heit : (*fitIn)->halfEdges()) {
        auto  eit     = (*heit)->edge();
        auto& eAttrib = mEdgeAttributes3D[(*eit)->id()];
        auto t        = typeSide2(eAttrib.mPred);
        if ( t == IF_INTERIOR) {
            return false;
        }
        if ( t == IF_FACE) {
            if (faceCnt == 0) {
                fitOut = eAttrib.mFit2;
                faceCnt++;
            }
            else if (fitOut != eAttrib.mFit2) {
                return false;
            }
        }
    }

    if (faceCnt > 0) {
        return true;
    }

    long     vertexCnt = 0;
    VertexIt vit1, vit2, vit3;
    for (auto heit : (*fitIn)->halfEdges()) {

        auto  src     = (*heit)->src();
        auto& vAttrib = mVertexAttributes[(*src)->id()];
        auto  t       = typeSide2(vAttrib.mPred);

        if ( t == IF_VERTEX ) {
            if (vertexCnt==0) {
                vit1 = vAttrib.mVit2;
                vertexCnt++;
            }
            else if (vertexCnt==1 && vit1 != vAttrib.mVit2) {
                vit2 = vAttrib.mVit2;
                vertexCnt++;
            }
            else if (vertexCnt==2 && vit1 != vAttrib.mVit2 && 
                                     vit2 != vAttrib.mVit2     ) {
                vit3 = vAttrib.mVit2;
                vertexCnt++;
            }
        }

        if (vertexCnt==3) {
            break;
        }
    }

    if (vertexCnt<3) {

        for (auto heit : (*fitIn)->halfEdges()) {
            auto  eit     = (*heit)->edge();
            auto& eAttrib = mEdgeAttributes3D[(*eit)->id()];
            auto t        = typeSide2(eAttrib.mPred);

            if ( t == IF_EDGE) {

                auto heit1 = (*(eAttrib.mEit2))->he1();
                auto src   = (*heit1)->src();
                auto dst   = (*heit1)->dst();
                
                if (vertexCnt==0) {
                    vit1 = src;
                    vertexCnt++;
                }
                else if (vertexCnt==1 && vit1 != src) {
                    vit2 = src;
                    vertexCnt++;
                }
                else if (vertexCnt==2 && vit1 != src && vit2 != src) {
                    vit3 = src;
                    vertexCnt++;
                }

                if (vertexCnt==0) {
                    vit1 = dst;
                    vertexCnt++;
                }
                else if (vertexCnt==1 && vit1 != dst) {
                    vit2 = dst;
                    vertexCnt++;
                }
                else if (vertexCnt==2 && vit1 != dst && vit2 != dst) {
                    vit3 = dst;
                    vertexCnt++;
                }
            }

            if (vertexCnt==3) {
                break;
            }
        }
    }

    if (vertexCnt==3) {

        fitOut = mM2.findFace(vit1, vit2, vit3);
        if (fitOut != mM2.faces().second) {
            return true;
        }    
    }

    return false;
}


bool IntersectionFinder::checkFaceIncidenceSide13D(
    FaceIt fitPoly1,
    FaceIt fitIntsec
) {

    for (auto heit : (*fitPoly1)->halfEdges()) {

        auto  eit     = (*heit)->edge();
        auto& eAttrib = mEdgeAttributes3D[(*eit)->id()];
        auto  te      = typeSide1(eAttrib.mPred);

        if ( te == IF_EDGE) {
            if ((*fitIntsec)->isIncident(eAttrib.mEit1)) {
                return false;
            }
        }

        auto  src     = (*heit)->src();
        auto& vAttrib = mVertexAttributes[(*src)->id()];
        auto  tv      = typeSide1(vAttrib.mPred);

        if (tv == IF_VERTEX) {
            if (!(*fitIntsec)->isIncident(vAttrib.mVit1)) {
                return false;
            }
        }
        if ( tv == IF_EDGE) {
            if (!(*fitIntsec)->isIncident(vAttrib.mEit1)) {
                return false;
            }
        }
    }

    return true;
}


bool IntersectionFinder::checkFaceIncidenceSide23D(
    FaceIt fitPoly2,
    FaceIt fitIntsec
) {

    for (auto heit : (*fitPoly2)->halfEdges()) {

        auto  eit     = (*heit)->edge();
        auto& eAttrib = mEdgeAttributes3D[(*eit)->id()];
        auto te        = typeSide2(eAttrib.mPred);

        if ( te == IF_EDGE) {
            if ((*fitIntsec)->isIncident(eAttrib.mEit2)) {
                return false;
            }
        }

        auto  src     = (*heit)->src();
        auto& vAttrib = mVertexAttributes[(*src)->id()];
        auto tv        = typeSide2(vAttrib.mPred);

        if (tv == IF_VERTEX) {
            if (!(*fitIntsec)->isIncident(vAttrib.mVit2)) {
                return false;
            }
        }
        if ( tv == IF_EDGE) {
            if (!(*fitIntsec)->isIncident(vAttrib.mEit2)) {
                return false;
            }
        }
    }

    return true;
}


}// namespace Makena


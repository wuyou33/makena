#ifndef _MAKENA_INTERSECTION_FINDER_HPP_
#define _MAKENA_INTERSECTION_FINDER_HPP_

#include <memory>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <exception>
#include <stdexcept>
#include <cmath>

#include <directed/di_base.hpp>
#include "loggable.hpp"
#include "primitives.hpp"
#include "manifold.hpp"
#include "intersection_convex_polygon_2d.hpp"
#include "convex_hull_2d.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file intersection_finder.hpp
 *
 * @brief Finds the intersection of two convex polytopes in 3D.
 *        The resuls is one of null, a point, an edge (line segment), 
 *        a convex polygon, or a convex polytope.
 */

namespace Makena {

using namespace std;

using IntSec2D = IntersectionFinderConvexPolygon2D;


class IntersectionFinder : public Loggable {

  public:

    class Attributes {
      public:

        Vec3           mP;
        enum predicate mPred;

        VertexIt       mVit1;
        EdgeIt         mEit1;
        FaceIt         mFit1;

        VertexIt       mVit2;
        EdgeIt         mEit2;
        FaceIt         mFit2;
    };

    /** @brief constructor
     *
     *  @param m1    (in): convex polytope 1 in its own LCS
     *
     *  @param Qmat1 (in): rotation matrix to transfer the points in m1
     *                     to GCS
     *  @param CoM1  (in): Center of mass of m1 in GCS.
     *
     *  @param m2    (in): convex polytope 2 in its own LCS
     *
     *  @param Qmat2 (in): rotation matrix to transfer the points in m2
     *                     to GCS
     *  @param CoM2  (in): Center of mass of m2 in GCS.
     *
     *  @param epsilonZero 
     *               (in): margin whithin in which numerical values are
     *                     considered zero.
     *
     *  @param epsilonZeroPCA
     *               (in): margin whithin in which numerical values are
     *                     considered zero for the spread (variance) in
     *                     the eigen value analysis.
     *
     *  @param epsilonAngle
     *               (in): margin whithin in which numerical values are
     *                     considered zero for dot product and squared norm 
     *                     of cross product of two normalized vectors.
     */
    IntersectionFinder( Manifold& m1, const Mat3x3& Qmat1, const Vec3& CoM1,
                        Manifold& m2, const Mat3x3& Qmat2, const Vec3& CoM2,
                        const double& epsilonZero,
                        const double& epsilonZeroPCA,
                        const double& epsilonAngle,
                        std::ostream& logStream = std::cerr
    );

    /** @brief dectructor */
    ~IntersectionFinder();


    /** @brief main function to construct the intersection and the 
     *         attributes
     *
     *  @param sepAxis (in): pseudo-separation axis of two convex polytopes
     *                       along which the two are separated with zero or
     *                       small translation.
     */
    void               find(const Vec3 sepAxis);
    void               find();


    /** @brief returns the dimension of the resultant intersection.
     *
     *  @return  3 : the result is a convex polytope (3D convex hull)
     *
     *  @return  2 : the result is a convex polygon  (2D convex hull)
     *
     *  @return  1 : the result is a line segment (edge)
     *
     *  @return  0 : the result is a point
     *
     *  @return -1 : the result is empty.
     */
    long               dimension() const;


    /** @brief returns the intersection convex polytope in GCS
     *         if the result is 3 dimensional.
     */
    Manifold&          hull3D();


    /** @brief returns the vertex attributes of the specified vertex in
     *         in the resultant convex polytope.
     *
     *  @param  vit  : pointer to a vertex in the resultant convex polytope
     *
     *  @return associated attributes in Attributes
     *
     *   Attributes::mType : Indicates the source of the vertex in the 
     *                       input convex polytope 1 and convex polytope 2
     *                       in this order. One of the following is set.
     *                  IF_VERTEX_VERTEX
     *                  IF_VERTEX_EDGE
     *                  IF_VERTEX_FACE          
     *                  IF_VERTEX_INTERIOR
     *                  IF_EDGE_VERTEX
     *                  IF_EDGE_EDGE
     *                  IF_EDGE_FACE
     *                  IF_EDGE_INTERIOR
     *                  IF_FACE_VERTEX
     *                  IF_FACE_EDGE
     *                  IF_FACE_FACE
     *                  IF_FACE_INTERIOR
     *                  IF_INTERIOR_VERTEX
     *                  IF_INTERIOR_EDGE
     *                  IF_INTERIOR_FACE
     *
     *   Attributes:mVit1 : Pointer to the original polytope 1, if the type is
     *                      IF_VERTEX_*
     *
     *   Attributes:mEit1 : Pointer to the original polytope 1, if the type is
     *                      IF_EDGE_*
     *
     *   Attributes:mFit1 : Pointer to the original polytope 1, if the type is
     *                      IF_FACE_*
     *
     *   Attributes:mVit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_VERTEX
     *
     *   Attributes:mEit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_EDGE
     *
     *   Attributes:mFit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_FACE
     */
    Attributes&        vertexAttributes3D(VertexIt vit);


    /** @brief returns the edge attributes of the specified edge in
     *         in the resultant convex polytope.
     *
     *  @param  eit  : pointer to an edge in the resultant convex polytope
     *
     *  @return associated attributes in Attributes
     *
     *   Attributes::mType : Indicates the source of the edge in the 
     *                       input convex polytope 1 and convex polytope 2
     *                       in this order. One of the following is set.
     *                  IF_EDGE_EDGE
     *                  IF_EDGE_FACE
     *                  IF_EDGE_INTERIOR
     *                  IF_FACE_EDGE
     *                  IF_FACE_FACE
     *                  IF_FACE_INTERIOR
     *                  IF_INTERIOR_EDGE
     *                  IF_INTERIOR_FACE
     *
     *   Attributes:mEit1 : Pointer to the original polytope 1, if the type is
     *                      IF_EDGE_*
     *
     *   Attributes:mFit1 : Pointer to the original polytope 1, if the type is
     *                      IF_FACE_*
     *
     *   Attributes:mEit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_EDGE
     *
     *   Attributes:mFit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_FACE
     */
    Attributes&        edgeAttributes3D(EdgeIt eit);


    /** @brief returns the face attributes of the specified edge in
     *         in the resultant convex polytope.
     *
     *  @param  fit  : pointer to a face in the resultant convex polytope
     *
     *  @return associated attributes in Attributes
     *
     *   Attributes::mType : Indicates the source of the face in the 
     *                       input convex polytope 1 and convex polytope 2
     *                       in this order. One of the following is set.
     *                  IF_FACE_FACE
     *                  IF_FACE_INTERIOR
     *                  IF_INTERIOR_FACE
     *
     *   Attributes:mFit1 : Pointer to the original polytope 1, if the type is
     *                      IF_FACE_*
     *
     *   Attributes:mFit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_FACE
     */
    Attributes&        faceAttributes3D(FaceIt fit);

    /** @brief returns the size (number of vertices)  of the resultant 
     *                 convex polygon.
     */
    long               size2D() const;


    /** @brief returns the vertex attributes of the specified vertex in
     *         in the resultant convex polygon
     *
     *  @param  index : index into the vertices of the convex polygon, which
     *                  is numbered from 0 to size2D()-1.
     *                  The orientation of the numbering (ccw or cw around
     *                  certain axis) is unspecified.
     *
     *  @return associated attributes in Attributes
     *
     *   Attributes::mType : Indicates the source of the vertex in the 
     *                       input convex polytope 1 and convex polytope 2
     *                       in this order. One of the following is set.
     *                  IF_VERTEX_VERTEX
     *                  IF_VERTEX_EDGE
     *                  IF_VERTEX_FACE          
     *                  IF_VERTEX_INTERIOR
     *                  IF_EDGE_VERTEX
     *                  IF_EDGE_EDGE
     *                  IF_EDGE_FACE
     *                  IF_EDGE_INTERIOR
     *                  IF_FACE_VERTEX
     *                  IF_FACE_EDGE
     *                  IF_FACE_FACE
     *                  IF_FACE_INTERIOR
     *                  IF_INTERIOR_VERTEX
     *                  IF_INTERIOR_EDGE
     *                  IF_INTERIOR_FACE
     *
     *   Attributes:mVit1 : Pointer to the original polytope 1, if the type is
     *                      IF_VERTEX_*
     *
     *   Attributes:mEit1 : Pointer to the original polytope 1, if the type is
     *                      IF_EDGE_*
     *
     *   Attributes:mFit1 : Pointer to the original polytope 1, if the type is
     *                      IF_FACE_*
     *
     *   Attributes:mVit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_VERTEX
     *
     *   Attributes:mEit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_EDGE
     *
     *   Attributes:mFit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_FACE
     */
    Attributes&        vertexAttributes2D(long index);


    /** @brief returns the edge attributes of the specified edge in
     *         in the resultant convex polygon.
     *
     *  @param  index : index into the vertices of the convex polygon, which
     *                  is numbered from 0 to size2D()-1.
     *                  The edge between index and (index + 1)%size2D()
     *                  is the target edge for this operation.
     *
     *  @return associated attributes in Attributes
     *
     *   Attributes::mType : Indicates the source of the edge in the 
     *                       input convex polytope 1 and convex polytope 2
     *                       in this order. One of the following is set.
     *                  IF_EDGE_EDGE
     *                  IF_EDGE_FACE
     *                  IF_EDGE_INTERIOR
     *                  IF_FACE_EDGE
     *                  IF_FACE_FACE
     *                  IF_FACE_INTERIOR
     *                  IF_INTERIOR_EDGE
     *                  IF_INTERIOR_FACE
     *
     *   Attributes:mEit1 : Pointer to the original polytope 1, if the type is
     *                      IF_EDGE_*
     *
     *   Attributes:mFit1 : Pointer to the original polytope 1, if the type is
     *                      IF_FACE_*
     *
     *   Attributes:mEit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_EDGE
     *
     *   Attributes:mFit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_FACE
     */
    Attributes&        edgeAttributes2D(long index);


    /** @brief returns the face attributes in the resultant convex polygon
     *
     *  @return associated attributes in Attributes
     *
     *   Attributes::mType : Indicates the source of the face in the 
     *                       input convex polytope 1 and convex polytope 2
     *                       in this order. One of the following is set.
     *                  IF_FACE_FACE
     *                  IF_FACE_INTERIOR
     *                  IF_INTERIOR_FACE
     *
     *   Attributes:mFit1 : Pointer to the original polytope 1, if the type is
     *                      IF_FACE_*
     *
     *   Attributes:mFit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_FACE
     */
    Attributes&        faceAttributes2D();


    /** @brief returns the vertex attributes of the specified vertex in
     *         in the resultant line segment
     *
     *  @param  index : vertex 0 or 1.
     *
     *  @return associated attributes in Attributes
     *
     *   Attributes::mType : Indicates the source of the vertex in the 
     *                       input convex polytope 1 and convex polytope 2
     *                       in this order. One of the following is set.
     *                  IF_VERTEX_VERTEX
     *                  IF_VERTEX_EDGE
     *                  IF_VERTEX_FACE          
     *                  IF_VERTEX_INTERIOR
     *                  IF_EDGE_VERTEX
     *                  IF_EDGE_EDGE
     *                  IF_EDGE_FACE
     *                  IF_EDGE_INTERIOR
     *                  IF_FACE_VERTEX
     *                  IF_FACE_EDGE
     *                  IF_FACE_FACE
     *                  IF_FACE_INTERIOR
     *                  IF_INTERIOR_VERTEX
     *                  IF_INTERIOR_EDGE
     *                  IF_INTERIOR_FACE
     *
     *   Attributes:mVit1 : Pointer to the original polytope 1, if the type is
     *                      IF_VERTEX_*
     *
     *   Attributes:mEit1 : Pointer to the original polytope 1, if the type is
     *                      IF_EDGE_*
     *
     *   Attributes:mFit1 : Pointer to the original polytope 1, if the type is
     *                      IF_FACE_*
     *
     *   Attributes:mVit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_VERTEX
     *
     *   Attributes:mEit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_EDGE
     *
     *   Attributes:mFit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_FACE
     */
    Attributes&        vertexAttributes1D(long index);


    /** @brief returns the edge attributes of the resultant line segment
     *
     *  @return associated attributes in Attributes
     *
     *   Attributes::mType : Indicates the source of the edge in the 
     *                       input convex polytope 1 and convex polytope 2
     *                       in this order. One of the following is set.
     *                  IF_EDGE_EDGE
     *                  IF_EDGE_FACE
     *                  IF_EDGE_INTERIOR
     *                  IF_FACE_EDGE
     *                  IF_FACE_FACE
     *                  IF_FACE_INTERIOR
     *                  IF_INTERIOR_EDGE
     *                  IF_INTERIOR_FACE
     *
     *   Attributes:mEit1 : Pointer to the original polytope 1, if the type is
     *                      IF_EDGE_*
     *
     *   Attributes:mFit1 : Pointer to the original polytope 1, if the type is
     *                      IF_FACE_*
     *
     *   Attributes:mEit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_EDGE
     *
     *   Attributes:mFit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_FACE
     */
    Attributes&        edgeAttributes1D();


    /** @brief returns the vertex attribute of the resultant point
     *
     *  @return associated attributes in Attributes
     *
     *   Attributes::mType : Indicates the source of the vertex in the 
     *                       input convex polytope 1 and convex polytope 2
     *                       in this order. One of the following is set.
     *                  IF_VERTEX_VERTEX
     *                  IF_VERTEX_EDGE
     *                  IF_VERTEX_FACE          
     *                  IF_VERTEX_INTERIOR
     *                  IF_EDGE_VERTEX
     *                  IF_EDGE_EDGE
     *                  IF_EDGE_FACE
     *                  IF_EDGE_INTERIOR
     *                  IF_FACE_VERTEX
     *                  IF_FACE_EDGE
     *                  IF_FACE_FACE
     *                  IF_FACE_INTERIOR
     *                  IF_INTERIOR_VERTEX
     *                  IF_INTERIOR_EDGE
     *                  IF_INTERIOR_FACE
     *
     *   Attributes:mVit1 : Pointer to the original polytope 1, if the type is
     *                      IF_VERTEX_*
     *
     *   Attributes:mEit1 : Pointer to the original polytope 1, if the type is
     *                      IF_EDGE_*
     *
     *   Attributes:mFit1 : Pointer to the original polytope 1, if the type is
     *                      IF_FACE_*
     *
     *   Attributes:mVit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_VERTEX
     *
     *   Attributes:mEit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_EDGE
     *
     *   Attributes:mFit2 : Pointer to the original polytope 2, if the type is
     *                      IF_*_FACE
     */
    Attributes&        vertexAttributes0D();

#ifdef UNIT_TESTS
  public:
#else
  private:  
#endif

    void prepareVertices(const Vec3 sepAxis);


    void prepareVertices();


    void prepareEdgesAndFaces();


    Mat3x3 rotMatAlignZDirToZAxis(const Vec3& zDir);


    void   rotateWorld(
        const Vec3& n,
        Mat3x3&     newQmat1,
        Vec3&       newCoM1,
        Mat3x3&     newQmat2,
        Vec3&       newCoM2,
        Mat3x3&     matAlignInv
    );


    void findEdgePointsToBeTested(
        const Vec3& pBase,
        const Vec3& pSrc,
        const Vec3& pDst,
        bool&       processSrc,
        bool&       processDst,
        bool&       processMid,
        Vec3&       pMid
    );


    void construct2DInputElementsForFace(
        FaceIt                       fit,
        const Mat3x3&                Qmat,
        const Vec3&                  CoM,
        vector<IntSec2D::InputElem>& elems,
        vector<VertexIt>&            vertices,
        vector<EdgeIt>&              edges,
        Vec3&                        pBase
    );


    void findIntersectionsEdgesOn1FacesOn2();
    void findIntersectionsFacesOn1EdgesOn2();


    void processIntersection2DMid(
        EdgeIt                       eit1,
        vector<IntSec2D::InputElem>& elems1,
        const Vec3&                  pGCS,
        FaceIt                       fit2,
        vector<IntSec2D::InputElem>& elems2,
        vector<VertexIt>&            vertices2,
        vector<EdgeIt>&              edges2
    );


    void processIntersection2DMid(
        FaceIt                       fit1,
        vector<IntSec2D::InputElem>& elems1,
        vector<VertexIt>&            vertices1,
        vector<EdgeIt>&              edges1,
        EdgeIt                       eit2,
        vector<IntSec2D::InputElem>& elems2,
        const Vec3&                  pGCS
    );


    void processIntersection2DEnd(
        FaceIt                       fit1,
        vector<IntSec2D::InputElem>& elems1,
        vector<VertexIt>&            vertices1,
        vector<EdgeIt>&              edges1,
        FaceIt                       fit2,
        vector<IntSec2D::InputElem>& elems2,
        vector<VertexIt>&            vertices2,
        vector<EdgeIt>&              edges2,
        const Vec3&                  pGCS
    );

    void processIntersection2DCoplanar(
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
    );

    Vec3 interpolatePoint(
        IntSec2D::InputElem&  ie1,
        IntSec2D::InputElem&  ie2,
        IntSec2D::OutputElem& oe,
        const Vec3&           p1_3D,
        const Vec3&           p2_3D
    );

    long edgeIndex(long i1, long i2, long numElems);


    void findRemainingInteriorVertices();


    enum predicate classifyVertexAgainstPlane(
        const Vec3& pTest,
        const Vec3& pPlane,
        const Vec3& nPlane
    );

    void constructVertexAttributes();


    void findDimension();


    long findDimensionByPCA();

    void processZeroDimensionalIntersection();

    void processOneDimensionalIntersection();

    void constructEdgeAttributes1D2D();

    void processEdgePredicatesSide1(
        Attributes&     srcAttr,
        Attributes&     dstAttr,
        enum predicate  predIn,
        EdgeIt&         eitOut,
        FaceIt&         fitOut,
        enum predicate& predOut
    );

    void processEdgePredicatesSide2(
        Attributes&     srcAttr,
        Attributes&     dstAttr,
        enum predicate  predIn,
        EdgeIt&         eitOut,
        FaceIt&         fitOut,
        enum predicate& predOut
    );

    void processTwoDimensionalIntersection();

    void constructVertexAttributes1D2D(vector<long>& indices);

    void constructFaceAttributes2D();

    bool checkFaceIncidenceSide12D(FaceIt fit);

    bool checkFaceIncidenceSide22D(FaceIt fit);

    bool findUnieuqFaceSide12D(FaceIt& fit);

    bool findUnieuqFaceSide22D(FaceIt& fit);

    void processThreeDimensionalIntersection();

    void constructVertexAttributes3D();

    enum predicate extractEdgePredicatesSide1(
        Attributes& srcAttr,
        Attributes& dstAttr
    );

    enum predicate extractEdgePredicatesSide2(
        Attributes& srcAttr,
        Attributes& dstAttr
    );

    void constructEdgeAttributes3D();

    enum predicate typeSide1(enum predicate p);
    enum predicate typeSide2(enum predicate p);

    void constructFaceAttributes3D();

    bool findUnieuqFaceSide13D(FaceIt fit, FaceIt& fitOut);

    bool findUnieuqFaceSide23D(FaceIt fit, FaceIt& fitOut);

    bool checkFaceIncidenceSide13D(FaceIt fitPoly1, FaceIt fitIntsec);

    bool checkFaceIncidenceSide23D(FaceIt fitPoly2, FaceIt fitIntsec);   

    Manifold&        mM1;
    Manifold&        mM2;
    const Mat3x3     mQmat1;
    const Mat3x3     mQmat2;
    const Vec3       mCoM1;
    const Vec3       mCoM2;
    const double     mEpsilonZero;
    const double     mEpsilonZeroPCA;
    const double     mEpsilonAngle;

    vector<VertexIt> mActiveVertices1;
    vector<VertexIt> mActiveVertices2;
    vector<EdgeIt>   mActiveEdges1;
    vector<EdgeIt>   mActiveEdges2;
    vector<FaceIt>   mActiveFaces1;
    vector<FaceIt>   mActiveFaces2;

    /** @brief following maps are used to keep track of the intersections */
    map< pair<long,             long             >, Vec3 > mMapVV;
    map< pair<long,             pair<long, long> >, Vec3 > mMapVE;
    map< pair<long,             long             >, Vec3 > mMapVF;
    map< pair<pair<long, long>, long             >, Vec3 > mMapEV;
    map< pair<pair<long, long>, pair<long, long> >, Vec3 > mMapEE;
    map< pair<pair<long, long>, long             >, Vec3 > mMapEF;
    map< pair<long,             long             >, Vec3 > mMapFV;
    map< pair<long,             pair<long, long> >, Vec3 > mMapFE;
    map< long,                                      Vec3 > mMapVI;
    map< long,                                      Vec3 > mMapIV;

    vector<Attributes>        mVertexAttributes;
    vector<Vec3>              mPoints;

    Mat3x3                    mEigenMatrix;
    Vec3                      mSpread;
    Vec3                      mMean;
    long                      mDimension;

    Manifold                  mResultMani;

    map<long,             Attributes> mVertexAttributes3D;
    map<pair<long, long>, Attributes> mEdgeAttributes3D;
    map<long,             Attributes> mFaceAttributes3D;

    vector<Attributes>        mVertexAttributes0D1D2D;
    vector<Attributes>        mEdgeAttributes1D2D;
    Attributes                mFaceAttributes2D;

};


}// namespace Makena

#endif/*_MAKENA_INTERSECTION_FINDER_HPP_*/

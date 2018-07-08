#ifndef _MAKENA_INTERSECTION_DECOMPOSER_HPP_
#define _MAKENA_INTERSECTION_DECOMPOSER_HPP_

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <exception>

#include "primitives.hpp"
#include "manifold.hpp"
#include "convex_rigid_body.hpp"
#include "contact_pair_info.hpp"
#include "broad_phase_aabb_collision_detector.hpp"
#include "obb_obb_test.hpp"
#include "gjk_origin_finder.hpp"
#include "bd_boundary_simplex_finder.hpp"
#include "contact_points_and_normal_generator.hpp"
#include "contact_updater.hpp"
#include "jacobian_constraint.hpp"
#include "constraint_manager.hpp"
#include "intersection_finder.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file intersection_decomposer.hpp
 *
 * @brief Decomposes the surface of the intersection of two convex polytopes
 *        into three categories: boundary, polytope 1, and polytope 2
 */
namespace Makena {


class IntersectionDecomposer :public Loggable {

  public:


    /** @brief constructor
     *
     *  @param logStream(in): Log output
     */
    IntersectionDecomposer(
        const double      epsilonZero,
        const double      epsilonAngle,
        std::ostream&     logStream
    );

    /** @brief descructor. */
    ~IntersectionDecomposer();



#ifdef UNIT_TESTS
  public:
#else
  private:
#endif


    class ConnectedComponent {
      public:
        long             mId;
        enum predicate   mPred;
        vector<VertexIt> mVertices;
        vector<FaceIt>   mFaces;
        vector<EdgeIt>   mEdges;
        vector<EdgeIt>   mDegenerateEdges;
        vector<VertexIt> mDegenerateVertices;
    };

    void constructComponents(IntersectionFinder& finder);

    void resetAndClassifyIndividualFeatures(
        IntersectionFinder& finder
    );

    void resetAndClassifyIndividualFace(
        IntersectionFinder& finder,
        FaceIt              fit
    );

    void resetAndClassifyIndividualEdge(
        IntersectionFinder& finder,
        EdgeIt              eit
    );

    void resetIndividualVertex(
        IntersectionFinder& finder,
        VertexIt            vit
    );

    void findConnectedFeaturesIntoComponents(
        IntersectionFinder&         finder,
        vector<ConnectedComponent>& components
    );

    void markUnvisitedDegenerateEdgesAndVertices(
        vector<ConnectedComponent>& components
    );

    ConnectedComponent findOneFeatureComponentOfPolytope(
        IntersectionFinder& finder,
        FaceIt              fit,
        long                componentId,
        bool                polytope1
    );

    void findOneFeatureComponentOfPolytopeVisitVertex(
        VertexIt             vit,
        const enum predicate THIS_FACE_TYPE,
        const enum predicate THIS_EDGE_TYPE,
        const enum predicate THIS_DEGENERATE_EDGE_TYPE,
        const long           componentId,
        vector<VertexIt>&    verticesToBeReset,
        list<FaceIt>&        BFSStackFaces,
        list<VertexIt>&      BFSStackVertices,
        ConnectedComponent&  comp
    );

    ConnectedComponent findOneFeatureComponentOfPolytope(
        IntersectionFinder& finder,
        EdgeIt              eit,
        long                componentId,
        bool                polytope1
    );

    void findOneFeatureComponentOfPolytopeVisitVertexDegenerate(
        VertexIt             vit,
        const enum predicate EDGE_TYPE,
        const long           componentId,
        vector<VertexIt>&    verticesToBeReset,
        list<VertexIt>&      BFSStack,
        ConnectedComponent&  comp
    );

    ConnectedComponent findOneFeatureComponentOfPolytope(
        IntersectionFinder& finder,
        VertexIt            vit ,
        long                componentId,
        bool                polytope1
    );

    ConnectedComponent findOneFeatureComponentOfBoundary(
        IntersectionFinder& finder,
        EdgeIt              eit,
        long                componentId
    );


    void findOneFeatureComponentOfBoundaryVisitVertex(
        VertexIt             vit,
        const long           componentId,
        list<FaceIt>&        BFSStackFaces,
        list<VertexIt>&      BFSStackVertices,
        ConnectedComponent&  comp
    );

    ConnectedComponent findOneFeatureComponentOfBoundary(
        IntersectionFinder& finder,
        VertexIt            vit,
        long                componentId
    );

    void linkComponents(
        IntersectionFinder&         finder,
        vector<ConnectedComponent>& components
    );

    void tryLinkingComponents(
        ConnectedComponent& comp,
        EdgeIt              eit,
        bool&               found1,
        bool&               found2
    );

    void tryLinkingComponentsByDegenerateEdge(
        ConnectedComponent& comp,
        EdgeIt              eit,
        bool&               found1,
        bool&               found2
    );

    void tryLinkingComponentsByDegenerateVertex(
        ConnectedComponent& comp,
        VertexIt            vit,
        bool&               found1,
        bool&               found2
    );

    bool isUnvisited         (FaceIt fit);
    void markVisited         (FaceIt fit);
    void markUnvisited       (FaceIt fit);
    void setComponentID      (FaceIt fit, long id);
    long componentID         (FaceIt fit);
    void setType             (FaceIt fit, enum predicate p);
    bool isType              (FaceIt fit, enum predicate p);
    enum predicate type      (FaceIt eit);
    bool isPolytopeType      (FaceIt fit);
    bool isPolytopeType1     (FaceIt fit);
    bool isPolytopeType2     (FaceIt fit);

    bool isUnvisited         (EdgeIt eit);
    void markVisited         (EdgeIt eit);
    void markUnvisited       (EdgeIt eit);
    void setComponentID      (EdgeIt eit, long id);
    void setComponentIDAUX   (EdgeIt eit, long id);
    void setType             (EdgeIt eit, enum predicate p);
    bool isType              (EdgeIt eit, enum predicate p);
    void incidentFaces       (EdgeIt eit, FaceIt& fit1, FaceIt& fit2);
    bool isDegenerate        (EdgeIt eit);
    bool isDegenerateEmptyPolytope1
                             (EdgeIt eit);
    void incidentVertices    (EdgeIt eit, VertexIt vit1, VertexIt vit2);
    enum predicate type      (EdgeIt eit);
    bool isBoundary          (EdgeIt eit);
    bool isProperBoundary    (EdgeIt eit);
    bool isDegenerateBoundary(EdgeIt eit);
 
    bool isUnvisited         (VertexIt vit);
    void markVisited         (VertexIt vit);
    void markUnvisited       (VertexIt vit);
    void setComponentID      (VertexIt vit, long id);
    void setComponentIDAUX   (VertexIt vit, long id);
    void setType             (VertexIt vit, enum predicate p);
    bool isType              (VertexIt vit, enum predicate p);
    bool isDegenerate        (VertexIt vit);
    bool isDegenerateEmptyPolytope1
                             (VertexIt vit);
    bool isDegenerateEmptyPolytope2
                             (VertexIt vit);
    FaceIt anIncidentFace    (VertexIt vit);

    const double             mEpsilonZero;
    const double             mEpsilonAngle;
};



}// namespace Makena


#endif /*_MAKENA_INTERSECTION_DECOMPOSER_HPP_*/

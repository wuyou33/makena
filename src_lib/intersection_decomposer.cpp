#include "intersection_decomposer.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file intersection_decomposer.hpp
 *
 * @brief
 */
namespace Makena {

IntersectionDecomposer::IntersectionDecomposer(
    const double      epsilonZero,
    const double      epsilonAngle,
    std::ostream&     logStream
):    
    Loggable(logStream),
    mEpsilonZero(epsilonZero),
    mEpsilonAngle(epsilonAngle)
    {;}


IntersectionDecomposer::~IntersectionDecomposer(){;}


bool IntersectionDecomposer::isUnvisited(FaceIt fit) {
    return (*fit)->mIFcnt==0;
}

void IntersectionDecomposer::markVisited(FaceIt fit) {
    (*fit)->mIFcnt=1;
}

void IntersectionDecomposer::markUnvisited(FaceIt fit) {
    (*fit)->mIFcnt=0;
}


void IntersectionDecomposer::setComponentID(FaceIt fit, long id) {
    (*fit)->mIFcomponentId = id;
}


long IntersectionDecomposer::componentID(FaceIt fit) {
    return (*fit)->mIFcomponentId;
}


void IntersectionDecomposer::setType(FaceIt fit, enum predicate p) {
    (*fit)->mIFflags = p;
}


bool IntersectionDecomposer::isType(FaceIt fit, enum predicate p) {
    return (*fit)->mIFflags == p;
}


enum predicate IntersectionDecomposer::type(FaceIt fit) {
    return (*fit)->mIFflags;
}


bool IntersectionDecomposer::isPolytopeType(FaceIt fit) {
    return (*fit)->mIFflags == IF_FACE_POLYTOPE_1 ||
           (*fit)->mIFflags == IF_FACE_POLYTOPE_2 ;
}


bool IntersectionDecomposer::isPolytopeType1(FaceIt fit) {
    return (*fit)->mIFflags == IF_FACE_POLYTOPE_1;
}


bool IntersectionDecomposer::isPolytopeType2(FaceIt fit) {
    return (*fit)->mIFflags == IF_FACE_POLYTOPE_2;
}


bool IntersectionDecomposer::isUnvisited(EdgeIt eit) {
    return (*eit)->mIFcnt==0;
}

void IntersectionDecomposer::markVisited(EdgeIt eit) {
    (*eit)->mIFcnt=1;
}

void IntersectionDecomposer::markUnvisited(EdgeIt eit) {
    (*eit)->mIFcnt=0;
}


void IntersectionDecomposer::setComponentID(EdgeIt eit, long id) {
    (*eit)->mIFcomponentId=id;
}


void IntersectionDecomposer::setComponentIDAUX(EdgeIt eit, long id) {
    (*eit)->mIFcomponentIdaux=id;
}


void IntersectionDecomposer::setType(EdgeIt eit, enum predicate p) {
    (*eit)->mIFflags = p;
}


bool IntersectionDecomposer::isType(EdgeIt eit, enum predicate p) {
    return (*eit)->mIFflags == p;
}


enum predicate IntersectionDecomposer::type(EdgeIt eit) {
    return (*eit)->mIFflags;
}


void IntersectionDecomposer::incidentFaces(EdgeIt eit, FaceIt& fit1, FaceIt& fit2) {
    auto he1  = (*eit)->he1();
    auto he2  = (*eit)->he2();
    fit1      = (*he1)->face();
    fit2      = (*he2)->face();
}


bool IntersectionDecomposer::isBoundary(EdgeIt eit) {

    switch ( type(eit)  ) {

      case IF_EDGE_BOUNDARY_INTERIOR:
      case IF_EDGE_BOUNDARY_PROPER:
      case IF_EDGE_BOUNDARY_ON_POLYTOPE_1:
      case IF_EDGE_BOUNDARY_ON_POLYTOPE_2:
      case IF_EDGE_DEGENERATE_EMPTY_POLYTOPE_1:
      case IF_EDGE_DEGENERATE_EMPTY_POLYTOPE_2:
        return true;

      default:
        return false;

    }
}


bool IntersectionDecomposer::isProperBoundary(EdgeIt eit) {

    switch ( type(eit)  ) {

      case IF_EDGE_BOUNDARY_PROPER:
      case IF_EDGE_BOUNDARY_ON_POLYTOPE_1:
      case IF_EDGE_BOUNDARY_ON_POLYTOPE_2:
        return true;

      default:
        return false;
    }
}


bool IntersectionDecomposer::isDegenerate(EdgeIt eit) {

    return (*eit)->mIFflags ==IF_EDGE_DEGENERATE_EMPTY_POLYTOPE_1 ||
           (*eit)->mIFflags ==IF_EDGE_DEGENERATE_EMPTY_POLYTOPE_2 ;

}


void IntersectionDecomposer::incidentVertices(
    EdgeIt   eit, 
    VertexIt vit1,
    VertexIt vit2
) {
    auto he = (*eit)->he1();
    vit1    = (*he)->src();
    vit2    = (*he)->dst();
}


bool IntersectionDecomposer::isDegenerateEmptyPolytope1(EdgeIt eit) {
    return (*eit)->mIFflags == IF_EDGE_DEGENERATE_EMPTY_POLYTOPE_1;
}


bool IntersectionDecomposer::isUnvisited(VertexIt vit) {
    return (*vit)->mIFcnt==0;
}

void IntersectionDecomposer::markVisited(VertexIt vit) {
    (*vit)->mIFcnt=1;
}

void IntersectionDecomposer::markUnvisited(VertexIt vit) {
    (*vit)->mIFcnt=0;
}


void IntersectionDecomposer::setComponentID(VertexIt vit, long id) {
    (*vit)->mIFcomponentId=id;
}


void IntersectionDecomposer::setComponentIDAUX(VertexIt vit, long id) {
    (*vit)->mIFcomponentIdaux=id;
}


void IntersectionDecomposer::setType(VertexIt vit, enum predicate p) {
    (*vit)->mIFflags = p;
}


bool IntersectionDecomposer::isType(VertexIt vit, enum predicate p) {
    return (*vit)->mIFflags == p;
}


bool IntersectionDecomposer::isDegenerate(VertexIt vit)
{

   long cnt_IF_FACE_POLYTOPE_1 = 0;
   long cnt_IF_FACE_POLYTOPE_2 = 0;

   long cnt_IF_EDGE_POLYTOPE_1_INTERIOR = 0;
   long cnt_IF_EDGE_POLYTOPE_2_INTERIOR = 0;

    for (auto he : (*vit)->halfEdges()) {

        if ((*he)->src()==vit) {

            auto eit = (*he)->edge();

            switch (type(eit)) {

              case IF_EDGE_POLYTOPE_1_INTERIOR:
                cnt_IF_EDGE_POLYTOPE_1_INTERIOR++;
                break;

              case IF_EDGE_POLYTOPE_2_INTERIOR:
                cnt_IF_EDGE_POLYTOPE_1_INTERIOR++;
                break;

              default:
                return false;
            }

            auto fit = (*he)->face();

            switch (type(fit)) {

              case IF_FACE_POLYTOPE_1:
                cnt_IF_FACE_POLYTOPE_1++;
                break;

              case IF_FACE_POLYTOPE_2:
                cnt_IF_FACE_POLYTOPE_2++;
                break;

              default:
                return false;
            }

        }

    }

    if ( ( cnt_IF_FACE_POLYTOPE_1          > 0 && 
           cnt_IF_FACE_POLYTOPE_2          > 0    ) ||
         ( cnt_IF_EDGE_POLYTOPE_1_INTERIOR > 0 && 
           cnt_IF_EDGE_POLYTOPE_2_INTERIOR > 0    ) ||  
         ( cnt_IF_FACE_POLYTOPE_1          > 0 && 
           cnt_IF_EDGE_POLYTOPE_2_INTERIOR > 0    ) ||
         ( cnt_IF_FACE_POLYTOPE_2          > 0 && 
           cnt_IF_EDGE_POLYTOPE_1_INTERIOR > 0    )   ) {

        return false;
    }
    else if (cnt_IF_FACE_POLYTOPE_1 > 0) {

        setType(vit, IF_VERTEX_DEGENERATE_EMPTY_POLYTOPE_2);
        return true;

    }
    else if (cnt_IF_FACE_POLYTOPE_2 > 0) {

        setType(vit, IF_VERTEX_DEGENERATE_EMPTY_POLYTOPE_1);
        return true;

    }
    return false;
}


bool IntersectionDecomposer::isDegenerateEmptyPolytope1(VertexIt vit)
{
    return (*vit)->mIFflags == IF_VERTEX_DEGENERATE_EMPTY_POLYTOPE_1;
}


bool IntersectionDecomposer::isDegenerateEmptyPolytope2(VertexIt vit)
{
    return (*vit)->mIFflags == IF_VERTEX_DEGENERATE_EMPTY_POLYTOPE_2;
}


FaceIt IntersectionDecomposer::anIncidentFace(VertexIt vit)
{
    auto he = *((*vit)->halfEdges()).begin();
    return (*he)->face();
}


void IntersectionDecomposer::constructComponents(IntersectionFinder& finder) {

    resetAndClassifyIndividualFeatures(finder);

    vector<ConnectedComponent> components;

    findConnectedFeaturesIntoComponents(finder, components);

}


void IntersectionDecomposer::resetAndClassifyIndividualFeatures(
    IntersectionFinder& finder
)
{
    auto& ch = finder.hull3D();

    auto fPair = ch.faces();
    for (auto fit = fPair.first; fit != fPair.second; fit++) {

        resetAndClassifyIndividualFace(finder, fit);

    }

    auto ePair = ch.edges();
    for (auto eit = ePair.first; eit != ePair.second; eit++) {

        resetAndClassifyIndividualEdge(finder, eit);

    }

    auto vPair = ch.vertices();
    for (auto vit = vPair.first; vit != vPair.second; vit++) {

        resetIndividualVertex(finder, vit);

    }
}


void IntersectionDecomposer::resetAndClassifyIndividualFace(
    IntersectionFinder& finder,
    FaceIt              fit
) {        
    markUnvisited(fit);
    setComponentID(fit, 0);

    auto& fatt = finder.faceAttributes3D(fit);
    switch (fatt.mPred) {

      case IF_FACE_FACE:

        setType(fit, IF_FACE_BOUNDARY );
        break;

      case IF_FACE_INTERIOR:

        setType(fit, IF_FACE_POLYTOPE_1 );
        break;

      case IF_INTERIOR_FACE:

        setType(fit, IF_FACE_POLYTOPE_2 );
        break;

      default:
        setType(fit, NONE);
        break;
    }
}


void IntersectionDecomposer::resetAndClassifyIndividualEdge(
    IntersectionFinder& finder,
    EdgeIt              eit
) {        
    
    markUnvisited(eit);
    setComponentID(eit, 0);
    setComponentIDAUX(eit, 0);

    auto& eatt = finder.edgeAttributes3D(eit);
    switch (eatt.mPred) {

      case IF_EDGE_EDGE:
      case IF_EDGE_FACE:
      case IF_FACE_EDGE:
      case IF_FACE_FACE:
        {
            FaceIt fit1, fit2;
            incidentFaces(eit, fit1, fit2);

            if ( isType(fit1, IF_FACE_BOUNDARY) ) {

                if ( isType(fit2, IF_FACE_BOUNDARY) ) {

                    setType(eit, IF_EDGE_BOUNDARY_INTERIOR);

                }
                else if ( isType(fit2, IF_FACE_POLYTOPE_1) ) {

                    setType(eit, IF_EDGE_BOUNDARY_ON_POLYTOPE_1);

                }
                else if ( isType(fit2, IF_FACE_POLYTOPE_2) ) {

                    setType(eit, IF_EDGE_BOUNDARY_ON_POLYTOPE_2);  

                }
                else {

                    setType(eit, NONE);

                }
            }
            else if ( isType(fit1, IF_FACE_POLYTOPE_1) ) {

                if ( isType(fit2, IF_FACE_BOUNDARY) ) {

                    setType(eit, IF_EDGE_BOUNDARY_ON_POLYTOPE_1);

                }
                else if ( isType(fit2, IF_FACE_POLYTOPE_1) ) {

                    setType(eit, IF_EDGE_DEGENERATE_EMPTY_POLYTOPE_1);

                }
                else if ( isType(fit2, IF_FACE_POLYTOPE_2) ) {

                    setType(eit, IF_EDGE_BOUNDARY_PROPER);

                }
                else {

                    setType(eit, NONE);

                }
            }
            else if ( isType(fit1, IF_FACE_POLYTOPE_2) ) {

                if ( isType(fit2, IF_FACE_BOUNDARY) ) {

                    setType(eit, IF_EDGE_BOUNDARY_ON_POLYTOPE_2);

                }
                else if ( isType(fit2, IF_FACE_POLYTOPE_1) ) {

                    setType(eit, IF_EDGE_BOUNDARY_PROPER);

                }
                else if ( isType(fit2, IF_FACE_POLYTOPE_2) ) {

                    setType(eit, IF_EDGE_DEGENERATE_EMPTY_POLYTOPE_2);

                }
                else {

                    setType(eit,  NONE);

                }
            }
            else {

                setType(eit,  NONE);

            }
        }
        break;

      case IF_EDGE_INTERIOR:
      case IF_FACE_INTERIOR:

        setType(eit, IF_EDGE_POLYTOPE_1_INTERIOR);
        break;

      case IF_INTERIOR_EDGE:
      case IF_INTERIOR_FACE:

        setType(eit, IF_EDGE_POLYTOPE_2_INTERIOR);
        break;

      default:
        setType(eit,  NONE);
        break;
    }
}


void IntersectionDecomposer::resetIndividualVertex(
    IntersectionFinder& finder,
    VertexIt            vit
) {        
    markUnvisited (vit);
    setComponentID(vit, 0);
    setComponentIDAUX(vit, 0);
    setType       (vit, NONE);
}


void IntersectionDecomposer::findConnectedFeaturesIntoComponents(
    IntersectionFinder&         finder,
    vector<ConnectedComponent>& components
) {
    long componentId = 0;

    auto& ch = finder.hull3D();


    // Starting from an unvisited component face, find all the connected 
    // features on either of the polytopes.

    auto fPair = ch.faces();
    for (auto fit = fPair.first; fit != fPair.second; fit++ ) {

        if ( isUnvisited(fit) && isPolytopeType(fit) ) {

            auto comp = findOneFeatureComponentOfPolytope(
                            finder,
                            fit,
                            componentId,
                            isPolytopeType1(fit)
                        );

            components.push_back(std::move(comp));
            componentId++;
        }        

    }


    // Starting from an unvisited degenerate edge, find all the connected 
    // features (other incident degenerate edges) on either of the polytopes.
    // The resultant component forms a connected tree.

    auto ePair = ch.edges();
    for (auto eit = ePair.first; eit != ePair.second; eit++) {

        if ( isUnvisited(eit) && isDegenerate(eit)) {

            auto comp = findOneFeatureComponentOfPolytope(
                            finder,
                            eit,
                            componentId,
                            isDegenerateEmptyPolytope1(eit)
                        );

            components.push_back(std::move(comp));
            componentId++;
        }        

    }


    // Find all the degenerate vertex, which are all isolated (i.e., surrounded
    // by the faces of either of the polytopes), make a component of polytope
    // type.
    auto vPair = ch.vertices();
    for (auto vit = vPair.first; vit != vPair.second; vit++) {

        if (isDegenerate(vit)) {
            auto comp = findOneFeatureComponentOfPolytope(
                            finder,
                            vit,
                            componentId,
                            isDegenerateEmptyPolytope1(vit)
                        );

            components.push_back(std::move(comp));
            componentId++;
        }        

    }


    // Mark all the degenerate edges and vertices unvisited, as they will have
    //  to be explored to find the boundary compopnents.

    markUnvisitedDegenerateEdgesAndVertices(components);


    // Starting from an unvisited degenerate or a boundary edge, find all the
    // connected features of the boundary type.

    for (auto eit = ePair.first; eit != ePair.second; eit++) {

        if ( isUnvisited(eit) && isBoundary(eit) ) {

            auto comp = findOneFeatureComponentOfBoundary(
                            finder,
                            eit,
                            componentId
                        );

            components.push_back(std::move(comp));
            componentId++;
        }        
    }

    // Find all the degenerate vertex, and make a component of boundary
    // of itself.

    for (auto vit = vPair.first; vit != vPair.second; vit++) {

        if ( isDegenerateEmptyPolytope1(vit) ||
             isDegenerateEmptyPolytope2(vit)    ) {

            auto comp = findOneFeatureComponentOfBoundary(
                            finder,
                            vit,
                            componentId
                        );

            components.push_back(std::move(comp));
            componentId++;
        }        
    }

}


void IntersectionDecomposer::markUnvisitedDegenerateEdgesAndVertices(
    vector<ConnectedComponent>& components
) {
    for (auto& comp : components) {
        for (auto eit : comp.mDegenerateEdges) {
            markUnvisited(eit);
        }
        for (auto vit : comp.mDegenerateVertices) {
            markUnvisited(vit);
        }
    }
}


IntersectionDecomposer::ConnectedComponent 
IntersectionDecomposer::findOneFeatureComponentOfPolytope(
    IntersectionFinder& finder,
    FaceIt              fit ,
    long                componentId,
    bool                polytope1
) {
    const enum predicate THIS_FACE_TYPE = polytope1 ?
                              IF_FACE_POLYTOPE_1 :
                              IF_FACE_POLYTOPE_2 ;

    const enum predicate THIS_EDGE_TYPE = polytope1 ?
                              IF_EDGE_POLYTOPE_1_INTERIOR :
                              IF_EDGE_POLYTOPE_2_INTERIOR ;  

    const enum predicate THIS_DEGENERATE_EDGE_TYPE = polytope1 ?
                              IF_EDGE_DEGENERATE_EMPTY_POLYTOPE_1 :
                              IF_EDGE_DEGENERATE_EMPTY_POLYTOPE_2 ;

    ConnectedComponent comp;
    comp.mId   = componentId;
    comp.mPred = polytope1 ? IF_POLYTOPE_1 : IF_POLYTOPE_2;


    vector<VertexIt> verticesToBeReset;
    list<FaceIt>     BFS_StackFaces;
    list<VertexIt>   BFS_StackVertices;

    BFS_StackFaces.push_back(fit);
    markVisited(fit);
    setComponentID(fit, componentId);
    comp.mFaces.push_back(fit);

    while ( !BFS_StackFaces.empty() || !BFS_StackVertices.empty() ) {

        if (!BFS_StackFaces.empty() ) {

            auto fit = (*BFS_StackFaces.begin());
            BFS_StackFaces.pop_front();

            for (auto he : (*fit)->halfEdges()) {

                findOneFeatureComponentOfPolytopeVisitVertex(
                    (*he)->src(),
                    THIS_FACE_TYPE,
                    THIS_EDGE_TYPE,
                    THIS_DEGENERATE_EDGE_TYPE,
                    componentId,
                    verticesToBeReset,
                    BFS_StackFaces,
                    BFS_StackVertices,
                    comp
                );
            }
        }
        else if (!BFS_StackVertices.empty() ) {

            auto vit = (*BFS_StackVertices.begin());
            BFS_StackVertices.pop_front();

            findOneFeatureComponentOfPolytopeVisitVertex(
                vit,
                THIS_FACE_TYPE,
                THIS_EDGE_TYPE,
                THIS_DEGENERATE_EDGE_TYPE,
                componentId,
                verticesToBeReset,
                BFS_StackFaces,
                BFS_StackVertices,
                comp
            );
        }
    }

    for (auto vit : verticesToBeReset) {
        markUnvisited(vit);
    }

    return comp;
}


void IntersectionDecomposer::findOneFeatureComponentOfPolytopeVisitVertex(
    VertexIt             vit,
    const enum predicate FACE_TYPE,
    const enum predicate EDGE_TYPE,
    const enum predicate DEGENERATE_EDGE_TYPE,
    const long           componentId,
    vector<VertexIt>&    verticesToBeReset,
    list<FaceIt>&        BFS_StackFaces,
    list<VertexIt>&      BFS_StackVertices,
    ConnectedComponent&  comp
) {

    for ( auto he : (*vit)->halfEdges() ) {

        if ( (*he)->src() == vit ) {

            auto fit = (*he)->face();

            if ( isUnvisited(fit) && isType(fit, FACE_TYPE) ) {

                 markVisited(fit);
                 setComponentID(fit, componentId);
                 comp.mFaces.push_back(fit);
                 BFS_StackFaces.push_back(fit);
             }

             auto eit = (*he)->edge();

             if ( isUnvisited(eit) && isType(eit, EDGE_TYPE) ) {

                 markVisited(eit);
                 setComponentID(eit, componentId);
                 comp.mEdges.push_back(eit);

             }
             else if ( isUnvisited(eit) && isType(eit, DEGENERATE_EDGE_TYPE) ) {

                 markVisited(eit);
                 setComponentIDAUX(eit, componentId);
                 comp.mDegenerateEdges.push_back(eit);

                 auto vit = (*he)->dst();
                 if ( isUnvisited(vit) ) {

                     markVisited(vit);
                     BFS_StackVertices.push_back(vit);
                     verticesToBeReset.push_back(vit);

                 }
             }
         }
    }
}



// Rarely occur. This is to detect a connected tree of edges
IntersectionDecomposer::ConnectedComponent 
IntersectionDecomposer::findOneFeatureComponentOfPolytope(
    IntersectionFinder& finder,
    EdgeIt              eit ,
    long                componentId,
    bool                polytope1
) {
    const enum predicate EDGE_TYPE = polytope1 ?
                                     IF_EDGE_DEGENERATE_EMPTY_POLYTOPE_1 :
                                     IF_EDGE_DEGENERATE_EMPTY_POLYTOPE_2 ;

    ConnectedComponent comp;
    comp.mId   = componentId;
    comp.mPred = polytope1 ? IF_POLYTOPE_1 : IF_POLYTOPE_2 ;

    vector<VertexIt> verticesToBeReset;
    list<VertexIt>   BFS_Stack;

    markVisited(eit);
    setComponentID(eit, componentId);
    comp.mDegenerateEdges.push_back(eit);

    VertexIt vit1, vit2;
    incidentVertices(eit, vit1, vit2);

    markVisited(vit1);
    markVisited(vit2);

    BFS_Stack.push_back(vit1);
    BFS_Stack.push_back(vit2);

    while ( !BFS_Stack.empty() ) {

        if (!BFS_Stack.empty() ) {

            auto vit = (*BFS_Stack.begin());
            BFS_Stack.pop_front();

            findOneFeatureComponentOfPolytopeVisitVertexDegenerate(
                vit,
                EDGE_TYPE,
                componentId,
                verticesToBeReset,
                BFS_Stack,
                comp
            );
        }
    }

    for (auto vit : verticesToBeReset) {
        markUnvisited(vit);
    }

    return comp;
}


void IntersectionDecomposer::findOneFeatureComponentOfPolytopeVisitVertexDegenerate(
    VertexIt             vit,
    const enum predicate EDGE_TYPE,
    const long           componentId,
    vector<VertexIt>&    verticesToBeReset,
    list<VertexIt>&      BFS_Stack,
    ConnectedComponent&  comp
) {

    for (auto he : (*vit)->halfEdges()) {

        if ((*he)->src()==vit) {

            auto eit = (*he)->edge();

            if ( isUnvisited(eit) && isType(eit, EDGE_TYPE) ) {

                markVisited(eit);
                setComponentIDAUX(eit, componentId);
                comp.mDegenerateEdges.push_back(eit);

                auto ait = (*he)->dst();
                if ( isUnvisited(ait) ) {

                    markVisited(ait);
                    BFS_Stack.push_back(ait);
                    verticesToBeReset.push_back(ait);

                }
            }
        }
    }
}


// Rarely occur. This is to handle an isolated vertex
IntersectionDecomposer::ConnectedComponent 
IntersectionDecomposer::findOneFeatureComponentOfPolytope(
    IntersectionFinder& finder,
    VertexIt            vit ,
    long                componentId,
    bool                polytope1
) {
    ConnectedComponent comp;
    comp.mId   = componentId;
    comp.mPred = polytope1 ? IF_POLYTOPE_1 : IF_POLYTOPE_2 ;

    setComponentIDAUX(vit, componentId);
    comp.mDegenerateVertices.push_back(vit);

    return comp;
}


IntersectionDecomposer::ConnectedComponent 
IntersectionDecomposer::findOneFeatureComponentOfBoundary(
    IntersectionFinder& finder,
    EdgeIt              eit,
    long                componentId
) {

    ConnectedComponent comp;
    comp.mId   = componentId;
    comp.mPred = IF_BOUNDARY;

    list<FaceIt>     BFS_StackFaces;
    list<VertexIt>   BFS_StackVertices;

    markVisited(eit);
    setComponentID(eit, componentId);

    if (isDegenerate(eit)) {
        comp.mDegenerateEdges.push_back(eit);
    }
    else {
        comp.mEdges.push_back(eit);
    }

    VertexIt vit1, vit2;
    incidentVertices(eit, vit1, vit2);

    markVisited(vit1);
    markVisited(vit2);

    BFS_StackVertices.push_back(vit1);
    BFS_StackVertices.push_back(vit2);

    while ( !BFS_StackFaces.empty() || !BFS_StackVertices.empty() ) {

        if (!BFS_StackFaces.empty() ) {

            auto fit = (*BFS_StackFaces.begin());
            BFS_StackFaces.pop_front();

            for (auto he : (*fit)->halfEdges()) {

                findOneFeatureComponentOfBoundaryVisitVertex(
                    (*he)->src(),
                    componentId,
                    BFS_StackFaces,
                    BFS_StackVertices,
                    comp
                );
            }
        }
        else if (!BFS_StackVertices.empty() ) {

            auto vit = (*BFS_StackVertices.begin());
            BFS_StackVertices.pop_front();

            findOneFeatureComponentOfBoundaryVisitVertex(
                vit,
                componentId,
                BFS_StackFaces,
                BFS_StackVertices,
                comp
            );
        }
    }

    return comp;
}


void IntersectionDecomposer::findOneFeatureComponentOfBoundaryVisitVertex(
    VertexIt             vit,
    const long           componentId,
    list<FaceIt>&        BFS_StackFaces,
    list<VertexIt>&      BFS_StackVertices,
    ConnectedComponent&  comp
) {

    for ( auto he : (*vit)->halfEdges() ) {

        if ( (*he)->src() == vit ) {

            auto fit = (*he)->face();

            if ( isUnvisited(fit) && isType(fit, IF_FACE_BOUNDARY) ) {

                markVisited(fit);
                setComponentID(fit, componentId);
                comp.mFaces.push_back(fit);
                BFS_StackFaces.push_back(fit);

            }

            auto eit = (*he)->edge();

            if ( isUnvisited(eit) && isBoundary(eit) ) {

                markVisited(eit);
                setComponentID(eit, componentId);

                if ( isDegenerate(eit) ) {

                    comp.mDegenerateEdges.push_back(eit);

                    auto ait = (*he)->dst();
                    if ( isUnvisited(ait) ){
                        markVisited(ait);
                        BFS_StackVertices.push_back(ait);
                    }

                }
                else {

                    comp.mEdges.push_back(eit);

                }
            }
        }
    }
}


// Rarely occur. This is to handle an isolated vertex
IntersectionDecomposer::ConnectedComponent 
IntersectionDecomposer::findOneFeatureComponentOfBoundary(
    IntersectionFinder& finder,
    VertexIt            vit ,
    long                componentId
) {

    ConnectedComponent comp;
    comp.mId   = componentId;
    comp.mPred = IF_BOUNDARY;

    setComponentID(vit, componentId);
    comp.mDegenerateVertices.push_back(vit);

    return comp;

}



}// namespace Makena

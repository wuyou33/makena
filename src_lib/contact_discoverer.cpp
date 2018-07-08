#include "contact_discoverer.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif

/**
 * @file contact_discoverer.cpp
 *
 * @brief
 */
namespace Makena {


ContactDiscoverer::ContactDiscoverer(
    ConvexRigidBody& body1,
    ConvexRigidBody& body2,
    const double     epsilonZero,
    const double     epsilonZeroPCA,
    const double     epsilonAngle,
    const double     obbScalingStep,
    const long       obbScalingNumExtra,
    const bool       useTempGeomConfig,
    std::ostream&    logStream
):
    Loggable           (logStream),
    mBody1             (body1),
    mBody2             (body2),
    mUseTempGeomConfig (useTempGeomConfig),
    mQmat1             (useTempGeomConfig ? body1.QmatTemp() : body1.Qmat()),
    mCoM1              (useTempGeomConfig ? body1.CoMTemp()  : body1.CoM() ),
    mQmat2             (useTempGeomConfig ? body2.QmatTemp() : body2.Qmat()),
    mCoM2              (useTempGeomConfig ? body2.CoMTemp()  : body2.CoM() ),
    mVlin1             (useTempGeomConfig ? body1.VlinTemp() : body1.Vlin()),
    mVang1             (useTempGeomConfig ? body1.VangTemp() : body1.Vang()),
    mVlin2             (useTempGeomConfig ? body2.VlinTemp() : body2.Vlin()),
    mVang2             (useTempGeomConfig ? body2.VangTemp() : body2.Vang()),
    mEpsilonZero       (epsilonZero),
    mEpsilonZeroPCA    (epsilonZeroPCA),
    mEpsilonAngle      (epsilonAngle),
    mObbScalingStep    ( std::min( std::max( obbScalingStep, 0.01 ), 0.99 ) ),
    mObbScalingNumExtra(obbScalingNumExtra),
    mISFinder          (

        body1.ConvexHull(), mQmat1, mCoM1,
        body2.ConvexHull(), mQmat2, mCoM2,
        epsilonZero,
        epsilonZero*100.0,
        epsilonZero,
        logStream

    )
    {;}


ContactDiscoverer::~ContactDiscoverer(){;}


ContactPairInfo ContactDiscoverer::discover()
{

    findSeparationAxesFromOBBscaling();

    if ( mSeparationAxes1To2.size() > 0 ) {

        mISFinder.find(mSeparationAxes1To2[0]);
    }
    else {

        mISFinder.find();
    }

    auto dimPM = mISFinder.dimension();

    switch (dimPM) {
      case -1:
        mInfo.mType1 = ContactPairInfo::FT_NONE;
        mInfo.mType2 = ContactPairInfo::FT_NONE;
        break;

      case 0:
        processPM_0D();
        break;

      case 1:
        processPM_1D();
        break;

      case 2:
        processPM_2D();
        break;

      case 3:
        {
            bool poly1inside, poly2inside;

            isOnePolytopeInteriorOfTheOther( poly1inside, poly2inside );

            if (poly1inside) {
                processPM_Poly1Inside();

            }
            else if (poly2inside) {
                processPM_Poly2Inside();

            }
            else {
                processPM_3D();

            }
        }
        break;

      default:
        break;
    }

    if (mInfo.mType1!=ContactPairInfo::FT_NONE) {
        updateContactNormal();
    }
    return mInfo;

}


void ContactDiscoverer::findSeparationAxesFromOBBscaling()
{

    double scaling = 1.0; 

    for (  ;
           scaling > mEpsilonZero; 
           scaling -= mObbScalingStep  ) {

        Vec3 sepAxis;
        auto doesIntersect = Makena::doIntersect(

            mBody1.OBBcenter(), mBody1.OBBaxes(), mBody1.OBBextents(),
            mBody2.OBBcenter(), mBody2.OBBaxes(), mBody2.OBBextents(),
            scaling, sepAxis                                          
        );

        if ( !doesIntersect ) {

            mSeparationAxes1To2.push_back(sepAxis);
            break;

        }
    }

    // Run 3 more times to find further separation axes until we find 3.
    for ( long i = 0 ; 

          i < mObbScalingNumExtra        &&
          mSeparationAxes1To2.size() < 3 && 
          scaling > mEpsilonZero         ;

          scaling -= mObbScalingStep, i++    ) {

        Makena::findSeparationAxes(

            mBody1.OBBcenter(), mBody1.OBBaxes(), mBody1.OBBextents(),
            mBody2.OBBcenter(), mBody2.OBBaxes(), mBody2.OBBextents(),
            scaling, mSeparationAxes1To2, mEpsilonZero
        );
    }
}


void ContactDiscoverer::processPM_0D()
{
    auto  att = mISFinder.vertexAttributes0D();

    switch (att.mPred) {
      case IF_VERTEX_VERTEX:

        mInfo.mType1 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit1  = att.mVit1;
        mInfo.mType2 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit2  = att.mVit2;
        break;

      case IF_VERTEX_EDGE:
        mInfo.mType1 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit1 = att.mVit1;
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2 = att.mEit2;
        break;

      case IF_VERTEX_FACE:
        mInfo.mType1 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit1 = att.mVit1;
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2 = att.mFit2;
        break;

      case IF_EDGE_VERTEX:
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1 = att.mEit1;
        mInfo.mType2 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit2 = att.mVit2;
        break;

      case IF_FACE_VERTEX:
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1 = att.mFit1;
        mInfo.mType2 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit2 = att.mVit2;
        break;

      case IF_EDGE_EDGE:
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1 = att.mEit1;
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2 = att.mEit2;
        break;

      case IF_EDGE_FACE:
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1 = att.mEit1;
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2 = att.mFit2;
        break;

      case IF_FACE_EDGE:
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1 = att.mFit1;
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2 = att.mEit2;
        break;

      case IF_FACE_FACE:
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1 = att.mFit1;
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2 = att.mFit2;
        break;

      default:
        log( Loggable::ERROR, __FILE__, __LINE__, 
             "Invalid 0-dim penetration manifold"  );
        break;
    }
}


void ContactDiscoverer::processPM_1D()
{
    auto  att = mISFinder.edgeAttributes1D();

    switch (att.mPred) {

      case IF_EDGE_EDGE:
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1 = att.mEit1;
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2 = att.mEit2;
        break;

      case IF_EDGE_FACE:
        mInfo.mType1 = ContactPairInfo::FT_EDGE;
        mInfo.mEit1 = att.mEit1;
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2 = att.mFit2;
        break;

      case IF_FACE_EDGE:
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1 = att.mFit1;
        mInfo.mType2 = ContactPairInfo::FT_EDGE;
        mInfo.mEit2 = att.mEit2;
        break;

      case IF_FACE_FACE:
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1 = att.mFit1;
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2 = att.mFit2;
        break;

      default:
        log( Loggable::ERROR, __FILE__, __LINE__, 
             "Invalid 1-dim penetration manifold"  );
        break;
    }
}


void ContactDiscoverer::processPM_2D()
{
    auto  att = mISFinder.faceAttributes2D();
    switch (att.mPred) {

      case IF_FACE_FACE:
        mInfo.mType1 = ContactPairInfo::FT_FACE;
        mInfo.mFit1 = att.mFit1;
        mInfo.mType2 = ContactPairInfo::FT_FACE;
        mInfo.mFit2 = att.mFit2;
        break;

      default:
        log( Loggable::ERROR, __FILE__, __LINE__, 
             "Invalid 2-dim penetration manifold"  );
        break;
    }
}


void ContactDiscoverer::isOnePolytopeInteriorOfTheOther(
    bool& polytope1Inside,
    bool& polytope2Inside
) {
    polytope1Inside = true;
    polytope2Inside = true;

    auto& ch    = mISFinder.hull3D();
    auto  vPair = ch.vertices();
    for ( auto vit = vPair.first; vit != vPair.second; vit++ ) {

        auto att  = mISFinder.vertexAttributes3D(vit);

        switch (att.mPred) {

          case IF_VERTEX_VERTEX:
          case IF_VERTEX_EDGE:
          case IF_VERTEX_FACE:
          case IF_EDGE_VERTEX:
          case IF_EDGE_EDGE:
          case IF_EDGE_FACE:
          case IF_FACE_VERTEX:
          case IF_FACE_EDGE:
          case IF_FACE_FACE:
            polytope1Inside = false;
            polytope2Inside = false;
            return;            

          case IF_VERTEX_INTERIOR:
          case IF_EDGE_INTERIOR:
          case IF_FACE_INTERIOR:
            polytope1Inside = false;
            if (!polytope2Inside) {
                return;
            }
            break;

          case IF_INTERIOR_VERTEX:
          case IF_INTERIOR_EDGE:
          case IF_INTERIOR_FACE:
            polytope2Inside = false;
            if (!polytope1Inside) {
                return;
            }
            break;
          default:
            break;
        }
    }
}


VertexIt ContactDiscoverer::findExtremePointOnPMInTheDirection(const Vec3& dir)
{
    auto& ch = mISFinder.hull3D();
    auto vPair = ch.vertices();

    double   dotMax = 0.0;
    VertexIt itMax  = vPair.second;

    for ( auto vit = vPair.first; vit != vPair.second; vit++ ) {

        auto dot = dir.dot((*vit)->pLCS());

        if ( itMax == vPair.second || dotMax < dot ) {

            itMax  = vit;
            dotMax = dot;       
        }
    }
    return itMax;
}


Mat3x3 ContactDiscoverer::rotMatAlignZDirToZAxis(const Vec3& zDir)
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


bool ContactDiscoverer::doesFaceIncludeOriginXY(
    FaceIt        fit,
    const Mat3x3& Q,
    const Vec3&   com
) {
    for ( auto he : (*fit)->halfEdges() ) {

        auto src = (*he)->src();
        auto dst = (*he)->dst();
        auto pSrc = (*src)->pGCS(Q, com);
        auto pDst = (*dst)->pGCS(Q, com);

        auto v12  = pDst - pSrc;
        auto vOri = pSrc * -1.0;

        if ( vOri.squaredNorm2() <= mEpsilonZero ) {
            return true;
        }
        else {

            v12. normalize();
            vOri.normalize();

            auto cr = v12.cross(vOri);
            if ( cr.z() <= -1.0 * mEpsilonZero ) {
                return false;
            }          
        }
    }
    return true;
}


void ContactDiscoverer::processPM_Poly1Inside()
{
    // Find the relative velocity at CoM 1.    
    const auto p21 = mCoM1 - mCoM2;
    auto       v21 = mVlin1 - ( mVlin2 - p21.cross(mVang2) );

    if ( v21.squaredNorm2() <= mEpsilonZero ) {

        log( Loggable::ERROR, __FILE__, __LINE__, 
                              "Relative velocity zero. Using diff of CoMs." );

        v21 = p21;

        if ( v21.squaredNorm2() <= mEpsilonZero ) {

            log( Loggable::ERROR, __FILE__, __LINE__, 
                        "Concentric objects. Using an arbitrary direction." );

            v21 = Vec3(1.0, 0.0, 0.0);
        }
    }

    v21.normalize();

    auto  vitPM = findExtremePointOnPMInTheDirection(v21);
    auto& attr  = mISFinder.vertexAttributes3D(vitPM);

    if ( attr.mPred == IF_VERTEX_INTERIOR ) {
        mInfo.mType1 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit1  = attr.mVit1;
    }
    else {

        log( Loggable::ERROR, __FILE__, __LINE__, 
             "Invalid 3-dim penetration interior manifold"  );

    }

    auto fit = 
          findIntersectingFaceOnBody2FromPointInDirection( mCoM1, v21 * -1.0 );
    mInfo.mType2 = ContactPairInfo::FT_FACE;
    mInfo.mFit2 = fit;
}


void ContactDiscoverer::processPM_Poly2Inside()
{
    // Find the relative velocity at CoM 2.
    const auto p12 = mCoM2 - mCoM1;
    auto       v12 = mVlin2 - ( mVlin1 - p12.cross(mVang1) );

    if ( v12.squaredNorm2() <= mEpsilonZero ) {

        log( Loggable::ERROR, __FILE__, __LINE__, 
                              "Relative velocity zero. Using diff of CoMs." );

        v12 = p12;

        if ( v12.squaredNorm2() <= mEpsilonZero ) {

            log( Loggable::ERROR, __FILE__, __LINE__, 
                        "Concentric objects. Using an arbitrary direction." );

            v12 = Vec3(1.0, 0.0, 0.0);
        }
    }

    v12.normalize();

    auto  vitPM = findExtremePointOnPMInTheDirection(v12);
    auto& attr  = mISFinder.vertexAttributes3D(vitPM);

    if ( attr.mPred == IF_INTERIOR_VERTEX ) {
        mInfo.mType2 = ContactPairInfo::FT_VERTEX;
        mInfo.mVit2  = attr.mVit2;
    }
    else {

        log( Loggable::ERROR, __FILE__, __LINE__, 
             "Invalid 3-dim penetration interior manifold"  );

    }

    auto fit = 
          findIntersectingFaceOnBody1FromPointInDirection( mCoM2, v12 * -1.0 );
    mInfo.mType1 = ContactPairInfo::FT_FACE;
    mInfo.mFit1 = fit;
}


FaceIt ContactDiscoverer::findIntersectingFaceOnBody2FromPointInDirection(
    const Vec3& pBase,
    const Vec3& vDir
) {

    auto M   = rotMatAlignZDirToZAxis( vDir );
    auto Q   = M * mQmat2;
    auto com = M * (mCoM2 - pBase);

    auto& ch    = mBody2.ConvexHull();
    auto  fPair = ch.faces();
    for ( auto fit = fPair.first; fit != fPair.second; fit++ ) {

        auto N = (*fit)->nGCS(Q);
        if (N.z() > 0.0) {
            if (doesFaceIncludeOriginXY(fit, Q, com)) {
                return fit;
            }
        }
    }

    log( Loggable::ERROR, __FILE__, __LINE__, "Could not find face." );
             
    return fPair.second;
}


FaceIt ContactDiscoverer::findIntersectingFaceOnBody1FromPointInDirection(
    const Vec3& pBase,
    const Vec3& vDir
) {

    auto M   = rotMatAlignZDirToZAxis( vDir );
    auto Q   = M * mQmat1;
    auto com = M * (mCoM1 - pBase);

    auto& ch    = mBody1.ConvexHull();
    auto  fPair = ch.faces();
    for ( auto fit = fPair.first; fit != fPair.second; fit++ ) {

        auto N = (*fit)->nGCS(Q);
        if (N.z() > 0.0) {
            if (doesFaceIncludeOriginXY(fit, Q, com)) {
                return fit;
            }
        }
    }

    log( Loggable::ERROR, __FILE__, __LINE__, "Could not find face." );
             
    return fPair.second;
}


void ContactDiscoverer::processPM_3D()
{
    Vec3   relVel12;
    long   relVelDimension;
    Vec3   primaryEigenVector;
    findRelativeVelocitiesOf1To2InPenetrationManifold(
                               relVel12, relVelDimension, primaryEigenVector );

    bool useRelVel = false;

    if (relVel12.squaredNorm2() > mEpsilonZero ) {

        relVel12.normalize();
        if ( isWithinSepAxesCone1To2(relVel12) ) {

            if (relVelDimension == 0) {
                useRelVel = true;
            }
            else if (relVelDimension == 1) {
                auto cr = relVel12.cross(primaryEigenVector);
                if (cr.squaredNorm2() <= mEpsilonZeroPCA) {
                    // Variation is in the direction of mean velocity.
                    useRelVel = true;
                }
            }
        }
    }

    Vec3 penDir1, penDir2;

    if ( useRelVel ) { 

        penDir1 = relVel12;
        penDir2 = relVel12 * -1.0;
    }
    else {
        penDir1 = findDirectionFromWeightedAverageSide1();
        penDir2 = findDirectionFromWeightedAverageSide2();
    }

    auto vit1 = findExtremeVertexSide1(penDir1);
    auto vit2 = findExtremeVertexSide2(penDir2);

    adjustFeatureAndSetContactInfoSide1(vit1, penDir1);
    adjustFeatureAndSetContactInfoSide2(vit2, penDir2);

}


bool ContactDiscoverer::isWithinSepAxesCone1To2(const Vec3& n)
{
    for (auto& ax : mSeparationAxes1To2) {
        if ( ax.dot(n) < mEpsilonZero ) {
            return false;
        }
    }
    return true;
}


bool ContactDiscoverer::isWithinSepAxesCone2To1(const Vec3& n)
{
    for (auto& ax : mSeparationAxes1To2) {
        if ( ax.dot(n) > -1.0 * mEpsilonZero ) {
            return false;
        }
    }
    return true;
}


void ContactDiscoverer::findRelativeVelocitiesOf1To2InPenetrationManifold(
    Vec3& mean,
    long& dimension,
    Vec3& primaryEigenVector
) {

    auto& ch = mISFinder.hull3D();
    auto vPair = ch.vertices();
    vector<Vec3> relVels1To2;

    for ( auto vit = vPair.first; vit != vPair.second; vit++ ) {

        auto pGCS = (*vit)->pLCS();
        auto r1   = pGCS - mCoM1;
        auto r2   = pGCS - mCoM2;
        auto v1   = mVlin1 - r1.cross(mVang1);
        auto v2   = mVlin2 - r2.cross(mVang2);

        relVels1To2.push_back(v1 - v2);

    }

    Vec3 spread;
    auto EV = findPrincipalComponents( relVels1To2, spread, mean );

    primaryEigenVector = EV.col(1);

    if (spread.z() <= mEpsilonZeroPCA) {

        if (spread.y() <= mEpsilonZeroPCA) {

            if (spread.x() <= mEpsilonZeroPCA) {

                dimension = 0;
            }
            else {

                dimension = 1;
            }
        }
        else {

            dimension = 2;
        }
    }
    else {

        dimension = 3;
    }


}


Vec3 ContactDiscoverer::findDirectionFromWeightedAverageSide1()
{

    Vec3 sum(0.0, 0.0, 0.0);

    auto& ch = mISFinder.hull3D();
    auto fPair = ch.faces();
    for ( auto fit = fPair.first; fit != fPair.second; fit++ ) {

        auto att  = mISFinder.faceAttributes3D(fit);

        if ( att.mPred == IF_FACE_INTERIOR    ) {

            auto nGCS = (*fit)->nLCS();

            if ( isWithinSepAxesCone1To2(nGCS) ) {

                auto area = (*fit)->areaIfConvex();
                sum += (nGCS * area);
            }
        }
    }

    if (sum.squaredNorm2() <= mEpsilonZero) {
        // Fallback in case there is no face.
        auto vPair = ch.vertices();
        for ( auto vit = vPair.first; vit != vPair.second; vit++ ) {

            auto att = mISFinder.vertexAttributes3D(vit);

            switch (att.mPred) {

              case IF_VERTEX_VERTEX:
              case IF_VERTEX_EDGE:
              case IF_VERTEX_FACE:
              case IF_VERTEX_INTERIOR:
              case IF_EDGE_VERTEX:
              case IF_EDGE_EDGE:
              case IF_EDGE_FACE:
              case IF_EDGE_INTERIOR:
              case IF_FACE_VERTEX:
              case IF_FACE_EDGE:
              case IF_FACE_FACE:
              case IF_FACE_INTERIOR:
                {
                    auto nGCS = (*vit)->nLCS();

                    if ( isWithinSepAxesCone1To2(nGCS) ) {
                        sum += nGCS;
                    }
                }
                break;
              default:
                break; 
            }
        }
    }

    sum.normalize();

    return sum;
}


Vec3 ContactDiscoverer::findDirectionFromWeightedAverageSide2()
{

    Vec3 sum(0.0, 0.0, 0.0);

    auto& ch = mISFinder.hull3D();
    auto fPair = ch.faces();
    for ( auto fit = fPair.first; fit != fPair.second; fit++ ) {

        auto att  = mISFinder.faceAttributes3D(fit);

        if ( att.mPred == IF_INTERIOR_FACE ) {

            auto nGCS = (*fit)->nLCS();

            if ( isWithinSepAxesCone2To1(nGCS) ) {

                auto area = (*fit)->areaIfConvex();
                sum += (nGCS * area);

            }
        }
    }

    if (sum.squaredNorm2() <= mEpsilonZero) {
        // Fallback in case there is no face.
        auto vPair = ch.vertices();
        for ( auto vit = vPair.first; vit != vPair.second; vit++ ) {

            auto att = mISFinder.vertexAttributes3D(vit);

            switch (att.mPred) {

              case IF_VERTEX_VERTEX:
              case IF_VERTEX_EDGE:
              case IF_VERTEX_FACE:
              case IF_EDGE_VERTEX:
              case IF_EDGE_EDGE:
              case IF_EDGE_FACE:
              case IF_FACE_VERTEX:
              case IF_FACE_EDGE:
              case IF_FACE_FACE:
              case IF_INTERIOR_VERTEX:
              case IF_INTERIOR_EDGE:
              case IF_INTERIOR_FACE:
                {
                    auto nGCS = (*vit)->nLCS();

                    if ( isWithinSepAxesCone2To1(nGCS) ) {
                        sum += nGCS;
                    }
                }
                break;
              default:
                break; 
            }
        }
    }

    sum.normalize();

    return sum;
}


VertexIt ContactDiscoverer::findExtremeVertexSide1(const Vec3& dir)
{

    auto&    ch     = mISFinder.hull3D();
    VertexIt maxV   = ch.vertices().second;
    double   maxDot;

    auto  vPair = ch.vertices();
    for ( auto vit = vPair.first ; vit != vPair.second ; vit++ ) {

        auto attr = mISFinder.vertexAttributes3D(vit);

        switch ( attr.mPred ) {

          case IF_VERTEX_VERTEX:
          case IF_VERTEX_EDGE:
          case IF_VERTEX_FACE:
          case IF_EDGE_VERTEX:
          case IF_EDGE_EDGE:
          case IF_EDGE_FACE:
          case IF_FACE_VERTEX:
          case IF_FACE_EDGE:
          case IF_FACE_FACE:
          case IF_VERTEX_INTERIOR:
          case IF_EDGE_INTERIOR:
          case IF_FACE_INTERIOR:
            {
                auto dot  = dir.dot( (*vit)->pLCS() );

                if ( maxV == ch.vertices().second || dot > maxDot ) {

                   maxV   = vit;
                   maxDot = dot;
                }
            }
            break;

          default:
            break;
        }
    }
    return maxV;
}


VertexIt ContactDiscoverer::findExtremeVertexSide2(const Vec3& dir)
{
    auto&    ch     = mISFinder.hull3D();
    VertexIt maxV   = ch.vertices().second;
    double   maxDot;

    auto  vPair = ch.vertices();
    for ( auto vit = vPair.first ; vit != vPair.second ; vit++ ) {

        auto attr = mISFinder.vertexAttributes3D(vit);

        switch ( attr.mPred ) {

          case IF_VERTEX_VERTEX:
          case IF_VERTEX_EDGE:
          case IF_VERTEX_FACE:
          case IF_EDGE_VERTEX:
          case IF_EDGE_EDGE:
          case IF_EDGE_FACE:
          case IF_FACE_VERTEX:
          case IF_FACE_EDGE:
          case IF_FACE_FACE:
          case IF_INTERIOR_VERTEX:
          case IF_INTERIOR_EDGE:
          case IF_INTERIOR_FACE:
            { 
                auto dot  = dir.dot( (*vit)->pLCS() );

                if ( maxV == ch.vertices().second || dot > maxDot ) {

                   maxV   = vit;
                   maxDot = dot;
                }
            }
            break;

          default:
            break;
        }
    }
    return maxV;
}


void ContactDiscoverer::adjustFeatureOnPMSide1(
    const Vec3& dir,
    VertexIt    vitBest,
    EdgeIt&     eitBest,
    FaceIt&     fitBest
) {
    auto& ch = mISFinder.hull3D();
    auto fitEnd = ch.faces().second;
    auto eitEnd = ch.edges().second;

    const Vec3 pBase = (*vitBest)->pLCS();

    fitBest = fitEnd;
    double maxFace;

    eitBest = eitEnd;
    double minEdge;

    for ( auto he : (*vitBest)->halfEdges() ) {

        if ((*he)->src() == vitBest) {

            auto dst = (*he)->dst();

            auto fit = (*he)->face();
            auto Nf  = (*fit)->nLCS();
            auto dot = dir.dot(Nf);

            auto attrf = mISFinder.faceAttributes3D(fit);
            if ( attrf.mPred == IF_FACE_INTERIOR ) {

                if ( dir.dot(Nf) >= (1.0 - mEpsilonAngle) ) {

                    if (fitBest == fitEnd || maxFace < dot) {

                        fitBest = fit;
                        maxFace = dot;

                    }
                }
            }

            auto eit = (*he)->edge();
            auto De  = (*dst)->pLCS() - pBase;
            De.normalize();
            auto fdot = fabs(dir.dot(De));

            auto attre = mISFinder.edgeAttributes3D(eit);

            switch (attre.mPred) {

              case IF_EDGE_INTERIOR:
              case IF_FACE_INTERIOR:
                if ( fdot <= mEpsilonAngle ) {

                   if (eitBest == eitEnd || minEdge > fdot) {

                        eitBest = eit;
                        minEdge = fdot;
                    }
                }
                break;
              default:
                break;

            }
        }
    }
}


void ContactDiscoverer::adjustFeatureOnPMSide2(
    const Vec3& dir,
    VertexIt    vitBest,
    EdgeIt&     eitBest,
    FaceIt&     fitBest
) {
    auto& ch = mISFinder.hull3D();
    auto fitEnd = ch.faces().second;
    auto eitEnd = ch.edges().second;

    const Vec3 pBase = (*vitBest)->pLCS();

    fitBest = fitEnd;
    double maxFace;

    eitBest = eitEnd;
    double minEdge;

    for ( auto he : (*vitBest)->halfEdges() ) {

        if ((*he)->src() == vitBest) {
            auto dst = (*he)->dst();

            auto fit = (*he)->face();
            auto Nf  = (*fit)->nLCS();
            auto dot = dir.dot(Nf);

            auto attrf = mISFinder.faceAttributes3D(fit);
            if ( attrf.mPred == IF_INTERIOR_FACE ) {
                if ( dir.dot(Nf) >= (1.0 - mEpsilonAngle) ) {

                    if (fitBest == fitEnd || maxFace < dot) {
                        fitBest = fit;
                        maxFace = dot;

                    }
                }
            }

            auto eit = (*he)->edge();
            auto De  = (*dst)->pLCS() - pBase;
            De.normalize();
            auto fdot = fabs(dir.dot(De));

            auto attre = mISFinder.edgeAttributes3D(eit);

            switch (attre.mPred) {

              case IF_INTERIOR_EDGE:
              case IF_INTERIOR_FACE:
                if ( fdot <= mEpsilonAngle ) {

                   if (eitBest == eitEnd || minEdge > fdot) {
                        eitBest = eit;
                        minEdge = fdot;
                    }
                }
                break;
              default:
                break;

            }
        }
    }
}


void ContactDiscoverer::adjustFeatureAndSetContactInfoSide1(
    VertexIt    vit,
    const Vec3& dir
) {

    EdgeIt eitBest;
    FaceIt fitBest;

    adjustFeatureOnPMSide1(dir, vit, eitBest, fitBest);

    auto& ch = mISFinder.hull3D();

    if (fitBest!=ch.faces().second) {

        auto& attr = mISFinder.faceAttributes3D(fitBest);
        switch (attr.mPred) {
          case IF_FACE_INTERIOR:
          case IF_FACE_EDGE:
          case IF_FACE_FACE:
            mInfo.mType1 = ContactPairInfo::FT_FACE;
            mInfo.mFit1  = attr.mFit1;
            break;

          default:
            log( Loggable::ERROR, __FILE__, __LINE__, "Invalid predicate");
            break;
        }

    }
    else if (eitBest!=ch.edges().second) {

        auto& attr = mISFinder.edgeAttributes3D(eitBest);

        switch (attr.mPred) {

          case IF_EDGE_INTERIOR:
          case IF_EDGE_EDGE:
          case IF_EDGE_FACE:
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = attr.mEit1;
            break;

          case IF_FACE_INTERIOR:
          case IF_FACE_EDGE:
          case IF_FACE_FACE:
            mInfo.mType1 = ContactPairInfo::FT_FACE;
            mInfo.mFit1  = attr.mFit1;
            break;

          default:
            log( Loggable::ERROR, __FILE__, __LINE__, "Invalid predicate");
            break;
        }

    }
    else {

        auto& attr = mISFinder.vertexAttributes3D(vit);

        switch (attr.mPred) {

          case IF_VERTEX_VERTEX:
          case IF_VERTEX_EDGE:
          case IF_VERTEX_FACE:
          case IF_VERTEX_INTERIOR:
            mInfo.mType1 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit1  = attr.mVit1;
            break;

          case IF_EDGE_VERTEX:
          case IF_EDGE_EDGE:
          case IF_EDGE_FACE:
          case IF_EDGE_INTERIOR:
            mInfo.mType1 = ContactPairInfo::FT_EDGE;
            mInfo.mEit1  = attr.mEit1;
            break;

          case IF_FACE_VERTEX:
          case IF_FACE_EDGE:
          case IF_FACE_FACE:
          case IF_FACE_INTERIOR:
            mInfo.mType1 = ContactPairInfo::FT_FACE;
            mInfo.mFit1  = attr.mFit1;
            break;

          default:
            log( Loggable::ERROR, __FILE__, __LINE__, "Invalid predicate");
            break;
        }
    }
}


void ContactDiscoverer::adjustFeatureAndSetContactInfoSide2(
    VertexIt    vit,
    const Vec3& dir
) {

    EdgeIt eitBest;
    FaceIt fitBest;

    adjustFeatureOnPMSide2(dir, vit, eitBest, fitBest);

    auto& ch = mISFinder.hull3D();

    if (fitBest!=ch.faces().second) {

        auto& attr = mISFinder.faceAttributes3D(fitBest);
        switch (attr.mPred) {
          case IF_INTERIOR_FACE:
          case IF_EDGE_FACE:
          case IF_FACE_FACE:
            mInfo.mType2 = ContactPairInfo::FT_FACE;
            mInfo.mFit2  = attr.mFit2;
            break;

          default:
            log( Loggable::ERROR, __FILE__, __LINE__, "Invalid predicate");
            break;
        }

    }
    else if (eitBest!=ch.edges().second) {

        auto& attr = mISFinder.edgeAttributes3D(eitBest);

        switch (attr.mPred) {

          case IF_INTERIOR_EDGE:
          case IF_EDGE_EDGE:
          case IF_FACE_EDGE:
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = attr.mEit2;
            break;

          case IF_INTERIOR_FACE:
          case IF_EDGE_FACE:
          case IF_FACE_FACE:
            mInfo.mType2 = ContactPairInfo::FT_FACE;
            mInfo.mFit2  = attr.mFit2;
            break;

          default:
            log( Loggable::ERROR, __FILE__, __LINE__, "Invalid predicate");
            break;
        }

    }
    else {

        auto& attr = mISFinder.vertexAttributes3D(vit);

        switch (attr.mPred) {

          case IF_VERTEX_VERTEX:
          case IF_EDGE_VERTEX:
          case IF_FACE_VERTEX:
          case IF_INTERIOR_VERTEX:
            mInfo.mType2 = ContactPairInfo::FT_VERTEX;
            mInfo.mVit2  = attr.mVit2;
            break;

          case IF_VERTEX_EDGE:
          case IF_EDGE_EDGE:
          case IF_FACE_EDGE:
          case IF_INTERIOR_EDGE:
            mInfo.mType2 = ContactPairInfo::FT_EDGE;
            mInfo.mEit2  = attr.mEit2;
            break;

          case IF_VERTEX_FACE:
          case IF_EDGE_FACE:
          case IF_FACE_FACE:
          case IF_INTERIOR_FACE:
            mInfo.mType2 = ContactPairInfo::FT_FACE;
            mInfo.mFit2  = attr.mFit2;
            break;

          default:
            log( Loggable::ERROR, __FILE__, __LINE__, "Invalid predicate");
            break;
        }
    }
}


void ContactDiscoverer::updateContactNormal()
{
    if ( mInfo.mType1 == ContactPairInfo::FT_FACE) {
        mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_1;
        mInfo.mContactNormal1To2 = (*mInfo.mFit1)->nGCS(mQmat1);
    }
    else if ( mInfo.mType1 == ContactPairInfo::FT_EDGE) {
        if ( mInfo.mType2 == ContactPairInfo::FT_FACE) {
            mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;
            mInfo.mContactNormal1To2 = 
                             (*mInfo.mFit2)->nGCS(mQmat2) * -1.0;
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE) {
            auto vit1src = (*((*mInfo.mEit1)->he1()))->src();
            auto vit1dst = (*((*mInfo.mEit1)->he1()))->dst();
            auto vit2src = (*((*mInfo.mEit2)->he1()))->src();
            auto vit2dst = (*((*mInfo.mEit2)->he1()))->dst();

            auto d1      = (*vit1dst)->pGCS(mQmat1) -
                           (*vit1src)->pGCS(mQmat1) ;

            auto d2      = (*vit2dst)->pGCS(mQmat2) -
                           (*vit2src)->pGCS(mQmat2) ;

            auto n1      =  (*mInfo.mEit1)->nGCS(mQmat1);
            auto n2      =  (*mInfo.mEit2)->nGCS(mQmat2);
            auto nAvg    = n1 - n2;
            nAvg.normalize();
            auto cr = d1.cross(d2);
            if (cr.squaredNorm2() <= mEpsilonAngle) {
                // two edges are parallel. Take average of the two edge normals
                mInfo.mContactNormal1To2 = nAvg;
                mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_EDGE_AVG;
            }
            else {
                cr.normalize();

                if (fabs(nAvg.dot(cr)) >= sqrt(2.0)/2.0) {

                    if (cr.dot(mCoM2 - mCoM1)< 0.0) {

                        cr.scale(-1.0);
                    }
                    mInfo.mContactNormal1To2 = cr;
                    mInfo.mContactNormalType = 
                                          ContactPairInfo::NT_EDGE_CROSS_EDGE;
                }
                else {

                    mInfo.mContactNormal1To2 = nAvg;
                    mInfo.mContactNormalType = 
                                          ContactPairInfo::NT_EDGE_EDGE_AVG;
                }
            }
        }
        else {
            mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_NORMAL_1;
            mInfo.mContactNormal1To2 = (*mInfo.mEit1)->nGCS(mQmat1);
        }
    }
    else {
        if ( mInfo.mType2 == ContactPairInfo::FT_FACE) {
            mInfo.mContactNormalType = ContactPairInfo::NT_FACE_NORMAL_2;
            mInfo.mContactNormal1To2 = 
                              (*mInfo.mFit2)->nGCS(mQmat2) * -1.0;
        }
        else if ( mInfo.mType2 == ContactPairInfo::FT_EDGE) {
            mInfo.mContactNormalType = ContactPairInfo::NT_EDGE_NORMAL_2;
            mInfo.mContactNormal1To2 = 
                              (*mInfo.mEit2)->nGCS(mQmat2) * -1.0;
        }
        else {
            mInfo.mContactNormal1To2 =
                (*mInfo.mVit1)->nGCS(mQmat1) + 
                (*mInfo.mVit2)->nGCS(mQmat2) * -1.0;
            mInfo.mContactNormal1To2.normalize();
            mInfo.mContactNormalType = ContactPairInfo::NT_VERTEX_VERTEX_AVG;
        }
    }
}



}// namespace Makena

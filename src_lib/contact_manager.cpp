#include "contact_manager.hpp"

/**
 * @file contact_manager.cpp
 */
namespace Makena {


const std::string ERR_LOCKED   = "Operation not allowed while locked";
const std::string ERR_UNLOCKED = "Operation not allowed while unlocked";

ContactManager::ContactManager(
    const double     KCorr,
    const double     velocityThreshold,
    const long       maxNumPointsPerContact,
    const double     epsilonZero,
    const double     epsilonAngle,
    const double     epsilonZeroGJK,
    const double     scalingGJK,
    const long       maxNumIter,
    const long       maxNumCycles,
    std::ostream&    logStream
):
    mKCorr                    (KCorr),
    mVelocityThresholdSquared (velocityThreshold),
    mMaxNumPointsPerContact   (maxNumPointsPerContact),
    mEpsilonZero              (epsilonZero),
    mEpsilonAngle             (epsilonAngle),
    mEpsilonZeroGJK           (epsilonZeroGJK),
    mScalingGJK               (scalingGJK),
    mMaxNumIter               (maxNumIter),
    mMaxNumCycles             (maxNumCycles),
    mLogStream                (logStream),
    mLocked                   (false)
    {;}


ContactManager::~ContactManager() 
{

    for (auto it = mActivePairs.begin(); it != mActivePairs.end(); it++) {

        auto* info = it->second;
        delete info;

    }

    for (auto it = mNewPairs.begin(); it != mNewPairs.end(); it++) {

        auto* info = it->second;
        delete info;

    }
}


void ContactManager::registerRigidBody(ConvexRigidBody* p)
{

    if (mLocked) {
        throw std::logic_error(ERR_LOCKED);
    }

    mAABBDetector.registerElem(&(p->aabb()));

}


void ContactManager::unregisterRigidBody(ConvexRigidBody* p)
{
    if (mLocked) {
        throw std::logic_error(ERR_LOCKED);
    }

    mAABBDetector.unregisterElem(&(p->aabb()));

    // Scan through the active pairs and remove elements.
    vector <CollisionMapItType>  toBeRemoved;
    for (auto it = mActivePairs.begin(); it != mActivePairs.end(); it++) {
        if (it->first.first == p || it->first.second == p) {
            toBeRemoved.push_back(it);
        }
    }

    for (auto it = mNewPairs.begin(); it != mNewPairs.end(); it++) {
        if (it->first.first == p || it->first.second == p) {
            toBeRemoved.push_back(it);
        }
    }

    for (auto it : toBeRemoved) {

        auto* info = it->second;
        delete info;
        mActivePairs.erase(it);
    }
}


void ContactManager::initializeAABB()
{

    mAABBDetector.update();

}


void ContactManager::lockForNextStep(const double deltaT)
{
    mLocked          = true;

    mDeltaTinv       = 1.0 / deltaT;

    mKCorrOverDeltaT = mKCorr * mDeltaTinv;
}


void ContactManager::unlockAfterStep(ConstraintManager& m)
{
    findContactPressure();

    mergeNewPairs();

    removeInactivePairs(m);

    cleanupJacobianConstraints(m);

    mLocked = false;
}


void ContactManager::removeInactivePairs(ConstraintManager& m)
{
    vector<CollisionMapItType> toBeRemoved;

    for (auto it = mActivePairs.begin(); it != mActivePairs.end(); it++) {

	auto* info      = it->second;

        bool  noContact = true;

        for (auto* cs : info->mUniConstraints) {

            if (cs->lambda() > 0.0) {
                noContact = false;
                break;
            }

        }

        if (noContact) {
            (info->mInactiveCount)++;
            if (info->mInactiveCount > 0) {
                toBeRemoved.push_back(it);
            }
        }
        else {
            info->mInactiveCount = 0;
        }
    }

    for (auto it : toBeRemoved) {
        auto* info = it->second;

        for (auto* cs : info->mUniConstraints) {
            m.unregisterConstraint(cs);
        }

        for (auto* cs : info->mBiConstraints) {
            m.unregisterConstraint(cs);
        }

        info->deleteConstraints();
        delete info;
        mActivePairs.erase(it);

    }
}


void ContactManager::mergeNewPairs()
{

    mActivePairs.insert(mNewPairs.begin(), mNewPairs.end());

    mNewPairs.clear();
}


void ContactManager::findContactPressure()
{
    for (auto it = mActivePairs.begin(); it != mActivePairs.end(); it++) {

	auto* info  = it->second;

        info->mSumLambdas = 0.0;

        for (auto* cs : info->mUniConstraints) {

            (info->mSumLambdas) += (cs->lambda());

        }
    }
}


void ContactManager::cleanupJacobianConstraints(ConstraintManager& m)
{

    for (auto it = mActivePairs.begin(); it != mActivePairs.end(); it++) {
	auto* info  = it->second;

        for (auto* cs : info->mUniConstraints) {
            m.unregisterConstraint(cs);
        }

        for (auto* cs : info->mBiConstraints) {
            m.unregisterConstraint(cs);
        }

        info->deleteConstraints();

    }

    for (auto it = mNewPairs.begin(); it != mNewPairs.end(); it++) {
	auto* info  = it->second;

        for (auto* cs : info->mUniConstraints) {
            m.unregisterConstraint(cs);
        }

        for (auto* cs : info->mBiConstraints) {
            m.unregisterConstraint(cs);
        }

        info->deleteConstraints();
    }
}


void ContactManager::updateActiveContacts(ConstraintManager& m)
{
    if (!mLocked) {
        throw std::logic_error(ERR_UNLOCKED);
    }

    vector<CollisionMapItType> toBeRemoved;

    for (auto it = mActivePairs.begin(); it != mActivePairs.end(); it++) {

        auto& body1 = *(it->first.first);
        auto& body2 = *(it->first.second);
	auto& info  = *(it->second);

	ContactUpdater updater( 
            body1,
            body2,
            info,
            mEpsilonZero,
            mEpsilonAngle,
            mEpsilonZeroGJK,
            mScalingGJK,
            mMaxNumPointsPerContact,
            mLogStream
        );

        auto noContact = updater.update();

        if (noContact) {
            toBeRemoved.push_back(it);
            
        }
        else {
            generateUnilateralConstraints(body1, body2, info);
            generateFriction(body1, body2, info);
        }
    }

    for (auto it : toBeRemoved) {

        auto* info = it->second;

        for (auto* cs : info->mUniConstraints) {
            m.unregisterConstraint(cs);
        }

        for (auto* cs : info->mBiConstraints) {
            m.unregisterConstraint(cs);
        }
        info->deleteConstraints();

        delete info;
        mActivePairs.erase(it);
    }
}


void ContactManager::registerActiveConstraints(ConstraintManager& m)
{
    for (auto it = mActivePairs.begin(); it != mActivePairs.end(); it++) {

	auto* info  = it->second;

        for (auto* c : info->mUniConstraints) {
            m.registerConstraint(c);
        }

        for (auto* c : info->mBiConstraints) {
            m.registerConstraint(c);
        }
    }
}
    

void ContactManager::registerNewConstraints(ConstraintManager& m)
{
    for (auto it = mNewPairs.begin(); it != mNewPairs.end(); it++) {

	auto* info  = it->second;

        for (auto* c : info->mUniConstraints) {
            m.registerConstraint(c);
        }

        for (auto* c : info->mBiConstraints) {
            m.registerConstraint(c);
        }
    }
}
    

bool ContactManager::discoverNewContacts(bool useGeomConfigTemp)
{
    if (!mLocked) {
        throw std::logic_error(ERR_UNLOCKED);
    }

    mAABBDetector.update();
    std::list< pair<ConvexRigidBody*const, 
                    ConvexRigidBody*const > >& collisionPairsAABB = 
                                                         mAABBDetector.pairs();
    for (auto& p : collisionPairsAABB) {
        
        if (mActivePairs.find(p) != mActivePairs.end()) {
            continue;
        }

        auto& body1 = *(p.first);
        auto& body2 = *(p.second);

        Vec3 sepAxisNotUsed;
        if (!Makena::doIntersect(
                body1.OBBcenter(), body1.OBBaxes(), body1.OBBextents(),
                body2.OBBcenter(), body2.OBBaxes(), body2.OBBextents(),
                1.0, sepAxisNotUsed                                     )) {
            continue;
        }

        ContactPairInfo* info = new ContactPairInfo();

        if (!discoverNewContactPair(useGeomConfigTemp, body1, body2, *info)) {
            delete info;
            continue;
        }
        generateUnilateralConstraints(body1, body2, *info);

        pair<ConvexRigidBody*const, ConvexRigidBody*const> 
                                                       newPair(&body1, &body2);
        mNewPairs[newPair] = info;
    }
    return mNewPairs.size() > 0;
}


bool ContactManager::findSeparationAxis(
    ConvexRigidBody& body1,
    ConvexRigidBody& body2,
    Vec3&            separationAxis
) {

    // Find the separation axis
    const double testStep = 0.05;

    for (double scaling = 1.0; scaling > mEpsilonZero; scaling -= testStep) {
        if (!Makena::doIntersect(
                body1.OBBcenter(), body1.OBBaxes(), body1.OBBextents(),
                body2.OBBcenter(), body2.OBBaxes(), body2.OBBextents(),
                scaling, separationAxis )) {
            return true;
            break;
        }
    }
    return false;
                                 
}



bool ContactManager::discoverNewContactPair(
    bool             useGeomConfigTemp,
    ConvexRigidBody& body1,
    ConvexRigidBody& body2,
    ContactPairInfo& info
) {
    GJKOriginFinder originFinder (
        body1,
        body2,
        info.mLastVertices,
        mMaxNumIter,
        mMaxNumCycles,
        mEpsilonZeroGJK,
        useGeomConfigTemp,
        mScalingGJK,
        mLogStream
    );

    if (!originFinder.findOrigin(true)) {
        return false;
    }
                  
    Makena::ContactDiscoverer discoverer(
        body1,
        body2, 
        EPSILON_LINEAR,
        EPSILON_LINEAR,
        0.1,
        0.05,
        0,
        true,
        std::cerr
    );

    info = discoverer.discover();

    if (info.mType1 ==ContactPairInfo::FT_NONE) {
        return false;
    }
        
    ContactPointsAndNormalGenerator pnGenerator(
        body1,
        body2,
        info,
        info.mLastVertices,
        mEpsilonZero,
        mEpsilonAngle,
        false,
        mMaxNumPointsPerContact,
        mLogStream
    );

    pnGenerator.generateContatctPointsAndNormalFromActiveFeatures();

    if (info.mIntsect1.size()==0) {
        return false;
    }

    return true;
}


void ContactManager::getPerpendicularVectors(const Vec3& n, Vec3& u ,Vec3& w)
{
    Vec3 a(0.0, 0.0, 0.0);
    double x = fabs(n.x());
    double y = fabs(n.y());
    double z = fabs(n.z());
    if (x <= y && x <= z ) {
        a.setX(1.0);
    }
    else if (y <= z && y <= x ) {
        a.setY(1.0);
    }
    else {
        a.setZ(1.0);
    }
    u = n.cross(a);
    u.normalize();
    w = n.cross(u);
    w.normalize();
}


void ContactManager::generateFriction(
    ConvexRigidBody& body1,
    ConvexRigidBody& body2,
    ContactPairInfo& info
) {

    auto muStatic  = std::min(body1.muStatic(),  body2.muStatic());
    auto muDynamic = std::min(body1.muDynamic(), body2.muDynamic());

    if (info.mSumLambdas > 0.0 && (muStatic > 0.0 || muDynamic > 0.0)) {
        Vec3 a1(0.0, 0.0, 0.0);
        for (auto& p : info.mIntsect1) {
            a1 += p;
        }
        a1.scale(1.0/info.mIntsect1.size());

        Vec3 a2(0.0, 0.0, 0.0);
        for (auto& p : info.mIntsect2) {
            a2 += p;
        }
        a2.scale(1.0/info.mIntsect2.size());

        a1 = body1.Qmat() * a1;
        a2 = body2.Qmat() * a2;
        auto r1 = a1 + body1.CoM();
        auto r2 = a2 + body2.CoM();

        auto relV = (body1.Vlin() + body1.Vang().cross(r1)) -
                    (body2.Vlin() + body2.Vang().cross(r2)) ;

        Vec3& n = info.mContactNormal1To2;

        double mu;

        if (fabs(relV.dot(n)) <= mVelocityThresholdSquared) {
            mu = muStatic;
        }
        else {
            mu = muDynamic;
        }
        if (mu==0.0) {
            return;
        }

        Vec3 u ,w, zero(0.0, 0.0, 0.0);
        getPerpendicularVectors(n, u ,w);

        auto* jc1 = new JacobianConstraint(
                            JacobianConstraint::BILATERAL_BOX, &body1, &body2);
        info.mBiConstraints.push_back(jc1);

        if (!body1.isFixed()){
            jc1->setLHSlin1(u*-1.0);
            jc1->setLHSang1((a1.crossMat() * u) *-1.0);
        }
        else {
            jc1->setLHSlin1(zero);
            jc1->setLHSang1(zero);
        }

        if (!body2.isFixed()){
            jc1->setLHSlin2(u);
            jc1->setLHSang2((a2.crossMat() * u)*  1.0);
        }
        else {
            jc1->setLHSlin2(zero);
            jc1->setLHSang2(zero);
        }
        jc1->setRHS(0.0);
        jc1->setBoxLimit(info.mSumLambdas * mu);

        auto* jc2 = new JacobianConstraint(
                            JacobianConstraint::BILATERAL_BOX, &body1, &body2);
        info.mBiConstraints.push_back(jc2);
        if (!body1.isFixed()){
            jc2->setLHSlin1(w*-1.0);
            jc2->setLHSang1((a1.crossMat() * w) *-1.0);
        }
        else {
            jc2->setLHSlin1(zero);
            jc2->setLHSang1(zero);
        }

        if (!body2.isFixed()){
            jc2->setLHSlin2(w);
            jc2->setLHSang2((a2.crossMat() * w)*  1.0);
        }
        else {
            jc2->setLHSlin2(zero);
            jc2->setLHSang2(zero);
        }
        jc2->setRHS(0.0);
        jc2->setBoxLimit(info.mSumLambdas * mu);
        if (isRotationalFrictionRequired(info)) {
            auto* jc3 = new JacobianConstraint(
                            JacobianConstraint::BILATERAL_BOX, &body1, &body2);
            info.mBiConstraints.push_back(jc3);

            if (!body1.isFixed()){
                jc3->setLHSlin1(zero);
                jc3->setLHSang1(n);
            }
            else {
                jc3->setLHSlin1(zero);
                jc3->setLHSang1(zero);
            }

            if (!body2.isFixed()){
                jc3->setLHSlin2(zero);
                jc3->setLHSang2(n*-1.0);
            }
            else {
                jc3->setLHSlin2(zero);
                jc3->setLHSang2(zero);
            }
            jc3->setRHS(0.0);
            jc3->setBoxLimit(info.mSumLambdas * mu);
        }
    }
}


bool ContactManager::isRotationalFrictionRequired(ContactPairInfo& info)
{
    return ( info.mType1 == ContactPairInfo::FT_FACE &&
             info.mType2 == ContactPairInfo::FT_FACE    ) ||
           ( info.mType1 == ContactPairInfo::FT_FACE &&
             info.mType2 == ContactPairInfo::FT_EDGE    ) ||
           ( info.mType1 == ContactPairInfo::FT_EDGE &&
             info.mType2 == ContactPairInfo::FT_FACE    ) ;
}


Vec3 ContactManager::findRelativeVelocityOfBody1RelativeToBody2(
    ConvexRigidBody& body1,
    ConvexRigidBody& body2,
    ContactPairInfo& info,
    bool             useGeomConfigTemp
) {
    vector<VertexIt> vertices1;
    vector<VertexIt> vertices2;

    decomposeQ(info.mLastVertices, vertices1, vertices2);

    FaceIt   fit1, fit2;
    EdgeIt   eit1, eit2;
    VertexIt vit1, vit2;
    enum ContactPairInfo::FeatureType type1, type2;

    findFeatureOfRigidBody(body1.ConvexHull(), vertices1,
                                     fit1, eit1, vit1, type1);

    findFeatureOfRigidBody(body2.ConvexHull(), vertices2,
                                     fit2, eit2, vit2, type2);
    Vec3 normal1(0.0, 0.0, 0.0);
    Vec3 normal2(0.0, 0.0, 0.0);

    if (type1 == ContactPairInfo::FT_FACE) {
        normal1 = useGeomConfigTemp?
               ((*fit1)->nGCS(body1.QmatTemp())):((*fit1)->nGCS(body1.Qmat()));
    }
    else if (type1 == ContactPairInfo::FT_EDGE) {
        normal1 = useGeomConfigTemp?
               ((*eit1)->nGCS(body1.QmatTemp())):((*eit1)->nGCS(body1.Qmat()));
    }
    else if (type1 == ContactPairInfo::FT_VERTEX) {
        normal1 = useGeomConfigTemp?
               ((*vit1)->nGCS(body1.QmatTemp())):((*vit1)->nGCS(body1.Qmat()));
    }

    if (type2 == ContactPairInfo::FT_FACE) {
        normal2 = useGeomConfigTemp?
               ((*fit2)->nGCS(body2.QmatTemp())):((*fit2)->nGCS(body2.Qmat()));
    }
    else if (type2 == ContactPairInfo::FT_EDGE) {
        normal2 = useGeomConfigTemp?
               ((*eit2)->nGCS(body2.QmatTemp())):((*eit2)->nGCS(body2.Qmat()));
    }
    else if (type2 == ContactPairInfo::FT_VERTEX) {
        normal2 = useGeomConfigTemp?
               ((*vit2)->nGCS(body2.QmatTemp())):((*vit2)->nGCS(body2.Qmat()));
    }

    if ( type1 != ContactPairInfo::FT_NONE ||
         type2 != ContactPairInfo::FT_NONE   ) {
        Vec3 N = normal1 - normal2;
        N.normalize();
        return N;
    }

    Vec3 com12 = body2.CoM() - body1.CoM();
    com12.normalize();

    Vec3 relV12(0.0, 0.0, 0.0);
    bool found = false;
    for (auto vit1 : vertices1) {
        Vec3 r1    = ((*vit1)->pLCS());
        Vec3 relV1 = useGeomConfigTemp?
                         (body1.VlinTemp() + body1.VangTemp().cross(r1)):
                         (body1.Vlin() + body1.Vang().cross(r1));

        for (auto vit2 : vertices2) {
            Vec3 r2    = ((*vit2)->pLCS());
            Vec3 relV2 = useGeomConfigTemp?
                             (body2.VlinTemp() + body2.VangTemp().cross(r2)):
                             (body2.Vlin() + body2.Vang().cross(r2));
            Vec3 relV12elem = relV1 - relV2;
            relV12elem.normalize();
            if (com12.dot(relV12elem) >= mEpsilonAngle) {
                relV12 += relV12elem;
                found = true;
            }
        }
    }
    if (found) {

        relV12.normalize();
        return relV12;
    }
    else {

        return com12;
    }

}


void ContactManager::findFeatureOfRigidBody(
    Manifold&                          ch,
    vector<VertexIt>&                  vits,
    FaceIt&                            fit,
    EdgeIt&                            eit,
    VertexIt&                          vit,
    enum ContactPairInfo::FeatureType& type
) {
    auto  noFace   = ch.faces().second;
    auto  noEdge   = ch.edges().second;
    auto  noVertex = ch.vertices().second;

    fit  = noFace;
    eit  = noEdge;
    vit  = noVertex;
    type = ContactPairInfo::FT_NONE;

    if (vits.size() == 3) {
        fit = ch.findFace(vits[0], vits[1], vits[2]);
        if (fit != noFace) {
            type = ContactPairInfo::FT_FACE;
        }
    }
    else if (vits.size() == 2) {
        eit = ch.findEdge(vits[0],vits[1]);
        if (eit != noEdge) {
            type = ContactPairInfo::FT_EDGE;
        }
        else {
            fit = ch.findFace(vits[0],vits[1]);
            if (fit != noFace) {
                type = ContactPairInfo::FT_FACE;
            }
        }
    }
    else if (vits.size() == 1) {
        vit = vits[0];
        type = ContactPairInfo::FT_VERTEX;
    }
}


void ContactManager::decomposeQ(
    vector<BDVertex>& Q,
    vector<VertexIt>& vertices1,
    vector<VertexIt>& vertices2
) {

    vertices1.clear();
    vertices2.clear();

    for (size_t i = 0; i < Q.size(); i++) {

        bool found = false;

        for (size_t j = 0; j < i; j++) {
            if (Q[j].v1() == Q[i].v1()) {
                found = true;
                break;
            }
        }

        if (!found) {
            vertices1.push_back(Q[i].v1());
        }

	found = false;
        for (size_t j = 0; j < i; j++) {
            if (Q[j].v2() == Q[i].v2()) {
	      found = true;
                break;
            }
        }

        if (!found) {
            vertices2.push_back(Q[i].v2());
        }
    }
}


void ContactManager::generateUnilateralConstraints(
    ConvexRigidBody& body1,
    ConvexRigidBody& body2,
    ContactPairInfo& info
) {
    auto numPairs  = info.mIntsect1reduced.size();
    auto&      n12 = info.mContactNormal1To2;
    const auto n21 = n12*-1.0;

    for (long i = 0 ; i < numPairs; i++) {

        auto* jc = new JacobianConstraint(
                               JacobianConstraint::UNILATERAL, &body1, &body2);
        info.mUniConstraints.push_back(jc);

        const auto& p1  = info.mIntsect1reduced[i];
        const auto& p2  = info.mIntsect2reduced[i];

        const auto a1 = body1.Qmat() * p1;
        const auto r1 = body1.CoM() + a1;
        const auto a2 = body2.Qmat() * p2;
        const auto r2 = body2.CoM() + a2;

        if (body1.isFixed()) {
            jc->setLHSlin1(0.0, 0.0, 0.0);
            jc->setLHSang1(0.0, 0.0, 0.0);
        }
        else {
            jc->setLHSlin1(n21);
            jc->setLHSang1(a1.crossMat() * n21);
        }

        if (body2.isFixed()) {
            jc->setLHSlin2(0.0, 0.0, 0.0);
            jc->setLHSang2(0.0, 0.0, 0.0);
        }
        else {
            jc->setLHSlin2(n12);
            jc->setLHSang2(a2.crossMat() * n12);
        }

        jc->setRHS(n21.dot(r1 - r2) * -1.0 * mKCorrOverDeltaT);
    }
}


}// namespace Makena

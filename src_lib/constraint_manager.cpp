#include "constraint_manager.hpp"

/**
 * @file constraint_manager.cpp
 */

namespace Makena {

const char* ConstraintManager::StateManager::EXCEPTION_MSG_01 =
                             "ConstraintManager::StateManager Invalid State";
const char* ConstraintManager::EXCEPTION_MSG_01 =
                             "ConstraintManager Invalid State";

void ConstraintManager::initializeStep(const double deltaT)
{

    mDeltaT = deltaT;

    long numCols = mStateMgr.numConstraintsNotCheckedIn();

    mMemMgr.setSizeAndInit(numCols);

    mState = IN_SIMULATION_STEP;
}


void ConstraintManager::terminateStep()
{

    mState = OUT_OF_SIMULATION_STEP;

    mStateMgr.resetConstraints();

    mColMgr.reset();
}


void ConstraintManager::update()
{

    list<ConvexRigidBody*> bodiesResetCopy(mStateMgr.bodiesReset());

    for (auto* p : bodiesResetCopy) {

        p->setU(mDeltaT);
        mStateMgr.moveBodyFromResetToNeedsUpdate(p);

    }

    for (auto* c : mStateMgr.constraintsNotCheckedIn()) {

        c->setMatIndex(mColMgr.allocateOne());

    }

    auto numCols = mColMgr.highestIndex();

    mMemMgr.adjustSize(numCols);

    SymMat M(numCols, false, false);
    M.receiveMemory(mMemMgr.M(), false);

    VarVec q(numCols, false, false);
    q.receiveMemory(mMemMgr.q(), false);

    VarVec z(numCols, false, false);
    z.receiveMemory(mMemMgr.z(), false);

    VarVec z_low(numCols, false, false);
    z_low.receiveMemory(mMemMgr.z_low(), false);

    VarVec z_high(numCols, false, false);
    z_high.receiveMemory(mMemMgr.z_high(), false);

    MLCP::WorkMemory wm(numCols, false, false);
    wm.receiveMemory(mMemMgr.WM(), true);

    for (auto* c : mStateMgr.constraintsNotCheckedIn()) {

        q.val(c->matIndex()) = -1.0 * c->RHS();

        if (c->type()==JacobianConstraint::UNILATERAL) {

            mMemMgr.setMLCPtype(c->matIndex(), 
                                    MLCP::UNILATERAL_NON_ACTIVE);
        }
        else if (c->type()==JacobianConstraint::BILATERAL_FREE) {

            mMemMgr.setMLCPtype(c->matIndex(), MLCP::BILATERAL_FREE);
        }
        else if (c->type()==JacobianConstraint::BILATERAL_BOX) {

            mMemMgr.setMLCPtype( c->matIndex(),
                                 MLCP::BILATERAL_BOXED_NOT_CLAMPED );
            z.val(c->matIndex())      = 0.0;
            z_low.val(c->matIndex())  = -1.0 * c->boxLimit();
            z_high.val(c->matIndex()) = c->boxLimit();

        }
    }

    // Extend M.
    for (auto* p : mStateMgr.bodiesNeedsUpdate()) {

        for (auto citit_i = p->constraintsNotCheckedIn().begin();
                  citit_i != p->constraintsNotCheckedIn().end(); citit_i++) {

            auto cit_i = *citit_i;
            auto i     = cit_i->matIndex();

            const Vec3& LHSlini = (cit_i->body1()==p)?
                                                 (cit_i->LHSlin1()):
                                                 (cit_i->LHSlin2());

            const Vec3& LHSangi = (cit_i->body1()==p)?
                                                 (cit_i->LHSang1()):
                                                 (cit_i->LHSang2());
            q.val(i) = q.const_val(i) + 
                       p->Ulin().dot(LHSlini) + p->Uang().dot(LHSangi);

            Vec3 MJ = p->Iinv() * LHSangi;

            for ( auto cit_j : p->constraintsCheckedIn() ) {

                auto j = cit_j->matIndex();

                const Vec3& LHSlinj = (cit_j->body1()==p)?
                                                     (cit_j->LHSlin1()):
                                                     (cit_j->LHSlin2());

                const Vec3& LHSangj = (cit_j->body1()==p)?
                                                     (cit_j->LHSang1()):
                                                     (cit_j->LHSang2());

                M.val(i,j) = M.const_val(i,j) +
                             LHSlini.dot(LHSlinj) * p->mInv() + 
                             MJ.dot(LHSangj);
            }

            for ( auto citit_j = citit_i; 
                  citit_j != p->constraintsNotCheckedIn().end(); citit_j++ ) {

                auto cit_j = *citit_j;
                auto j     = cit_j->matIndex();

                const Vec3& LHSlinj = (cit_j->body1()==p)?
                                                     (cit_j->LHSlin1()):
                                                     (cit_j->LHSlin2());

                const Vec3& LHSangj = (cit_j->body1()==p)?
                                                     (cit_j->LHSang1()):
                                                     (cit_j->LHSang2());

                M.val(i,j) = M.const_val(i,j) +
                             LHSlini.dot(LHSlinj) * p->mInv() + 
                             MJ.dot(LHSangj);
            }
        }
    }

    // Add CFM to diagonal elements
    for (auto* c : mStateMgr.constraintsNotCheckedIn()) {
        auto i = c->matIndex();
        M.val(i,i) += mMLCP_cfm;
    }

    mStateMgr.checkInConstraints();

    // Perform MLCP
    MLCP::solve(  M, q, z, z_low, z_high, numCols, mMemMgr.MLCPtypes(),
                  mMLCP_maxIter, mMLCP_kgs, mMLCP_ksm, mMLCP_tolerance, wm );

    
    // Put lambda to the constraints and 
    // calculate induced forces and torques.
    for (auto* body : mStateMgr.bodiesManaged()) {

        body->resetInducedForceTorque();

    }

    double maxZ = 0.0;
    for (auto* c : mStateMgr.constraintsCheckedIn()) {
        maxZ = std::max(maxZ, fabs(z.val(c->matIndex())));
    }

//    cerr << "lamdas: \n";
    for (auto* c : mStateMgr.constraintsCheckedIn()) {
/*
        if (c->type()==JacobianConstraint::BILATERAL_BOX) {
           cerr << "box lambda[" << c->matIndex() << "]: " <<
                                   z.val(c->matIndex()) << "\n";
        }

        if (c->type()==JacobianConstraint::BILATERAL_FREE) {
           cerr << "free lambda[" << c->matIndex() << "]: " <<
                                   z.val(c->matIndex()) << "\n";
        }

        if (c->type()==JacobianConstraint::UNILATERAL) {
           cerr << "uni lambda[" << c->matIndex() << "]: " <<
                                   z.val(c->matIndex()) << "\n";
        }
*/
//        if (z.val(c->matIndex())>10.0) {
//            z.val(c->matIndex()) = 10.0;
//        }

        c->setLambda(z.val(c->matIndex()));

        if (!isBody1Fixed(c)) {
            c->body1()->addInducedForceTorque(
                             c->LHSlin1(), c->LHSang1(), c->lambda());
        }

        if (!isBody2Fixed(c)) {
            c->body2()->addInducedForceTorque(
                             c->LHSlin2(), c->LHSang2(), c->lambda());
        }
    }


    // Update Geometric configuration, and determine normal force.
    for (auto* body : mStateMgr.bodiesManaged()) {

        body->updateGeomConfig(mDeltaT);
    }

    for (auto* body : mStateMgr.bodiesFree()) {

        body->resetInducedForceTorque();            
        body->updateGeomConfig(mDeltaT);
    }

//    printDebug();

}


void ConstraintManager::commitUpdate() {

    for (auto* body : mStateMgr.bodiesManaged()) {

        body->commitGeomConfig();
    }

    for (auto* body : mStateMgr.bodiesFree()) {

        body->commitGeomConfig();
    }

}


#ifdef UNIT_TESTS
void ConstraintManager::printDebug() 
{
    std::cerr << "UNILATERAL\n";
    for (auto* c : mStateMgr.constraintsCheckedIn()) {
        if (c->type()==JacobianConstraint::UNILATERAL) {
            std::cerr <<
                "Index:["    << c->matIndex()    << "]\t"
                 << "Lin1: " << c->LHSlin1().x() << "\t"
                             << c->LHSlin1().y() << "\t"
                             << c->LHSlin1().z() << "\t"
                 << "Ang1: " << c->LHSang1().x() << "\t"
                             << c->LHSang1().y() << "\t"
                             << c->LHSang1().z() << "\t"
                 << "Lin2: " << c->LHSlin2().x() << "\t"
                             << c->LHSlin2().y() << "\t"
                             << c->LHSlin2().z() << "\t"
                 << "Ang2: " << c->LHSang2().x() << "\t"
                             << c->LHSang2().y() << "\t"
                             << c->LHSang2().z() << "\t"
                 << "RHS: "  << c->RHS() << "\t"
                 << "Lambda: " << c->lambda() << "\n";
        }
    }

    std::cerr << "BILATERAL BOX\n";
    for (auto* c : mStateMgr.constraintsCheckedIn()) {

        if (c->type()==JacobianConstraint::BILATERAL_BOX) {
            std::cerr <<
                "Index:["    << c->matIndex()    << "]\t"
                 << "Lin1: " << c->LHSlin1().x() << "\t"
                             << c->LHSlin1().y() << "\t"
                             << c->LHSlin1().z() << "\t"
                 << "Ang1: " << c->LHSang1().x() << "\t"
                             << c->LHSang1().y() << "\t"
                             << c->LHSang1().z() << "\t"
                 << "Lin2: " << c->LHSlin2().x() << "\t"
                             << c->LHSlin2().y() << "\t"
                             << c->LHSlin2().z() << "\t"
                 << "Ang2: " << c->LHSang2().x() << "\t"
                             << c->LHSang2().y() << "\t"
                             << c->LHSang2().z() << "\t"
                 << "RHS: "  << c->RHS() << "\t"
                 << "Lambda: " << c->lambda() << "\n";
        }
    }

    std::cerr << "BILATERAL FREE\n";
    for (auto* c : mStateMgr.constraintsCheckedIn()) {

        if (c->type()==JacobianConstraint::BILATERAL_FREE) {
            std::cerr <<
                "Index:["    << c->matIndex()    << "]\t"
                 << "Lin1: " << c->LHSlin1().x() << "\t"
                             << c->LHSlin1().y() << "\t"
                             << c->LHSlin1().z() << "\t"
                 << "Ang1: " << c->LHSang1().x() << "\t"
                             << c->LHSang1().y() << "\t"
                             << c->LHSang1().z() << "\t"
                 << "Lin2: " << c->LHSlin2().x() << "\t"
                             << c->LHSlin2().y() << "\t"
                             << c->LHSlin2().z() << "\t"
                 << "Ang2: " << c->LHSang2().x() << "\t"
                             << c->LHSang2().y() << "\t"
                             << c->LHSang2().z() << "\t"
                 << "RHS: "  << c->RHS() << "\t"
                 << "Lambda: " << c->lambda() << "\n";
        }
    }

    std::cerr << "BODIES\n";    
    for (auto* body : mStateMgr.bodiesManaged()) {
        std::cerr << 
               "Index:[" << body->id() << "]\n"
            << "Fext: "  << body->Fext().x()  << "\t"
                         << body->Fext().y()  << "\t"
                         << body->Fext().z()  << "\t"
            << "Text: "  << body->Text().x()  << "\t"
                         << body->Text().y()  << "\t"
                         << body->Text().z()  << "\n"
            << "Fint: "  << body->FintDt().x()/mDeltaT  << "\t"
                         << body->FintDt().y()/mDeltaT  << "\t"
                         << body->FintDt().z()/mDeltaT  << "\t"
            << "Tint: "  << body->TintDt().x()/mDeltaT  << "\t"
                         << body->TintDt().y()/mDeltaT  << "\t"
                         << body->TintDt().z()/mDeltaT  << "\n"
            << "Vlin: "  << body->Vlin().x() << "\t"
                         << body->Vlin().y()  << "\t"
                         << body->Vlin().z()  << "\t"
            << "Vang: "  << body->Vang().x() << "\t"
                         << body->Vang().y()  << "\t"
                         << body->Vang().z()  << "\n"
            << "CoM:"     << body->CoM().x()   << "\t"
                         << body->CoM().y()   << "\t"
                         << body->CoM().z()   << "\t"
            << "q:"      << body->Q().s()   << "\t"
                         << body->Q().x()   << "\t"
                         << body->Q().y()   << "\t"
                         << body->Q().z()   << "\t"
            << "\n\n";
    }


}
#endif



}// namespace Makena

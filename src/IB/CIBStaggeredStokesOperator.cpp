// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/CIBStaggeredStokesOperator.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/PETScSAMRAIVectorReal.h"

#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Timers.
static Timer* t_apply;
static Timer* t_initialize_operator_state;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CIBStaggeredStokesOperator::CIBStaggeredStokesOperator(std::string object_name,
                                                       Pointer<CIBStrategy> cib_strategy,
                                                       bool homogeneous_bc)
    : StaggeredStokesOperator(std::move(object_name), homogeneous_bc), d_cib_strategy(cib_strategy)
{
    // Setup Timers.
    IBAMR_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBAMR::CIBStaggeredStokesOperator::apply()");
                  t_initialize_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::CIBStaggeredStokesOperator::initializeOperatorState()"););
} // CIBStaggeredStokesOperator

void
CIBStaggeredStokesOperator::setInterpScaleFactor(const double beta)
{
    d_scale_interp = beta;
} // setInterpScaleFactor

void
CIBStaggeredStokesOperator::setSpreadScaleFactor(const double gamma)
{
    d_scale_spread = gamma;
} // setSpreadScaleFactor

void
CIBStaggeredStokesOperator::setRegularizeMobilityFactor(const double delta)
{
    d_reg_mob_factor = delta;
} // setRegularizeMobilityFactor

void
CIBStaggeredStokesOperator::setNormalizeSpreadForce(const bool normalize_force)
{
    d_normalize_spread_force = normalize_force;
} // setNormalizeSpreadForce

CIBStaggeredStokesOperator::~CIBStaggeredStokesOperator()
{
    deallocateOperatorState();
} // ~CIBStaggeredStokesOperator

void
CIBStaggeredStokesOperator::apply(Vec x, Vec y)
{
    IBAMR_TIMER_START(t_apply);

    const double half_time = 0.5 * (d_new_time + d_current_time);
    Pointer<IBStrategy> ib_method_ops = d_cib_strategy;

    // Get some vectors and unpack them.
    Vec *vx, *vy;
    VecNestGetSubVecs(x, nullptr, &vx);
    VecNestGetSubVecs(y, nullptr, &vy);
    Pointer<SAMRAIVectorReal<NDIM, double> > vx0, vy0;
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVectorRead(vx[0], &vx0);
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vy[0], &vy0);
    SAMRAIVectorReal<NDIM, double>& u_p = *vx0;
    Vec L = vx[1];
    Vec U = vx[2];
    SAMRAIVectorReal<NDIM, double>& g_f = *vy0;
    Vec V = vy[1];
    Vec F = vy[2];

    // Temporary vectors.
    Vec Vrigid;
    VecDuplicate(V, &Vrigid);

    // Get the Eulerian vector components.
    const int U_idx = u_p.getComponentDescriptorIndex(0);
    const int P_idx = u_p.getComponentDescriptorIndex(1);
    const int A_U_idx = g_f.getComponentDescriptorIndex(0);
    const int A_P_idx = g_f.getComponentDescriptorIndex(1);
    const int U_scratch_idx = d_x->getComponentDescriptorIndex(0);

    Pointer<SideVariable<NDIM, double> > U_sc_var = u_p.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > P_cc_var = u_p.getComponentVariable(1);
    Pointer<SideVariable<NDIM, double> > A_U_sc_var = g_f.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > A_P_cc_var = g_f.getComponentVariable(1);

    // Simultaneously fill ghost cell values for u and p.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(2);
    transaction_comps[0] = InterpolationTransactionComponent(U_scratch_idx,
                                                             U_idx,
                                                             DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_U_bc_coefs,
                                                             d_U_fill_pattern);
    transaction_comps[1] = InterpolationTransactionComponent(P_idx,
                                                             DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_P_bc_coef,
                                                             d_P_fill_pattern);
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
        d_U_bc_coefs, d_P_bc_coef, U_scratch_idx, P_idx, d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    // Compute the action of the operator:
    // A*[u;p;U;L] :=
    //     [A_u;A_p;A_U;A_L] = [(C*I+D*L)*u + Grad P - gamma*S L; -Div u; T L;
    //                          -beta*J u + beta*T^{*} U
    //                          -beta*delta*Reg*L]

    // (a) Momentum equation.
    d_hier_math_ops->grad(A_U_idx, A_U_sc_var, /*cf_bdry_synch*/ false, 1.0, P_idx, P_cc_var, d_no_fill, half_time);
    d_hier_math_ops->laplace(A_U_idx,
                             A_U_sc_var,
                             d_U_problem_coefs,
                             U_scratch_idx,
                             U_sc_var,
                             d_no_fill,
                             half_time,
                             1.0,
                             A_U_idx,
                             A_U_sc_var);

    d_cib_strategy->setConstraintForce(L, half_time, -1.0 * d_scale_spread);
    ib_method_ops->spreadForce(A_U_idx, nullptr, std::vector<Pointer<RefineSchedule<NDIM> > >(), half_time);
    if (d_normalize_spread_force)
    {
        d_cib_strategy->subtractMeanConstraintForce(L, A_U_idx, -1 * d_scale_spread);
    }

    // (b) Divergence-free constraint.
    d_hier_math_ops->div(A_P_idx,
                         A_P_cc_var,
                         -1.0,
                         U_scratch_idx,
                         U_sc_var,
                         d_no_fill,
                         half_time,
                         /*cf_bdry_synch*/ true);
    d_bc_helper->copyDataAtDirichletBoundaries(A_U_idx, U_scratch_idx);

    // (c) Rigid body velocity constraint.
    d_cib_strategy->setInterpolatedVelocityVector(V, half_time);
    ib_method_ops->interpolateVelocity(U_scratch_idx,
                                       std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
                                       std::vector<Pointer<RefineSchedule<NDIM> > >(),
                                       half_time);

    d_cib_strategy->getInterpolatedVelocity(V, half_time, d_scale_interp);
    VecSet(Vrigid, 0.0);
    d_cib_strategy->setRigidBodyVelocity(U,
                                         Vrigid,
                                         /*only_free_dofs*/ true,
                                         /*only_imposed_dofs*/ false);

    VecScale(Vrigid, d_scale_interp);
    VecAYPX(V, -1.0, Vrigid);
    if (!MathUtilities<double>::equalEps(d_reg_mob_factor, 0.0))
    {
        d_cib_strategy->computeMobilityRegularization(Vrigid, L);
        VecAXPY(V, -1.0 * d_scale_interp * d_reg_mob_factor, Vrigid);
    }

    // (d) Force and torque constraint.
    d_cib_strategy->computeNetRigidGeneralizedForce(L,
                                                    F,
                                                    /*only_free_dofs*/ true,
                                                    /*only_imposed_dofs*/ false);
    // Delete temporary vectors.
    VecDestroy(&Vrigid);

    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(vx[0], &vx0);
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVector(vy[0], &vy0);

    IBAMR_TIMER_STOP(t_apply);
} // apply

void
CIBStaggeredStokesOperator::modifyRhsForBcs(Vec y)
{
    const double half_time = 0.5 * (d_new_time + d_current_time);
    Pointer<IBStrategy> ib_method_ops = d_cib_strategy;

    // Get vectors corresponding to fluid and Lagrangian velocity.
    Vec* vy;
    VecNestGetSubVecs(y, nullptr, &vy);
    Pointer<SAMRAIVectorReal<NDIM, double> > vy0;
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vy[0], &vy0);
    SAMRAIVectorReal<NDIM, double>& b = *vy0;
    Vec W = vy[1];

    // Modify RHS for fluid Bcs.
    modifyRhsForBcs(b);

    // Modify RHS for the action of interpolation operator J on 0.
    if (!d_homogeneous_bc)
    {
        Vec V;
        VecDuplicate(W, &V);
        Pointer<SAMRAIVectorReal<NDIM, double> > x = b.cloneVector("");
        x->allocateVectorData();
        x->setToScalar(0.0);
        const int U_idx = x->getComponentDescriptorIndex(0);
        const int P_idx = x->getComponentDescriptorIndex(1);
        StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
            d_U_bc_coefs, d_P_bc_coef, U_idx, P_idx, d_homogeneous_bc);
        d_bc_helper->enforceNormalVelocityBoundaryConditions(U_idx, P_idx, d_U_bc_coefs, d_new_time, d_homogeneous_bc);
        StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);

        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> U_transaction_comps(1);
        U_transaction_comps[0] = InterpolationTransactionComponent(U_idx,
                                                                   DATA_REFINE_TYPE,
                                                                   USE_CF_INTERPOLATION,
                                                                   DATA_COARSEN_TYPE,
                                                                   BDRY_EXTRAP_TYPE,
                                                                   CONSISTENT_TYPE_2_BDRY,
                                                                   d_U_bc_coefs,
                                                                   d_U_fill_pattern);
        Pointer<HierarchyGhostCellInterpolation> U_bdry_fill = new IBTK::HierarchyGhostCellInterpolation();
        U_bdry_fill->initializeOperatorState(U_transaction_comps, x->getPatchHierarchy());
        U_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
        U_bdry_fill->fillData(d_solution_time);

        d_cib_strategy->setInterpolatedVelocityVector(V, half_time);
        ib_method_ops->interpolateVelocity(U_idx,
                                           std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
                                           std::vector<Pointer<RefineSchedule<NDIM> > >(),
                                           half_time);

        d_cib_strategy->getInterpolatedVelocity(V, half_time, -1.0 * d_scale_interp);
        VecAXPY(W, -1.0, V);

        // Deallocate scratch data.
        x->freeVectorComponents();
        VecDestroy(&V);
    }
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVector(vy[0], &vy0);
} // modifyRhsForBcs

void
CIBStaggeredStokesOperator::imposeSolBcs(Vec x)
{
    Vec* vx;
    VecNestGetSubVecs(x, nullptr, &vx);
    Pointer<SAMRAIVectorReal<NDIM, double> > vx0;
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vx[0], &vx0);
    SAMRAIVectorReal<NDIM, double>& u_p = *vx0;
    imposeSolBcs(u_p);
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVector(vx[0], &vx0);
} // imposeSolBcs

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

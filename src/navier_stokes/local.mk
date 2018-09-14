srcs-core.cpp += $(call thisdir, \
	INSCollocatedCenteredConvectiveOperator.cpp \
	INSCollocatedConvectiveOperatorManager.cpp \
	INSCollocatedHierarchyIntegrator.cpp \
	INSCollocatedPPMConvectiveOperator.cpp \
	INSCollocatedVelocityBcCoef.cpp \
	INSCollocatedWavePropConvectiveOperator.cpp \
	INSHierarchyIntegrator.cpp \
	INSIntermediateVelocityBcCoef.cpp \
	INSProjectionBcCoef.cpp \
	INSStaggeredCenteredConvectiveOperator.cpp \
	INSStaggeredConvectiveOperatorManager.cpp \
	INSStaggeredHierarchyIntegrator.cpp \
	INSStaggeredPPMConvectiveOperator.cpp \
	INSStaggeredPressureBcCoef.cpp \
	INSStaggeredStabilizedPPMConvectiveOperator.cpp \
	INSStaggeredStochasticForcing.cpp \
	INSStaggeredUpwindConvectiveOperator.cpp \
	INSStaggeredVelocityBcCoef.cpp \
	INSStaggeredWavePropConvectiveOperator.cpp \
	INSVCStaggeredConservativeHierarchyIntegrator.cpp \
	INSVCStaggeredConservativeMassMomentumIntegrator.cpp \
	INSVCStaggeredHierarchyIntegrator.cpp \
	INSVCStaggeredNonConservativeHierarchyIntegrator.cpp \
	INSVCStaggeredPressureBcCoef.cpp \
	INSVCStaggeredVelocityBcCoef.cpp \
	KrylovLinearSolverStaggeredStokesSolverInterface.cpp \
	PETScKrylovStaggeredStokesSolver.cpp \
	SpongeLayerForceFunction.cpp \
	StaggeredStokesBlockFactorizationPreconditioner.cpp \
	StaggeredStokesBlockPreconditioner.cpp \
	StaggeredStokesBoxRelaxationFACOperator.cpp \
	StaggeredStokesFACPreconditioner.cpp \
	StaggeredStokesFACPreconditionerStrategy.cpp \
	StaggeredStokesLevelRelaxationFACOperator.cpp \
	StaggeredStokesOpenBoundaryStabilizer.cpp \
	StaggeredStokesOperator.cpp \
	StaggeredStokesPETScLevelSolver.cpp \
	StaggeredStokesPETScMatUtilities.cpp \
	StaggeredStokesPETScVecUtilities.cpp \
	StaggeredStokesPhysicalBoundaryHelper.cpp \
	StaggeredStokesProjectionPreconditioner.cpp \
	StaggeredStokesSolver.cpp \
	StaggeredStokesSolverManager.cpp \
	StokesBcCoefStrategy.cpp \
	StokesSpecifications.cpp \
	SurfaceTensionForceFunction.cpp \
	VCStaggeredStokesOperator.cpp \
	VCStaggeredStokesProjectionPreconditioner.cpp \
	)

include $(call incsubdirs,fortran)

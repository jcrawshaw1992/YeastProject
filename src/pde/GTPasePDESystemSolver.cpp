#include "GTPasePDESystemSolver.hpp"
#include "DistributedVectorFactory.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::VELOCITY[][SPACE_DIM] = {{ActinGrowthRate * cos(ANGLES[0]), ActinGrowthRate * sin(ANGLES[0])},
																				     {ActinGrowthRate * cos(ANGLES[1]), ActinGrowthRate * sin(ANGLES[1])},
																				     {ActinGrowthRate * cos(ANGLES[2]), ActinGrowthRate * sin(ANGLES[2])},
																				     {ActinGrowthRate * cos(ANGLES[3]), ActinGrowthRate * sin(ANGLES[3])},
																				     {ActinGrowthRate * cos(ANGLES[4]), ActinGrowthRate * sin(ANGLES[4])}};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActiveCDC42SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double active_cdc42 = rU[CDC42];
    double inactive_cdc42 = rU[INACTIVE_CDC42];

    return CDC42ActivationRate(rU, rX, t, pNode) * (inactive_cdc42 / TotalCDC42) - CDC42DecayRate * active_cdc42;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActiveRacSourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double active_rac = rU[RAC];
    double inactive_rac = rU[INACTIVE_RAC];

    return RacActivationRate(rU) * (inactive_rac / TotalRac) - RacDecayRate * active_rac;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActiveRhoSourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double active_rho = rU[RHO];
    double inactive_rho = rU[INACTIVE_RHO];

    return RhoActivationRate(rU) * (inactive_rho / TotalRho) - RhoDecayRate * active_rho;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::InactiveCDC42SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double active_cdc42 = rU[CDC42];
    double inactive_cdc42 = rU[INACTIVE_CDC42];

    return - CDC42ActivationRate(rU, rX, t, pNode) * (inactive_cdc42 / TotalCDC42) + CDC42DecayRate * active_cdc42;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::InactiveRacSourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double active_rac = rU[RAC];
    double inactive_rac = rU[INACTIVE_RAC];

    return - RacActivationRate(rU) * (inactive_rac / TotalRac) + RacDecayRate * active_rac;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::InactiveRhoSourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double active_rho = rU[RHO];
    double inactive_rho = rU[INACTIVE_RHO];

    return - RhoActivationRate(rU) * (inactive_rho / TotalRho) + RhoDecayRate * active_rho;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActivePIP1SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double active_PIP1 = rU[PIP1];
    double active_PIP2 = rU[PIP2];
    double active_rac= rU[RAC];

    return  PIP1ImputCoefficient - (PIP1DecayCoefficient * active_PIP1)  + (PIP2toPIP1BaseConversionCoefficient * active_PIP2) - PIP1toPIP2BaseConversionCoefficient/2 * ( 1 + active_rac / BasalRac ) * active_PIP1 ;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActivePIP2SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double active_PIP1 = rU[PIP1];
    double active_PIP2 = rU[PIP2];
    double active_PIP3 = rU[PIP3];
    double active_rac  = rU[RAC];
    double active_rho  = rU[RHO];

    return -PIP2toPIP1BaseConversionCoefficient * active_PIP2 + PIP1toPIP2BaseConversionCoefficient/2 * (1 + active_rac / BasalRac ) * active_PIP1 - PIP2toPIP3BaseConversionCoefficient /2 * (1+ active_rac / BasalRac) * active_PIP2 + PIP3toPIP2BaseConversionCoefficient / 2 * (1+ active_rho / BasalRho) * active_PIP3 ;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActivePIP3SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double active_PIP2 = rU[PIP2];
    double active_PIP3 = rU[PIP3];
    double active_rac= rU[RAC];
    double active_rho= rU[RHO];

    return PIP2toPIP3BaseConversionCoefficient / 2 * (1 + active_rac / BasalRac) * active_PIP2 - PIP3toPIP2BaseConversionCoefficient / 2 * (1+ active_rho / BasalRho ) * active_PIP3;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ArpDynamicsSourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double active_Arp = rU[ARP];

    return ArpActivationRate(rU, rX, t) - ArpDependentBranching(rU, rX, t) * TotalFilamentDensity(rU, rX, t) - ArpDecayRate*active_Arp;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActinFilamentDensity_M0_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double barbedEndDensity_M0 = rU[B0];
    double actinFilamentDensity_M0 = rU[F0];

    return ActinGrowthRate * barbedEndDensity_M0 - ActinTurnoverRate * actinFilamentDensity_M0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActinFilamentDensity_M1_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double barbedEndDensity_M1 = rU[B1];
    double actinFilamentDensity_M1 = rU[F1];

    return ActinGrowthRate * barbedEndDensity_M1 - ActinTurnoverRate * actinFilamentDensity_M1;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActinFilamentDensity_M2_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double barbedEndDensity_M2 = rU[B2];
    double actinFilamentDensity_M2 = rU[F2];

    return ActinGrowthRate * barbedEndDensity_M2 - ActinTurnoverRate * actinFilamentDensity_M2;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActinFilamentDensity_M3_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double barbedEndDensity_M3 = rU[B3];
    double actinFilamentDensity_M3 = rU[F3];

    return ActinGrowthRate * barbedEndDensity_M3 - ActinTurnoverRate * actinFilamentDensity_M3;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActinFilamentDensity_M4_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double barbedEndDensity_M4 = rU[B4];
    double actinFilamentDensity_M4 = rU[F4];

    return ActinGrowthRate * barbedEndDensity_M4 - ActinTurnoverRate * actinFilamentDensity_M4;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::BarbedEndDensity_M0_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double barbedEndDensity_M0 = rU[B0];
    double actinFilamentDensity_M1 = rU[F1];
    double actinFilamentDensity_M4 = rU[F4];

    return 1.0/2 * ScaleFactorConctoB * ArpDependentBranching(rU, rX, t) * (actinFilamentDensity_M4 + actinFilamentDensity_M1) - CappingBarbedEnd(rU, rX, t) * barbedEndDensity_M0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::BarbedEndDensity_M1_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double barbedEndDensity_M1 = rU[B1];
    double actinFilamentDensity_M2 = rU[F2];
    double actinFilamentDensity_M0 = rU[F0];

    return 1.0/2 * ScaleFactorConctoB * ArpDependentBranching(rU, rX, t) * (actinFilamentDensity_M0 + actinFilamentDensity_M2) - CappingBarbedEnd(rU, rX, t) * barbedEndDensity_M1;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::BarbedEndDensity_M2_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double barbedEndDensity_M2 = rU[B2];
    double actinFilamentDensity_M3 = rU[F3];
    double actinFilamentDensity_M1 = rU[F1];

    return 1.0/2 * ScaleFactorConctoB * ArpDependentBranching(rU, rX, t) * (actinFilamentDensity_M1 + actinFilamentDensity_M3) - CappingBarbedEnd(rU, rX, t) * barbedEndDensity_M2;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::BarbedEndDensity_M3_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double barbedEndDensity_M3 = rU[B3];
    double actinFilamentDensity_M4 = rU[F4];
    double actinFilamentDensity_M2 = rU[F2];

    return 1.0/2 * ScaleFactorConctoB * ArpDependentBranching(rU, rX, t) * (actinFilamentDensity_M2 + actinFilamentDensity_M4) - CappingBarbedEnd(rU, rX, t) * barbedEndDensity_M3;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::BarbedEndDensity_M4_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double barbedEndDensity_M4 = rU[B4];
    double actinFilamentDensity_M3 = rU[F3];
    double actinFilamentDensity_M0 = rU[F0];

    return 1.0/2 * ScaleFactorConctoB * ArpDependentBranching(rU, rX, t) * (actinFilamentDensity_M3 + actinFilamentDensity_M0) - CappingBarbedEnd(rU, rX, t) * barbedEndDensity_M4;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ForceBearingBarbedEnd_M0_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double forceBearingBarbedEnd_M0 = rU[P0];
    double barbedEndDensity_M0 = rU[B0];

    const c_vector<double,SPACE_DIM>& edgeNormal = mEdgeNormalForNodesOnBoundary.at(pNode);
    double velocityDotEdgeNormal = VELOCITY[0][0]*edgeNormal[0] + VELOCITY[0][1]*edgeNormal[1];

    return barbedEndDensity_M0 * velocityDotEdgeNormal - ReductionCappingLeadingEdge * CappingBarbedEnd(rU, rX, t) * forceBearingBarbedEnd_M0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ForceBearingBarbedEnd_M1_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double forceBearingBarbedEnd_M1 = rU[P1];
    double barbedEndDensity_M1 = rU[B1];

    const c_vector<double,SPACE_DIM>& edgeNormal = mEdgeNormalForNodesOnBoundary.at(pNode);
    double velocityDotEdgeNormal = VELOCITY[1][0]*edgeNormal[0] + VELOCITY[1][1]*edgeNormal[1];

    return barbedEndDensity_M1 * velocityDotEdgeNormal - ReductionCappingLeadingEdge * CappingBarbedEnd(rU, rX, t) * forceBearingBarbedEnd_M1;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ForceBearingBarbedEnd_M2_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double forceBearingBarbedEnd_M2 = rU[P2];
    double barbedEndDensity_M2 = rU[B2];

    const c_vector<double,SPACE_DIM>& edgeNormal = mEdgeNormalForNodesOnBoundary.at(pNode);
    double velocityDotEdgeNormal = VELOCITY[2][0]*edgeNormal[0] + VELOCITY[2][1]*edgeNormal[1];

    return barbedEndDensity_M2 * velocityDotEdgeNormal - ReductionCappingLeadingEdge * CappingBarbedEnd(rU, rX, t) * forceBearingBarbedEnd_M2;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ForceBearingBarbedEnd_M3_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double forceBearingBarbedEnd_M3 = rU[P3];
    double barbedEndDensity_M3 = rU[B3];

    const c_vector<double,SPACE_DIM>& edgeNormal = mEdgeNormalForNodesOnBoundary.at(pNode);
    double velocityDotEdgeNormal = VELOCITY[3][0]*edgeNormal[0] + VELOCITY[3][1]*edgeNormal[1];

    return barbedEndDensity_M3 * velocityDotEdgeNormal - ReductionCappingLeadingEdge * CappingBarbedEnd(rU, rX, t) * forceBearingBarbedEnd_M3;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ForceBearingBarbedEnd_M4_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double forceBearingBarbedEnd_M4 = rU[P4];
    double barbedEndDensity_M4 = rU[B4];

    const c_vector<double,SPACE_DIM>& edgeNormal = mEdgeNormalForNodesOnBoundary.at(pNode);
    double velocityDotEdgeNormal = VELOCITY[4][0]*edgeNormal[0] + VELOCITY[4][1]*edgeNormal[1];

    return barbedEndDensity_M4 * velocityDotEdgeNormal - ReductionCappingLeadingEdge * CappingBarbedEnd(rU, rX, t) * forceBearingBarbedEnd_M4;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::CDC42ActivationRate(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const
{
    double active_rho = rU[RHO];
    double active_PIP3 = rU[PIP3];

    double stimulus = mpStimulusProtocol->GetCDC42ActivationRate(pNode, t);

    return ((stimulus +  CDC42BaselineActivationRate) / (1 + pow((active_rho / RhoHalfMaximalDropCDC42), HillCoeffMutualInhibition))) * ( (1- feedback) + feedback * active_PIP3  / BasalPIP3);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::RacActivationRate(const c_vector<double,GTPASE_PROBLEM_DIM>& rU) const
{
    double active_cdc42 = rU[CDC42];
    double active_PIP3 = rU[PIP3];

    return ( RacBaselineActivationRate + CDC42AmplificationOfRac * active_cdc42) * ((1- feedback) + feedback * active_PIP3  / BasalPIP3);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::RhoActivationRate(const c_vector<double,GTPASE_PROBLEM_DIM>& rU) const
{
    double active_cdc42 = rU[CDC42];
    double active_rac = rU[RAC];

    return (RhoBaselineActivationRate + RacEnhancedRhoActivation * active_rac) /
            (1 + pow(active_cdc42 / CDC42HalfMaximalDropCDC42, HillCoeffMutualInhibition));
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ArpActivationRate(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t) const
{
    double active_cdc42 = rU[CDC42];
    double active_PIP2 = rU[PIP2];

    return ArpActivationRatePIP2dep/2 * (pow(active_PIP2, HillCoeffPIP2MediatedArpActivation)/ (pow(ThresholdPIP2forArpActivation,HillCoeffPIP2MediatedArpActivation) + pow(active_PIP2,HillCoeffPIP2MediatedArpActivation))) * (1 + active_cdc42/BasalCDC42);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ArpDependentBranching(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t) const
{
    double active_Arp = rU[ARP];

    return ArpNucleationRate* (active_Arp / (SaturationConstantArpNucleation + active_Arp + ScaleFactorFtoConc * TotalFilamentDensity(rU, rX, t)));
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::CappingBarbedEnd(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t) const
{
   double active_PIP2 = rU[PIP2];

   return BarbedEndCappingRate - MaxReductionCappingPIP2 * (pow(active_PIP2, HillCoeffPIP2MediatedArpActivation)/ (pow(ThresholdPIP2forArpActivation,HillCoeffPIP2MediatedArpActivation) + pow(active_PIP2,HillCoeffPIP2MediatedArpActivation))) ;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::TotalFilamentDensity(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t) const
{
    double total_filaments = 0.0;

    for (unsigned variable=F0; variable<=F4; ++variable)
    {
        total_filaments += rU[variable];
    }

    return total_filaments;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::TotalBarbedEndDensity(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t) const
{
    double total_barbed = 0.0;

    for (unsigned variable=B0; variable<=B4; ++variable)
    {
        total_barbed += rU[variable];
    }

    return total_barbed;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::TotalForceBearingBarbedEnd(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t) const
{
    double total_force_bearing = 0.0;

    for (unsigned variable=P0; variable<=P4; ++variable)
    {
        total_force_bearing += rU[variable];
    }

    return total_force_bearing;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::SetupLinearSystem(Vec currentSolution, bool computeMatrix)
{
    if (computeMatrix)
    {
        GTPasePDESystemLHSAssembler<ELEMENT_DIM, SPACE_DIM> lhs_matrix_assembler(this->mpMesh);
        GTPasePDESystemMassMatrixAssembler<ELEMENT_DIM, SPACE_DIM> mass_matrix_assembler(this->mpMesh);

        lhs_matrix_assembler.SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
        lhs_matrix_assembler.AssembleMatrix();

        mass_matrix_assembler.SetMatrixToAssemble(mMassMatrix);
        mass_matrix_assembler.AssembleMatrix();

        this->mpLinearSystem->FinaliseLhsMatrix();
        PetscMatTools::Finalise(mMassMatrix);
    }

    //////////////////////////////////////////
    // Set up z in b=Mz
    //////////////////////////////////////////
    double dt = PdeSimulationTime::GetPdeTimeStep();
    double t = PdeSimulationTime::GetTime();

    // Create distributed vectors from the current solution and the vector used for RHS assemble to facilitate parallel access
    DistributedVectorFactory* p_vector_factory = this->mpMesh->GetDistributedVectorFactory();
    DistributedVector distributed_current_solution = p_vector_factory->CreateDistributedVector(currentSolution, true);
    DistributedVector distributed_vec_for_rhs = p_vector_factory->CreateDistributedVector(mVecForConstructingRhs);

    // The stripe objects will allow access to each variable independently using an iterator
    std::vector<DistributedVector::Stripe*> current_solution_variables(GTPASE_PROBLEM_DIM);
    std::vector<DistributedVector::Stripe*> vec_for_rhs_variables(GTPASE_PROBLEM_DIM);
    for (unsigned variable = 0; variable < GTPASE_PROBLEM_DIM; ++variable)
    {
        current_solution_variables[variable] = new DistributedVector::Stripe(distributed_current_solution, variable);
        vec_for_rhs_variables[variable] = new DistributedVector::Stripe(distributed_vec_for_rhs, variable);
    }

    // Iterate over the vector used for RHS assemble and populate it by calling the PDE source terms
    // Each iteration corresponds to one mesh node with the GTPASE_PROBLEM_DIM dof accessible via the stride objects
    for (DistributedVector::Iterator index = distributed_vec_for_rhs.Begin();
         index!= distributed_current_solution.End();
         ++index)
    {
        Node<SPACE_DIM>* p_node = this->mpMesh->GetNode(index.Global);

        // Put together an auxiliary vector with all the dof corresponding to the current mesh node
        c_vector<double,GTPASE_PROBLEM_DIM> node_dofs;
        for (unsigned variable = 0; variable < GTPASE_PROBLEM_DIM; ++variable)
        {
            node_dofs[variable] = (*current_solution_variables[variable])[index];
        }

        // Calculate the vector entry for each of the dof associated with the current node
        for (unsigned variable = 0; variable < GTPASE_PROBLEM_DIM; ++variable)
        {
            if (EQUATION_TYPE[variable] != ODE)
            {
                double current_solution = node_dofs[variable];
                double rhs_value = current_solution / dt + (this->*mSourceTerms[variable])(node_dofs, p_node->GetPoint(), t+dt, p_node);
                (*vec_for_rhs_variables[variable])[index] = rhs_value;
            }
        }

    }

    // Flush element accessed via stride objects back into the original Vec
    distributed_vec_for_rhs.Restore();

    // Clean up stride objects
    for (unsigned variable = 0; variable < GTPASE_PROBLEM_DIM; ++variable)
    {
        delete current_solution_variables[variable];
        delete vec_for_rhs_variables[variable];
    }

    //////////////////////////////////////////
    // b = Mz
    //////////////////////////////////////////
    MatMult(mMassMatrix, mVecForConstructingRhs, this->mpLinearSystem->rGetRhsVector());

    /* Apply the dirichlet BCs from the BCC to the linear system */
    mBoundaryConditions.ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);

//        NaturalNeumannSurfaceTermAssembler<DIM,DIM,1> surface_integral_assembler(this->mpMesh, mpBoundaryConditions);
//        surface_integral_assembler.SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false /*don't zero vector before assembling!*/);
//        surface_integral_assembler.Assemble();
//
//        /* Some necessary PETSc communication before applying Dirichet BCs */
//        this->mpLinearSystem->FinaliseRhsVector();         // (Petsc communication)
//        this->mpLinearSystem->SwitchWriteModeLhsMatrix();  // (Petsc communication - needs to called when going from adding entries to inserting entries)

    FinaliseODEEquations(currentSolution, computeMatrix);

    // Reassembling needed after setting individual elements
    this->mpLinearSystem->AssembleFinalLinearSystem();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::FinaliseODEEquations(Vec currentSolution, bool computeMatrix)
{
    // When assembling the matrix, set to 1 the diagonal entries corresponding to the ODE equations
    if (computeMatrix)
    {
        unsigned num_equations = this->mpLinearSystem->GetSize();
        for (unsigned equation_num=0; equation_num<num_equations; equation_num++)
        {
            if (EQUATION_TYPE[equation_num % GTPASE_PROBLEM_DIM] == ODE)
            {
                PetscMatTools::SetElement(this->mpLinearSystem->rGetLhsMatrix(), equation_num, equation_num, 1.0);
            }
        }
    }

    double dt = PdeSimulationTime::GetPdeTimeStep();
    double t = PdeSimulationTime::GetTime();

    // Create distributed vectors from the current solution and the system RHS to facilitate parallel access
    DistributedVectorFactory* p_vector_factory = this->mpMesh->GetDistributedVectorFactory();
    DistributedVector distributed_current_solution = p_vector_factory->CreateDistributedVector(currentSolution, true);
    DistributedVector distributed_rhs = p_vector_factory->CreateDistributedVector(this->mpLinearSystem->rGetRhsVector());

    // The stripe objects will allow access to each variable independently using an iterator
    std::vector<DistributedVector::Stripe*> current_solution_variables(GTPASE_PROBLEM_DIM);
    std::vector<DistributedVector::Stripe*> rhs_variables(GTPASE_PROBLEM_DIM);
    for (unsigned variable = 0; variable < GTPASE_PROBLEM_DIM; ++variable)
    {
        current_solution_variables[variable] = new DistributedVector::Stripe(distributed_current_solution, variable);
        rhs_variables[variable] = new DistributedVector::Stripe(distributed_rhs, variable);
    }

    // Iterate over the system RHS and update the entries corresponding to ODEs
    // Each iteration corresponds to one mesh node with the GTPASE_PROBLEM_DIM dof accessible via the stride objects
    for (DistributedVector::Iterator index = distributed_rhs.Begin();
         index!= distributed_rhs.End();
         ++index)
    {
        Node<SPACE_DIM>* p_node = this->mpMesh->GetNode(index.Global);

        // Put together an auxiliary vector with all the dof corresponding to the current mesh node
        c_vector<double,GTPASE_PROBLEM_DIM> node_dofs;
        for (unsigned variable = 0; variable < GTPASE_PROBLEM_DIM; ++variable)
        {
            node_dofs[variable] = (*current_solution_variables[variable])[index];
        }

        // Calculate the vector entries corresponding to ODEs
        for (unsigned f_variable = F0; f_variable <= F4; ++f_variable)
        {
            double current_solution = node_dofs[f_variable];
            double rhs_value = current_solution + dt * (this->*mSourceTerms[f_variable])(node_dofs, p_node->GetPoint(), t+dt, p_node);
            (*rhs_variables[f_variable])[index] = rhs_value;
        }
        for (unsigned p_variable = P0; p_variable <= P4; ++p_variable)
        {
            double current_solution = node_dofs[p_variable];
            // Set to zero the RHS for pushing barbed ends equations defined in internal nodes
            double rhs_value = (!p_node->IsBoundaryNode()) ? 0.0 : current_solution + dt * (this->*mSourceTerms[p_variable])(node_dofs, p_node->GetPoint(), t+dt, p_node);
            (*rhs_variables[p_variable])[index] = rhs_value;
        }
    }

    distributed_rhs.Restore();

    for (unsigned variable = 0; variable < GTPASE_PROBLEM_DIM; ++variable)
    {
        delete current_solution_variables[variable];
        delete rhs_variables[variable];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    // First create the linear system
    AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,GTPASE_PROBLEM_DIM>::InitialiseForSolve(initialSolution);

    Vec& r_template = this->mpLinearSystem->rGetRhsVector();
    VecDuplicate(r_template, &mVecForConstructingRhs);
    PetscInt ownership_range_lo;
    PetscInt ownership_range_hi;
    VecGetOwnershipRange(r_template, &ownership_range_lo, &ownership_range_hi);
    PetscInt local_size = ownership_range_hi - ownership_range_lo;
    unsigned dofs = GTPASE_PROBLEM_DIM*this->mpMesh->GetNumNodes();
    PetscTools::SetupMat(mMassMatrix, dofs, dofs, GTPASE_PROBLEM_DIM*this->mpMesh->CalculateMaximumNodeConnectivityPerProcess(), local_size, local_size);

    // System is not symmetric (due to discretisation of convective operator in hyperbolic equations), use GMRES.
    this->mpLinearSystem->SetKspType("gmres");

    /// \todo Use PETSc's default preconditioner for now. AMG seemed to slow down simulations. This needs further investigation.
//    this->mpLinearSystem->SetPcType("hypre");
//    PetscTools::SetOption("-pc_hypre_type", "boomeramg");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::PopulateBoundaryConditionsContainer(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    assert(SPACE_DIM == 2);

    // Specify no-flux boundary conditions (i.e. zero Neumann) for reaction-diffusion equations (parabolic) and zero concentration
    // of B at the inlets of the reaction-convection equations (hyperbolic)

    // Object gets freed inside BoundaryConditionsContainer
    ConstBoundaryCondition<2>* p_dirichlet_bc = new ConstBoundaryCondition<2>(0.0);

    // TODO: Velocity is now a static member variable that we should be able to access from here without recalculating it.
    c_matrix<double, 5, 2> Velocity;
    for (unsigned orientation_index=0; orientation_index<5; orientation_index++)
    {
        Velocity(orientation_index,0) = ActinGrowthRate * cos(ANGLES[orientation_index]);
        Velocity(orientation_index,1) = ActinGrowthRate * sin(ANGLES[orientation_index]);
    }

    for (unsigned unknown_index=0; unknown_index<GTPASE_PROBLEM_DIM; unknown_index++)
    {
        switch(EQUATION_TYPE[unknown_index])
        {
            case PARABOLIC_PDE:
                mBoundaryConditions.DefineZeroNeumannOnMeshBoundary(pMesh, unknown_index);
                break;
            case HYPERBOLIC_PDE:
            {
                for (typename TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator elem_iter = pMesh->GetBoundaryElementIteratorBegin();
                     elem_iter != pMesh->GetBoundaryElementIteratorEnd();
                     ++elem_iter)
                {
                    c_vector<double, SPACE_DIM> normal = (*elem_iter)->CalculateNormal();
                    c_vector<double, SPACE_DIM> velocity = matrix_row<c_matrix<double, 5, SPACE_DIM> >(Velocity, unknown_index-B0);

                    // dot product between velocity and edge normal is negative
                    if (inner_prod(normal, velocity) < 0)
                    {
                        mBoundaryConditions.AddDirichletBoundaryCondition((*elem_iter)->GetNode(0), p_dirichlet_bc, unknown_index);
                        mBoundaryConditions.AddDirichletBoundaryCondition((*elem_iter)->GetNode(1), p_dirichlet_bc, unknown_index);
                    }
                }

                break;
            }
            case ODE:
                break;
            default:
                NEVER_REACHED;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::PopulateEdgeNormalForNodesOnBoundary(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    typedef std::map<unsigned, std::vector<c_vector<double,ELEMENT_DIM> > > NormalsList;
    NormalsList all_normals_found_for_node_on_boundary;

    for (TetrahedralMesh<2,2>::BoundaryElementIterator it = pMesh->GetBoundaryElementIteratorBegin();
         it != pMesh->GetBoundaryElementIteratorEnd();
         ++it)
    {
        assert((*it)->GetNumNodes() == 2);
        for (unsigned i=0; i<(*it)->GetNumNodes(); i++)
        {
            all_normals_found_for_node_on_boundary[(*it)->GetNode(i)->GetIndex()].push_back((*it)->CalculateNormal());
        }
    }

    for (typename NormalsList::const_iterator iter = all_normals_found_for_node_on_boundary.begin();
         iter != all_normals_found_for_node_on_boundary.end(); ++iter)
    {
        c_vector<double,ELEMENT_DIM> mean_vector = std::accumulate(iter->second.begin(), iter->second.end(), Create_c_vector(0,0));
        if(norm_2(mean_vector) == 0.0)
        {
            // This happens in some degenerate (yet still correct) cases, specially on regular grids.
            std::cout << "WARNING: zero normal computed for node " << iter->first << std::endl;
        }
        else
        {
            mean_vector /= norm_2(mean_vector);
        }
        mEdgeNormalForNodesOnBoundary[pMesh->GetNode(iter->first)] = mean_vector;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::GTPasePDESystemSolver(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh, const AbstractCDC42StimulusProtocol<ELEMENT_DIM,SPACE_DIM>* stimulusProtocol, Vec initialCondition)
    : AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,GTPASE_PROBLEM_DIM>(pMesh),
      mFreeInitialCondition(false),
      mpStimulusProtocol(stimulusProtocol)
{
    // At the moment VELOCITY is hardcoded as SPACE_DIM=2
    assert(SPACE_DIM==2);

    PopulateBoundaryConditionsContainer(pMesh);

    /*
     * By default the dynamic solvers will reassemble the matrix each timestep. In this problem the matrix is constant
     * and only needs to be assembled once. Make sure we tell the solver this, otherwise performance will be destroyed.
     */
    this->mMatrixIsConstant = true;

    // Register the callback functions used to evaluate each source term
    mSourceTerms[0] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActiveCDC42SourceTerm;
    mSourceTerms[1] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActiveRacSourceTerm;
    mSourceTerms[2] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActiveRhoSourceTerm;
    mSourceTerms[3] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::InactiveCDC42SourceTerm;
    mSourceTerms[4] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::InactiveRacSourceTerm;
    mSourceTerms[5] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::InactiveRhoSourceTerm;
    mSourceTerms[6] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActivePIP1SourceTerm;
    mSourceTerms[7] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActivePIP2SourceTerm;
    mSourceTerms[8] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActivePIP3SourceTerm;
    mSourceTerms[9] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ArpDynamicsSourceTerm;
    mSourceTerms[10] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActinFilamentDensity_M0_SourceTerm;
    mSourceTerms[11] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActinFilamentDensity_M1_SourceTerm;
    mSourceTerms[12] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActinFilamentDensity_M2_SourceTerm;
    mSourceTerms[13] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActinFilamentDensity_M3_SourceTerm;
    mSourceTerms[14] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ActinFilamentDensity_M4_SourceTerm;
    mSourceTerms[15] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::BarbedEndDensity_M0_SourceTerm;
    mSourceTerms[16] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::BarbedEndDensity_M1_SourceTerm;
    mSourceTerms[17] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::BarbedEndDensity_M2_SourceTerm;
    mSourceTerms[18] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::BarbedEndDensity_M3_SourceTerm;
    mSourceTerms[19] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::BarbedEndDensity_M4_SourceTerm;
    mSourceTerms[20] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ForceBearingBarbedEnd_M0_SourceTerm;
    mSourceTerms[21] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ForceBearingBarbedEnd_M1_SourceTerm;
    mSourceTerms[22] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ForceBearingBarbedEnd_M2_SourceTerm;
    mSourceTerms[23] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ForceBearingBarbedEnd_M3_SourceTerm;
    mSourceTerms[24] = &GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ForceBearingBarbedEnd_M4_SourceTerm;

    if (initialCondition == NULL)
    {
        this->SetInitialCondition(CreateInitialCondition(pMesh));
        mFreeInitialCondition = true;
    }
    else
    {
        this->SetInitialCondition(initialCondition);
    }

    pMesh->CheckOutwardNormals();
    PopulateEdgeNormalForNodesOnBoundary(pMesh);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Vec GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::CreateInitialCondition(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    std::vector<double> init_cond(GTPASE_PROBLEM_DIM*pMesh->GetNumNodes(), 0.0);
    for (unsigned node_index=0; node_index<pMesh->GetNumNodes(); ++node_index)
    {
        init_cond[GTPASE_PROBLEM_DIM*node_index] = BasalCDC42;
        init_cond[GTPASE_PROBLEM_DIM*node_index + 1] = BasalRac;
        init_cond[GTPASE_PROBLEM_DIM*node_index + 2] = BasalRho;
        init_cond[GTPASE_PROBLEM_DIM*node_index + 3] = TotalCDC42 - BasalCDC42;
        init_cond[GTPASE_PROBLEM_DIM*node_index + 4] = TotalRac - BasalRac;
        init_cond[GTPASE_PROBLEM_DIM*node_index + 5] = TotalRho - BasalRho;
        init_cond[GTPASE_PROBLEM_DIM*node_index + 6] = BasalPIP1;
        init_cond[GTPASE_PROBLEM_DIM*node_index + 7] = BasalPIP2;
        init_cond[GTPASE_PROBLEM_DIM*node_index + 8] = BasalPIP3;
        init_cond[GTPASE_PROBLEM_DIM*node_index + 9] = BasalARP;

        for (unsigned angle_index=F0; angle_index<=F4; angle_index++)
        {
            init_cond[GTPASE_PROBLEM_DIM*node_index + angle_index] = BasalFilamentDensity;
        }

        for (unsigned angle_index=B0; angle_index<=B4; angle_index++)
        {
            init_cond[GTPASE_PROBLEM_DIM*node_index + angle_index] = BasalBarbedEndDensity;
        }

        for (unsigned angle_index=P0; angle_index<=P4; angle_index++)
        {
            init_cond[GTPASE_PROBLEM_DIM*node_index + angle_index] = BasalEdgeBarbedEndDensity;
        }
    }
    return PetscTools::CreateVec(init_cond);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::~GTPasePDESystemSolver()
{
    if (mFreeInitialCondition)
    {
        PetscTools::Destroy(this->mInitialCondition);
    }

    PetscTools::Destroy(mMassMatrix);
    PetscTools::Destroy(mVecForConstructingRhs);
}

// Explicit instantiation
template class GTPasePDESystemSolver<2,2>;

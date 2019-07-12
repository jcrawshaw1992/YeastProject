#ifndef GTPASEPDESYSTEMSOLVER_HPP_
#define GTPASEPDESYSTEMSOLVER_HPP_

#include <numeric>

// We need to include the following two headers if we are going to use a combination
// of (element_dim, space_dim, problem_dim) that isn't explicitly instantiated
// in BoundaryConditionsContainer.cpp
#include "BoundaryConditionsContainerImplementation.hpp"
#include "AbstractBoundaryConditionsContainerImplementation.hpp"

#include "AbstractDynamicLinearPdeSolver.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractCDC42StimulusProtocol.hpp"
#include "GTPasePDESystemLHSAssembler.hpp"
#include "GTPasePDESystemMassMatrixAssembler.hpp"
#include "GTPasePDESystemParameters.hpp"

/*
 * This class implements a solver for a system of 6 reaction-diffusion equations modelling GTPase signalling
 * involved in cell polarisation and motility as described in Maree et al, PLOS Comp. Bio. (2012).
 *
 * Let V=[u_0, u_1, ..., u_5] be the vector of unknowns, the system takes the form:
 *
 *    V_t = a*Laplacian(V) + f(V,x,t)
 *
 * where a is the diffusion coefficient and f the source term. No-flux boundary conditions dV/dn = 0 are applied
 * to the whole domain boundary. The code implements a semi-implicit time discretisation such that
 *
 *    (u_i^{n+1} - u_i^{n})/dt = a*Laplacian(u_i^{n+1}) + f(V^{n},x,t^{n+1}).
 *
 * This leads to a block diagonal linear system with each of the 6 blocks being of the form [M/dt+a*K] U_i^{n+1}
 * with the associated right-hand-side being M/dt U^{n} + b, with b_i = integral( f(V^{n},x_i,t^{n+1}) phi_i dV ).
 *
 * The unknown are ordered: CDC42, Rac, Rho, Inactive CDC42, Inactive Rac, Inactive Rho (indexed from 0 to GTPASE_PROBLEM_DIM-1)
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class GTPasePDESystemSolver
    : public AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM, GTPASE_PROBLEM_DIM>
{
private:

    static const double VELOCITY[][SPACE_DIM];

    double ActiveCDC42SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ActiveRacSourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ActiveRhoSourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double InactiveCDC42SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double InactiveRacSourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double InactiveRhoSourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ActivePIP1SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ActivePIP2SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ActivePIP3SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ArpDynamicsSourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ActinFilamentDensity_M0_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ActinFilamentDensity_M1_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ActinFilamentDensity_M2_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ActinFilamentDensity_M3_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ActinFilamentDensity_M4_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double BarbedEndDensity_M0_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double BarbedEndDensity_M1_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double BarbedEndDensity_M2_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double BarbedEndDensity_M3_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double BarbedEndDensity_M4_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ForceBearingBarbedEnd_M0_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ForceBearingBarbedEnd_M1_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ForceBearingBarbedEnd_M2_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ForceBearingBarbedEnd_M3_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double ForceBearingBarbedEnd_M4_SourceTerm(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;

    double CDC42ActivationRate(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;
    double RacActivationRate(const c_vector<double,GTPASE_PROBLEM_DIM>& rU) const;
    double RhoActivationRate(const c_vector<double,GTPASE_PROBLEM_DIM>& rU) const;
    double ArpActivationRate(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t) const;
    double ArpDependentBranching(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t) const;
    double CappingBarbedEnd(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t) const;

    double TotalFilamentDensity(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t) const;
    double TotalBarbedEndDensity(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t) const;
    double TotalForceBearingBarbedEnd(const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t) const;

    void SetupLinearSystem(Vec currentSolution, bool computeMatrix);
    void FinaliseODEEquations(Vec currentSolution, bool computeMatrix);
    void InitialiseForSolve(Vec initialSolution);
    void PopulateBoundaryConditionsContainer(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);
    void PopulateEdgeNormalForNodesOnBoundary(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    // Mass matrix, used to computing the RHS vector
    Mat mMassMatrix;


    // The vector z to be multiplied by the mass matrix M to compute the system RHS, i.e. b=Mz.
    Vec mVecForConstructingRhs;

    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,GTPASE_PROBLEM_DIM> mBoundaryConditions;

    std::map<const Node<SPACE_DIM>*, c_vector<double,ELEMENT_DIM> > mEdgeNormalForNodesOnBoundary;

    // Array that stores the pointers to the functions implementing the source terms of each PDE
    double (GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::* mSourceTerms[GTPASE_PROBLEM_DIM]) (const c_vector<double,GTPASE_PROBLEM_DIM>& rU, const ChastePoint<SPACE_DIM>& rX, double t, const Node<SPACE_DIM>* pNode) const;

    // If the initial condition was allocated within this object, it needs freeing in the destructor
    bool mFreeInitialCondition;

    const AbstractCDC42StimulusProtocol<ELEMENT_DIM,SPACE_DIM>* mpStimulusProtocol;

public:
    GTPasePDESystemSolver(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh, const AbstractCDC42StimulusProtocol<ELEMENT_DIM,SPACE_DIM>* stimulusProtocol, Vec initialCondition=NULL);

    static Vec CreateInitialCondition(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    ~GTPasePDESystemSolver();
};

#endif /* GTPASEPDESYSTEMSOLVER_HPP_ */

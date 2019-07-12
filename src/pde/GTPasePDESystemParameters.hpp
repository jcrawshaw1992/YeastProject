#ifndef GTPASEPDESYSTEMPARAMETERS_HPP_
#define GTPASEPDESYSTEMPARAMETERS_HPP_

#define GTPASE_PROBLEM_DIM 25

// Unknown names, this enumerated type should contain GTPASE_PROBLEM_DIM entries
enum UnknownNames
{
    CDC42 = 0,
    RAC,
    RHO,
    INACTIVE_CDC42,
    INACTIVE_RAC,
    INACTIVE_RHO,
    PIP1,
    PIP2,
    PIP3,
    ARP,
    F0,
    F1,
    F2,
    F3,
    F4,
    B0,
    B1,
    B2,
    B3,
    B4,
    P0,
    P1,
    P2,
    P3,
    P4
};

// We need to handle three different types of differential equations, with different discretisations
enum EquationType
{
    PARABOLIC_PDE,
    ODE,
    HYPERBOLIC_PDE
};

// Label each equation with its type
const EquationType EQUATION_TYPE[GTPASE_PROBLEM_DIM] = {PARABOLIC_PDE, PARABOLIC_PDE, PARABOLIC_PDE, PARABOLIC_PDE, PARABOLIC_PDE,
                                                        PARABOLIC_PDE, PARABOLIC_PDE, PARABOLIC_PDE, PARABOLIC_PDE, PARABOLIC_PDE,
                                                        ODE, ODE, ODE, ODE, ODE,
                                                        HYPERBOLIC_PDE, HYPERBOLIC_PDE, HYPERBOLIC_PDE, HYPERBOLIC_PDE, HYPERBOLIC_PDE,
                                                        ODE, ODE, ODE, ODE, ODE};

// Orientation angle for each family of actin filaments
const double ANGLES[] = {0., 2*M_PI/5, 4*M_PI/5, 6*M_PI/5, 8*M_PI/5};

const double CDC42BaselineActivationRate = 2.95 ; // I_C
const double RacBaselineActivationRate = 0.5; // I_r
const double RhoBaselineActivationRate = 3.3; // I_rho
const unsigned HillCoeffMutualInhibition = 3; // n
const double RhoHalfMaximalDropCDC42 = 1.25; // a_1
const double CDC42HalfMaximalDropCDC42 = 1.0; // a_2
const double CDC42AmplificationOfRac = 4.5; // alpha
const double RacEnhancedRhoActivation = 0.3; // beta
const double TotalCDC42 = 2.4; // C_tot
const double TotalRac = 7.5; // R_tot
const double TotalRho = 3.1; // rho_tot
const double BasalCDC42 = 1.0; // C_b
const double BasalRac = 3.0; // R_b
const double BasalRho = 1.25; // rho_b
const double CDC42DecayRate = 1.0; // d_C
const double RacDecayRate = 1.0; // d_R
const double RhoDecayRate = 1.0; // d_rho
const double ActiveDiffusionCoefficient = 0.1e-12; // D_m
const double InactiveDiffusionCoefficient = 50e-12; // D_mc
const double PIDiffusionCoefficient = 5e-12; // D_P
const double PIP1ImputCoefficient = 0.0;//10.5; // I_p1
const double PIP1DecayCoefficient = 0.0;//0.21; // d_p1
const double PIP1toPIP2BaseConversionCoefficient = 0.084; // k_PI5K
const double PIP2toPIP1BaseConversionCoefficient = 0.14; // 0.021 is the value given in page 14 of Maree 2012. The value in Table 1 is 0.014 // k_K21
const double PIP2toPIP3BaseConversionCoefficient = 0.00072; // k_PI3K
const double PIP3toPIP2BaseConversionCoefficient = 0.432; // k_PTEN
const double BasalPIP1 = 50.0; // PIP1_b
const double BasalPIP2 = 30.0; // PIP2_b
const double BasalPIP3 = 0.05; // PIP3_b
const double feedback = 0.4; // 0<=feedback<=1
const double ArpActivationRatePIP2dep = 1.0; // u_arp_PIP2dep
const unsigned HillCoeffPIP2MediatedArpActivation = 3; // n_arp
const double ThresholdPIP2forArpActivation = 35.0; // P_2half
const double ArpDecayRate = 0.1; // d_Arp
const double ArpDiffusionCoefficient = 1e-12; // D_Arp
const double ArpNucleationRate = 60e-9; // n_arp
const double SaturationConstantArpNucleation = 2.0; // K_m
const double ScaleFactorFtoConc = 0.255e-6; // l
const double ScaleFactorConctoB = 106e+12; // k
const double ActinGrowthRate = 0.5e-6; // v_o
const double ActinTurnoverRate = 0.03; // d_F
const double BarbedEndCappingRate = 2.8; // k_max
const double MaxReductionCappingPIP2 = 2.8; // k_P2
const double ReductionCappingLeadingEdge = 0.14; // r
const double CouplingEnergyperboundarysite = 0.75e6; // J_cm
const double CellInelasticity = 4e18; // lambda_a
const double TargetAreaofCell = 300e-12; // A
const double MembraneInelasticity = 0.016e18; // lambda_p
const double TargetPerimeterofCell = 150e-12; // P
const double MembraneYield = 0.046e9; // H_b
const double Temperature = 0.008e9; // T
const double RhoEffectonContraction = 0.0025e9; // E  *Correction Maree
const double RhoContractionThreshold = 1.25; // P_th  *Correction Maree
const double BasalARP = 2; // A*
const double BasalFilamentDensity = 0.278e9; // F*
const double BasalBarbedEndDensity = 1.7e13; // B*
const double BasalEdgeBarbedEndDensity = 0.05e9; // P*

#endif

#ifndef GTPASEPDESYSTEMLHSASSEMBLER_HPP_
#define GTPASEPDESYSTEMLHSASSEMBLER_HPP_

#include "AbstractFeVolumeIntegralAssembler.hpp"
#include "PdeSimulationTime.hpp"
#include "GTPasePDESystemParameters.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class GTPasePDESystemLHSAssembler
    : public AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, GTPASE_PROBLEM_DIM, false, true, NORMAL>
{
private:

    static const double VELOCITY[][SPACE_DIM];

    /* Provide the (elemental contribution to the) LHS matrix. */
    c_matrix<double,GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1),GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeMatrixTerm(c_vector<double,ELEMENT_DIM+1>& rPhi,
                                                                                                             c_matrix<double,SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
                                                                                                             ChastePoint<SPACE_DIM>& rX,
                                                                                                             c_vector<double,GTPASE_PROBLEM_DIM>& rU,
                                                                                                             c_matrix<double,GTPASE_PROBLEM_DIM,SPACE_DIM>& rGradU,
                                                                                                             Element<ELEMENT_DIM,SPACE_DIM>* pElement);

    // Diffusion coefficient in each PDE.
    double mDiffusionCoefficients[GTPASE_PROBLEM_DIM];

public:

    GTPasePDESystemLHSAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);
};

#endif

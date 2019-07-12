#ifndef GTPASEPDESYSTEMMASSMATRIXASSEMBLER_HPP_
#define GTPASEPDESYSTEMMASSMATRIXASSEMBLER_HPP_

#include "AbstractFeVolumeIntegralAssembler.hpp"
#include "GTPasePDESystemParameters.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class GTPasePDESystemMassMatrixAssembler
    : public AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, GTPASE_PROBLEM_DIM, false, true, NORMAL>
{
private:

    c_matrix<double,GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1),GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeMatrixTerm(c_vector<double,ELEMENT_DIM+1>& rPhi,
                                                                                                             c_matrix<double,SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
                                                                                                             ChastePoint<SPACE_DIM>& rX,
                                                                                                             c_vector<double,GTPASE_PROBLEM_DIM>& rU,
                                                                                                             c_matrix<double,GTPASE_PROBLEM_DIM,SPACE_DIM>& rGradU,
                                                                                                             Element<ELEMENT_DIM,SPACE_DIM>* pElement);

public:

    GTPasePDESystemMassMatrixAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);
};

#endif /* GTPASEPDESYSTEMMASSMATRIXASSEMBLER_HPP_ */

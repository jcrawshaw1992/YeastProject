#include "GTPasePDESystemMassMatrixAssembler.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1),GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1)>
    GTPasePDESystemMassMatrixAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeMatrixTerm(c_vector<double,ELEMENT_DIM+1>& rPhi,
                                                                                  c_matrix<double,SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
                                                                                  ChastePoint<SPACE_DIM>& rX,
                                                                                  c_vector<double,GTPASE_PROBLEM_DIM>& rU,
                                                                                  c_matrix<double,GTPASE_PROBLEM_DIM,SPACE_DIM>& rGradU,
                                                                                  Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    c_matrix<double,GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1),GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1)> ret = zero_matrix<double>(GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1),GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1));

    for (unsigned i=0; i<ELEMENT_DIM+1; i++)
    {
        for (unsigned j=0; j<ELEMENT_DIM+1; j++)
        {
            for (unsigned k=0; k<GTPASE_PROBLEM_DIM; k++)
            {
                if (EQUATION_TYPE[k] != ODE)
                {
                    // mass matrix
                    ret(GTPASE_PROBLEM_DIM*i+k, GTPASE_PROBLEM_DIM*j+k) = rPhi(i)*rPhi(j);
            	}
            }
        }
    }

    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GTPasePDESystemMassMatrixAssembler<ELEMENT_DIM, SPACE_DIM>::GTPasePDESystemMassMatrixAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    : AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, GTPASE_PROBLEM_DIM, false, true, NORMAL>(pMesh)
{
}

// Explicit instantiation
template class GTPasePDESystemMassMatrixAssembler<2,2>;

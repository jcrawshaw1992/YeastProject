#include "GTPasePDESystemLHSAssembler.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const double GTPasePDESystemLHSAssembler<ELEMENT_DIM, SPACE_DIM>::VELOCITY[][SPACE_DIM] = {{ActinGrowthRate * cos(ANGLES[0]), ActinGrowthRate * sin(ANGLES[0])},
                                                                                           {ActinGrowthRate * cos(ANGLES[1]), ActinGrowthRate * sin(ANGLES[1])},
                                                                                           {ActinGrowthRate * cos(ANGLES[2]), ActinGrowthRate * sin(ANGLES[2])},
                                                                                           {ActinGrowthRate * cos(ANGLES[3]), ActinGrowthRate * sin(ANGLES[3])},
                                                                                           {ActinGrowthRate * cos(ANGLES[4]), ActinGrowthRate * sin(ANGLES[4])}};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1),GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1)>
    GTPasePDESystemLHSAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeMatrixTerm(c_vector<double,ELEMENT_DIM+1>& rPhi,
                                                                           c_matrix<double,SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
                                                                           ChastePoint<SPACE_DIM>& rX,
                                                                           c_vector<double,GTPASE_PROBLEM_DIM>& rU,
                                                                           c_matrix<double,GTPASE_PROBLEM_DIM,SPACE_DIM>& rGradU,
                                                                           Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    c_matrix<double,GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1),GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1)> ret = zero_matrix<double>(GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1),GTPASE_PROBLEM_DIM*(ELEMENT_DIM+1));

    double dt = PdeSimulationTime::GetPdeTimeStep();

    for (unsigned i=0; i<ELEMENT_DIM+1; i++)
    {
        for (unsigned j=0; j<ELEMENT_DIM+1; j++)
        {
            for (unsigned k=0; k<GTPASE_PROBLEM_DIM; k++)
            {
                if (EQUATION_TYPE[k] != ODE)
                {
                    // mass matrix
                    ret(GTPASE_PROBLEM_DIM*i+k, GTPASE_PROBLEM_DIM*j+k) = rPhi(i)*rPhi(j)/dt;
                }
            }

            // stiffness matrix
            for (unsigned dim=0; dim<SPACE_DIM; dim++)
            {
                for (unsigned k=0; k<GTPASE_PROBLEM_DIM; k++)
                {
                    switch(EQUATION_TYPE[k])
                    {
                    case PARABOLIC_PDE:
                        ret(GTPASE_PROBLEM_DIM*i+k, GTPASE_PROBLEM_DIM*j+k) += mDiffusionCoefficients[k]*rGradPhi(dim,i)*rGradPhi(dim,j);
                        break;
                    case HYPERBOLIC_PDE:
                    {
                        int orientation_index = k-B0;
                        assert(orientation_index >= 0);
                        assert(orientation_index < 5);
                        double velocity = VELOCITY[orientation_index][dim];
                        ret(GTPASE_PROBLEM_DIM*i+k, GTPASE_PROBLEM_DIM*j+k) += velocity*rPhi(i)*rGradPhi(dim,j);
                        break;
                    }
                    default:
                        break;
                    }
                }
            }
        }
    }

    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GTPasePDESystemLHSAssembler<ELEMENT_DIM, SPACE_DIM>::GTPasePDESystemLHSAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    : AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, GTPASE_PROBLEM_DIM, false, true, NORMAL>(pMesh)
{
    // Initialise the array containing the diffusion coefficient of each PDE
    mDiffusionCoefficients[0] = ActiveDiffusionCoefficient;
    mDiffusionCoefficients[1] = ActiveDiffusionCoefficient;
    mDiffusionCoefficients[2] = ActiveDiffusionCoefficient;
    mDiffusionCoefficients[3] = InactiveDiffusionCoefficient;
    mDiffusionCoefficients[4] = InactiveDiffusionCoefficient;
    mDiffusionCoefficients[5] = InactiveDiffusionCoefficient;
    mDiffusionCoefficients[6] = PIDiffusionCoefficient;
    mDiffusionCoefficients[7] = PIDiffusionCoefficient;
    mDiffusionCoefficients[8] = PIDiffusionCoefficient;
    mDiffusionCoefficients[9] = ArpDiffusionCoefficient;
}

// Explicit instantiation
template class GTPasePDESystemLHSAssembler<2,2>;

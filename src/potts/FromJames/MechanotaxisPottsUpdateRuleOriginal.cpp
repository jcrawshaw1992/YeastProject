/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "MechanotaxisPottsUpdateRule.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"

template<unsigned DIM>
MechanotaxisPottsUpdateRule<DIM>::MechanotaxisPottsUpdateRule()
    : AbstractPottsUpdateRule<DIM>(),
      mTractionCorrelationParameter(0.01) // @todo Made up defaults

{
}

template<unsigned DIM>
MechanotaxisPottsUpdateRule<DIM>::~MechanotaxisPottsUpdateRule()
{
}

template<unsigned DIM>
double MechanotaxisPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                         unsigned targetNodeIndex,
                                                                         PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    assert(DIM==3);

    assert(!rCellPopulation.GetNode(targetNodeIndex)->IsBoundaryNode());
    assert(!rCellPopulation.GetNode(currentNodeIndex)->IsBoundaryNode());

    c_vector<double, DIM> displacement_direction = rCellPopulation.GetNode(targetNodeIndex)->rGetLocation() - rCellPopulation.GetNode(currentNodeIndex)->rGetLocation();
    displacement_direction /= norm_2(displacement_direction);

    PottsArbitrarySurfaceIn3DMesh* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh*> (&(rCellPopulation.rGetMesh()));
    CellData* lattice_cell_data = p_static_cast_potts_mesh->GetPottsLatticeSite(currentNodeIndex)->GetLatticeSiteData();

    c_vector<double,DIM> traction_direction;
    traction_direction(0) = lattice_cell_data->GetItem("applied_shear_stress_x");
    traction_direction(1) = lattice_cell_data->GetItem("applied_shear_stress_y");
    traction_direction(2) = lattice_cell_data->GetItem("applied_shear_stress_z");
    traction_direction /= norm_2(traction_direction);

    // Favour moves along the shear stress vector (3rd term ranges [0,2] hence the 0.5 factor)
    double delta_H = - mTractionCorrelationParameter * 0.5 * (1 + inner_prod(displacement_direction, traction_direction));

    return delta_H;
}

template<unsigned DIM>
double MechanotaxisPottsUpdateRule<DIM>::GetTractionCorrelationParameter()
{
    return mTractionCorrelationParameter;
}

template<unsigned DIM>
void MechanotaxisPottsUpdateRule<DIM>::SetTractionCorrelationParameter(double tractionCorrelationParameter)
{
    mTractionCorrelationParameter = tractionCorrelationParameter;
}


template<unsigned DIM>
void MechanotaxisPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<TractionCorrelationParameter>" << mTractionCorrelationParameter << "</TractionCorrelationParameter>\n";

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MechanotaxisPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MechanotaxisPottsUpdateRule<3>)

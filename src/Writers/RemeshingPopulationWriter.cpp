/*

Copyright (c) 2005-2020, University of Oxford.
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


#include "RemeshingPopulationWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"

RemeshingPopulationWriter::RemeshingPopulationWriter()
    : AbstractCellPopulationWriter<2, 3>("RemeshingPopulationData.dat")
{
}


void RemeshingPopulationWriter::Visit(HistoryDepMeshBasedCellPopulation<2, 3>* pCellPopulation)
{

    TRACE("Write is set up right")
    // HistoryDepMeshBasedCellPopulation<2, 3>* pCellPopulation = static_cast<HistoryDepMeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
   

//     std::map<unsigned,double > NewNodeToOldElementMap =  pCellPopulation->GetNodeToOldElementMap();
//     std::map<unsigned,double> NewNodeToOldElementDistanceMap = pCellPopulation->GetNewNodeToOldElementDistanceMap();
//     std::map<unsigned,c_vector<double,2>> MappingVariables_a_b = pCellPopulation->GetMappingVariables_a_b( );
//     std::map<unsigned, double> MappingVariables_alpha = pCellPopulation->GetMappingVariables_alpha( );

//     std::map<unsigned,c_vector<double,3>> MappingVariables_z_basis = pCellPopulation->GetMappingVariables_z_basis( );
//     std::map<unsigned,c_vector<double,3>> MappingVariables_PointInNewRef = pCellPopulation->GetMappingVariables_PointInNewRef( );
//     std::map<unsigned,c_vector<double,3>> MappingVariables_Difference = pCellPopulation->GetMappingVariables_Difference( );
//     std::map<unsigned,c_vector<double,3>> MappingVariables_P_Translated = pCellPopulation->GetMappingVariables_P_Translated( );
//     std::map<unsigned,c_vector<double,3>> MappingVariables_Cs = pCellPopulation->GetMappingVariables_Cs( );

//     std::map<std::vector<int>, std::vector<unsigned> > BinMap = pCellPopulation->GetAllBins();

//     // c_vector<c_vector<double, 2>, 3>  InitalPosition = GetInitalVectors(unsigned elem_index);
// // c_vector<c_vector<double, 3>, 2> HistoryDepMeshBasedCellPopulation<2, 3>::GetInitalShapeFunction(unsigned elem_index)
// // (*cell_iter)->GetCellData()->SetItem("Boundary", 0);
// // unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);
// // std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
// // Want everything in the bins
// // std::vector<int> HistoryDepMeshBasedCellPopulation<2, 3>::GetBin(c_vector<double, 3> Location)






//     for (map<unsigned, double>::iterator it= NewNodeToOldElementMap.begin(); it != NewNodeToOldElementMap.end(); it++)
//     {
//         std::cout << it->first    // string (key)
//                 << ':'
//                 << it->second   // string's value 
//                 << std::endl;
//     }

    
    
    //  for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin(); cell_iter != this->mCells.end(); ++cell_iter)
    // {
    //     //find the distance to the nearest neighbour
    //     (*cell_iter)->GetCellData()->SetItem("Curvature", 0);
    // }
    
    //
    // std::vector< c_vector<double, 3> > t1_swap_locations = pCellPopulation->rGetMesh().GetLocationsOfT1Swaps();

    // *this->mpOutStream << t1_swap_locations.size() << "\t";
    // for (unsigned index = 0;  index < t1_swap_locations.size(); index++)
    // {
    //     for (unsigned i=0; i<3; i++)
    //     {
    //         *this->mpOutStream << t1_swap_locations[index][i] << "\t";
    //     }
    // }

}



// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RemeshingPopulationWriter)
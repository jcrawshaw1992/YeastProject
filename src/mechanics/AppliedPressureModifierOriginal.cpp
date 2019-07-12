/*

Copyright (c) 2005-2014, University of Oxford.
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


#include "AppliedPressureModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "UblasCustomFunctions.hpp"


#include "Debug.hpp"
#include <math.h>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AppliedPressureModifier<ELEMENT_DIM,SPACE_DIM>::AppliedPressureModifier()
    : AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>(),
      mResetPressureOnCells(true),
      mPressureFile(""),
      mEdgeDivisionThreshold(DBL_MAX) //Defaults to no division
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AppliedPressureModifier<ELEMENT_DIM,SPACE_DIM>::~AppliedPressureModifier()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedPressureModifier<ELEMENT_DIM,SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
	TRACE("SetupSolve");
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	if(mResetPressureOnCells && mPressureFile !="")
	{
		LoadPressureFromFile();
		UpdateCellData(rCellPopulation);
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedPressureModifier<ELEMENT_DIM,SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
	assert(SPACE_DIM==3);

	// See if any edges are too long and if so divide them
	double num_cells = rCellPopulation.GetNumRealCells();

// TRACE("UpdateAtEndOfTimeStep");
	MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population = dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(rCellPopulation));
	assert(p_cell_population != NULL);

	p_cell_population->DivideLongSprings(mEdgeDivisionThreshold);

    //If there has been any remeshing recalculate the applied force from the Pressure file.
    if (num_cells != rCellPopulation.GetNumRealCells())
    {
        PRINT_2_VARIABLES(num_cells,rCellPopulation.GetNumRealCells());
        UpdateCellData(rCellPopulation);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedPressureModifier<ELEMENT_DIM,SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
	assert(SPACE_DIM==3); // Currently assumes that SPACE_DIM = 3
	std::map<unsigned, c_vector<unsigned, 2>  > LatticeToNodeMap;


	MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
	

	for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
		 cell_iter != rCellPopulation.End();
		 ++cell_iter)
	{
		
		
		c_vector<double, SPACE_DIM> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
		 unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
		//  PRINT_VARIABLE(node_index);
		 Node<SPACE_DIM>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);

		unsigned nearest_fluid_site = UNSIGNED_UNSET;
		long double distance_to_fluid_site = DBL_MAX;
		long double distance ;
		for (unsigned fluid_site_index = 0; fluid_site_index <  mAppliedPosition.size(); fluid_site_index++)
		{
			// Find the closest fluid site 
			distance = norm_2(location - mAppliedPosition[fluid_site_index]);
			if (distance < distance_to_fluid_site)
			{
				distance_to_fluid_site = distance;
				
				nearest_fluid_site = fluid_site_index;

			}

		}






		assert(distance_to_fluid_site<1e-3);
		assert(nearest_fluid_site != UNSIGNED_UNSET);

		long double Pressure =mAppliedPressure[nearest_fluid_site];
		cell_iter->GetCellData()->SetItem("PressureScalar", Pressure);
	

	}


}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AppliedPressureModifier<ELEMENT_DIM,SPACE_DIM>::GetResetPressureOnCells()
{
    return mResetPressureOnCells;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedPressureModifier<ELEMENT_DIM,SPACE_DIM>::SetResetPressureOnCells(bool resetPressureOnCells, std::string pressureFile)
{
	assert(!resetPressureOnCells || pressureFile!="");
	mResetPressureOnCells = resetPressureOnCells;
	mPressureFile = pressureFile;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AppliedPressureModifier<ELEMENT_DIM,SPACE_DIM>::GetEdgeDivisionThreshold()
{
    return mEdgeDivisionThreshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedPressureModifier<ELEMENT_DIM,SPACE_DIM>::SetEdgeDivisionThreshold(double edgeDivisionThreshold)
{
	mEdgeDivisionThreshold = edgeDivisionThreshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedPressureModifier<ELEMENT_DIM,SPACE_DIM>::LoadPressureFromFile()
{

		double InletX =  -5.95e-3;
		double MaxPressure = 0.03;

		double OutletX =  5.95e-3;
		double MinPressure =  0.0225;


	// PRINT_VARIABLE(mPressureFile);
	FILE* pressure_file = fopen((char*)mPressureFile.c_str(), "r");
	assert(pressure_file != NULL);
	
	hemelb::io::writers::xdr::XdrFileReader reader(pressure_file);
	

	// File format described in http://pauli.chem.ucl.ac.uk/trac/hemelb/wiki/ExpressureFiles

    unsigned hemelb_magic_number, expressure_magic_number, expressure_version_number;
	TRACE("Reading the Pressure file");
    reader.readUnsignedInt(hemelb_magic_number);
    reader.readUnsignedInt(expressure_magic_number);
    reader.readUnsignedInt(expressure_version_number);
	

    assert(hemelb_magic_number == 0x686c6221);
    assert(expressure_magic_number == 0x78747204);
    assert(expressure_version_number == 4);

    double voxel_size;
    reader.readDouble(voxel_size);

    c_vector<double,3> origin;
    reader.readDouble(origin[0]);
    reader.readDouble(origin[1]);
    reader.readDouble(origin[2]); // Dont know what this is 


    unsigned long long number_fluid_sites;
    reader.readUnsignedLong(number_fluid_sites);
    unsigned field_count;
    reader.readUnsignedInt(field_count);
    assert(field_count == 1); // Traction and tangetial component of pressure
    unsigned header_length;
    reader.readUnsignedInt(header_length);

    // Traction field header
    std::string field_name;
    unsigned number_floats;
    double pressure_offset;

    reader.readString(field_name, header_length);
    reader.readUnsignedInt(number_floats);
	// PRINT_VARIABLE(number_floats);
    assert(number_floats == 1); // dont think this will be a problem 
    reader.readDouble(pressure_offset);
	

    // Data section (we are reading a single timestep)
    unsigned long long timestep_num;
    reader.readUnsignedLong(timestep_num);

    mAppliedPosition.clear();
    mAppliedPressure.clear();
	
	c_vector<double, 3> FluidLatticeLocation;
    for (unsigned fluid_site_index = 0; fluid_site_index <  number_fluid_sites; fluid_site_index++)
    {
    	{
    		c_vector<unsigned,3> coords;
    		reader.readUnsignedInt(coords[0]);
    		reader.readUnsignedInt(coords[1]);
    		reader.readUnsignedInt(coords[2]);

			FluidLatticeLocation = origin+voxel_size*coords;
    		mAppliedPosition.push_back(FluidLatticeLocation);
    	}

    	{ 
			
			float pressure; 	
			reader.readFloat(pressure);
			pressure += pressure_offset;
	
			
			assert(fabs(pressure)<1e10);

			mAppliedPressure[fluid_site_index] = pressure;

    	}
    }

    assert(mAppliedPosition.size() == number_fluid_sites);
    assert(mAppliedPressure.size() == number_fluid_sites);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedPressureModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}


// Explicit instantiation
template class AppliedPressureModifier<1,1>;
template class AppliedPressureModifier<1,2>;
template class AppliedPressureModifier<2,2>;
template class AppliedPressureModifier<1,3>;
template class AppliedPressureModifier<2,3>;
template class AppliedPressureModifier<3,3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AppliedPressureModifier)


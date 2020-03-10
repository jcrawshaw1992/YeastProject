#include "HemodynamicPressure.hpp"
#include "UblasCustomFunctions.hpp"
#include "Debug.hpp"
// Area constrant on a membrane 


HemodynamicPressure::HemodynamicPressure()
   : AbstractForce<2,3>()
{
}


void HemodynamicPressure::LoadPressureFromFile(std::string PressureFile)
{
    mPressureFile = PressureFile;
    PRINT_VARIABLE(mPressureFile);
	FILE* Pressure_file = fopen((char*)mPressureFile.c_str(), "r");
	assert(Pressure_file != NULL);
    
    hemelb::io::writers::xdr::XdrFileReader reader(Pressure_file);

// File format described in http://pauli.chem.ucl.ac.uk/trac/hemelb/wiki/ExtractionFiles

    unsigned hemelb_magic_number, extraction_magic_number, extraction_version_number;
	TRACE("Reading the Pressure file");
    reader.readUnsignedInt(hemelb_magic_number);

    reader.readUnsignedInt(extraction_magic_number);

    reader.readUnsignedInt(extraction_version_number);
	

    assert(hemelb_magic_number == 0x686c6221);
    assert(extraction_magic_number == 0x78747204);
    assert(extraction_version_number == 4);

    double voxel_size;
    reader.readDouble(voxel_size);

    c_vector<double,3> origin;
    PRINT_VECTOR(origin);
    reader.readDouble(origin[0]);
    reader.readDouble(origin[1]);
    reader.readDouble(origin[2]);


    unsigned long long number_fluid_sites;
    reader.readUnsignedLong(number_fluid_sites);
    unsigned field_count;
    reader.readUnsignedInt(field_count);
    PRINT_VARIABLE(field_count)
    assert(field_count == 1); // Pressure is a scalar
    unsigned header_length;

    reader.readUnsignedInt(header_length);

    // Traction field header
    std::string field_name;
    unsigned number_floats;
    double pressure_offset;

    reader.readString(field_name, header_length);
    reader.readUnsignedInt(number_floats);
    assert(number_floats == 1);
    reader.readDouble(pressure_offset);



    PRINT_2_VARIABLES(pressure_offset, number_fluid_sites);

    // Tangential traction field header
    // double tanget_traction_offset;

    // reader.readString(field_name, header_length);
    // reader.readUnsignedInt(number_floats);
    // assert(number_floats == 3);
    // reader.readDouble(tanget_traction_offset);

    // Data section (we are reading a single timestep)
    unsigned long long timestep_num;
    reader.readUnsignedLong(timestep_num);

    mAppliedPosition.clear();
    mAppliedPressures.clear();
    // mAppliedTangentTractions.clear();

    for (unsigned fluid_site_index = 0; fluid_site_index <  number_fluid_sites; fluid_site_index++)
    {
    	{
    		c_vector<unsigned,3> coords;
    		reader.readUnsignedInt(coords[0]);
    		reader.readUnsignedInt(coords[1]);
    		reader.readUnsignedInt(coords[2]);

    		mAppliedPosition.push_back(origin+voxel_size*coords);
    	}

    	{
			c_vector<float,3> pressure;
			reader.readFloat(pressure[0]);
			pressure[0] += pressure_offset;
			reader.readFloat(pressure[1]);
			pressure[1] += pressure_offset;
			reader.readFloat(pressure[2]);
			pressure[2] += pressure_offset;

			assert(fabs(pressure[0])<1e10);
			assert(fabs(pressure[1])<1e10);
			assert(fabs(pressure[2])<1e10);

			mAppliedPressures.push_back(pressure);
    	}

    }

    assert(mAppliedPosition.size() == number_fluid_sites);
    assert(mAppliedPressures.size() == number_fluid_sites);

}










void HemodynamicPressure::AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation)
{
/*
 * Itterate over all the cells, get the force, then add the force .
*/
    MeshBasedCellPopulation<2,3>* p_cell_population = static_cast<MeshBasedCellPopulation<2,3>*>(&rCellPopulation);


    for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
	{
		c_vector<double, 3> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);

		unsigned nearest_fluid_site = UNSIGNED_UNSET;
		double distance_to_fluid_site = DBL_MAX;
        // PRINT_2_VARIABLES(nearest_fluid_site, distance_to_fluid_site );
 
        // Iterate over all the LB positions
		for (unsigned fluid_site_index = 0; fluid_site_index <  mAppliedPosition.size(); fluid_site_index++)
		{

            // LOOKS LIKE YOU ITERATE OVER ALL THE FLUID SITES AND ACCEPT THE ONE WITH THE LEST DISTANCE
			double distance = norm_2(location - mAppliedPosition[fluid_site_index]);
            // PRINT_VARIABLE(distance);
			if (distance < distance_to_fluid_site)
			{
				distance_to_fluid_site = distance;
				nearest_fluid_site = fluid_site_index;
			}

		}
		 assert(nearest_fluid_site != UNSIGNED_UNSET);


		// Calculate the approximate area of the voronoi region around the cell by including a third of the area
		// of all surrounding triangles. Useful for turning stresses into forces.
		double voronoi_cell_area = rCellPopulation.GetVolumeOfCell(*cell_iter);
		// PRINT_VARIABLE(voronoi_cell_area);
        c_vector<double,3> appliedpressure = mAppliedPressures[nearest_fluid_site];
        // PRINT_VECTOR(appliedpressure);
		c_vector<double,3> force = 1* mAppliedPressures[nearest_fluid_site]*voronoi_cell_area;
		// PRINT_VECTOR(force);
        PRINT_VECTOR(force);

		// assert(fabs(force[0])<1e10);
		// assert(fabs(force[1])<1e10);
		// assert(fabs(force[2])<1e10);
		cell_iter->GetCellData()->SetItem("force_x", force[0]);
        cell_iter->GetCellData()->SetItem("force_y", force[1]);
        cell_iter->GetCellData()->SetItem("force_z", force[2]);
        

        double normal_force = norm_2(force);
        // PRINT_VARIABLE(normal_force);

        cell_iter->GetCellData()->SetItem("NormalForce", normal_force);
        pNode->AddAppliedForceContribution(force);
        // if (node_index == 100)
        // {
        //     double normal_force = norm_2(force);
        //     PRINT_VECTOR(force);

		// 	PRINT_VARIABLE(normal_force);
        // }

	}

} 




void HemodynamicPressure::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2,3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(HemodynamicPressure)



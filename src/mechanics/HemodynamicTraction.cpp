#include "HemodynamicTraction.hpp"
#include "UblasCustomFunctions.hpp"
#include "Debug.hpp"
// Area constrant on a membrane 



HemodynamicTraction::HemodynamicTraction()
   : AbstractForce<2,3>()
{
}


void HemodynamicTraction::LoadTractionFromFile(std::string TractionFile)
{
     mTractionFile = TractionFile;
    PRINT_VARIABLE(mTractionFile);
	FILE* Traction_file = fopen((char*)mTractionFile.c_str(), "r");
	assert(Traction_file != NULL);
    
    hemelb::io::writers::xdr::XdrFileReader reader(Traction_file);

// File format described in http://pauli.chem.ucl.ac.uk/trac/hemelb/wiki/ExtractionFiles

    unsigned hemelb_magic_number, extraction_magic_number, extraction_version_number;
	TRACE("Reading the Traction file");
    reader.readUnsignedInt(hemelb_magic_number);

    reader.readUnsignedInt(extraction_magic_number);

    reader.readUnsignedInt(extraction_version_number);
	

    assert(hemelb_magic_number == 0x686c6221);
    assert(extraction_magic_number == 0x78747204);
    assert(extraction_version_number == 4);

  
    double voxel_size;
    reader.readDouble(voxel_size);

    c_vector<double,3> origin;
    reader.readDouble(origin[0]);
    reader.readDouble(origin[1]);
    reader.readDouble(origin[2]);


    unsigned long long number_fluid_sites;
    reader.readUnsignedLong(number_fluid_sites);
    unsigned field_count;
    reader.readUnsignedInt(field_count);
    assert(field_count == 2); // Traction and tangetial component of traction
    unsigned header_length;
    reader.readUnsignedInt(header_length);

    // Traction field header
    std::string field_name;
    unsigned number_floats;
    double traction_offset;

    reader.readString(field_name, header_length);
    reader.readUnsignedInt(number_floats);
    assert(number_floats == 3);
    reader.readDouble(traction_offset);

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
    mAppliedTractions.clear();
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
			c_vector<float,3> traction;
			reader.readFloat(traction[0]);
			traction[0] += traction_offset;
			reader.readFloat(traction[1]);
			traction[1] += traction_offset;
			reader.readFloat(traction[2]);
			traction[2] += traction_offset;

			// assert(fabs(traction[0])<1e10);
			// assert(fabs(traction[1])<1e10);
			// assert(fabs(traction[2])<1e10);

			mAppliedTractions.push_back(traction);
    	}

    	// {
		// 	c_vector<float,3> tangent_traction;
		// 	reader.readFloat(tangent_traction[0]);
		// 	tangent_traction[0] += tanget_traction_offset;
		// 	reader.readFloat(tangent_traction[1]);
		// 	tangent_traction[1] += tanget_traction_offset;
		// 	reader.readFloat(tangent_traction[2]);
		// 	tangent_traction[2] += tanget_traction_offset;

		// 	assert(fabs(tangent_traction[0])<1e10);
		// 	assert(fabs(tangent_traction[1])<1e10);
		// 	assert(fabs(tangent_traction[2])<1e10);

		// 	mAppliedTangentTractions.push_back(tangent_traction);
    	// }
    }

    // assert(mAppliedPosition.size() == number_fluid_sites);
    // assert(mAppliedTractions.size() == number_fluid_sites);
    // assert(mAppliedTangentTractions.size() == number_fluid_sites);

}










void HemodynamicTraction::AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation)
{
/*
 * Itterate over all the cells, get the force, then add the force .
*/

        TRACE("Finding tangental and normal shear stress");
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

            // LOOKS LIKE YOU ITERATE OVER ALL THE FLUID SITES AND ACCEPT THE ONE WITH THE LEAST DISTANCE
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
		// double voronoi_cell_area = rCellPopulation.GetVolumeOfCell(*cell_iter);
		// PRINT_VARIABLE(voronoi_cell_area);

		c_vector<double,3> TractionForce =  mAppliedTractions[nearest_fluid_site] * rCellPopulation.GetVolumeOfCell(*cell_iter) ;
        // c_vector<double,3> ShearFroce = mAppliedTangentTractions[nearest_fluid_site];
        // c_vector<double,3> Force =  mAppliedTractions[nearest_fluid_site]
        
		PRINT_VECTOR(TractionForce);

		// assert(fabs(force[0])<1e10);
		// assert(fabs(force[1])<1e10);
		// assert(fabs(force[2])<1e10);
		cell_iter->GetCellData()->SetItem("Traction_x", TractionForce[0]);
        cell_iter->GetCellData()->SetItem("Traction_y", TractionForce[1]);
        cell_iter->GetCellData()->SetItem("Traction_z", TractionForce[2]);


        // cell_iter->GetCellData()->SetItem("Shear_x", mAppliedTangentTractions[nearest_fluid_site][0]);
        // cell_iter->GetCellData()->SetItem("Shear_y", mAppliedTangentTractions[nearest_fluid_site][1]);
        // cell_iter->GetCellData()->SetItem("Shear_z", mAppliedTangentTractions[nearest_fluid_site][2]);




        double normal_TractionForce = norm_2(TractionForce);

        cell_iter->GetCellData()->SetItem("ForceNormal", normal_TractionForce);



        pNode->AddAppliedForceContribution(TractionForce);
        // if (node_index == 100)
        // {
        //     double normal_TractionForce = norm_2(TractionForce);
        //     PRINT_VECTOR(TractionForce);

		// 	PRINT_VARIABLE(normal_TractionForce);
        // }

	}





//     MeshBasedCellPopulation<2,3>* p_cell_population = static_cast<MeshBasedCellPopulation<2,3>*>(&rCellPopulation);

//    for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
//          cell_iter != rCellPopulation.End();
//          ++cell_iter)
//     {
//         unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
//         Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);
//         // c_vector<double,3> force = 



//     //    force[0] =  4 * cell_iter->GetCellData()->GetItem("applied_force_x");
// 	//    force[1] =  4 * cell_iter->GetCellData()->GetItem("applied_force_y");
// 	//    force[2] =  4 * cell_iter->GetCellData()->GetItem("applied_force_z");
//     // //    PRINT_VECTOR(force);
//     //    TRACE("AM I EVEN HERE");

//         // pNode->AddAppliedForceContribution(force);
//         // TRACE("Drag corrected");
//         //  PRINT_2_VARIABLES(NormOldForce , NormNewForce);
//      }
} 




void HemodynamicTraction::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2,3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(HemodynamicTraction)



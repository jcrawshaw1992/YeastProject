
#include "TractionDataLoader.hpp"
#include "XdrFileReader.hpp"
#include "Debug.hpp"

TractionDataLoader::TractionDataLoader(const std::string& tractionFilename)
{
	FILE* traction_file = fopen((char*)tractionFilename.c_str(), "r");
	assert(traction_file != NULL);
	hemelb::io::writers::xdr::XdrFileReader reader(traction_file);

	// File format described in http://pauli.chem.ucl.ac.uk/trac/hemelb/wiki/ExtractionFiles

    unsigned hemelb_magic_number, extraction_magic_number, extraction_version_number;
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
    double tanget_traction_offset;

    reader.readString(field_name, header_length);
    reader.readUnsignedInt(number_floats);
    assert(number_floats == 3);
    reader.readDouble(tanget_traction_offset);

    // Data section (we are reading a single timestep)
    unsigned long long timestep_num;
    reader.readUnsignedLong(timestep_num);

    mAppliedPosition.clear();
    mAppliedTractions.clear();

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
			// PRINT_VECTOR(traction);
			assert(fabs(traction[0])<1e10);
			assert(fabs(traction[1])<1e10);
			assert(fabs(traction[2])<1e10);

			mAppliedTractions.push_back(traction);
    	}

    	{
			c_vector<float,3> tangent_traction;
			reader.readFloat(tangent_traction[0]);
			tangent_traction[0] += tanget_traction_offset;
			reader.readFloat(tangent_traction[1]);
			tangent_traction[1] += tanget_traction_offset;
			reader.readFloat(tangent_traction[2]);
			tangent_traction[2] += tanget_traction_offset;

			assert(fabs(tangent_traction[0])<1e10);
			assert(fabs(tangent_traction[1])<1e10);
			assert(fabs(tangent_traction[2])<1e10);

			mAppliedTangentTractions.push_back(tangent_traction);
    	}
    }

    assert(mAppliedPosition.size() == number_fluid_sites);
    assert(mAppliedTractions.size() == number_fluid_sites);
    assert(mAppliedTangentTractions.size() == number_fluid_sites);
}

void TractionDataLoader::UpdateLatticeSiteData(PottsArbitrarySurfaceIn3DMesh<3> * pPottsMesh)
{
    for (PottsArbitrarySurfaceIn3DMesh<3>::NodeIterator iter = pPottsMesh->GetNodeIteratorBegin();
         iter != pPottsMesh->GetNodeIteratorEnd();
         ++iter)
    {
    	if (!iter->IsBoundaryNode())
    	{
			c_vector<double, 3> location = iter->rGetLocation();

			unsigned nearest_fluid_site = UNSIGNED_UNSET;
			double distance_to_fluid_site = DBL_MAX;
			for (unsigned fluid_site_index = 0; fluid_site_index <  mAppliedPosition.size(); fluid_site_index++)
			{
				double distance = norm_2(location - mAppliedPosition[fluid_site_index]);
				if (distance < distance_to_fluid_site)
				{
					distance_to_fluid_site = distance;
					nearest_fluid_site = fluid_site_index;
				}

			}
			assert(nearest_fluid_site != UNSIGNED_UNSET);


			// Calculate the approximate area of the voronoi region around the cell by including a third of the area
			// of all surrounding triangles. Useful for turning stresses into forces.
			double lattice_cell_area = pPottsMesh->GetVolumeOfLatticeSite(iter->GetIndex());

			c_vector<double,3> force = mAppliedTractions[nearest_fluid_site]*lattice_cell_area;
			c_vector<double,3> shear_stress = mAppliedTangentTractions[nearest_fluid_site];

			assert(fabs(force[0])<1e10);
			assert(fabs(force[1])<1e10);
			assert(fabs(force[2])<1e10);
			assert(fabs(shear_stress[0])<1e10);
			assert(fabs(shear_stress[1])<1e10);
			assert(fabs(shear_stress[2])<1e10);
			PRINT_VARIABLE(iter->GetIndex());

			// Store the force in CellData -- i think this might need a static cast 
			CellData* p_lattice_site_data = pPottsMesh->GetPottsLatticeSite(iter->GetIndex())->GetLatticeSiteData();

			// p_lattice_site_data->SetItem("applied_force_x", force[0]);
			// p_lattice_site_data->SetItem("applied_force_y", force[1]);
			// p_lattice_site_data->SetItem("applied_force_z", force[2]);
			// p_lattice_site_data->SetItem("applied_force_mag", norm_2(force));
			// p_lattice_site_data->SetItem("lattice_cell_area", lattice_cell_area);
			// p_lattice_site_data->SetItem("applied_shear_stress_x", shear_stress[0]);
			// p_lattice_site_data->SetItem("applied_shear_stress_y", shear_stress[1]);
			// p_lattice_site_data->SetItem("applied_shear_stress_z", shear_stress[2]);
			// p_lattice_site_data->SetItem("applied_shear_stress_mag", norm_2(shear_stress));
    	}
    }

}



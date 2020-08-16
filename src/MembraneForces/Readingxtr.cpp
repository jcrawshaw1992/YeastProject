template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::LoadTractionFromFile()
{
    std::string TractionFile = mHemeLBDirectory + "results/Extracted/surface-tractions.xtr";
    // std::string TractionFile = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexus/results/Extracted/surface-tractions.xtr";
	TRACE("Load tracrtion file");
	// PRINT_VARIABLE(TractionFile); 
	FILE* traction_file = fopen((char*)TractionFile.c_str(), "r");
	assert(traction_file != NULL);
	
	hemelb::io::writers::xdr::XdrFileReader reader(traction_file);
	
	// File format described in http://pauli.chem.ucl.ac.uk/trac/hemelb/wiki/ExtractionFiles

    unsigned hemelb_magic_number, extraction_magic_number, extraction_version_number;
	// TRACE("Reading the traction file");
    reader.readUnsignedInt(hemelb_magic_number);
    reader.readUnsignedInt(extraction_magic_number);
    reader.readUnsignedInt(extraction_version_number);
	

    assert(hemelb_magic_number == 0x686c6221);
    assert(extraction_magic_number == 0x78747204);
    assert(extraction_version_number == 4);

    double voxel_size;
    reader.readDouble(voxel_size);
    // PRINT_VARIABLE(voxel_size)

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
    // PRINT_VARIABLE(timestep_num)

    mAppliedPosition.clear();
    mAppliedTractions.clear();
    mAppliedTangentTractions.clear();

    for (unsigned fluid_site_index = 0; fluid_site_index <  number_fluid_sites; fluid_site_index++)
    {
    	{
    		c_vector<unsigned,3> coords;
    		reader.readUnsignedInt(coords[0]);
    		reader.readUnsignedInt(coords[1]);
    		reader.readUnsignedInt(coords[2]);
    	
            mAppliedPosition.push_back((origin+coords*voxel_size)*1e3); // mAppliedPosition.push_back(origin+voxel_size*coords); <---- this was the original line, its not good, give completluy wrong results, it looks like a time based problem, may need to revist 
            // c_vector<double, 3> C00rds =voxel_size*coords*1e3;
            // PRINT_VECTOR(origin)
            // PRINT_VECTOR(C00rds)

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

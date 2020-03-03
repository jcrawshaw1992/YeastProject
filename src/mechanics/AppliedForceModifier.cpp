
#include "AppliedForceModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "UblasCustomFunctions.hpp"

#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"


//  #include "MembraneShearForce.hpp"

#include "Debug.hpp"
#include <math.h>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::AppliedForceModifier()
    : AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>(),
      mResetTractionsOnCells(true),
      mTractionFile(""),
      mEdgeDivisionThreshold(DBL_MAX) //Defaults to no division
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::~AppliedForceModifier()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{

	TRACE("In the set up solve");
	// This is an empty function because the mutatant elements need to be identified external to the simulation so that the membrane properties 
	// can then also be updated before the simulation. 
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	if(mResetTractionsOnCells && mTractionFile !="")
	{
		LoadTractionFromFile();
		UpdateCellData(rCellPopulation);
	}
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::SetupVessel(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	// if(mResetTractionsOnCells && mTractionFile !="")
	// {
		
		LoadTractionFromFile();
		TRACE("Update the applied forces. This should once at the start of each iteration");
		UpdateCellData(rCellPopulation);
		// TRACE("Hi Jess");
	// }
}





template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
	// This part is redundant because there is no remeshing 
	assert(SPACE_DIM==3);
		// See if any edges are too long and if so divide them
	double num_cells = rCellPopulation.GetNumRealCells();
	
	// TRACE("UpdateAtEndOfTimeStep");
	MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population = dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(rCellPopulation));
	assert(p_cell_population != NULL);
	

	p_cell_population->DivideLongSprings(mEdgeDivisionThreshold);

    //If there has been any remeshing recalculate the applied force from the traction file.
    if (num_cells != rCellPopulation.GetNumRealCells())
    {
        PRINT_2_VARIABLES(num_cells,rCellPopulation.GetNumRealCells());
        UpdateCellData(rCellPopulation);
		TRACE("I think this should only happen if there is remeshing");
		// UpdateCellDataBeforeSimulation(rCellPopulation);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
	
	assert(SPACE_DIM==3); // Currently assumes that SPACE_DIM = 3
	
	std::map<unsigned, c_vector<unsigned, 2>  > LatticeToNodeMap;

	// std::map<unsigned, c_vector<double, 3>  > ForceMap;
	
	
	MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);

	// c_vector<double, SPACE_DIM> centroid = zero_vector<double>(3);
	// for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
	// 		cell_iter != rCellPopulation.End();
	// 		++cell_iter)
	// {
	// 	unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
	// 	centroid += rCellPopulation.GetNode(node_index)->rGetLocation();
	// 	// Area += rCellPopulation.GetVolumeOfCell(*cell_iter);
	// }
	// centroid /= rCellPopulation.GetNumRealCells();


	for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
		 cell_iter != rCellPopulation.End();
		 ++cell_iter)
	{
		c_vector<double, SPACE_DIM> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
		unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
		Node<SPACE_DIM>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);
		unsigned nearest_fluid_site = UNSIGNED_UNSET;
		double distance_to_fluid_site = DBL_MAX;
	
	
		for (unsigned fluid_site_index = 0; fluid_site_index <  mAppliedPosition.size(); fluid_site_index++)
		{
			// Find the closest fluid site 
			double distance = norm_2(location - mAppliedPosition[fluid_site_index]);
			if (distance < distance_to_fluid_site)
			{
				distance_to_fluid_site = distance;	
				nearest_fluid_site = fluid_site_index;
			}
		}
		// PRINT_VARIABLE(distance_to_fluid_site);
		// assert(distance_to_fluid_site<1e-3);
		LatticeToNodeMap[node_index] = Create_c_vector(nearest_fluid_site, 0 );
		assert(nearest_fluid_site != UNSIGNED_UNSET);


		/* --------------
			 Normals
		   --------------

		   We are defining the normal to the node as the distance from the center line. 
		   This is being done as the pressure gradient distorts a basic cylindrical shape 
		   into a conical shape. Once this has started to happen, the surface normal now has 
		   some compoent in the axial (or z) direction causing the nodes to be pushed downwards 
		   and distorting the geometry in an unrealistic manner. 

		   Biological justification -- Vessel (at least the arteries and veins) are held in an elongated state 
		   by the surrounding tissue. The must be some tethering to the tissue for this to occur, 
		   as such the axial force from the  fluid flow will be opposed by an equal and opposite force
		   from the tissue. I might consider how to code this up, or just leave it as it is 


		   ----------------
		    Jess needs to 
		   ----------------

		   I need to adjust this code to read in the center lines file, find the nearest point on the center line 
		   and take the radial normal as the tangent to this line. 
		*/

		// c_vector<double, 3> cell_location = pNode->rGetLocation() - centroid;
		// cell_location(2) = 0.0;
		// cell_location /= norm_2(cell_location);

		// c_vector<long double, 3> NormalVector = cell_location; 

		c_vector<long double, 3> NormalVector = Create_c_vector(0,0,0);
		std::set<unsigned>& containing_elements = pNode->rGetContainingElementIndices();
        assert(containing_elements.size() > 0);
		// Finding the normal to the node -- need to check this 
        for (std::set<unsigned>::iterator iter = containing_elements.begin();
            iter != containing_elements.end();
            ++iter)
        {
            Node<SPACE_DIM>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
            Node<SPACE_DIM>* pNode1 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
            Node<SPACE_DIM>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

            c_vector<long double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
            c_vector<long double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

            NormalVector  += VectorProduct(vector_12, vector_13);
        }
		NormalVector /=norm_2(NormalVector);

		// Get the HemeLB force at the closest lattice site 
		c_vector<long double,3> force = mAppliedTractions[nearest_fluid_site]/133.3223874;//*voronoi_cell_area;
		long double Pressure = norm_2(force);     


		// Check I have maintained outward pointing norm by getting the direction of the force 
		// Vector and making sure the normal and the force vector are goin in the same direcvtion 
		// DO this by checking the angle betweem these two vectors is below a certain value -- basic dot proudct thing

		// c_vector<long double,3> forceDirection = force / Pressure;

		// double Angle = abs(acos(inner_prod(forceDirection, NormalVector) )) ;

		// if (Angle > M_PI/2)
		// {
		// 	// TRACE(" Normal was the wrong way");
		// 	NormalVector = -NormalVector;

		// }



		
		
		// Calculate the approximate area of the voronoi region around the cell by including a third of the area
		// of all surrounding triangles. Useful for turning stresses into forces.
		double voronoi_cell_area = rCellPopulation.GetVolumeOfCell(*cell_iter);
		

		// PRINT_VARIABLE(Pressure);
		// // XXXX  -Replaced the direction of the force with the normal -- did this because the force acts normal to the lattice site it was selected from, which is not the identical position to the node, but is slightly off
		// location = location - centroid;
		// location /= norm_2(location);
		
		c_vector<long double,3> Force = Pressure * NormalVector; 

		c_vector<double,3> shear_stress = mAppliedTangentTractions[nearest_fluid_site];

		assert(fabs(Force[0])<1e10);
		assert(fabs(Force[1])<1e10);
		assert(fabs(Force[2])<1e10);
		assert(fabs(shear_stress[0])<1e10);
		assert(fabs(shear_stress[1])<1e10);
		assert(fabs(shear_stress[2])<1e10);


		// Store the force in CellData

		cell_iter->GetCellData()->SetItem("Pressure", Pressure);
		// PRINT_VARIABLE(Pressure);


		// TRACE("Setting Force");
		cell_iter->GetCellData()->SetItem("applied_force_x", Force[0]);
		// PRINT_VECTOR(Force);
		cell_iter->GetCellData()->SetItem("applied_force_y", Force[1]);
		cell_iter->GetCellData()->SetItem("applied_force_z", Force[2]);
		cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(Force));
		cell_iter->GetCellData()->SetItem("voronoi_cell_area", voronoi_cell_area);
		cell_iter->GetCellData()->SetItem("applied_shear_stress_x", shear_stress[0]);
		cell_iter->GetCellData()->SetItem("applied_shear_stress_y", shear_stress[1]);
		cell_iter->GetCellData()->SetItem("applied_shear_stress_z", shear_stress[2]);
		cell_iter->GetCellData()->SetItem("applied_shear_stress_mag", norm_2(shear_stress));
	}


}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::GetResetTractionsOnCells()
{
    return mResetTractionsOnCells;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::SetResetTractionsOnCells(bool resetTractionsOnCells, std::string tractionFile)
{
	assert(!resetTractionsOnCells || tractionFile!="");
	mResetTractionsOnCells = resetTractionsOnCells;
	mTractionFile = tractionFile;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::GetEdgeDivisionThreshold()
{
    return mEdgeDivisionThreshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::SetEdgeDivisionThreshold(double edgeDivisionThreshold)
{
	mEdgeDivisionThreshold = edgeDivisionThreshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::LoadTractionFromFile()
{
	TRACE("Load tracrtion file");
	// PRINT_VARIABLE(mTractionFile);
	FILE* traction_file = fopen((char*)mTractionFile.c_str(), "r");
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
    mAppliedTangentTractions.clear();

    for (unsigned fluid_site_index = 0; fluid_site_index <  number_fluid_sites; fluid_site_index++)
    {
    	{
    		c_vector<unsigned,3> coords;
    		reader.readUnsignedInt(coords[0]);
    		reader.readUnsignedInt(coords[1]);
    		reader.readUnsignedInt(coords[2]);

    		mAppliedPosition.push_back(origin+voxel_size*coords);
			// c_vector<long double,3> LatticeSite= origin+voxel_size*coords;
		// PRINT_VECTOR(LatticeSite);
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}


// Explicit instantiation
template class AppliedForceModifier<1,1>;
template class AppliedForceModifier<1,2>;
template class AppliedForceModifier<2,2>;
template class AppliedForceModifier<1,3>;
template class AppliedForceModifier<2,3>;
template class AppliedForceModifier<3,3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AppliedForceModifier)



// if (cell_iter->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>() )
        // {
        //   //Dont set the maps, use cell data insead  // mElasticShearModulusMap[node_index] = mElasticShearModulusB;
            
        // }
        // else
        // {
        //     mElasticShearModulusMap[node_index] = mElasticShearModulus;
        //     mAreaDilationModulusMap[node_index] = mAreaDilationModulus;
        // }




	// for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
	// 	 cell_iter != rCellPopulation.End();
	// 	 ++cell_iter)
	// {

	// 	 unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

	// 	 for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter2 = rCellPopulation.Begin();
	// 	 cell_iter2 != rCellPopulation.End();
	// 	 ++cell_iter2)
	// 			{
	// 				unsigned Othernode_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter2);

	// 					unsigned LatticeSite0  = LatticeToNodeMap[node_index][0];
	// 					unsigned LatticeSite1  = LatticeToNodeMap[Othernode_index][0];

	// 				if (LatticeSite0 == LatticeSite1)
	// 				{

	// 					LatticeToNodeMap[node_index][1] = LatticeToNodeMap[node_index][1]+1;

	// 				}
					
	// 				if (Othernode_index == node_index)
	// 				{
	// 						LatticeToNodeMap[node_index][1] = LatticeToNodeMap[node_index][1]-1;
	// 				}

	// 			}

	// 	PRINT_3_VARIABLES(node_index, LatticeToNodeMap[node_index], ForceMap[node_index]);


	// 	cell_iter->GetCellData()->SetItem("SharingLattice", LatticeToNodeMap[node_index][1]);

	// }





	// Normals that didnt work



		
		// // Finding the normal to the node -- need to check this 
        // assert(containing_elements.size() > 0);
        // for (std::set<unsigned>::iterator iter = containing_elements.begin();
        //     iter != containing_elements.end();
        //     ++iter)
        // {
        //     Node<SPACE_DIM>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
        //     Node<SPACE_DIM>* pNode1 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
        //     Node<SPACE_DIM>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

        //     c_vector<long double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
        //     c_vector<long double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

        //     NormalVector  += VectorProduct(vector_12, vector_13);
        // }
		// NormalVector /=norm_2(NormalVector);
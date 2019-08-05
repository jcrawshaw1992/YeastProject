#include "AppliedForceOffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "Debug.hpp"
#include <math.h>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AppliedForceOffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::AppliedForceOffLatticeSimulation(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                      std::string tractionFile,
                                                      double edgeDivisionThreshold,
                                                      double shearThreshold,
                                                      bool deleteCellPopulationInDestructor,
                                                      bool initialiseCells)
    : OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells),
      mResetTractionsOnCells(true),
      mTractionFile(tractionFile),
      mEdgeDivisionThreshold(edgeDivisionThreshold),
      mShearThreshold(shearThreshold),
      mMaxDivisions(500) // if to large can cause elements to become degenerate
{
	if(tractionFile !="") // if loading from archive don't want to reload tractions here.
	{
		LoadTractionFromFile();
	}

	MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population = dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation));
    // Loop over all edges and work out whether rest length needs adjusting based on shear stress
    for (typename MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator spring_iterator = p_cell_population->SpringsBegin();
         spring_iterator != p_cell_population->SpringsEnd();
         ++spring_iterator)
    {
        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();



        std::pair<unsigned, unsigned> node_pair =  this->mrCellPopulation.CreateOrderedPair(nodeA_global_index, nodeB_global_index);

        mDivisionsApplied[node_pair] = 0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AppliedForceOffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::~AppliedForceOffLatticeSimulation()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceOffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::SetupSolve()
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
	if(mResetTractionsOnCells)
	{
		LoadTractionFromFile();
		UpdateCellData();
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceOffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::UpdateAtEndOfTimeStep()
{
	assert(SPACE_DIM==3);

	// See if any edges are too long and if so divide them
	double num_cells = this->mrCellPopulation.GetNumRealCells();
    
	MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population = dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation));
	assert(p_cell_population != NULL);

	p_cell_population->DivideLongSprings(mEdgeDivisionThreshold);

    //If there has been any remeshing recalculate the applied force from the traction file.
    if (num_cells != this->mrCellPopulation.GetNumRealCells())
    {
        PRINT_2_VARIABLES(num_cells,this->mrCellPopulation.GetNumRealCells());
        UpdateCellData();
    }
        
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    if (p_simulation_time->GetTimeStepsElapsed()%this->mSamplingTimestepMultiple == 0)
	{
    	PRINT_VARIABLE(p_simulation_time->GetTime());
	}

//    // Loop over all edges and work out whether rest length needs adjusting based on shear stress
//    for (typename MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator spring_iterator = p_cell_population->SpringsBegin();
//         spring_iterator != p_cell_population->SpringsEnd();
//         ++spring_iterator)
//    {
//        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
//        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
//
//        double nodeA_shear_stress = p_cell_population->GetCellUsingLocationIndex(nodeA_global_index)->GetCellData()->GetItem("applied_shear_stress_mag");
//        double nodeB_shear_stress = p_cell_population->GetCellUsingLocationIndex(nodeB_global_index)->GetCellData()->GetItem("applied_shear_stress_mag");
//
//        // TODO: To be made configurable
//        double reduction_factor = 0.0005;
//        if ((nodeA_shear_stress<mShearThreshold) && (nodeB_shear_stress<mShearThreshold))
//        {
//        	unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
//            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
//
//            std::pair<unsigned, unsigned> node_pair =  this->mrCellPopulation.CreateOrderedPair(nodeA_global_index, nodeB_global_index);
//
//            if (mDivisionsApplied.at(node_pair) < mMaxDivisions)
//            {
//            	double previous_rest_length = p_cell_population->GetRestLength(nodeA_global_index, nodeB_global_index);
//            	double new_rest_length = (1 - reduction_factor) * previous_rest_length;
//            	p_cell_population->SetRestLength(nodeA_global_index, nodeB_global_index, new_rest_length);
//            	mDivisionsApplied[node_pair]++;
//            }
//        }
//    }


	// TODO: To be made configurable
	double reduction_factor = 0.0005;

    // Loop over elements and calculate the new rest lengths
    for (typename MutableMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
             elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
             ++elem_iter)
	{
    	//TODO get this from cell data note this is unit vector.
    	c_vector<double, SPACE_DIM> tangent_vector;
    	tangent_vector[0] = 0.0;
    	tangent_vector[1] = 0.0;
    	tangent_vector[2] = 1.0;

    	tangent_vector/=norm_2(tangent_vector);

    	// TODO get from cell data and only shrink if below thereshold!
		//cell_iter->GetCellData()->GetItem("applied_shear_stress_mag");
		double shear_stress_magnitude = 0.0;

    	// Make sure tangent_vector lies in plane
    	c_vector<double,SPACE_DIM> element_normal = elem_iter->CalculateNormal();

    	tangent_vector -= element_normal * inner_prod(tangent_vector,element_normal);
    	tangent_vector/=norm_2(tangent_vector);


    	// Get nodes of element
    	assert(elem_iter->GetNumNodes() == 3);
    	Node<SPACE_DIM>* node_zero = elem_iter->GetNode(0);
    	Node<SPACE_DIM>* node_one = elem_iter->GetNode(1);
    	Node<SPACE_DIM>* node_two = elem_iter->GetNode(2);

    	// Get direction of edges. Edge i links nodes i and i+1 mod 3
    	c_vector<double,SPACE_DIM> edge_zero = node_one->GetPoint().rGetLocation() - node_zero->GetPoint().rGetLocation();
    	edge_zero /= norm_2(edge_zero);
    	c_vector<double,SPACE_DIM> edge_one = node_two->GetPoint().rGetLocation() - node_one->GetPoint().rGetLocation();
    	edge_one /= norm_2(edge_one);
    	c_vector<double,SPACE_DIM> edge_two = node_zero->GetPoint().rGetLocation() - node_two->GetPoint().rGetLocation();
    	edge_two /= norm_2(edge_two);

    	// Calculate which edge is most parallel to tangent_vector
    	double t_dot_zero = fabs(inner_prod(tangent_vector, edge_zero));
    	double t_dot_one = fabs(inner_prod(tangent_vector, edge_one));
    	double t_dot_two = fabs(inner_prod(tangent_vector, edge_two));

    	unsigned most_parallel_edge;
    	if (t_dot_zero >= t_dot_one && t_dot_zero >= t_dot_two)
    	{
    		most_parallel_edge = 0;
    	}
    	else if (t_dot_one >= t_dot_zero && t_dot_one >= t_dot_two)
    	{
    		most_parallel_edge = 1;
    	}
    	else if (t_dot_two >= t_dot_zero && t_dot_two >= t_dot_one)
    	{
    		most_parallel_edge = 2;
    	}
    	else
    	{
    		NEVER_REACHED;
    	}

    	// Calculate the orthogonal vector to the tangent_vector which lies in the plane of the element

    	// Vertices a and b are on Fixed edge and vertex C moves to get correct scaling
    	Node<SPACE_DIM>* vertex_a = elem_iter->GetNode(most_parallel_edge);
    	Node<SPACE_DIM>* vertex_b = elem_iter->GetNode((most_parallel_edge+1)%3);
    	Node<SPACE_DIM>* vertex_c = elem_iter->GetNode((most_parallel_edge+2)%3);

    	c_vector<double,SPACE_DIM> vertex_a_location = vertex_a->GetPoint().rGetLocation();
    	c_vector<double,SPACE_DIM> vertex_b_location = vertex_b->GetPoint().rGetLocation();
    	c_vector<double,SPACE_DIM> vertex_c_location = vertex_c->GetPoint().rGetLocation();

    	// Get direction and magnitude of edges. Edge i links nodes i and i+1 mod 3
    	c_vector<double,SPACE_DIM> edge_a = vertex_b_location - vertex_a_location;
    	c_vector<double,SPACE_DIM> edge_b = vertex_c_location - vertex_b_location;
    	c_vector<double,SPACE_DIM> edge_c = vertex_a_location - vertex_c_location;

    	// Caclulate rest length scalings for the other two edges.
    	c_vector<double,SPACE_DIM> normal_in_plane = edge_b - tangent_vector*inner_prod(edge_b,tangent_vector);
    	normal_in_plane /= norm_2(normal_in_plane);

    	double height_b = inner_prod(edge_b,normal_in_plane);
    	double height_c = inner_prod(-edge_c,normal_in_plane);

    	if ( (height_c<=0) ||(height_b<=0))
    	{
    		    	// PRINT_3_VARIABLES(node_zero->GetPoint().rGetLocation(),node_one->GetPoint().rGetLocation(),node_two->GetPoint().rGetLocation());
    		    	// PRINT_3_VARIABLES(tangent_vector,normal_in_plane,element_normal);
    		    	PRINT_VARIABLE(most_parallel_edge);
    		    	// PRINT_3_VARIABLES(edge_zero,edge_one,edge_two);
    		    	// PRINT_3_VARIABLES(t_dot_zero,t_dot_one,t_dot_two)
    		    	// PRINT_3_VARIABLES(vertex_a_location,vertex_b_location,vertex_c_location);
//    		    	PRINT_4_VARIABLES(vertex_a_location,vertex_b_location,vertex_c_location,vertex_c_new_location);

    		    	// PRINT_3_VARIABLES(edge_a,edge_b,edge_c);
    		    	// PRINT_2_VARIABLES(height_b,height_c);
    		    	//PRINT_3_VARIABLES(edge_a_reduction_ratio,edge_b_reduction_ratio,edge_c_reduction_ratio)

    	}

    	assert(height_b>0);
    	assert(height_c>0);

    	double height = fmax(height_b,height_c);

    	c_vector<double,SPACE_DIM> vertex_c_new_location = vertex_c_location - reduction_factor*height*normal_in_plane;

    	double edge_a_reduction_ratio = 1.0;
    	double edge_b_reduction_ratio = norm_2(vertex_b_location - vertex_c_new_location)/norm_2(vertex_b_location - vertex_c_location);
    	double edge_c_reduction_ratio = norm_2(vertex_a_location - vertex_c_new_location)/norm_2(vertex_a_location - vertex_c_location);


    	// Reduce the rest length according to the above ratio
    	//TODO make this addative not compund

    	// Edge b
    	std::pair<unsigned, unsigned> node_pair = p_cell_population->CreateOrderedPair(vertex_b->GetIndex(),vertex_c->GetIndex());

		if (mDivisionsApplied.at(node_pair) < mMaxDivisions)
		{
			double previous_rest_length = p_cell_population->GetRestLength(node_pair.first, node_pair.second);
			double new_rest_length = (1 - reduction_factor) * previous_rest_length;
			p_cell_population->SetRestLength(node_pair.first, node_pair.second, new_rest_length);
			mDivisionsApplied[node_pair]++;
		}

		// Edge c
		node_pair = p_cell_population->CreateOrderedPair(vertex_a->GetIndex(),vertex_c->GetIndex());

		if (mDivisionsApplied.at(node_pair) < mMaxDivisions)
		{
			double previous_rest_length = p_cell_population->GetRestLength(node_pair.first, node_pair.second);
			double new_rest_length = (1 - reduction_factor) * previous_rest_length;
			p_cell_population->SetRestLength(node_pair.first, node_pair.second, new_rest_length);
			mDivisionsApplied[node_pair]++;
		}
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceOffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::UpdateCellData()
{
	assert(SPACE_DIM==3); // Currently assumes that SPACE_DIM = 3



    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
         cell_iter != this->mrCellPopulation.End();
         ++cell_iter)
    {
        c_vector<double, SPACE_DIM> location = this->mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

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
        double voronoi_cell_area = this->mrCellPopulation.GetVolumeOfCell(*cell_iter);

        c_vector<double,3> force = mAppliedTractions[nearest_fluid_site]*voronoi_cell_area;
        c_vector<double,3> shear_stress = mAppliedTangentTractions[nearest_fluid_site];

        assert(fabs(force[0])<1e10);
        assert(fabs(force[1])<1e10);
        assert(fabs(force[2])<1e10);
        assert(fabs(shear_stress[0])<1e10);
        assert(fabs(shear_stress[1])<1e10);
        assert(fabs(shear_stress[2])<1e10);

        // Store the force in CellData
        cell_iter->GetCellData()->SetItem("applied_force_x", force[0]);
        cell_iter->GetCellData()->SetItem("applied_force_y", force[1]);
        cell_iter->GetCellData()->SetItem("applied_force_z", force[2]);
        cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(force));
        cell_iter->GetCellData()->SetItem("voronoi_cell_area", voronoi_cell_area);
        cell_iter->GetCellData()->SetItem("applied_shear_stress_x", shear_stress[0]);
        cell_iter->GetCellData()->SetItem("applied_shear_stress_y", shear_stress[1]);
        cell_iter->GetCellData()->SetItem("applied_shear_stress_z", shear_stress[2]);
        cell_iter->GetCellData()->SetItem("applied_shear_stress_mag", norm_2(shear_stress));
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AppliedForceOffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::GetAppliedForce(c_vector<double, SPACE_DIM> location)
{
	c_vector<double, SPACE_DIM> force = zero_vector<double>(SPACE_DIM);

	///////// Pulse going along cylindrical pipe /////////////////////////////////
//	// Currently creates a simple force field based on position but will eventually load this data from file.
//    SimulationTime* p_simulation_time = SimulationTime::Instance();
//    double current_time = p_simulation_time->GetTime();
//
//	c_vector<double, SPACE_DIM> radial_location = location;
//	radial_location[SPACE_DIM-1] = 0;
//
//	double height = location[SPACE_DIM-1];
//	double radial_distance = norm_2(radial_location);
//
//	double magnitude = 2.0*(tanh(height+0.5-(current_time-7.5)) - tanh(height-0.5-(current_time-7.5)));
//
//	if (radial_distance > 1e-4)
//	{
//		force = magnitude*radial_location/radial_distance;
//	}
	//////////////////////////////////////////////////////////////////////////////

	return force;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AppliedForceOffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::GetResetTractionsOnCells()
{
    return mResetTractionsOnCells;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceOffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::SetResetTractionsOnCells(bool resetTractionsOnCells, std::string tractionFile)
{
	assert(!resetTractionsOnCells || tractionFile!="");
	mResetTractionsOnCells = resetTractionsOnCells;
	mTractionFile = tractionFile;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceOffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::LoadTractionFromFile()
{
	PRINT_VARIABLE(mTractionFile);
	FILE* traction_file = fopen((char*)mTractionFile.c_str(), "r");
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

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AppliedForceOffLatticeSimulation<1,1>;
template class AppliedForceOffLatticeSimulation<1,2>;
template class AppliedForceOffLatticeSimulation<2,2>;
template class AppliedForceOffLatticeSimulation<1,3>;
template class AppliedForceOffLatticeSimulation<2,3>;
template class AppliedForceOffLatticeSimulation<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AppliedForceOffLatticeSimulation)

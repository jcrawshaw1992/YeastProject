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

#include "SpringLengthModifier.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "Debug.hpp"
#include <math.h>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SpringLengthModifier<ELEMENT_DIM,SPACE_DIM>::SpringLengthModifier()
    : AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>(),
      mShearThreshold(DOUBLE_UNSET), // Needs to be set!
      mReductionFactor(DOUBLE_UNSET), // Needs to be set!
      mMaxDivisions(UNSIGNED_UNSET), // Needs to be set!
      mFirstRun(true)
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SpringLengthModifier<ELEMENT_DIM,SPACE_DIM>::~SpringLengthModifier()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SpringLengthModifier<ELEMENT_DIM,SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */

	// on the first run store the steady state rest lengths
	if(mFirstRun)
	{
		MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population = dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);
		for (typename MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator spring_iterator = p_cell_population->SpringsBegin();
			 spring_iterator != p_cell_population->SpringsEnd();
			 ++spring_iterator)
		{
			unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
			unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

			std::pair<unsigned, unsigned> node_pair =  rCellPopulation.CreateOrderedPair(nodeA_global_index, nodeB_global_index);

			mDivisionsApplied[node_pair] = 0;
			mOriginalRestLength[node_pair] = p_cell_population->GetRestLength(node_pair.first, node_pair.second);
			mTargetRestLengthProportion[node_pair] = 1.0;
		}
		mFirstRun = false;
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SpringLengthModifier<ELEMENT_DIM,SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // These member variables need to be set to work
    assert(mShearThreshold!=DOUBLE_UNSET);
    assert(mReductionFactor!=DOUBLE_UNSET);
    assert(mMaxDivisions!=UNSIGNED_UNSET);
    //assert(mReductionFactor*mMaxDivisions<1); // If not then get could get degenerate elements.


	assert(SPACE_DIM==3);

	MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population = dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(rCellPopulation));
	assert(p_cell_population != NULL);


    // Loop over elements and calculate the new rest lengths
    for (typename MutableMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
             elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
             ++elem_iter)
	{
    	CellPtr p_cell_0 = rCellPopulation.GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(0));
    	CellPtr p_cell_1 = rCellPopulation.GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(1));
    	CellPtr p_cell_2 = rCellPopulation.GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(2));

    	// Get the applied shear stress from cell data and only shrink if below thereshold!
    	c_vector<double, SPACE_DIM> tangent_vector;
    	tangent_vector[0] = ( p_cell_0->GetCellData()->GetItem("applied_shear_stress_x")
    	    		        + p_cell_1->GetCellData()->GetItem("applied_shear_stress_x")
    	    			    + p_cell_2->GetCellData()->GetItem("applied_shear_stress_x") ) /3.0;
    	tangent_vector[1] = ( p_cell_0->GetCellData()->GetItem("applied_shear_stress_y")
    						+ p_cell_1->GetCellData()->GetItem("applied_shear_stress_y")
    						+ p_cell_2->GetCellData()->GetItem("applied_shear_stress_y") ) /3.0;
    	tangent_vector[2] = ( p_cell_0->GetCellData()->GetItem("applied_shear_stress_z")
    						+ p_cell_1->GetCellData()->GetItem("applied_shear_stress_z")
    						+ p_cell_2->GetCellData()->GetItem("applied_shear_stress_z") ) /3.0;

    	double element_shear= norm_2(tangent_vector);
    	tangent_vector/=element_shear;

    	if (element_shear<mShearThreshold)
    	{
			// Make sure tangent_vector lies in plane. Remove component in normal direction
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

            // Calculate the orthogonal vector to the tangent_vector which lies in the plane of the element
			c_vector<double,SPACE_DIM> normal_in_plane = edge_b - tangent_vector*inner_prod(edge_b,tangent_vector);
			normal_in_plane /= norm_2(normal_in_plane);

//		    // Caclulate rest length scalings for the other two edges.
//			double n_component_of_b = inner_prod(edge_b,normal_in_plane);
//			double n_component_of_c = inner_prod(-edge_c,normal_in_plane);
//
//			double t_component_of_b = inner_prod(edge_b,tangent_vector);
//			double t_component_of_c = inner_prod(edge_c,tangent_vector);
//
//			PRINT_4_VARIABLES(n_component_of_b,n_component_of_c,t_component_of_b,t_component_of_c);
//
//
	//    	if ( (height_c<=0) ||(height_b<=0))
	//    	{
//	//    		    	PRINT_3_VARIABLES(node_zero->GetPoint().rGetLocation(),node_one->GetPoint().rGetLocation(),node_two->GetPoint().rGetLocation());
//	    		    	PRINT_3_VARIABLES(tangent_vector[0],tangent_vector[1],tangent_vector[2]);
//	    		    	PRINT_3_VARIABLES(normal_in_plane[0],normal_in_plane[1],normal_in_plane[2]);
//
//	    		    	//PRINT_3_VARIABLES(tangent_vector,normal_in_plane,element_normal);
//
//	//    		    	PRINT_VARIABLE(most_parallel_edge);
//	//    		    	PRINT_3_VARIABLES(edge_zero,edge_one,edge_two);
//	//    		    	PRINT_3_VARIABLES(t_dot_zero,t_dot_one,t_dot_two)
//	    	//	    	PRINT_3_VARIABLES(vertex_a_location,vertex_b_location,vertex_c_location);
//	    		    	PRINT_3_VARIABLES(vertex_a_location[0],vertex_a_location[1],vertex_a_location[2]);
//	    		    	PRINT_3_VARIABLES(vertex_b_location[0],vertex_b_location[1],vertex_b_location[2]);
//	    		    	PRINT_3_VARIABLES(vertex_c_location[0],vertex_c_location[1],vertex_c_location[2]);
	//    		    	//PRINT_4_VARIABLES(vertex_a_location,vertex_b_location,vertex_c_location,vertex_c_new_location);
	//    		    	PRINT_3_VARIABLES(edge_a,edge_b,edge_c);
	//    		    	PRINT_2_VARIABLES(height_b,height_c);
	//    		    	//PRINT_3_VARIABLES(edge_a_reduction_ratio,edge_b_reduction_ratio,edge_c_reduction_ratio)
	//    	}
//
//			assert(n_component_of_b>0);
//			assert(n_component_of_c>0);
//			assert(t_component_of_b*t_component_of_c>=0);
//
//			//double height = fmin(height_b,height_c);
//			double height = (t_component_of_c * n_component_of_b + t_component_of_b * n_component_of_c) /(t_component_of_b +t_component_of_c);
//
//			c_vector<double,SPACE_DIM> vertex_c_new_location = vertex_c_location - mReductionFactor*height*normal_in_plane;
//
//
			c_vector<double,SPACE_DIM> normal_to_edge_b = edge_b - inner_prod(edge_b,edge_a)/inner_prod(edge_a,edge_a)*edge_a;

			// assert(fabs(inner_prod(edge_b,normal_to_edge_b))<1e-10);
			assert(fabs(inner_prod(edge_b,normal_to_edge_b))<1.0);

	        double alpha = inner_prod(edge_c, normal_to_edge_b)/inner_prod(normal_in_plane,normal_to_edge_b);

	        c_vector<double,SPACE_DIM> vertex_c_new_location = vertex_c_location + mReductionFactor*alpha*normal_in_plane;

           // PRINT_3_VARIABLES(vertex_c_new_location[0],vertex_c_new_location[1],vertex_c_new_location[2]);

			double edge_a_reduction_ratio = 1.0;
			double edge_b_reduction_ratio = norm_2(vertex_b_location - vertex_c_new_location)/norm_2(vertex_b_location - vertex_c_location);
			double edge_c_reduction_ratio = norm_2(vertex_a_location - vertex_c_new_location)/norm_2(vertex_a_location - vertex_c_location);

			assert(edge_b_reduction_ratio <=1.0);
			assert(edge_c_reduction_ratio <=1.0);

			//Stop edges shrinking more than mReductionFactor each step
            if(edge_b_reduction_ratio < 1.0 - mReductionFactor)
            {
                edge_b_reduction_ratio = 1.0 - mReductionFactor;
            }
            if(edge_c_reduction_ratio <= 1.0 - mReductionFactor)
            {
                edge_c_reduction_ratio = 1.0 - mReductionFactor;
            }
            assert(edge_b_reduction_ratio >= 1.0 - mReductionFactor);
            assert(edge_c_reduction_ratio >= 1.0 - mReductionFactor);


			//PRINT_3_VARIABLES(edge_a_reduction_ratio,edge_b_reduction_ratio,edge_c_reduction_ratio);

			// Reduce the rest length according to the above ratios

            // Edge a
            std::pair<unsigned, unsigned> node_pair = p_cell_population->CreateOrderedPair(vertex_b->GetIndex(),vertex_c->GetIndex());
            mDivisionsApplied[node_pair]++;

			// Edge b
			node_pair = p_cell_population->CreateOrderedPair(vertex_b->GetIndex(),vertex_c->GetIndex());

			if (mDivisionsApplied.at(node_pair) < mMaxDivisions)
			{
				double previous_rest_length = p_cell_population->GetRestLength(node_pair.first, node_pair.second);
				double new_rest_length = previous_rest_length - (1.0-edge_b_reduction_ratio) * mOriginalRestLength.at(node_pair);
				assert(new_rest_length>0);
				p_cell_population->SetRestLength(node_pair.first, node_pair.second, new_rest_length);
				mDivisionsApplied[node_pair]++;
			}

			// Edge c
			node_pair = p_cell_population->CreateOrderedPair(vertex_a->GetIndex(),vertex_c->GetIndex());

			if (mDivisionsApplied.at(node_pair) < mMaxDivisions)
			{
				double previous_rest_length = p_cell_population->GetRestLength(node_pair.first, node_pair.second);
				double new_rest_length = previous_rest_length - (1.0-edge_c_reduction_ratio) * mOriginalRestLength.at(node_pair);
				assert(new_rest_length>0);
				p_cell_population->SetRestLength(node_pair.first, node_pair.second, new_rest_length);
				mDivisionsApplied[node_pair]++;
			}
    	}
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SpringLengthModifier<ELEMENT_DIM,SPACE_DIM>::GetShearThreshold()
{
    return mShearThreshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SpringLengthModifier<ELEMENT_DIM,SPACE_DIM>::SetShearThreshold(double shearThreshold)
{
	mShearThreshold = shearThreshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SpringLengthModifier<ELEMENT_DIM,SPACE_DIM>::GetReductionFactor()
{
    return mReductionFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SpringLengthModifier<ELEMENT_DIM,SPACE_DIM>::SetReductionFactor(double reductionFactor)
{
    mReductionFactor = reductionFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SpringLengthModifier<ELEMENT_DIM,SPACE_DIM>::GetMaxDivisions()
{
    return mMaxDivisions;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SpringLengthModifier<ELEMENT_DIM,SPACE_DIM>::SetMaxDivisions(double maxDivisions)
{
	mMaxDivisions = maxDivisions;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SpringLengthModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class SpringLengthModifier<1,1>;
template class SpringLengthModifier<1,2>;
template class SpringLengthModifier<2,2>;
template class SpringLengthModifier<1,3>;
template class SpringLengthModifier<2,3>;
template class SpringLengthModifier<3,3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SpringLengthModifier)


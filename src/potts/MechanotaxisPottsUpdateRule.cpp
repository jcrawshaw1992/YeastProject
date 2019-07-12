

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

// template<unsigned DIM>
// double MechanotaxisPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
//                                                                          unsigned targetNodeIndex,
//                                                                          PottsBasedCellPopulation<DIM>& rCellPopulation)
// {
//     assert(DIM==3);
//     double delta_H = 0;

//   //  TRACE("Getting Hamiltonian for shear stress");
//     // assert(!rCellPopulation.GetNode(targetNodeIndex)->IsBoundaryNode());
//     // assert(!rCellPopulation.GetNode(currentNodeIndex)->IsBoundaryNode()); // this causes it to break if one or the other is a neighbour, which is stupid 

//     PottsArbitrarySurfaceIn3DMesh<DIM>* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*> (&(rCellPopulation.rGetMesh())); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>
    
//     c_vector<double, DIM> displacement_direction = p_static_cast_potts_mesh->GetNode(targetNodeIndex)->rGetLocation() - p_static_cast_potts_mesh->GetNode(currentNodeIndex)->rGetLocation();
//     displacement_direction /= norm_2(displacement_direction);

//     std::set<unsigned> containing_elements = p_static_cast_potts_mesh->GetNode(currentNodeIndex)->rGetContainingElementIndices();
//     std::set<unsigned> new_location_containing_elements = p_static_cast_potts_mesh->GetNode(targetNodeIndex)->rGetContainingElementIndices();

//     bool current_node_contained = !containing_elements.empty();
//     bool target_node_contained = !new_location_containing_elements.empty();

//      if (!current_node_contained && !target_node_contained)
//     {
//         EXCEPTION("At least one of the current node or target node must be in an element.");
//     }

//     if (current_node_contained && target_node_contained)
//     {
//         if (*(new_location_containing_elements.begin()) == *(containing_elements.begin()))
//         {
//             EXCEPTION("The current node and target node must not be in the same element.");
//         }
//     }

    
//     // if (current_node_contained) // current node is in an element
//     // {
//     // unsigned current_element = (*containing_elements.begin());

//     // double TheApsectRatio = p_static_cast_potts_mesh->GetAspectRatio(current_element);
//     // }

//     // PRINT_VARIABLE(TheApsectRatio );


  

//     // // Note that we define these vectors before setting them as otherwise the profiling build will break (see #2367)
//     // c_vector<double, DIM> current_location;
//     // current_location = p_static_cast_potts_mesh->GetNode(currentNodeIndex)->rGetLocation();
//     // c_vector<double, DIM> target_location;
//     // target_location = p_static_cast_potts_mesh->GetNode(targetNodeIndex)->rGetLocation();

   
//     //     if (target_location[2] > current_location[2])
//     //     {
//     //         delta_H -= 2;
//     //     }
//     //     else if (target_location[2] < current_location[2])
//     //     {
//     //         delta_H += 2;
//     //     }


//     // return delta_H;



// // //  // Just silenced this to see if I can get any migration at all 




//     //  if (current_node_contained) // current node is in an element
//     // {
//     //     unsigned current_element = (*containing_elements.begin());
//     //     PottsElement<DIM>* pCurrentElement = p_static_cast_potts_mesh->GetElement(current_element);


//     //     c_vector<double, DIM> CurrentShearStress = p_static_cast_potts_mesh->GetTractionOnElement(current_element);

//     //     double TheApsectRatio = p_static_cast_potts_mesh->GetAspectRatio();
   

//     //     // Mock change of spin to call GetSurfaceAreaOfElement with new configuration
//     //     Node<DIM>* pTargetNode = p_static_cast_potts_mesh->GetNode(targetNodeIndex);
//     //     pCurrentElement->AddNode(pTargetNode);

//     //     // Calculate the energy change from the source cell
//     //     if (target_node_contained) 
//     //     {
//     //       // if the target node is already in a cell, it needs to be removed from the cell and put it in the other cell 
//     //         unsigned target_element = (*new_location_containing_elements.begin());
//     //         PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);
//     //         pTargetElement->DeleteNode(pTargetElement->GetNodeLocalIndex(pTargetNode->GetIndex()));
//     //     }

//     //     c_vector<double, DIM> ShearStressAfterSwitch  = p_static_cast_potts_mesh->GetTractionOnElement(current_element);
        
        
//     //     // Undo mocked change
//     //     pCurrentElement->DeleteNode(pCurrentElement->GetNodeLocalIndex(pTargetNode->GetIndex()));
//     //     if (target_node_contained)
//     //     {
//     //         unsigned target_element = (*new_location_containing_elements.begin());
//     //         PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);
//     //         pTargetElement->AddNode(pTargetNode);
//     //     }

//     //      c_vector<double, DIM> ShearStressGrad = ShearStressAfterSwitch-CurrentShearStress ;
         
//     //      PRINT_VARIABLE(norm_2(ShearStressGrad)); //isnan(0.0)

//     //      delta_H +=  mTractionCorrelationParameter * norm_2(ShearStressGrad); //  0.5 * (1 + inner_prod(displacement_direction, ShearStressGrad));
//     //    }




//     // // Calculate the energy change from the target cell
//     // if (target_node_contained) // target node is in an element
//     // {
//     //     unsigned target_element = (*new_location_containing_elements.begin());
//     //     PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);

//     //     c_vector<double, DIM> TargetShearStress = p_static_cast_potts_mesh->GetTractionOnElement(target_element);
      
//     //     Node<DIM>* pTargetNode = p_static_cast_potts_mesh->GetNode(targetNodeIndex);
//     //     pTargetElement->DeleteNode(pTargetElement->GetNodeLocalIndex(pTargetNode->GetIndex()));

//     //     c_vector<double, DIM> TargetShearStressAfterSwitch = p_static_cast_potts_mesh->GetTractionOnElement(target_element);
       
//     //     // GetSurfaceAreaOfElement(target_element, pTargetElement); 

//     //     c_vector<double, DIM> TargetShearStressGradient = TargetShearStressAfterSwitch -  TargetShearStress;
//     //      PRINT_VECTOR(TargetShearStress);
//     //       PRINT_VECTOR(TargetShearStressAfterSwitch);
//     //     PRINT_VECTOR(TargetShearStressGradient);
  
//     //     // Add node back after calculating aspect ratio
//     //     pTargetElement->AddNode(pTargetNode);
//     //     PRINT_VARIABLE(norm_2(TargetShearStressGradient));
        
        
//     //     delta_H +=  mTractionCorrelationParameter * norm_2(TargetShearStressGradient); // * 0.5 * (1 + inner_prod(displacement_direction, TargetShearStressGradient));
//     //   }




//     return delta_H;
// }




    // several cases, if they are in the same element, then dont bother, if they are in different elements then do? 
    // if (current_node_contained) // current node is in an element
    // {
    //    unsigned current_element = (*containing_elements.begin());
    //   PRINT_VARIABLE(current_element);

    //   TRACE("get here");
    //   c_vector<double,DIM> ShearStress = p_static_cast_potts_mesh->GetTractionOnElement(current_element);
    //   TRACE("B");
    //   double ShearStressMagnitiude = norm_2(ShearStress);
    //   ShearStress /= ShearStressMagnitiude ;
    //   PRINT_VECTOR(ShearStress);
      

    //   // Favour moves along the shear stress vector (3rd term ranges [0,2] hence the 0.5 factor)
    //   delta_H = - mTractionCorrelationParameter * 0.5 * (1 + inner_prod(displacement_direction, ShearStress));
    // //   }
       
    //    c_vector<double,DIM> ShearStress = p_static_cast_potts_mesh->GetTractionOnLattice(currentNodeIndex);
    //    c_vector<double,DIM> ShearStressTarget = p_static_cast_potts_mesh->GetTractionOnLattice(targetNodeIndex);
    //     c_vector<double,DIM> ShearStressGrad = ShearStress-ShearStressTarget;
    //     double absGrad = norm_2(ShearStressGrad);
    //     PRINT_VARIABLE(absGrad );

    // delta_H = - mTractionCorrelationParameter * 0.5 * (1 + inner_prod(displacement_direction, ShearStressGrad));
        

    





// template<unsigned DIM>
// MechanotaxisPottsUpdateRule<DIM>::~MechanotaxisPottsUpdateRule()
// {
// }

template<unsigned DIM>
double MechanotaxisPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                         unsigned targetNodeIndex,
                                                                         PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    assert(DIM==3);

    if (!rCellPopulation.GetNode(targetNodeIndex)->IsBoundaryNode() )
    {
        return 0;
    }
    else if (!rCellPopulation.GetNode(currentNodeIndex)->IsBoundaryNode() )
    {
        return 0;
    }
    // assert(!rCellPopulation.GetNode(targetNodeIndex)->IsBoundaryNode());
    // assert(!rCellPopulation.GetNode(currentNodeIndex)->IsBoundaryNode());

    PottsArbitrarySurfaceIn3DMesh<DIM>* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*> (&(rCellPopulation.rGetMesh())); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>
    
    c_vector<double, DIM> displacement_direction = p_static_cast_potts_mesh->GetNode(targetNodeIndex)->rGetLocation() - p_static_cast_potts_mesh->GetNode(currentNodeIndex)->rGetLocation();
    displacement_direction /= norm_2(displacement_direction);


    c_vector<double,DIM> traction_direction =  p_static_cast_potts_mesh->GetTractionOnLattice(currentNodeIndex);
     // traction_direction /= norm_2(traction_direction);

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

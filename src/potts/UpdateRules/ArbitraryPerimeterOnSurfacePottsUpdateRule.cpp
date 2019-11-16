

#include "ArbitraryPerimeterOnSurfacePottsUpdateRule.hpp"
#include "CPMGTPaseEventHandler.hpp"
#include "Debug.hpp"

template <unsigned DIM>
ArbitraryPerimeterOnSurfacePottsUpdateRule<DIM>::ArbitraryPerimeterOnSurfacePottsUpdateRule()
        : AbstractWrappedPottsUpdateRule<DIM>(),
          mSurfaceAreaEnergyParameter(0.5), // Educated guess
          mTargetSurfaceArea(16.0) // Defaults to a 4*4 cell size

{
    /// \todo Default values don't apply in 3D.
}

template <unsigned DIM>
ArbitraryPerimeterOnSurfacePottsUpdateRule<DIM>::~ArbitraryPerimeterOnSurfacePottsUpdateRule()
{
}
//********



// 



template <unsigned DIM>
double ArbitraryPerimeterOnSurfacePottsUpdateRule<DIM>::DetermineAReasonableTargetPerimeter(double x, double y)
{
   // 1) Decide how many lattice sites I want in my cell, which is x wide and y long. Im thinking 5*10 (5 accross and 10 up)

    unsigned lattice_site_index = 22; 
    std::vector<std::pair<unsigned, unsigned> > ElementPairsAssociatedWithLatticeSite = mMapNodesToAssociateElementPairs[lattice_site_index];
    std::vector<std::pair<unsigned, unsigned> > ElementPairsAssociatedWithPottsCell;
   
    std::vector<unsigned> EdgeLenghts ;
    for (int j = 0; j < ElementPairsAssociatedWithLatticeSite.size(); j++)
    {
        EdgeLenghts.push_back(mDistanceBetweenElements[ElementPairsAssociatedWithLatticeSite[j] ]);
        
    }
    PRINT_VARIABLE(ElementPairsAssociatedWithLatticeSite.size());
    PRINT_VECTOR(EdgeLenghts); // These are giving zeros, dont know why 


    unsigned up = *std::min_element(EdgeLenghts.begin(), EdgeLenghts.end());
    unsigned diagonal = *std::max_element(EdgeLenghts.begin(), EdgeLenghts.end());

    // This is calculated based on a rectangular Potts cell 
    double UpwardsLength = y * up + y * 2 * diagonal;
    double HorizontalLenght = x * 2 * diagonal;

    double TargetPerimeter = 2*(UpwardsLength  + HorizontalLenght  );
    PRINT_VARIABLE(TargetPerimeter);
    return TargetPerimeter;


}




//---------------------

template <unsigned DIM>
double ArbitraryPerimeterOnSurfacePottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                                           unsigned targetNodeIndex,
                                                                                           WrappedPottsBasedCellPopulation<DIM>& rCellPopulation)
{
    CPMGTPaseEventHandler::BeginEvent(CPMGTPaseEventHandler::SURFACE);
  
//   TRACE("Evaluate Surface area Hamiltonian");
   PottsArbitrarySurfaceIn3DMesh<DIM>* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*> (&(rCellPopulation.rGetMesh())); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>
 
 
// note the other way to do this is do the static cast each time you need it 
//static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&rCellPopulation.rGetMesh())->GetVolumeOfLatticeSite(targetNodeIndex);


    double delta_H = 0.0;
    std::set<unsigned> containing_elements = p_static_cast_potts_mesh->GetNode(currentNodeIndex)->rGetContainingElementIndices();
    // TRACE("B");
    std::set<unsigned> new_location_containing_elements = p_static_cast_potts_mesh->GetNode(targetNodeIndex)->rGetContainingElementIndices();
    // TRACE("C");
    bool current_node_contained = !containing_elements.empty();
    bool target_node_contained = !new_location_containing_elements.empty();

    // Every node must each be in at most one element
    assert(new_location_containing_elements.size() < 2);

    if (!current_node_contained && !target_node_contained)
    {
        EXCEPTION("At least one of the current node or target node must be in an element.");
    }

    if (current_node_contained && target_node_contained)
    {
        if (*(new_location_containing_elements.begin()) == *(containing_elements.begin()))
        {
            EXCEPTION("The current node and target node must not be in the same element.");
        }
    }



// IF the target element is also the sister WHEN the sister only has one node, the code will break because you need to do a trial switch, which will leave the sister empty!!!!
      

    //double surface_area_target_lattice_site = mLatticePerimeters[targetNodeIndex] ;//static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&rCellPopulation.rGetMesh())->GetSurfaceAreaOfLatticeSite(targetNodeIndex);
    
    if (current_node_contained) // current node is in an element
    {
        unsigned current_element = (*containing_elements.begin());
        PottsElement<DIM>* pCurrentElement = p_static_cast_potts_mesh->GetElement(current_element);

        double current_surface_area = 0;
  
        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(current_element);

          if(rCellPopulation.IsPottsSimulationPeriodic())
            {
                double Center = p_cell->GetCellData()->GetItem("Center");
                unsigned Sister_Element = rCellPopulation.GetSister(current_element);
                current_surface_area = p_static_cast_potts_mesh->GetPerimeterOfCoupledElements(current_element, Sister_Element, Center );
            }else
            {
                current_surface_area = p_static_cast_potts_mesh->GetPerimeterOfElement(current_element);
            }
        
        double current_surface_area_difference = current_surface_area - mTargetSurfaceArea;
      
      
        // Mock change of spin to call GetPerimeterOfElement with new configuration
        Node<DIM>* pTargetNode = p_static_cast_potts_mesh->GetNode(targetNodeIndex);
        pCurrentElement->AddNode(pTargetNode);
        if (target_node_contained)
        {
       
            unsigned target_element = (*new_location_containing_elements.begin());
            PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);
            pTargetElement->DeleteNode(pTargetElement->GetNodeLocalIndex(pTargetNode->GetIndex()));

        }


        double current_surface_area_after_switch;
        if(rCellPopulation.IsPottsSimulationPeriodic())
            {
                double Center = p_cell->GetCellData()->GetItem("Center");
                unsigned Sister_Element = rCellPopulation.GetSister(current_element);
                current_surface_area_after_switch = p_static_cast_potts_mesh->GetPerimeterOfCoupledElements(current_element, Sister_Element, Center );
            }else
            {
                current_surface_area_after_switch = p_static_cast_potts_mesh->GetPerimeterOfElement(current_element);
            }
        double current_surface_area_difference_after_switch = current_surface_area_after_switch - mTargetSurfaceArea;
   
        // Undo mocked change
        pCurrentElement->DeleteNode(pCurrentElement->GetNodeLocalIndex(pTargetNode->GetIndex()));
        if (target_node_contained)
        {
            unsigned target_element = (*new_location_containing_elements.begin());
            PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);
            pTargetElement->AddNode(pTargetNode);
        }

        delta_H += mSurfaceAreaEnergyParameter * (current_surface_area_difference_after_switch * current_surface_area_difference_after_switch - current_surface_area_difference * current_surface_area_difference);
    }


    
    if (target_node_contained) // target node is in an element
    {
        unsigned target_element = (*new_location_containing_elements.begin());
        PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);


        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(target_element);
          double target_surface_area;
         if(rCellPopulation.IsPottsSimulationPeriodic())
            {
                double Center = p_cell->GetCellData()->GetItem("Center");
                unsigned Sister_Element = rCellPopulation.GetSister(target_element);
                target_surface_area = p_static_cast_potts_mesh->GetPerimeterOfCoupledElements(target_element, Sister_Element, Center );
            }else
            {
                target_surface_area = p_static_cast_potts_mesh->GetPerimeterOfElement(target_element);
            }
        double target_surface_area_difference = target_surface_area - mTargetSurfaceArea;

        Node<DIM>* pTargetNode = p_static_cast_potts_mesh->GetNode(targetNodeIndex);
        pTargetElement->DeleteNode(pTargetElement->GetNodeLocalIndex(pTargetNode->GetIndex()));

        // Make switch 
        double target_surface_area_after_switch;
         if(rCellPopulation.IsPottsSimulationPeriodic())
            {
                double Center = p_cell->GetCellData()->GetItem("Center");
                unsigned Sister_Element = rCellPopulation.GetSister(target_element);
                target_surface_area_after_switch = p_static_cast_potts_mesh->GetPerimeterOfCoupledElements(target_element, Sister_Element, Center );
            }else
            {
                target_surface_area_after_switch = p_static_cast_potts_mesh->GetPerimeterOfElement(target_element);
            }
        
        double target_surface_area_difference_after_switch = target_surface_area_after_switch - mTargetSurfaceArea;

        // Add node back after calculating aspect ratio
        pTargetElement->AddNode(pTargetNode);

        delta_H += mSurfaceAreaEnergyParameter * (target_surface_area_difference_after_switch * target_surface_area_difference_after_switch - target_surface_area_difference * target_surface_area_difference);
    }

    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::SURFACE);

    return delta_H;

}

template <unsigned DIM>
double ArbitraryPerimeterOnSurfacePottsUpdateRule<DIM>::GetSurfaceAreaEnergyParameter()
{
    return mSurfaceAreaEnergyParameter;
}

template <unsigned DIM>
void ArbitraryPerimeterOnSurfacePottsUpdateRule<DIM>::SetSurfaceAreaEnergyParameter(double surfaceAreaEnergyParameter)
{
    mSurfaceAreaEnergyParameter = surfaceAreaEnergyParameter;
}


template <unsigned DIM>
double ArbitraryPerimeterOnSurfacePottsUpdateRule<DIM>::GetTargetSurfaceArea() const
{
    return mTargetSurfaceArea;
}

template <unsigned DIM>
void ArbitraryPerimeterOnSurfacePottsUpdateRule<DIM>::SetTargetSurfaceArea(double targetSurfaceArea)
{
    assert(targetSurfaceArea >= 0.0);
    mTargetSurfaceArea = targetSurfaceArea;
}

   
template <unsigned DIM>
void ArbitraryPerimeterOnSurfacePottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SurfaceAreaEnergyParameter >" << mSurfaceAreaEnergyParameter << "</SurfaceAreaEnergyParameter >\n";
    *rParamsFile << "\t\t\t<TargetSurfaceArea>" << mTargetSurfaceArea << "</TargetSurfaceArea>\n";

    // Call method on direct parent class
    AbstractWrappedPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class ArbitraryVolumeConstraintPottsUpdateRule<1>;
template class ArbitraryPerimeterOnSurfacePottsUpdateRule<2>;
template class ArbitraryPerimeterOnSurfacePottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(ArbitraryVolumeConstraintPottsUpdateRule)
EXPORT_TEMPLATE_CLASS1(ArbitraryPerimeterOnSurfacePottsUpdateRule, 2)
EXPORT_TEMPLATE_CLASS1(ArbitraryPerimeterOnSurfacePottsUpdateRule, 3)

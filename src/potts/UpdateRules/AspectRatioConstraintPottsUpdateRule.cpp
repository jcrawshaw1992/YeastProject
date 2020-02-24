
#include "AspectRatioConstraintPottsUpdateRule.hpp"
#include "Debug.hpp"

template <unsigned DIM>
AspectRatioConstraintPottsUpdateRule<DIM>::AspectRatioConstraintPottsUpdateRule()
        : AbstractWrappedPottsUpdateRule<DIM>(),
          mAspectRatioEnergyParameter(0.5), // @todo Made up defaults
          mTargetAspectRatio(40.0)
{
}

template <unsigned DIM>
AspectRatioConstraintPottsUpdateRule<DIM>::~AspectRatioConstraintPottsUpdateRule()
{
}

template <unsigned DIM>
double AspectRatioConstraintPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                                  unsigned targetNodeIndex,
                                                                                  WrappedPottsBasedCellPopulation<DIM>& rCellPopulation)
{
    assert(DIM == 2 || DIM == 3);
    // DIM == 2;
    PottsArbitrarySurfaceIn3DMesh<DIM>* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&(rCellPopulation.rGetMesh())); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>

    double delta_H = 0.0;
    double delta_E = 0.0;
    double delta_O = 0.0;

    std::set<unsigned> containing_elements = rCellPopulation.GetNode(currentNodeIndex)->rGetContainingElementIndices();
    std::set<unsigned> new_location_containing_elements = rCellPopulation.GetNode(targetNodeIndex)->rGetContainingElementIndices();
    bool current_node_contained = !containing_elements.empty();
    bool target_node_contained = !new_location_containing_elements.empty();

    // Every node must each be in at most one element
    assert(new_location_containing_elements.size() < 2);
    if (!current_node_contained && !target_node_contained)
    {
        TRACE("!current_node_contained && !target_node_contained");
        EXCEPTION("At least one of the current node or target node must be in an element.");
    }
    if (current_node_contained && target_node_contained)
    {
        TRACE("current_node_contained && target_node_contained");
        if (*(new_location_containing_elements.begin()) == *(containing_elements.begin()))
        {
            EXCEPTION("The current node and target node must not be in the same element.");
        }
    }

    if (current_node_contained) // current node is in an element
    {
        unsigned current_element = (*containing_elements.begin());
        PottsElement<DIM>* pCurrentElement = rCellPopulation.rGetMesh().GetElement(current_element);
        double Aspect_Ratio_i0;
        c_vector<double, 2> ShearStress_i0;
        c_vector<double, 2> MajorAxis_i0;

        if (rCellPopulation.IsPottsSimulationPeriodic())
        {
            unsigned Sister_Element = rCellPopulation.GetSister(current_element);
            PottsElement<DIM>* pSisterElement = p_static_cast_potts_mesh->GetElement(Sister_Element);

            Aspect_Ratio_i0 = p_static_cast_potts_mesh->GetAspectRatio(current_element, Sister_Element);      
            // Get shear stress alignment with major axis 
            MajorAxis_i0 = p_static_cast_potts_mesh->GetMajorAxisVector(current_element);
            c_vector<double, 3> ShearStress = (p_static_cast_potts_mesh->GetTractionOnElement(current_element) + p_static_cast_potts_mesh->GetTractionOnElement(Sister_Element))/(pCurrentElement->GetNumNodes()+ pSisterElement->GetNumNodes()); 
            ShearStress_i0 = Create_c_vector(ShearStress[0],ShearStress[1]);
            // ShearStressProjection_0 = norm_2(Average_ShearStress) * cos(inner_prod(Average_ShearStress, MajorAxis)/ norm_2(Average_ShearStress) );
        }
        else
        {
            Aspect_Ratio_i0 = p_static_cast_potts_mesh->GetAspectRatio(current_element);
           
            // Get shear stress alignment with major axis 
            MajorAxis_i0 = p_static_cast_potts_mesh->GetMajorAxisVector(current_element);
            c_vector<double, 3> ShearStress = p_static_cast_potts_mesh->GetTractionOnElement(current_element)/pCurrentElement->GetNumNodes();      
            ShearStress_i0 = Create_c_vector(ShearStress[0],ShearStress[1]);

        }
        double Aspect_Ratio_Diff_i0 = Aspect_Ratio_i0 - mTargetAspectRatio;

        // PRINT_3_VARIABLES(current_aspect_ratio, Aspect_Ratio_Diff_i0, mTargetAspectRatio);
        // AverageShearStress= norm_2(p_static_cast_potts_mesh->GetTractionOnElement(current_element) + p_static_cast_potts_mesh->GetTractionOnElement(Sister_Element))/(pCurrentElement->GetNumNodes()+ pSisterElement->GetNumNodes()); 
           
        // Add node to element here
        Node<DIM>* pTargetNode = rCellPopulation.rGetMesh().GetNode(targetNodeIndex);
        pCurrentElement->AddNode(pTargetNode);

        double Aspect_Ratio_i1 ;
        c_vector<double, 2> ShearStress_i1;
        c_vector<double, 2> MajorAxis_i1;

        if (rCellPopulation.IsPottsSimulationPeriodic())
        {   
            unsigned Sister_Element = rCellPopulation.GetSister(current_element);
            PottsElement<DIM>* pSisterElement = p_static_cast_potts_mesh->GetElement(Sister_Element);

            Aspect_Ratio_i1 = p_static_cast_potts_mesh->GetAspectRatio(current_element, Sister_Element);      
            // Get shear stress alignment with major axis 
            MajorAxis_i1 = p_static_cast_potts_mesh->GetMajorAxisVector(current_element);
            c_vector<double, 3> ShearStress = (p_static_cast_potts_mesh->GetTractionOnElement(current_element) + p_static_cast_potts_mesh->GetTractionOnElement(Sister_Element))/(pCurrentElement->GetNumNodes()+ pSisterElement->GetNumNodes()); 
            ShearStress_i1 = Create_c_vector(ShearStress[0],ShearStress[1]);
            // ShearStressProjection_0 = norm_2(Average_ShearStress) * cos(inner_prod(Average_ShearStress, MajorAxis)/ norm_2(Average_ShearStress) );
        }
        else
        {
            Aspect_Ratio_i1 = p_static_cast_potts_mesh->GetAspectRatio(current_element);
           
            // Get shear stress alignment with major axis 
            MajorAxis_i1 = p_static_cast_potts_mesh->GetMajorAxisVector(current_element);
            c_vector<double, 3> ShearStress = p_static_cast_potts_mesh->GetTractionOnElement(current_element)/pCurrentElement->GetNumNodes();      
            ShearStress_i1 = Create_c_vector(ShearStress[0],ShearStress[1]);
        }

        double AverageShearStress = norm_2(ShearStress_i1 -  ShearStress_i0)/2;
        double Aspect_Ratio_Diff_i1 = Aspect_Ratio_i1 - mTargetAspectRatio;

        //  PRINT_3_VARIABLES(Aspect_Ratio_Diff_i1, Aspect_Ratio_Diff_i0_after_switch, mTargetAspectRatio);

        // Remove node after calculating aspect ratio
        pCurrentElement->DeleteNode(pCurrentElement->GetNodeLocalIndex(pTargetNode->GetIndex()));

        // delta_H += mAspectRatioEnergyParameter * AverageShearStress*(Aspect_Ratio_Diff_i0_after_switch * Aspect_Ratio_Diff_i0_after_switch - Aspect_Ratio_Diff_i0 * Aspect_Ratio_Diff_i0);
        // delta_H += mAspectRatioEnergyParameter *( ShearStressProjection_1 *  Aspect_Ratio_Diff_i1 * Aspect_Ratio_Diff_i1   -ShearStressProjection_0 * Aspect_Ratio_Diff_i0 * Aspect_Ratio_Diff_i0);
       
        delta_E += mAspectRatioEnergyParameter *( Aspect_Ratio_Diff_i1 * Aspect_Ratio_Diff_i1 - Aspect_Ratio_Diff_i0 * Aspect_Ratio_Diff_i0);
        delta_O += mOrientationParameter *AverageShearStress *( acos(inner_prod(ShearStress_i1, MajorAxis_i1)/ norm_2(ShearStress_i1)) - acos(inner_prod(ShearStress_i0, MajorAxis_i0)/ norm_2(ShearStress_i0)));
        
        // PRINT_3_VARIABLES(ShearStressProjection_0, ShearStressProjection_1, delta_H)
       }


        // 
       if (target_node_contained) // target node is in an element
       {

        unsigned target_element = (*new_location_containing_elements.begin());
        PottsElement<DIM>* pTargetElement = rCellPopulation.rGetMesh().GetElement(target_element);
        double Aspect_Ratio_j0;
        c_vector<double, 2> ShearStress_j0;
        c_vector<double, 2> MajorAxis_j0;

        if (rCellPopulation.IsPottsSimulationPeriodic())
        {
            unsigned Sister_Element = rCellPopulation.GetSister(target_element);
            PottsElement<DIM>* pSisterElement = p_static_cast_potts_mesh->GetElement(Sister_Element);

            Aspect_Ratio_j0 = p_static_cast_potts_mesh->GetAspectRatio(target_element, Sister_Element);
            // Get shear stress alignment with major axis
            MajorAxis_j0 = p_static_cast_potts_mesh->GetMajorAxisVector(target_element);
            c_vector<double, 3> ShearStress = (p_static_cast_potts_mesh->GetTractionOnElement(target_element) + p_static_cast_potts_mesh->GetTractionOnElement(Sister_Element))/(pTargetElement->GetNumNodes()+ pSisterElement->GetNumNodes());
            ShearStress_j0 = Create_c_vector(ShearStress[0],ShearStress[1]);
        }
        else
        {
        Aspect_Ratio_j0 = p_static_cast_potts_mesh->GetAspectRatio(target_element);
        // Get shear stress alignment with major axis
        MajorAxis_j0 = p_static_cast_potts_mesh->GetMajorAxisVector(target_element);
        c_vector<double, 3> ShearStress = p_static_cast_potts_mesh->GetTractionOnElement(target_element)/pTargetElement->GetNumNodes();
        ShearStress_j0 = Create_c_vector(ShearStress[0],ShearStress[1]);

        }
        double Aspect_Ratio_Diff_j0 = Aspect_Ratio_j0 - mTargetAspectRatio;

        Node<DIM>* pTargetNode = rCellPopulation.rGetMesh().GetNode(targetNodeIndex);
        pTargetElement->DeleteNode(pTargetElement->GetNodeLocalIndex(pTargetNode->GetIndex()));

        double Aspect_Ratio_j1;
        c_vector<double, 2> ShearStress_j1;
        c_vector<double, 2> MajorAxis_j1;

        if (rCellPopulation.IsPottsSimulationPeriodic())
        {
        unsigned Sister_Element = rCellPopulation.GetSister(target_element);
        PottsElement<DIM>* pSisterElement = p_static_cast_potts_mesh->GetElement(Sister_Element);

        Aspect_Ratio_j1 = p_static_cast_potts_mesh->GetAspectRatio(target_element, Sister_Element);
        // Get shear stress alignment with major axis
        MajorAxis_j1 = p_static_cast_potts_mesh->GetMajorAxisVector(target_element);
        c_vector<double, 3> ShearStress = (p_static_cast_potts_mesh->GetTractionOnElement(target_element) + p_static_cast_potts_mesh->GetTractionOnElement(Sister_Element))/(pTargetElement->GetNumNodes()+ pSisterElement->GetNumNodes());
        ShearStress_j1 = Create_c_vector(ShearStress[0],ShearStress[1]);
        // ShearStressProjection_0 = norm_2(Average_ShearStress) * cos(inner_prod(Average_ShearStress, MajorAxis)/ norm_2(Average_ShearStress) );
        }
        else
        {
        Aspect_Ratio_j1 = p_static_cast_potts_mesh->GetAspectRatio(target_element);
        // Get shear stress alignment with major axis
        MajorAxis_j1 = p_static_cast_potts_mesh->GetMajorAxisVector(target_element);
        c_vector<double, 3> ShearStress = p_static_cast_potts_mesh->GetTractionOnElement(target_element)/pTargetElement->GetNumNodes();
        ShearStress_j1 = Create_c_vector(ShearStress[0],ShearStress[1]);
        }

        double AverageShearStress = norm_2(ShearStress_j1 - ShearStress_j0)/2;
        double Aspect_Ratio_Diff_j1 = Aspect_Ratio_j1 - mTargetAspectRatio;

        // PRINT_3_VARIABLES(Aspect_Ratio_Diff_j1, Aspect_Ratio_Diff_j0_after_switch, mTargetAspectRatio);

        // Remove node after calculating aspect ratio
        pTargetElement->AddNode(pTargetNode);

        // delta_H += mAspectRatioEnergyParameter * AverageShearStress*(Aspect_Ratio_Diff_j0_after_switch * Aspect_Ratio_Diff_j0_after_switch - Aspect_Ratio_Diff_j0 * Aspect_Ratio_Diff_j0);
        // delta_H += mAspectRatioEnergyParameter *( ShearStressProjection_1 * Aspect_Ratio_Diff_j1 * Aspect_Ratio_Diff_j1 -ShearStressProjection_0 * Aspect_Ratio_Diff_j0 * Aspect_Ratio_Diff_j0);
        delta_E += mAspectRatioEnergyParameter  *( Aspect_Ratio_Diff_j1 * Aspect_Ratio_Diff_j1 - Aspect_Ratio_Diff_j0 * Aspect_Ratio_Diff_j0);
        delta_O += mOrientationParameter *AverageShearStress *( acos(inner_prod(ShearStress_j1, MajorAxis_j1)/ norm_2(ShearStress_j1)) - acos(inner_prod(ShearStress_j0, MajorAxis_j0)/ norm_2(ShearStress_j0)));
        // PRINT_3_VARIABLES(ShearStressProjection_0, ShearStressProjection_1, delta_H)

        }


    delta_H = delta_O +delta_E;

   
    return delta_H;
}

template <unsigned DIM>
double AspectRatioConstraintPottsUpdateRule<DIM>::GetAspectRatioEnergyParameter()
{
    return mAspectRatioEnergyParameter;
}

template <unsigned DIM>
void AspectRatioConstraintPottsUpdateRule<DIM>::SetAspectRatioEnergyParameter(double aspectRatioEnergyParameter)
{
    mAspectRatioEnergyParameter = aspectRatioEnergyParameter;
}

template <unsigned DIM>
double AspectRatioConstraintPottsUpdateRule<DIM>::GetOrientationParameter()
{
    return mOrientationParameter;
}

template <unsigned DIM>
void AspectRatioConstraintPottsUpdateRule<DIM>::SetOrientationParameter(double OrientationParameter)
{
    mOrientationParameter = OrientationParameter;
}

template <unsigned DIM>
double AspectRatioConstraintPottsUpdateRule<DIM>::GetTargetAspectRatio() const
{
    return mTargetAspectRatio;
}

template <unsigned DIM>
void AspectRatioConstraintPottsUpdateRule<DIM>::SetTargetAspectRatio(double targetAspectRatio)
{
    assert(targetAspectRatio >= 0.0);
    mTargetAspectRatio = targetAspectRatio;
}

template <unsigned DIM>
void AspectRatioConstraintPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AspectRatioEnergyParameter>" << mAspectRatioEnergyParameter << "</AspectRatioEnergyParameter>\n";
    *rParamsFile << "\t\t\t<TargetAspectRatio>" << mTargetAspectRatio << "</TargetAspectRatio>\n";

    // Call method on direct parent class
    AbstractWrappedPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AspectRatioConstraintPottsUpdateRule<1>;
template class AspectRatioConstraintPottsUpdateRule<2>;
template class AspectRatioConstraintPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AspectRatioConstraintPottsUpdateRule)

/*

Shear Stress update rule 

*/

#include "MechanotaxisPottsWrappedUpdateRule.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"

template <unsigned DIM>
MechanotaxisPottsWrappedUpdateRule<DIM>::MechanotaxisPottsWrappedUpdateRule()
        : AbstractWrappedPottsUpdateRule<DIM>(),
          mTractionCorrelationParameter(0.01) // @todo Made up defaults

{
    // TRACE("Shear constructor")
}

template <unsigned DIM>
MechanotaxisPottsWrappedUpdateRule<DIM>::~MechanotaxisPottsWrappedUpdateRule()
{
}

template <unsigned DIM>
double MechanotaxisPottsWrappedUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                                unsigned targetNodeIndex,
                                                                                WrappedPottsBasedCellPopulation<DIM>& rCellPopulation)
{

    // TRACE(" Migration shear stress --- EvaluateHamiltonianContribution ")
    // Here We assess the difference in the average shear stress on the cell

    c_vector<double, DIM> displacement_direction = rCellPopulation.GetNode(targetNodeIndex)->rGetLocation() - rCellPopulation.GetNode(currentNodeIndex)->rGetLocation();

    

    PottsArbitrarySurfaceIn3DMesh<DIM>* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&(rCellPopulation.rGetMesh())); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>

    // Favour moves along the shear stress vector (3rd term ranges [0,2] hence the 0.5 factor)
    c_vector<double, DIM> WSS = p_static_cast_potts_mesh->GetTractionOnLattice(currentNodeIndex); //
    
    // // This if pretty much Gradf.u For the directional derivative
    // double delta_H = mTractionCorrelationParameter * inner_prod(displacement_direction, WSS);
    // return delta_H;



    // New way 

    c_vector<double, DIM> DeltaSS =  p_static_cast_potts_mesh->GetTractionOnLattice(targetNodeIndex) -p_static_cast_potts_mesh->GetTractionOnLattice(currentNodeIndex) / norm_2(displacement_direction) ;

    displacement_direction /= norm_2(displacement_direction);
    // PRINT_VECTOR(DeltaSS); 
    // // PRINT_VECTOR(displacement_direction); 
    // PRINT_VARIABLE(inner_prod(displacement_direction, DeltaSS));
    // TRACE("Should this be divided by the distace??????")
    double delta_H = -mTractionCorrelationParameter * inner_prod(displacement_direction, DeltaSS);
    return delta_H ;



}

/*

// Try as James suggested with the basic gradient average 

  



    // PRINT_2_VARIABLES(delta_H ,WSS_grad );





With this method I was using the average shear stress, thinking this would be maximised as the cell went 
up the shear stress gradient -- didnt work, the cell stretched out very annoying 

 // THis is the difference in the average shear stress over the element before and after the switch 
        unsigned current_element = (*containing_elements.begin());
        PottsElement<DIM>* pCurrentElement = p_static_cast_potts_mesh->GetElement(current_element);
        N = pCurrentElement->GetNumNodes();
        T_Cell_i = p_static_cast_potts_mesh->GetSumShearStressOnElement(current_element);


        // unsigned Sister_Element = rCellPopulation.GetSister(current_element);
        // PottsElement<DIM>* pSisterElement = p_static_cast_potts_mesh->GetElement(Sister_Element);
        //  N += pSisterElement->GetNumNodes();

        // T_Cell_i +=p_static_cast_potts_mesh->GetSumShearStressOnElement(Sister_Element);

        //  double T_j = norm_2(p_static_cast_potts_mesh->GetTractionOnLattice(targetNodeIndex));
        // double T_grad = (T_j * N - T_Cell_i) / (N * (N + 1));
        // delta_H = mTractionCorrelationParameter* T_grad ;

        // PRINT_4_VARIABLES(T_grad,T_j , T_Cell_i, N)

















    James' way 
    Favour moves along the shear stress vector (3rd term ranges [0,2] hence the 0.5 factor)
    double delta_H = - mTractionCorrelationParameter * 0.5 * (1 + inner_prod(displacement_direction, WSS));


    From contact inhibition Update rule in 
    Contact-Inhibited Chemotaxis in De Novo and Sprouting Blood-Vessel Growth ---  Merks and Perryn 2008
    [not sure how this will help]
    delta_H = - mTractionCorrelationParameter* (ShearStressGrad/(1+3* ShearStressGrad) - norm_2(WSS)/ (1+3* norm_2(WSS));



    From traction with substract durotaxis 
    Cells Contractility Facilitates Alignment of cells and tissue to static uniaxial stretch ----   Ren and Merks 2017
    double WSS_thresh = 0;
    double delta_H = - mTractionCorrelationParameter* 100/( 1+exp(-1e-3*(ShearStressGrad-WSS_thresh))) *inner_prod(displacement_direction, WSS) * inner_prod(displacement_direction, WSS);
     
    My second idea 

     WSS_grad = WWS_j - Mag_WSS;
     delta_H = - mTractionCorrelationParameter *0.5* (1 + inner_prod(displacement_direction, WSS* WSS_grad));

      double WSS_grad = -(WWS_j - CurrentShearStress) / (N * (N + 1));

    
    // double delta_H = - mTractionCorrelationParameter *0.5* (1 + inner_prod(displacement_direction, WSS* WSS_grad));
  



*/

template <unsigned DIM>
double MechanotaxisPottsWrappedUpdateRule<DIM>::GetTractionCorrelationParameter()
{
    // TRACE("GetTractionCorrelationParameter ")
    return mTractionCorrelationParameter;
}

template <unsigned DIM>
void MechanotaxisPottsWrappedUpdateRule<DIM>::SetTractionCorrelationParameter(double tractionCorrelationParameter)
{
    // TRACE("SetTractionCorrelationParameter")
    mTractionCorrelationParameter = tractionCorrelationParameter;
}

template <unsigned DIM>
void MechanotaxisPottsWrappedUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<TractionCorrelationParameter>" << mTractionCorrelationParameter << "</TractionCorrelationParameter>\n";

    // Call method on direct parent class
    AbstractWrappedPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MechanotaxisPottsWrappedUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MechanotaxisPottsWrappedUpdateRule<3>)

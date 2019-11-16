/*

Shear Stress update rule 

*/

#include "MechanotaxisPottsUpdateRule.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"

template <unsigned DIM>
MechanotaxisPottsUpdateRule<DIM>::MechanotaxisPottsUpdateRule()
        : AbstractPottsUpdateRule<DIM>(),
          mTractionCorrelationParameter(0.01) // @todo Made up defaults

{
    // TRACE("Shear constructor")
}

template <unsigned DIM>
MechanotaxisPottsUpdateRule<DIM>::~MechanotaxisPottsUpdateRule()
{
}

template <unsigned DIM>
double MechanotaxisPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                         unsigned targetNodeIndex,
                                                                         PottsBasedCellPopulation<DIM>& rCellPopulation)
{
 
    // TRACE(" Migration shear stress --- EvaluateHamiltonianContribution ")
    assert(DIM == 3);

    c_vector<double, DIM> displacement_direction = rCellPopulation.GetNode(targetNodeIndex)->rGetLocation() - rCellPopulation.GetNode(currentNodeIndex)->rGetLocation();

    displacement_direction /= norm_2(displacement_direction);

    PottsArbitrarySurfaceIn3DMesh<DIM>* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&(rCellPopulation.rGetMesh())); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>

    std::set<unsigned> containing_elements = p_static_cast_potts_mesh->GetNode(currentNodeIndex)->rGetContainingElementIndices();
    unsigned current_element = (*containing_elements.begin());

    // double SS = containing_elements.size();

    // unsigned N = p_static_cast_potts_mesh->GetElement(current_element)->GetNumNodes();
    double N;
   bool current_node_contained = !containing_elements.empty();
   double T_Cell_i ;
   double delta_H;
     if (current_node_contained) // current node is in an element
    {

        unsigned current_element = (*containing_elements.begin());
        PottsElement<DIM>* pCurrentElement = p_static_cast_potts_mesh->GetElement(current_element);
        N = pCurrentElement->GetNumNodes();
        T_Cell_i = p_static_cast_potts_mesh->GetSumShearStressOnElement(current_element);

         double T_j = norm_2(p_static_cast_potts_mesh->GetTractionOnLattice(targetNodeIndex));
        double T_grad = (T_j * N - T_Cell_i) / (N * (N + 1));
        delta_H = -mTractionCorrelationParameter* T_grad ;


        c_vector<double, DIM> WSS = p_static_cast_potts_mesh->GetTractionOnLattice(currentNodeIndex);
        double Mag_WSS = norm_2(WSS);
        WSS /= norm_2(WSS);


        // double Direc = inner_prod(displacement_direction, WSS);
        // PRINT_2_VARIABLES(T_grad , Direc);
    }
    else{
    //   N=1;
    //  T_Cell_i =0;
     delta_H  = 0 ;
    }


   

    // double ShearStressGrad = norm_2(p_static_cast_potts_mesh->GetTractionOnLattice(targetNodeIndex));

    // Jess's idea

    // double delta_H = - mTractionCorrelationParameter *0.5 * (1 + inner_prod(displacement_direction, WSS* Mag_WSS));
    // double WWS_j = norm_2(p_static_cast_potts_mesh->GetTractionOnLattice(targetNodeIndex));

    // double CurrentShearStress = N * N * norm_2(p_static_cast_potts_mesh->GetTractionOnElement(current_element));
    
    // From traction with substract durotaxis  --  Cells Contractility Facilitates Alignment of cells and tissue to static uniaxial stretch ----   Ren and Merks 2017
    // double WSS_thresh = 0;
    // double delta_H = - mTractionCorrelationParameter* 100/( 1+exp(-1e-3*(ShearStressGrad-WSS_thresh))) *inner_prod(displacement_direction, WSS * Mag_WSS) * inner_prod(displacement_direction, WSS * Mag_WSS);
  




// Try as James suggested with the basic gradient average 

  



    // PRINT_2_VARIABLES(delta_H ,WSS_grad );
    return delta_H; 
}

/*
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
double MechanotaxisPottsUpdateRule<DIM>::GetTractionCorrelationParameter()
{
    // TRACE("GetTractionCorrelationParameter ")
    return mTractionCorrelationParameter;
}

template <unsigned DIM>
void MechanotaxisPottsUpdateRule<DIM>::SetTractionCorrelationParameter(double tractionCorrelationParameter)
{
    // TRACE("SetTractionCorrelationParameter")
    mTractionCorrelationParameter = tractionCorrelationParameter;
}

template <unsigned DIM>
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

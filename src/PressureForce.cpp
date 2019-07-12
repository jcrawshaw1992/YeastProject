#include "PressureForce.hpp"
#include "UblasCustomFunctions.hpp"
#include "Debug.hpp"
// Area constrant on a membrane 


PressureForce::PressureForce()
   : AbstractForce<2,3>()
{
}



void PressureForce::AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation)
{
/*
 * Itterate over all the cells, get the force, then add the force .
*/


    MeshBasedCellPopulation<2,3>* p_cell_population = static_cast<MeshBasedCellPopulation<2,3>*>(&rCellPopulation);

   for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);
        c_vector<double,3> force;


       force[0] =   cell_iter->GetCellData()->GetItem("applied_force_x");
	   force[1] =   cell_iter->GetCellData()->GetItem("applied_force_y");
	   force[2] =   cell_iter->GetCellData()->GetItem("applied_force_z");
    //    PRINT_VECTOR(force);
    //    force =100 *force;
       PRINT_VECTOR(force);
       TRACE("AM I EVEN HERE");

        pNode->AddAppliedForceContribution(force);
        // TRACE("Drag corrected");
        //  PRINT_2_VARIABLES(NormOldForce , NormNewForce);
     }
} 




void PressureForce::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2,3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PressureForce)



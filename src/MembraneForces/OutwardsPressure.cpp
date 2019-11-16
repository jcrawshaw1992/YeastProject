#include "OutwardsPressure.hpp"

#include "UblasCustomFunctions.hpp"
#include "MathsFunctions.hpp"


// #include "MathsFunctions.hpp"

OutwardsPressure::OutwardsPressure()
        : AbstractForce<2, 3>()
{
}

void OutwardsPressure::SetPressure(double Pressure)
{
      mStrength =  Pressure;
}




void OutwardsPressure::AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
{
   MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
     

       // Calculate midpoint
        // c_vector<double, 3> centroid = zero_vector<double>(3);
        // for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
        //      cell_iter != rCellPopulation.End();
        //      ++cell_iter)
        // {
        //     unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        //     centroid += rCellPopulation.GetNode(node_index)->rGetLocation();
        //    // Area += rCellPopulation.GetVolumeOfCell(*cell_iter);
        // }
        // centroid /= rCellPopulation.GetNumRealCells();
        
        for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {

            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);
            // c_vector<double, 3> cell_location = p_node->rGetLocation() - centroid;
            // cell_location(2) = 0.0;
            // cell_location /= norm_2(cell_location);

            c_vector<double, 3> Normal = zero_vector<double>(3);
            double Area = 0;
            std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
            assert(containing_elements.size() > 0);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                 iter != containing_elements.end();
                 ++iter)
            {
                Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
                Node<3>* pNode1= p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
                Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

                c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
                c_vector<double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

                c_vector<double, 3> normalVector = VectorProduct(vector_12, vector_13);
                
                // Area+= norm_2(normalVector)/6;
                Normal += normalVector;///norm_2(normalVector);

            }
            Normal /=norm_2(Normal);
             c_vector<double, 3> force = mStrength *Normal; // / norm_2(cell_location);
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);
        }

}
 

void OutwardsPressure::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(OutwardsPressure)

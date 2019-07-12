#include "AppliedForce.hpp"
#include "Debug.hpp"
#include "MeshBasedCellPopulation.hpp"
#include <cmath>
#include "UblasCustomFunctions.hpp"


#include "Debug.hpp"
#include <math.h>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AppliedForce<ELEMENT_DIM, SPACE_DIM>::AppliedForce()
    : AbstractForce<ELEMENT_DIM, SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AppliedForce<ELEMENT_DIM, SPACE_DIM>::~AppliedForce()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForce<ELEMENT_DIM, SPACE_DIM>::AddForceContribution(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    //  TRACE("Applied force");
	assert(SPACE_DIM==3); // Currently assumes that SPACE_DIM = 3

    MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);



 for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<SPACE_DIM>* pNode = p_cell_population->rGetMesh().GetNode(node_index);
 
        // double ScalarPressure = cell_iter->GetCellData()->GetItem("Pressure");
        // c_vector<long double,3> Force = ScalarPressure  * NormalVector;
        // c_vector<long double,3> Force = ScalarPressure  * cell_location;
        //       c_vector<double, 3> cell_location = pNode->rGetLocation() ;
        //     cell_location(2) = 0.0;
        //     cell_location /= norm_2(cell_location);
        // c_vector<long double,3> Force =0.01 *  cell_location;

        c_vector<long double,3> Force ;

        // TRACE("Getting Force");
        Force[0] = cell_iter->GetCellData()->GetItem("applied_force_x");
        Force[1] = cell_iter->GetCellData()->GetItem("applied_force_y");
        Force[2] = cell_iter->GetCellData()->GetItem("applied_force_z");
        double Pressure  = cell_iter->GetCellData()->GetItem("Pressure");
        // PRINT_VARIABLE(Pressure);
        
        //  Force= Create_c_vector(0.001,0.001,0);


        //  Force *= 1e-6;
        // Force *= 1e-1;
        // // if (node_index ==200)
        // {
        //     PRINT_VECTOR(Force);
        // }
        // 
      
        // double Pressure = norm_2(Force) ;
        // PRINT_VARIABLE(cell_iter->GetCellData()->GetItem("Pressure"));
        rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(Force);
    //  PRINT_VECTOR(Force);
        // applied_force[0]=cell_iter->GetCellData()->GetItem("applied_force_x");
        // applied_force[1]=cell_iter->GetCellData()->GetItem("applied_force_y");
        // applied_force[2]=cell_iter->GetCellData()->GetItem("applied_force_z");
        // applied_force = applied_force/20;
        // PRINT_VECTOR(applied_force);


        // assert(fabs(Force[0])<1e10);
        // assert(fabs(Force[1])<1e10);
        // assert(fabs(Force[2])<1e10);

    }





   
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include

    // Call method on direct parent class
    AbstractForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////
template class AppliedForce<2,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AppliedForce)



 // double MaxNodes = 20*20*1.5;
        // double Nc =20;

        //   if (node_index < Nc || node_index > MaxNodes -Nc)
        //     {
        //         TRACE("pressure too small");
        //         PRINT_VARIABLE(ScalarPressure );
        //     }

        //   if (node_index == Nc-1 || node_index == MaxNodes- Nc -2)
        //     {
        //         TRACE("Expected Pressure");
        //         PRINT_VARIABLE(ScalarPressure );
        //     }







// for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
//                      cell_iter != rCellPopulation.End();
//                      ++cell_iter)
//         {
//             TRACE("Have an abstract cell populaiton ");
//             unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
//             Node<SPACE_DIM>* p_node = rCellPopulation.GetNode(node_index);
           
//             c_vector<double, SPACE_DIM> force = zero_vector<double>(SPACE_DIM);

//             //Calculate cell normal (average of element normals)
//             c_vector<double,SPACE_DIM> normal = zero_vector<double>(SPACE_DIM);

//             std::set<unsigned>&  containing_elements = p_node->rGetContainingElementIndices();
//             assert(containing_elements.size()>0);
//             for (std::set<unsigned>::iterator iter = containing_elements.begin();
//                     iter != containing_elements.end();
//                     ++iter)
//             {
//                 // Negative as normals point inwards for these surface meshes
//                 normal +=  -p_cell_population->rGetMesh().GetElement(*iter)->CalculateNormal();
//             }
//             normal /= norm_2(normal);
//             double Pressure  = cell_iter->GetCellData()->GetItem("Pressure");
//             force = Pressure  * normal; // cell_location / norm_2(cell_location);
            
//             rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);

//             cell_iter->GetCellData()->SetItem("applied_force_x", force[0]);
//             cell_iter->GetCellData()->SetItem("applied_force_y", force[1]);
//             cell_iter->GetCellData()->SetItem("applied_force_z", force[2]);
//             cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(force));

//             cell_iter->GetCellData()->SetItem("norm_x", normal[0]);
//             cell_iter->GetCellData()->SetItem("norm_y", normal[1]);
//             cell_iter->GetCellData()->SetItem("norm_z", normal[2]);
//         }
    
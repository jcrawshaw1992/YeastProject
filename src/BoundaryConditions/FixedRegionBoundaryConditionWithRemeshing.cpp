#include "FixedRegionBoundaryConditionWithRemeshing.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "CellLabel.hpp"
#include "Debug.hpp"
#include "VertexBasedCellPopulation.hpp"

// Normal needs to be pointing out, point needs be be close to the mesh/region, and the radius will give you how far away you need to be :)

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FixedRegionBoundaryConditionWithRemeshing<ELEMENT_DIM, SPACE_DIM>::FixedRegionBoundaryConditionWithRemeshing(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
                                                                                                             c_vector<double, SPACE_DIM> point,
                                                                                                             c_vector<double, SPACE_DIM> normal,
                                                                                                             double radius)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(pCellPopulation),
          mPointOnPlane(point),
          mRadius(radius)
{
    assert(norm_2(normal) > 0.0);
    assert(radius > 0.0);
    mNormalToPlane = normal / norm_2(normal);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& FixedRegionBoundaryConditionWithRemeshing<ELEMENT_DIM, SPACE_DIM>::rGetPointOnPlane() const
{
    return mPointOnPlane;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& FixedRegionBoundaryConditionWithRemeshing<ELEMENT_DIM, SPACE_DIM>::rGetNormalToPlane() const
{
    return mNormalToPlane;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const double& FixedRegionBoundaryConditionWithRemeshing<ELEMENT_DIM, SPACE_DIM>::rGetRadius() const
{
    return mRadius;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FixedRegionBoundaryConditionWithRemeshing<ELEMENT_DIM, SPACE_DIM>::ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations)
{

    ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation) == NULL)
    {
        EXCEPTION("FixedRegionBoundaryConditionWithRemeshing requires a subclass of AbstractOffLatticeCellPopulation.");
    }
    if (SPACE_DIM != 1)
    {

        //assert(dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation));
        // AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>* pStaticCastCellPopulation = static_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation);
        HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pStaticCastCellPopulation = static_cast<HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation);


        // WHen RemeshedStatis ==0 then there has just been a remeshing event, and there is no last condition to call, and this one needs to be skipped 
        bool RemeshedStatis = pStaticCastCellPopulation->GetUpdateBoundaryConditions();
 
        if (RemeshedStatis == 1)
        {
            if (dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation))
            {
                for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
                    cell_iter != this->mpCellPopulation->End();
                    ++cell_iter)
                {
                    unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
                    Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

                    c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();
                    c_vector<double, SPACE_DIM> node_old_location = rOldLocations.find(p_node)->second;

                    double signed_distance = inner_prod(node_location - mPointOnPlane, mNormalToPlane);
                    if (signed_distance > 0.0)
                    {
                        c_vector<double,SPACE_DIM> displacement = node_location - node_old_location;
                        p_node->rGetModifiableLocation() = node_old_location + displacement - inner_prod(displacement,mNormalToPlane)*mNormalToPlane;
                        cell_iter->GetCellData()->SetItem("FixedBoundary", 1);
                    } 
                   
    
                ////
                }
            }


        }else
        {

            pStaticCastCellPopulation->UpdateBoundaryConditions(); // this will reset RemeshedStatis to 1. if RemeshedStatis ==0 then the remeshing step happened very recently
        }
        //  bool RemeshedStatis = pStaticCastCellPopulation->GetUpdateBoundaryConditions();
    }
}





// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void FixedRegionBoundaryConditionWithRemeshing<ELEMENT_DIM, SPACE_DIM>::ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations)
// {

//     ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
//     if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation) == NULL)
//     {
//         EXCEPTION("FixedRegionBoundaryConditionWithRemeshing requires a subclass of AbstractOffLatticeCellPopulation.");
//     }
//     if (SPACE_DIM != 1)
//     {

//         //assert(dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation));
//         // AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>* pStaticCastCellPopulation = static_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation);
//         HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pStaticCastCellPopulation = static_cast<HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation);


//         // WHen RemeshedStatis ==0 then there has just been a remeshing event, and there is no last condition to call, and this one needs to be skipped 
//         bool RemeshedStatis = pStaticCastCellPopulation->GetUpdateBoundaryConditions();
//         // PRINT_VARIABLE(RemeshedStatis) 
// // RemeshedStatis = 0
//         if (RemeshedStatis == 1)
//         {
//             unsigned num_nodes = pStaticCastCellPopulation->GetNumNodes();
//             for (unsigned node_index = 0; node_index < num_nodes; node_index++)
//             {
//                 Node<SPACE_DIM>* p_node = pStaticCastCellPopulation->GetNode(node_index);
//                 c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();
//                 // Find the old location of the node
//                 c_vector<double, SPACE_DIM> node_old_location = rOldLocations.find(p_node)->second;

//                 // calculate the distance of the node from the central point if its to far then dont impose boundary condition
//                 double distance_from_point = norm_2(node_old_location - mPointOnPlane);

//                 if (distance_from_point <= mRadius)
//                 {
//                     double signed_distance = inner_prod(node_location - mPointOnPlane, mNormalToPlane);
//                     if (signed_distance > 0.0)
//                     {

//                         // //  // For the closest point on the plane we travel from node_location the signed_distance in the direction of -mNormalToPlane
//                         // // c_vector<double, SPACE_DIM> nearest_point = node_location - signed_distance*mNormalToPlane;

//                         // // p_node->rGetModifiableLocation() = nearest_point;
//                         // // // Only allow movement perpendicular to the normal -- for some reason this is not working
//                         // c_vector<double, SPACE_DIM> displacement = node_location - node_old_location;
//                         // p_node->rGetModifiableLocation() = node_old_location + displacement - inner_prod(displacement, mNormalToPlane) * mNormalToPlane;
//                         // // p_node->rGetModifiableLocation() = node_location - inner_prod(displacement,mNormalToPlane)*mNormalToPlane;

//                         // CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(node_index);
//                         // //  p_cell->RemoveCellProperty<CellLabel>();

//                         // // Add label to cell so can see them in Visualizer
//                         // // boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
//                         // // this->mpCellPopulation->GetCellUsingLocationIndex(node_index)->AddCellProperty(p_label); // else


//                     // Only allow movement perpendicular to the normal
//                     c_vector<double,SPACE_DIM> displacement = node_location - node_old_location;

//                     p_node->rGetModifiableLocation() = node_old_location + displacement - inner_prod(displacement,mNormalToPlane)*mNormalToPlane;

//                     // Uncomment this code if you want iolets labelled in Paraview for viz purposes
                    
//                     // Remove any CellLabels as can only have one label
//                     CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(node_index);
//                     p_cell->RemoveCellProperty<CellLabel>();
    
    
//                     // Add label to cell so can see them in Visualizer
//                     boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
//                     this->mpCellPopulation->GetCellUsingLocationIndex(node_index)->AddCellProperty(p_label);



//                         p_cell->GetCellData()->SetItem("FixedBoundary", 1);
//                         //  TRACE("A")
//                     }
//                 }
//             }
//         }else
//         {

//             pStaticCastCellPopulation->UpdateBoundaryConditions(); // this will reset RemeshedStatis to 1. if RemeshedStatis ==0 then the remeshing step happened very recently
//         }
//         //  bool RemeshedStatis = pStaticCastCellPopulation->GetUpdateBoundaryConditions();
//     }
// }

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool FixedRegionBoundaryConditionWithRemeshing<ELEMENT_DIM, SPACE_DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;
    // TODO implement check that the fixed points dont move

    return condition_satisfied;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FixedRegionBoundaryConditionWithRemeshing<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PointOnPlane>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mPointOnPlane[index] << ",";
    }
    *rParamsFile << mPointOnPlane[SPACE_DIM - 1] << "</PointOnPlane>\n";

    *rParamsFile << "\t\t\t<NormalToPlane>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mNormalToPlane[index] << ",";
    }
    *rParamsFile << mNormalToPlane[SPACE_DIM - 1] << "</NormalToPlane>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class FixedRegionBoundaryConditionWithRemeshing<1, 1>;
template class FixedRegionBoundaryConditionWithRemeshing<1, 2>;
template class FixedRegionBoundaryConditionWithRemeshing<2, 2>;
template class FixedRegionBoundaryConditionWithRemeshing<1, 3>;
template class FixedRegionBoundaryConditionWithRemeshing<2, 3>;
template class FixedRegionBoundaryConditionWithRemeshing<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(FixedRegionBoundaryConditionWithRemeshing)

// Uncomment this code if you want iolets labelled in Paraview for viz purposes

// // Remove any CellLabels as can only have one label
// CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(node_index);
// //  if (mIntialTime==1)
// //  {
//         p_cell->GetCellData()->SetItem("BoundarySet", 1);
// TRACE("new boundary node")
//  }

// p_cell->RemoveCellProperty<CellLabel>();
// // TRACE("Stuffremoved here")

// // Add label to cell so can see them in Visualizer
// boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
// this->mpCellPopulation->GetCellUsingLocationIndex(node_index)->AddCellProperty(p_label); // else
//      {
//          CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(node_index);
//          p_cell->GetCellData()->SetItem("Boundary", 0);
//      }
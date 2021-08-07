#include "EnclosedRegionBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellLabel.hpp"
// #include "Debug.hpp"

// Normal needs to be pointing out, point needs be be close to the mesh/region, and the radius will give you how far away you need to be :) 

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
EnclosedRegionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::EnclosedRegionBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* pCellPopulation,
                                                    c_vector<double, SPACE_DIM> point,
                                                    c_vector<double, SPACE_DIM> normal,
                                                    c_vector<double, SPACE_DIM> point2,
                                                    c_vector<double, SPACE_DIM> normal2,
                                                    double radius)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>(pCellPopulation),
          mPointOnPlane(point),
          mRadius(radius)
{
    assert(norm_2(normal) > 0.0);
    assert(radius>0.0);
    mNormalToPlane = normal/norm_2(normal);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& EnclosedRegionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::rGetPointOnPlane() const
{
    return mPointOnPlane;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& EnclosedRegionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::rGetNormalToPlane() const
{
    return mNormalToPlane;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& EnclosedRegionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::rGetPointOnPlane2() const
{
    return mPointOnPlane2;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& EnclosedRegionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::rGetNormalToPlane2() const
{
    return mNormalToPlane2;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const double& EnclosedRegionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::rGetRadius() const
{
    return mRadius;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void EnclosedRegionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations)
{
    ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation)==NULL)
    {
        EXCEPTION("EnclosedRegionBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }
    if (SPACE_DIM != 1)
    {
		//assert(dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation));
		AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>* pStaticCastCellPopulation = static_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation);
		// Iterate over all nodes and update their positions according to the boundary conditions
		unsigned num_nodes = pStaticCastCellPopulation->GetNumNodes();
		for (unsigned node_index=0; node_index<num_nodes; node_index++)
		{
			Node<SPACE_DIM>* p_node = pStaticCastCellPopulation->GetNode(node_index);
			c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();          
			// Find the old location of the node
			c_vector<double, SPACE_DIM> node_old_location = rOldLocations.find(p_node)->second;

            // calculate the distance of the node from the central point if its to far then dont impose boundary condition
            double distance_from_point = norm_2(node_old_location - mPointOnPlane);
            
            if (distance_from_point <= mRadius)
            {    
                double signed_distance = inner_prod(node_location - mPointOnPlane, mNormalToPlane);
                double signed_distance2 = inner_prod(node_location - mPointOnPlane2, mNormalToPlane2);
                if (signed_distance >= 0.0 && signed_distance<=0.0)
                {

                    //  // For the closest point on the plane we travel from node_location the signed_distance in the direction of -mNormalToPlane
                    // c_vector<double, SPACE_DIM> nearest_point = node_location - signed_distance*mNormalToPlane;
    
                    // p_node->rGetModifiableLocation() = nearest_point;
                    // // Only allow movement perpendicular to the normal -- for some reason this is not working 
                    c_vector<double,SPACE_DIM> displacement = node_location - node_old_location;
                    p_node->rGetModifiableLocation() = node_old_location + displacement - inner_prod(displacement,mNormalToPlane)*mNormalToPlane;
                    // p_node->rGetModifiableLocation() = node_location - inner_prod(displacement,mNormalToPlane)*mNormalToPlane
                    CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(node_index);
                    p_cell->GetCellData()->SetItem("FixedBoundary", 1);
                    //  TRACE("A")
                 }
              
            } 
           
		}
    }

}




template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool EnclosedRegionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;
    // TODO implement check that the fixed points dont move

    return condition_satisfied;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void EnclosedRegionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PointOnPlane>";
    for (unsigned index=0; index != SPACE_DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mPointOnPlane[index] << ",";
    }
    *rParamsFile << mPointOnPlane[SPACE_DIM-1] << "</PointOnPlane>\n";

    *rParamsFile << "\t\t\t<NormalToPlane>";
     for (unsigned index=0; index != SPACE_DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
     {
         *rParamsFile << mNormalToPlane[index] << ",";
     }
     *rParamsFile << mNormalToPlane[SPACE_DIM-1] << "</NormalToPlane>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class EnclosedRegionBoundaryCondition<1,1>;
template class EnclosedRegionBoundaryCondition<1,2>;
template class EnclosedRegionBoundaryCondition<2,2>;
template class EnclosedRegionBoundaryCondition<1,3>;
template class EnclosedRegionBoundaryCondition<2,3>;
template class EnclosedRegionBoundaryCondition<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(EnclosedRegionBoundaryCondition)


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
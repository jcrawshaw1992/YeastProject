#include "PottsArbitrarySurfaceIn3DMesh.hpp"


template<unsigned SPACE_DIM>
PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::PottsArbitrarySurfaceIn3DMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                                                        std::vector<PottsElement<SPACE_DIM>*> pottsElements,
                                                                        std::vector< std::set<unsigned> > vonNeumannNeighbouringNodeIndices,
                                                                        std::vector< std::set<unsigned> > mooreNeighbouringNodeIndices,
                                                                        MutableMesh<2,SPACE_DIM>* pDelaunayMesh) :
    PottsMesh<SPACE_DIM>(nodes, pottsElements, vonNeumannNeighbouringNodeIndices, mooreNeighbouringNodeIndices),
    mpDelaunayMesh(pDelaunayMesh)
{
    assert(SPACE_DIM==2 || SPACE_DIM==3);

    mMeshElementMidPoints.clear();
}

template<unsigned SPACE_DIM>
PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::~PottsArbitrarySurfaceIn3DMesh()
{
}


template<unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetVolumeOfLatticeSite(unsigned latticeSiteIndex)
{
    // TRACE(" UPDATED ------  GetVolumeOfLatticeSite");
    // TRACE("GetVolumeOfLatticeSite has been updated for Potts nodes on Mesh nodes");
     // TODO: This method's output can be computed in the constructor and cached. Cache will need updating after a call to UpdatePottsNodeLocationFromDelaunay
    // Element<2,SPACE_DIM>* mesh_element = mpDelaunayMesh->GetElement(latticeSiteIndex);
    
    // Jess editied this so that the area of the potts lattice site was caclulated based on the area of 
    // the mesh node, rather than the area of the elements 
    long double Node_Area = 0;
    Node<SPACE_DIM> * p_node =  mpDelaunayMesh->GetNode(latticeSiteIndex);
    std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
    assert(containing_elements.size() > 0);

    for (std::set<unsigned>::iterator iter = containing_elements.begin();
        iter != containing_elements.end();
        ++iter)

    {   unsigned element_index = *(iter);
        Element<2,SPACE_DIM>* mesh_element = mpDelaunayMesh->GetElement(element_index );

        // Loops over the elements and calc the area
        c_matrix<double, SPACE_DIM, 2> jacobian;
        double determinant;
        mesh_element->CalculateJacobian(jacobian, determinant);
        double VolumeOfele = mesh_element->GetVolume(determinant);
        Node_Area +=VolumeOfele/3;       
    }

    return Node_Area;
}



template<unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetVolumeOfElement(unsigned pottsElementIndex)
{
    // TRACE(" UPDATED ------  GetVolumeOfElement");
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);

    double potts_element_volume = 0.0;

    // An element is made of a number of lattice sites, which are centered around a number of nodes in the Delaunay mesh
    for (unsigned node_index = 0; node_index < p_potts_element->GetNumNodes(); ++node_index)
    {
        potts_element_volume += GetVolumeOfLatticeSite(p_potts_element->GetNodeGlobalIndex(node_index));
    }

    return potts_element_volume;
}

template<unsigned SPACE_DIM>
inline double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::ComputeTriangleArea(const c_vector<double, SPACE_DIM>& vertexA, const c_vector<double, SPACE_DIM>& vertexB, const c_vector<double, SPACE_DIM>& vertexC) const
{
    const c_vector<double, SPACE_DIM> AC = vertexC-vertexA;
    const c_vector<double, SPACE_DIM> AB = vertexB-vertexA;
    return 0.5 * norm_2(VectorProduct(AC, AB));
}

template<unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetPerimeterOfLatticeSite(unsigned latticeSiteIndex)
{
        
    // TRACE("Updating now ");
    // TODO: This method's output can be computed in the constructor and cached. Cache will need updating after a call to UpdatePottsNodeLocationFromDelaunay
    Element<2,SPACE_DIM>* mesh_element = mpDelaunayMesh->GetElement(latticeSiteIndex);
    // TRACE("A");

    /// TODO make this a template parameter and code up the ELEMENT_DIM=3 case
    const unsigned ELEMENT_DIM = 2;
    // TRACE("B");

    double lattice_site_area = 0.0;


    for (unsigned local_node_index=0; local_node_index<ELEMENT_DIM+1; ++local_node_index)
    {
        // TRACE("C");
        unsigned node_a_local_index = local_node_index;
        unsigned node_b_local_index = (local_node_index+1) % (ELEMENT_DIM+1);

        const c_vector<double, SPACE_DIM>& node_a_location = mesh_element->GetNode(node_a_local_index)->rGetLocation();
        const c_vector<double, SPACE_DIM>& node_b_location = mesh_element->GetNode(node_b_local_index)->rGetLocation();

        lattice_site_area += norm_2(node_a_location - node_b_location); // what??? shouldnt this be the cross XXXX

    }

// TRACE("D");
//     PRINT_VARIABLE(lattice_site_area);


   long double Lattice_Perimiter = 0;
    Node<SPACE_DIM> * p_node =  mpDelaunayMesh->GetNode(latticeSiteIndex);

    c_vector<double , 3> Position = p_node->rGetLocation();
    // PRINT_VECTOR(Position );
    std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
    std::vector<unsigned> Elements;
    Elements.assign(containing_elements.begin(), containing_elements.end());
    // PRINT_VECTOR(Elements);

    // NeighboursToNode.erase(NeighboursToNode.begin()+i);

    assert(containing_elements.size() > 0);
    // unsigned NumberOfElements = Elements.size();
    // const unsigned p = NumberOfElements;
    // PRINT_2_VARIABLES(NumberOfElements, p  );
   
    c_vector<c_vector<double , 3> , 6 > MidPoints; 
    c_vector<c_vector<double , 3> , 3 > LocationVector; 
    // // c_vector<c_vector<Element<2,SPACE_DIM>* , 2>  , 6 > NeighbourELement;
    // std::map<unsigned, std::vector<Element<2,SPACE_DIM>*>   > NeighbourELementMap;
    
    // TRACE("Here");
    double j=0;
    std::map<unsigned, std::vector<unsigned>   > NeighbourELementMap;
    double LatticePerimiter =0;
// TRACE("E");

    // Loop over the containing elements 
    for (std::set<unsigned>::iterator iter = containing_elements.begin();
        iter != containing_elements.end();
        ++iter)

    {   
        unsigned element_index = *(iter);
        //  unsigned My_New_Element_index = *(iter)->GetIndex();
        //  PRINT_2_VARIABLE(element_index , My_New_Element_index);
        Element<2,SPACE_DIM>* mesh_element = mpDelaunayMesh->GetElement(element_index);
        // Need to get the nodes of this element 
        c_vector<c_vector<double, SPACE_DIM>, 3> Locaiton;
        c_vector<double , 3> CenterPoint = Create_c_vector(0,0,0);
        c_vector<double , 3> NodeIncides;
        double LocalCenterIndex;
// TRACE("F");
        for (int i=0; i<3; i++)
        {
            Node<SPACE_DIM>* pNode = mpDelaunayMesh->GetNode(mesh_element->GetNodeGlobalIndex(i));
            CenterPoint += pNode->rGetLocation();
            LocationVector[i]=pNode->rGetLocation();
            NodeIncides[i] = pNode->GetIndex();
            if ( NodeIncides[i] == latticeSiteIndex)
            {
                // PRINT_VECTOR(LocationVector[i])
                LocalCenterIndex = i;
            }
        }
        // Now we can get the center point 
        MidPoints[j]=CenterPoint/3;

        // Find which index is the index of interest 
        c_vector<double , 3> N1= Create_c_vector(0,0,0);
        c_vector<double , 3> N2= Create_c_vector(0,0,0);
        c_vector<double , 3> MidPoint = MidPoints[j] - LocationVector[LocalCenterIndex];
        bool HaveN1 =0;
        // Translate everything such that the initial lattice site is the origin
        for (int i=0; i<3; i++)
        {
            LocationVector[i] -= LocationVector[LocalCenterIndex];
            if ( i != LocalCenterIndex )
            {
                if (HaveN1 ==0)
                {
                 N1 = LocationVector[i];
                 HaveN1 =1;
                }
                else 
                {
                 N2 = LocationVector[i];
                }
            }
        }
        j+=1;   

        // Vector projections 

        c_vector<double , 3> b1 =(inner_prod(N1,MidPoint)/ inner_prod(N1,N1) * N1)- MidPoint;
        
        c_vector<double , 3> b2 =(inner_prod(N2,MidPoint)/ inner_prod(N2,N2) * N2)- MidPoint;
        // PRINT_2_VARIABLES(norm_2(b1), norm_2(b2));
        // LatticePerimiter += norm_2(b1)+norm_2(b2);
    }




// PRINT_2_VARIABLES(LatticePerimiter, lattice_site_area);
    
    return lattice_site_area;
}

template<unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetPerimeterOfElement(unsigned pottsElementIndex)
{
    // TRACE("A");
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);

    double potts_element_surface = 0.0;

    
    // PRINT_VARIABLE(p_potts_element->GetNumNodes());
    for (unsigned local_lattice_index = 0; local_lattice_index < p_potts_element->GetNumNodes(); ++local_lattice_index)
    {
        unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(local_lattice_index);
// TRACE("A");
        // I think the calculation for the surface area is a) based on the element centered Potts model and b) wrong -> Needs a cross product in there 
        // potts_element_surface += GetSurfaceAreaOfLatticeSite(lattice_site_index);
        potts_element_surface += GetVolumeOfLatticeSite(lattice_site_index);

// TRACE("B");
        // Im not sure if this is necessary XXXXX
        // Why would we subtract the overlapping but --- what overlapping bits? 

        // for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[lattice_site_index].begin();
        //      neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[lattice_site_index].end();
        //      ++neighbour_lattice_index)
        // {
        //     if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index))
        //     {
        //         potts_element_surface -= GetContactAreaBetweenLatticeSite(lattice_site_index,
        //                                                                   *neighbour_lattice_index);
        //     }
        // }
    }

    return potts_element_surface;
}

template<unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetSurfaceAreaOfLatticeSite(unsigned latticeSiteIndex)
{
    TRACE(" NOT  UPDATED ------  GetSurfaceAreaOfLatticeSite");
    // TODO: This method's output can be computed in the constructor and cached. Cache will need updating after a call to UpdatePottsNodeLocationFromDelaunay
    Element<2,SPACE_DIM>* mesh_element = mpDelaunayMesh->GetElement(latticeSiteIndex);

    /// TODO make this a template parameter and code up the ELEMENT_DIM=3 case
    const unsigned ELEMENT_DIM = 2;

    double lattice_site_area = 0.0;

    for (unsigned local_node_index=0; local_node_index<ELEMENT_DIM+1; ++local_node_index)
    {
        unsigned node_a_local_index = local_node_index;
        unsigned node_b_local_index = (local_node_index+1) % (ELEMENT_DIM+1);

        const c_vector<double, SPACE_DIM>& node_a_location = mesh_element->GetNode(node_a_local_index)->rGetLocation();
        const c_vector<double, SPACE_DIM>& node_b_location = mesh_element->GetNode(node_b_local_index)->rGetLocation();

        lattice_site_area += norm_2(node_a_location - node_b_location); // what??? shouldnt this be the cross XXXX

    }

    return lattice_site_area;
}




// Havent fixed because I dont use, this will cause problems later 
template<unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetSurfaceAreaOfElement(unsigned pottsElementIndex)
{
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);

    double potts_element_surface = 0.0;

    // TRACE(" NOT  UPDATED ------  GetSurfaceAreaOfElement");

    for (unsigned local_lattice_index = 0; local_lattice_index < p_potts_element->GetNumNodes(); ++local_lattice_index)
    {
        unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(local_lattice_index);

        // I think the calculation for the surface area is a) based on the element centered Potts model and b) wrong -> Needs a cross product in there 
        // potts_element_surface += GetSurfaceAreaOfLatticeSite(lattice_site_index);
        potts_element_surface += GetPerimeterOfLatticeSite(lattice_site_index);


        // Im not sure if this is necessary XXXXX
        // Why would we subtract the overlapping but --- what overlapping bits? 

        for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[lattice_site_index].begin();
             neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[lattice_site_index].end();
             ++neighbour_lattice_index)
        {
            if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index))
            {
                potts_element_surface -= GetContactAreaBetweenLatticeSite(lattice_site_index,
                                                                          *neighbour_lattice_index);
            }
        }
    }

    return potts_element_surface;
}









template<unsigned SPACE_DIM>
bool PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::DoNodesShareElement(unsigned indexNodeA, unsigned indexNodeB)
{
    std::set<unsigned> node_a_elements = this->GetNode(indexNodeA)->rGetContainingElementIndices();
    std::set<unsigned> node_b_elements = this->GetNode(indexNodeB)->rGetContainingElementIndices();

    bool is_a_medium = node_a_elements.empty();
    bool is_b_medium = node_b_elements.empty();

    if (is_a_medium xor is_b_medium)
    {
        return false;
    }

    if (is_a_medium and is_b_medium)
    {
        return true;
    }
    else
    {
        assert(node_a_elements.size() == 1);
        assert(node_b_elements.size() == 1);
        return *node_a_elements.begin() == *node_b_elements.begin();
    }
}




template<unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetContactAreaBetweenLatticeSite(unsigned index_A, unsigned index_B)
{
    assert(index_A != index_B);
    // TRACE("Need to check this");

    /// TODO make this a template parameter and code up the ELEMENT_DIM=3 case
    const unsigned ELEMENT_DIM = 2;


    // Get the elements for node a and b 
    Element<2,SPACE_DIM>* mesh_element_a = mpDelaunayMesh->GetElement(index_A);
    Element<2,SPACE_DIM>* mesh_element_b = mpDelaunayMesh->GetElement(index_B);

    
    std::vector<Node<SPACE_DIM>*> node_set_a(ELEMENT_DIM+1);
    std::vector<Node<SPACE_DIM>*> node_set_b(ELEMENT_DIM+1);

    
    // Get the pointers for the three nodes containted in element A and B
    for (unsigned local_node_index=0; local_node_index<ELEMENT_DIM+1; ++local_node_index)
    {
        node_set_a[local_node_index] = mesh_element_a->GetNode(local_node_index);
        node_set_b[local_node_index] = mesh_element_b->GetNode(local_node_index);
    }
    
    // Order the nodes in each element
    std::sort(node_set_a.begin(), node_set_a.end());
    std::sort(node_set_b.begin(), node_set_b.end());
   
    // PRINT_VECTOR(node_set_a);
    // PRINT_VECTOR(node_set_b);


    // Find the common nodes
    std::vector<Node<SPACE_DIM>*> node_set_intersecion(ELEMENT_DIM);
    typename std::vector<Node<SPACE_DIM>*>::iterator node_set_intersection_end;
    node_set_intersection_end = std::set_intersection(node_set_a.begin(), node_set_a.end(),
                                                      node_set_b.begin(), node_set_b.end(),
                                                      node_set_intersecion.begin());

    node_set_intersecion.resize(node_set_intersection_end - node_set_intersecion.begin());
    // PRINT_VECTOR(node_set_intersecion);
    


    switch (node_set_intersecion.size())
    {
        case 3:
            /// TODO code up for the ELEMENT_DIM=3 case. Use ComputeTriangleArea helper method
            NEVER_REACHED;
        case 2:
        {
            const c_vector<double, SPACE_DIM>& node_a_location = node_set_intersecion[0]->rGetLocation();
            const c_vector<double, SPACE_DIM>& node_b_location = node_set_intersecion[1]->rGetLocation();
            // TRACE("Two nodes shared between thesments");
            return norm_2(node_a_location - node_b_location);
        }
        case 1:
        // TRACE("One nodes shared between the elements -- which is defult becuase there neighbours ");
        case 0:
        // TRACE("No nodes shared between the elements --- this is bad??");
            return 0.;
        default:
            NEVER_REACHED;
    }
}




template<unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::UpdatePottsNodeLocationFromDelaunay()
{
    TRACE("Should not be here");

    error

    // mMyMember = 7;
    // This code was edited by Jess. 
    // This function updates the node locations in the potts mesh based on where the nodes have been moved to in the Delaunay mesh 
    typename MutableMesh<2,SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
    typename PottsMesh<SPACE_DIM>::NodeIterator potts_node_iter = this->GetNodeIteratorBegin();

    // At some point, I want to come through there and make some artifical changes in the 
    // node position so that I can check things work 

// TRACE("HERE")

    for(;
        node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
        ++node_iter, ++potts_node_iter)
    {
            
            // PRINT_2_VARIABLES(potts_node_iter->GetIndex(), node_iter->GetIndex());
            assert(potts_node_iter->GetIndex() == node_iter->GetIndex());
            potts_node_iter->rGetModifiableLocation() = node_iter->rGetLocation() ;
    }


}

// template<unsigned SPACE_DIM>
// void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::UpdatePottsVolumeAndArea(MutableMesh<2,SPACE_DIM>& rMesh)
// {
  
//     // This function updates the node locations in the potts mesh based on where the nodes have been moved to in the Delaunay mesh 
//     typename MutableMesh<2,SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
//     typename PottsMesh<SPACE_DIM>::NodeIterator potts_node_iter = this->GetNodeIteratorBegin();


// // Iterate over the Mesh elements to start with 

// // for (typename MutableMesh<2,SPACE_DIM>::ElementIterator elem_iter = mpDelaunayMesh->GetElementIteratorBegin();
// //             elem_iter != rMesh.GetElementIteratorEnd();
// //             ++elem_iter)
// //             {
// //                 // TRACE("hey");
// //                 // Calculate the center point for each node 
// //                 unsigned Element_index = elem_iter->GetIndex();
// //                 c_vector<double , 3> CenterPoint = Create_c_vector(0,0,0);
// //                 // c_vector<double , 3> NodeIncides;

// //                 for (int i=0; i<3; i++)
// //                 {
// //                     Node<SPACE_DIM>* pNode = rMesh.GetNode(rMesh.GetElement(Element_index)->GetNodeGlobalIndex(i));
// //                     CenterPoint += pNode->rGetLocation();
// //                     // NodeIncides[i] = pNode->GetIndex();  
// //                 }
// //                 mMeshElementMidPoints[Element_index] = CenterPoint /3;
// //             }


// TRACE("Jess is the best")



// }



template<unsigned SPACE_DIM>
MutableMesh<2,SPACE_DIM>* PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetDelaunayMesh()
{
	return mpDelaunayMesh;
}


template class PottsArbitrarySurfaceIn3DMesh<2>;
template class PottsArbitrarySurfaceIn3DMesh<3>;

//#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS1(PottsArbitrarySurfaceIn3DMesh, 2)
//EXPORT_TEMPLATE_CLASS1(PottsArbitrarySurfaceIn3DMesh, 3)











//  Code graveyard 

    //     // Figure out which elements are the nighbours to this one -- should be one or two neighbours 
    //     double a= 0;
    //     c_vector<std::vector<unsigned> , 6 > NeighbourELement;
    //     for (std::vector<unsigned>::iterator  Neighbour_iter = Elements.begin();
    //      Neighbour_iter != Elements.end();
    //     ++ Neighbour_iter)

    //     {   
    //         c_vector<double , 3> NeighbourElementNodeIncides;
    //         std::vector<unsigned> CommonNodes;
            
    //         // TRACE("Iterating over the nodes in each neighbouring mesh element ");
    //         Element<2,SPACE_DIM>* Neighbour_element = mpDelaunayMesh->GetElement(*  Neighbour_iter);
    //         unsigned N_Element_index = Neighbour_element->GetIndex();
  
    //         for (int i=0; i<3; i++)
    //         {
    //             Node<SPACE_DIM>* pNode = mpDelaunayMesh->GetNode(Neighbour_element->GetNodeGlobalIndex(i));
    //             NeighbourElementNodeIncides[i]  = pNode->GetIndex();
    //         }
            
    //         // See if there are any nodes from element iter intersecting with the original nodes
    //          std::set_intersection(NeighbourElementNodeIncides.begin(), NeighbourElementNodeIncides.end(),
    //                           NodeIncides.begin(), NodeIncides.end(),
    //                          std::back_inserter(CommonNodes));

    //         if(CommonNodes.size()>1 & CommonNodes.size()!=3) // IF the two elements share 2 nodes, they are neighbours 
    //         {
    //             TRACE("Common nodes greater than 1");  
    //            NeighbourELementMap[element_index].push_back(N_Element_index );
    //            NeighbourELementMap[N_Element_index ].push_back(element_index);

    //         }
    //         // Elements.erase(Neighbour_iter);
    //         // PRINT_VECTOR(NeighbourElementNodeIncides);
    //         // PRINT_VECTOR(NodeIncides);
    //         // PRINT_VECTOR(CommonNodes);
    //     }
    //     TRACE("Done looping");

    // }

    // for (std::set<unsigned>::iterator iter = containing_elements.begin();
    //     iter != containing_elements.end();
    //     ++iter)

    // {    
    //         Element<2,SPACE_DIM>* Neighbour_element = mpDelaunayMesh->GetElement(* iter);
    //         unsigned N_Element_index = Neighbour_element->GetIndex();
    //        PRINT_VARIABLE(N_Element_index);
    //         sort(NeighbourELementMap[N_Element_index].begin(), NeighbourELementMap[N_Element_index].end());
    //         NeighbourELementMap[N_Element_index].erase(unique(NeighbourELementMap[N_Element_index].begin(), NeighbourELementMap[N_Element_index].end()), NeighbourELementMap[N_Element_index].end());
    //         PRINT_VECTOR(NeighbourELementMap[N_Element_index]);

          
    //     }



/*
This code remeshes the vascular geomertey when trigger has been called 
*/

#include "RemeshingModifier.hpp"


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RemeshingModifier<ELEMENT_DIM, SPACE_DIM>::RemeshingModifier()
        : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RemeshingModifier<ELEMENT_DIM, SPACE_DIM>::~RemeshingModifier()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{

    TRACE("In the set up solve for the Remeshing -- dont know if its worth doing much here, record the inital conditions ")
    assert(ELEMENT_DIM == 2);
    assert(SPACE_DIM == 3);
    // std::map<unsigned, c_vector<unsigned, 5> > mNearestNodesMap;

    // for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
    //      cell_iter != rCellPopulation.End();
    //      ++cell_iter)
    // {
      
       
    //     if (cell_iter->GetCellData()->GetItem("Boundary") == 1)
    //     {

    //         unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //         Node<SPACE_DIM>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);
    //         c_vector<long double, SPACE_DIM> CellLocation = pNode->rGetLocation();

        
    //         for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter_Search = rCellPopulation.Begin();
    //              cell_iter_Search != rCellPopulation.End();
    //              ++cell_iter_Search)
    //         {
                
    //             unsigned node_index_Search = rCellPopulation.GetLocationIndexUsingCell(*cell_iter_Search);

    //             if (cell_iter_Search->GetCellData()->GetItem("Boundary") == 0 && node_index_Search != node_index)
    //             {
    //                 // Need location of each

    //                 Node<SPACE_DIM>* pNode1 = rCellPopulation.rGetMesh().GetNode(node_index_Search);

    //                 c_vector<long double, SPACE_DIM> LocationOfNode = pNode1->rGetLocation();
    //                 //  PRINT_VECTOR(LocationOfNode)

    //                 double Distance = norm_2(LocationOfNode - CellLocation);
        
    //             }
    //         }
    //         mNearestNodesMap[node_index] = NearestNodes;
    //     }
    // }
}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingModifier<ELEMENT_DIM, SPACE_DIM>::SetChasteOutputDirectory(std::string ChasteOutputDirectory, double startime)
{
    // Note std::string is easier to add together, but in the end what I will be putting into system will be a char * 
    char * C = getenv("CHASTE_TEST_OUTPUT");
    std::string directory;
    directory+=C;

    if (startime == (int) startime)
    {
        directory +="/"+ChasteOutputDirectory + "/results_from_time_"+std::to_string((int) startime)+ "/";
    }
    else
    {
        std::string TimeStamp = std::to_string(startime);
        TimeStamp.erase(TimeStamp.find_last_not_of('0') + 1, std::string::npos );
        std::cout<<  TimeStamp << std::endl;
        directory +="/"+ChasteOutputDirectory + "/results_from_time_"+TimeStamp+ "/";
    }
    mChasteOutputDirectory = directory;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingModifier<ELEMENT_DIM, SPACE_DIM>::RemeshGeometry()
{

    /*
	 * Here I will remesh the current chaste geometery
     * 1) Take the latest .vtu file convert it into an .stl 
     * 2) Remesh this stl
     * 3) Convert this stl into a .vtu so I can read it into chaste and use it  
	 */

      // Convert the .vtu into an .stl

      std::string stlfile = mChasteOutputDirectory + "CurrentMesh.stl";
      std::string vtu2stlCommand ="python projects/VascularRemodelling/apps/ConvertCurrentGeometryToSTL.py -ChasteOutput " +mChasteOutputDirectory +"  -stlOutput " + stlfile;
      std::system(vtu2stlCommand.c_str()); // system only takes char * 
      
      // Remesh the .stl
      std::string Remeshedstl =  mChasteOutputDirectory + "RemeshedGeometry.stl";
      std::string RemeshCommand ="vmtksurfaceremeshing -ifile " + stlfile+ " -iterations 10 -area 1e-7 -ofile " + Remeshedstl;
      std::system(RemeshCommand.c_str());


     //Finally conver the Remeshed stl into a vtu -- meshio is your friend
     std::string Remeshedvtu =  mChasteOutputDirectory + "RemeshedGeometry.vtu";
     std::string stl2vtuCommand =" meshio-convert " + Remeshedstl + " " +Remeshedvtu;
     std::system(stl2vtuCommand.c_str());

}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingModifier<ELEMENT_DIM, SPACE_DIM>::SaveInitalConditions(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
    // TRACE("Find the centroids for the original mesh")
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator iter = p_cell_population->rGetMesh().GetNodeIteratorBegin();
         iter != p_cell_population->rGetMesh().GetNodeIteratorBegin();
         ++iter)
    {
        unsigned Node_index =iter->GetIndex();
        mOriginalNodePositions[Node_index] = iter->rGetLocation();  
    }

}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::map<unsigned, c_vector<double, SPACE_DIM> > RemeshingModifier<ELEMENT_DIM, SPACE_DIM>::GetInitalConditions()
{
    return mInitalPositionOfRemeshedNodes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingModifier<ELEMENT_DIM, SPACE_DIM>::MappingAdaptedMeshToInitalGeometry(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, MutableMesh<2, 3>& rMesh )
{
    assert(SPACE_DIM == 3);  assert(ELEMENT_DIM == 2); 
    std::map<unsigned, c_vector<double, SPACE_DIM> > CentroidMap; // Save the centroid for each element
    TRACE("Now to do the mapping thing")



    
    MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
    // TRACE("Find the centroids for the original mesh")
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        Node<SPACE_DIM>* pNode0 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0));
        Node<SPACE_DIM>* pNode1 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(1));
        Node<SPACE_DIM>* pNode2 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));
        // Centroid is the average of the three triangle nodes 
        CentroidMap[elem_index] = (pNode0->rGetLocation() +pNode1->rGetLocation()+pNode2->rGetLocation())/3;        
    }

    // Read in the new mesh 
    VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mChasteOutputDirectory + "RemeshedGeometry.vtu");
    MutableMesh<ELEMENT_DIM, SPACE_DIM> New_mesh;
    New_mesh.ConstructFromMeshReader(mesh_reader);
    TRACE("Have the new mesh ");

    // Now I need to figure out how to iterate over the nodes -- the worst 
     for (typename MutableMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator iter = New_mesh.GetNodeIteratorBegin();
         iter != New_mesh.GetNodeIteratorEnd();
         ++iter)
    {
        unsigned node_index = iter->GetIndex();
        c_vector<double,SPACE_DIM> Location = iter->rGetLocation();
        double distance = 10000;
        unsigned ClosestElement;

        // Now I need to get the nearest centroid 
        for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = CentroidMap.begin();
            Centroid_iter != CentroidMap.end();
            ++Centroid_iter)
        {
            if ( std::abs(norm_2(Location - Centroid_iter->second )) < distance)
            {
                ClosestElement = Centroid_iter->first;
                distance  =  std::abs(norm_2(Location - Centroid_iter->second ));
            }
        }

        // SO i have looked over all the centorids and found the closest element 
        // PRINT_VARIABLE(ClosestElement)
        Element<ELEMENT_DIM,SPACE_DIM>* p_element = p_cell_population->rGetMesh().GetElement(ClosestElement);


        //Collect the intial configuration of the nodes, and the deformed configu, send them into the mapping function 
        c_vector<c_vector<double, SPACE_DIM> , SPACE_DIM> Element_0;
        c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> DeformedElement;

        for (int i =0; i<3; ++i)
        {
            Node<SPACE_DIM>* pNode = p_cell_population->rGetMesh().GetNode(p_element->GetNodeGlobalIndex(i));
            Element_0[i] = mOriginalNodePositions[pNode->GetIndex()]; 
            DeformedElement[i] = pNode->rGetLocation();
        }

        mInitalPositionOfRemeshedNodes[node_index]  = NewNodeInInitalConfigurationFromChangeOfBasis(Element_0, DeformedElement, Location);        
    }
 
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> RemeshingModifier<ELEMENT_DIM, SPACE_DIM>::NewNodeInInitalConfigurationFromChangeOfBasis(c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> Element_0, c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> DeformedElement, c_vector<double, SPACE_DIM> NewNode)
{
    // % Translate each to origin
    // % InitalConfiguration

    c_vector<double,SPACE_DIM> vector_12_0 = Element_0[1] - Element_0[0];
    c_vector<double,SPACE_DIM> vector_13_0 = Element_0[2] - Element_0[0];
    double a = norm_2(vector_12_0); //%% Lenght a -> edge connecting P1 and P2
    double b = norm_2(vector_13_0); //%% Lenght b -> edge connecting P1 and P3

    double alpha = acos(inner_prod(vector_12_0, vector_13_0) / (a*b));

    c_vector<double,SPACE_DIM> x1_0 = Create_c_vector(0,0,0);
    c_vector<double,SPACE_DIM> x2_0 = Create_c_vector(a,0,0);
    c_vector<double,SPACE_DIM> x3_0 = Create_c_vector(b * cos(alpha),b* sin(alpha),0);

    // Nodes at time t

    c_vector<double,SPACE_DIM> vector_12 = DeformedElement[1] - DeformedElement[0];
    c_vector<double,SPACE_DIM> vector_13 = DeformedElement[2] - DeformedElement[0];

    a = norm_2(vector_12);
    b = norm_2(vector_13);

    alpha = acos(inner_prod(vector_12, vector_13) / (a*b));

    // c_vector<double,SPACE_DIM>  x1 = Create_c_vector(0,0,0);
    // c_vector<double,SPACE_DIM>  x2 = Create_c_vector(a,0,0);
    // c_vector<double,SPACE_DIM>  x3 = Create_c_vector(b * cos(alpha),b* sin(alpha),0);
   
    // Need to put P into a corrdinate system with the basis vectors
    // vector_12, and vector_13. First need to translate by -Nodes(1,:)
    // beacuse all the points in this triangle have ( including the
    // basis vectors)

    c_vector<double,SPACE_DIM> z_basis = VectorProduct(vector_12, vector_13);
    c_vector<double,SPACE_DIM> P_translated = NewNode - DeformedElement[1];

    // Change of basis, simple algebra/rearanging gives us 
    //P_translated  == C1*vector_12 + C2*vector_13 + C3*z_basis 

    double x1 = vector_12[0]; double y1 = vector_12[1]; double z1 = vector_12[2];
    double x2 = vector_13[0]; double y2 = vector_13[1]; double z2 = vector_13[2];
    double x3 = z_basis[0];   double y3 = z_basis[1];   double z3 = z_basis[2];
    double p1 = P_translated[0]; double p2 = P_translated[1]; double p3 = P_translated[2]; 

    // Bunch of constants to make the similtaneous equations easier
    double E = 1/(y2*z3-y3*z2); double D = p2*z3-p3*y3; double F = z1*y3-y1*z3;
    double h = E*D*z2*x3/z3;
    double g = x1+F*E*x2-z1*x3/z3-E*F*z2*x3/z3;

    double C1=(p1-E*D*x2-p3*x3/z3-h)/g;
    double C2 = E*(D+a*F);
    double C3 =(p3-a*z1-b*z2)/z3;

    // so the point in the coordinate system described by the local element is
    // P_translated  == C1*vector_12 + C2*vector_13 + C3*z_basis 
    // Now take these basis vectors and replace them with the vector's original configuration 

    z_basis = VectorProduct(vector_12_0, vector_13_0);

    c_vector<double,SPACE_DIM> InitalPoint_Translated =  C1*vector_12_0 + C2*vector_13_0 + C3*z_basis;
    c_vector<double,SPACE_DIM> P_0 = InitalPoint_Translated + Element_0[1];

    return P_0;
}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{

    bool Trigger =0;
    if(Trigger == 1)
    {
        TRACE("Not Doing this yet");
    }
    // 
    // UpdateCellData( rCellPopulation);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{

    // assert(ELEMENT_DIM ==2);
    // TRACE("Switharoo")
    // assert(SPACE_DIM == 3);
    // // std::map<unsigned, c_vector<unsigned, 5> > mNearestNodesMap;

    //      for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
    //             cell_iter != rCellPopulation.End();
    //             ++cell_iter)
    //     {
    //         unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //         Node<SPACE_DIM>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);
    //         c_vector<unsigned, 5>  NearestNodes = mNearestNodesMap[node_index];
    //         // p_new_node->ClearAppliedForce();
    //         // pNode->AddAppliedForceContribution(MembraneForceMap[node_index] ); // Add the new force
    //        cell_iter->GetCellData()->SetItem("MembraneForce", -350 );

    //         //             Node<3>* pReferenceNode = p_cell_population->rGetMesh().GetNode(ReferenceNode);
    //       c_vector<long double, 3> ForceOnNode = pNode->rGetAppliedForce();
    //       pNode->ClearAppliedForce(); // remove any applied force, this stops there begin two expanding forces at this node
    //       pNode->AddAppliedForceContribution(Create_c_vector(-1,-1,0));
    //     }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class RemeshingModifier<2, 2>;
template class RemeshingModifier<2, 3>;
template class RemeshingModifier<3, 3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RemeshingModifier)

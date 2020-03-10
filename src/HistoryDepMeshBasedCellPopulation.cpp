
#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "VtkMeshWriter.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::HistoryDepMeshBasedCellPopulation(MutableMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                                                                             std::vector<CellPtr>& rCells,
                                                                                             const std::vector<unsigned> locationIndices,
                                                                                             bool deleteMesh,
                                                                                             bool validate)
        : MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>(rMesh, rCells, locationIndices, deleteMesh, validate)
{
    this->SaveInitalConditions();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::HistoryDepMeshBasedCellPopulation(MutableMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
        : MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>(rMesh)
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::~HistoryDepMeshBasedCellPopulation()
{
    TRACE("Is there an issue in my destructors -- HistoryDepMeshBasedCellPopulation")
    // mNew_mesh.~HistoryDepMutableMesh();
    // TRACE("Got rid of the new mesh class too ")
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SaveInitalConditions()
{
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        mOriginalNodePositions[node_index] = node_iter->rGetLocation();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetChasteOutputDirectory(std::string ChasteOutputDirectory, double startime)
{
    // Note std::string is easier to add together, but in the end what I will be putting into system will be a char *
    char* C = getenv("CHASTE_TEST_OUTPUT");
    std::string directory;
    directory += C;

    if (startime == (int)startime)
    {
        directory += "/" + ChasteOutputDirectory + "/results_from_time_" + std::to_string((int)startime) + "/";
    }
    else
    {
        std::string TimeStamp = std::to_string(startime);
        TimeStamp.erase(TimeStamp.find_last_not_of('0') + 1, std::string::npos);
        std::cout << TimeStamp << std::endl;
        directory += "/" + ChasteOutputDirectory + "/results_from_time_" + TimeStamp + "/";
    }
    mChasteOutputDirectory = directory;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::ExecuteHistoryDependentRemeshing()
{


    // VtkMeshReader<2,3> mesh_reader2(mChasteOutputDirectory + "RemeshedGeometry.vtu");
    // MutableMesh<2,3> New_mesh2;
    // New_mesh2.ConstructFromMeshReader(mesh_reader2);
    // New_mesh2.~MutableMesh();


    // step by step
    // 1) remesh geometry
    // 2) Map
    // 3) Delete old mesh and reset with new mesh
    // 4) detele old cells and relpace with the cells for the current mesh

    //  this->RemeshGeometry();
    //  this->MappingAdaptedMeshToInitalGeometry();

    // VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mChasteOutputDirectory + "RemeshedGeometry.vtu");
    // MutableMesh<ELEMENT_DIM, SPACE_DIM> New_mesh;
    // New_mesh.ConstructFromMeshReader(mesh_reader);
    // New_mesh.~MutableMesh();
    // TRACE("Distrucotr executed")

    VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mChasteOutputDirectory + "RemeshedGeometry.vtu");
    mNew_mesh.ConstructFromMeshReader(mesh_reader);

    static_cast<HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>&>((this->mrMesh)).DeleteMesh();
    // MutableMesh<ELEMENT_DIM, SPACE_DIM>* pNew_mesh = &mNew_mesh; 
    static_cast<HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>&>((this->mrMesh)).AssignNewMesh(&mNew_mesh);
   

//  TRACE("History dependent remeshing")

    // static_cast<HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>&>((this->mrMesh)).AddANewNodeBehindBoundary();
    //     VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer("CheckingremeshsedMesh", "config", false);
    // mesh_writer.WriteFilesUsingMesh(static_cast<AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>&>((this->mrMesh)));

//     // Now look at adding a new cell for the new node 
//     // Option 1, follow the line of constructors to fine where I need to clear/delete cells and add new cells on 

 

    TRACE("New nodes")
     MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

    std::vector<CellPtr> cells;
    CellsGenerator<FixedG1GenerationalCellCycleModel, ELEMENT_DIM> cells_generator;
    cells_generator.GenerateBasicRandom(cells, this->mrMesh.GetNumNodes(), p_differentiated_type);


//1)
    // From the  AbstractCellPopulation Constructor

     /** List of cells. */
    // std::list<CellPtr> mCells; -- Protected

    // std::list<CellPtr>& rGetCells(); @return reference to mCells.
     this->mCells.clear();
     PRINT_VARIABLE(this->mCells.size())
    for (std::vector<CellPtr>::iterator i=cells.begin(); i!= cells.end(); ++i)
    {
           this->mCells.push_back(*i);
    }
     /*
     * To avoid double-counting problems, clear the passed-in cells vector.
     * We force a reallocation of memory so that subsequent usage of the
     * vector is more likely to give an error.
     */
   std::vector<CellPtr>().swap(cells);
   

    // There must be a one-one correspondence between cells and location indices

    if (this->mCells.size() != this->mrMesh.GetNumNodes())
    {
        EXCEPTION("There is not a one-one correspondence between cells and location indices");
    }
    // Set up the map between location indices and cells
    this->mLocationCellMap.clear();
    this->mCellLocationMap.clear();


    std::list<CellPtr>::iterator it = this->mCells.begin();
    for (unsigned i=0; it != this->mCells.end(); ++it, ++i)
    {
        // These are cell things
        // Give each cell a pointer to the property registry (we have taken ownership in this constructor)

        (*it)->rGetCellPropertyCollection().SetCellPropertyRegistry(this->mpCellPropertyRegistry.get());
    }
    
// ////-----------------------------------------
// 2)  From the AbstractOffLatticeCellPopulation Constructor 
// 3)  From the AbstractCentreBasedCellPopulation Constructor 

// If no location indices are specified, associate with nodes from the mesh.
    std::list<CellPtr>::iterator iter = this->mCells.begin();
    PRINT_2_VARIABLES(this->mCells.size() , this->mrMesh.GetNumNodes())
    typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();

    for (unsigned i=0; iter != this->mCells.end(); ++iter, ++i, ++node_iter)
    {
        unsigned index =  node_iter->GetIndex(); // assume that the ordering matches
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::AddCellUsingLocationIndex(index,*iter);
        
    }

    this->mpCentreBasedDivisionRule.reset(new RandomDirectionCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>());


// 4)
// From the MeshBasedCellPopultion Constructor 
    PRINT_2_VARIABLES(this->mCells.size(), this->mrMesh.GetNumNodes())
    assert(this->mCells.size() == this->mrMesh.GetNumNodes());
    bool validate =1;
    if (validate)
    {
      
        this->Validate();
    }

    // Initialise the applied force at each node to zero
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }

// //  this->Update(0);
// //  this->InitialiseCells();
    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::RemeshGeometry()
{

    /*
	 * Here I will remesh the current chaste geometery
     * 1) Take the latest .vtu file convert it into an .stl 
     * 2) Remesh this stl
     * 3) Convert this stl into a .vtu so I can read it into chaste and use it  
	 */

    // Convert the .vtu into an .stl

    std::string stlfile = mChasteOutputDirectory + "CurrentMesh.stl";
    std::string vtu2stlCommand = "python projects/VascularRemodelling/apps/ConvertCurrentGeometryToSTL.py -ChasteOutput " + mChasteOutputDirectory + "  -stlOutput " + stlfile;
    std::system(vtu2stlCommand.c_str()); // system only takes char *

    // Remesh the .stl
    std::string Remeshedstl = mChasteOutputDirectory + "RemeshedGeometry.stl";
    std::string RemeshCommand = "vmtksurfaceremeshing -ifile " + stlfile + " -iterations 10 -area 1e-7 -ofile " + Remeshedstl;
    std::system(RemeshCommand.c_str());

    //Finally conver the Remeshed stl into a vtu -- meshio is your friend
    std::string Remeshedvtu = mChasteOutputDirectory + "RemeshedGeometry.vtu";
    std::string stl2vtuCommand = " meshio-convert " + Remeshedstl + " " + Remeshedvtu;
    std::system(stl2vtuCommand.c_str());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::MappingAdaptedMeshToInitalGeometry()
{
    assert(SPACE_DIM == 3);
    assert(ELEMENT_DIM == 2);
    std::map<unsigned, c_vector<double, SPACE_DIM> > CentroidMap; // Save the centroid for each element
    TRACE("Now to do the mapping thing")

    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
         elem_iter != this->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        Node<SPACE_DIM>* pNode0 = this->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0));
        Node<SPACE_DIM>* pNode1 = this->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(1));
        Node<SPACE_DIM>* pNode2 = this->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));
        CentroidMap[elem_index] = (pNode0->rGetLocation() + pNode1->rGetLocation() + pNode2->rGetLocation()) / 3;
    }

    // Read in the new mesh
    VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mChasteOutputDirectory + "RemeshedGeometry.vtu");
    mNew_mesh.ConstructFromMeshReader(mesh_reader);
    TRACE("Have the new mesh ;) ");

    // Now I need to figure out how to iterate over the nodes -- the worst
    for (typename MutableMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator iter = mNew_mesh.GetNodeIteratorBegin();
         iter != mNew_mesh.GetNodeIteratorEnd();
         ++iter)
    {
        unsigned node_index = iter->GetIndex();
        c_vector<double, SPACE_DIM> NewNodeLocation = iter->rGetLocation();
        double distance = 10000;
        unsigned ClosestElement;

        // Now I need to get the nearest centroid
        for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = CentroidMap.begin();
             Centroid_iter != CentroidMap.end();
             ++Centroid_iter)
        {
            if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance)
            {
                ClosestElement = Centroid_iter->first;
                distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
            }
        }

        // SO i have looked over all the centorids and found the closest element
        // PRINT_VARIABLE(ClosestElement)
        Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->rGetMesh().GetElement(ClosestElement);
        mInitalPositionOfRemeshedNodes[node_index] = NewNodeInInitalConfigurationFromChangeOfBasis(ClosestElement, NewNodeLocation);
    }
}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::NewNodeInInitalConfigurationFromChangeOfBasis(unsigned ClosestElement_OldMeshIndex, c_vector<double, SPACE_DIM> NewNode)
{
    Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->rGetMesh().GetElement(ClosestElement_OldMeshIndex);

    // //Collect the intial configuration of the nodes, and the deformed configu, send them into the mapping function
    c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> Element_0;
    c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> DeformedElement;

    for (int i = 0; i < 3; ++i)
    {
        Node<SPACE_DIM>* pNode = this->rGetMesh().GetNode(p_element->GetNodeGlobalIndex(i));
        Element_0[i] = mOriginalNodePositions[pNode->GetIndex()];
        DeformedElement[i] = pNode->rGetLocation();
    }

    c_vector<double, SPACE_DIM> vector_12_0 = Element_0[1] - Element_0[0];
    c_vector<double, SPACE_DIM> vector_13_0 = Element_0[2] - Element_0[0];
    double a = norm_2(vector_12_0); //%% Lenght a -> edge connecting P1 and P2
    double b = norm_2(vector_13_0); //%% Lenght b -> edge connecting P1 and P3

    double alpha = acos(inner_prod(vector_12_0, vector_13_0) / (a * b));

    c_vector<double, SPACE_DIM> x1_0 = Create_c_vector(0, 0, 0);
    c_vector<double, SPACE_DIM> x2_0 = Create_c_vector(a, 0, 0);
    c_vector<double, SPACE_DIM> x3_0 = Create_c_vector(b * cos(alpha), b * sin(alpha), 0);

    // Nodes at time t

    c_vector<double, SPACE_DIM> vector_12 = DeformedElement[1] - DeformedElement[0];
    c_vector<double, SPACE_DIM> vector_13 = DeformedElement[2] - DeformedElement[0];

    a = norm_2(vector_12);
    b = norm_2(vector_13);

    alpha = acos(inner_prod(vector_12, vector_13) / (a * b));

    // Need to put P into a corrdinate system with the basis vectors
    // vector_12, and vector_13. First need to translate by -Nodes(1,:)
    // beacuse all the points in this triangle have ( including the
    // basis vectors)

    c_vector<double, SPACE_DIM> z_basis = VectorProduct(vector_12, vector_13);
    c_vector<double, SPACE_DIM> P_translated = NewNode - DeformedElement[1];

    // Change of basis, simple algebra/rearanging gives us
    //P_translated  == C1*vector_12 + C2*vector_13 + C3*z_basis

    double x1 = vector_12[0];
    double y1 = vector_12[1];
    double z1 = vector_12[2];
    double x2 = vector_13[0];
    double y2 = vector_13[1];
    double z2 = vector_13[2];
    double x3 = z_basis[0];
    double y3 = z_basis[1];
    double z3 = z_basis[2];
    double p1 = P_translated[0];
    double p2 = P_translated[1];
    double p3 = P_translated[2];

    // Bunch of constants to make the similtaneous equations easier
    double E = 1 / (y2 * z3 - y3 * z2);
    double D = p2 * z3 - p3 * y3;
    double F = z1 * y3 - y1 * z3;
    double h = E * D * z2 * x3 / z3;
    double g = x1 + F * E * x2 - z1 * x3 / z3 - E * F * z2 * x3 / z3;

    double C1 = (p1 - E * D * x2 - p3 * x3 / z3 - h) / g;
    double C2 = E * (D + a * F);
    double C3 = (p3 - a * z1 - b * z2) / z3;

    // so the point in the coordinate system described by the local element is
    // P_translated  == C1*vector_12 + C2*vector_13 + C3*z_basis
    // Now take these basis vectors and replace them with the vector's original configuration

    z_basis = VectorProduct(vector_12_0, vector_13_0);

    c_vector<double, SPACE_DIM> InitalPoint_Translated = C1 * vector_12_0 + C2 * vector_13_0 + C3 * z_basis;
    c_vector<double, SPACE_DIM> P_0 = InitalPoint_Translated + Element_0[1];

    return P_0;
}

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::DeleteCellPopulation()
// {
//    // Need to be able to delete the cells and mesh
// }
// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::DeleteMesh(out_stream& rParamsFile)
// {
//    // Need to be able to delete the cells and mesh
// }

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationParameters(rParamsFile);
}

// Explicit instantiation
// template class HistoryDepMeshBasedCellPopulation<2, 2>;
template class HistoryDepMeshBasedCellPopulation<2, 3>;
// template class HistoryDepMeshBasedCellPopulation<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(HistoryDepMeshBasedCellPopulation)

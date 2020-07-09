
#include "HistoryDepMeshBasedCellPopulation.hpp"
#include <ctime>
#include "MathsFunctions.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"
#include "VtkMeshWriter.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::HistoryDepMeshBasedCellPopulation(MutableMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                                                                             std::vector<CellPtr>& rCells,
                                                                                             const std::vector<unsigned> locationIndices,
                                                                                             bool deleteMesh,
                                                                                             bool validate)
        : MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>(rMesh, rCells, locationIndices, deleteMesh, validate)
{
    // TRACE("SAVE initial config")
    SaveInitalConditions();
    // TRACE("SAVE initial config at Time 0 ")
    SaveInitalConditionsAtTime0();
    // TRACE("Mark Boudary ")
    DetermineDimensions();
    SetBinningIntervals(1,1,1);
    mFudgeScalling = 0.05;
    SetBinningRegions2();


    bool SetUpInitialConfigurations =0;
    if (SetUpInitialConfigurations ==1)
    {
        SaveInitalConditions();
        SaveInitalConditionsAtTime0();
        MarkBoundaryNodes();
        // // Define all the inital configurations ... this needed to be done and the start, and whenever there is a remeshing
        SetupMembraneConfiguration();
        // SetInitialAnlgesAcrossMembrane();
        SetMaxEdgelength();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::HistoryDepMeshBasedCellPopulation(MutableMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
        : MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>(rMesh)
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::~HistoryDepMeshBasedCellPopulation()
{

    // mNew_mesh.~HistoryDepMutableMesh();
    // TRACE("Got rid of the new mesh class too ")
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::ExecuteHistoryDependentRemeshing()
{
    TRACE("ExecuteHistoryDependentRemeshing")
    // step by step
    // 1) remesh geometry
    // 2) Map
    // 3) Delete old mesh and reset with new mesh
    // 4) detele old cells and relpace with the cells for the current mesh

    if (mRemeshingSoftwear == "VMTK")
    {
        this->RemeshGeometryWithVMTK();
    }
    else if (mRemeshingSoftwear == "CGAL")
    {
        this->RemeshGeometry();
    }else if 
    (mRemeshingSoftwear == "PreAllocatedMatlabMesh")
    {
        this->TakeInPreAllocatedRemeshGeometry();
    }
    this->SetBinningRegions2();
    this->MappingAdaptedMeshToInitalGeometry();
    if (mPrintRemeshedIC)
    {
        WriteOutMappedInitalConfig();
    }

    /*
        --------------------------------
              The Old switcharoo 
        
        Now Switch the old and new meshes
        ---------------------------------
     */

    // Delete all the nodes and elements in the old mesh
    static_cast<HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>&>((this->mrMesh)).DeleteMesh();
    // Assign the nodes and elements of the new mesh to the simulation's mesh
    static_cast<HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>&>((this->mrMesh)).AssignNewMesh(&mNew_mesh);

    MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
    std::vector<CellPtr> cells;
    CellsGenerator<FixedG1GenerationalCellCycleModel, ELEMENT_DIM> cells_generator;
    cells_generator.GenerateBasicRandom(cells, this->mrMesh.GetNumNodes(), p_differentiated_type);

    // 1) From the  AbstractCellPopulation Constructor  std::list<CellPtr> mCells; -- Protected
    this->mCells.clear();
    for (std::vector<CellPtr>::iterator i = cells.begin(); i != cells.end(); ++i)
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
    for (unsigned i = 0; it != this->mCells.end(); ++it, ++i)
    {
        // Give each cell a pointer to the property registry (we have taken ownership in this constructor)
        (*it)->rGetCellPropertyCollection().SetCellPropertyRegistry(this->mpCellPropertyRegistry.get());
    }

    // 3)  From the AbstractCentreBasedCellPopulation Constructor -- If no location indices are specified, associate with nodes from the mesh.

    std::list<CellPtr>::iterator iter = this->mCells.begin();
    typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();

    for (unsigned i = 0; iter != this->mCells.end(); ++iter, ++i, ++node_iter)
    {
        unsigned index = node_iter->GetIndex(); // assume that the ordering matches
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::AddCellUsingLocationIndex(index, *iter);
    }

    this->mpCentreBasedDivisionRule.reset(new RandomDirectionCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>());

    // 4) From the MeshBasedCellPopultion Constructor

    assert(this->mCells.size() == this->mrMesh.GetNumNodes());
    this->Validate();

    // Now redefine all the inital configurations ... this needed to be done and the start, and whenever there is a remeshing
    MarkBoundaryNodes();
    SetupMembraneConfiguration();
    // SetInitialAnlgesAcrossMembrane();

    // Initialise the applied force at each node to zero
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        node_iter->ClearAppliedForce();
        CellPtr p_cell = this->GetCellUsingLocationIndex(node_iter->GetIndex());
        p_cell->GetCellData()->SetItem("MappingMethod", mMapOfProbNodes[node_iter->GetIndex()]);
    }

    this->SetBinningRegions2();
    // TRACE("Binning dome")

    // I want to iterate around and set the problem nodes from the mapping as cell data
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::WriteOutMappedInitalConfig()
{
    // Read in the new mesh and have it as its own thjing
    //TRACE("Try write out the inital mesh");

    // Read in the new mesh
    VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mChasteOutputDirectory + "RemeshedGeometry.vtu");

    MutableMesh<ELEMENT_DIM, SPACE_DIM> MappedICmesh;

    MappedICmesh.ConstructFromMeshReader(mesh_reader);

    //TRACE("Have the adapted mesh, which now I can change to the mapped inital mesh and write out to see what is going on ;) ");
    //PRINT_VARIABLE(MappedICmesh.GetNumNodes());

    for (int i = 0; i < MappedICmesh.GetNumNodes(); i++)
    {
        Node<SPACE_DIM>* p_node = MappedICmesh.GetNode(i);
        p_node->rGetModifiableLocation() = mInitalPositionOfRemeshedNodes[i];
    }
    // TRACE("Adapted nodes, now write");
    std::string OutputFile = "NewInitalConfiguration" + std::to_string(mNumberOfChanges);
    VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(mRelativePath, OutputFile, false);
    // TRACE("Creater writer");

    mesh_writer.WriteFilesUsingMesh(MappedICmesh);
    // TRACE("Written")
    mNumberOfChanges += 1;
    // MappedICmesh.~MutableMesh();
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::RemeshGeometry()
{

    // TRACE("Remesh geometry here ")
    // Outpust the mesh in .vtu format for HemeLB setup tool to pick up (first converted to stl, though).
    // std::string stlfile = mChasteOutputDirectory;
    VtkMeshWriter<2, 3> mesh_writer(mRelativePath, "config", false);
    mesh_writer.WriteFilesUsingMesh(this->rGetMesh());

    // For using .off and CGAL
    /*
	 * Here I will remesh the current chaste geometery
     * 1) Take the latest .vtu file convert it into an .off
     * 2) Remesh this stl
     * 3) Convert this stl into a .vtu so I can read it into chaste and use it  
	 */

    // Convert the .vtu into an .off
    std::string offfile = mChasteOutputDirectory + "CurrentMesh.off";
    // TRACE("Create .off")
    PRINT_VARIABLE(mChasteOutputDirectory)
    std::string vtu2offCommand = "meshio-convert " + mChasteOutputDirectory + "config.vtu " + offfile;
    std::system(vtu2offCommand.c_str()); // system only takes char *.

    // Now excute the CGAL command to remesh the current geometry - not the input and output within this file have to be pre-set. I will explore if I can make this more neat later, should care.... dont care
    std::string CGALRemeshingCommand = "(cd  ~/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/;./isotropic_remeshing_ForChaste -input " + offfile + " -output " + mChasteOutputDirectory + "CurrentPlexusRemeshed.off -target_edge_length " + std::to_string(mTargetRemeshingEdgeLength) + " -iterations " + std::to_string(mIterations) + " )";
    std::system(CGALRemeshingCommand.c_str()); // system only takes char *

    // Now ned to convert from .off back to a .vtu

    std::string Remeshedvtu = mChasteOutputDirectory + "RemeshedGeometry.vtu";
    std::string off2vtuCommand = "meshio-convert " + mChasteOutputDirectory + "CurrentPlexusRemeshed.off " + Remeshedvtu;
    std::system(off2vtuCommand.c_str()); // system only takes char *

    // /// Delete the orignal .off file -- else if the CurrentMesh.off fails to transfer, the old transfer will be the one taken
    // std::string DeleteOldCurrentConfigOff = "rm ~/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/data/CurrentMesh.off";
    // std::system( DeleteOldCurrentConfigOff.c_str()); // system only takes char *

    // Read in the new remeshsed mesh
    VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mChasteOutputDirectory + "RemeshedGeometry.vtu");
    mNew_mesh.ConstructFromMeshReader(mesh_reader);
    // TRACE("Have the new mesh ;) ")
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::RemeshGeometryWithVMTK()
{
    // old, when using vmtk and stls
    /*
    * Here I will remesh the current chaste geometery
    * 1) Take the latest .vtu file convert it into an .stl 
    * 2) Remesh this stl
    * 3) Convert this stl into a .vtu so I can read it into chaste and use it  
    */

    // Outpust the mesh in .vtu format for HemeLB setup tool to pick up (first converted to stl, though).
    // std::string stlfile = mChasteOutputDirectory;
    VtkMeshWriter<2, 3> mesh_writer(mRelativePath, "config", false);
    mesh_writer.WriteFilesUsingMesh(this->rGetMesh());

    // Convert the .vtu into an .off
    std::string stlfile = mChasteOutputDirectory + "CurrentMesh.stl";
    // TRACE("Create .stl")
    std::string vtu2stlCommand = "meshio-convert " + mChasteOutputDirectory + "config.vtu " + stlfile;
    std::system(vtu2stlCommand.c_str()); // system only takes char *.

    // Remesh the .stl
    std::string Remeshedstl = mChasteOutputDirectory + "RemeshedGeometry.stl";
    // TRACE("Remeshing geometry")

    // More info for the remeshing options found at https://sourceforge.net/p/vmtk/mailman/message/29391820/ and  http://www.vmtk.org/vmtkscripts/vmtksurfaceremeshing.html
    // std::string RemeshCommand = "vmtksurfaceremeshing -ifile " + stlfile + " -iterations " + std::to_string(mIterations) + " -area " + std::to_string(mTargetRemeshingElementArea) + " -maxarea "+ std::to_string(1.2*mTargetRemeshingElementArea) +" -ofile " + Remeshedstl;
    // std::system(RemeshCommand.c_str());
    // PRINT_VARIABLE(mIterations)
    std::string RemeshCommand = "vmtksurfaceremeshing -ifile " + stlfile + " -iterations " + std::to_string(mIterations) + " -elementsizemode 'edgelength' -edgelength " + std::to_string(mTargetRemeshingEdgeLength) + " -ofile " + Remeshedstl + " -elementsizemode edgelength";
    std::system(RemeshCommand.c_str());

    //Finally conver the Remeshed stl into a vtu -- meshio is your friend
    std::string Remeshedvtu = mChasteOutputDirectory + "RemeshedGeometry.vtu";
    std::string stl2vtuCommand = " meshio-convert " + Remeshedstl + " " + Remeshedvtu;
    std::system(stl2vtuCommand.c_str());

    // Read in the new mesh
    VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mChasteOutputDirectory + "RemeshedGeometry.vtu");
    mNew_mesh.ConstructFromMeshReader(mesh_reader);
    // TRACE("Have the new mesh ;) ");
}

// Need to load up the mesh 
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::TakeInPreAllocatedRemeshGeometry()
{
    std::string Remeshedvtu = mChasteOutputDirectory + "RemeshedGeometry.vtu";
    PRINT_VARIABLE(mChasteOutputDirectory)
    copyfile(mPreAllocatedRemeshedMesh.c_str(), Remeshedvtu.c_str(), NULL, COPYFILE_DATA | COPYFILE_XATTR);
     TRACE("Read in the new mesh")
    VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mPreAllocatedRemeshedMesh);
    mNew_mesh.ConstructFromMeshReader(mesh_reader);
    // TRACE("Have the new mesh ;) ");
}





template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetCentroidMap()
{
    mCentroidMap.clear();
    assert(SPACE_DIM == 3);
    assert(ELEMENT_DIM == 2);
    // Loop over the old map and get the centroids of the old map
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
         elem_iter != this->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        Node<SPACE_DIM>* pNode0 = this->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0));
        Node<SPACE_DIM>* pNode1 = this->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(1));
        Node<SPACE_DIM>* pNode2 = this->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));
        mCentroidMap[elem_index] = (pNode0->rGetLocation() + pNode1->rGetLocation() + pNode2->rGetLocation()) / 3;
    }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::DetermineDimensions()
{
    mDIM=2;
    double Z = 0;//mNew_mesh.GetNodeIteratorBegin()->rGetLocation()[2];
    for (typename MutableMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator iter = mNew_mesh.GetNodeIteratorBegin();
         iter != mNew_mesh.GetNodeIteratorEnd();
         ++iter)
    {
        unsigned node_index = iter->GetIndex();
        c_vector<double, SPACE_DIM> NewNodeLocation = iter->rGetLocation();
        if (NewNodeLocation[2]!=Z)
        {
            mDIM =3;
            break;
        }
    }
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetCommonNodes(unsigned elem_1, unsigned elem_2)
{
    std::vector<unsigned> NodeSet1Vec;
    std::vector<unsigned> NodeSet2Vec;

    std::vector<unsigned> CommonNodes;
    Element<ELEMENT_DIM, SPACE_DIM>* pElement1 = this->rGetMesh().GetElement(elem_1);
    Element<ELEMENT_DIM, SPACE_DIM>* pElement2 = this->rGetMesh().GetElement(elem_2);

    // Loop over the old map and get the centroids of the old map
    for (int i = 0; i < 3; ++i)
    {
        Node<SPACE_DIM>* pNode1 = this->rGetMesh().GetNode(pElement1->GetNodeGlobalIndex(i));
        Node<SPACE_DIM>* pNode2 = this->rGetMesh().GetNode(pElement2->GetNodeGlobalIndex(i));
        unsigned node_index1 = pNode1->GetIndex();
        unsigned node_index2 = pNode2->GetIndex();
        NodeSet1Vec.push_back(node_index1);
        NodeSet2Vec.push_back(node_index2);
    }

    std::set_intersection(NodeSet1Vec.begin(), NodeSet1Vec.end(),
                          NodeSet2Vec.begin(), NodeSet2Vec.end(),
                          back_inserter(CommonNodes));

    return CommonNodes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<int> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::DetermineBin(c_vector<double, SPACE_DIM> Point)
{
    bool FoundXRegion = 0;
    bool FoundYRegion = 0;
    bool FoundZRegion = 0;
    int xStep = 1;
    int yStep = 1;
    int zStep = 1;
    while (FoundXRegion == 0 && xStep <= mNx)
    {  double  itervalx = (mMaxX - mMinX)/mNx;
        if (Point[0] <mMinX + itervalx * xStep) // Look for x bin
        {

        // if (Point[0] - mMinX <= (mMaxX - mMinX) * xStep / mNx) // Look for x bin
        // {
            FoundXRegion = 1; // SO we know it should be in this spacial step
            while (FoundYRegion == 0 && yStep <= mNy)
            {
                if (Point[1] <= (mMaxY - mMinY) * yStep / mNy) // Look for y bin
                {
                    
                // if (Point[1] - mMinY <= (mMaxY - mMinY) * yStep / mNy) // Look for y bin
                // {
                    FoundYRegion = 1; // We now also have the y bin
                    while (FoundZRegion == 0 && zStep <= mNz) // Look for z bin
                    {
                        if (Point[2] <= (mMaxZ - mMinZ) * zStep / mNz)
                        {

                        // if (Point[2] - mMinZ <= (mMaxZ - mMinZ) * zStep / mNz)
                        // {
                            FoundZRegion = 1; // We now also have the y bin
                        }
                        else
                        {
                            zStep += 1;
                        }
                    }
                }
                else
                {
                    yStep += 1;
                }
            }
        }
        else
        {
            xStep += 1;
        }
    }

    xStep = std::min(xStep, mNx); // if there are needed my discrtisation is probably crap
    yStep = std::min(yStep, mNy);
    zStep = std::min(zStep, mNz);
    assert(xStep <= mNx && yStep <= mNy & zStep <= mNz); // If this is called, a bin has,t been found, and I need to figure out why
    std::vector<int> Bin ={xStep, yStep, zStep};
    // c_vector<double,3> Bin =Create_c_vector(xStep, yStep, zStep);
    return Bin;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringElements(unsigned elem_index)
{
    Element<ELEMENT_DIM, SPACE_DIM>* pElement = this->rGetMesh().GetElement(elem_index);

    std::set<unsigned> neighbouring_elements_indices;
    std::set<unsigned> neighbouring_elements;

    // Form a set of neighbouring elements via the nodes
    for (unsigned i = 0; i < 3; i++)
    {
        Node<SPACE_DIM>* p_node = pElement->GetNode(i);
        neighbouring_elements_indices = p_node->rGetContainingElementIndices();

        std::set_union(neighbouring_elements.begin(), neighbouring_elements.end(),
                       neighbouring_elements_indices.begin(), neighbouring_elements_indices.end(),
                       std::inserter(neighbouring_elements, neighbouring_elements.begin()));
    }

    return neighbouring_elements;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetElementNormal(unsigned elem_index)
{
    assert(SPACE_DIM == 3);
    Element<ELEMENT_DIM, SPACE_DIM>* pElement = this->rGetMesh().GetElement(elem_index);
    Node<SPACE_DIM>* p_node1 = pElement->GetNode(0);
    Node<SPACE_DIM>* p_node2 = pElement->GetNode(1);
    Node<SPACE_DIM>* p_node3 = pElement->GetNode(2);

    c_vector<double, SPACE_DIM> Vector1 = p_node1->rGetLocation();
    c_vector<double, SPACE_DIM> Vector2 = p_node2->rGetLocation();
    c_vector<double, SPACE_DIM> Vector3 = p_node3->rGetLocation();

    c_vector<double, SPACE_DIM> Vector12 = Vector2 - Vector1;
    c_vector<double, SPACE_DIM> Vector13 = Vector3 - Vector1;

    c_vector<double, SPACE_DIM> Normal = VectorProduct(Vector13, Vector12);
    Normal /= norm_2(Normal);

    return Normal;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::MappingAdaptedMeshToInitalGeometry()
{
    std::map<unsigned, c_vector<double, SPACE_DIM> > InitalPositionOfRemeshedNodes;
    SetCentroidMap();
    assert(SPACE_DIM == 3);
    assert(ELEMENT_DIM == 2);

    // Now iterate over the nodes of the new mesh, and find the cloesest element from the old mesh -
    // To accuratly determine the closest, I need to find the two closest, and then determine which one it is closer to. Sometimes might not be super clear

    for (typename MutableMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator iter = mNew_mesh.GetNodeIteratorBegin();
         iter != mNew_mesh.GetNodeIteratorEnd();
         ++iter)
    {
        double distance = 1e10;
        unsigned node_index = iter->GetIndex();
        c_vector<double, SPACE_DIM> NewNodeLocation = iter->rGetLocation();

        // Now I need to get the nearest centroid
        // double ClosestElement =-10;
        // for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
        //      Centroid_iter != mCentroidMap.end();
        //      ++Centroid_iter)
        // {
        //     if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance)
        //     {
        //         ClosestElement = Centroid_iter->first;
        //         distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
        //     }
        // }

        //  unsigned num_elements = this->rGetMesh().GetNumElements();
         
        // // PRINT_3_VARIABLES(ClosestElement,num_elements, mCentroidMap.size());
        // assert(ClosestElement <= num_elements);
        // assert( ClosestElement != -10 );



        unsigned ClosestElement =  GetClosestElementInOldMeshMethod2( node_index,   NewNodeLocation);
        InitalPositionOfRemeshedNodes[node_index] = NewNodeInInitalConfigurationFromChangeOfBasis(ClosestElement, NewNodeLocation);
        mInitalPositionOfRemeshedNodes = InitalPositionOfRemeshedNodes;
    }

    mOriginalNodePositions = InitalPositionOfRemeshedNodes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetClosestElementInOldMeshMethod2(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation)
{

    // This method is super simple. -- Just find the closest element -- it isnt perfect, 
    unsigned ClosestElement;
    double distance = 100;
    bool HaveElement =0;
    
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
            elem_iter != this->rGetMesh().GetElementIteratorEnd();
            ++elem_iter)
    {
        // Now need to find the cloeset element in this bin
        unsigned elem_index = elem_iter->GetIndex();
        c_vector<double, SPACE_DIM> Centroid = mCentroidMap[elem_index];

        if (abs(norm_2(NewNodeLocation - Centroid)) <= distance)
        {
            distance = abs(norm_2(NewNodeLocation - Centroid));
            ClosestElement = elem_index;  
            // if (PointInTriangle3D(NewNodeLocation, ClosestElement))
            //     {
            //         break;
            //     } 
        }
    }
//   HaveElement =1;
    if (PointInTriangle3D(NewNodeLocation, ClosestElement))
    {
      HaveElement =1;
      mMapOfProbNodes[node_index] = 1; 
    }
    
    unsigned num_elements = this->rGetMesh().GetNumElements();
/// ----------------------- ////////
    if (HaveElement == 0)
    {
        unsigned OrginalClosest = ClosestElement;
        distance = 100;
        std::pair<unsigned,unsigned> ClosestEdge;
        c_vector<double, SPACE_DIM> ClosestEdgeNorm ;
        c_vector<double, SPACE_DIM> ClosestMidPoint;
         for (typename MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SpringIterator spring_iterator = this->SpringsBegin();
            spring_iterator != this->SpringsEnd();
            ++spring_iterator)
            {
                unsigned NodeIndexA = spring_iterator.GetNodeA()->GetIndex();
                unsigned NodeIndexB = spring_iterator.GetNodeB()->GetIndex();

                c_vector<double, SPACE_DIM> P1 = this->GetNode(NodeIndexA)->rGetLocation();
                c_vector<double, SPACE_DIM> P2 = this->GetNode(NodeIndexB)->rGetLocation();
                c_vector<double, SPACE_DIM> EdgeVector = P2- P1;
                EdgeVector/=norm_2(EdgeVector);
                c_vector<double, SPACE_DIM> V = NewNodeLocation - P1;
                 // Normal to the plane containing the edge and the new node
                c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(EdgeVector, V)/norm_2(V);

                c_vector<double, SPACE_DIM> Norm;

                PlaneNormal/=norm_2(PlaneNormal);
            
            
                // Now can find the normal to the edge that is contaned in the plane defined by the three points
                c_vector<double, SPACE_DIM> EdgeNormal = VectorProduct(EdgeVector, PlaneNormal);
                EdgeNormal/=norm_2(EdgeNormal);
                double Dist2Edge = abs(inner_prod(EdgeNormal, V)); 
        
                if ( Dist2Edge<= distance )
                {
                    distance = Dist2Edge;
                    ClosestEdge = std::pair<unsigned,unsigned>(NodeIndexA, NodeIndexB);
                    ClosestEdgeNorm = EdgeNormal;
                }
            }
        c_vector<double, SPACE_DIM> P1 = this->GetNode(ClosestEdge.first)->rGetLocation();
        c_vector<double, SPACE_DIM> P2 = this->GetNode(ClosestEdge.second)->rGetLocation();

        std::set<unsigned> elements_containing_node1 = this->GetNode(ClosestEdge.first)->rGetContainingElementIndices();
        std::set<unsigned> elements_containing_node3 = this->GetNode(ClosestEdge.second)->rGetContainingElementIndices();

        std::set<unsigned> shared_elements;
        std::set_intersection(elements_containing_node1.begin(),  elements_containing_node1.end(), elements_containing_node3.begin(),  elements_containing_node3.end(),  std::inserter(shared_elements, shared_elements.begin()));
        typename std::set<unsigned>::iterator CommonElements = shared_elements.begin();

        unsigned Element1 = *CommonElements;
        if (shared_elements.size() ==1)
        {
                ClosestElement =Element1; HaveElement=1;
                mMapOfProbNodes[node_index] = 2; 
        }
        else 
        {   
            std::advance(CommonElements, 1); 
            unsigned Element2 = *CommonElements;
            if (PointInTriangle3D(NewNodeLocation, Element1)==1)
            {
                ClosestElement =Element1; HaveElement=1;
                mMapOfProbNodes[node_index] = 3; 
                TRACE("Option 1")
            }else if( PointInTriangle3D(NewNodeLocation, Element2) ==1)
            {
                ClosestElement =Element2;  HaveElement=1;

                mMapOfProbNodes[node_index] = 3; 
                TRACE("Option 2")
            }
        }
        if (HaveElement ==0)
        {
            distance =100;
            ClosestElement = OrginalClosest;
            // for (typename MutableMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator iter = mNew_mesh.GetNodeIteratorBegin();
            // iter != mNew_mesh.GetNodeIteratorEnd();
            // ++iter)
            // {
            //     unsigned node_index = iter->GetIndex();
            //     c_vector<double, SPACE_DIM> NodeLocation = iter->rGetLocation();
            //     double DistToNode = norm_2(NodeLocation  - NewNodeLocation);

            //     if (DistToNode < distance)
            //     {
            //         distance = DistToNode;
            //         ClosestNode = node_index;
            //     }
            // }
            // PRINT_VARIABLE(ClosestNode)
            //  std::set<unsigned> ContainingElements = this->GetNode(ClosestNode)->rGetContainingElementIndices();
            //  -- need to iterate over the containing elements and get the one containing 
            //  for (std::set<unsigned>::iterator iter = ContainingElements.begin();
                //  iter != ContainingElements.end();
                //  ++iter)
                // {

                //     if (PointInTriangle3D(NewNodeLocation, *iter))
                //     {
                //             ClosestElement = *iter;
                //     }
                // }
        }


    }

    double ClosestCentroidDistance = abs(norm_2(NewNodeLocation - mCentroidMap[ClosestElement]));
    // PRINT_VARIABLE(ClosestCentroidDistance);
    if(ClosestCentroidDistance >0.85)
    {
        TRACE("Broken node here")
    }
// check closest point is within certain region 


    assert(ClosestElement !=-10);
    // mMapOfProbNodes[node_index] = 3;
    return ClosestElement;
}

// else if ( ClosestPointInTriangle(NewNodeLocation, Element1) <= ClosestPointInTriangle(NewNodeLocation, Element2))
                // {
                //     ClosestElement =Element1;  HaveElement=1;
                //     mMapOfProbNodes[node_index] = 4; 
                //     TRACE("L")
                    
                // }else
                // {
                //     ClosestElement =Element2;  HaveElement=1;
                //     mMapOfProbNodes[node_index] = 4; 
                //     TRACE("K")
                // }

            // unsigned ClosestElement = WhichElement2(P1, P2, NewNodeLocation,  Element1,  Element2);
            // // if (ClosestElement == WinningElement2)
            // // {
            //     mMapOfProbNodes[node_index] = 5; 
            // }

            // unsigned WinningElement = WhichElement(P1, P2, NewNodeLocation, Element1,  Element2);
        
            // PRINT_4_VARIABLES(ClosestElement, WinningElement,Element1, Element2 )
            // assert(WinningElement==WinningElement2);
            // if (( ClosestElement == Element1  || ClosestElement == Element2))
            // {
            //     PRINT_3_VARIABLES(ClosestElement, Element1, Element2)
            // }
            // assert( ClosestElement == Element1  || ClosestElement == Element2);

            // ClosestElement = WinningElement2;
            // assert( ClosestElement == Element1  || ClosestElement == Element2 );



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::WhichElement(c_vector<double, SPACE_DIM> P1, c_vector<double, SPACE_DIM> P2, c_vector<double, SPACE_DIM> NewNodeLocation,unsigned Element1, unsigned Element2)
{
    c_vector<double, 3> Edge = P2 - P1;
    c_vector<double, 3> EdgeMidpoint = (P2 + P1) / 2;

    c_vector<double, SPACE_DIM> Normal1 = GetElementNormal(Element1);
    c_vector<double, SPACE_DIM> Normal2 = GetElementNormal(Element2);

    // Average normal of the two elements
    c_vector<double, SPACE_DIM> Normal = (Normal1 + Normal2) / 2;

    // The normal should be pointing outwards, so the direction of the plane normal will depend on the direction of the edge vector.
    // The edge vector is totally depenedent on the order of the nodes pair, which should be sorted numerically
    c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(Normal, Edge);
    PlaneNormal /= norm_2(PlaneNormal);

    // Need to determine the direction of the normal

    c_vector<double, SPACE_DIM> NormalArm = PlaneNormal + EdgeMidpoint;
    double DistanceToA = norm_2(NormalArm -mCentroidMap[Element1]); // A being the center point of element 1
    double DistanceToB = norm_2(NormalArm -mCentroidMap[Element2]); // A being the center point of element 2
    assert(DistanceToB != DistanceToA);

    if (DistanceToA < DistanceToB)
    {

        // TRACE("The normal is pointing towards element 1 ")
        // I now have a plane, Dependeing on what side of the plane the point is, will give me the closest element, but I dont know which way the plane is pointing
        double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
        if (side>0)
        {
            return Element1;
            // TRACE("The new node is on the element ONE side of the midline plane ")
        }
        else
        {
                // TRACE("The new node is on the element TWO side of the midline plane")
            return Element2;
        }

    }else if  (DistanceToB < DistanceToA)
    {
        // TRACE("The normal is pointing towards element 2 ")
        double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
        if (side>0)
        {
            // TRACE("The new node is on the element TWO side of the midline plane")
            return Element2;
        }
        else
        {
            return Element1;
            //   TRACE("The new node is on the element ONE side of the midline plane")
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::WhichElement2(c_vector<double, SPACE_DIM> P1, c_vector<double, SPACE_DIM> P2, c_vector<double, SPACE_DIM> NewNodeLocation, unsigned Element1, unsigned Element2)
{

    //    FIgure out where it is the same way I used as above in method 2 -- the current way I have doesnt work!!!!!!
            c_vector<double, 3> Edge = P2 - P1;
            c_vector<double, 3> EdgeMidpoint = (P2 + P1) / 2;

            c_vector<double, SPACE_DIM> Normal1 = GetElementNormal(Element1);
            c_vector<double, SPACE_DIM> Normal2 = GetElementNormal(Element2);

            // Average normal of the two elements
            c_vector<double, SPACE_DIM> Normal = (Normal1 + Normal2) / 2;

            // The normal should be pointing outwards, so the direction of the plane normal will depend on the direction of the edge vector.
            // The edge vector is totally depenedent on the order of the nodes pair, which should be sorted numerically
            c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(Normal, Edge);
            PlaneNormal /= norm_2(PlaneNormal);

            // Need to determine the direction of the normal

            c_vector<double, SPACE_DIM> NormalArm = PlaneNormal + EdgeMidpoint;
            double DistanceToA = abs(norm_2(NormalArm - mCentroidMap[Element1])); // A being the center point of element 1
            double DistanceToB = abs(norm_2(NormalArm - mCentroidMap[Element2])); // A being the center point of element 2
            assert(DistanceToB != DistanceToA);
            if (DistanceToA < DistanceToB)
            {
                // TRACE("The normal is pointing towards element 1 ")
                // I now have a plane, Dependeing on what side of the plane the point is, will give me the closest element, but I dont know which way the plane is pointing
                double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
                if (side > 0)
                {
                    // TRACE("The new node is on the element ONE side of the midline plane ")
                    return Element1;
                }
                else
                {
                    // TRACE("The new node is on the element TWO side of the midline plane")
                    return Element2;
                }

            }  else if (DistanceToB > DistanceToA)
            {
            // TRACE("The normal is pointing towards element 2 ")
            double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
            if (side > 0)
            {
                // TRACE("The new node is on the element TWO side of the midline plane")
                return Element2;
            }
            else
            {
                //   TRACE("The new node is on the element ONE side of the midline plane")
                return Element1;
            }
    }
}



// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetClosestElementInOldMeshMethod4(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation)
// {

//     unsigned ClosestElement = -10;
//     double distance =  std::max( (mMaxX - mMinX), (mMaxY - mMinY))/(2*mNx);

//     // Determine which bin I am in, then look in that bin for the element containing this node
//     c_vector<double,3> Bin = DetermineBin(NewNodeLocation);
    
//      for (std::vector<double>::iterator ele_iter = mBinMap[Bin].begin(); ele_iter != mBinMap[Bin].end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
//     {
//         // Now need to find the cloeset element in this bin
//         c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
//         if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
//         {
//             distance = abs(norm_2(NewNodeLocation - Centroid));
//             ClosestElement = *ele_iter;
//             if ( PointInTriangle3D(NewNodeLocation, ClosestElement))
//             {

//                 mMapOfProbNodes[node_index] = 1; 
//                 break;
//             }
//         }
//     }
   


    
//     assert(ClosestElement !=-10);
//     mMapOfProbNodes[node_index] = 3;
//     return ClosestElement;mMinX
// }


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetClosestElementInOldMeshMethod6(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation)
{
    

    unsigned num_elements = this->rGetMesh().GetNumElements();
    mMapOfProbNodes[node_index] = 5; // -5 meaning I cant find the needed element
    double ClosestElement = -10;
    double distance =  1000*std::max( (mMaxX - mMinX), (mMaxY - mMinY))/(2*mNx);

    // Determine which bin I am in, then look in that bin for the element containing this node
    std::vector<int> Bin2 = DetermineBin(NewNodeLocation);
    PRINT_VECTOR(Bin2)
    std::vector<int> Bin = {1,1,1};

    // 

    bool HaveElement = 0;
     for (std::vector<unsigned>::iterator ele_iter = mBinMap[Bin].begin(); ele_iter != mBinMap[Bin].end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
    {
        // Now need to find the cloeset element in this bin
        c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
        if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
        {
            distance = abs(norm_2(NewNodeLocation - Centroid));
            ClosestElement = *ele_iter;
           
        }
    }
    // if ( PointInTriangle3D(NewNodeLocation, ClosestElement))
    // {
    //     HaveElement =1;
    //     mMapOfProbNodes[node_index] = 1; 
    // }
    if (HaveElement ==0)
    {
       typename std::vector<int>::iterator iter = Bin.begin();
        int XBin = *iter;
        std::advance(iter, 1);
        int YBin = *iter;
        std::advance(iter, 1);
        int ZBin = *iter;
        

        std::vector<int> NeighbouringBins; 
        for (int i =-1; i<=1; ++i)
        { 
            if (XBin+i > mNx || XBin+i<1)
            { continue; }
            if (HaveElement ==1 )
            { break;}
            for (int j =-1; j<=1; ++j)
            {
                if (YBin+j > mNy || YBin+j<1|| HaveElement ==1 )
                {  continue;  }
                for (int k =-1; k<=1; ++k)
                {
                    if (ZBin+k > mNz || ZBin+k<1 || (k  ==0 && j ==0 && i==0) || HaveElement ==1)
                    { continue; }
                    // 
                        std::vector<int> NeighBinAddress = {XBin+i, YBin+j, ZBin+k };
                        std::vector<unsigned> NeighBin = mBinMap[NeighBinAddress];
                        // NeighbouringBins.insert(NeighbouringBins.end(), NeighBin.begin(), NeighBin.end());
                        // PRINT_3_VARIABLES(ZBin+k , YBin+j, XBin+i)

                        if (NeighBin.size() <1)
                        {
                            continue;
                        }
                        for (std::vector<unsigned>::iterator ele_iter = NeighBin.begin(); ele_iter != NeighBin.end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
                        {
                            // Now need to find the cloeset element in this bin
                            c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
                            if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
                            {
                                distance = abs(norm_2(NewNodeLocation - Centroid));
                                ClosestElement = *ele_iter;
                                if ( PointInTriangle3D(NewNodeLocation, ClosestElement))
                                {
                                    HaveElement =1;
                                    mMapOfProbNodes[node_index] = 1; 
                                    break;

                                }
                            }
                        }
                }    
            }
        }
    }    
    assert(ClosestElement !=-10 && ClosestElement <= num_elements);
    return ClosestElement;
// 

}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::NewNodeInInitalConfigurationFromChangeOfBasis(unsigned ClosestElement_OldMeshIndex, c_vector<double, SPACE_DIM> NewNode)
{
    Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->rGetMesh().GetElement(ClosestElement_OldMeshIndex);
    // Collect the intial configuration of the nodes, and the deformed config, send them into the mapping function
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
    c_vector<double, SPACE_DIM> z_basis_0 = VectorProduct(vector_12_0, vector_13_0);

    // Nodes at time t

    c_vector<double, SPACE_DIM> vector_12 = DeformedElement[1] - DeformedElement[0];
    c_vector<double, SPACE_DIM> vector_13 = DeformedElement[2] - DeformedElement[0];

    double a = norm_2(vector_12);
    double b = norm_2(vector_13);

    double alpha = acos(inner_prod(vector_12, vector_13) / (a * b));

    // Need to put P into a corrdinate system with the basis vectors
    // vector_12, and vector_13. First need to translate by -Nodes(1,:)
    // beacuse all the points in this triangle have ( including the
    // basis vectors)

    c_vector<double, SPACE_DIM> z_basis = VectorProduct(vector_12, vector_13);
    c_vector<double, SPACE_DIM> P_translated = NewNode - DeformedElement[0];

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

    double C1 = -((p1 * y2 * z3 - p1 * y3 * z2 - p2 * x2 * z3 + p2 * x3 * z2 + p3 * x2 * y3 - p3 * x3 * y2) / (-x1 * y2 * z3 + x1 * y3 * z2 + x2 * y1 * z3 - x2 * y3 * z1 - x3 * y1 * z2 + x3 * y2 * z1));
    double C2 = -((-p1 * y1 * z3 + p1 * y3 * z1 + p2 * x1 * z3 - p2 * x3 * z1 - p3 * x1 * y3 + p3 * x3 * y1) / (-x1 * y2 * z3 + x1 * y3 * z2 + x2 * y1 * z3 - x2 * y3 * z1 - x3 * y1 * z2 + x3 * y2 * z1));
    double C3 = -((-p1 * y1 * z2 + p1 * y2 * z1 + p2 * x1 * z2 - p2 * x2 * z1 - p3 * x1 * y2 + p3 * x2 * y1) / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1));

    // so the point in the coordinate system described by the local element is
    // P_translated  == C1*vector_12 + C2*vector_13 + C3*z_basis
    // Now take these basis vectors and replace them with the vector's original configuration

    c_vector<double, SPACE_DIM> PointInNewRef = C1 * vector_12 + C2 * vector_13 + C3 * z_basis;
    c_vector<double, SPACE_DIM> Difference = PointInNewRef - P_translated;

    assert(norm_2(Difference) < 0.00001);

    c_vector<double, SPACE_DIM> InitalPoint_Translated = C1 * vector_12_0 + C2 * vector_13_0 + C3 * z_basis_0;
    c_vector<double, SPACE_DIM> P_0 = InitalPoint_Translated + Element_0[0];

    return P_0;
}


// For the shear and area forces
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetupMembraneConfiguration()
{
    MathsFunctions M;

    double Kalpha = pow(10, -6.8124);
    double KA = 0;
    double KS = pow(10, -7);

    std::map<unsigned, c_vector<double, 3> > MembraneForceMap;
    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin(); cell_iter != this->mCells.end(); ++cell_iter)
    {
        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);
        MembraneForceMap[node_index] = Create_c_vector(0, 0, 0);
    }

    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
         elem_iter != this->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        Node<SPACE_DIM>* pNode0 = this->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0));
        Node<SPACE_DIM>* pNode1 = this->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(1));
        Node<SPACE_DIM>* pNode2 = this->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));

        c_vector<double, SPACE_DIM> Vec0 = mInitalPositionOfRemeshedNodes[elem_iter->GetNodeGlobalIndex(0)];
        c_vector<double, SPACE_DIM> Vec1 = mInitalPositionOfRemeshedNodes[elem_iter->GetNodeGlobalIndex(1)];
        c_vector<double, SPACE_DIM> Vec2 = mInitalPositionOfRemeshedNodes[elem_iter->GetNodeGlobalIndex(2)];

        // Vectors connecting the nodes
        c_vector<double, SPACE_DIM> vector_12 = Vec1 - Vec0; // Vector 1 to 2
        c_vector<double, SPACE_DIM> vector_13 = Vec2 - Vec0; // Vector 1 to 3

        // define the necessary objects to be used in the loop
        c_vector<double, SPACE_DIM> aVector;
        c_vector<double, SPACE_DIM> bVector;

        // Find the intial area, A0 = 0.5*norm(normal)
        c_vector<double, SPACE_DIM> normalVector = VectorProduct(vector_12, vector_13);
        // PRINT_VECTOR(normalVector);
        double Area = 0.5 * norm_2(normalVector);

        // Save intial area in a map with the element as the key
        mArea0[elem_index] = Area;
        double a = norm_2(vector_12); // Lenght a -> edge connecting P1 and P2
        double b = norm_2(vector_13); // Lenght b -> edge connecting P1 and P3

        double alpha = acos(inner_prod(vector_12, vector_13) / (a * b));

        c_vector<double, 2> x1 = Create_c_vector(0, 0);
        c_vector<double, 2> x2 = Create_c_vector(a, 0);
        c_vector<double, 2> x3 = Create_c_vector(b * cos(alpha), b * sin(alpha));

        //Save the intial position vectors
        mInitalVectors[elem_index][0] = x1;
        mInitalVectors[elem_index][1] = x2;
        mInitalVectors[elem_index][2] = x3;

        // Find the 6 shape function constants
        aVector[0] = (x2[1] - x3[1]) / (2 * Area);
        aVector[1] = (x3[1] - x1[1]) / (2 * Area);
        aVector[2] = (x1[1] - x2[1]) / (2 * Area);

        bVector[0] = (x3[0] - x2[0]) / (2 * Area);
        bVector[1] = (x1[0] - x3[0]) / (2 * Area);
        bVector[2] = (x2[0] - x1[0]) / (2 * Area);

        mACoefficients[elem_index] = aVector;
        mBCoefficients[elem_index] = bVector;

        // Now go through and set the inital strain over the mesh
        vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
        vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 33);

        unsigned node_index;
        // CellPtr p_cell;

        // Get side lengths and angle to create an equal triangle at the origin
        a = norm_2(vector_12); // Lenght a -> edge connecting P0 and P1
        b = norm_2(vector_13); // Lenght b -> edge connecting P0 and P2
        alpha = acos(inner_prod(vector_12, vector_13) / (a * b));

        // This will create an equal triangle at the origin
        x1 = Create_c_vector(0, 0);
        x2 = Create_c_vector(a, 0);
        x3 = Create_c_vector(b * cos(alpha), b * sin(alpha));

        // Get the original triangle
        //Displacement vectors.
        c_vector<double, 2> V1 = x1 - mInitalVectors[elem_index][0];
        c_vector<double, 2> V2 = x2 - mInitalVectors[elem_index][1];
        c_vector<double, 2> V3 = x3 - mInitalVectors[elem_index][2];

        // Get the shape function coefficents for this element
        double Area0 = mArea0[elem_index];

        // Deformation tensor
        double Dxx = 1 + aVector[0] * V1[0] + aVector[1] * V2[0] + aVector[2] * V3[0];
        double Dxy = bVector[0] * V1[0] + bVector[1] * V2[0] + bVector[2] * V3[0];
        double Dyx = aVector[0] * V1[1] + aVector[1] * V2[1] + aVector[2] * V3[1];
        double Dyy = 1 + bVector[0] * V1[1] + bVector[1] * V2[1] + bVector[2] * V3[1];

        c_vector<c_vector<double, 2>, 2> G;

        // G =DTD  -- Caughy green
        G[0][0] = Dxx * Dxx + Dyx * Dyx;
        G[0][1] = Dxx * Dxy + Dyx * Dyy;
        G[1][0] = Dxx * Dxy + Dyy * Dyx;
        G[1][1] = Dxy * Dxy + Dyy * Dyy;

        // Strain invarients
        double I1 = M.tr(G) - 2;
        double I2 = M.det(G) - 1;

        c_vector<c_vector<double, 3>, 3> RotatedForceOnNode;
        c_vector<double, 3> RotatedMag;

        double dedvX;
        double dedvY;

        c_vector<c_vector<double, 3>, 3> X;

        X[0] = Create_c_vector(a, 0, 0);
        X[1] = Create_c_vector(b * cos(alpha), b * sin(alpha), 0);
        X[2] = Create_c_vector(0, 0, 1);
        X = M.MatrixTranspose(X);

        c_vector<c_vector<double, 3>, 3> F_rp;
        c_vector<c_vector<double, 3>, 3> ForceOnNode;

        for (int i = 0; i < 3; i++)
        {
            dedvX = KS * (I1 + 1) * (Dxx * aVector[i] + Dxy * bVector[i]) + (-KS + Kalpha * I2) * ((Dxx * Dyy * Dyy - Dxy * Dyx * Dyy) * aVector[i] + (Dxy * Dyx * Dyx - Dxx * Dyx * Dyy) * bVector[i]);
            dedvY = KS * (I1 + 1) * (Dyx * aVector[i] + Dyy * bVector[i]) + (-KS + Kalpha * I2) * ((Dxx * Dxx * Dyy - Dxx * Dxy * Dyx) * bVector[i] + (Dyx * Dxy * Dxy - Dxx * Dxy * Dyy) * aVector[i]);
            RotatedForceOnNode[i] = -Area0 / 3 * Create_c_vector(dedvX, dedvY, 0);
            F_rp[i] = M.MatrixMultiplication(M.Inverse(X), RotatedForceOnNode[i]);

            // Rotate the force to the original configuretion
            ForceOnNode[i] = F_rp[i][0] * vector_12 + F_rp[i][1] * vector_13 + F_rp[i][2] * X[2];
            // This is the force for each node  for the triangle situated at the origin
            // Matrix with each row containing the force on the corresponding element

            node_index = elem_iter->GetNodeGlobalIndex(i);
            double CellArea = this->GetVolumeOfCell(this->GetCellUsingLocationIndex(node_index));
            // ForceOnNode[i] /= CellArea;
            MembraneForceMap[node_index] += ForceOnNode[i] / CellArea;
        }

        // // Rotate the force to the original configuretion
        // ForceOnNode[0] = F_rp[0][0] * vector_12 + F_rp[0][1] * vector_13 + F_rp[0][2] * X[2];
        // ForceOnNode[1] = F_rp[1][0] * vector_12 + F_rp[1][1] * vector_13 + F_rp[1][2] * X[2];
        // ForceOnNode[2] = F_rp[2][0] * vector_12 + F_rp[2][1] * vector_13 + F_rp[2][2] * X[2];
    }

    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin(); cell_iter != this->mCells.end(); ++cell_iter)
    {
        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);
        Node<3>* pNode = this->rGetMesh().GetNode(node_index);
        // (*cell_iter)->GetCellData()->SetItem("AreaDilationModulus" , 1);

        if ((*cell_iter)->GetCellData()->GetItem("Boundary") == 1)
        {
            c_vector<double, 3> AverageForce = Create_c_vector(0, 0, 0);
            c_vector<unsigned, 3> NearestNodes = this->GetNearestInternalNodes(node_index);
            for (int i = 0; i < 3; i++)
            {
                AverageForce += MembraneForceMap[NearestNodes[i]];
            }

            (*cell_iter)->GetCellData()->SetItem("MembraneForce", norm_2(AverageForce / 3));
        }
        else
        {
            (*cell_iter)->GetCellData()->SetItem("MembraneForce", norm_2(MembraneForceMap[node_index]));
        }
    }
}

// For the bending force
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetInitialAnlgesAcrossMembrane()
{

    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin(); cell_iter != this->mCells.end(); ++cell_iter)
    {
        //find the distance to the nearest neighbour
        (*cell_iter)->GetCellData()->SetItem("Curvature", 0);
    }
    double MeanAngle = 0;
    double counter = 0;
    for (typename MutableMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator edge_iterator = this->rGetMesh().EdgesBegin();
         edge_iterator != this->rGetMesh().EdgesEnd();
         ++edge_iterator)
    {
        std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > UnitNormals;
        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> otherNodes;

        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edge = std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>(edge_iterator.GetNodeA(), edge_iterator.GetNodeB());

        bool boundary_edge_found = CalculateElementNormals(edge, UnitNormals, otherNodes);

        if (boundary_edge_found)
        {
            mOriginalAngles[std::pair<unsigned, unsigned>(edge.first->GetIndex(), edge.second->GetIndex())] = DOUBLE_UNSET;
        }
        else
        {
            // I need this edge to be the inital edge
            Node<SPACE_DIM>* pNode1 = edge_iterator.GetNodeA(); // unsigned Node_1_index = pNode1->GetIndex();
            Node<SPACE_DIM>* pNode3 = edge_iterator.GetNodeB(); // unsigned Node_2_index = pNode3->GetIndex();
            Node<SPACE_DIM>* pNode2 = otherNodes.first; // unsigned Shared_Node_1_index = pNode2->GetIndex();
            Node<SPACE_DIM>* pNode4 = otherNodes.second; // unsigned Shared_Node_2_index = pNode4->GetIndex();

            c_vector<c_vector<double, SPACE_DIM>, 4> PositionVectors;
            // Here is where all the position stuff happens and I can jump in
            PositionVectors[0] = mInitalPositionOfRemeshedNodes[pNode1->GetIndex()];
            PositionVectors[1] = mInitalPositionOfRemeshedNodes[pNode2->GetIndex()];
            PositionVectors[2] = mInitalPositionOfRemeshedNodes[pNode3->GetIndex()];
            PositionVectors[3] = mInitalPositionOfRemeshedNodes[pNode4->GetIndex()];

            // Element 1
            c_vector<double, 3> Element_1_vector_12 = PositionVectors[0] - PositionVectors[3]; // Vector 1 to 2
            c_vector<double, 3> Element_1_vector_13 = PositionVectors[2] - PositionVectors[3]; // Vector 1 to 3

            c_vector<double, 3> normal_1 = VectorProduct(Element_1_vector_12, Element_1_vector_13);

            // Element 2
            c_vector<double, 3> Element_2_vector_12 = PositionVectors[1] - PositionVectors[0]; // Vector 1 to 2
            c_vector<double, 3> Element_2_vector_13 = PositionVectors[2] - PositionVectors[0]; // Vector 1 to 3

            c_vector<double, 3> normal_2 = VectorProduct(Element_2_vector_12, Element_2_vector_13);

            c_vector<double, 2> NodeLocation = Create_c_vector(PositionVectors[2][0], PositionVectors[2][2]);
            c_vector<double, 2> Projection = Create_c_vector(normal_1[0], normal_1[1]);

            normal_2 = normal_2 / norm_2(normal_2);
            normal_1 = normal_1 / norm_2(normal_1);
            double NormalsDot = inner_prod(normal_1, normal_2);
            double Angle = acos(NormalsDot);
            if (NormalsDot == 1)
            {
                Angle = 0;
            }
            mOriginalAngles[std::pair<unsigned, unsigned>(edge.first->GetIndex(), edge.second->GetIndex())] = Angle;

            if (Angle > 10)
            {
                CellPtr p_cell1 = this->GetCellUsingLocationIndex(pNode1->GetIndex());
                CellPtr p_cell2 = this->GetCellUsingLocationIndex(pNode2->GetIndex());
                CellPtr p_cell3 = this->GetCellUsingLocationIndex(pNode3->GetIndex());
                CellPtr p_cell4 = this->GetCellUsingLocationIndex(pNode4->GetIndex());

                p_cell1->GetCellData()->SetItem("Curvature", 1);
                p_cell2->GetCellData()->SetItem("Curvature", 1);
                p_cell3->GetCellData()->SetItem("Curvature", 1);
                p_cell4->GetCellData()->SetItem("Curvature", 1);
            }

            MeanAngle += Angle;
            counter += 1;
        }
    }
    MeanAngle /= counter;
    double variance = 0;

    for (typename MutableMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator edge_iterator = this->rGetMesh().EdgesBegin();
         edge_iterator != this->rGetMesh().EdgesEnd();
         ++edge_iterator)
    {
        std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > UnitNormals;
        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> otherNodes;

        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edge = std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>(edge_iterator.GetNodeA(), edge_iterator.GetNodeB());

        bool boundary_edge_found = CalculateElementNormals(edge, UnitNormals, otherNodes);
        if (!boundary_edge_found)
        {
            variance += pow(mOriginalAngles[std::pair<unsigned, unsigned>(edge.first->GetIndex(), edge.second->GetIndex())] - MeanAngle, 2);
        }
    }
    variance /= counter;

    double stdev = sqrt(variance);
    double Threshold1 = MeanAngle - 1.5 * stdev;
    double Threshold2 = MeanAngle + 1.5 * stdev;
    // Now go through and see what is greater than the threshold

    std::vector<Element<ELEMENT_DIM, SPACE_DIM>*> ElementsToRefine;

    for (typename MutableMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator edge_iterator = this->rGetMesh().EdgesBegin();
         edge_iterator != this->rGetMesh().EdgesEnd();
         ++edge_iterator)
    {
        std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > nonUnitNormals;
        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> otherNodes;

        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edge = std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>(edge_iterator.GetNodeA(), edge_iterator.GetNodeB());

        bool boundary_edge_found = CalculateElementNormals(edge, nonUnitNormals, otherNodes);

        if (!boundary_edge_found)
        {
            double Angle = mOriginalAngles[std::pair<unsigned, unsigned>(edge.first->GetIndex(), edge.second->GetIndex())];

            if (Angle < Threshold1 || Angle > Threshold2)
            {
                // TRACE("High curveature")
                // I need this edge to be the inital edge
                Node<SPACE_DIM>* pNode1 = edge_iterator.GetNodeA(); // unsigned Node_1_index = pNode1->GetIndex();
                Node<SPACE_DIM>* pNode3 = edge_iterator.GetNodeB(); // unsigned Node_2_index = pNode3->GetIndex();
                Node<SPACE_DIM>* pNode2 = otherNodes.first; // unsigned Shared_Node_1_index = pNode2->GetIndex();
                Node<SPACE_DIM>* pNode4 = otherNodes.second; // unsigned Shared_Node_2_index = pNode4->GetIndex();

                // get the cell here and mark it

                CellPtr p_cell1 = this->GetCellUsingLocationIndex(pNode1->GetIndex());
                CellPtr p_cell2 = this->GetCellUsingLocationIndex(pNode2->GetIndex());
                CellPtr p_cell3 = this->GetCellUsingLocationIndex(pNode3->GetIndex());
                CellPtr p_cell4 = this->GetCellUsingLocationIndex(pNode4->GetIndex());
                p_cell1->GetCellData()->SetItem("Curvature", 1);
                p_cell2->GetCellData()->SetItem("Curvature", 1);
                p_cell3->GetCellData()->SetItem("Curvature", 1);
                p_cell4->GetCellData()->SetItem("Curvature", 1);

                // // Need to mark the element

                std::set<unsigned> elements_containing_node1 = pNode1->rGetContainingElementIndices();
                std::set<unsigned> elements_containing_node3 = pNode3->rGetContainingElementIndices();

                // Find common elements
                std::set<unsigned> shared_elements;
                std::set_intersection(elements_containing_node1.begin(),
                                      elements_containing_node1.end(),
                                      elements_containing_node3.begin(),
                                      elements_containing_node3.end(),
                                      std::inserter(shared_elements, shared_elements.begin()));

                std::set<unsigned>::iterator set_iter = shared_elements.begin();
                Element<ELEMENT_DIM, SPACE_DIM>* pElement1 = this->rGetMesh().GetElement(*set_iter);
                ++set_iter;
                Element<ELEMENT_DIM, SPACE_DIM>* pElement2 = this->rGetMesh().GetElement(*set_iter);

                ElementsToRefine.push_back(pElement1);
                ElementsToRefine.push_back(pElement2);
            }
        }
    }

    // for (typename std::vector<Element<ELEMENT_DIM, SPACE_DIM>*>::iterator iter =  ElementsToRefine.begin(); iter!= ElementsToRefine.end(); ++iter)
    // {
    //     // PRINT_VARIABLE(iter)
    // //        for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
    // //      elem_iter != this->rGetMesh().GetElementIteratorEnd();
    // //      ++elem_iter)
    // // {
    //     //    unsigned elem_index = iter->sGetIndex();
    //     // unsigned elem_index = this->rGetMesh().GetIndex(iter);
    //     unsigned temp = iter->GetNodeGlobalIndex(0);
    //         // Node<SPACE_DIM>* pNode0 = this->rGetMesh().GetNode(iter->GetNodeGlobalIndex(0));
    //     Node<SPACE_DIM>* pNode1 = this->rGetMesh().GetNode(iter->GetNodeGlobalIndex(1));
    //     Node<SPACE_DIM>* pNode2 = this->rGetMesh().GetNode(iter->GetNodeGlobalIndex(2));
    // c_vector<double, 3> Centroid = (pNode0->rGetLocation() + pNode1->rGetLocation() + pNode2->rGetLocation()) / 3;

    // Add this as a new node
    //     Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->rGetMesh().GetElement(elem_index);

    // }

    // }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetOriginalAngle(std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edge)
{
    return mOriginalAngles.at(std::pair<unsigned, unsigned>(edge.first->GetIndex(), edge.second->GetIndex()));
}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::map<unsigned, c_vector<double, SPACE_DIM> > HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetInitalNodePositions()
{
    return mOriginalNodePositions;
}

///=====
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<c_vector<double, 2>, 3> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetInitalVectors(unsigned elem_index)
{
    return mInitalVectors[elem_index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<c_vector<double, 3>, 2> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetInitalShapeFunction(unsigned elem_index)
{

    c_vector<double, 3> aVector = mACoefficients[elem_index];
    c_vector<double, 3> bVector = mBCoefficients[elem_index];
    c_vector<c_vector<double, 3>, 2> ShapeFunction;
    ShapeFunction[0] = aVector;
    ShapeFunction[1] = bVector;
    return ShapeFunction;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetOriginalArea(unsigned elem_index)
{
    return mArea0[elem_index];
}





template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetBinningRegions2()
{
    SetMaxDomainDimensions();
    SetCentroidMap();
    // Clear the map so I dont have the old elements and the new elements messing the binning up
    mBinMap.clear();
    // Need to iterate over the elements and determine which bin each element centroid is in. I look for the cloesest old centeroid for each new node, so it fits that I will have the centroids sorted in the bins
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
         elem_iter != this->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        std::vector<int> Bin = {1,1,1};
        mBinMap[Bin].push_back(elem_index);

        // c_vector<double, 3> Centroid = mCentroidMap[elem_index];
        // bool FoundXRegion = 0;
        // bool FoundYRegion = 0;
        // bool FoundZRegion = 0;
        // int xStep = 1;
        // int yStep = 1;
        // int zStep = 1;
        // while (FoundXRegion == 0 && xStep <= mNx)
        // {
        //     if (Centroid[0] - mMinX <= (mMaxX - mMinX) * xStep / mNx) // Look for x bins
        //     {
        //         FoundXRegion = 1; // SO we know it should be in this spacial step
        //         while (FoundYRegion == 0 && yStep <= mNy)
        //         {
        //             if (Centroid[1] - mMinY <= (mMaxY - mMinY) * yStep / mNy) // Look for y bins
        //             {
        //                 FoundYRegion = 1; // We now also have the y bin
        //                 while (FoundZRegion == 0 && zStep <= mNz) // Look for z bins
        //                 {
        //                     if (Centroid[2] - mMinZ <= (mMaxZ - mMinZ) * zStep / mNz)
        //                     {
        //                         FoundZRegion = 1; // We now also have the y bin
        //                     }
        //                     else
        //                     {
        //                         zStep += 1;
        //                     }
        //                 }
        //             }
        //             else
        //             {
        //                 yStep += 1;
        //             }
        //         }
        //     }
        //     else
        //     {
        //         xStep += 1;
        //     }
        // }

        // // assert(xStep ==1 && yStep ==1&& zStep==1);
        // std::vector<int> Bin = {xStep, yStep, zStep};
        // mBinMap[Bin].push_back(elem_index);
        // // // assert(elem_index < this->mrMesh->mElements.size());
        // // // I want a visulisation so I can easily check if my binning is looking okay
        // for (int i = 0; i < 3; i++)
        // {
        //     CellPtr p_cell = this->GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(i));
        //     p_cell->GetCellData()->SetItem("BinX", xStep);
        //     p_cell->GetCellData()->SetItem("BinY", yStep);
        //     p_cell->GetCellData()->SetItem("BinZ", zStep);
        // }
    }
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetMaxDomainDimensions()
{
    // std::map < std::pair<int, int>, std::vector<double>> BinMap;

    // I need the maximal dimensions
    c_vector<double, 3> Location_0 = (this->rGetMesh().GetNodeIteratorBegin())->rGetLocation();
    mMaxX = Location_0[0];
    mMinX = Location_0[0];

    mMaxY = Location_0[1];
    mMinY = Location_0[1];

    mMaxZ = Location_0[2];
    mMinZ = Location_0[2];

    // Start by finding the max x,y,z locations
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        c_vector<double, 3> Location = node_iter->rGetLocation();

        if (Location[0] < mMinX)
        {
            mMinX = Location[0];
        }
        else if (Location[0] > mMaxX)
        {
            mMaxX = Location[0];
        }
        if (Location[1] < mMinY)
        {
            mMinY = Location[1];
        }
        else if (Location[1] > mMaxY)
        {
            mMaxY = Location[1];
        }

        if (Location[2] < mMinZ)
        {
            mMinZ = Location[2];
        }
        else if (Location[2] > mMaxZ)
        {
            mMaxZ = Location[2];
        }
    }

    double Xdom = (mMaxX - mMinX)/10;
    double Ydom = (mMaxY - mMinY)/10;
    double Zdom = (mMaxZ - mMinZ)/10;
    mMaxX +=Xdom;
    mMaxY +=Ydom;
    mMaxZ +=Zdom;

    mMinX -=Xdom;
    mMinY -=Ydom;
    mMinZ -=Zdom;

}

// Set boundaries
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::MarkBoundaryNodes()
{

    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin(); cell_iter != this->mCells.end(); ++cell_iter)
    {
        (*cell_iter)->GetCellData()->SetItem("Boundary", 0);
        (*cell_iter)->GetCellData()->SetItem("FixedBoundary", 0);
    }
    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin(); cell_iter != this->mCells.end(); ++cell_iter)
    {
        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);
        Node<SPACE_DIM>* p_node = this->GetNode(node_index);
        // (*cell_iter)->GetCellData()->SetItem("Boundary", 0);
        std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();

        if (containing_elements.size() < 5)
        {
            (*cell_iter)->GetCellData()->SetItem("Boundary", 1);
            std::set<unsigned> neighbouring_node_indices = this->GetNeighbouringNodeIndices(node_index);
            for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
                 iter != neighbouring_node_indices.end();
                 ++iter)
            {
                CellPtr p_cell = this->GetCellUsingLocationIndex(*iter);
                p_cell->GetCellData()->SetItem("Boundary", 1);
            }
        }
    }

    // TRACE("In the set up solve for the boundaries")
    assert(ELEMENT_DIM == 2);
    assert(SPACE_DIM == 3);
    // std::map<unsigned, c_vector<unsigned, 2> > mNearestNodesMap;
    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin(); cell_iter != this->mCells.end(); ++cell_iter)
    {
        //find the distance to the nearest neighbour
        if ((*cell_iter)->GetCellData()->GetItem("Boundary") == 1)
        {

            unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);
            Node<SPACE_DIM>* pNode = this->rGetMesh().GetNode(node_index);
            c_vector<long double, SPACE_DIM> CellLocation = pNode->rGetLocation();

            c_vector<unsigned, 3> NearestNodes;

            double Distance3 = 0.5;
            double Distance2 = 0.5;
            double Distance1 = 0.5;
            // TRACE("A")
            for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter_Search = this->Begin();
                 cell_iter_Search != this->End();
                 ++cell_iter_Search)
            {
                unsigned node_index_Search = this->GetLocationIndexUsingCell(*cell_iter_Search);
                if (cell_iter_Search->GetCellData()->GetItem("Boundary") == 0 && node_index_Search != node_index)
                {
                    Node<SPACE_DIM>* pNode1 = this->rGetMesh().GetNode(node_index_Search);
                    c_vector<long double, SPACE_DIM> LocationOfNode = pNode1->rGetLocation();
                    double Distance = norm_2(LocationOfNode - CellLocation);

                    if (Distance <= Distance1)
                    {
                        NearestNodes[0] = node_index_Search;
                        Distance1 = Distance;
                    }
                    else if (Distance <= Distance2)
                    {
                        NearestNodes[1] = node_index_Search;
                        Distance2 = Distance;
                    }
                    else if (Distance <= Distance3)
                    {
                        NearestNodes[2] = node_index_Search;
                        Distance3 = Distance;
                    }
                    // else if (Distance <= Distance4)
                    // {
                    //     NearestNodes[3] = node_index_Search;
                    //     Distance4 = Distance;
                    //     // TRACE("D4")
                    // }
                    // else if (Distance <= Distance5)
                    // {
                    //     NearestNodes[4] = node_index_Search;
                    //     Distance5 = Distance;
                    // }
                }
            }
            mNearestNodesMap[node_index] = NearestNodes;
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<unsigned, 3> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetNearestInternalNodes(unsigned node_index)
{
    return mNearestNodesMap[node_index];
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
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetChasteOutputDirectory(std::string ChasteOutputDirectory)
{

    // Note std::string is easier to add together, but in the end what I will be putting into system will be a char *
    char* C = getenv("CHASTE_TEST_OUTPUT");
    std::string directory;
    directory += C;

    directory += "/" + ChasteOutputDirectory;
    mChasteOutputDirectory = directory;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetRelativePath(std::string ChasteOutputDirectory, double startime)
{
    // Note std::string is easier to add together, but in the end what I will be putting into system will be a char *

    std::string directory;
    if (startime == (int)startime)
    {
        directory = ChasteOutputDirectory + "/results_from_time_" + std::to_string((int)startime) + "/";
    }
    else
    {
        std::string TimeStamp = std::to_string(startime);
        TimeStamp.erase(TimeStamp.find_last_not_of('0') + 1, std::string::npos);
        std::cout << TimeStamp << std::endl;
        directory = ChasteOutputDirectory + "/results_from_time_" + TimeStamp + "/";
    }
    mRelativePath = directory;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetRelativePath(std::string ChasteOutputDirectory)
{
    mRelativePath = ChasteOutputDirectory;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetPaths(std::string ChasteOutputDirectory)
{
    mRelativePath = ChasteOutputDirectory;

     // Note std::string is easier to add together, but in the end what I will be putting into system will be a char *
    char* C = getenv("CHASTE_TEST_OUTPUT");
    std::string directory;
    directory += C;

    directory += "/" + ChasteOutputDirectory;
    mChasteOutputDirectory = directory;


}





//============




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
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SaveInitalConditionsAtTime0()
{
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        mInitalPositionOfRemeshedNodes[node_index] = node_iter->rGetLocation();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetMaxEdgelength()
{
    mMaxEdgelength = 0;
    for (typename MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SpringIterator spring_iterator = this->SpringsBegin();
         spring_iterator != this->SpringsEnd();
         ++spring_iterator)
    {
        // mInitalPositionOfRemeshedNodes
        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

        c_vector<double, SPACE_DIM> A_location = this->GetNode(nodeA_global_index)->rGetLocation();
        c_vector<double, SPACE_DIM> B_location = this->GetNode(nodeB_global_index)->rGetLocation();

        if (abs(norm_2(A_location - B_location)) > mMaxEdgelength)
        {
            mMaxEdgelength = abs(norm_2(A_location - B_location));
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetTargetRemeshingEdgeLength(double TargetRemeshingEdgeLength)
{
    mTargetRemeshingEdgeLength = TargetRemeshingEdgeLength;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetTargetRemeshingIterations(int Iterations)
{
    mIterations = Iterations;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetPrintRemeshedIC(bool PrintRemeshedIC)
{
    mPrintRemeshedIC = PrintRemeshedIC;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetRemeshingSoftwear(std::string RemeshingSoftwear)
{
    if (RemeshingSoftwear != "CGAL" && RemeshingSoftwear != "VMTK" && RemeshingSoftwear !="PreAllocatedMatlabMesh")
    {
        EXCEPTION("RemeshingSoftwear must be CGAL or VMTK or PreAllocatedMatlabMesh");
    }
    mRemeshingSoftwear = RemeshingSoftwear;
}



// Set the pathway to the mesh I need to read in
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::PreAllocatedRemeshedMesh(std::string RemeshedMesh)
{
    mPreAllocatedRemeshedMesh = RemeshedMesh;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::JiggleNodes()
{

    clock_t t = clock();
    /// Loop over new and old mesh and find the nodes in the new mesh that are the cloest to the old mesh and put these nodes here --- this should be done in the mesh
    // Loop over old mesh first, then new mesh
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        c_vector<double, SPACE_DIM> OldNodeLocation = node_iter->rGetLocation();

        // Now want to find the cloest node in the new mesh

        double MinDistance = 1000000000000;
        unsigned ClosestNewNode;
        c_vector<double, SPACE_DIM> NewNodeLocation;

        // Now want the location of the nodes in the new mesh

        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator NewN_iter = mNew_mesh.GetNodeIteratorBegin();
             NewN_iter != mNew_mesh.GetNodeIteratorEnd();
             ++NewN_iter)
        {
            NewNodeLocation = NewN_iter->rGetLocation();
            if (norm_2(NewNodeLocation - OldNodeLocation) < MinDistance)
            {
                unsigned Newnode_index = NewN_iter->GetIndex();
                MinDistance = norm_2(NewNodeLocation - OldNodeLocation);
                ClosestNewNode = Newnode_index;
            }
        }
        if (MinDistance < 0.005)
        {

            Node<SPACE_DIM>* p_node = mNew_mesh.GetNode(ClosestNewNode);
            p_node->rGetModifiableLocation() = OldNodeLocation;
        }
    }

    VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(mRelativePath, "RemeshedGeometry2", false);
    mesh_writer.WriteFilesUsingMesh(mNew_mesh);

    t = clock() - t;
    std::cout << "simulation time: " << t / CLOCKS_PER_SEC << " seconds" << std::endl;
}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::CalculateElementNormals(std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edge,
                                                                                        std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> >& nonUnitNormals,
                                                                                        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>& otherNodes)
{
    Node<SPACE_DIM>* pNode1 = edge.first;
    Node<SPACE_DIM>* pNode3 = edge.second;

    /*
     *  Find common triangles
     
      */
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_node1 = pNode1->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_node3 = pNode3->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_node1.begin(),
                          elements_containing_node1.end(),
                          elements_containing_node3.begin(),
                          elements_containing_node3.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    switch (shared_elements.size())
    {
        case 1:
            // We found a boundary edge and we finish here
            return true;
        case 2:
            break;
        default:
            NEVER_REACHED;
    }

    std::set<unsigned>::iterator set_iter = shared_elements.begin();
    Element<2, SPACE_DIM>* pElement1 = this->rGetMesh().GetElement(*set_iter);
    ++set_iter;
    Element<2, SPACE_DIM>* pElement2 = this->rGetMesh().GetElement(*set_iter);

    // Find additional nodes
    Node<SPACE_DIM>* pNode2 = NULL;
    Node<SPACE_DIM>* pNode4 = NULL;
    for (unsigned local_index = 0; local_index < 3; ++local_index)
    {
        unsigned index_for_node2 = pElement1->GetNodeGlobalIndex(local_index);
        unsigned index_for_node4 = pElement2->GetNodeGlobalIndex(local_index);

        if ((index_for_node2 != pNode1->GetIndex()) && (index_for_node2 != pNode3->GetIndex()))
        {
            pNode2 = pElement1->GetNode(local_index);
        }

        if ((index_for_node4 != pNode1->GetIndex()) && (index_for_node4 != pNode3->GetIndex()))
        {
            pNode4 = pElement2->GetNode(local_index);
        }
    }
    assert(pNode2 != NULL);
    assert(pNode4 != NULL);

    // Calculate the force acting on each node
    c_vector<double, SPACE_DIM> vector_A = pNode1->rGetLocation() - pNode3->rGetLocation();
    c_vector<double, SPACE_DIM> vector_B = pNode2->rGetLocation() - pNode3->rGetLocation();
    c_vector<double, SPACE_DIM> normal_1 = VectorProduct(vector_A, vector_B);

    vector_A = pNode4->rGetLocation() - pNode3->rGetLocation();
    vector_B = pNode1->rGetLocation() - pNode3->rGetLocation();

    c_vector<double, SPACE_DIM> normal_2 = VectorProduct(vector_A, vector_B);
    double Area1 = 0.5 * norm_2(normal_1);
    // double Area2 = 0.5 * norm_2(normal_2);
    normal_2 /= norm_2(normal_2);
    normal_1 /= norm_2(normal_1);
    nonUnitNormals = std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> >(normal_1, normal_2);
    otherNodes = std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>(pNode2, pNode4);

    return false;
}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::PointInTriangle3D(c_vector<double, SPACE_DIM> Point, unsigned ClosestElement)
{

    Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->rGetMesh().GetElement(ClosestElement);
    // Collect the intial configuration of the nodes, and the deformed config, send them into the mapping function
    c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> Triangle;
    for (int i = 0; i < 3; ++i)
    {
        Triangle[i] = this->rGetMesh().GetNode(p_element->GetNodeGlobalIndex(i))->rGetLocation();
    }

    // Get normal

    c_vector<double, SPACE_DIM> CenterPoint = (Triangle[1] + Triangle[0] + Triangle[2]) / 3;

    c_vector<double, SPACE_DIM> Line1 = Triangle[1] - Triangle[0]; // Vector 1 to 2
    c_vector<double, SPACE_DIM> Line2 = Triangle[2] - Triangle[0]; // Vector 1 to 3
    c_vector<double, SPACE_DIM> Line3 = Triangle[2] - Triangle[1]; // Vector 2 to 3

    c_vector<double, SPACE_DIM> Midpoint1 = (Triangle[1] + Triangle[0]) / 2;
    c_vector<double, SPACE_DIM> Midpoint2 = (Triangle[2] + Triangle[0]) / 2;
    c_vector<double, SPACE_DIM> Midpoint3 = (Triangle[2] + Triangle[1]) / 2;

    c_vector<double, SPACE_DIM> NormalToElement = VectorProduct(Line1, Line2);
    NormalToElement /= norm_2(NormalToElement);

    c_vector<double, SPACE_DIM> Normal1 = VectorProduct(Line1, NormalToElement);
    Normal1 /= norm_2(Normal1);

    c_vector<double, SPACE_DIM> Normal2 = VectorProduct(Line2, NormalToElement);
    Normal2 /= norm_2(Normal2);

    c_vector<double, SPACE_DIM> Normal3 = VectorProduct(Line3, NormalToElement);
    Normal3 /= norm_2(Normal3);

    // Think I will need to do the anle things to check if the normal is pointing towards the center or not.....

    c_vector<double, SPACE_DIM> Center2MidPoint1 = CenterPoint - Midpoint1;
    c_vector<double, SPACE_DIM> Center2MidPoint2 = CenterPoint - Midpoint2;
    c_vector<double, SPACE_DIM> Center2MidPoint3 = CenterPoint - Midpoint3;
    // Needs to be below plane
    double sign1 = inner_prod(Normal1, (CenterPoint - Midpoint1));
    if (sign1 > 0)
    {
        Normal1 *= -1;
    }

    double sign2 = inner_prod(Normal2, (CenterPoint - Midpoint2));
    if (sign2 > 0)
    {
        Normal2 *= -1;
        sign2 = inner_prod(Normal2, (CenterPoint - Midpoint2));
    }
    double sign3 = inner_prod(Normal3, (CenterPoint - Midpoint3));
    if (sign3 > 0)
    {
        Normal3 *= -1;
        sign3 = inner_prod(Normal3, (CenterPoint - Midpoint3));
    }
    sign1 = inner_prod(Normal1, (Point - Midpoint1));
    sign2 = inner_prod(Normal2, (Point - Midpoint2));
    sign3 = inner_prod(Normal3, (Point - Midpoint3));
    bool Answer;
    if (sign1 <= 0 && sign2 <= 0 && sign3 <= 0)
    {
        Answer = 1;
    }
    else
    {
        Answer = 0;
    }

    // if (SameSideOfPlane(Point,Triangle[0], Triangle[1],Triangle[2]) & SameSideOfPlane(Point,Triangle[1], Triangle[0],Triangle[2])& SameSideOfPlane(Point,Triangle[2], Triangle[0],Triangle[1]))
    // {
    //     Answer = 1;
    // }
    // else
    // {
    //     Answer = 0;
    // }
    return Answer;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::PointInTriangle2D(c_vector<double, SPACE_DIM> Point, unsigned ClosestElement)
{
    bool Answer =0;
    Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->rGetMesh().GetElement(ClosestElement);
    // Collect the intial configuration of the nodes, and the deformed config, send them into the mapping function
    c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> Triangle;
    for (int i = 0; i < 3; ++i)
    {
        Triangle[i] = this->rGetMesh().GetNode(p_element->GetNodeGlobalIndex(i))->rGetLocation();
    }

    if (SameSideOfPlane(Point,Triangle[0], Triangle[1],Triangle[2]) & SameSideOfPlane(Point,Triangle[1], Triangle[0],Triangle[2])& SameSideOfPlane(Point,Triangle[2], Triangle[0],Triangle[1]))
    {
        Answer = 1;
    }
    else
    {
        Answer = 0;
    }
    return Answer;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::ClosestPointInTriangle(c_vector<double, SPACE_DIM> Point, unsigned ClosestElement)
{

    Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->rGetMesh().GetElement(ClosestElement);
    // Collect the intial configuration of the nodes, and the deformed config, send them into the mapping function
    c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> Triangle;
    for (int i = 0; i < 3; ++i)
    {
        Triangle[i] = this->rGetMesh().GetNode(p_element->GetNodeGlobalIndex(i))->rGetLocation();
    }

    // Get normal
    // c_vector<double, SPACE_DIM> CenterPoint = (Triangle[1] + Triangle[0] + Triangle[2]) / 3;

    c_vector<double, SPACE_DIM> Line1 = Triangle[1] - Triangle[0]; // Vector 1 to 2
    c_vector<double, SPACE_DIM> Line2 = Triangle[2] - Triangle[0]; // Vector 1 to 3
    c_vector<double, SPACE_DIM> Line3 = Triangle[2] - Triangle[1]; // Vector 2 to 3

    c_vector<double, SPACE_DIM> ElementNormal =  GetElementNormal(ClosestElement);

    // Need to create a plane containing the triangle. 
    // Let  Triangle[0] be the origin
    c_vector<double, SPACE_DIM> v = Point-Triangle[0];
    double dist = inner_prod(ElementNormal, v);
    c_vector<double, SPACE_DIM> projected_point = Point - dist*ElementNormal;
    // Check that my normal was facing the right way - 
    assert( inner_prod(ElementNormal, projected_point-Triangle[0]) < 1e-13 );

    if (PointInTriangle3D(projected_point, ClosestElement) ==1)
    {
        return dist;
    }
    else{
    // Get the closest point on each line --- 
    
        c_vector<double, SPACE_DIM> Normal1 = VectorProduct(Line1, ElementNormal);
        Normal1 /= norm_2(Normal1);

        c_vector<double, SPACE_DIM> Normal2 = VectorProduct(Line2, ElementNormal);
        Normal2 /= norm_2(Normal2);

        c_vector<double, SPACE_DIM> Normal3 = VectorProduct(Line3, ElementNormal);
        Normal3 /= norm_2(Normal3);

        double dist1 = abs(inner_prod(Normal1, (projected_point-Triangle[1])));
        double dist2 = abs(inner_prod(Normal2, (projected_point-Triangle[2])));
        double dist3 = abs(inner_prod(Normal3, (projected_point-Triangle[2])));
        // Find the smallest distance 

        double MinDistane = std::min(dist1, std::min(dist2, dist3));

        double PointDistanceToElement = sqrt(MinDistane*MinDistane +dist*dist );
        return PointDistanceToElement;
    }

 
}





template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::PointInTriangle3D(c_vector<double, SPACE_DIM> Point, c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> Triangle)
{

    // Get normals
    c_vector<double, SPACE_DIM> CenterPoint = (Triangle[1] + Triangle[0] + Triangle[2]) / 3;

    c_vector<double, SPACE_DIM> Line1 = Triangle[1] - Triangle[0]; // Vector 1 to 2
    c_vector<double, SPACE_DIM> Line2 = Triangle[2] - Triangle[0]; // Vector 1 to 3
    c_vector<double, SPACE_DIM> Line3 = Triangle[2] - Triangle[1]; // Vector 2 to 3

    c_vector<double, SPACE_DIM> Midpoint1 = (Triangle[1] + Triangle[0]) / 2;
    c_vector<double, SPACE_DIM> Midpoint2 = (Triangle[2] + Triangle[0]) / 2;
    c_vector<double, SPACE_DIM> Midpoint3 = (Triangle[2] + Triangle[1]) / 2;

    c_vector<double, SPACE_DIM> NormalToElement = VectorProduct(Line1, Line2);
    NormalToElement /= norm_2(NormalToElement);

    c_vector<double, SPACE_DIM> Normal1 = VectorProduct(Line1, NormalToElement);
    Normal1 /= norm_2(Normal1);

    c_vector<double, SPACE_DIM> Normal2 = VectorProduct(Line2, NormalToElement);
    Normal2 /= norm_2(Normal2);

    c_vector<double, SPACE_DIM> Normal3 = VectorProduct(Line3, NormalToElement);
    Normal3 /= norm_2(Normal3);

    // Think I will need to do the anle things to check if the normal is pointing towards the center or not.....

    c_vector<double, SPACE_DIM> Center2MidPoint1 = CenterPoint - Midpoint1;
    c_vector<double, SPACE_DIM> Center2MidPoint2 = CenterPoint - Midpoint2;
    c_vector<double, SPACE_DIM> Center2MidPoint3 = CenterPoint - Midpoint3;
    // Needs to be below plane
    double sign1 = inner_prod(Normal1, (CenterPoint - Midpoint1));
    if (sign1 > 0)
    {
        Normal1 *= -1;
    }

    double sign2 = inner_prod(Normal2, (CenterPoint - Midpoint2));
    if (sign2 > 0)
    {
        Normal2 *= -1;
        sign2 = inner_prod(Normal2, (CenterPoint - Midpoint2));
    }
    double sign3 = inner_prod(Normal3, (CenterPoint - Midpoint3));
    if (sign3 > 0)
    {
        Normal3 *= -1;
        sign3 = inner_prod(Normal3, (CenterPoint - Midpoint3));
    }
    sign1 = inner_prod(Normal1, (Point - Midpoint1));
    sign2 = inner_prod(Normal2, (Point - Midpoint2));
    sign3 = inner_prod(Normal3, (Point - Midpoint3));
    bool Answer;
    if (sign1 <= 0 && sign2 <= 0 && sign3 <= 0)
    {
        Answer = 1;
    }
    else
    {
        Answer = 0;
    }
    bool Answer1;
    if (SameSideOfPlane(Point,Triangle[0], Triangle[1],Triangle[2]) & SameSideOfPlane(Point,Triangle[1], Triangle[0],Triangle[2])& SameSideOfPlane(Point,Triangle[2], Triangle[0],Triangle[1]))
    {
        Answer1 = 1;
    }
    else
    {
        Answer1 = 0;
    }
    if (Answer1 !=1)
    {
        PRINT_2_VARIABLES(Answer, Answer1)
    }
    return Answer;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::PointInTriangle3D2(c_vector<double, SPACE_DIM> Point, unsigned ClosestElement)
{

    Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->rGetMesh().GetElement(ClosestElement);
    // Collect the intial configuration of the nodes, and the deformed config, send them into the mapping function
    c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> Triangle;
    for (int i = 0; i < 3; ++i)
    {
        Triangle[i] = this->rGetMesh().GetNode(p_element->GetNodeGlobalIndex(i))->rGetLocation();
    }

    bool Answer;
    if (SameSide(Point, Triangle[0], Triangle[1], Triangle[2]) & SameSide(Point, Triangle[1], Triangle[0], Triangle[2]) & SameSide(Point, Triangle[2], Triangle[0], Triangle[1]))
    {
        Answer = 1;
    }
    else
    {
        Answer = 0;
    }
    return Answer;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SameSideOfPlane(c_vector<double, SPACE_DIM> P1, c_vector<double, SPACE_DIM> P2, c_vector<double, SPACE_DIM> a, c_vector<double, SPACE_DIM> b)
{
    bool Answer;
    ;
    // % Need to create plane from one of the lines and the opposite point -- here
    // % order of points is very important... dam
    c_vector<double, SPACE_DIM> Line = b - a;
    c_vector<double, SPACE_DIM> Midpoint = (a + b) / 2;

    c_vector<double, SPACE_DIM> vector_12 = a - P2; // Vector 1 to 2
    c_vector<double, SPACE_DIM> vector_13 = b - P2; // Vector 1 to 3

    c_vector<double, SPACE_DIM> NormalToElement = -VectorProduct(vector_12, vector_13);
    NormalToElement /= norm_2(NormalToElement);
    c_vector<double, SPACE_DIM> NormalToPlane = VectorProduct(NormalToElement, Line);
    NormalToPlane /= -norm_2(NormalToPlane);

    // Need to check which way the normal is pointing.... need the normal
    if (norm_2(P2 - (Midpoint + 0.000001 * NormalToPlane)) < norm_2(P2 - Midpoint))
    {
        // The normal is pointing the wrong way, flip it around
        NormalToPlane *= -1;
    }

    // Now determine if this is above or below plane
    double sign = inner_prod(NormalToPlane, (P1 - Midpoint));
    // If sign is neg, point is below plane, i.e on the triangle side of the
    // plane, is sign is pos, the point is away from the triangle

    if (sign <= 0)
    {
        Answer = 1;
    }
    else
    {
        Answer = 0;
    }
    return Answer;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SameSide(c_vector<double, SPACE_DIM> P1, c_vector<double, SPACE_DIM> P2, c_vector<double, SPACE_DIM> a, c_vector<double, SPACE_DIM> b)
{
    bool Answer;
    // % Need to create plane from one of the lines and the opposite point -- here
    c_vector<double, SPACE_DIM> A = b - a;
    c_vector<double, SPACE_DIM> B = P1 - a;
    c_vector<double, SPACE_DIM> C = P2 - a;

    c_vector<double, SPACE_DIM> cp1 = VectorProduct(A, B);
    c_vector<double, SPACE_DIM> cp2 = VectorProduct(A, C);
    if (inner_prod(cp1, cp2) >= 0)
    {
        Answer = 1;
    }
    else
    {
        Answer = 0;
    }
    return Answer;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetBinningIntervals(int nX, int nY, int nZ)
{
    mNx = nX;
    mNy = nY;
    mNz = nZ;
}



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

// VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mChasteOutputDirectory + "RemeshedGeometry.vtu");
// mNew_mesh.ConstructFromMeshReader(mesh_reader);

//  TRACE("History dependent remeshing")

// static_cast<HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>&>((this->mrMesh)).AddANewNodeBehindBoundary();
//     VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer("CheckingremeshsedMesh", "config", false);
// mesh_writer.WriteFilesUsingMesh(static_cast<AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>&>((this->mrMesh)));

//     // Now look at adding a new cell for the new node
//     // Option 1, follow the line of constructors to fine where I need to clear/delete cells and add new cells on

// old, when using vmtk and stls
/*
	 * Here I will remesh the current chaste geometery
     * 1) Take the latest .vtu file convert it into an .stl 
     * 2) Remesh this stl
     * 3) Convert this stl into a .vtu so I can read it into chaste and use it  
	 */

// // Convert the .vtu into an .stl

// std::string stlfile = mChasteOutputDirectory + "CurrentMesh.stl";
// TRACE("Create stl")
// PRINT_VARIABLE(mChasteOutputDirectory)
// std::string vtu2stlCommand = "python projects/VascularRemodelling/apps/ConvertCurrentGeometryToSTL.py -ChasteOutput " + mChasteOutputDirectory + "  -stlOutput " + stlfile;
// std::system(vtu2stlCommand.c_str()); // system only takes char *

// // Remesh the .stl
// std::string Remeshedstl = mChasteOutputDirectory + "RemeshedGeometry.stl";
// TRACE("Remeshing geometry")

// // More info for the remeshing options found at https://sourceforge.net/p/vmtk/mailman/message/29391820/ and  http://www.vmtk.org/vmtkscripts/vmtksurfaceremeshing.html
// std::string RemeshCommand = "vmtksurfaceremeshing -ifile " + stlfile + " -iterations " + std::to_string(mIterations) + " -area " + std::to_string(mTargetRemeshingElementArea) + " -maxarea "+ std::to_string(1.2*mTargetRemeshingElementArea) +" -ofile " + Remeshedstl;
// std::system(RemeshCommand.c_str());

// //Finally conver the Remeshed stl into a vtu -- meshio is your friend
// std::string Remeshedvtu = mChasteOutputDirectory + "RemeshedGeometry.vtu";
// std::string stl2vtuCommand = " meshio-convert " + Remeshedstl + " " + Remeshedvtu;
// std::system(stl2vtuCommand.c_str());

//         // Read in the new mesh
// VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mChasteOutputDirectory + "RemeshedGeometry.vtu");
// mNew_mesh.ConstructFromMeshReader(mesh_reader);
// TRACE("Have the new mesh ;) ");

// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::CheckCurvature()
// {
//     // Loop over all springs and check the angles between each
//     for (typename MeshBasedCellPopulation<2, 3>::SpringIterator spring_iterator = p_static_cast_cell_population->SpringsBegin();
//          spring_iterator != p_static_cast_cell_population->SpringsEnd();
//          ++spring_iterator)
//     {

//         Node<3>* pNode1 = spring_iterator.GetNodeA();
//         Node<3>* pNode3 = spring_iterator.GetNodeB();

//         std::pair<Node<3>*, Node<3>*> edge = std::pair<Node<3>*, Node<3>*>(pNode1, pNode3);
//         double OriginalAngle = p_static_cast_cell_population->GetOriginalAngle(edge);

//         std::pair<c_vector<double, 3>, c_vector<double, 3> > nonUnitNormals;

//         std::pair<Node<3>*, Node<3>*> otherNodes;

//         bool boundary_edge_found = p_static_cast_cell_population->CalculateElementNormals( edge, nonUnitNormals, otherNodes);

// }

// ----------------------------------------------------------------------

/*

        //  c_vector<double,3> Bin(xStep, yStep, zStep); //This is the key in the map, then the value is some vector that this is added to
       

        // // Now I need to get the nearest centroid
        // for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
        //      Centroid_iter != mCentroidMap.end();
        //      ++Centroid_iter)
        // {
        //     // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
        //     // if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance && PointInTriangle3D(NewNodeLocation,  Centroid_iter->first))
        //     // {
        //      if (PointInTriangle3D(NewNodeLocation,  Centroid_iter->first) ==1)
        //     {
        //         // now test if it is in the element %
        //         ClosestElement = Centroid_iter->first;
        //         break;

        //     }
        // }

        //  if (ClosestElement == -10)  //
        // {
        //     TRACE("If this fails the new mesh might not be a good representation of the old mesh, will need to loop over everything again")
        //     TRACE("first try 2D edge checking .....")
        //     for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
        //      Centroid_iter != mCentroidMap.end();
        //      ++Centroid_iter)
        //     {
        //         // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
        //           if (PointInTriangle3D2(NewNodeLocation,  Centroid_iter->first) ==1)
        //             {
        //                 // now test if it is in the element %
        //                 ClosestElement = Centroid_iter->first;
        //                 break;

        //             }
        //     }
        // }

        // if (ClosestElement == -10)  //
        // {
        // TRACE(" 2D edge checking didnt work, now need to try cloeset center point -- there is something wrong")
        for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
             Centroid_iter != mCentroidMap.end();
             ++Centroid_iter)
        {
            // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
            if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance)
            {
                // now test if it is in the element %
                ClosestElement = Centroid_iter->first;
                distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
            }
        }
       

        // PRINT_VARIABLE(InTriangle)
        distance *= 10;
        unsigned WrongClosestElement = ClosestElement;
        if (InTriangle == 0)
        {

            // this is the wrong triangle, try the next one
            for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
                 Centroid_iter != mCentroidMap.end();
                 ++Centroid_iter)
            {
                // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
                if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance && Centroid_iter->first != WrongClosestElement)
                {
                    // now test if it is in the element %
                    ClosestElement = Centroid_iter->first;
                    distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
                }
            }
        }
        InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
        distance *= 10;
        unsigned WrongClosestElement2 = ClosestElement;
        if (InTriangle == 0)
        {
            TRACE("This is the third try")
            // this is the wrong triangle, try the next one
            for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
                 Centroid_iter != mCentroidMap.end();
                 ++Centroid_iter)
            {
                // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
                if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance && Centroid_iter->first != WrongClosestElement && Centroid_iter->first != WrongClosestElement2)
                {
                    // now test if it is in the element %
                    ClosestElement = Centroid_iter->first;
                    distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
                }
            }
        }
        InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
        distance *= 10;
        unsigned WrongClosestElement3 = ClosestElement;

        if (InTriangle == 0)
        {
            TRACE("This is the fourth try")
            // this is the wrong triangle, try the next one
            for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
                 Centroid_iter != mCentroidMap.end();
                 ++Centroid_iter)
            {
                // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
                if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance && Centroid_iter->first != WrongClosestElement && Centroid_iter->first != WrongClosestElement2 && Centroid_iter->first != WrongClosestElement3)
                {
                    // now test if it is in the element %
                    ClosestElement = Centroid_iter->first;
                    distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
                }
            }
        }

        InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
        distance *= 10;
        unsigned WrongClosestElement4 = ClosestElement;

        if (InTriangle == 0)
        {
            TRACE("This is the fith try")
            // this is the wrong triangle, try the next one
            for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
                 Centroid_iter != mCentroidMap.end();
                 ++Centroid_iter)
            {
                // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
                if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance && Centroid_iter->first != WrongClosestElement4 && Centroid_iter->first != WrongClosestElement && Centroid_iter->first != WrongClosestElement2 && Centroid_iter->first != WrongClosestElement3)
                {
                    // now test if it is in the element %
                    ClosestElement = Centroid_iter->first;
                    distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
                }
            }
        }
        InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
        distance *= 10;
        unsigned WrongClosestElement5 = ClosestElement;

        if (InTriangle == 0)
        {
            TRACE("This is the sixth try")
            // this is the wrong triangle, try the next one
            for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
                 Centroid_iter != mCentroidMap.end();
                 ++Centroid_iter)
            {
                // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
                if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance && Centroid_iter->first != WrongClosestElement5 && Centroid_iter->first != WrongClosestElement4 && Centroid_iter->first != WrongClosestElement && Centroid_iter->first != WrongClosestElement2 && Centroid_iter->first != WrongClosestElement3)
                {
                    // now test if it is in the element %
                    ClosestElement = Centroid_iter->first;
                    distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
                }
            }
        }
        InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);

        if (InTriangle == 0)
        {
            ClosestElement = WrongClosestElement;
        }

        // PRINT_VARIABLE(InTriangle)
        // I want to know what this says about being in the same ele
        // }
        */

// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::UpdateCellData()
// {
//      // When outputting any CellData, we assume that the first cell is representative of all cells
//     unsigned num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
//     std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();

//     std::vector<std::vector<double> > cell_data;
//     for (unsigned var=0; var<num_cell_data_items; var++)
//     {
//         std::vector<double> cell_data_var(num_cells_from_mesh);
//         cell_data.push_back(cell_data_var);
//     }

//     if (mOutputMeshInVtk)
//     {
//         // Create mesh writer for VTK output
//         VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(rDirectory, "mesh_"+time.str(), false);
//         mesh_writer.WriteFilesUsingMesh(rGetMesh());
//     }

// }




// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetBinningRegions()
// {
//     SetMaxDomainDimensions();

//     SetCentroidMap();
//     // Clear the map so I dont have the old elements and the new elements messing the binning up

//     mBinMap.clear();
//     // Need to iterate over the elements and determine which bin each element centroid is in. I look for the cloesest old centeroid for each new node, so it fits that I will have the centroids sorted in the bins
//     for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
//          elem_iter != this->rGetMesh().GetElementIteratorEnd();
//          ++elem_iter)
//     {
//         unsigned elem_index = elem_iter->GetIndex();
//         c_vector<double, 3> Centroid = mCentroidMap[elem_index];
//         bool FoundXRegion = 0;
//         bool FoundYRegion = 0;
//         bool FoundZRegion = 0;
//         int xStep = 1;
//         int yStep = 1;
//         int zStep = 1;
//         double FudgeX = mFudgeScalling * (mMaxX - mMinX);
//         double FudgeY = mFudgeScalling * (mMaxY - mMinY);
//         double FudgeZ = mFudgeScalling * (mMaxZ - mMinZ);

//         bool FudgeInX = 0;
//         bool FudgeInY = 0;
//         bool FudgeInZ = 0;
//         bool lowerFudgeX = 0;
//         bool lowerFudgeY = 0;
//         bool lowerFudgeZ = 0;
//         while (FoundXRegion == 0 && xStep <= mNx)
//         {
//             if (Centroid[0] - mMinX <= (mMaxX - mMinX) * xStep / mNx) // Look for x bins
//             {
//                 FoundXRegion = 1; // SO we know it should be in this spacial step
//                 while (FoundYRegion == 0 && yStep <= mNy)
//                 {
//                     if (Centroid[1] - mMinY <= (mMaxY - mMinY) * yStep / mNy) // Look for y bins
//                     {
//                         FoundYRegion = 1; // We now also have the y bin
//                         while (FoundZRegion == 0 && zStep <= mNz) // Look for z bins
//                         {
//                             if (Centroid[2] - mMinZ <= (mMaxZ - mMinZ) * zStep / mNz)
//                             {
//                                 FoundZRegion = 1; // We now also have the y bin

//                                 if (Centroid[2] - mMinZ >= (mMaxZ - mMinZ) * zStep / mNz - FudgeZ && zStep + 1 <= mNz) // Look for x bins
//                                 {
//                                     FudgeInZ = 1;
//                                     // SHould also be in the next X bin across
//                                     c_vector<double,3> Bin= Create_c_vector(xStep, yStep, zStep + 1); //This is the key in the map, then the value is some vector that this is added to
//                                     mBinMap[Bin].push_back(elem_index);
//                                 }
//                                 else if (Centroid[2] - mMinZ <= (mMaxZ - mMinZ) * (zStep - 1) / mNz + FudgeZ && zStep > 1) // // What about the region behind?????
//                                 {
//                                     lowerFudgeZ = 1;
//                                     c_vector<double,3> Bin= Create_c_vector(xStep, yStep, zStep - 1); // SHould also be in the next X bin across
//                                     mBinMap[Bin].push_back(elem_index);
//                                 }
//                             }
//                             else
//                             {
//                                 zStep += 1;
//                             }
//                         }

//                         if (Centroid[1] - mMinY >= (mMaxY - mMinY) * yStep / mNy - FudgeY && yStep + 1 <= mNy) // Look for x bins
//                         {
//                             FudgeInY = 1;
//                             // SHould also be in the next X bin across
//                             c_vector<double,3> Bin= Create_c_vector(xStep, yStep + 1, zStep); //This is the key in the map, then the value is some vector that this is added to
//                             mBinMap[Bin].push_back(elem_index);

//                             if (FudgeInZ)
//                             {
//                                 c_vector<double,3> Bin= Create_c_vector(xStep, yStep + 1, zStep + 1); //This is the key in the map, then the value is some vector that this is added to
//                                 mBinMap[Bin].push_back(elem_index);
//                             }
//                         }
//                         else if (Centroid[1] - mMinY <= (mMaxY - mMinY) * (yStep - 1) / mNy + FudgeY && yStep > 1) // // What about the region behind?????
//                         {
//                             lowerFudgeY = 1;
//                             c_vector<double,3> Bin= Create_c_vector(xStep, yStep - 1, zStep); // SHould also be in the next X bin across
//                             mBinMap[Bin].push_back(elem_index);

//                             if (lowerFudgeZ)
//                             {
//                                 c_vector<double,3> Bin= Create_c_vector(xStep, yStep - 1, zStep - 1); //This is the key in the map, then the value is some vector that this is added to
//                                 mBinMap[Bin].push_back(elem_index);
//                             }
//                         }
//                     }
//                     else
//                     {
//                         yStep += 1;
//                     }
//                 }
//                 // Is this element also in the fudge region of the next bin across?
//                 if (Centroid[0] - mMinX >= (mMaxX - mMinX) * xStep / mNx - FudgeX && xStep + 1 <= mNx) // Look for x bins
//                 {
//                     FudgeInX = 1;
//                     c_vector<double,3> Bin= Create_c_vector(xStep + 1, yStep, zStep); // SHould also be in the next X bin across
//                     mBinMap[Bin].push_back(elem_index);
//                     // Also want to account for the fudge in the other directions :)
//                     if (FudgeInY)
//                     {
//                         c_vector<double,3> Bin= Create_c_vector(xStep + 1, yStep + 1, zStep); //This is the key in the map, then the value is some vector that this is added to
//                         mBinMap[Bin].push_back(elem_index);
//                         if (FudgeInZ)
//                         {
//                             c_vector<double,3> Bin= Create_c_vector(xStep + 1, yStep + 1, zStep + 1); //This is the key in the map, then the value is some vector that this is added to
//                             mBinMap[Bin].push_back(elem_index);
//                         }
//                     }
//                     if (FudgeInZ)
//                     {
//                         c_vector<double,3> Bin= Create_c_vector(xStep + 1, yStep, zStep + 1); //This is the key in the map, then the value is some vector that this is added to
//                         mBinMap[Bin].push_back(elem_index);
//                     }
//                 }
//                 else if (Centroid[0] - mMinX <= (mMaxX - mMinX) * (xStep - 1) / mNx + FudgeX && xStep > 1) // What about the region behind?????
//                 {
//                     lowerFudgeX = 1;
//                     c_vector<double,3> Bin= Create_c_vector(xStep - 1, yStep, zStep); // SHould also be in the next X bin across
//                     mBinMap[Bin].push_back(elem_index);

//                     // Also want to account for the fudge in the other directions :)
//                     if (lowerFudgeY)
//                     {
//                         c_vector<double,3> Bin= Create_c_vector(xStep - 1, yStep - 1, zStep); //This is the key in the map, then the value is some vector that this is added to
//                         mBinMap[Bin].push_back(elem_index);
//                         if (lowerFudgeZ)
//                         {
//                             c_vector<double,3> Bin= Create_c_vector(xStep - 1, yStep - 1, zStep - 1); //This is the key in the map, then the value is some vector that this is added to
//                             mBinMap[Bin].push_back(elem_index);
//                         }
//                     }
//                     if (lowerFudgeZ)
//                     {
//                         c_vector<double,3> Bin =Create_c_vector(xStep - 1, yStep, zStep - 1); //This is the key in the map, then the value is some vector that this is added to
//                         mBinMap[Bin].push_back(elem_index);
//                     }
//                 }
//             }
//             else
//             {
//                 xStep += 1;
//             }
//         }

//         c_vector<double,3> Bin= Create_c_vector(xStep, yStep, zStep); //This is the key in the map, then the value is some vector that this is added to

//         // I think the fugde region needs to be implemented else where -- there are still some problems with my nearest element search too --- this could be causing me pain
//         mBinMap[Bin].push_back(elem_index);
//         // TRACE(" In bin mapping")
//         // PRINT_VARIABLE(elem_index)
//         // assert(elem_index < this->mrMesh->mElements.size());

//         // I want a visulisation so I can easily check if my binning is looking okay
//         for (int i = 0; i < 3; i++)
//         {
//             CellPtr p_cell = this->GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(i));
//             p_cell->GetCellData()->SetItem("BinX", xStep);
//             p_cell->GetCellData()->SetItem("BinY", yStep);
//             p_cell->GetCellData()->SetItem("BinZ", zStep);

//             if (lowerFudgeZ == 1 || lowerFudgeX == 1 || lowerFudgeY == 1)
//             {
//                 p_cell->GetCellData()->SetItem("FudgeL", 1);
//             }
//             else
//             {
//                 p_cell->GetCellData()->SetItem("FudgeL", 0);
//             }
//             if (FudgeInZ == 1 || FudgeInX == 1 || FudgeInY == 1)
//             {
//                 p_cell->GetCellData()->SetItem("Fudge", 1);
//             }
//             else
//             {
//                 p_cell->GetCellData()->SetItem("Fudge", 0);
//             }
//         }
//     }
// }





// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetClosestElementInOldMeshMethod3(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation)
// {

//     mMapOfProbNodes[node_index] = 5; // -5 meaning I cant find the needed element
//     unsigned ClosestElement = -10;
//     double distance =  std::max( (mMaxX - mMinX), (mMaxY - mMinY))/(2*mNx);

//     // Determine which bin I am in, then look in that bin for the element containing this node
//     c_vector<double,3> Bin = DetermineBin(NewNodeLocation);
// std::vector<double> CurrentBin = mBinMap[Bin];
//     bool HaveElement = 0;
//      for (std::vector<double>::iterator ele_iter = CurrentBin.begin(); ele_iter != CurrentBin.end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
//     {
//         // Now need to find the cloeset element in this bin
//         c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
//         if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
//         {
//             distance = abs(norm_2(NewNodeLocation - Centroid));
//             ClosestElement = *ele_iter;
//             if ( PointInTriangle3D(NewNodeLocation, ClosestElement))
//             {
//                 HaveElement =1;
//                 mMapOfProbNodes[node_index] = 1; 
//                 break;
//             }
//         }
//     }
   
//     // unsigned num_elements = this->rGetMesh().GetNumElements();
    

//     // if (HaveElement == 0)
//     // {
//     //     distance = mMaxEdgelength;
//     //     std::pair<unsigned,unsigned> ClosestEdge;
//     //     c_vector<double, SPACE_DIM> ClosestEdgeNorm ;
//     //     c_vector<double, SPACE_DIM> ClosestMidPoint;
//     //      for (typename MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SpringIterator spring_iterator = this->SpringsBegin();
//     //      spring_iterator != this->SpringsEnd();
//     //      ++spring_iterator)
//     //      {
//     //         // mInitalPositionOfRemeshedNodes
//     //         unsigned NodeIndexA = spring_iterator.GetNodeA()->GetIndex();
//     //         unsigned NodeIndexB = spring_iterator.GetNodeB()->GetIndex();

//     //         c_vector<double, SPACE_DIM> P1 = this->GetNode(NodeIndexA)->rGetLocation();
//     //         c_vector<double, SPACE_DIM> P2 = this->GetNode(NodeIndexB)->rGetLocation();
//     //         c_vector<double, SPACE_DIM> EdgeVector = P2- P1;
//     //         c_vector<double, SPACE_DIM> MidPoint = (P2+ P1)/2;
//     //         EdgeVector/=norm_2(EdgeVector);
//     //         c_vector<double, SPACE_DIM> V = NewNodeLocation - P1;
 
//     //         // Normal to the plane containing the edge and the new node
//     //         c_vector<double, SPACE_DIM> N = VectorProduct(EdgeVector, V)/norm_2(V);
//     //         N/=norm_2(N);
//     //         // Now can find the normal to the edge that is contaned in the plane defined by the three points
            
//     //         c_vector<double, SPACE_DIM> normal = VectorProduct(EdgeVector, N);
//     //         normal/=norm_2(normal);
//     //         if (norm_2(MidPoint-normal -NewNodeLocation) < norm_2(MidPoint+normal -NewNodeLocation) )
//     //         {
//     //             normal*=-1;
//     //         }
//     //         // Need to check the direction of the normal :) 

//     //         // normal to edge will be normal to the edge and the element normal -- se get the cross product of
//     //         double theta = acos(inner_prod(normal, V/norm_2(V)) );
//     //         double d =  norm_2(V) *cos(theta);
//     //         if (d< distance )
//     //         {
//     //             distance =d;
//     //             // need to think about what to do here, 
//     //             ClosestEdge = std::pair<unsigned,unsigned>(NodeIndexA, NodeIndexB);
//     //             ClosestEdgeNorm = normal;
//     //             ClosestMidPoint = MidPoint;
//     //         }
            
//     //     }

//     //     c_vector<double, SPACE_DIM> P1 = this->GetNode(ClosestEdge.first)->rGetLocation();
//     //     c_vector<double, SPACE_DIM> P2 = this->GetNode(ClosestEdge.second)->rGetLocation();

//     //     //Check that this is close enough // how ?? 
//     //     if (distance < norm_2(P1 - NewNodeLocation )  || distance < norm_2(P2 - NewNodeLocation )  )
//     //     {
//     //         std::set<unsigned> elements_containing_node1 = this->GetNode(ClosestEdge.first)->rGetContainingElementIndices();
//     //         std::set<unsigned> elements_containing_node3 = this->GetNode(ClosestEdge.second)->rGetContainingElementIndices();

//     //         std::set<unsigned> shared_elements;
//     //         std::set_intersection(elements_containing_node1.begin(),
//     //                                 elements_containing_node1.end(),
//     //                                 elements_containing_node3.begin(),
//     //                                 elements_containing_node3.end(),
//     //                                 std::inserter(shared_elements, shared_elements.begin()));
//     //         typename std::set<unsigned>::iterator CommonElements = shared_elements.begin();
//     //         unsigned Element1 = *CommonElements;
//     //         if (shared_elements.size() ==1)
//     //         {
//     //             // This is a boundary edge 
//     //              ClosestElement =Element1;
//     //              HaveElement=1;
//     //         }else 
//     //         {   
//     //             std::advance(CommonElements, 1);
//     //             unsigned Element2 = *CommonElements;

//     //             if (PointInTriangle3D(NewNodeLocation, Element1))
//     //             {
//     //                 ClosestElement =Element1;
//     //                 HaveElement=1;
//     //                 // std::cout<<" Inside element 1\n";
//     //                 // TRACE(" Inside element 1\n")
//     //             } else if( PointInTriangle3D(NewNodeLocation, Element2) ==1)
//     //             {
//     //                 ClosestElement =Element2;
//     //                 HaveElement=1;
    
//     //                 // TRACE(" Inside element 2\n")
//     //             }else
//     //             {  
//     //                 // TRACE("  Not inside the elements connected to the edge. Dam it\n")
//     //                  // unsigned Node1 = ClosestEdge.first;
//     //             // unsigned Node2 = ClosestEdge.second;
                    
//     //                 c_vector<double, SPACE_DIM> Centroid1 = mCentroidMap[Element1];
//     //                 c_vector<double, SPACE_DIM> Centroid2 = mCentroidMap[Element2];

//     //                 double Distance1 = norm_2(Centroid1 - NewNodeLocation);
//     //                 double Distance2 = norm_2(Centroid2 - NewNodeLocation);
//     //                 if (Distance1 <Distance2)
//     //                 {
//     //                     ClosestElement =Element1;
//     //                 }
//     //                 else{
//     //                     ClosestElement =Element2;

//     //                 }

//     //                 // // FIgure out where it is the same way I used as above in method 2 -- the current way I have doesnt work!!!!!!
//     //                 // c_vector<double, 3> Edge = P2 - P1;
//     //                 // c_vector<double, 3> EdgeMidpoint = (P2 + P1) / 2;

//     //                 // c_vector<double, SPACE_DIM> Normal1 = GetElementNormal(Element1);
//     //                 // c_vector<double, SPACE_DIM> Normal2 = GetElementNormal(Element2);

//     //                 // // Average normal of the two elements
//     //                 // c_vector<double, SPACE_DIM> Normal = (Normal1 + Normal2) / 2;

//     //                 // // The normal should be pointing outwards, so the direction of the plane normal will depend on the direction of the edge vector.
//     //                 // // The edge vector is totally depenedent on the order of the nodes pair, which should be sorted numerically
//     //                 // c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(Normal, Edge);
//     //                 // PlaneNormal /= norm_2(PlaneNormal);

//     //                 // // Need to determine the direction of the normal

//     //                 // c_vector<double, SPACE_DIM> NormalArm = PlaneNormal + EdgeMidpoint;
//     //                 // double DistanceToA = norm_2(NormalArm - mCentroidMap[Element1]); // A being the center point of element 1
//     //                 // double DistanceToB = norm_2(NormalArm - mCentroidMap[Element2]); // A being the center point of element 2
//     //                 // assert(DistanceToB != DistanceToA);
//     //                 // if (DistanceToA < DistanceToB)
//     //                 // {
//     //                 //     // TRACE("The normal is pointing towards element 1 ")
//     //                 //     // I now have a plane, Dependeing on what side of the plane the point is, will give me the closest element, but I dont know which way the plane is pointing
//     //                 //     double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//     //                 //     if (side > 0)
//     //                 //     {
//     //                 //         // TRACE("The new node is on the element ONE side of the midline plane ")
//     //                 //         ClosestElement = Element1;
//     //                 //         HaveElement=1;
//     //                 //     }
//     //                 //     else
//     //                 //     {
//     //                 //         // TRACE("The new node is on the element TWO side of the midline plane")
//     //                 //         ClosestElement = Element2;
//     //                 //         HaveElement=1;
//     //                 //     }

//     //                 // }  else if (DistanceToB < DistanceToA)
//     //                 // {
//     //                 // // TRACE("The normal is pointing towards element 2 ")
//     //                 // double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//     //                 // if (side > 0)
//     //                 // {
//     //                 //     // TRACE("The new node is on the element TWO side of the midline plane")
//     //                 //     ClosestElement = Element2;
//     //                 //     HaveElement=1;
//     //                 // }
//     //                 // else
//     //                 // {
//     //                 //     //   TRACE("The new node is on the element ONE side of the midline plane")
//     //                 //     ClosestElement = Element1;
//     //                 //     HaveElement=1;
//     //                 // }
//     //             }

                
               
//     //             // }
              
//     //         // //////////////////////////////v//////
//     //         //     // Edge has two elements
//     //         //    // Now determine which way the normal is faceing? -- Away from new point or towards? 
//     //         //     std::advance(CommonElements, 1);
//     //         //     unsigned Element2 = *CommonElements;
//     //         //     c_vector<double, SPACE_DIM> Centroid1 = mCentroidMap[Element1];
//     //         //     c_vector<double, SPACE_DIM> Centroid2 = mCentroidMap[Element2];

//     //         //      if (norm_2((ClosestMidPoint+ ClosestEdgeNorm) - Centroid1 )< norm_2((ClosestMidPoint+ ClosestEdgeNorm) - Centroid2)  && num_elements<=Element1)
//     //         //     {
//     //         //         // Element 1 is closer//
//     //         //         ClosestElement =Element1;
//     //         //         HaveElement=1;
//     //         //     }else 
//     //         //     {
//     //         //         // Element 2 is closer//
//     //         //         ClosestElement =Element2;
//     //         //         HaveElement=1;
//     //         //      }
//     //         }
//     //       mMapOfProbNodes[node_index] =3;
//     //     }
//     //     //So I Have the element this node is closest to. 
//     // }
//     // if (HaveElement == 0)
//     // {
//     //      mMapOfProbNodes[node_index] =4;
//     // }
//     assert(ClosestElement !=-10);
//     return ClosestElement;
// }



// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetClosestElementInOldMeshMethod5(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation)
// {

//     mMapOfProbNodes[node_index] = 5; // -5 meaning I cant find the needed element
//     unsigned ClosestElement = -10;
//     double distance =  std::max( (mMaxX - mMinX), (mMaxY - mMinY))/(2*mNx);

//     // Determine which bin I am in, then look in that bin for the element containing this node
//     c_vector<double,3> Bin = DetermineBin(NewNodeLocation);

// std::vector<double> CurrentBin = mBinMap[Bin];
    
//     bool HaveElement = 0;
//      for (std::vector<double>::iterator ele_iter = CurrentBin.begin(); ele_iter != CurrentBin.end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
//     {
//         // Now need to find the cloeset element in this bin
//         c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
//         if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
//         {
//             distance = abs(norm_2(NewNodeLocation - Centroid));
//             ClosestElement = *ele_iter;
//             if ( PointInTriangle3D(NewNodeLocation, ClosestElement))
//             {
//                 HaveElement =1;
//                 mMapOfProbNodes[node_index] = 1; 
//                 break;
//             }
//         }
//     }
   
//     unsigned num_elements = this->rGetMesh().GetNumElements();
    

//     if (HaveElement == 0)
//     {
//         distance = mMaxEdgelength;
//         std::pair<unsigned,unsigned> ClosestEdge;
//         c_vector<double, SPACE_DIM> ClosestEdgeNorm ;
//         c_vector<double, SPACE_DIM> ClosestMidPoint;
//          for (typename MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SpringIterator spring_iterator = this->SpringsBegin();
//          spring_iterator != this->SpringsEnd();
//          ++spring_iterator)
//          {
//             // mInitalPositionOfRemeshedNodes
//             unsigned NodeIndexA = spring_iterator.GetNodeA()->GetIndex();
//             unsigned NodeIndexB = spring_iterator.GetNodeB()->GetIndex();

//             c_vector<double, SPACE_DIM> P1 = this->GetNode(NodeIndexA)->rGetLocation();
//             c_vector<double, SPACE_DIM> P2 = this->GetNode(NodeIndexB)->rGetLocation();
//             c_vector<double, SPACE_DIM> EdgeVector = P2- P1;
//             c_vector<double, SPACE_DIM> MidPoint = (P2+ P1)/2;
//             EdgeVector/=norm_2(EdgeVector);
//             c_vector<double, SPACE_DIM> V = NewNodeLocation - P1;
 
//             // Normal to the plane containing the edge and the new node
//             c_vector<double, SPACE_DIM> N = VectorProduct(EdgeVector, V)/norm_2(V);
//             N/=norm_2(N);
//             // Now can find the normal to the edge that is contaned in the plane defined by the three points
            
//             c_vector<double, SPACE_DIM> normal = VectorProduct(EdgeVector, N);
//             normal/=norm_2(normal);
//             if (norm_2(MidPoint-normal -NewNodeLocation) < norm_2(MidPoint+normal -NewNodeLocation) )
//             {
//                 normal*=-1;
//             }
//             // Need to check the direction of the normal :) 

//             // normal to edge will be normal to the edge and the element normal -- se get the cross product of
//             double theta = acos(inner_prod(normal, V/norm_2(V)) );
//             double d =  norm_2(V) *cos(theta);
//             if (d< distance )
//             {
//                 distance =d;
//                 // need to think about what to do here, 
//                 ClosestEdge = std::pair<unsigned,unsigned>(NodeIndexA, NodeIndexB);
//                 ClosestEdgeNorm = normal;
//                 ClosestMidPoint = MidPoint;
//             }
            
//         }

//         c_vector<double, SPACE_DIM> P1 = this->GetNode(ClosestEdge.first)->rGetLocation();
//         c_vector<double, SPACE_DIM> P2 = this->GetNode(ClosestEdge.second)->rGetLocation();

//         //Check that this is close enough // how ??  -- distance to edge is smaller than distance to element center
//         if (distance < norm_2(P1 - NewNodeLocation )  || distance < norm_2(P2 - NewNodeLocation )  )
//         {
//             std::set<unsigned> elements_containing_node1 = this->GetNode(ClosestEdge.first)->rGetContainingElementIndices();
//             std::set<unsigned> elements_containing_node3 = this->GetNode(ClosestEdge.second)->rGetContainingElementIndices();

//             std::set<unsigned> shared_elements;
//             std::set_intersection(elements_containing_node1.begin(),
//                                     elements_containing_node1.end(),
//                                     elements_containing_node3.begin(),
//                                     elements_containing_node3.end(),
//                                     std::inserter(shared_elements, shared_elements.begin()));
//             typename std::set<unsigned>::iterator CommonElements = shared_elements.begin();
//             unsigned Element1 = *CommonElements;
//             if (shared_elements.size() ==1)
//             {
//                 // This is a boundary edge 
//                  ClosestElement =Element1;
//                  HaveElement=1;
//             }else 
//             {   
//                 std::advance(CommonElements, 1);
//                 unsigned Element2 = *CommonElements;

//                 if (PointInTriangle3D(NewNodeLocation, Element1)==1)
//                 {
//                     ClosestElement =Element1;
//                     HaveElement=1;
//                     // std::cout<<" Inside element 1\n";
//                     // TRACE(" Inside element 1\n")
//                 } else if( PointInTriangle3D(NewNodeLocation, Element2) ==1)
//                 {
//                     ClosestElement =Element2;
//                     HaveElement=1;
//                     // TRACE(" Inside element 2\n")
//                 }
//                 // The new node is now not in either of the elements, need to find which one it is closest to 
//                 else if ( ClosestPointInTriangle(NewNodeLocation, Element1) <= ClosestPointInTriangle(NewNodeLocation, Element2))
//                 {
//                     ClosestElement =Element1;
//                     HaveElement=1;
                    
//                 }else
//                 {
//                     ClosestElement =Element2;
//                     HaveElement=1;

//                 }

//             }
//         }
//     }

//                 // else
//                 // {  
//                 //     // TRACE("  Not inside the elements connected to the edge. Dam it\n")
//                 //      // unsigned Node1 = ClosestEdge.first;
//                 // // unsigned Node2 = ClosestEdge.second;
                    
//                 //     c_vector<double, SPACE_DIM> Centroid1 = mCentroidMap[Element1];
//                 //     c_vector<double, SPACE_DIM> Centroid2 = mCentroidMap[Element2];

//                 //     double Distance1 = norm_2(Centroid1 - NewNodeLocation);
//                 //     double Distance2 = norm_2(Centroid2 - NewNodeLocation);
//                 //     if (Distance1 <Distance2)
//                 //     {
//                 //         ClosestElement =Element1;
//                 //     }
//                 //     else{
//                 //         ClosestElement =Element2;

//                 //     }

//                     // // FIgure out where it is the same way I used as above in method 2 -- the current way I have doesnt work!!!!!!
//                     // c_vector<double, 3> Edge = P2 - P1;
//                     // c_vector<double, 3> EdgeMidpoint = (P2 + P1) / 2;

//                     // c_vector<double, SPACE_DIM> Normal1 = GetElementNormal(Element1);
//                     // c_vector<double, SPACE_DIM> Normal2 = GetElementNormal(Element2);

//                     // // Average normal of the two elements
//                     // c_vector<double, SPACE_DIM> Normal = (Normal1 + Normal2) / 2;

//                     // // The normal should be pointing outwards, so the direction of the plane normal will depend on the direction of the edge vector.
//                     // // The edge vector is totally depenedent on the order of the nodes pair, which should be sorted numerically
//                     // c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(Normal, Edge);
//                     // PlaneNormal /= norm_2(PlaneNormal);

//                     // // Need to determine the direction of the normal

//                     // c_vector<double, SPACE_DIM> NormalArm = PlaneNormal + EdgeMidpoint;
//                     // double DistanceToA = norm_2(NormalArm - mCentroidMap[Element1]); // A being the center point of element 1
//                     // double DistanceToB = norm_2(NormalArm - mCentroidMap[Element2]); // A being the center point of element 2
//                     // assert(DistanceToB != DistanceToA);
//                     // if (DistanceToA < DistanceToB)
//                     // {
//                     //     // TRACE("The normal is pointing towards element 1 ")
//                     //     // I now have a plane, Dependeing on what side of the plane the point is, will give me the closest element, but I dont know which way the plane is pointing
//                     //     double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//                     //     if (side > 0)
//                     //     {
//                     //         // TRACE("The new node is on the element ONE side of the midline plane ")
//                     //         ClosestElement = Element1;
//                     //         HaveElement=1;
//                     //     }
//                     //     else
//                     //     {
//                     //         // TRACE("The new node is on the element TWO side of the midline plane")
//                     //         ClosestElement = Element2;
//                     //         HaveElement=1;
//                     //     }

//                     // }  else if (DistanceToB < DistanceToA)
//                     // {
//                     // // TRACE("The normal is pointing towards element 2 ")
//                     // double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//                     // if (side > 0)
//                     // {
//                     //     // TRACE("The new node is on the element TWO side of the midline plane")
//                     //     ClosestElement = Element2;
//                     //     HaveElement=1;
//                     // }
//                     // else
//                     // {
//                     //     //   TRACE("The new node is on the element ONE side of the midline plane")
//                     //     ClosestElement = Element1;
//                     //     HaveElement=1;
//                     // }
//             // }

                
               
//                 // }
              
//             // //////////////////////////////v//////
//             //     // Edge has two elements
//             //    // Now determine which way the normal is faceing? -- Away from new point or towards? 
//             //     std::advance(CommonElements, 1);
//             //     unsigned Element2 = *CommonElements;
            
//             //     c_vector<double, SPACE_DIM> Centroid1 = mCentroidMap[Element1];
//             //     c_vector<double, SPACE_DIM> Centroid2 = mCentroidMap[Element2];

//             //      if (norm_2((ClosestMidPoint+ ClosestEdgeNorm) - Centroid1 )< norm_2((ClosestMidPoint+ ClosestEdgeNorm) - Centroid2)  && num_elements<=Element1)
//             //     {
//             //         // Element 1 is closer//
//             //         ClosestElement =Element1;
//             //         HaveElement=1;
            
//             //     }else 
//             //     {
//             //         // Element 2 is closer//
//             //         ClosestElement =Element2;
//             //         HaveElement=1;
//             //      }
//         //     // }
//         //   mMapOfProbNodes[node_index] =3;
//         // }
//         //So I Have the element this node is closest to. 
//     // }
    
//     if (HaveElement == 0)
//     {
//          mMapOfProbNodes[node_index] =4;
//     }
//     assert(ClosestElement !=-10);
//     return ClosestElement;
// }



// std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > UnitNormals;
// std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> otherNodes;

// std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edgey = std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>(this->GetNode(ClosestEdge.first), this->GetNode(ClosestEdge.second));

// bool boundary_edge_found = CalculateElementNormals(edgey, UnitNormals, otherNodes);
// if (boundary_edge_found ==1)
// {
//     // PRINT_VARIABLE(boundary_edge_found)
//     // PRINT_4_VARIABLES(num_elements, Element1, Element2, shared_elements.size())
//     ClosestElement =Element1;
//     HaveElement=1;
    
// }




// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetClosestElementInOldMesh(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation)
// {
//     mMapOfProbNodes[node_index] = 5; // -5 meaning I cant find the needed element
//     unsigned ClosestElement = -10;

//     // Determine which bin I am in, then look in that bin for the element containing this node
//     c_vector<double,3> Bin = DetermineBin(NewNodeLocation);
//     //  std::cout << Bin[0] << " " << Bin[1] << " " << Bin[2] << std::endl;
//     unsigned num_elements = this->rGetMesh().GetNumElements();
//     assert(mNx >= Bin[0] && mNy >= Bin[1] && mNz >= Bin[2]);

//     bool HaveElement = 0;

//      std::vector<double> CurrentBin = mBinMap[Bin];
//     for (std::vector<double>::iterator ele_iter = CurrentBin.begin(); ele_iter != CurrentBin.end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
//     {
//         // Now need to find the element containing this node in this bin
//         bool InTriangle = PointInTriangle3D(NewNodeLocation, *ele_iter);
//         if (InTriangle)
//         {
//             // TRACE("found nearest element by looping over bin and checking if in each element element")
//             ClosestElement = *ele_iter;
//             HaveElement = 1; // If I have the element I can end this search here :) --- I might wrap this search up to its own funciton

//             mMapOfProbNodes[node_index] = 0;
//             break;
//         }
//     }
//     // If I havent found the element by searching the closest bin, then look across whole mesh 

//     if (HaveElement == 0)
//     {
//         for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
//              elem_iter != this->rGetMesh().GetElementIteratorEnd();
//              ++elem_iter)
//         {
//             bool InTriangle = PointInTriangle3D(NewNodeLocation, elem_iter->GetIndex());
//             if (InTriangle)
//             {
//                 ClosestElement = elem_iter->GetIndex(); //TRACE(" Have found the nearest element by looping over every element and checking if it is in the element")
               
//                 HaveElement = 1;
//                 mMapOfProbNodes[node_index] = 1;
//                 break;
//             }
//         }
//     }

//     if (HaveElement == 0)
//     {

//         double distance = 1e3; //
//         double distance2 = 1e3; //

// std::vector<double> CurrentBin = mBinMap[Bin];
//         unsigned ClosestElement2 = -10;
//         // Now need to loop over this bin to find the nearest element centroid  -- iterate over vector -- getting element centorid indices out -- can go straight to the centroid here
//         for (std::vector<double>::iterator ele_iter = CurrentBin.begin(); ele_iter != CurrentBin.end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
//         {
//             // Now need to find the cloeset element in this bin
//             c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
//             if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
//             {
//                 distance = abs(norm_2(NewNodeLocation - Centroid));
//                 ClosestElement = *ele_iter;
//             }
//         }
//         bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
//         if (InTriangle)
//         {
//             mMapOfProbNodes[node_index] = 2;
//         }
//         else if (InTriangle == 0)
//         {
//             // TRACE("Not in triganlge")

//             // Now I have the 'closest element', get the neighbours to this element, and check the node isnt actually in one of them
//             std::set<unsigned> NeighbouringElements = GetNeighbouringElements(ClosestElement);
//             for (typename std::set<unsigned>::iterator neigh_iter = NeighbouringElements.begin();
//                  neigh_iter != NeighbouringElements.end();
//                  ++neigh_iter)
//             {
//                 if (*neigh_iter != ClosestElement)
//                 {
//                     c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*neigh_iter];

//                     if (abs(norm_2(NewNodeLocation - Centroid)) < distance2)
//                     {
//                         distance2 = abs(norm_2(NewNodeLocation - Centroid));
//                         ClosestElement2 = *neigh_iter;
//                     }
//                 }
//             }
//             // TRACE("Not in the second cloeset element either, need to figure out which side of the plane containing the edge")
//             std::vector<unsigned> CommonNodes = GetCommonNodes(ClosestElement, ClosestElement2);
//             assert(CommonNodes.size() < 3);
//             if (CommonNodes.size() == 2)
//             {
                
//                 typename std::vector<unsigned>::iterator CommonNodeIter = CommonNodes.begin();

//                 Node<SPACE_DIM>* pNode1 = this->rGetMesh().GetNode(*CommonNodeIter);
//                 c_vector<double, SPACE_DIM> P1 = pNode1->rGetLocation();
//                 std::advance(CommonNodeIter, 1);

//                 Node<SPACE_DIM>* pNode2 = this->rGetMesh().GetNode(*CommonNodeIter);
//                 c_vector<double, SPACE_DIM> P2 = pNode2->rGetLocation();

//                 c_vector<double, 3> Edge = P2 - P1;
//                 c_vector<double, 3> EdgeMidpoint = (P2 + P1) / 2;

//                 c_vector<double, SPACE_DIM> Normal1 = GetElementNormal(ClosestElement);
//                 c_vector<double, SPACE_DIM> Normal2 = GetElementNormal(ClosestElement2);

//                 // Average normal of the two elements
//                 c_vector<double, SPACE_DIM> Normal = (Normal1 + Normal2) / 2;

//                 // The normal should be pointing outwards, so the direction of the plane normal will depend on the direction of the edge vector.
//                 // The edge vector is totally depenedent on the order of the nodes pair, which should be sorted numerically
//                 c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(Normal, Edge);
//                 PlaneNormal /= norm_2(PlaneNormal);

//                 // Need to determine the direction of the normal

//                 c_vector<double, SPACE_DIM> NormalArm = PlaneNormal + EdgeMidpoint;
//                 double DistanceToA = norm_2(NormalArm - mCentroidMap[ClosestElement]); // A being the center point of element 1
//                 double DistanceToB = norm_2(NormalArm - mCentroidMap[ClosestElement2]); // A being the center point of element 2
//                 assert(DistanceToB != DistanceToA);

//                 if (DistanceToA < DistanceToB)
//                 {

//                     // TRACE("The normal is pointing towards element 1 ")
//                     // I now have a plane, Dependeing on what side of the plane the point is, will give me the closest element, but I dont know which way the plane is pointing
//                     double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//                     if (side > 0)
//                     {
//                         // TRACE("The new node is on the element ONE side of the midline plane ")
//                     }
//                     else
//                     {
//                         // TRACE("The new node is on the element TWO side of the midline plane")
//                         ClosestElement = ClosestElement2;
//                     }
                 
//                 }
//                 else if (DistanceToB < DistanceToA)
//                 {
//                     // TRACE("The normal is pointing towards element 2 ")
//                     double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//                     if (side > 0)
//                     {
//                         // TRACE("The new node is on the element TWO side of the midline plane")
//                         ClosestElement = ClosestElement2;
//                     }
//                     else
//                     {
//                         //   TRACE("The new node is on the element ONE side of the midline plane")
//                     }
               
//                 }
//                 mMapOfProbNodes[node_index] = 3;
//                 HaveElement = 1;
//             }
//         }
//     }
//     assert(ClosestElement != -10);
//     if (HaveElement == 0)
//     {
//         // I think this is the cases I cant find
//         mMapOfProbNodes[node_index] = 5; // -5 meaning I cant find the needed element
//         // If I havent been able to find the containing or the closest element, then the defult is the closest
//     }
//     return ClosestElement;
// }




// // Now need to loop over this bin to find the nearest element centroid  -- iterate over vector -- getting element centorid indices out -- can go straight to the centroid here
// for (std::vector<double>::iterator ele_iter = CurrentBin.begin(); ele_iter != mBinMap[Bin].end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
// {

//     // Now need to find the cloeset element in this bin
//     c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
//     if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
//     {
//         distance = abs(norm_2(NewNodeLocation - Centroid));
//         ClosestElement = *ele_iter;
//     }
// }
// assert(num_elements >= ClosestElement && ClosestElement != -10);
// bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
// if (InTriangle ==1 )
// {
//     A+=1;
// }
// if (InTriangle ==0 )
// {
//     // TRACE("Not in triganlge")

//     // Now I have the 'closest element', get the neighbours to this element, and check the node isnt actually in one of them
//     std::set<unsigned> NeighbouringElements = GetNeighbouringElements(ClosestElement);
//     for (typename std::set<unsigned>::iterator neigh_iter = NeighbouringElements.begin();
//         neigh_iter != NeighbouringElements.end();
//         ++neigh_iter)
//     {
//         if (*neigh_iter != ClosestElement)
//         {
//             c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*neigh_iter];

//             if (abs(norm_2(NewNodeLocation - Centroid)) < distance2)
//             {
//                 distance2 = abs(norm_2(NewNodeLocation - Centroid));
//                 ClosestElement2 = *neigh_iter;
//             }
//         }
//     }

//     bool InTriangle2 = PointInTriangle3D(NewNodeLocation, ClosestElement);
//     if  (InTriangle2==1)
//     {
//         // TRACE("In the second triangle")// Looks like this never happens?
//         B+=1;
//     }
//     if (InTriangle2==0)
//     {
//         // TRACE("Not in the second cloeset element either, need to figure out which side of the plane containing the edge")
//         std::vector<unsigned> CommonNodes = GetCommonNodes(ClosestElement, ClosestElement2);
//         assert(CommonNodes.size() < 3);
//         if (CommonNodes.size() == 2)
//         {
//             // TRACE("Two common nodes")
//             typename std::vector<unsigned>::iterator CommonNodeIter = CommonNodes.begin();

//             Node<SPACE_DIM>* pNode1 = this->rGetMesh().GetNode(*CommonNodeIter);
//             c_vector<double, SPACE_DIM> P1 = pNode1->rGetLocation();
//             std::advance(CommonNodeIter, 1);

//             Node<SPACE_DIM>* pNode2 = this->rGetMesh().GetNode(*CommonNodeIter);
//             c_vector<double, SPACE_DIM> P2 = pNode2->rGetLocation();

//             c_vector<double, 3> Edge = P2 - P1;
//             c_vector<double, 3> EdgeMidpoint = (P2 + P1) / 2;

//             c_vector<double, SPACE_DIM> Normal1 = GetElementNormal(ClosestElement);
//             c_vector<double, SPACE_DIM> Normal2 = GetElementNormal(ClosestElement2);

//             // Average normal of the two elements
//             c_vector<double, SPACE_DIM> Normal = (Normal1 + Normal2) / 2;

//             // The normal should be pointing outwards, so the direction of the plane normal will depend on the direction of the edge vector.
//             // The edge vector is totally depenedent on the order of the nodes pair, which should be sorted numerically
//             c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(Normal, Edge);
//             PlaneNormal /= norm_2(PlaneNormal);

//             // Need to determine the direction of the normal

//             c_vector<double, SPACE_DIM> NormalArm = PlaneNormal + EdgeMidpoint;
//             double DistanceToA = norm_2(NormalArm -mCentroidMap[ClosestElement]); // A being the center point of element 1
//             double DistanceToB = norm_2(NormalArm -mCentroidMap[ClosestElement2]); // A being the center point of element 2
//             assert(DistanceToB != DistanceToA);

//             if (DistanceToA < DistanceToB)
//             {

//                 // TRACE("The normal is pointing towards element 1 ")
//                 // I now have a plane, Dependeing on what side of the plane the point is, will give me the closest element, but I dont know which way the plane is pointing
//                 double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//                 if (side>0)
//                 {
//                     // TRACE("The new node is on the element ONE side of the midline plane ")
//                 }
//                 else
//                 {
//                         // TRACE("The new node is on the element TWO side of the midline plane")
//                         ClosestElement = ClosestElement2;
//                 }
//             C+=1;

//             }else if  (DistanceToB < DistanceToA)
//             {
//                 // TRACE("The normal is pointing towards element 2 ")
//                 double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//                 if (side>0)
//                 {
//                     // TRACE("The new node is on the element TWO side of the midline plane")
//                     ClosestElement = ClosestElement2;
//                 }
//                 else
//                 {
//                     //   TRACE("The new node is on the element ONE side of the midline plane")
//                 }
//             C+=1;

//             }

//         } else if (CommonNodes.size() <2)
//         {


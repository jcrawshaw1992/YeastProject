#include "HistoryDepMeshBasedCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::HistoryDepMeshBasedCellPopulation(MutableMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                                                                             std::vector<CellPtr>& rCells,
                                                                                             const std::vector<unsigned> locationIndices,
                                                                                             bool deleteMesh,
                                                                                             bool validate)
        : MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>(rMesh, rCells, locationIndices, deleteMesh, validate)
{
    SaveInitalConditions();
    this->SetBinningRegions();
    // SetBinningIntervals(1, 1, 1);
    bool InitialRemesh =0;
    if (InitialRemesh)
    {
        RemeshGeometry();
    }
  
    if (mSetUpInitialConfigurations == 1)
    {
        MarkBoundaryNodes();
        // Define all the inital configurations ... this needed to be done and the start, and whenever there is a remeshing
        SetupMembraneConfiguration();
        SetInitialAnlgesAcrossMembrane();
        // SetMaxEdgelength();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetUpInitialConfig(bool SetUpInitialConfigurations)
{
    mSetUpInitialConfigurations = SetUpInitialConfigurations;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::HistoryDepMeshBasedCellPopulation(MutableMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
        : MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>(rMesh)
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::~HistoryDepMeshBasedCellPopulation()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetBoundaries(bool SetB)
{
    mSetBoundaries = SetB;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::ExecuteHistoryDependentRemeshing()
{

    TRACE("ExecuteHistoryDependentRemeshing")
    /*
     * 1) Remesh geometry
     * 2) Map
     * 3) Delete old mesh and reset with new mesh
     * 4) detele old cells and relpace with the cells for the current mesh  
	 */
    TRACE("Remesh Geometry")
    if (mRemeshingSoftwear == "VMTK")
    {
        this->RemeshGeometryWithVMTK();
    }
    else if (mRemeshingSoftwear == "CGAL")
    {
        TRACE("REmeshing")
        this->RemeshGeometry();
    }
    else if (mRemeshingSoftwear == "PreAllocatedMatlabMesh")
    {
        this->TakeInPreAllocatedRemeshGeometry();
    }
    TRACE("SetBinningRegions")
    this->SetBinningRegions();
    TRACE("MappingAdaptedMeshToInitalGeometry")
    this->MappingAdaptedMeshToInitalGeometry();
    if (mPrintRemeshedIC)
    {
        TRACE("WriteOutMappedInitalConfig")
        WriteOutMappedInitalConfig();
    }
    


    VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(mRelativePath, "RemeshedGeometryAdjusted", false);
    mesh_writer.WriteFilesUsingMesh(mNew_mesh);//mNew_mesh

    TRACE("TempTurnOff")
    /*
        --------------------------------
              The Old switcharoo 
        
        Now Switch the old and new meshes
        ---------------------------------
     */

    assert(ELEMENT_DIM == 2 && SPACE_DIM == 3);
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
    SetInitialAnlgesAcrossMembrane();

    // Initialise the applied force at each node to zero
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        node_iter->ClearAppliedForce();
        CellPtr p_cell = this->GetCellUsingLocationIndex(node_iter->GetIndex());
        p_cell->GetCellData()->SetItem("MappingMethod", mMapOfProbNodes[node_iter->GetIndex()]);
    }
    this->SetBinningRegions();//Remove this later :) 

    mUpdateComplete = 0;
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

    for (unsigned i = 0; i < MappedICmesh.GetNumNodes(); i++)
    {
        Node<SPACE_DIM>* p_node = MappedICmesh.GetNode(i);
        p_node->rGetModifiableLocation() = mInitalPositionOfRemeshedNodes[i];
        // PRINT_VECTOR(mInitalPositionOfRemeshedNodes[i])
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
    int SystemOutput;
    // TRACE("Remesh geometry here ")
    // Outpust the mesh in .vtu format for HemeLB setup tool to pick up (first converted to stl, though).
    VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(mRelativePath, "config", false);
    mesh_writer.WriteFilesUsingMesh(this->rGetMesh());


    if (mVariableEdgeLength ==1)
    {
        if (mTargetRemeshingEdgeLength> 0.9*1e-2)
        {
            mVariableEdgeLength =0;
        }
        else{
            mTargetRemeshingEdgeLength *= mEdgeLengthMultiple;
        }
         
    }
  
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
    // PRINT_VARIABLE(mChasteOutputDirectory)
    std::string vtu2offCommand = "meshio-convert " + mChasteOutputDirectory + "config.vtu " + offfile;
    SystemOutput = std::system(vtu2offCommand.c_str()); // system only takes char *.
    // cmake -DCGAL_DIR=$HOME/CGAL-5.2.1 -DCMAKE_BUILD_TYPE=Release .
    // make isotropic_remeshing_ForChaste
    // Now excute the CGAL command to remesh the current geometry - not the input and output within this file have to be pre-set. I will explore if I can make this more neat later, should care.... dont care
    std::string CGALRemeshingCommand;
    if (mServer ==1)
    {
        CGALRemeshingCommand = "(cd  /home/vascrem/CGAL-5.0.2/Polygon_mesh_processing/examples/Polygon_mesh_processing/;./isotropic_remeshing_ForChaste -input " + offfile + " -output " + mChasteOutputDirectory + "CurrentPlexusRemeshed.off -target_edge_length " + std::to_string(mTargetRemeshingEdgeLength) + " -iterations " + std::to_string(mIterations) + " > null)";
        SystemOutput = std::system(CGALRemeshingCommand.c_str()); // system only takes char *
    } 
    else if (mServer ==0)
    {
        CGALRemeshingCommand = "(cd  ~/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/;./isotropic_remeshing_ForChaste -input " + offfile + " -output " + mChasteOutputDirectory + "CurrentPlexusRemeshed.off -target_edge_length " + std::to_string(mTargetRemeshingEdgeLength) + " -iterations " + std::to_string(mIterations) + " > null)";
        SystemOutput = std::system(CGALRemeshingCommand.c_str()); // system only takes char *
    } 

    // Now ned to convert from .off back to a .vtu
    // TRACE("Now ned to convert from .off back to a .vtu")
    std::string Remeshedvtu = mChasteOutputDirectory + "RemeshedGeometry.vtu";
    std::string off2vtuCommand = "meshio-convert " + mChasteOutputDirectory + "CurrentPlexusRemeshed.off " + Remeshedvtu;
    SystemOutput = std::system(off2vtuCommand.c_str()); // system only takes char *

    // /// Delete the orignal .off file -- else if the CurrentMesh.off fails to transfer, the old transfer will be the one taken
    // std::string DeleteOldCurrentConfigOff = "rm ~/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/data/CurrentMesh.off";
    // SystemOutput = std::system( DeleteOldCurrentConfigOff.c_str()); // system only takes char *
    // TRACE("Read in the new mesh")
    // Read in the new remeshsed mesh
    VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mChasteOutputDirectory + "RemeshedGeometry.vtu");
    mNew_mesh.ConstructFromMeshReader(mesh_reader);



    // double LowerAspectRatio = GetAspectRatioFromMesh();
    // double NewIterations = mIterations;
    // double NewTargetRemeshingEdgeLength = mTargetRemeshingEdgeLength;
    // double counter = 1;
    // PRINT_VARIABLE(LowerAspectRatio)
    // while( LowerAspectRatio <0.7 )
    // {
    //     TRACE("Need to remesh again!")

    //     NewTargetRemeshingEdgeLength*=0.95;
    //     NewIterations +=5;

    //     offfile = mChasteOutputDirectory +"CurrentPlexusRemeshed.off";
    //     if (mServer ==1)
    //     {
    //         CGALRemeshingCommand = "(cd  /home/vascrem/CGAL-5.0.2/Polygon_mesh_processing/examples/Polygon_mesh_processing/;./isotropic_remeshing_ForChaste -input " + offfile + " -output " + mChasteOutputDirectory + "CurrentPlexusRemeshed.off -target_edge_length " + std::to_string(NewTargetRemeshingEdgeLength) + " -iterations " + std::to_string(NewIterations) + " ) > null";
    //         SystemOutput = std::system(CGALRemeshingCommand.c_str()); // system only takes char *
    //     } 
    //     else if (mServer ==0)
    //     {
    //         CGALRemeshingCommand = "(cd  ~/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/;./isotropic_remeshing_ForChaste -input " + offfile + " -output " + mChasteOutputDirectory + "CurrentPlexusRemeshed.off -target_edge_length " + std::to_string(NewTargetRemeshingEdgeLength) + " -iterations " + std::to_string(NewIterations) + " > null) > null";
    //         SystemOutput = std::system(CGALRemeshingCommand.c_str()); // system only takes char *
    //     } 

    //     // Now ned to convert from .off back to a .vtu
    //     // TRACE("Now ned to convert from .off back to a .vtu")
    //     std::string Remeshedvtu = mChasteOutputDirectory + "RemeshedGeometry.vtu";
    //     std::string off2vtuCommand = "meshio-convert " + mChasteOutputDirectory + "CurrentPlexusRemeshed.off " + Remeshedvtu;
    //     SystemOutput = std::system(off2vtuCommand.c_str()); // system only takes char *

    //     std::stringstream outN;
    //     outN <<  counter;
    //     std::string MeshNumber = outN.str();

    //     std::string RemeshedvtuAgain = mChasteOutputDirectory + "RemeshedGeometry"+MeshNumber+".vtu";
    //     off2vtuCommand = "meshio-convert " + mChasteOutputDirectory + "CurrentPlexusRemeshed.off " + RemeshedvtuAgain;
    //     SystemOutput = std::system(off2vtuCommand.c_str()); // system only takes char *
    //     counter+=1;


    //     LowerAspectRatio = GetAspectRatioFromMesh();

    //     PRINT_2_VARIABLES(LowerAspectRatio,NewIterations)

    // }
    // TRACE("New Mesh Made")
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetAspectRatioFromMesh()
{
   
    assert(SPACE_DIM == 3);
    assert(ELEMENT_DIM == 2);
    std::vector<double> AspectRatioVector;
    double MinAspectRatio = 100;
    // Loop over the old map and get the centroids of the old map
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = mNew_mesh.GetElementIteratorBegin();
         elem_iter != mNew_mesh.GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        Node<SPACE_DIM>* pNode0 = mNew_mesh.GetNode(elem_iter->GetNodeGlobalIndex(0));
        Node<SPACE_DIM>* pNode1 = mNew_mesh.GetNode(elem_iter->GetNodeGlobalIndex(1));
        Node<SPACE_DIM>* pNode2 = mNew_mesh.GetNode(elem_iter->GetNodeGlobalIndex(2));

        double AspectRatio = CalculateAspectRatio(pNode0->rGetLocation(), pNode1->rGetLocation(), pNode2->rGetLocation() );
        // PRINT_VARIABLE(AspectRatio)
        if (AspectRatio < MinAspectRatio)
        {
            MinAspectRatio = AspectRatio;
            PRINT_VARIABLE(AspectRatio)
        }

    }

    return MinAspectRatio;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::RemeshGeometryWithVMTK()
{
    int SystemOutput;
    // old, when using vmtk and stls
    /*
    * Here I will remesh the current chaste geometery
    * 1) Take the latest .vtu file convert it into an .stl 
    * 2) Remesh this stl
    * 3) Convert this stl into a .vtu so I can read it into chaste and use it  
    */

    // Outpust the mesh in .vtu format for HemeLB setup tool to pick up (first converted to stl, though).
    // std::string stlfile = mChasteOutputDirectory;
    VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(mRelativePath, "config", false);
    mesh_writer.WriteFilesUsingMesh(this->rGetMesh());

    // Convert the .vtu into an .off
    std::string stlfile = mChasteOutputDirectory + "CurrentMesh.stl";
    // TRACE("Create .stl")
    std::string vtu2stlCommand = "meshio-convert " + mChasteOutputDirectory + "config.vtu " + stlfile;
    SystemOutput = std::system(vtu2stlCommand.c_str()); // system only takes char *.

    // Remesh the .stl
    std::string Remeshedstl = mChasteOutputDirectory + "RemeshedGeometry.stl";
    // TRACE("Remeshing geometry")

    // More info for the remeshing options found at https://sourceforge.net/p/vmtk/mailman/message/29391820/ and  http://www.vmtk.org/vmtkscripts/vmtksurfaceremeshing.html
    // std::string RemeshCommand = "vmtksurfaceremeshing -ifile " + stlfile + " -iterations " + std::to_string(mIterations) + " -area " + std::to_string(mTargetRemeshingElementArea) + " -maxarea "+ std::to_string(1.2*mTargetRemeshingElementArea) +" -ofile " + Remeshedstl;
    // SystemOutput = std::system(RemeshCommand.c_str());


    std::string RemeshCommand = "vmtksurfaceremeshing -ifile " + stlfile + " -iterations " + std::to_string(mIterations) + " -edgelength " + std::to_string(mTargetRemeshingEdgeLength) + " -ofile " + Remeshedstl + " -elementsizemode edgelength -internalangletolerance 0.3 -relaxation 0.6";
    SystemOutput = std::system(RemeshCommand.c_str());

    PRINT_VARIABLE(RemeshCommand)

    //Finally conver the Remeshed stl into a vtu -- meshio is your friend
    std::string Remeshedvtu = mChasteOutputDirectory + "RemeshedGeometry.vtu";
    std::string stl2vtuCommand = " meshio-convert " + Remeshedstl + " " + Remeshedvtu;
    SystemOutput = std::system(stl2vtuCommand.c_str());

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
    // PRINT_VARIABLE(mChasteOutputDirectory)
    // bool copy_file(mPreAllocatedRemeshedMesh.c_str(), Remeshedvtu.c_str(), NULL, COPYFILE_DATA | COPYFILE_XATTR);
    boost::filesystem::copy_file(mPreAllocatedRemeshedMesh.c_str(), Remeshedvtu.c_str()); // This was changed so the code could work on linux and mac
    // TRACE("Read in the new mesh")
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





// Binning methods  -- this needs to be simplified

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetBinningIntervals(int nX, int nY, int nZ)
{
    mNx = nX; mNy = nY; mNz = nZ;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetMeshSize()
{
    // Want to get the max and min dims in x, y and z
    // Get some first place, then find what is smaller and larger in each direction 
    c_vector<double, 3> Location_0 = (this->rGetMesh().GetNodeIteratorBegin())->rGetLocation();
    mMaxX = Location_0[0]; mMinX = Location_0[0];
    mMaxY = Location_0[1]; mMinY = Location_0[1];
    mMaxZ = Location_0[2]; mMinZ = Location_0[2];

    // Start by finding the max x, y, z locations
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        c_vector<double, 3> Location = node_iter->rGetLocation();
        // X directions
        if (Location[0] < mMinX)
        {
         mMinX = Location[0];
        } else if (Location[0] > mMaxX)
        {
         mMaxX = Location[0];
        }
        // Y directions
        if (Location[1] < mMinY)
        {
         mMinY = Location[1];
        } else if (Location[1] > mMaxY)
        {
            mMaxY = Location[1];
        }
        // Z directions
        if (Location[2] < mMinZ)
        {
            mMinZ = Location[2];
        } else if (Location[2] > mMaxZ)
        {
            mMaxZ = Location[2];
        }
    }

    double Xdom = (mMaxX - mMinX)/20;
    double Ydom = (mMaxY - mMinY)/20;
    double Zdom = (mMaxZ - mMinZ)/20;
    mMaxX += Xdom; mMaxY += Ydom;  mMaxZ += Zdom;
    mMinX -= Xdom; mMinY -= Ydom;  mMinZ -= Zdom;

}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetBinningWidth()
{
    // I need the maximal dimensions   
    double X = mMaxX - mMinX;
    double Y = mMaxY - mMinY;
    double Z = mMaxZ - mMinZ;
    
     double NumberOfBins = mNz*mNy*mNz;
     for (int i=0; i<mNx+1;++i)
     {
          for (int j=0; j<mNy+1;++j)
            {
                for (int k=0; k<mNz+1;++k)
                {
                    std::vector<int> BinIdentifier = { i,j,k };
                    double Xinterval_L = mMinX +(X*(i)/mNx);
                  
                    double Xinterval_U = mMinX +(X*(i+1)/mNx) ;

                    double Yinterval_L = mMinY +(Y*j/mNy);
                    double Yinterval_U = mMinY +(Y*(j+1)/mNy);

                    double Zinterval_L = mMinZ +(Z*k/mNz);
                    double Zinterval_U = mMinZ +(Z*(k+1)/mNz);
                    mBinCoords[BinIdentifier] = { Xinterval_L,Xinterval_U, Yinterval_L,Yinterval_U , Zinterval_L , Zinterval_U};
                    // PRINT_VECTOR(mBinCoords[BinIdentifier]);
                }
            }
         
     }
     
    
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetBinningRegions()
{

    assert(SPACE_DIM == 3 && ELEMENT_DIM == 2 );

    // Setting the binning regions of the old map here 
    this->SetMeshSize();
    this->SetBinningWidth();
    this->SetCentroidMap();
    // Clear the map so I dont have the old elements and the new elements messing the binning up
    mBin.clear();
    mEdgeBin.clear();
    double BlurryRegion = 10;
    double X = (mMaxX - mMinX)/BlurryRegion;
    double Y = (mMaxY - mMinY)/BlurryRegion;
    double Z = (mMaxZ - mMinZ)/BlurryRegion;
    

    // Need to iterate over the elements and determine which bin each element centroid is in. I look for the cloesest old centeroid for each new node, so it fits that I will have the centroids sorted in the bins
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
         elem_iter != this->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        c_vector<double, SPACE_DIM> Centroid = mCentroidMap[elem_index];

        for (int i=0; i<mNx+1;++i)
        {
            for (int j=0; j<mNy+1;++j)
            {
                for (int k=0; k<mNz+1;++k)
                {
                    // std::vector<double> BinLimits =mBinCoords[{ i,j,k }];
                    double MinX = GetBinLowerX(mBinCoords[{ i,j,k }])-X ; double MaxX = GetBinUpperX(mBinCoords[{ i,j,k }])+X ;
                    double MinY = GetBinLowerY(mBinCoords[{ i,j,k }])-Y; double MaxY = GetBinUpperY(mBinCoords[{ i,j,k }])+Y;
                    double MinZ = GetBinLowerZ(mBinCoords[{ i,j,k }])-Z; double MaxZ = GetBinUpperZ(mBinCoords[{ i,j,k }])+Z;

                    // Check if this element is in this bin (can be in multiple bins )
                    if (Centroid[0]>= MinX &&  Centroid[0]<= MaxX )
                      {if (Centroid[1]>= MinY &&  Centroid[1]<= MaxY )
                        { if (Centroid[2]>= MinZ &&  Centroid[2]<= MaxZ )
                            {
                                mBin[{ i,j,k }].push_back(elem_index);
                                // Label everything for visulisation sake
                                 for (int l= 0; l < 3; l++)
                                    {
                                        CellPtr p_cell = this->GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(l));
                                        p_cell->GetCellData()->SetItem("i", i);
                                        p_cell->GetCellData()->SetItem("j", j);
                                        p_cell->GetCellData()->SetItem("k", k);
                                        p_cell->GetCellData()->SetItem("BIN", i+j+k);
                                    }
                            }
                        }
                    }
                }
            }
        }
    }


    // Also need to get bin for edges

    for (typename MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SpringIterator spring_iterator = this->SpringsBegin();
            spring_iterator != this->SpringsEnd();
            ++spring_iterator)
    {
        unsigned NodeIndexA = spring_iterator.GetNodeA()->GetIndex();
        unsigned NodeIndexB = spring_iterator.GetNodeB()->GetIndex();

        c_vector<double, SPACE_DIM> MidPoint = 0.5*(this->GetNode(NodeIndexA)->rGetLocation()+this->GetNode(NodeIndexB)->rGetLocation());

        for (int i=0; i<mNx+1;++i)
        {
            for (int j=0; j<mNy+1;++j)
            {
                for (int k=0; k<mNz+1;++k)
                {
                    double MinX = GetBinLowerX(mBinCoords[{ i,j,k }])-X ; double MaxX = GetBinUpperX(mBinCoords[{ i,j,k }])+X ;
                    double MinY = GetBinLowerY(mBinCoords[{ i,j,k }])-Y; double MaxY = GetBinUpperY(mBinCoords[{ i,j,k }])+Y;
                    double MinZ = GetBinLowerZ(mBinCoords[{ i,j,k }])-Z; double MaxZ = GetBinUpperZ(mBinCoords[{ i,j,k }])+Z;

                    // Check if this element is in this bin (can be in multiple bins )
                    if (MidPoint[0]>= MinX &&  MidPoint[0]<= MaxX )
                      { if (MidPoint[1]>= MinY &&  MidPoint[1]<= MaxY )
                         { if (MidPoint[2]>= MinZ &&  MidPoint[2]<= MaxZ )
                            {
                                mEdgeBin[{ i,j,k }].push_back(std::make_pair(NodeIndexA, NodeIndexB));  
                                    CellPtr p_cell = this->GetCellUsingLocationIndex(NodeIndexA);
                                    p_cell->GetCellData()->SetItem("Edgei", i);
                                    p_cell->GetCellData()->SetItem("Edgej", j);
                                    p_cell->GetCellData()->SetItem("Edgek", k);
                                    p_cell->GetCellData()->SetItem("EdgeBIN", i+j+k);   

                                    CellPtr p_cell2 = this->GetCellUsingLocationIndex(NodeIndexB);
                                    p_cell2->GetCellData()->SetItem("Edgei", i);
                                    p_cell2->GetCellData()->SetItem("Edgej", j);
                                    p_cell2->GetCellData()->SetItem("Edgek", k);
                                    p_cell2->GetCellData()->SetItem("EdgeBIN", i+j+k);  
                            }
                         }
                    }
                }
            }
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<int> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetBin(c_vector<double, SPACE_DIM> Location)
{
    assert(SPACE_DIM == 3 &&ELEMENT_DIM == 2 );

    std::vector<int> Bin;

    // Need to iterate over the elements and determine which bin each element centroid is in. I look for the cloesest old centeroid for each new node, so it fits that I will have the centroids sorted in the bins

        for (int i=0; i<mNx+1;++i)
        {
            for (int j=0; j<mNy+1;++j)
            {
                for (int k=0; k<mNz+1;++k)
                {
                    double MinX = GetBinLowerX(mBinCoords[{ i,j,k }]); double MaxX = GetBinUpperX(mBinCoords[{ i,j,k }]);
                    double MinY = GetBinLowerY(mBinCoords[{ i,j,k }]); double MaxY = GetBinUpperY(mBinCoords[{ i,j,k }]);
                    double MinZ = GetBinLowerZ(mBinCoords[{ i,j,k }]); double MaxZ = GetBinUpperZ(mBinCoords[{ i,j,k }]);
                    // Check if this element is in this bin (can be in multiple bins )
                    if (Location[0]>= MinX &&  Location[0]<= MaxX )
                      {if (Location[1]>= MinY &&  Location[1]<= MaxY )
                        { if (Location[2]>= MinZ &&  Location[2]<= MaxZ )
                            {
                                Bin = { i,j,k };
                            }
                        }
                    }
                }
            }
        }
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
    TRACE("MappingAdaptedMeshToInitalGeometry")
    std::map<unsigned, c_vector<double, SPACE_DIM> > InitalPositionOfRemeshedNodes;
    mInitalPositionOfRemeshedNodes.clear();
    
    SetCentroidMap();
    TRACE("A")
    assert(SPACE_DIM == 3);
    assert(ELEMENT_DIM == 2);

    // Now iterate over the nodes of the new mesh, and find the cloesest element from the old mesh -
    // To accuratly determine the closest, I need to find the two closest, and then determine which one it is closer to. Sometimes might not be super clear

    for (typename MutableMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator iter = mNew_mesh.GetNodeIteratorBegin();
         iter != mNew_mesh.GetNodeIteratorEnd();
         ++iter)
    {
        unsigned node_index = iter->GetIndex();
        c_vector<double, SPACE_DIM> NewNodeLocation = iter->rGetLocation();

        c_vector<double, 3> ClosestElementOrEdge = GetClosestElementInOldMesh(node_index, NewNodeLocation);
        if (ClosestElementOrEdge[0]==0)/* Closest thing is a element */
        {
            double ClosestElement = ClosestElementOrEdge[1];
            InitalPositionOfRemeshedNodes[node_index] = NewNodeInInitalConfigurationFromChangeOfBasis(ClosestElement, NewNodeLocation);
        }else if (ClosestElementOrEdge[0]==1)/* Closest thing is a edge */
        {
            double EdgeNode1 = ClosestElementOrEdge[1];
            double EdgeNode2 = ClosestElementOrEdge[2];
            InitalPositionOfRemeshedNodes[node_index] = NewNodeInInitalConfigurationFromClosestEdge(EdgeNode1, EdgeNode2, NewNodeLocation, node_index);
        }
        else
        {
            /* Should Not be triggered */
            assert(ClosestElementOrEdge[0]==1 || ClosestElementOrEdge[0]==0);

        }
        
    }
    mInitalPositionOfRemeshedNodes = InitalPositionOfRemeshedNodes;
    mOriginalNodePositions.clear();
    mOriginalNodePositions = InitalPositionOfRemeshedNodes;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetClosestElementInOldMesh(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation)
{

    
    assert(SPACE_DIM == 3);
    // This method is super simple. -- Just find the closest element -- it isnt perfect,
    int ClosestElement;
    int ContainingElement;
    bool Accept = 0;
    
    double ClosestElementDistance = 10;
    double ContainingElementDistance =20;

    std::vector<int> Bin = GetBin(NewNodeLocation);
    std::vector<unsigned> ElementsInDaBin= mBin[Bin];

    // ELement is 0 and Edge is 1
    double ElementIdentifier = 0;
    double EdgeIdentifier = 1;
    c_vector<double, 3>  LocalElementOrEdge;
    
    double ContainedInElements =0;
    for (std::vector<unsigned>::iterator elem_index = ElementsInDaBin.begin(); elem_index != ElementsInDaBin.end(); ++elem_index)
    {
        c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*elem_index];
        if (norm_2(NewNodeLocation - Centroid) <= ClosestElementDistance)
        {
            ClosestElementDistance = norm_2(NewNodeLocation - Centroid);
            ClosestElement = *elem_index;
        }
        ContainedInElements +=PointInTriangle3D(NewNodeLocation, *elem_index);
        if (PointInTriangle3D(NewNodeLocation, *elem_index)==1)
        {
            ContainingElementDistance = norm_2(NewNodeLocation - Centroid);
            ContainingElement = *elem_index;
        }
    }
   
    bool ClosetElementIsContained = PointInTriangle2D(NewNodeLocation, ClosestElement) + PointInTriangle3D(NewNodeLocation, ClosestElement);
    
    if( ClosetElementIsContained>0)
    {
        /*
         ---------------------------------------------------------------------
                    The closest element is the containing element 

                 We accept this, and we are very happy. 
         ----------------------------------------------------------------------
        */

        LocalElementOrEdge = Create_c_vector(ElementIdentifier,ClosestElement,0);
        Accept = 1;

        double ClosestCentroidDistance = abs(norm_2(NewNodeLocation - mCentroidMap[ClosestElement]));
        assert(ClosestCentroidDistance < 2);


    }
    else if ( ClosetElementIsContained==0)
    {
        if (ContainedInElements!=0 && ContainingElementDistance < 2*ClosestElementDistance)  // Not in the closest, but is in ContainingElement, and the containing element is closish, just accept i
        { 
            /*
            ---------------------------------------------------------------------------------------------------
                The containing element is pretty close, close than two times the distance of the closes

                            We accept this, I am happy for now, but might revisit. 
            ----------------------------------------------------------------------------------------------------
            */

            LocalElementOrEdge = Create_c_vector(ElementIdentifier,ContainingElement,0);
            Accept = 1;

            double ClosestCentroidDistance = abs(norm_2(NewNodeLocation - mCentroidMap[ClosestElement]));
            assert(ClosestCentroidDistance < 2);
            
        }
        else 
        {
            /*
            -----------------------------------------------------------------
                        The closest thing is probably an edge. 

                The containing element is not anywhere near the new node 

                       Going to need to iterate over all the edges
            ------------------------------------------------------------------
            */
            double EdgeDistance = 10;
            std::pair<unsigned, unsigned> ClosestEdge;

            /* Iterate over the local bin --- Might remove this, dependeing on how things go ... I think things are okay with this */
            std::vector< std::pair<unsigned, unsigned> > EdgesInDaBin= mEdgeBin[Bin];
    
            for (std::vector< std::pair<unsigned, unsigned> >::iterator Edge_iter = EdgesInDaBin.begin(); Edge_iter != EdgesInDaBin.end(); ++Edge_iter)
            {
                std::pair<unsigned, unsigned> edgeIndex = (*Edge_iter);
                unsigned NodeIndexA = (*Edge_iter).first;  unsigned NodeIndexB = (*Edge_iter).second;

                c_vector<double, SPACE_DIM> P1 = this->GetNode(NodeIndexA)->rGetLocation(); c_vector<double, SPACE_DIM> P2 = this->GetNode(NodeIndexB)->rGetLocation();
                // double ProjectionToEdge = ProjectPointToLine(P1, P2,  NewNodeLocation);
                double ProjectionToEdge = ProjectPointToLine(edgeIndex,  NewNodeLocation);
                
                c_vector<double, SPACE_DIM> EdgeVector = (P2 - P1)/norm_2(P2 - P1);

                /* Distance between New node and point 1 on the line */
                c_vector<double, SPACE_DIM> V = NewNodeLocation - P1;
                
                /*  Normal to the plane containing the edge and the new node */
                c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(EdgeVector, V)/norm_2(V); // Not sure why need /norm_2(V)
                PlaneNormal /= norm_2(PlaneNormal);

                /*  Now can find the normal to the edge contaned in the plane defined by the three points */
                c_vector<double, SPACE_DIM> EdgeNormal = VectorProduct(EdgeVector, PlaneNormal);
                EdgeNormal/=norm_2(EdgeNormal);
                
                double DistanceToEdge = abs(inner_prod(EdgeNormal, V));
                PRINT_3_VARIABLES(ProjectionToEdge, DistanceToEdge, ProjectionToEdge- DistanceToEdge )
                assert(ProjectionToEdge == DistanceToEdge );
                if (DistanceToEdge <= EdgeDistance)
                {
                    EdgeDistance = DistanceToEdge;
                    ClosestEdge = *Edge_iter;
                }
            }


            /*  Get the nodes from the edge */
            c_vector<double, SPACE_DIM> P1 = this->GetNode(ClosestEdge.first)->rGetLocation(); c_vector<double, SPACE_DIM> P2 = this->GetNode(ClosestEdge.second)->rGetLocation();
            std::set<unsigned> elements_containing_node1 = this->GetNode(ClosestEdge.first)->rGetContainingElementIndices();   std::set<unsigned> elements_containing_node3 = this->GetNode(ClosestEdge.second)->rGetContainingElementIndices();


            /* Get containing elements */
            std::set<unsigned> shared_elements;
            std::set_intersection(elements_containing_node1.begin(), elements_containing_node1.end(), elements_containing_node3.begin(), elements_containing_node3.end(), std::inserter(shared_elements, shared_elements.begin()));

            
            if (shared_elements.size() == 1)
            {
                 /*
                -----------------------------------------------------------------
                    We have a boundary node -- We need to mark that this is a boundary, 
                    and ensure that the new node is on the boundary 
                ------------------------------------------------------------------
                */

                // LocalElementOrEdge = Create_c_vector(EdgeIdentifier,ClosestEdge);
                LocalElementOrEdge = Create_c_vector(EdgeIdentifier,ClosestEdge.first, ClosestEdge.second);
                Accept=1; 
            }
            else  // In case there are two connected elememtns
            {
                typename std::set<unsigned>::iterator CommonElements = shared_elements.begin();
                unsigned Element1 = *CommonElements;


                
                std::advance(CommonElements, 1);
                unsigned Element2 = *CommonElements;
                if (PointInTriangle3D(NewNodeLocation, Element1) == 1)
                {
                    TRACE("A")
                    ClosestElement = Element1;
                    LocalElementOrEdge = Create_c_vector(ElementIdentifier,ClosestElement,0);
                    Accept=1;
                }
                else if (PointInTriangle3D(NewNodeLocation, Element2) == 1)
                {
                    TRACE("B")
                    ClosestElement = Element2;
                    LocalElementOrEdge = Create_c_vector(ElementIdentifier,ClosestElement,0);
                    Accept=1;
                    
                }else
                {
                    TRACE("I dont want to get here -- this is BAD")
                    TRACE("HERE WE HAVE PROBLEMMMMMMMM!!!!!!!")
                    // assert(1==0);
                    // ClosestElement = WhichElement(P1, P2, NewNodeLocation, Element1, Element2);


                    double DistanceToElement1 = DistanceBetweenPointAndElement(NewNodeLocation, Element1, EdgeDistance);
                    double DistanceToElement2 = DistanceBetweenPointAndElement(NewNodeLocation, Element2, EdgeDistance);
                    if (DistanceToElement1 < DistanceToElement2)
                    {
                        // AcceptElement1 
                        ClosestElement = Element1;

                    }
                    else if (DistanceToElement2 < DistanceToElement1)
                    {
                        // AcceptElement2
                        ClosestElement = Element2;
                        
                    }else if (DistanceToElement2 == DistanceToElement1)
                    {
                        TRACE("didnt expect to get here, accept the edge??? or maybe either??? ")
                        ClosestElement = Element2;
                        LocalElementOrEdge = Create_c_vector(EdgeIdentifier,ClosestEdge.first, ClosestEdge.second);
                        Accept=1; 
                        
                    }



                    // PRINT_VARIABLE(ClosestElement)
                    // LocalElementOrEdge = Create_c_vector(ElementIdentifier,ClosestElement,0);
                    // PRINT_VECTOR(LocalElementOrEdge)
                    Accept=1;

                }
            }
        }
    
        
    }

    assert(Accept ==1);

    

    
    return LocalElementOrEdge;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double  HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::DistanceBetweenPointAndElement( c_vector<double, SPACE_DIM>  NewPoint, unsigned OldElement, double DistanceToNearestLine)
{
    // Orthogal projection onto plane, if the point is in the plane, then we find the distance between this point and the new point, that will give us the distacce
    // If this point is not in the plane, then the closest point is on one of the edges, so just test the edges?? 
    // If these points are the same distance apart then lets say the closes this is the edge and go with that?>>? 

    double ClosestPoint;
    std::pair<double, c_vector<double, SPACE_DIM> >  ProjectionToTheElement = ProjectPointToPlane(NewPoint, OldElement);

    bool IsInElement = PointInTriangle2D(ProjectionToTheElement.second, OldElement);

    if (IsInElement ==0)// The closest point isnt in an element, so the closest point is on an edge :S This will be the closest edge. Have this ordered right 
    {
        ClosestPoint = DistanceToNearestLine;
    }else
    {
        ClosestPoint = ProjectionToTheElement.first;
    }
    

    return ClosestPoint;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<double, c_vector<double, SPACE_DIM> >  HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::ProjectPointToPlane( c_vector<double, SPACE_DIM>  NewPoint, unsigned OldElement)
{

    assert(SPACE_DIM == 3);
    Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->rGetMesh().GetElement(OldElement);
    // Collect the intial configuration of the nodes, and the deformed config, send them into the mapping function
    c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> x;
    for (int i = 0; i < 3; ++i)
    {
        x[i] = this->rGetMesh().GetNode(p_element->GetNodeGlobalIndex(i))->rGetLocation();
    }

    // u1 and u2 are two edge vectors of the element
    c_vector<double, SPACE_DIM> u1 = x[1]-x[0];
    c_vector<double, SPACE_DIM> u2 = x[2]-x[0];

    c_vector<double, SPACE_DIM> NormalToPlane = VectorProduct(u1, u2)/norm_2(VectorProduct(u1, u2));
    // a is the vector between the first point of the element and the new point
    c_vector<double, SPACE_DIM> a = NewPoint - x[0];

    double Distance = inner_prod(NormalToPlane, a);

    c_vector<double, SPACE_DIM> NearestPointInThePlane = NewPoint -abs(Distance) * NormalToPlane;
    std::pair<double, c_vector<double, SPACE_DIM> > ProjectionToTheElement = std::pair<double , c_vector<double, SPACE_DIM> >(Distance, NearestPointInThePlane);
  
    return ProjectionToTheElement;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::ProjectPointToLine(std::pair<unsigned, unsigned> edgeIndex , c_vector<double, SPACE_DIM>  NewPoint)
{

    unsigned NodeIndexA = edgeIndex.first;  unsigned NodeIndexB = edgeIndex.second;

    c_vector<double, SPACE_DIM> x1 = this->GetNode(NodeIndexA)->rGetLocation(); c_vector<double, SPACE_DIM> x2 = this->GetNode(NodeIndexB)->rGetLocation();

    c_vector<double, SPACE_DIM> a = NewPoint - x1;

    c_vector<double, SPACE_DIM> unitEdgeVector = x2-x1;
    unitEdgeVector/=norm_2(unitEdgeVector);

    // double Projection1 = inner_prod(mCentroidMap[Element1] - EdgeMidpoint, PlaneNormal);
    double d1 = inner_prod(a, unitEdgeVector);
    c_vector<double, SPACE_DIM> NearestPointOnLine = d1 * unitEdgeVector + x1;
    double DistanceToLine = norm_2(NearestPointOnLine - NewPoint);
    return DistanceToLine;

}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::ProjectPointToLine(c_vector<double, SPACE_DIM> x1, c_vector<double, SPACE_DIM> x2, c_vector<double, SPACE_DIM>  NewPoint)
{

    c_vector<double, SPACE_DIM> a = NewPoint - x1;

    c_vector<double, SPACE_DIM> unitEdgeVector = x2-x1;
    unitEdgeVector/=norm_2(unitEdgeVector);

    // double Projection1 = inner_prod(mCentroidMap[Element1] - EdgeMidpoint, PlaneNormal);
    double d1 = inner_prod(a, unitEdgeVector);
    c_vector<double, SPACE_DIM> NearestPointOnLine = d1 * unitEdgeVector + x1;
    double DistanceToLine = norm_2(NearestPointOnLine - NewPoint);
    return DistanceToLine;

}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::WhichElement(c_vector<double, SPACE_DIM> P1, c_vector<double, SPACE_DIM> P2, c_vector<double, SPACE_DIM> NewNodeLocation, unsigned Element1, unsigned Element2)
{
    // PRINT_  VARIABLE
    PRINT_2_VARIABLES(Element1,Element2)


    assert(SPACE_DIM == 3);
    unsigned ChosenElement;
    //    FIgure out where it is the same way I used as above in method 2 -- the current way I have doesnt work!!!!!!
    c_vector<double, 3> Edge = P2 - P1;
    c_vector<double, 3> EdgeMidpoint = (P2 + P1) / 2;
    PRINT_VECTOR(Edge)

    PRINT_VECTOR(EdgeMidpoint)

    c_vector<double, 3> Normal1 = GetElementNormal(Element1);
    c_vector<double, 3> Normal2 = GetElementNormal(Element2);
    
PRINT_VECTOR(Normal1)
PRINT_VECTOR(Normal2)
    // Average normal of the two elements
    c_vector<double, 3> Normal = (Normal1 + Normal2) / 2;

PRINT_VECTOR(Normal)
    // The normal should be pointing outwards, so the direction of the plane normal will depend on the direction of the edge vector.
    // The edge vector is depenedent on the order of the nodes pair, which should be sorted numerically
    c_vector<double, 3> PlaneNormal = VectorProduct(Normal, Edge);
    PlaneNormal /= norm_2(PlaneNormal);

PRINT_VECTOR(PlaneNormal)

    // Need to determine the direction of the normal

    c_vector<double, SPACE_DIM> NormalArm = PlaneNormal;
    PRINT_VECTOR(NormalArm)

    double Projection1 = inner_prod(mCentroidMap[Element1] - EdgeMidpoint, PlaneNormal);
    double Projection2 = inner_prod(mCentroidMap[Element2] - EdgeMidpoint, PlaneNormal);
    if (Projection1 == Projection2)
    {
        PRINT_2_VARIABLES(Projection1, Projection2)
    }
    PRINT_2_VARIABLES(Projection1, Projection2)
    assert(Projection1 != Projection2);
    // if on the same side as element 2, then projection 2 will be the shorter
    if (Projection1 > Projection2)
    {
        double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
        // If side is positive, then the new node will be on the same side as the plane normal is pointing, which is element 2. 
        if (side > 0)
        {
            // TRACE("The new node is on the element ONE side of the midline plane ")
            ChosenElement = Element2;
        }else
        {
            // TRACE("The new node is on the element ONE side of the midline plane ")
            ChosenElement = Element1;
        }

    }else if (Projection1 < Projection2)
     if (Projection1 > Projection2)
    {
        double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
        // If side is positive, then the new node will be on the same side as the plane normal is pointing, which is element 2. 
        if (side > 0)
        {
            // TRACE("The new node is on the element ONE side of the midline plane ")
            ChosenElement = Element1;
        }else
        {
            // TRACE("The new node is on the element ONE side of the midline plane ")
            ChosenElement = Element2;
        }

    }
    PRINT_VARIABLE(ChosenElement)
    
    return ChosenElement;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::NewNodeInInitalConfigurationFromClosestEdge(unsigned EdgeNode1, unsigned EdgeNode2, c_vector<double, SPACE_DIM> NewNode, unsigned NodeIndex)
{
        assert(SPACE_DIM == 3);
        // PRINT_VARIABLE(ClosestEdge_OldMeshIndex)
        // /*  Get the nodes from the edge */
        TRACE("Jess is here")
        c_vector<double, SPACE_DIM> P1 = this->GetNode(EdgeNode1)->rGetLocation(); c_vector<double, SPACE_DIM> P2 = this->GetNode(EdgeNode2)->rGetLocation();
            // // I want the distance from the edge

        c_vector<double, SPACE_DIM> EdgeVector = (P2 - P1);

            /* Distance between New node and point 1 on the line */
            c_vector<double, SPACE_DIM> V = NewNode - P1;

            /*  Normal to the plane containing the edge and the new node */
            c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(EdgeVector, V)/norm_2(V); // Not sure why need /norm_2(V)
            PlaneNormal /= norm_2(PlaneNormal);

            /*  Now can find the normal to the edge contaned in the plane defined by the three points */
            c_vector<double, SPACE_DIM> EdgeNormal = VectorProduct(EdgeVector, PlaneNormal);
            EdgeNormal/=norm_2(EdgeNormal);
    
            double DistanceToEdge = abs(inner_prod(EdgeNormal, V));   
            c_vector<double, SPACE_DIM> ProjectToEdge= DistanceToEdge * EdgeNormal; // DOnt know if this will be pushing in the right direction 


            // PRINT_VECTOR(ProjectToEdge)
            c_vector<double, 3> ProjectedPosition = NewNode + ProjectToEdge; 

            // Chekc this new point is closer to line
            double NewDistance = GetDistanceToLine(ProjectedPosition, P1, P2);


            if (NewDistance > DistanceToEdge)
            {
                TRACE("Projected wrong way")
                ProjectedPosition = NewNode - ProjectToEdge; 
                NewDistance = GetDistanceToLine(ProjectedPosition, P1, P2);
                if (NewDistance > DistanceToEdge)
                {
                    TRACE("Not going to work")
                }
            }

            // DIstance from node 1 and node 2
            double DistanceFromNode1 = norm_2(P1 - ProjectedPosition);
            // TRACE("D")
            double DistanceFromNode2 = norm_2(P2 - ProjectedPosition);
            double EdgeLength = norm_2(P2 - P1);

            double ProportionalDistanceFromNode1 = DistanceFromNode1/EdgeLength;
            // PRINT_4_VARIABLES(DistanceFromNode1,DistanceFromNode2, EdgeLength, ProportionalDistanceFromNode1)

            // // Get orginalNode Positions
            // PRINT_2_VARIABLES(EdgeNode1,EdgeNode2)
            c_vector<double, 3> OrginialPositonOfNode1 = mOriginalNodePositions[EdgeNode1];
            c_vector<double, 3> OrginialPositonOfNode2 = mOriginalNodePositions[EdgeNode2];
     
     
            c_vector<double, 3> OriginalEdgeVector = OrginialPositonOfNode2- OrginialPositonOfNode1;
            double originalLength = norm_2(OriginalEdgeVector);

            c_vector<double, 3> OriginalPosition = OrginialPositonOfNode1 + ProportionalDistanceFromNode1*OriginalEdgeVector;///originalLength;



            // Now need point 
            // PRINT_VECTOR(OriginalPosition)
            assert(norm_2(OriginalPosition) <10);
            return OriginalPosition;

}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::NewNodeInInitalConfigurationFromChangeOfBasis(unsigned ClosestElement_OldMeshIndex, c_vector<double, SPACE_DIM> NewNode)
{
    
    assert(SPACE_DIM == 3);
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
   
    assert(norm_2(P_0) <10);
    return P_0;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetDistanceToLine( c_vector<double, SPACE_DIM> NewNode, c_vector<double, SPACE_DIM> P1, c_vector<double, SPACE_DIM>  P2)
{

        c_vector<double, SPACE_DIM> EdgeVector = (P2 - P1);

        /* Distance between New node and point 1 on the line */
        c_vector<double, SPACE_DIM> V = NewNode - P1;

        /*  Normal to the plane containing the edge and the new node */
        c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(EdgeVector, V)/norm_2(V); // Not sure why need /norm_2(V)
        PlaneNormal /= norm_2(PlaneNormal);

        /*  Now can find the normal to the edge contaned in the plane defined by the three points */
        c_vector<double, SPACE_DIM> EdgeNormal = VectorProduct(EdgeVector, PlaneNormal);
        EdgeNormal/=norm_2(EdgeNormal);

        double DistanceToEdge = abs(inner_prod(EdgeNormal, V));
     
  
        return DistanceToEdge;


}

// For the shear and area forces
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetupMembraneConfiguration()
{
    MathsFunctions M;
    assert(SPACE_DIM == 3);
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
        Node<SPACE_DIM>* pNode = this->rGetMesh().GetNode(node_index);
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
    assert(SPACE_DIM == 3);
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
    double Threshold1a = MeanAngle - 1.5 * stdev;
    double Threshold2a = MeanAngle + 1.5 * stdev;

    double Threshold1b = MeanAngle - 2.5 * stdev;
    double Threshold2b = MeanAngle + 2.5 * stdev;
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

            if (Angle < Threshold1b || Angle > Threshold2b)
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
                p_cell1->GetCellData()->SetItem("Curvature", 2);
                p_cell2->GetCellData()->SetItem("Curvature", 2);
                p_cell3->GetCellData()->SetItem("Curvature", 2);
                p_cell4->GetCellData()->SetItem("Curvature", 2);

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
            }else if(Angle < Threshold1a || Angle > Threshold2a)
            {

                // TRACE("Medium curveature")
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


            }
        }
    }





    
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


// Set boundaries
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::MarkBoundaryNodes()
{

    assert(ELEMENT_DIM == 2 & SPACE_DIM == 3);
    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin(); cell_iter != this->mCells.end(); ++cell_iter)
    {
        (*cell_iter)->GetCellData()->SetItem("Boundary", 0);
        (*cell_iter)->GetCellData()->SetItem("FixedBoundary", 0);
    }

    // IF we are checking for boundaries
    if (mSetBoundaries == 1)
    {
        for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin(); cell_iter != this->mCells.end(); ++cell_iter)
        {
            unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);
            Node<SPACE_DIM>* p_node = this->GetNode(node_index);
            // (*cell_iter)->GetCellData()->SetItem("Boundary", 0);
            std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();

            if (containing_elements.size() < 5)
            {
                (*cell_iter)->GetCellData()->SetItem("Boundary", 1);
                bool SetSurroundingNodesAsEdges =0; // TODO make this a member variable I can adapt with a function 
                if (SetSurroundingNodesAsEdges)
                {
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
        }
    }

    // TRACE("In the set up solve for the boundaries")

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

    SetRelativePath(ChasteOutputDirectory, startime);

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
    SetRelativePath(ChasteOutputDirectory);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetStartTime(double StartTime)
{
    mStartTime = StartTime;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetStartTime()
{
    return mStartTime;
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
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::EdgeLengthVariable(double EdgeLengthMultiple)
{
    mEdgeLengthMultiple = EdgeLengthMultiple;
    mVariableEdgeLength = 1;
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
    if (RemeshingSoftwear != "CGAL" && RemeshingSoftwear != "VMTK" && RemeshingSoftwear != "PreAllocatedMatlabMesh")
    {
        EXCEPTION("RemeshingSoftwear must be CGAL or VMTK or PreAllocatedMatlabMesh");
    }
    mRemeshingSoftwear = RemeshingSoftwear;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetOperatingSystem(std::string OperatingSystem)
{
    if (OperatingSystem == "server" || OperatingSystem == "Server" || OperatingSystem == "linux" || OperatingSystem == "Linux")
    {
        mServer = 1;
    }
    else if (OperatingSystem == "mac" || OperatingSystem == "Mac" || OperatingSystem == "local" || OperatingSystem == "Local")
    {
        mServer = 0;
    }
    else 
    {
        EXCEPTION("OperatingSystem must be Mac or Linux");
    }
    
}


    void SetOperatingSystem( std::string OperatingSystem);
    bool mServer =1;



// Set the pathway to the mesh I need to read in
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::PreAllocatedRemeshedMesh(std::string RemeshedMesh)
{
    mPreAllocatedRemeshedMesh = RemeshedMesh;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::CalculateElementNormals(std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edge,
                                                                                        std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> >& nonUnitNormals,
                                                                                        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>& otherNodes)
{
    assert(SPACE_DIM == 3 && ELEMENT_DIM == 2);
    Node<SPACE_DIM>* pNode1 = edge.first;
    Node<SPACE_DIM>* pNode3 = edge.second;

    /*
     *  Find common triangles
    */
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
        case 1: // We found a boundary edge and we finish here
            return true;
        case 2:
            break;
        default:
            PRINT_VARIABLE(shared_elements.size())
            NEVER_REACHED;
    }

    std::set<unsigned>::iterator set_iter = shared_elements.begin();
    Element<ELEMENT_DIM, SPACE_DIM>* pElement1 = this->rGetMesh().GetElement(*set_iter);
    ++set_iter;
    Element<ELEMENT_DIM, SPACE_DIM>* pElement2 = this->rGetMesh().GetElement(*set_iter);

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
    assert(SPACE_DIM == 3);
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

    return Answer;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::PointInTriangle2D(c_vector<double, SPACE_DIM> Point, unsigned ClosestElement)
{
    bool Answer = 0;
    Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->rGetMesh().GetElement(ClosestElement);
    // Collect the intial configuration of the nodes, and the deformed config, send them into the mapping function
    c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> Triangle;
    for (int i = 0; i < 3; ++i)
    {
        Triangle[i] = this->rGetMesh().GetNode(p_element->GetNodeGlobalIndex(i))->rGetLocation();
    }

    if (SameSideOfPlane(Point, Triangle[0], Triangle[1], Triangle[2]) & SameSideOfPlane(Point, Triangle[1], Triangle[0], Triangle[2]) & SameSideOfPlane(Point, Triangle[2], Triangle[0], Triangle[1]))
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
    assert(SPACE_DIM == 3);
    Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->rGetMesh().GetElement(ClosestElement);
    // Collect the intial configuration of the nodes, and the deformed config, send them into the mapping function
    c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> Triangle;
    for (int i = 0; i < 3; ++i)
    {
        Triangle[i] = this->rGetMesh().GetNode(p_element->GetNodeGlobalIndex(i))->rGetLocation();
    }

    c_vector<double, SPACE_DIM> Line1 = Triangle[1] - Triangle[0]; // Vector 1 to 2
    c_vector<double, SPACE_DIM> Line2 = Triangle[2] - Triangle[0]; // Vector 1 to 3
    c_vector<double, SPACE_DIM> Line3 = Triangle[2] - Triangle[1]; // Vector 2 to 3

    c_vector<double, SPACE_DIM> ElementNormal = GetElementNormal(ClosestElement);

    // Need to create a plane containing the triangle.
    // Let  Triangle[0] be the origin
    c_vector<double, SPACE_DIM> v = Point - Triangle[0];
    double dist = inner_prod(ElementNormal, v);
    c_vector<double, SPACE_DIM> projected_point = Point - dist * ElementNormal;
    // Check that my normal was facing the right way -
    assert(inner_prod(ElementNormal, projected_point - Triangle[0]) < 1e-13);

    if (PointInTriangle3D(projected_point, ClosestElement) == 1)
    {
        return dist;
    }
    else
    {
        // Get the closest point on each line ---

        c_vector<double, SPACE_DIM> Normal1 = VectorProduct(Line1, ElementNormal);
        Normal1 /= norm_2(Normal1);

        c_vector<double, SPACE_DIM> Normal2 = VectorProduct(Line2, ElementNormal);
        Normal2 /= norm_2(Normal2);

        c_vector<double, SPACE_DIM> Normal3 = VectorProduct(Line3, ElementNormal);
        Normal3 /= norm_2(Normal3);

        double dist1 = abs(inner_prod(Normal1, (projected_point - Triangle[1])));
        double dist2 = abs(inner_prod(Normal2, (projected_point - Triangle[2])));
        double dist3 = abs(inner_prod(Normal3, (projected_point - Triangle[2])));
        // Find the smallest distance

        double MinDistane = std::min(dist1, std::min(dist2, dist3));

        double PointDistanceToElement = sqrt(MinDistane * MinDistane + dist * dist);
        return PointDistanceToElement;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::PointInTriangle3D(c_vector<double, SPACE_DIM> Point, c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> Triangle)
{
    assert(SPACE_DIM == 3);
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
    if (SameSideOfPlane(Point, Triangle[0], Triangle[1], Triangle[2]) & SameSideOfPlane(Point, Triangle[1], Triangle[0], Triangle[2]) & SameSideOfPlane(Point, Triangle[2], Triangle[0], Triangle[1]))
    {
        Answer1 = 1;
    }
    else
    {
        Answer1 = 0;
    }
    if (Answer1 != 1)
    {
        PRINT_2_VARIABLES(Answer, Answer1)
    }
    return Answer;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SameSideOfPlane(c_vector<double, SPACE_DIM> P1, c_vector<double, SPACE_DIM> P2, c_vector<double, SPACE_DIM> a, c_vector<double, SPACE_DIM> b)
{
    bool Answer;
    assert(SPACE_DIM == 3);
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
    assert(SPACE_DIM == 3);
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
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetBinLowerX(std::vector<double> Bin)
{
    typename std::vector<double>::iterator iter = Bin.begin();
    return *iter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetBinUpperX(std::vector<double> Bin)
{
    typename std::vector<double>::iterator iter = Bin.begin();
    std::advance(iter, 1);
    return *iter;
}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetBinLowerY(std::vector<double> Bin)
{
    typename std::vector<double>::iterator iter = Bin.begin();
    std::advance(iter, 2);
    return *iter;
}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetBinUpperY(std::vector<double> Bin)
{
    typename std::vector<double>::iterator iter = Bin.begin();
    std::advance(iter, 3);
    return *iter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetBinLowerZ(std::vector<double> Bin)
{
    typename std::vector<double>::iterator iter = Bin.begin();
    std::advance(iter, 4);
    return *iter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetBinUpperZ(std::vector<double> Bin)
{
    typename std::vector<double>::iterator iter = Bin.begin();
    std::advance(iter, 5);
    double value = *iter;
    return *iter;
}


// this needs fixing too 

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::CalculateAspectRatio(c_vector<double, SPACE_DIM> Node1, c_vector<double, SPACE_DIM> Node2,c_vector<double, SPACE_DIM> Node3 )
{
    assert(SPACE_DIM ==3);
    c_vector<double, SPACE_DIM> vector_12 = Node2 - Node1; 
    c_vector<double, SPACE_DIM> vector_13 = Node3 - Node1; 
    c_vector<double, SPACE_DIM> vector_23 = Node3 - Node2; 

    double l1 = norm_2(vector_12); 
    double l2 = norm_2(vector_13); 
    double l3 = norm_2(vector_23);

    double ElementArea = 0.5*norm_2(VectorProduct(vector_12, vector_13));
    double MaxEdgeLength = std::max( std::max(l1,l3) , std::max(l2,l3)) ; 

    double AspectRatio = 4/sqrt(3) * ElementArea/(MaxEdgeLength * MaxEdgeLength);



    // The second option 
    // double MinEdgeLength = std::min( std::min(l1,l3) , std::min(l2,l3)) ; 
    // double AR = MinEdgeLength/MaxEdgeLength;

    // Thrid option
    // c_vector<double, SPACE_DIM> MaxEdge;
    // c_vector<double, SPACE_DIM>P1; 
    // c_vector<double, SPACE_DIM> P2 ;
    // c_vector<double, SPACE_DIM> P0;
    // if (l1 == MaxEdgeLength)
    // {
    //     MaxEdge = vector_12;
    //     P1 = Node1; 
    //     P2 = Node2;
    //     P0 = Node3;

    // }else if(l2 == MaxEdgeLength)
    // {
    //     MaxEdge = vector_13;
    //     P1 = Node1;
    //     P2 = Node3;
    //     P0 = Node2;
    // }else if (l3 == MaxEdgeLength)
    // {
    //     MaxEdge = vector_23;
    //     P1 = Node2;
    //     P2 = Node3;
    //     P0 = Node1;
    // }else
    // {
    //     assert(0==1);
    // }
    // double A = P1[1]-P0[1];
    // double B = P1[0]-P0[0];
    // double Height = abs((P2[0]-P1[0])*A - B*(P2[1]-P1[1]) )/MaxEdgeLength;
    // double ThirdAspectRatio = 2*ElementArea/MaxEdgeLength;

    // PRINT_2_VARIABLES(ElementArea, MaxEdgeLength)
    // PRINT_3_VARIABLES(AspectRatio, AR,ThirdAspectRatio)
    return AspectRatio;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::MinimumElementAspectRatio()
{

     for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        CellPtr p_cell = this->GetCellUsingLocationIndex(node_iter->GetIndex());
        p_cell->GetCellData()->SetItem("AspectRatio", 1.1);
    }

   
    double MinimumAspectRatio = 100;
    assert(SPACE_DIM == 3);
    assert(ELEMENT_DIM == 2);
    std::vector<double> AspectRatioVector;
    // Loop over the old map and get the centroids of the old map
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
         elem_iter != this->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {

        // I want to exclude the edge region 
        unsigned elem_index = elem_iter->GetIndex();
        Node<SPACE_DIM>* pNode0 = this->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0));
        Node<SPACE_DIM>* pNode1 = this->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(1));
        Node<SPACE_DIM>* pNode2 = this->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));


        CellPtr p_cell0 = this->GetCellUsingLocationIndex(pNode0->GetIndex());
        CellPtr p_cell1 = this->GetCellUsingLocationIndex(pNode1->GetIndex());
        CellPtr p_cell2 = this->GetCellUsingLocationIndex(pNode2->GetIndex());       
        
        
        if( p_cell0->GetCellData()->GetItem("Boundary") == 0 && p_cell1->GetCellData()->GetItem("Boundary") == 0 && p_cell2->GetCellData()->GetItem("Boundary") == 0 )
        {
            // I am also going to record the element aspect ratios at this point, because might as well while I am iterating over the elements, multitask
            double AspectRatio = CalculateAspectRatio(pNode0->rGetLocation(), pNode1->rGetLocation(), pNode2->rGetLocation() );

            AspectRatioVector.push_back(AspectRatio);

            double AR2 = p_cell2->GetCellData()->GetItem("AspectRatio");
            double AR0 = p_cell0->GetCellData()->GetItem("AspectRatio");
            double AR1 = p_cell1->GetCellData()->GetItem("AspectRatio");
            if (AspectRatio < AR1)
            {
                p_cell1->GetCellData()->SetItem("AspectRatio", AspectRatio);
            }
            if (AspectRatio < AR0)
            {
                p_cell0->GetCellData()->SetItem("AspectRatio", AspectRatio);
            }
            if (AspectRatio < AR2)
            {
                p_cell2->GetCellData()->SetItem("AspectRatio", AspectRatio);
            }
        }
    }
    return AspectRatioVector;
    // std::vector<double> quartilesNeeded = { 0.25, 0.5, 0.75};
    // std::vector<double> quartiles = Quantile(AspectRatioVector, quartilesNeeded);
    // PRINT_VECTOR(quartiles)
    // double Minimum = *std::min_element( std::begin(AspectRatioVector), std::end(AspectRatioVector) );
    // // typename std::vector<double>::iterator Iterator = quartiles.begin();
    // // double Quartile1 = *Iterator;
    // // std::advance(Iterator, 1);
    // // double Med = *Iterator;
    // // std::advance(Iterator, 1);
    // // double Quartile3 = *Iterator;

    // // c_vector<double, SPACE_DIM> AspectRatioQuatiles = Create_c_vector(Quartile1, Med, Quartile3);

    // // PRINT_VECTOR(AspectRatioQuatiles)
    // return Minimum;
}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::Quantile(std::vector<double>& data, std::vector<double>& probs)
{
    if (data.empty())
    {
        return std::vector<double>();
    }

    if (1 == data.size())
    {
        return std::vector<double>(1, data[0]);
    }

    std::sort(data.begin(), data.end());
    std::vector<double> quantiles;

    for (size_t i = 0; i < probs.size(); ++i)
    {
        // double poi = Lerp<double>(-0.5, data.size() - 0.5, probs[i]);
        double poi = (1 - probs[i]) * -0.5 + probs[i] * (data.size() - 0.5);

        size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
        size_t right = std::min(int64_t(std::ceil(poi)), int64_t(data.size() - 1));

        double datLeft = data.at(left);
        double datRight = data.at(right);

        // double quantile = Lerp<double>(datLeft, datRight, poi - left);
        double quantile =  (1 - (poi - left)) * datLeft + (poi - left) * datRight;

        quantiles.push_back(quantile);
    }
    return quantiles;
}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::UpdateBoundaryConditions()
{
    assert(ELEMENT_DIM == 2 && SPACE_DIM == 3);
   mUpdateComplete =1;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetUpdateBoundaryConditions()
{
    assert(ELEMENT_DIM == 2 && SPACE_DIM == 3);
   return mUpdateComplete;
}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationParameters(rParamsFile);
}


// Explicit instantiation
template class HistoryDepMeshBasedCellPopulation<1, 2>;
template class HistoryDepMeshBasedCellPopulation<1, 3>;
template class HistoryDepMeshBasedCellPopulation<2, 2>;
template class HistoryDepMeshBasedCellPopulation<2, 3>;
template class HistoryDepMeshBasedCellPopulation<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS2(HistoryDepMeshBasedCellPopulation, 1, 3)
EXPORT_TEMPLATE_CLASS2(HistoryDepMeshBasedCellPopulation, 2, 3)
EXPORT_TEMPLATE_CLASS2(HistoryDepMeshBasedCellPopulation, 2, 2)
EXPORT_TEMPLATE_CLASS2(HistoryDepMeshBasedCellPopulation, 3, 3)
/*

*/

#ifndef HISTORYDEPMESHBASEDCELLPOPULATION_HPP_
#define HISTORYDEPMESHBASEDCELLPOPULATION_HPP_

#include <algorithm>
#include <map>
#include "AbstractCentreBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MutableMesh.hpp"
// #include "VertexMesh.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"

#include "CellBasedEventHandler.hpp"
#include "CellId.hpp"

#include <time.h>
#include "Debug.hpp"
#include "HistoryDepMutableMesh.hpp"
#include "NodeVelocityWriter.hpp"
#include "NodesOnlyMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "VoronoiDataWriter.hpp"
#include "VtkMeshReader.hpp"
#include "VtkMeshWriter.hpp"
// #include <filesystem>

#include <boost/filesystem.hpp>

#include <ctime>
#include "MathsFunctions.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"


#include <math.h>
#include "Debug.hpp"



/**
 * A facade class encapsulating a mesh-based 'cell population'.
 *
 * Contains a group of cells and maintains the associations between cells and
 * nodes in the mesh.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class HistoryDepMeshBasedCellPopulation : public MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>
{
    // friend class TestMeshBasedCellPopulation;
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM> >(*this);
        // archive & boost::serialization::base_object<AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM> >(*this);



        archive& mOriginalNodePositions;
        archive& mInitalPositionOfRemeshedNodes;
        archive& mNew_mesh;
        archive& mInitalVectors;
        archive& mACoefficients;
        archive& mBCoefficients;
        archive& mArea0;
        archive& mOriginalAngles;
        archive& mTargetRemeshingEdgeLength;
        archive& mIterations;
        archive& mRelativePath;
        archive& mPrintRemeshedIC;
        // archive& mMaxEdgelength;
        // archive& mMapOfProbNodes;
        archive& mNumberOfChanges;//////
        archive& mRemeshingSoftwear;
        archive& mNearestNodesMap;
        archive& mNx;
        archive& mNy;
        archive& mNz;
        archive& mCentroidMap;
        archive& mStartTime;
        archive& mServer;
        archive& mSetUpInitialConfigurations;
        //
        archive& mChasteOutputDirectory;
        archive& mSetBoundaries;
        archive& mPreAllocatedRemeshedMesh;
        archive& mCounter;
        archive& mNewNodeToOldElementMap;
        archive& mNewNodeToOldElementDistanceMap;
        archive& mMappingVariables_a_b;
        archive& mMappingVariables_alpha;
        archive& mMappingVariables_z_basis;
        archive& mMappingVariables_PointInNewRef;
        archive& mMappingVariables_Difference;
        archive& mMappingVariables_P_Translated;
        archive& mMappingVariables_Cs;
        archive& mBinMap;
        archive& mBin;
        archive& mEdgeBin;
        archive& mBinCoords;
        archive& mDIM;
        archive& mMaxX;
        archive& mMinX;
        archive& mMaxY;
        archive& mMinY;
        archive& mMaxZ;
        archive& mMinZ;
        archive& mVariableEdgeLength;
        archive& mEdgeLengthMultiple;
        archive& mUpdateComplete;

        // archive& mOutputMeshInVtk;
        // archive& mWriteVtkAsPoints;
        // this->Validate();
    }





protected:
    /** Static cast of the mesh from AbstractCellPopulation */
    HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM> mNew_mesh;

public:
        // bool mOutputMeshInVtk =1;
        // bool mWriteVtkAsPoints = 0;

    /**
     * Create a new cell population facade from a mesh and collection of cells.
     *
     * There must be precisely 1 cell for each node of the mesh.
     *
     * @param rMesh a mutable tetrahedral mesh
     * @param rCells cells corresponding to the nodes of the mesh
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     * @param validate whether to validate the cell population
     */
    HistoryDepMeshBasedCellPopulation(MutableMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices = std::vector<unsigned>(),
                                      bool deleteMesh = false,
                                      bool validate = true);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a mutable tetrahedral mesh.
     */
    HistoryDepMeshBasedCellPopulation(MutableMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

    /**
     * Destructor.
     */
    virtual ~HistoryDepMeshBasedCellPopulation();

    
    // Set bool to turn off and on the initial config setup -- sometimes i will need to skip this
    void SetUpInitialConfig(bool SetUpInitialConfigurations);
    bool mSetUpInitialConfigurations =1;




    /*
     * Need a direction to the folder everything is savinging to  
     */

    
    void SetChasteOutputDirectory(std::string ChasteOutputDirectory, double startime);
    void SetChasteOutputDirectory(std::string ChasteOutputDirectory);

    std::string mChasteOutputDirectory;

    void SetRelativePath(std::string ChasteOutputDirectory, double startime);
    void SetRelativePath(std::string ChasteOutputDirectory);
    void SetPaths(std::string ChasteOutputDirectory);
    std::string mRelativePath;


    void SetStartTime(double StartTime);
    double GetStartTime();
    double mStartTime=0;

    void ExecuteHistoryDependentRemeshing();
    void MappingAdaptedMeshToInitalGeometry();

    void SetBoundaries(bool SetB);
    bool mSetBoundaries = 1;

    // Function to find the cloeset element in last mesh
    
    // c_vector<double, 3> GetClosestElementInOldMesh(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation);
     double GetClosestElementInOldMesh(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation);



    // std::map<unsigned, unsigned> mMapOfProbNodes;

    // When I get a population of new cells, I need to be able to give these new cells the same CellData as the Precious cells, this is particullary important for the
    // void UpdateCellData();

    /**
     * Map the new node to the old element   Element_0 has the node position from the inital condition 
     * Element_0 has the node position from the inital condition 
     * DeformedElement has the node position from the deformed position 
     * NewNodes is the position of the new node in the deformed configuration 
     * 
     * Goal is to get the inital position for this remeshed node
     */

    
    // c_vector<double, 3> NewNodeInInitalConfigurationFromClosestEdge(unsigned EdgeNode1, unsigned EdgeNode2, c_vector<double, SPACE_DIM> NewNode,unsigned NodeIndex);
    c_vector<double, 3> NewNodeInInitalConfigurationFromChangeOfBasis(unsigned ClosestElement_OldMeshIndex, c_vector<double, SPACE_DIM> NewNode, unsigned NodeIndex);
    std::map<unsigned, c_vector<double, SPACE_DIM> > mInitalPositionOfRemeshedNodes;
    int mNumberOfChanges = 1;

    /**
     * Writing out the inital node locations in to their own vtu so I can see what is going on 
     */

    void WriteOutMappedInitalConfig();

    /**
     * Calling the python code to do the remeshing  
     */
    void RemeshGeometry();
    void RemeshGeometryWithVMTK();
    void TakeInPreAllocatedRemeshGeometry();
    void PreAllocatedRemeshedMesh(std::string RemeshedMesh);
    std::string mPreAllocatedRemeshedMesh ="/data/vascrem/RemeshedGeometry" ;


    double mCounter = 0;

    void CheckCurvature();
    void SetRemeshingSoftwear(std::string RemeshingSoftwear);
    std::string mRemeshingSoftwear = "CGAL";
    void SetOperatingSystem( std::string OperatingSystem);
    bool mServer =1;

    void SaveInitalConditions();
    std::map<unsigned, c_vector<double, SPACE_DIM> > mOriginalNodePositions;
    std::map<unsigned, c_vector<double, SPACE_DIM> > GetInitalNodePositions();
    // void SetMaxEdgelength();
    bool PointInTriangle3D(c_vector<double, SPACE_DIM> Point, unsigned ClosestElement);
    bool PointInTriangle3D(c_vector<double, SPACE_DIM> Point, c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> Triangle);
    bool PointInTriangle2D(c_vector<double, SPACE_DIM> Point, unsigned ClosestElement);
    double ClosestPointInTriangle(c_vector<double, SPACE_DIM> Point, unsigned ClosestElement);
    unsigned WhichElement(c_vector<double, SPACE_DIM> P1, c_vector<double, SPACE_DIM> P2, c_vector<double, SPACE_DIM> NewNodeLocation, unsigned Element1, unsigned Element2);

    bool SameSideOfPlane(c_vector<double, SPACE_DIM> P1, c_vector<double, SPACE_DIM> P2, c_vector<double, SPACE_DIM> a, c_vector<double, SPACE_DIM> b);
    bool SameSide(c_vector<double, SPACE_DIM> P1, c_vector<double, SPACE_DIM> P2, c_vector<double, SPACE_DIM> a, c_vector<double, SPACE_DIM> b);

    // double mMaxEdgelength;

    void SetMinArea();

    /**
     * THis method will loop over the edges in the remeshed mesh and make sure the edge lengths in the IC are not longer than the edge lenght in the original remesh's original config  
    */

    void CheckRemeshedIC();

    /**
     * Set up the intial conditions for for the intial configuration for the forces 
     */

    void SetupMembraneConfiguration();
    std::map<unsigned, c_vector<c_vector<double, 2>, 3> > mInitalVectors;


    // Mapping the new nodes to the old elements; am recording this for debugging purposes
    std::map<unsigned,double > mNewNodeToOldElementMap;

    // Bunch of members and methods to store and later save data for debugging

    std::map<unsigned,double >  GetNodeToOldElementMap();
    std::map<unsigned,double > GetNewNodeToOldElementDistanceMap();

    std::map<unsigned,double >  mNewNodeToOldElementDistanceMap;


    std::map<unsigned,c_vector<double,2> > GetMappingVariables_a_b( );
    std::map<unsigned,c_vector<double,2> > mMappingVariables_a_b;

    std::map<unsigned, double > GetMappingVariables_alpha( );
    std::map<unsigned, double > mMappingVariables_alpha;

    std::map<unsigned,c_vector<double,3> > GetMappingVariables_z_basis( );
    std::map<unsigned,c_vector<double,3> > mMappingVariables_z_basis;


    std::map<unsigned,c_vector<double,3> > GetMappingVariables_PointInNewRef( );
    std::map<unsigned,c_vector<double,3> > mMappingVariables_PointInNewRef;

    std::map<unsigned,c_vector<double,3> > GetMappingVariables_Difference( );
    std::map<unsigned,c_vector<double,3> > mMappingVariables_Difference;


    std::map<unsigned,c_vector<double,3> > GetMappingVariables_P_Translated( );
    std::map<unsigned,c_vector<double,3> > mMappingVariables_P_Translated;

    std::map<unsigned,c_vector<double,3> > GetMappingVariables_Cs( );
    std::map<unsigned,c_vector<double,3> > mMappingVariables_Cs;







    ///////


    // Map to the inital aCoefficients for each element
    std::map<unsigned, c_vector<double, 3> > mACoefficients;
    // Map to the inital aCoefficients for each element
    std::map<unsigned, c_vector<double, 3> > mBCoefficients;
    std::map<unsigned, double> mArea0;
    std::map<unsigned, c_vector<unsigned, 2> > mNearestNodesMap;

    // Make the intial conditions accessable from other classes
    // Access the inital angle for the bending force from other classes
    double GetOriginalAngle(std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edge);

    // Get the inital shape functions and areas and vectors for the shear and area forces

    ///=====

    c_vector<c_vector<double, 2>, 3> GetInitalVectors(unsigned elem_index);
    c_vector<c_vector<double, 3>, 2> GetInitalShapeFunction(unsigned elem_index);
    double GetOriginalArea(unsigned elem_index);

    /*
     * Mark the boundary nodes as boundary nodes -- I have decidede that will be when a node has less than 5 containing elements  
     */
    void MarkBoundaryNodes();
    std::set<unsigned> GetNeighbouringElements(unsigned elem_index);
    std::vector<unsigned> GetCommonNodes(unsigned elem_1, unsigned elem_2);
    c_vector<double, SPACE_DIM> GetElementNormal(unsigned elem_index);

    // Can manually change how many bins in each dim there are from the test
    void SetBinningIntervals(int nX, int nY, int nZ);
    void SetBinningRegions();


    std::map<std::vector<int>, std::vector<unsigned> > mBinMap;

    std::map<std::vector<int>, std::vector<unsigned> > GetAllBins();

    // std::map<unsigned, c_vector<double, SPACE_DIM> > mBin;

    // THis map provides a list of the elements in each bin, the bin being the key 
    std::map<std::vector<int>, std::vector<unsigned > > mBin;
    // Map for the edges in each bin -- each bin has a vector of edges, the edges are descirbed by the nodes on the edge
    std::map<std::vector<int>, std::vector<  std::pair<unsigned, unsigned>  > > mEdgeBin;
    std::map<std::vector<int>, std::vector<double> > mBinCoords;
    // std::map<std::vector<int>,  std::vector<c_vector<double, SPACE_DIM> >> mBinCoords1;

    //Finds the bin of any given point
    std::vector<int> GetBin(c_vector<double, SPACE_DIM> Location);

    // This function is called in SetBinningRegions to set the maximal dimesions of the domain, letting me determine whitch nth of the domain each point is in
    void SetBinningWidth();
    void SetMeshSize();



    // Makes it easier to access each of the boundaries for the bin
    double GetBinLowerX(std::vector<double> Bin);
    double GetBinUpperX(std::vector<double> Bin);
    double GetBinLowerY(std::vector<double> Bin);
    double GetBinUpperY(std::vector<double> Bin);
    double GetBinLowerZ(std::vector<double> Bin);
    double GetBinUpperZ(std::vector<double> Bin);
        

    // Binning member variables
    // Set the number of intervals  -- default is 1, but can be changed with SetBinningIntervals
    int mNx = 1;
    int mNy = 1;
    int mNz = 2;

    // If the geometry is 2d, need to record it  with this member in
    int mDIM = 3;

    // Need max dimensions of the domain
    double mMaxX;
    double mMinX;

    double mMaxY;
    double mMinY;

    double mMaxZ;
    double mMinZ;


    // Get the nearest nodes
    c_vector<unsigned, 2> GetNearestInternalNodes(unsigned node_index);

    // Several methods need the centroids, so here is a method to create a map of the centoroids
    void SetCentroidMap();

    std::map<unsigned, c_vector<double, SPACE_DIM> > mCentroidMap; // Save the centroid for each element

    /**
     * Set up the angles for the bending force -- useing the inital conditions 
     */

    void SetInitialAnlgesAcrossMembrane();
    std::map<std::pair<unsigned, unsigned>, double> mOriginalAngles;

    bool CalculateElementNormals(std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edge,
                                 std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> >& nonUnitNormals,
                                 std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>& otherNodes);

    // Set the element are to set the remeshing element area to be -- determines discretisation of the geometry

    void SetTargetRemeshingEdgeLength(double TargetRemeshingEdgeLength);
    double mTargetRemeshingEdgeLength = 1e-7;
    void SetTargetRemeshingIterations(int Iterations);
    int mIterations = 10;
    void EdgeLengthVariable(double EdgeLengthMultiple);
    bool mVariableEdgeLength =0;
    double mEdgeLengthMultiple =1;

    void SetPrintRemeshedIC(bool PrintRemeshedIC);
    bool mPrintRemeshedIC = 0;


    double CalculateAspectRatio(c_vector<double, SPACE_DIM> Node1, c_vector<double, SPACE_DIM> Node2,c_vector<double, SPACE_DIM> Node3 );
    std::vector<double> MinimumElementAspectRatio();

    std::vector<double> Quantile(std::vector<double>& inData, std::vector<double>& probs);



    void UpdateBoundaryConditions();
    bool GetUpdateBoundaryConditions();
    bool mUpdateComplete =1;
    double GetAspectRatioFromMesh();

    // double GetDistanceToLine( c_vector<double, SPACE_DIM> NewNode, c_vector<double, SPACE_DIM> P1, c_vector<double, SPACE_DIM>  P2);

    double DistanceBetweenPointAndElement( c_vector<double, SPACE_DIM>  NewPoint, unsigned Element, double DistanceToNearestLine);
    double DistanceBetweenPointAndElement( c_vector<double, SPACE_DIM>  NewPoint, unsigned OldElement);

    // Returns the distance and the nearest point 
    std::pair<double, c_vector<double, SPACE_DIM> >  ProjectPointToPlane( c_vector<double, SPACE_DIM>  NewPoint, unsigned Element);
    double ProjectPointToLine(c_vector<double, SPACE_DIM> x1, c_vector<double, SPACE_DIM> x2, c_vector<double, SPACE_DIM>  NewPoint);
    double ProjectPointToLine(std::pair<unsigned, unsigned> edgeIndex , c_vector<double, SPACE_DIM>  NewPoint);


    // Added to solve my mesh problem  .....mChasteOutputDirectory
    // void WriteVtkResultsToFile(const std::string& rDirectory);


    /**
     * Overridden WriteResultsToFiles() method.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     */
    // virtual void WriteResultsToFiles(const std::string& rDirectory);

    // void OutputCellPopulationParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(HistoryDepMeshBasedCellPopulation)

namespace boost
{
    namespace serialization
    {
        /**
         * Serialize information required to construct a HistoryDepMeshBasedCellPopulation.
         */
        template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
        inline void save_construct_data(
            Archive & ar, const HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
        {
            // Save data required to construct instance
            const MutableMesh<ELEMENT_DIM, SPACE_DIM>* p_mesh = &(t->rGetMesh());
            ar & p_mesh;
        }

    /**
     * De-serialize constructor parameters and initialise a HistoryDepMeshBasedCellPopulation.
     * Loads the mesh from separate files.
     */
        template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
        inline void load_construct_data(
            Archive& ar, HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
        {
            // Retrieve data from archive required to construct new instance
            MutableMesh<ELEMENT_DIM, SPACE_DIM>* p_mesh;
            ar >> p_mesh;

            // Invoke inplace constructor to initialise instance
            ::new(t)HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>(*p_mesh);
        }
    } // namespace serialization
} // namespace boost

// #include "SerializationExportWrapper.hpp"
// EXPORT_TEMPLATE_CLASS_ALL_DIMS(HistoryDepMeshBasedCellPopulation)
#endif /*HistoryDepMeshBasedCellPopulation_HPP_*/
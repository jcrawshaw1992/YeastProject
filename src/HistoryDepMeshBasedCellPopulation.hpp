/*

*/

#ifndef HISTORYDEPMESHBASEDCELLPOPULATION_HPP_
#define HISTORYDEPMESHBASEDCELLPOPULATION_HPP_

#include <map>
#include <algorithm>
#include "AbstractCentreBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MutableMesh.hpp"
// #include "VertexMesh.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellsGenerator.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

#include "CellBasedEventHandler.hpp"
#include "CellId.hpp"

#include "NodeVelocityWriter.hpp"
#include "NodesOnlyMesh.hpp"
#include "VoronoiDataWriter.hpp"
#include "VtkMeshWriter.hpp"
#include "UblasCustomFunctions.hpp"
#include "Debug.hpp"
#include "HistoryDepMutableMesh.hpp"
#include "VtkMeshReader.hpp"



/**
 * A facade class encapsulating a mesh-based 'cell population'.
 *
 * Contains a group of cells and maintains the associations between cells and
 * nodes in the mesh.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class HistoryDepMeshBasedCellPopulation : public MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>
{
    // friend class TestMeshBasedCellPopulation;
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    // friend class boost::AbstractCentreBasedCellPopulation::access;
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
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM> >(*this);

        archive & mOriginalNodePositions;
        archive & mInitalPositionOfRemeshedNodes;
        // archive & mpMutableMesh;
        archive & mNew_mesh;

    }

protected:
   
   
    /** Static cast of the mesh from AbstractCellPopulation */
    // HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>* mpMutableMesh;
    HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM> mNew_mesh;

public:
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
                    const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                    bool deleteMesh=false,
                    bool validate=true);

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

     /*
     * Need a direction to the folder everything is savinging to  
     */
    void SetChasteOutputDirectory(std::string ChasteOutputDirectory, double startime);
    std::string mChasteOutputDirectory;

    void ExecuteHistoryDependentRemeshing();

    void MappingAdaptedMeshToInitalGeometry();

     /**
     * Map the new node to the old element   Element_0 has the node position from the inital condition 
     * Element_0 has the node position from the inital condition 
     * DeformedElement has the node position from the deformed position 
     * NewNodes is the position of the new node in the deformed configuration 
     * 
     * Goal is to get the inital position for this remeshed node
     */

    c_vector<double, 3> NewNodeInInitalConfigurationFromChangeOfBasis(unsigned ClosestElement_OldMeshIndex, c_vector<double, SPACE_DIM> NewNode);
    std::map<unsigned, c_vector<double, SPACE_DIM> > mInitalPositionOfRemeshedNodes;


    /**
     * Calling the python code to do the remeshing  
     */
    void RemeshGeometry();

    void SaveInitalConditions();
    std::map<unsigned, c_vector<double, SPACE_DIM> > mOriginalNodePositions;




    /**
     * Overridden WriteResultsToFiles() method.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     */
    // virtual void WriteResultsToFiles(const std::string& rDirectory);

    void OutputCellPopulationParameters(out_stream& rParamsFile);
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
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
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
template<class Archive,  unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    MutableMesh<ELEMENT_DIM, SPACE_DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*HistoryDepMeshBasedCellPopulation_HPP_*/

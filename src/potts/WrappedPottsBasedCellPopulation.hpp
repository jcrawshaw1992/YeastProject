/*

Jess's WrappedWrappedPottsBasedCellPopulation

*/

#ifndef WRAPPEDPOTTSBASEDCELLPOPULATION_HPP_
#define WRAPPEDPOTTSBASEDCELLPOPULATION_HPP_

#include "AbstractOnLatticeCellPopulation.hpp"
// #include "PottsMesh.hpp"
#include "VertexMesh.hpp"
#include "AbstractUpdateRule.hpp"
// #include "AbstractWrappedPottsUpdateRule.hpp"
#include "MutableMesh.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include "UblasCustomFunctions.hpp"
#include <utility>
#include <math.h>
#include <algorithm>

// Added by Jess
#include "PottsArbitrarySurfaceIn3DMesh.hpp"

#include "FixedG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
// #include "projects/VascularRemodelling/src/potts/PottsArbitrarySurfaceIn3DMesh.hpp"


/**
 * A facade class encapsulating a cell population under the Cellular
 * Potts Model framework.
 *
 * Contains a group of cells and maintains the associations
 * between CellPtrs and elements in a specialised PottsMesh class.
 *
 * The code currently requires the PottsMesh object to be fixed,
 * in the sense that no new nodes or elements can be added.
 */
template<unsigned DIM>
class WrappedPottsBasedCellPopulation : public PottsBasedCellPopulation<DIM>
{
    // friend class TestWrappedPottsBasedCellPopulation;

private:

    /** Set the threshold after which the cell breaks in two. Jess set this for her Potts wrapped around the vessel */
    

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
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<PottsBasedCellPopulation<DIM> >(*this);

        /*
         * In its current form the code does not allow the direct serialization
         * of the PottsMesh class, so instead we delete mpVoronoiTessellation.
         */
        // delete mpElementTessellation;
        // mpElementTessellation = nullptr;
        // archive & mTemperature;
         archive & mElementPairing;
         archive & mElementPairingVector;
         archive & mNumSweepsPerTimestep;
         archive & mIsPeriodic;

    }

public:
    


    /**
     * Create a new cell population facade from a mesh and collection of cells.
     *
     * There must be precisely one CellPtr for each PottsElement in
     * the mesh.
     *
     * @param rMesh reference to a PottsMesh
     * @param rCells reference to a vector of CellPtrs
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     *                   (defaults to false)
     * @param validate whether to validate the cell population when it is created (defaults to true)
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    WrappedPottsBasedCellPopulation(PottsMesh<DIM>& rMesh,
                             std::vector<CellPtr>& rCells,
                             std::map< unsigned, unsigned > ElementPairing,
                             bool deleteMesh=false,
                             bool validate=true,
                             const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    WrappedPottsBasedCellPopulation(PottsMesh<DIM>& rMesh,
                             std::vector<CellPtr>& rCells,
                             std::vector<unsigned> BoundaryVector,
                             bool deleteMesh=false,
                             bool validate=true,
                             const std::vector<unsigned> locationIndices=std::vector<unsigned>());                

    std::vector<unsigned> mBoundaryVector;    

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a vertex mesh.
     */
    WrappedPottsBasedCellPopulation(PottsMesh<DIM>& rMesh);

     std::map< unsigned, unsigned > mElementPairing;
     std::vector<c_vector<unsigned,2> > mElementPairingVector;

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~WrappedPottsBasedCellPopulation();

    
    double GetNumSweepsPerTimestep();

    void SetNumSweepsPerTimestep(double numSweepsPerTimestep);

    double mNumSweepsPerTimestep;

    bool mIsPeriodic;

    bool IsPottsSimulationPeriodic();


  



   
    // Jess wrote this to divide the cell. I know there is a divide element function, but it was randomly splitting the cells, so Im writing my own
    //
    CellPtr DivideCell(CellPtr pNewCell, CellPtr pParentCell=CellPtr());
    void DividePopulation();
    void LabelEdges( double NodesWide, double NodesLong);

     // Jess added this function
     // Returns the sister element to any element number 
 
    unsigned GetSister(unsigned Element_Index);

    // Returns the mapping between the paired elements
    std::map< unsigned, unsigned > GetElementPairingMap();
    std::vector<c_vector<unsigned,2> > GetElementPairingVector();

    unsigned GetDivisionStatus(unsigned Element_Index);

    
    void UpdateCellLocations(double dt);

    /**
     * Overridden AddUpdateRule() method.
     *
     * @param pUpdateRule pointer to an update rule
     */

     virtual void AddUpdateRule(boost::shared_ptr<AbstractUpdateRule<DIM> > pUpdateRule);


    //  void AddUpdateRule(boost::shared_ptr<AbstractWrappedUpdateRule<DIM> > pUpdateRule);



    /**
     * Outputs CellPopulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    // void OutputCellPopulationParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(WrappedPottsBasedCellPopulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a WrappedPottsBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const WrappedPottsBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const PottsMesh<DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a WrappedPottsBasedCellPopulation.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, WrappedPottsBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    PottsMesh<DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)WrappedPottsBasedCellPopulation<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*WrappedPottsBasedCellPopulation_HPP_*/

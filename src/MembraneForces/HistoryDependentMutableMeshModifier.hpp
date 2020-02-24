/*

    History dependence in the mutable mesh -- For every node in the new mesh we find the 3 local/closest nodes in the old mesh.
    Make a triangle out of these old nodes. 
    Translate and rotate/map this triangle back to the origin
    Make two of the edges the new basis vectors 
    Describe the Node in the new Mesh P using the basis vectors (needs to be translated toward the origina as the triangle was)
    Find the shape funcitons of this triangle 
    Figure out where these nodes where in the initial configuaration and make an inital triangle 
    project/map this triangle back to the orign
    Find the displacement vector
    Using the dispacement equation V4 = P - po = N1*V1 + N2*V2 + N3*V3, find po, the inital position of the new node in the inital configuarion
    Use two edges of the inital triangle as the basis vectors to discribe po, then unmap the inital triangle to its proper/unmapped configuration

*/

#ifndef HistoryDependentMutableMeshModifier_HPP_
#define HistoryDependentMutableMeshModifier_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <math.h>
#include "AbstractCellBasedSimulationModifier.hpp"
#include "PetscTools.hpp"
#include "UblasCustomFunctions.hpp"
#include "XdrFileReader.hpp"

#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"

/**
 * A modifier class which at each simulation time step
 *
 * To be used in conjunction with the AppliedForce Class
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class HistoryDependentMutableMeshModifier : public AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>
{
private:
    std::map<double, c_vector<long double, 4> > mGrowthMaps;
 

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> >(*this);
        // archive & mSetupSolve;
    }

public:
    /**
	 * Default constructor.
	 */
    HistoryDependentMutableMeshModifier();

    /**
     * Destructor.
     */
    virtual ~HistoryDependentMutableMeshModifier();

    void RemeshingWithHistoryDepenance(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, MutableMesh<ELEMENT_DIM, SPACE_DIM5>& rOldMesh,MutableMesh<ELEMENT_DIM, SPACE_DIM5>& rNewMesh, MutableMesh<ELEMENT_DIM, SPACE_DIM5>& rInitalOldMesh,)
    void NewNodeInInitalConfiguration(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)








    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Boundaries for the basement regions. the normals point into the center of the basement 
     *  
     */
    
    /**
   * Helper method to store the applied tractions in CellData.
   */
    void UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);
  
    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(HistoryDependentMutableMeshModifier);

#endif /*HistoryDependentMutableMeshModifier_HPP_*/

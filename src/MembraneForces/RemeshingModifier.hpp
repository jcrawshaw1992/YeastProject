/*

This code controles the membrame properties. 

At the start of a simulation, the intial conditions delievered at 20(?)X less than original starting configuration of the vessel 
The system then develops to an equlibirum (that will be the realisitc size of a network). After this point the memnbrane properties will 
be variable.

For the inital impletmentation of this code, the membrane properties in the center will be changed to decrease the vessel size to 2X the intial 
condition, with some smoothing on the vessel properties just outside the region. 

After this, the membrane constants will be dependent on the number of cells in the region and the denoted mutation state. 

*/

#ifndef RemeshingModifier_HPP_
#define RemeshingModifier_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <string> 
#include <map>
#include <math.h>
#include "AbstractCellBasedSimulationModifier.hpp"
#include "PetscTools.hpp"
#include "UblasCustomFunctions.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

#include "Debug.hpp"



/**
 * A modifier class which at each simulation time step
 *
 * To be used in conjunction with the AppliedForce Class
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class RemeshingModifier : public AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>
{
private:
 

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
        archive& mOriginalNodePositions;
        archive& mInitalPositionOfRemeshedNodes;
    }

public:
    /**
	 * Default constructor.
	 */
    RemeshingModifier();

    /**
     * Destructor.
     */
    virtual ~RemeshingModifier();

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
     * Finding the inital postion of the adapted mesh   
     */
    void MappingAdaptedMeshToInitalGeometry(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, MutableMesh<2, 3>& rMesh);
    std::map<unsigned, c_vector<double, SPACE_DIM> > mInitalPositionOfRemeshedNodes;
    /**
     * Map the new node to the old element   Element_0 has the node position from the inital condition 
     * Element_0 has the node position from the inital condition 
     * DeformedElement has the node position from the deformed position 
     * NewNodes is the position of the new node in the deformed configuration 
     * 
     * Goal is to get the inital position for this remeshed node
     */
    c_vector<double, 3> NewNodeInInitalConfigurationFromChangeOfBasis(c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> Element_0, c_vector<c_vector<double, SPACE_DIM>, SPACE_DIM> DeformedElement, c_vector<double, SPACE_DIM> NewNode);

    /**
     * Need a direction to the folder everything is savinging to  
     */

    void SetChasteOutputDirectory(std::string ChasteOutputDirectory, double starttime);
    std::string mChasteOutputDirectory;

    /**
     * Calling the python code to do the remeshing  
     */
    void RemeshGeometry();

    void SaveInitalConditions(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation );

    // Map saving Original Node position 
    std::map<unsigned, c_vector<double, SPACE_DIM> > mOriginalNodePositions;
    std::map<unsigned, c_vector<double, SPACE_DIM> > GetInitalConditions();



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
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RemeshingModifier);

#endif /*RemeshingModifier_HPP_*/

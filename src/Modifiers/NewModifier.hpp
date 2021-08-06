
#ifndef NewModifier_HPP_
#define NewModifier_HPP_

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


#include "HistoryDepMeshBasedCellPopulation.hpp"

/**
 * A modifier class which at each simulation time step
 *
 * To be used in conjunction with the AppliedForce Class
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class NewModifier : public AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>
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
    }

public:
    /**
	 * Default constructor.
	 */
    NewModifier();

    /**
     * Destructor.
     */
    virtual ~NewModifier();


void SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory);


void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

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
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NewModifier);

#endif /*APPLIEDFORCEMODIFIER_HPP_*/

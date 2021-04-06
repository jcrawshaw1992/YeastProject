#ifndef RemeshingPopulationWriter_HPP_
#define RemeshingPopulationWriter_HPP_


#include "Debug.hpp"
#include "AbstractCellPopulationWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A class written using the visitor pattern for writing node location from a cell population to file.
 *
 * The output file is called T1SwapLocations.dat by default.
 */

class RemeshingPopulationWriter : public AbstractCellPopulationWriter<2, 3>
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationWriter<2, 3> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    RemeshingPopulationWriter();


    void Visit(HistoryDepMeshBasedCellPopulation<2, 3>* pCellPopulation);

    
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(RemeshingPopulationWriter)


#endif /* RemeshingPopulationWriter_HPP_ */

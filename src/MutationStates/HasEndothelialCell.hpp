

#ifndef HasEndothelialCell_HPP_
#define HasEndothelialCell_HPP_

#include "AbstractCellMutationState.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Subclass of AbstractCellMutationState defining a Beta-catenin with
 * a change at residue 45 mutation state.
 */
class HasEndothelialCell : public AbstractCellMutationState
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell mutation state.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellMutationState>(*this);
    }

public:
    /**
     * Constructor.
     */
    HasEndothelialCell();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(HasEndothelialCell)

#endif /* HasEndothelialCell_HPP_ */

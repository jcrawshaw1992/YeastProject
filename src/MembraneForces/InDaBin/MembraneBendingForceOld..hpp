
#ifndef MembraneBendingForceOld_HPP_
#define MembraneBendingForceOld_HPP_

#include "AbstractForce.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"

#include "Debug.hpp"
#include "UblasCustomFunctions.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"

/**
 * Membrane Stiffness Force!
 */
class MembraneBendingForceOld : public AbstractForce<2, 3>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<2,3> >(*this);
        archive & mMembraneStiffness;
        
    }

public:

    /**
     * Constructor.
     */
    MembraneBendingForceOld();

    void SetMembraneStiffness(double membraneStiffness, double Nc, double Nz);

    double mNc;
    double mNz;
    double mMembraneStiffness;
    std::map<unsigned, c_vector<unsigned, 5> > mNearestNodesMap;

    void SetNearestNodesForBoundaryNodesBending(std::map<unsigned, c_vector<unsigned, 5> > NearestNodesMap);



    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};


// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MembraneBendingForceOld)

#endif /*MembraneBendingForceOld_HPP_*/

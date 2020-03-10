

#ifndef HemodynamicPressure_HPP_
#define HemodynamicPressure_HPP_

#include "AbstractForce.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "MeshBasedCellPopulation.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "PetscTools.hpp"
#include "XdrFileReader.hpp"
#include <map>




/**
 * Membrane Surface Force
 * Force tyring to minimising surface area to the relaxed state 
 */
class HemodynamicPressure : public AbstractForce<2, 3>

{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
   // friend class AbstractTetrahedralMesh;
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
    }


public:

    /**
     * Constructor.
     */
    HemodynamicPressure();


/*
 * This method is called at the begining, specifically from the primary simulation to calculate the areas of each element. This 
 * is then saved as a protected member variable (a map) with the element index as the key, and the area as the vaule. This member 
 * can be easily accessed from the other methods in this class. 
*/

 
    void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation);
    void LoadPressureFromFile(std::string PressureFile);
    std::string mPressureFile;
    double mAppliedPressure; 
    std::vector < c_vector<double,3> > mAppliedPosition;

    std::vector < c_vector<double,3> > mAppliedPressures;

    

     
 
    
    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};


// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(HemodynamicPressure)

#endif /*HemodynamicPressure_HPP_*/

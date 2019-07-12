

#ifndef HemodynamicTraction_HPP_
#define HemodynamicTraction_HPP_

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
class HemodynamicTraction : public AbstractForce<2, 3>

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
    HemodynamicTraction();


/*
 * This method is called at the begining, specifically from the primary simulation to calculate the areas of each element. This 
 * is then saved as a protected member variable (a map) with the element index as the key, and the area as the vaule. This member 
 * can be easily accessed from the other methods in this class. 
*/

 
    void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation);
    void LoadTractionFromFile(std::string TractionFile);

	/**
	 * 	Store the loaded traction data; this is the applied traction at a given point.
	 */
	std::vector < c_vector<double,3> > mAppliedTractions;

	/**
	 * 	Store the loaded tangential component of the traction data; this is the applied tangential traction at a given point.
	 */
	std::vector < c_vector<double,3> > mAppliedTangentTractions;

	/**
	 * 	Store the loaded traction data; this is the location that the traction is defined at.
	 */
	std::vector < c_vector<double,3> > mAppliedPosition;

	bool mResetTractionsOnCells;

	std::string mTractionFile;


    

     
 
    
    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};


// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(HemodynamicTraction)

#endif /*HemodynamicTraction_HPP_*/

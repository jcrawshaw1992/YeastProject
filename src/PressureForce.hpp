

#ifndef PressureForce_HPP_
#define PressureForce_HPP_

#include "AbstractForce.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "MutableVertexMesh.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include <string.h>


/**
 * Membrane Surface Force
 * Force tyring to minimising surface area to the relaxed state 
 */
class PressureForce : public AbstractForce<2, 3>

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
    PressureForce();


/*
 * This method is called at the begining, specifically from the primary simulation to calculate the areas of each element. This 
 * is then saved as a protected member variable (a map) with the element index as the key, and the area as the vaule. This member 
 * can be easily accessed from the other methods in this class. 
*/

    // Add the forces to the nodes 
    void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation);
    void GetConcentrationProfile(AbstractCellPopulation<2,3>& rCellPopulation);
    void SetReactionDiffusionResultsFile(std::string mesh_file);
    void UpdateReactionDiffusionProfile();
    void SetOutputDirectory(std::string outputDirectory);
    void WriteOutVtuFile(std::string OutputDirectory);



/*
 * Member variables. 
*/
    double mLoadPreviousResult = 0;
    std::string mMesh_file;
    std::map<int, double >  mConcentrationAcrossYeast;
    std::vector<double> mConcentrationVector;
    MutableMesh<2, 3> *mMesh;
    std::string mOutputDirectory;

    double mFRMIterationCounter =50;
    double mCounter =1;

    double FEMTimeStep = 1e-3;
    double mCurrentTime = 0;
    double mCurrentEndTime = 1e-2;



     
 
    
    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};


// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(PressureForce)

#endif /*PressureForce_HPP_*/

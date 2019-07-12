

#ifndef MembraneSurfaceForce_HPP_
#define MembraneSurfaceForce_HPP_

#include "AbstractForce.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "MutableVertexMesh.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

// #include "EmptyBasementMatrix.hpp"
// #include "LostEndothelialCell.hpp"
// #include "HasEndothelialCell.hpp"

#include "projects/VascularRemodelling/src/MutationStates/EmptyBasementMatrix.hpp"
#include "projects/VascularRemodelling/src/MutationStates/LostEndothelialCell.hpp"
#include "projects/VascularRemodelling/src/MutationStates/HasEndothelialCell.hpp"




/**
 * Membrane Surface Force
 * Force tyring to minimising surface area to the relaxed state 
 */
class MembraneSurfaceForce : public AbstractForce<2, 3>

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
        archive & mOriginalAreas;
        archive & mMembraneSurface;
        archive & mMembraneSurfaceMaps;
        archive & mScalling;

    }

protected:

    double mMembraneSurface;
    // bool  mUpdateSwitch;
    std::map<unsigned, double> mOriginalAreas;
    std::map<unsigned, double> mMembraneSurfaceMaps;
public:

    /**
     * Constructor.
     */
    MembraneSurfaceForce();


/*
 * This method is called at the begining, specifically from the primary simulation to calculate the areas of each element. This 
 * is then saved as a protected member variable (a map) with the element index as the key, and the area as the vaule. This member 
 * can be easily accessed from the other methods in this class. 
*/
    double mScalling;
   void SetupInitialAreas(AbstractCellPopulation<2,3>& rCellPopulation);

    void UpdateMembraneSurfaceForceProperties(AbstractCellPopulation<2, 3>& rCellPopulation);

   void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation);

   void SetMembraneStiffness(double membrane_constant);

    void SetScallingArea(double Scalling);
//    void SetAreaUpdate(bool number );
    
    // void SetNc(double Nc);


    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};


// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MembraneSurfaceForce)

#endif /*MembraneSurfaceForce_HPP_*/



#ifndef CylinderOutwardsPressure_HPP_
#define CylinderOutwardsPressure_HPP_

#include "AbstractForce.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MutableVertexMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "VertexBasedCellPopulation.hpp"

#include <math.h>
#include "CellLabel.hpp"
#include "Debug.hpp"
// #include "EdgeNodeMutationState.hpp"
#include "UblasCustomFunctions.hpp"


/**
 * Membrane Surface Force
 * Force tyring to minimising surface area to the relaxed state 
 */
class CylinderOutwardsPressure : public AbstractForce<2, 3>

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

    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractForce<2, 3> >(*this);
        archive & mStrength;
    }


       

public:

    double mStrength;
    /**
     * Constructor.
     */
    CylinderOutwardsPressure();

     // Map of the inital areas
    // Might need these later, so just silence 
    void SetPressure(double Pressure);

    void AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CylinderOutwardsPressure)

#endif  /*CylinderOutwardsPressure_HPP_*/

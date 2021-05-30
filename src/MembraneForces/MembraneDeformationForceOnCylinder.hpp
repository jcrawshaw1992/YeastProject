

#ifndef MembraneDeformationForceOnCylinder_HPP_
#define MembraneDeformationForceOnCylinder_HPP_

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
#include "UblasCustomFunctions.hpp"


/**
 * Membrane Surface Force
 * Force tyring to minimising surface area to the relaxed state 
 */
class MembraneDeformationForceOnCylinder : public AbstractForce<2, 3>

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

protected:  
    // double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
    // double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)
    double mStrength =  0.002133152 - 0.001466542;
   

public:
    /**
     * Constructor.
     */
    MembraneDeformationForceOnCylinder();

    

/*
 * Loop over all elements and save the inital membrane confuguration as map, element index as the key, to a vector of vector
 * between each of the nodes 
*/

    c_vector<c_vector<double, 3>, 3> Inverse(c_vector<c_vector<double, 3>, 3> Matrix);
    c_vector<c_vector<double, 3>, 3> Inverse(c_vector<c_vector<double, 3>, 3> Matrix, double elem);

    c_vector<c_vector<double, 3>, 3> MatrixMultiplication(c_vector<c_vector<double, 3>, 3> Matrix1, c_vector<c_vector<double, 3>, 3> Matrix2);
    c_vector<double, 3> MatrixMultiplication(c_vector<c_vector<double, 3>, 3> Matrix, c_vector<double, 3> Vector); 
    c_vector<c_vector<double, 3>, 3> MatrixTranspose(c_vector<c_vector<double, 3>, 3> Matrix1);
    double det(c_vector<c_vector<double, 2>, 2> Matrix);
    double tr(c_vector<c_vector<double, 2>, 2> Matrix);


/*
 * Cast the cell population to the MeshBasedCellPopulation so that I can loop over all of the elements, calulcate the area
 * of each element, then caclulate the area difference and add the area force contribution of this element to each of the 
 * of this element.
*/

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
CHASTE_CLASS_EXPORT(MembraneDeformationForceOnCylinder)

#endif /*MembraneDeformationForceOnCylinder_HPP_*/

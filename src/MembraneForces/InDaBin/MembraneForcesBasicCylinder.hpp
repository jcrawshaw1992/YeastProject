

#ifndef MembraneForcesBasicCylinder_HPP_
#define MembraneForcesBasicCylinder_HPP_

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

#include "projects/VascularRemodelling/src/MutationStates/EmptyBasementMatrix.hpp"
#include "projects/VascularRemodelling/src/MutationStates/LostEndothelialCell.hpp"
#include "projects/VascularRemodelling/src/MutationStates/HasEndothelialCell.hpp"


/**
 * Membrane Surface Force
 * Force tyring to minimising surface area to the relaxed state 
 */
class MembraneForcesBasicCylinder : public AbstractForce<2, 3>

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
        archive & mInitalVectors;
        archive & mACoefficients;
        archive & mBCoefficients;
        archive & mArea0;
    }

protected:


    // double mAreaDilationModulus;
    // double mElasticShearModulus;
    // double mMembraneSurface;
  
    // Map to the inital aCoefficients for each element
    std::map<unsigned, c_vector<double, 3> > mACoefficients;

    // Map to the inital aCoefficients for each element
    std::map<unsigned, c_vector<double, 3> > mBCoefficients;

    std::map<unsigned, double> mArea0;
    
    std::map<unsigned, std::vector<unsigned>  > mNeighbours;
      

    // Map of node1 Iital locations
    //std::map<unsigned,  c_vector<double, 3> > mNode0_IC;
    std::map<unsigned, c_vector<c_vector<double, 2>, 3> > mInitalVectors;
    


public:
    /**
     * Constructor.
     */
    MembraneForcesBasicCylinder();

     // Map of the inital areas
    // Might need these later, so just silence 
    void SetAreaDilationModulus(double AreaDilationModulus);
    void SetElasticShearModulus(double ElasticShearModulus);
    void SetMembraneStiffness(double MembraneSurface);

/*
 * Loop over all elements and save the inital membrane confuguration as map, element index as the key, to a vector of vector
 * between each of the nodes 
*/

    void SetupMembraneConfiguration(AbstractCellPopulation<2, 3>& rCellPopulation);

    c_vector<c_vector<long double, 3>, 3> Inverse(c_vector<c_vector<long double, 3>, 3> Matrix);
    c_vector<c_vector<long double, 3>, 3> Inverse(c_vector<c_vector<long double, 3>, 3> Matrix, double elem);

    c_vector<c_vector<long double, 3>, 3> MatrixMultiplication(c_vector<c_vector<long double, 3>, 3> Matrix1, c_vector<c_vector<long double, 3>, 3> Matrix2);
    c_vector<long double, 3> MatrixMultiplication(c_vector<c_vector<long double, 3>, 3> Matrix, c_vector<long double, 3> Vector); 
    c_vector<c_vector<long double, 3>, 3> MatrixTranspose(c_vector<c_vector<long double, 3>, 3> Matrix1);
    c_vector<c_vector<long double, 3>, 3> RowReduction(c_vector<c_vector<long double, 3>, 3> Matrix);
    c_vector<c_vector<long double, 3>, 3> MappingMatrix(MeshBasedCellPopulation<2, 3>* p_cell_population, typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter, double a, double b, double theta);
    long double det(c_vector<c_vector<long double, 2>, 2> Matrix);
    long double tr(c_vector<c_vector<long double, 2>, 2> Matrix);


/*
 * Cast the cell population to the MeshBasedCellPopulation so that I can loop over all of the elements, calulcate the area
 * of each element, then caclulate the area difference and add the area force contribution of this element to each of the 
 * of this element.
*/

    void AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation);


/*
 * Function that saves a Map of all of the neighbours for all of the nodes 
*/
    void FindNeighbours(AbstractCellPopulation<2, 3>& rCellPopulation);



    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MembraneForcesBasicCylinder)

#endif /*MembraneForcesBasicCylinder_HPP_*/
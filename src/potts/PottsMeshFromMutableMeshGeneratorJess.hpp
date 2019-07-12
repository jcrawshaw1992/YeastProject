
#ifndef PottsMeshFromMutableMeshGeneratorJess_HPP_
#define PottsMeshFromMutableMeshGeneratorJess_HPP_

#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "MutableMesh.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"


// #include "PottsBasedCellPopulationOnVessel.hpp"

/**
 * Generator of arbitrary Potts meshes, used as starting points for many simulations.
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
template<unsigned SPACE_DIM>
class PottsMeshFromMutableMeshGeneratorJess
{
private:
	/** We are only going to consider this case for the time being */
	static const unsigned ELEMENT_DIM = 2;

    /** A pointer to the mesh this class creates */
    PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>* mpMesh;

public:

    /**
     * Constructor to generate a Potts mesh from an existing MutableMesh. The approach is
     * to consider each of the triangles in the mesh as a lattice site in the Potts simulation.
     *
     * @param rMesh the mutable mesh to be used as template
     */
    PottsMeshFromMutableMeshGeneratorJess(MutableMesh<ELEMENT_DIM,SPACE_DIM>& rMesh);

    /**
     * Destructor - deletes the mesh object and pointer
     */
    virtual ~PottsMeshFromMutableMeshGeneratorJess();

    /**
     * @return the generated Potts mesh.
     */
    PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>* GetMesh();
};

//#include "SerializationExportWrapper.hpp"
//EXPORT_TEMPLATE_CLASS1(PottsMeshFromMutableMeshGeneratorJess, 2)
//EXPORT_TEMPLATE_CLASS1(PottsMeshFromMutableMeshGeneratorJess, 3)

#endif /*POTTSMESHFROMMUTABLEMESHGENERATOR_HPP_*/

/*

*/

#ifndef _TETRAHEDRALMESHJESS_HPP_
#define _TETRAHEDRALMESHJESS_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "UblasVectorInclude.hpp"
#include "UblasMatrixInclude.hpp"

#include <vector>
#include <string>
#include <set>

#include "AbstractTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractMeshReader.hpp"
#include "ChastePoint.hpp"


struct triangulateio; /**< Forward declaration for triangle helper methods (used in MutableMesh QuadraticMesh)*/

//////////////////////////////////////////////////////////////////////////
//   DECLARATION
//////////////////////////////////////////////////////////////////////////

/**
 * A concrete tetrahedral mesh class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TetrahedralMeshJess : public TetrahedralMesh< ELEMENT_DIM, SPACE_DIM>
{
   private:
    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the mesh.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<TetrahedralMesh<ELEMENT_DIM, SPACE_DIM> >(*this);
    }


   
public:

    /**
     * Constructor.
     */
    TetrahedralMeshJess();

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader);


    
    };
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(TetrahedralMeshJess)
s
#endif //_TETRAHEDRALMESH_HPP_

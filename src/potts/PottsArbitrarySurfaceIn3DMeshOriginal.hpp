#ifndef POTTSARBITRARYSURFACEIN3DMESH_HPP_
#define POTTSARBITRARYSURFACEIN3DMESH_HPP_

#include "PottsMesh.hpp"
#include "MutableMesh.hpp"

#include <algorithm>
#include "UblasCustomFunctions.hpp"
#include "Debug.hpp"
#include <utility>

#include <math.h>


template<unsigned SPACE_DIM>
class PottsArbitrarySurfaceIn3DMesh : public PottsMesh<SPACE_DIM>
{

public:


    /**
     * Constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param pottsElements vector of pointers to PottsElements
     * @param vonNeumannNeighbouringNodeIndices vector of set of Moore neighbours for each node
     * @param mooreNeighbouringNodeIndices vector of set of Von Neumann neighbours for each node
     * @param pDelaunayMesh pointer to associated underlying Delaunay mesh
     */
    PottsArbitrarySurfaceIn3DMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                  std::vector<PottsElement<SPACE_DIM>*> pottsElements,
                                  std::vector< std::set<unsigned> > vonNeumannNeighbouringNodeIndices,
                                  std::vector< std::set<unsigned> > mooreNeighbouringNodeIndices,
                                  MutableMesh<2,SPACE_DIM>* pDelaunayMesh);

    virtual ~PottsArbitrarySurfaceIn3DMesh();

    /**
     * Get the volume (or area in 2D, or length in 1D) of a lattice site.
     *
     * @param index the global index of a specified lattice site
     * @return the volume of the lattice site
     */
    double GetVolumeOfLatticeSite(unsigned index);

    // std::map<double, double> mMyMap;

    /**
     * Get the volume (or area in 2D, or length in 1D) of a PottsElement.
     *
     * @param index the global index of a specified PottsElement element
     * @return the volume of the element
     */
    double GetVolumeOfElement(unsigned index);

    /**
     * Get the surface area (or perimeter in 2D) of a lattice site.
     *
     * @param latticeSiteIndex the global index of a specified lattice site
     * @return the surface area of the lattice site
     */
    double GetSurfaceAreaOfLatticeSite(unsigned latticeSiteIndex);
    double GetPerimeterOfLatticeSite(unsigned latticeSiteIndex);
    

    /**
     * Compute the surface area (or perimeter in 2D) of a PottsElement.
     *
     * @param index the global index of a specified PottsElement
     * @return the surface area of the element
     */
    double GetSurfaceAreaOfElement(unsigned index);
    double GetPerimeterOfElement(unsigned index);

    /**
     * Compute the contact area (or edge in 2D) between 2 lattice sites.
     *
     * @param index the index of the first lattice site
     * @param index the index of the second lattice site
     * @return the contact area between the sites
     */
    double GetContactAreaBetweenLatticeSite(unsigned index_A, unsigned index_B);

    /**
     * Update the location of the Potts nodes when the underlying Delaunay mesh changes
     */
    void UpdatePottsNodeLocationFromDelaunay();

    /**
     * Get the original triangular mesh used to define the object
     */
    MutableMesh<2,SPACE_DIM>* GetDelaunayMesh();



// protected:



private:

    inline double ComputeTriangleArea(const c_vector<double, SPACE_DIM>& vertexA, const c_vector<double, SPACE_DIM>& vertexB, const c_vector<double, SPACE_DIM>& vertexC) const;

    /**
     * Test if two nodes (i.e. two lattices sites) are contained in the same element (i.e. have the same spin).
     * Also true if both are part of the medium.
     *
     * @param indexNodeA Global index for node a.
     * @param indexNodeB Global index for node b.
     * @return true if both lattice sites are contained in the same element (including the medium)
     */
    inline bool DoNodesShareElement(unsigned indexNodeA, unsigned indexNodeB);

    /** Pointer to original delaunay surface mesh */

    MutableMesh<2,SPACE_DIM>* mpDelaunayMesh;

    std::map<unsigned, double> mMeshElementMidPoints;
};

//#include "SerializationExportWrapper.hpp"
//EXPORT_TEMPLATE_CLASS1(PottsArbitrarySurfaceIn3DMesh, 2)
//EXPORT_TEMPLATE_CLASS1(PottsArbitrarySurfaceIn3DMesh, 3)

#endif /* POTTSARBITRARYSURFACEIN3DMESH_HPP_ */







// template<unsigned SPACE_DIM>
// class PottsArbitrarySurfaceIn3DMesh : public PottsMesh<SPACE_DIM>
// {

//   // protected:
//   // double mMyMember;
//   //     std::map<double, double> mMyMap;

// public:
 
//     /**
//      * Constructor.
//      *
//      * @param nodes vector of pointers to nodes
//      * @param pottsElements vector of pointers to PottsElements
//      * @param vonNeumannNeighbouringNodeIndices vector of set of Moore neighbours for each node
//      * @param mooreNeighbouringNodeIndices vector of set of Von Neumann neighbours for each node
//      * @param pDelaunayMesh pointer to associated underlying Delaunay mesh
//      */
//     PottsArbitrarySurfaceIn3DMesh(std::vector<Node<SPACE_DIM>*> nodes,
//                                   std::vector<PottsElement<SPACE_DIM>*> pottsElements,
//                                   std::vector< std::set<unsigned> > vonNeumannNeighbouringNodeIndices,
//                                   std::vector< std::set<unsigned> > mooreNeighbouringNodeIndices,
//                                   MutableMesh<2,SPACE_DIM>* pDelaunayMesh);

//     virtual ~PottsArbitrarySurfaceIn3DMesh();


//     void SetElasticShearModulus(double ElasticShearModulus);

//     double mElasticShearModulus;
//     std::map<double, double > mACoefficients;


//     /**
//      * Get the volume (or area in 2D, or length in 1D) of a lattice site.
//      *
//      * @param index the global index of a specified lattice site
//      * @return the volume of the lattice site
//      */
//     double GetVolumeOfLatticeSite(unsigned index);
    

//     /**
//      * Get the volume (or area in 2D, or length in 1D) of a PottsElement.
//      *
//      * @param index the global index of a specified PottsElement element
//      * @return the volume of the element
//      */
//     double GetVolumeOfElement(unsigned index);

//     /**
//      * Get the surface area (or perimeter in 2D) of a lattice site.
//      *
//      * @param latticeSiteIndex the global index of a specified lattice site
//      * @return the surface area of the lattice site
//      */
//     double GetSurfaceAreaOfLatticeSite(unsigned latticeSiteIndex);
//     double GetPerimeterOfLatticeSite(unsigned latticeSiteIndex);
    

//     /**
//      * Compute the surface area (or perimeter in 2D) of a PottsElement.
//      *
//      * @param index the global index of a specified PottsElement
//      * @return the surface area of the element
//      */
//     double GetSurfaceAreaOfElement(unsigned index);
//     double GetPerimeterOfElement(unsigned index);

//     /**
//      * Compute the contact area (or edge in 2D) between 2 lattice sites.
//      *
//      * @param index the index of the first lattice site
//      * @param index the index of the second lattice site
//      * @return the contact area between the sites
//      */
//     double GetContactAreaBetweenLatticeSite(unsigned index_A, unsigned index_B);

//     /**
//      * Update the location of the Potts nodes when the underlying Delaunay mesh changes
//      */
//     void UpdatePottsNodeLocationFromDelaunay();

//      /**
//      * Update the stored area and perimeter of the Potts elements when the underlying Delaunay mesh changes
//      */

//     void UpdatePottsVolumeAndArea(MutableMesh<2,SPACE_DIM>& rMesh);

//     /**
//      * Get the original triangular mesh used to define the object
//      */
//     MutableMesh<2,SPACE_DIM>* GetDelaunayMesh();
      
    



//     // std::map<unsigned, c_vector<double , 3>  > mMeshElementMidPoints;
//     // std::map<std::pair<unsigned, unsigned>, double> mOriginalAngles;
    
   
    

// private:
//     // 
//     inline double ComputeTriangleArea(const c_vector<double, SPACE_DIM>& vertexA, const c_vector<double, SPACE_DIM>& vertexB, const c_vector<double, SPACE_DIM>& vertexC) const;

//     /**
//      * Test if two nodes (i.e. two lattices sites) are contained in the same element (i.e. have the same spin).
//      * Also true if both are part of the medium.
//      *
//      * @param indexNodeA Global index for node a.
//      * @param indexNodeB Global index for node b.
//      * @return true if both lattice sites are contained in the same element (including the medium)
//      */
//     inline bool DoNodesShareElement(unsigned indexNodeA, unsigned indexNodeB);

//     /** Pointer to original delaunay surface mesh */
//     MutableMesh<2,SPACE_DIM>* mpDelaunayMesh;


//     /** Needed for serialization. */
//     friend class boost::serialization::access;
//     /**
//      * Archive the object and its member variables.
//      *
//      * @param archive the archive
//      * @param version the current version of this class
//      */
//     template<class Archive>
//     void serialize(Archive & archive, const unsigned int version)
//     {
//         archive & boost::serialization::base_object<PottsMesh<SPACE_DIM> >(*this);
//         // archive & mOriginalAngles;
//     }


// };

// //#include "SerializationExportWrapper.hpp"
// //EXPORT_TEMPLATE_CLASS1(PottsArbitrarySurfaceIn3DMesh, 2)
// //EXPORT_TEMPLATE_CLASS1(PottsArbitrarySurfaceIn3DMesh, 3)

// #endif /* POTTSARBITRARYSURFACEIN3DMESH_HPP_ */

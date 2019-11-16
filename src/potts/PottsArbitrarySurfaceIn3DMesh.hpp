#ifndef POTTSARBITRARYSURFACEIN3DMESH_HPP_
#define POTTSARBITRARYSURFACEIN3DMESH_HPP_



#include "PottsMesh.hpp"
#include "MutableMesh.hpp"
#include "RandomNumberGenerator.hpp"
 #include "WrappedPottsBasedCellPopulation.hpp"

#include "AbstractPottsUpdateRule.hpp"

#include <algorithm>
#include "UblasCustomFunctions.hpp"

#include "Debug.hpp"
#include <utility>
#include <math.h>


// HEADERS NEEDED TO READ IN TRACTION FILE

#include "XdrFileReader.hpp"
#include <vector>
#include <map>
#include "CellLabel.hpp"



#include "MutableElement.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellPopulation.hpp"
        
#include "Cell.hpp"

#include <boost/shared_ptr.hpp>

#include <iterator>


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
     * Overridden UpdateCellLocations() method.
     *
     * @param dt time step
     */



    void SetBoundaries(std::vector<unsigned > BoundaryVector );

    std::vector<unsigned > mBoundaryVector ;
    std::map< std::pair<unsigned,unsigned>  ,unsigned > mIsEdgeBoundaryMap;


    /*
     * Map of the edge lattice sites for each element 
     */
    // std::map< unsigned , std::vector<bool> > mEdgeLatticeSites;


    /*
     * Loops over all the mesh elements and determins the neighboring element pairs. This is to 
     * be run once at the start of the simulation. This member saves the following
     * a) mMapNodesToAssociateElementPairs;
     * b) mMapElementPairsToNodes; 
     * c) mMeshElementPairs; 
     */
    void FindElementNeighbours();
    void GetEdges();


    /*
     * Loops over the edges of a lattice site and calulcates 
     * a) the midpoint of each PottsElement 
     *    This is saved in mMeshElementMidPoints;
     * b) the lenght of each edge between neighbouring mesh elements (the edge of the potts element)
     *    This is saved in mDistanceBetweenElements; 
     * c) calculated the perimeter of the lattice site, this last one might be useless
     * Needs to be ran at the start of each PottsLoop of the coupled model 
     */
    void CalculateEdgeLenghts();
    void EdgeLenghts();



    // // Calculate the curviture between each pair of lattice neigbours ( in vorinno region)
    // void CalculateCurvature();


    // Calculate the traction at each lattice site
    void CalculateTraction();

    // Get the traction data from the files
    void TractionDataLoader(const std::string& tractionFilename);

    // Set a constant wall shear stress -- used to test and validate the Shear stress migration on a 2D lattice 
    void SetConstantWallShearStress(double WallShearStress);

    /*
     * Returns the average traction on a CPM element (has been multipled by the area )
     */
    c_vector<double,3> GetTractionOnElement(unsigned pottsElementIndex);


     /*
     * Returns the sum shear stress mag acting on a Potts element
     */
     double GetSumShearStressOnElement(unsigned pottsElementIndex);

    
    /*
     * Returns the traction on a CPM lattice site  (has been multipled by the area )
     */
    c_vector<double,3> GetTractionOnLattice(unsigned latticeSiteIndex);

    

    




    /*
     * Loops over the enternal edges and sums the lenght to give the perimeter of the element
     */

    double GetPerimeterOfElement(unsigned pottsElementIndex);

    double GetPerimeterOfCoupledElements(unsigned pottsElementIndex, unsigned pottsSisterElementIndex, double Center);

    bool IsNodeNearCenter(unsigned LatticeIndex, double Center );

    void CalculateLatticeVolumes();

 
    double GetVolumeOfElement(unsigned pottsElementIndex);
    double GetNumberOfLatticesInElement(unsigned pottsElementIndex);

    double GetVolumeOfLatticeSite(unsigned NodeIndex);
    bool IsLatticeInElement(unsigned Latticeindex, unsigned pottsElementIndex);




    /*
     * Returns the average mean curvature over the Potts element
     */
    double GetCurvatureOfElement(unsigned pottsElementIndex);


    /*
     * Check the AspectRatio
     */

    double GetAspectRatio(unsigned pottsElementIndex);


     /*
     * Check the Angle between the major axis of the cell and the circumferncial direction 
     */
    double GetMajorAxisAngle(unsigned pottsElementIndex);
  
    void GetRadius(double Radius);

    double mRadius;
    double mMaxC;
    double mMinC;
    double mLength;
    double mPeriodic;
    double mLatticeSpaceing ;
    
    


    std::map< unsigned , c_vector<double,2> > MapCellTo2DPlane(unsigned pottsElementIndex);
    double CellLength(unsigned pottsElementIndex);


    c_vector<double, SPACE_DIM> GetCellCentre(unsigned pottsElementIndex);
    c_vector<double, SPACE_DIM> GetCellCentre(std::vector<unsigned> Edges ,unsigned pottsElementIndex);



//    double min_value(std::vector<double> Vector);
//    double max_value(std::vector<double> Vector);


// bool TestRandomWalk(unsigned pottsElementIndex);

   double mAngleRange;




 // Returns the number of lattice sites in a Potts element 
unsigned GetSizeOfElement(unsigned pottsElementIndex);


    




     /*
     * Returns pointer to node for a given lattice site 
     */
     Node<SPACE_DIM> * GetPottsLatticeSite(unsigned latticeSiteIndex);

    /*
     * Maps the cylinder to the plane so I can gind the aspect ratio in 2D
     */

    void MapCylinderToPlane();

     



    // /*
    //  * ------   MEMBER VARIABLES STORING TRACTION DATA
    //  */ 

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


    /**
	 * 	Store the force on each lattice
	 */

    std::map< unsigned , c_vector<double,3> > mForceOnLattice;

    /**
	 * 	Store the traction on each lattice
	 */

    std::map< unsigned , c_vector<double,3> > mTractionOnLattice;



    // /*
    //  * ------   MEMBER VARIABLES
    //  * 
    //  */ 

    // /* 
    //  * Saves the location of each node as mapped to a 2D plane -- an unwrapped cylinder
    //  */

    std::map<unsigned, c_vector<double, 2> > mMappedLocation;


    // /* 
    //  * Saves the location of each node as mapped to a 2D plane -- an unwrapped cylinder
    //  */

    std::map<unsigned, c_vector<double, 2> > mMajorAxis;



    // /* 
    //  * Maps element to the location of the midpoint  -- Absolute, not relative 
    //  */
    std::map<unsigned, c_vector<double, 3>> mMeshElementMidPoints;

    //  Maps element to all its neighbouring elements 
    // std::map<unsigned, std::vector< unsigned>  > mMeshElementNeighbours;

    /* 
     * Saves the lattice perimeter for each lattice site
     */ 
    std::map<unsigned, double> mLatticePerimeters;
    
    /* 
     * A vector containing all of the element Pairs -- this is set once at the start  
     */ 
      
    std::vector< std::pair <unsigned,unsigned> > mMeshElementPairs; 

    /* 
     * A map mapping the each element pair to the two common nodes (aka the edge connecting them) -- Only needs to be set once at the start
     */ 
      
    std::map< std::pair <unsigned,unsigned> , std::vector<unsigned>  > mMapElementPairsToNodes; 
    std::map< unsigned , std::pair <unsigned,unsigned> > mEdges; 
    std::map< unsigned ,  std::vector<unsigned> > mGetElementsFromEdge;
    std::map< unsigned , double > mIsEdgeBoundary;

    // Will return if a pair of nodes (making an edge) is on the boundary or not
    std::map< std::pair <unsigned,unsigned> , bool  >mGetNodePairOnBoundary;


        double mBoundaryEdgeUnitLength;// the length between two nodes along the top and bottom boundary in a regular mesh 


    /* 
     * Saves the distance between the center points of each paired element. Remember this is a distance in 3d along a surface. 
     */ 
    
    std::map< std::pair <unsigned,unsigned> , double > mDistanceBetweenElements;

    std::map< std::vector<unsigned> , double > mPerimeterBetweenLatticeSites;

/// THe lattice edge is basically perpendicular to the mesh edge, so this is the length between the center of the two elements bisecting the mesh edges
    std::map< std::pair <unsigned,unsigned> , double > mLengthOfLatticeEdge;


     /* 
     * Maps two neighbouring Nodes (pair) to the two associated elements (pair). 
     */ 
    std::map< std::pair <unsigned,unsigned> , std::pair <unsigned,unsigned>  >  mMapNodesPairsToElementPairs;
    std::map< std::pair <unsigned,unsigned> , std::vector<unsigned>  >  mMapNodesPairsToElementVec;


     /* 
     * Edge lattice sites of the cell  
     */ 
    // std::map< unsigned , std::vector <unsigned>  >  mEdgeLatticeSites;




    /* 
     * Maps the neighbouring nodes to the instantaneous curvature between them . 
     */ 
    
    std::map< std::pair <unsigned,unsigned> , double > mInstantaneousCurvature;

    /* 
     * Angle between neighbouring elements. 
     */ 



    std::map< std::pair <unsigned,unsigned> , double >  mAngleBetweenNodes;


    /* 
     * Maps the neighbouring nodes to the instantaneous curvature between them . 
     */ 
    
    std::map< unsigned , double > mMeanCurvature;


    /* 
     * Map of the lattice sites to the lattice volumes 
     */ 

    std::map<unsigned, double> mLatticeVolume;

    std::map<unsigned, c_vector<double, 3>> mLatticeNormal;
    

    /* 
     * Maps each node index will all of the asscoaited element pairs. Kinda the inverse mapping to  mMapElementPairsToNodes.
     */ 
      
    std::map< unsigned, std::vector< std::pair <unsigned,unsigned>> > mMapNodesToAssociateElementPairs; 

     /* 
     * Maps for each potts element - saves a vector list of each lattice site, labelling them as 1 or 0 if they are an edge lattice or not respectivly.
     */ 
      
    std::map< unsigned, std::map< unsigned, bool > > mEdgeLatticeSitesList; 



     /* 
     * Adds angles together -- periodic 
     */ 

    // double AddAngles(double alpha, double beta);
   
    double MiddleAngle(std::vector<double> Vector);

    double mQuaters;

// Average edge length;
    double mAverageEdgeLength;

    // Average area of the elements
    double mAverageArea;







    // Ensure that the unit normal is outwards pointing 

    // double MaintainOutwardsPointingNormal( c_vector<double, 3> Normal,  c_vector<double, 3> x1);

     /* 
     * Create a pair, Order is important, so orders the smallest first 
     */ 

    // std::pair <unsigned,unsigned> Create_pair(unsigned x, unsigned y );
    // std::pair <double,double> Create_pair(double x, double y );


    // Axciallary functions 

    // void PRINT_PAIR(std::pair<unsigned, unsigned> Pair);


         /* 
     * Returns a vector with the common elements from each vector 
     */ 
     
    // std::vector<unsigned> Intersection(std::vector<unsigned> Vector1, std::vector<unsigned> Vector2);

    /* 
     * Checks if a pair is already in the vector  
     */ 
    
    // bool IsPairInVector(std::vector<std::pair<unsigned, unsigned >> Vector, std::pair<unsigned, unsigned > Pair);

    /* 
     * Checks is an unsigned or double is already in the vector 
     */ 
     
    // bool IsNumberInVector(std::vector< unsigned > Vector, unsigned number);
    // bool IsNumberInVector(std::vector< double > Vector, double number);
    bool IsVectorInVector(std::vector< c_vector<double,SPACE_DIM>  > Vector, c_vector<double,SPACE_DIM> Locaiton );

    /* 
     * Check if the vectors are the same 
     */ 
      
    double AreVectorsSame(std::vector<unsigned> Vector1, std::vector<unsigned> Vector2);
    double AreVectorsSame(std::vector<std::pair<unsigned, unsigned >> Vector1, std::vector<std::pair<unsigned, unsigned >> Vector2);

    /* 
     * Removes the number from the vector  
     */ 
    // std::vector<double> RemoveElement(std::vector<double> Vector1, double number);
    // std::vector<unsigned> RemoveElement(std::vector<unsigned> Vector1, unsigned number);



    /* 
     * Function to remove the interal edges of a cell, leaving only the external edges to then iterate over and calulate the perimeter 
     */ 
     std::vector<std::pair<unsigned, unsigned> > RemoveInternalEdges(std::vector<std::pair<unsigned, unsigned> > ElementPairsAssociatedWithPottsCell);

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
     * Unwrap a cylinder so i have a 2D rectangle easily 
     */
    void UnwrapCylinder(double radius);



    /**
     * This corrects the volume and edge length of the mesh elements near the edges of the mesh for the lazy unwrapped cylinder version of a periodic rectangle
     */
    void CorrectedPeriodicEdges();




    /**
     * Get the original triangular mesh used to define the object
     */
    MutableMesh<2,SPACE_DIM>* GetDelaunayMesh();



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
    
    

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<PottsMesh<SPACE_DIM> >(*this);
        archive & mAngleRange;
        archive & mRadius;
        archive & mMaxC;
        archive & mMinC;
        archive & mLength;
        archive & mMappedLocation;
        archive & mMajorAxis;
        archive & mMeshElementMidPoints;
        archive & mLatticePerimeters;
        archive & mMeshElementPairs; 
        archive & mMapElementPairsToNodes; 
        archive & mDistanceBetweenElements;
        archive & mPerimeterBetweenLatticeSites;
        archive &  mMapNodesPairsToElementPairs; 
        archive & mAngleBetweenNodes;
        archive & mMeanCurvature;
        archive & mLatticeVolume;
        archive & mLatticeNormal;
        archive & mMapNodesToAssociateElementPairs; 
        archive & mEdgeLatticeSitesList; 
        archive & mpDelaunayMesh;
        archive & mLatticeSpaceing;
        archive & mEdges;
    }

};


// #include "SerializationExportWrapper.hpp"
// EXPORT_TEMPLATE_CLASS1(PottsArbitrarySurfaceIn3DMesh, 2)
// EXPORT_TEMPLATE_CLASS1(PottsArbitrarySurfaceIn3DMesh, 3)

#endif /* POTTSARBITRARYSURFACEIN3DMESH_HPP_ */






// Functions that I have removed 

// TestRandomWalkAlongEdges
// 


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
    
   


// };

// //#include "SerializationExportWrapper.hpp"
// //EXPORT_TEMPLATE_CLASS1(PottsArbitrarySurfaceIn3DMesh, 2)
// //EXPORT_TEMPLATE_CLASS1(PottsArbitrarySurfaceIn3DMesh, 3)

// #endif /* POTTSARBITRARYSURFACEIN3DMESH_HPP_ */

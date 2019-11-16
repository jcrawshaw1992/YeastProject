
#ifndef ARBITRARYSURFACEAREACONSTRAINTUPDATERULE_HPP_
#define ARBITRARYSURFACEAREACONSTRAINTUPDATERULE_HPP_


#include "PottsArbitrarySurfaceIn3DMesh.hpp"


#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractPottsUpdateRule.hpp"
#include "PottsBasedCellPopulation.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "CellLabel.hpp"

/**
 * An area surface constraint update rule class for use in Potts based simulations.
 *
 * Note this currently assumes cells don't grow, i.e the target surface area is constant
 * for each cell over time.
 */
template<unsigned DIM>
class ArbitrarySurfaceAreaConstraintPottsUpdateRule : public AbstractPottsUpdateRule<DIM>
{
friend class TestPottsUpdateRules;

private:

    /**
     * Cell deformation energy parameter.
     * Set to the default value 0.5 in the constructor.
     * \todo provide units
     */
    double mSurfaceAreaEnergyParameter;

    /**
     * Non-dimensional target volume of a mature (fully-grown) cell,
     * given in number of lattice sites.
     * Set to the default value 16 in the constructor.
     */
    double mTargetSurfaceArea;

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractPottsUpdateRule<DIM> >(*this);

        archive & mSurfaceAreaEnergyParameter;
        archive & mTargetSurfaceArea;

    }

public:

    /**
     * Constructor.
     */
    ArbitrarySurfaceAreaConstraintPottsUpdateRule();

    /**
     * Destructor.
     */
    ~ArbitrarySurfaceAreaConstraintPottsUpdateRule();

    /**
     * Overridden EvaluateHamiltonianContribution() method
     *
     * Uses sum_elements alpha (V_i - V_i^T)^2.
     *
     * @param currentNodeIndex The index of the current node/lattice site
     * @param targetNodeIndex The index of the target node/lattice site
     * @param rCellPopulation The cell population
     *
     * @return The difference in the Hamiltonian with the configuration of the target node
     * having the same spin as the current node with the current configuration. i.e H_1-H_0
     */
    double EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                           unsigned targetNodeIndex,
                                           PottsBasedCellPopulation<DIM>& rCellPopulation);


    /*
     * Loops over the enternal edges and sums the lenght to give the perimeter of the element
     */

    // double GetSurfaceAreaOfElement(unsigned pottsElementIndex, PottsElement<DIM>* pCurrentElement );
       

    /*
     * Loops over the edges of a lattice site and calulcates 
     * a) the midpoint of each PottsElement 
     *    This is saved in mMeshElementMidPoints;
     * b) the lenght of each edge between neighbouring mesh elements (the edge of the potts element)
     *    This is saved in mDistanceBetweenElements; 
     * c) calculated the perimeter of the lattice site, this last one might be useless
     * Needs to be ran at the start of each PottsLoop of the coupled model 
     */
    void CalculateEdgeLenghts(MutableMesh<2,DIM>& rMesh);

    /*
     * Loops over all the mesh elements and determins the neighboring element pairs. This is to 
     * be run once at the start of the simulation. This member saves the following
     * a) mMapNodesToAssociateElementPairs;
     * b) mMapElementPairsToNodes; 
     * c) mMeshElementPairs; 
     */
    void FindElementNeighbours(MutableMesh<2,DIM>& rMesh);




    /**
     * Set mDeformationEnergyParameter.
     *
     * @param deformationEnergyParameter the new value of mDeformationEnergyParameter
     */
    void SetSurfaceAreaEnergyParameter(double surfaceAreaEnergyParameter);

    
    
    /**
     * @return mMatureCellTargetSurfaceArea
     */
    double GetTargetSurfaceArea() const;

    /**
     * Set mMatureCellTargetSurfaceArea.
     *
     * @param matureCellTargetSurfaceArea the new value of mMatureCellTargetSurfaceArea
     */
    void SetTargetSurfaceArea(double TargetSurfaceArea);





                                           
     /* 
     * Returns a vector with the common elements from each vector 
     */ 
     
    std::vector<unsigned> Intersection(std::vector<unsigned> Vector1, std::vector<unsigned> Vector2);

    /* 
     * Checks if a pair is already in the vector  
     */ 
    
    bool IsPairInVector(std::vector<std::pair<unsigned, unsigned >> Vector, std::pair<unsigned, unsigned > Pair);

    /* 
     * Checks is an unsigned or double is already in the vector 
     */ 
     
    bool IsNumberInVector(std::vector< unsigned > Vector, unsigned number);
    bool IsNumberInVector(std::vector< double > Vector, double number);

    /* 
     * Check if the vectors are the same 
     */ 
      
    double AreVectorsSame(std::vector<unsigned> Vector1, std::vector<unsigned> Vector2);
    double AreVectorsSame(std::vector<std::pair<unsigned, unsigned >> Vector1, std::vector<std::pair<unsigned, unsigned >> Vector2);

    /* 
     * Function to remove the interal edges of a cell, leaving only the external edges to then iterate over and calulate the perimeter 
     */ 
     std::vector<std::pair<unsigned, unsigned> > RemoveInternalEdges(std::vector<std::pair<unsigned, unsigned> > ElementPairsAssociatedWithPottsCell);

    /* 
     * A function testing the RemoveInternalEdges function. Will throw an error if RemoveInternalEdges dosnt give the expected result 
     */ 
    void TestRemoveInternalEdgesFunction();
    
    /* 
     * This function is supposed to determine the target perimeter for a cell x wide and y long
     */ 

    double DetermineAReasonableTargetPerimeter(double x, double y);



    /**
     * @return mDeformationEnergyParameter
     */
    double GetSurfaceAreaEnergyParameter();

    /* 
     * Maps element to the location of the midpoint  -- Absolute, not relative 
     */
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


    /* 
     * Saves the distance between the center points of each paired element. Remember this is a distance in 3d along a surface. 
     */ 
    
    std::map< std::pair <unsigned,unsigned> , double > mDistanceBetweenElements;
    

    /* 
     * Maps each node index will all of the asscoaited element pairs. Kinda the inverse mapping to  mMapElementPairsToNodes.
     */ 
      
    std::map< unsigned, std::vector< std::pair <unsigned,unsigned>> > mMapNodesToAssociateElementPairs; 


    
    // Prints the two parts of a pair
    void PRINT_PAIR(std::pair<unsigned, unsigned> Pair);

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ArbitrarySurfaceAreaConstraintPottsUpdateRule)

#endif /*ARBITRARYSURFACEAREACONSTRAINTUPDATERULE_HPP_*/

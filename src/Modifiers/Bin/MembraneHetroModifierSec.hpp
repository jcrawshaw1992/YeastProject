/*

This code controles the membrame properties. 

At the start of a simulation, the intial conditions delievered at 20(?)X less than original starting configuration of the vessel 
The system then develops to an equlibirum (that will be the realisitc size of a network). After this point the memnbrane properties will 
be variable.

For the inital impletmentation of this code, the membrane properties in the center will be changed to decrease the vessel size to 2X the intial 
condition, with some smoothing on the vessel properties just outside the region. 

After this, the membrane constants will be dependent on the number of cells in the region and the denoted mutation state. 

*/

#ifndef MembraneHetroModifierSec_HPP_
#define MembraneHetroModifierSec_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "PetscTools.hpp"
#include "XdrFileReader.hpp"
#include <map>
#include <math.h>
#include <algorithm>
#include "UblasCustomFunctions.hpp"
#include <iostream>


#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"


/**
 * A modifier class which at each simulation time step
 *
 * To be used in conjunction with the AppliedForce Class
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MembraneHetroModifierSec : public AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>
{
private:

  std::map<double, c_vector<long double, 4> > mGrowthMaps;  
  bool mOn;
  double mMaxZ;
  double mMinZ;
  double mKe;
  unsigned mSamplebasementNode;
  unsigned mSampleECNode;
  bool mAchievedTargetK;

  double mMaxBasementZ;
  double mMinBasementZ;

   double mStrength;
   bool mHetro;

   c_vector<double, 6> mEdgeOfBasementPlane; //A(x-a)+B(y-b)+C(z-z)=0
   c_vector<double, 6> mEdgeOfECPlane;

   std::vector<unsigned> mNodesNextToBasement;
   std::vector<unsigned> mBasementNodes;
   double mDistanceBetweenPlanes;

   c_vector<double, 3> mCentroidEC;
   c_vector<double, 3> mCentroidBasement;
   c_vector<double, 3> mNormalToPlane;

   c_vector<double, 3> mBoundaryPoint;
   c_vector<double, 3> mBoundaryNormal;


   std::vector<unsigned> mStoredBoundaryDefiningNodes;




   double mCounter;


	/** Needed for serialization. */
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
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mGrowthMaps;
        archive & mOn;
        archive & mMaxZ;
        archive & mMinZ;
        archive & mKe;
        archive & mSamplebasementNode;
        archive & mSampleECNode;
        archive & mAchievedTargetK;

        archive & mMaxBasementZ;
        archive & mMinBasementZ;

        archive & mNodesNextToBasement;
        archive & mBasementNodes;

        archive & mDistanceBetweenPlanes;
        archive & mEdgeOfBasementPlane; //A(x-a)+B(y-b)+C(z-z)=0
        archive & mEdgeOfECPlane;
        archive & mCentroidEC;
        archive & mCentroidBasement;
        archive & mNormalToPlane;
        archive & mCounter;

        archive & mStoredBoundaryDefiningNodes;
    }

public:

    /**
	 * Default constructor.
	 */
    MembraneHetroModifierSec();

    /**
     * Destructor.
     */
    virtual ~MembraneHetroModifierSec();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolveBasedOnGivenRegionOrNodes(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory);
    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory);
    void SetAConstantHetrogenaity(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory);





	/**
	 * Set mResetTractionsOnCells.
	 *
	 * Set the member variable mGrowthMaps to be the paramters sent in 
	 */

  void  SetMembranePropeties( std::map<double, c_vector<long double, 4> > GrowthMaps, double Strength, bool Hetrogeneous);

  void SetBoundaryNode(unsigned Node_Index);
  void SetBoundary( c_vector<long double, 3>  Point, c_vector<long double, 3>  normal);

  /**
   * Helper method to store the applied tractions in CellData.
   */
  void UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MembraneHetroModifierSec);

#endif /*MembraneHetroModifierSec_HPP_*/

/*

This code controles the membrame properties. 

At the start of a simulation, the intial conditions delievered at 20(?)X less than original starting configuration of the vessel 
The system then develops to an equlibirum (that will be the realisitc size of a network). After this point the memnbrane properties will 
be variable.

For the inital impletmentation of this code, the membrane properties in the center will be changed to decrease the vessel size to 2X the intial 
condition, with some smoothing on the vessel properties just outside the region. 

After this, the membrane constants will be dependent on the number of cells in the region and the denoted mutation state. 

*/

#ifndef MembranePropertiesModifier_HPP_
#define MembranePropertiesModifier_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <math.h>
#include "AbstractCellBasedSimulationModifier.hpp"
#include "PetscTools.hpp"
#include "UblasCustomFunctions.hpp"
#include "XdrFileReader.hpp"

#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"

/**
 * A modifier class which at each simulation time step
 *
 * To be used in conjunction with the AppliedForce Class
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MembranePropertiesModifier : public AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>
{
private:
    std::map<double, c_vector<long double, 4> > mGrowthMaps;
    // mGrowthMaps[Strength] = [ KA, Kalpha, Ks, kb];
    mGrowthMaps[10] = Create_c_vector(pow(10,-6.9), pow(10,-8.2459),pow(10, -9), 1e-14 );
    mGrowthMaps[8] = Create_c_vector(pow(10,-6.9), pow(10,-8.0160),pow(10, -9) , 1e-14 );
    mGrowthMaps[6] = Create_c_vector(pow(10,-6.9), pow(10,-7.7300),pow(10, -9) , 1e-14 );


    // mGrowthMaps[5] = Create_c_vector(pow(10,-6.9341)*1e-3, pow(10,-7.7)*1e-2,pow(10, -8) *1e-2, 1e-14 );

    mGrowthMaps[5] = Create_c_vector(pow(10,-6.9341), pow(10,-7.7),pow(10, -8), 1e-14 );
    mGrowthMaps[4] = Create_c_vector(pow(10,-6.9), pow(10,-7.4224),pow(10, 8), 1e-14  );

    mGrowthMaps[2] = Create_c_vector(pow(10,-6.8), pow(10,-6.8124),pow(10, -7) , 1e-14 );
    mGrowthMaps[1.5] = Create_c_vector(pow(10,-6.5), pow(10,-6.3491),pow(10, -7) , 1e-14 );
    mGrowthMaps[1.2] =  Create_c_vector(pow(10,-6.2), pow(10,-5.8360),pow(10, -7), 1e-14  );



    bool mOn =0;
    double mMaxZ;
    double mMinZ;

    unsigned mSamplebasementNode;
    unsigned mSampleECNode;
    bool mAchievedTargetK;

    double mStrength;
    bool mHetro;
    double mStepSize;
    double mCounter;
    double mThreshold;
    double mSetupSolve;

    std::vector<unsigned> mNodesNextToBasement;
    std::vector<unsigned> mBasementNodes;

    std::vector<std::vector<  c_vector<double, 3> > > mBoundaries;



    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive& mGrowthMaps;
        archive& mOn;
        archive& mMaxZ;
        archive& mMinZ;

        archive& mSamplebasementNode;
        archive& mSampleECNode;
        archive& mAchievedTargetK;

        archive& mNodesNextToBasement;
        archive& mBasementNodes;
        archive& mStepSize;
        archive& mCounter;
        archive& mBoundaries;
        // archive & mSetupSolve;
    }

public:
    /**
	 * Default constructor.
	 */
    MembranePropertiesModifier();

    /**
     * Destructor.
     */
    virtual ~MembranePropertiesModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Boundaries for the basement regions. the normals point into the center of the basement 
     *  
     */
    
    void Boundaries(c_vector<double, 3> UpperPlaneNormal, c_vector<double, 3> UpperPlanePoint, c_vector<double, 3> LowerPlaneNormal, c_vector<double, 3> PlanePoint);


    /**
	 * Set mResetTractionsOnCells.
	 *
	 * Set the member variable mGrowthMaps to be the paramters sent in 
	 */

    void SetMembranePropeties(std::map<double, c_vector<long double, 4> > GrowthMaps, double Strength, bool Hetrogeneous, double StepSize, double SetupSolve);

    void SetBendingForce(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, double BendingConstant);

    /**
   * Helper method to store the applied tractions in CellData.
   */
    void UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);
    void StepChange(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MembranePropertiesModifier);

#endif /*MembranePropertiesModifier_HPP_*/

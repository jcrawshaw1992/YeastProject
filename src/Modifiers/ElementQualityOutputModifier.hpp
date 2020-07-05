/*

This code controles the membrame properties. 

At the start of a simulation, the intial conditions delievered at 20(?)X less than original starting configuration of the vessel 
The system then develops to an equlibirum (that will be the realisitc size of a network). After this point the memnbrane properties will 
be variable.

For the inital impletmentation of this code, the membrane properties in the center will be changed to decrease the vessel size to 2X the intial 
condition, with some smoothing on the vessel properties just outside the region. 

After this, the membrane constants will be dependent on the number of cells in the region and the denoted mutation state. 

*/

#ifndef ElementQualityOutputModifier_HPP_
#define ElementQualityOutputModifier_HPP_

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

#include <iostream>
#include <fstream>

/**
 * A modifier class which at each simulation time step writes out dimensions of elements that I want 
 *
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ElementQualityOutputModifier : public AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>
{
private:

 double mExecute;
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
    }

public:
    /**
	 * Default constructor.
	 */
    ElementQualityOutputModifier();

    /**
     * Destructor.
     */
    virtual ~ElementQualityOutputModifier();


    ofstream mOutputFile;
    std::string mFileName;
    std::string mAreaFileName;
    std::string mInternalAnglesFileName;
    std::string mAspectRatioFileName;

    std::string mRadiusRatioFileName;
    std::string mEdgeLengthsFileName;
    std::string mMeanChangeInLocalCurvatureFileName;
    std::string mElementRadiFileName;
    std::string mReadMe;

    std::string mEdgeAnglefile;
    std::string mLocalNodesForEdege;
    std::string mContainingElementsForEdges;
    std::string mRegionofContainingEdgesForEdges;
    std::string mTimeFiles;
    double mOutPutCounter;
    double mTimeCounter;



    double mWriteoutInterval = 200;
    double mWriteoutIntervalCounter = 0;

    void SetElementMetricsFileName(std::string FileName);
    void SetWriteoutInterval(double WriteoutInterval);

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);
    void WriteOutDataToFiles(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

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
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ElementQualityOutputModifier);

#endif /*ElementQualityOutputModifier_HPP_*/

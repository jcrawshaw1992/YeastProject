#ifndef HemeLBForce_HPP_
#define HemeLBForce_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractMesh.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "CellLabel.hpp"

#include <boost/filesystem.hpp>
#include <ctime>
#include <boost/algorithm/string/predicate.hpp>
#include <cstdlib>
#include <filesystem>
// namespace fs = std::filesystem;

// #include <direct.h>
#include <boost/filesystem.hpp>
#include <mpi.h>
#include <unistd.h>  
#include <fstream>

#include <iostream>
#include <copyfile.h>

#include <cstdio>
#include <memory>
#include <string>
#include <array>
#include <utility>

#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h> // Need to read the centerlines file

#include <algorithm>
#include <vector>
#include <cmath>

#include "HistoryDepMutableMesh.hpp"
#include "VtkMeshWriter.hpp" // Need a mesh writer -- writing out the current mesh as a vtu so HemeLB can read
#include <sstream>
#include "MutableMesh.hpp"

// #include "PetscTools.hpp"
#include "XdrFileReader.hpp"
#include <map>



#include "vtkXMLUnstructuredGridReader.h"

#include <vtkPolyData.h>


#include "AppliedForceModifier.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "Debug.hpp"
#include <math.h>



#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"

#include "AbstractCellBasedSimulationModifier.hpp"
#include "PetscTools.hpp"
#include "XdrFileReader.hpp"
#include <map>

#include "EmptyBasementMatrix.hpp"
#include "LostEndothelialCell.hpp"
#include "HasEndothelialCell.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"


/**
 * A force class to be used with an HemeLBForceOffLatticeSimulation to impose an externally defined force.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class HemeLBForce  : public AbstractForce<ELEMENT_DIM, SPACE_DIM>
{
friend class TestForces;

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

    
public:

    /**
     * Constructor.
     */
    HemeLBForce();

    /**
     * Destructor.
     */
    ~HemeLBForce();

    /**
     * Overridden AddForceContribution() method. Which uses the applied force stored in CellData.
     *
     * @param rCellPopulation reference to the cell population
     *
     */
    void AddForceContribution(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);


    void Writepr2File(std::string Directory, double SimulationDuration);
    std::string mRadiusdataName = "Radius";
    double mRadius;


    void WriteHemeLBBashScript();
    void CopyFile(std::string InputDirectory, std::string OutputDirectory);
    void UpdateCurrentyFlowVtuCount();
    int mCenterlinesNumber = 1;
    bool mRunHemeLB = 1;
    bool mSetupHemeLB = 1;
    

    void WriteOutVtuFile(std::string outputDirectory);
    HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM> *mMesh;


    void SetUpHemeLBConfiguration(std::string outputDirectory, HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>& Mesh, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);
    
    void ExecuteHemeLB();
    bool CheckIfSteadyStateAchieved();
    void ReRunHemeLB();
    void LoadTractionFromFile();
    void UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);



    void LoadTractionFromVTKFile();
    void UpdateCellData2(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);


    double mLatestFinialHemeLBVTU = -1;
    double mExecuteHemeLBCounter = 0;
    std::map<unsigned, c_vector<double, 3>> mForceMap;


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
	 * 	Store the loaded traction data; this is the location that the traction is defined at.
	 */
	std::vector <double > mAppliedPressure;

    


    // Inlets and out lets -- this will let me write the pr2 file 
    void Inlets(c_vector<double, 3> PlaneNormal, c_vector<double, 3> Point, double pressure, std::string FlowDirection);

    std::vector<std::vector<c_vector<double, 3>>> mIolets; // Each vector has one boundary, this is a vecotr of two vectors, none being the normal and the other is the point defining the plane

    std::vector<double> mPressure;
    std::vector<std::string> mType;
    double mHemeLBScalling = 1;//1e3;
    

    std::string mOutputDirectory;
    std::string mHemeLBDirectory;
    std::string mHemeLB_output;
    std::string mhemelb_setup_exe = "env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui";


    // Function to help return system commands 
    
    std::pair<std::string, int> exec(const char* cmd);


    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(HemeLBForce)

#endif /*HemeLBForce_HPP_*/

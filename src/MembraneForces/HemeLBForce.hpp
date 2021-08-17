#ifndef HemeLBForce_HPP_
#define HemeLBForce_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractMesh.hpp"
#include "CellLabel.hpp"
#include <boost/filesystem.hpp>// Needed here to avoid serialization errors (on Boost<1.37)
#include <ctime>
#include <boost/algorithm/string/predicate.hpp>
#include <cstdlib>
// #include <direct.h>
#include <boost/filesystem.hpp>
#include <mpi.h>
#include <unistd.h>  
#include <fstream>

#include <iostream>


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
#include <iomanip>
#include "MutableMesh.hpp"

// #include "PetscTools.hpp"
#include <map>
#include "vtkXMLUnstructuredGridReader.h"
#include <vtkPolyData.h>
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
#include "UblasCustomFunctions.hpp"


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
    int mCenterlinesNumber = 164;
    bool mRunHemeLB = 1;
    bool mSetupHemeLB = 1;
    double mConstantPressure =0;
    bool mNewInlets =1;
    double mRemeshingCounter =0;

    c_vector<c_vector<double, 3> , 4> mCollapsedRegion;
    

    void SetConstantPressure(double Pressure);

    void WriteOutVtuFile(std::string outputDirectory); 
    void WriteOpenVtus(int Period, int mCenterlinesNumber);


    HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM> *mMesh;

    void SetUpHemeLBConfiguration(std::string outputDirectory,  AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);
    void SetUpHemeLBConfiguration(std::string outputDirectory,  AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, bool RunInitalHemeLB);
    void SetUpFilePaths(std::string outputDirectory, bool CreateFiles, bool RenamePriorResults);

    void ExecuteHemeLB();
    bool CheckIfSteadyStateAchieved();
    void ReRunHemeLB();
    void LoadTractionFromFile();
    void UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    void SetStartTime(double StartTime);

    double mStartTime ;
    double GetStartTime();




    // void LoadTractionFromVTKFile();
    // void UpdateCellData2(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);


    double mLatestFinialHemeLBVTU = -1;
    double mExecuteHemeLBCounter = 0;
    double mTriggerHemeLB = 1500;
    
    std::string mResultsDirectory;
    std::string mRemoveResultsDirectory;
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
  


   void  Network(std::string Network);
   std::string  mNetwork = "Honeycomb";
   double mMinSS; 
   double mMaxSS; 
   double mRegionOfForceCollection = 0.0015;

    


    // Inlets and out lets -- this will let me write the pr2 file 
    void Inlets(c_vector<double, 3> PlaneNormal, c_vector<double, 3> Point, double pressure, std::string FlowDirection);
    void CollapsedRegions(c_vector<double, 3> UpperPlaneNormal, c_vector<double, 3> UpperPoint, c_vector<double, 3> LowerPlaneNormal, c_vector<double, 3> LowerPoint);

    void SetFluidSolidIterations(double Iterations);
    std::vector<std::vector<c_vector<double, 3>>> mIolets; // Each vector has one boundary, this is a vecotr of two vectors, none being the normal and the other is the point defining the plane

    std::vector<double> mPressure;
    double mEstimatedIC;
    double mExpectedVelocity;
    std::vector<std::string> mType;
    double mHemeLBScalling = 1;//1e3;
    
    std::string mChasteOutputDirectory;
    std::string mOutputDirectory;
    std::string mHemeLBDirectory;
    std::string mHemeLB_output;
    // std::string mhemelb_setup_exe = "env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui";
    std::string mhemelb_setup_exe = "env PYTHONPATH=/home/vascrem/hemelb-dev/Tools/setuptool:$PYTHONPATH /home/vascrem/hemelb-dev/Tools/setuptool/scripts/hemelb-setup-nogui";
    std::string mHemeLBPath = "/home/vascrem/hemelb-dev/";
    void SetHemeLBPath(std::string HemeLBPath);

    void SetGenerateFlowVtus( bool FlowVtus);
    bool mFlowVtus=0;

    void SetMachine(std::string Machine);
    std::string mMachine ="server"; // Machine can be mac or Linux server, will make this better soon 

    // Function to help return system commands     
    std::pair<std::string, int> exec(const char* cmd);
    
    // Dont lose precision when going number to string   
    std::string double_to_string(double Number, long double precision);


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









/*
********************
   Code graveyard 
********************
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::UpdateCellData2(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{

	assert(SPACE_DIM==3); // Currently assumes that SPACE_DIM = 3
	std::map<unsigned, c_vector<unsigned, 2>  > LatticeToNodeMap;
	MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);

	for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
		 cell_iter != rCellPopulation.End();
		 ++cell_iter)
	{
		c_vector<double, SPACE_DIM> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
		unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
		Node<SPACE_DIM>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);
		unsigned nearest_fluid_site = UNSIGNED_UNSET;
		double distance_to_fluid_site = DBL_MAX;
	
		for (unsigned fluid_site_index = 0; fluid_site_index <  mAppliedPosition.size(); fluid_site_index++)
		{
			// Find the closest fluid site 
			double distance = norm_2(location - mAppliedPosition[fluid_site_index]);
			if (distance < distance_to_fluid_site)
			{
				distance_to_fluid_site = distance;	
				nearest_fluid_site = fluid_site_index;
			}
		}
		PRINT_2_VARIABLES(distance_to_fluid_site, nearest_fluid_site);
		LatticeToNodeMap[node_index] = Create_c_vector(nearest_fluid_site, 0 );
		assert(nearest_fluid_site != UNSIGNED_UNSET);

		c_vector<double, 3> NormalVector = Create_c_vector(0,0,0);
		std::set<unsigned>& containing_elements = pNode->rGetContainingElementIndices();
        assert(containing_elements.size() > 0);
		// Finding the normal to the node -- need to check this 
        for (std::set<unsigned>::iterator iter = containing_elements.begin();
            iter != containing_elements.end();
            ++iter)
        {
            Node<SPACE_DIM>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
            Node<SPACE_DIM>* pNode1 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
            Node<SPACE_DIM>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

            c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
            c_vector<double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

            NormalVector  += VectorProduct(vector_13, vector_12);
        }
		NormalVector /=norm_2(NormalVector);
		double Pressure = mAppliedPressure[nearest_fluid_site];    
		c_vector<long double,3> Force = Pressure * NormalVector; 
    
		assert(fabs(Force[0])<1e10);
		assert(fabs(Force[1])<1e10);
		assert(fabs(Force[2])<1e10);
		

		// Store the force in CellData
		cell_iter->GetCellData()->SetItem("Pressure", Pressure);
		cell_iter->GetCellData()->SetItem("applied_force_x", Force[0]);
		cell_iter->GetCellData()->SetItem("applied_force_y", Force[1]);
		cell_iter->GetCellData()->SetItem("applied_force_z", Force[2]);
    	}
}





template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::LoadTractionFromVTKFile()
{
    //   std::string TractionFile = mHemeLBDirectory + "results/Extracted/surface-tractions.xtr";
      std::string file = mHemeLBDirectory + "results/Extracted/surface-pressure_3348.vtu";
    
        
      vtkSmartPointer<vtkXMLUnstructuredGridReader> Reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
      Reader->SetFileName(file.c_str());
      Reader->Update();

      vtkUnstructuredGrid* Output = Reader->GetOutput();
      vtkPointData* point_data = Reader->GetOutput()->GetPointData();

    //   double NumberOfDataPoints = Reader->GetOutput()->GetNumberOfPoints();

    //   // Write all of the coordinates of the points in the vtkPolyData to the console.
    //   for(int i = 0; i < Output->GetNumberOfPoints(); i++)
    //   {
    //     double p[3];
    //     Output->GetPoint(i,p);
    //     c_vector<double,3> Point;
    //     Point[0]= p[0]; Point[1]= p[1]; Point[2]= p[2];
    //     mAppliedPosition.push_back(Point*1e3);
    //   }



      // PRINT_VARIABLE(NumberOfDataPoints)
    //   std::cout << *point_data << std::endl;
    //   std::cout << *Output << std::endl;
      // vtkPointData* point_data = Reader->GetOutput()->GetPointData();


    //   This will get the fluid property at each point -- still need the corrds 
      vtkCellData *cellData = Output->GetCellData();

    //   vtkDataSet *DataSet = Reader->GetOutputAsDataSet();

    //   std::cout << *cellData << std::endl;
    //   std::cout << *DataSet << std::endl;

      
    //   for (int i = 0; i < cellData->GetNumberOfArrays(); i++)
    //   {
          vtkDataArray* data = cellData->GetArray(0);
        //   cout << "name " << data->GetName() << endl;
          for (int j = 0; j < data->GetNumberOfTuples(); j++)
          {
              double value = data->GetTuple1(j);
            //   cout << "  value " << j << "th is " << value << endl;
              mAppliedPressure.push_back(value);

          }
    //   }
    PRINT_2_VARIABLES(data->GetNumberOfTuples(), Output->GetNumberOfPoints())
    // assert(mAppliedPosition.size() == Output->GetNumberOfPoints());
    // assert(mAppliedTractions.size() == number_fluid_sites);
    // assert(mAppliedTangentTractions.size() == number_fluid_sites);


}

*/
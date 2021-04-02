#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"

#include "Debug.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "SmartPointers.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"

#include "CommandLineArguments.hpp"
#include "FixedRegionBoundaryCondition.hpp"

#include "MembraneDeformationForce.hpp"
#include "OutwardsPressure.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:
    void TestParametersOverCylinder() throw(Exception)
    {


        std::string ArchivedDirectory = "ParameterSweepWithRemeshing/Cylinder/Param_10_DilationParam_10_DeformationParam_5.5/archive";
        

                // Load and fix any settings in the simulator
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(ArchivedDirectory, 10);



        double AreaParameter = 10; //CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-AreaParameter");
        double DilationParameter = 10;//CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-DilationParameter");

        double dt= 0.01;

        double startime = 0;
        double EndTime = 10;

        double SamplingTimestepMultiple = 500;
  
        std::string OperatingSystem = "Mac";

        // Now I have done the first one, going to iterate over all the others



        // ///////////
        // // std::string output_dir = ArchivedDirectory;
        // double DeformationParamter[9] = {6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10};
        // for (unsigned D_index = 0; D_index < 9; D_index++)
        // {

        //     std::stringstream outN;
        //     outN << "Param_" << AreaParameter << "_DilationParam_" << DilationParameter << "_DeformationParam_" << DeformationParamter[D_index];
        //     std::string ParameterSetN = outN.str();
        //     std::string output_dir = "ParameterSweepWithRemeshing/Cylinder/"+ParameterSetN;


        //     // Load and fix any settings in the simulator
        //     PRINT_VARIABLE(ArchivedDirectory)
        //     PRINT_VARIABLE(EndTime)
        //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(ArchivedDirectory, EndTime);

               
        //     ArchivedDirectory = output_dir;
        //     startime = EndTime;
        //     EndTime = EndTime +10;

            
        //     /* Update the ouput directory for the population  */
        //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, startime);
        //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(startime);

         

        //     p_simulator->SetEndTime(EndTime);
        //     p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        //     p_simulator->SetDt(dt);
        //     p_simulator->SetOutputDirectory(output_dir);

        
        //     /* 
        //     -----------------------------
        //     Update membrane properties
        //     ----------------------------
        //     */
        //     std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        //     boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnHeteroMeshModifier<2, 3> >(*iter);
        //     std::map<double, c_vector<long double, 4> > GrowthMaps;

        //     GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter[D_index]), 0);
        //     //                                 Strength,hetro,stepsize, setupsolve
        //     p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 0, 100, 1);
                

        //     p_simulator->Solve();
        //     CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);


        // }

    
    }
};

#endif /*TESTRELAXATION_HPP_*/

#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"

// #include <cstdio>
// #include <ctime>
// #include <cmath>
// #include <vector>

#include "Debug.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "SmartPointers.hpp"
#include "VtkMeshReader.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"


#include "FixedRegionBoundaryCondition.hpp"
#include "HemeLBForce.hpp"
#include "OutwardsPressureWithBreaks.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"
#include "CylinderOutwardsPressure.hpp"

#include "MembraneDeformationForceOnCylinder.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:
    void offTestRunningArchieve() throw(Exception)
    {

        TRACE("Jess is best")
     
        double AreaParameter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-AreaParameter");
        double DilationParameter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-DilationParameter");
        double DeformationParamter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-DeformationParamter");
       

        PRINT_3_VARIABLES(AreaParameter, DilationParameter, DeformationParamter)
        double dt= 0.0005;
         if (CommandLineArguments::Instance()->OptionExists("-dt"))
        {
            dt= CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-dt");
        }
        double NewEndTime = 30;
        if (CommandLineArguments::Instance()->OptionExists("-NewEndTime"))
        {
            NewEndTime = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-NewEndTime");
        }
        double EndTime = 10;
        if (CommandLineArguments::Instance()->OptionExists("-EndTime"))
        {
            EndTime = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-EndTime");
        }
        double SamplingTimestepMultiple = 5000;
        if (CommandLineArguments::Instance()->OptionExists("-SamplingTimestepMultiple"))
        {
            SamplingTimestepMultiple = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-SamplingTimestepMultiple");
        }

        std::string ArchivedDirectory = "MembraneParameterSweep/Cylinder";
        if (CommandLineArguments::Instance()->OptionExists("-ArchivedDirectory"))
        {
            ArchivedDirectory  = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-ArchivedDirectory");
        }
        
         double StartingParameterForSlowIncrease = 1e-9;
        if (CommandLineArguments::Instance()->OptionExists("-StartingParameterForSlowIncrease"))
        {
            StartingParameterForSlowIncrease= CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-StartingParameterForSlowIncrease");
        }
        // PRINT_4_VARIABLES(dt, NewEndTime,EndTime, SamplingTimestepMultiple)
        std::stringstream out;
        out << "Param_" << AreaParameter << "_DilationParam_" << DilationParameter << "_DeformationParam_" << DeformationParamter;
        std::string ParameterSet = out.str();
        std::string output_dir = ArchivedDirectory + "/Parameters/";

        // Load and fix any settings in the simulator
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(ArchivedDirectory, 10);

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        p_simulator->RemoveAllForces();
        p_simulator->SetEndTime(EndTime + NewEndTime);
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir + ParameterSet);

        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForceOnCylinder> p_shear_force(new MembraneDeformationForceOnCylinder());
        p_simulator->AddForce(p_shear_force);

        /* 
        -----------------------------
        Update membrane properties
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnHeteroMeshModifier<2, 3> >(*iter);
        std::map<double, c_vector<long double, 4> > GrowthMaps;
        GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
        // Strength,hetro,stepsize, setupsolve
        // GrowthMaps, Strength, Hetrogeneous,  StepSize,SetupSolve
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 0, 100, 1);

        if (AreaParameter < (double)7 || DilationParameter < (double)7 || DeformationParamter < (double)7)
        {
            p_Mesh_modifier->SetStartingParameterForSlowIncrease(StartingParameterForSlowIncrease);
            p_Mesh_modifier->SetSlowIncreaseInMembraneStrength(1, 1);
        }
        p_simulator->Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

    }



    void TestRunningArchieveExtended() throw(Exception)
    {

        TRACE("Jess is best")
        double DilationParameter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-DilationParameter");
        double AreaParameter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-AreaParameter");
        double DeformationParamter =CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-DeformationParamter");
    
        PRINT_3_VARIABLES(AreaParameter, DilationParameter, DeformationParamter)
        double dt= CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-dt");
        double NewEndTime = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-NewEndTime");;
        double EndTime = 50;//CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-EndTime");
      
        double SamplingTimestepMultiple = 1000;
        if (CommandLineArguments::Instance()->OptionExists("-SamplingTimestepMultiple"))
        {
            SamplingTimestepMultiple = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-SamplingTimestepMultiple");
        }

        std::stringstream out;
        out << "Param_" << AreaParameter << "_DilationParam_" << DilationParameter << "_DeformationParam_" << DeformationParamter;
        std::string ParameterSet = out.str();
        std::string output_dir = "MembraneParameterSweep/Cylinder/Parameters/" + ParameterSet ;

        // Load and fix any settings in the simulator
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(output_dir, EndTime);

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        p_simulator->SetEndTime(EndTime + NewEndTime);
      
        p_simulator->Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        // SimulationTime::Instance()->Destroy();
        // SimulationTime::Instance()->SetStartTime(0.0);

    }

};

#endif /*TESTRELAXATION_HPP_*/

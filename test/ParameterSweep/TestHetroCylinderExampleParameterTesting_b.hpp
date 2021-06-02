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
// #include "VtkMeshReader.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"


#include "FixedRegionBoundaryCondition.hpp"
#include "MembraneDeformationForceOnCylinder.hpp"
#include "OutwardsPressureWithBreaks.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:
 
    void TestIntroduceHetro1() throw(Exception)
    {

        std::string Archieved = "HetroCylinderExampleThird/";
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, 25);
        // Load and fix any settings in the simulator

        double scale = 1e3;
        c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(0, 0, 1);
        // c_vector<double, 3> UpperPlanePoint = Create_c_vector(12e-6 * scale, 0, 0);
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0, 0, 15e-6 * scale);
        c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0, 0, -1);
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0, 0,25e-6 * scale);

        double DilationParameter=8.4;
        double AreaParameter=7.5;
        double DeformationParamter=8.3;

     
        double dt= 0.001;
        double NewEndTime = 100;
        double EndTime = 25;
        
        double SamplingTimestepMultiple = 1000;
        std::string output_dir = "HetroCylinder_b/1/";
        
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        p_simulator->SetEndTime(EndTime + NewEndTime);
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);

        /* 
        -----------------------------
        Update membrane properties
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnHeteroMeshModifier<2, 3> >(*iter);
        std::map<double, c_vector<long double, 4> > GrowthMaps;
        GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
        GrowthMaps[0] = Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -5), 0);
        // Strength,hetro,stepsize, setupsolve
        // GrowthMaps, Strength, Hetrogeneous,  StepSize,SetupSolve
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 1, 1e-10, 1);
        p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
        p_Mesh_modifier->SetBasementMembraneStrength(0);
        p_Mesh_modifier->SetPlateauParameters(8, 1);

        p_simulator->Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);


    }

     void TestIntroduceHetro2() throw(Exception)
    {

        std::string Archieved = "HetroCylinderExampleThird/";
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, 25);
        // Load and fix any settings in the simulator

        double scale = 1e3;
        c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(0, 0, 1);
        // c_vector<double, 3> UpperPlanePoint = Create_c_vector(12e-6 * scale, 0, 0);
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0, 0, 15e-6 * scale);
        c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0, 0, -1);
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0, 0,25e-6 * scale);

        double DilationParameter=8.4;
        double AreaParameter=7.5;
        double DeformationParamter=8.3;

     
        double dt= 0.001;
        double NewEndTime = 100;
        double EndTime = 25;
        
        double SamplingTimestepMultiple = 1000;
        std::string output_dir = "HetroCylinder_b/4/";
        
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        p_simulator->SetEndTime(EndTime + NewEndTime);
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);

        /* 
        -----------------------------
        Update membrane properties
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnHeteroMeshModifier<2, 3> >(*iter);
        std::map<double, c_vector<long double, 4> > GrowthMaps;
        GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
        GrowthMaps[0] = Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -5), 0);
        // Strength,hetro,stepsize, setupsolve
        // GrowthMaps, Strength, Hetrogeneous,  StepSize,SetupSolve
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 1, 1e-10, 1);
        p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
        p_Mesh_modifier->SetBasementMembraneStrength(0);
        p_Mesh_modifier->SetPlateauParameters(8, 4);

        p_simulator->Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);


    }







};

#endif /*TESTRELAXATION_HPP_*/

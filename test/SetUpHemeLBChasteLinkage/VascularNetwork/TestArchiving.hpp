#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"

#include "Debug.hpp"

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


// #include "FixedRegionBoundaryCondition.hpp"
// #include "EnclosedRegionBoundaryCondition.hpp"
// #include "HemeLBForce.hpp"
// #include "MembraneDeformationForce.hpp"
// #include "OutwardsPressureWithBreaks.hpp"
// #include "OutwardsPressure.hpp"
// #include "RemeshingTriggerOnStepHeteroModifier.hpp"
// #include "MembraneBendingForce.hpp"




class TestRemeshing : public AbstractCellBasedTestSuite
{
public:

 void TestWithConstantForce() throw(Exception)
   {

        double AreaParameter = -5;  double DilationParameter = -5.5; double DeformationParamter = -5; double BendingParameter = -7;
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), pow(10, BendingParameter)) }, {0,  Create_c_vector(pow(10, -4), pow(10, -4), pow(10, -4),pow(10, BendingParameter))} };

        TRACE("Jess is good")
        double EndTime = 0;
        double scale = 0.0011; 

        double SamplingStep = 100;
        double dt = 0.02;
        double RemeshingTime = 500;
        double EdgeLength =1.5*scale;
        
        /////////////////////////////////////////////////////////////////////////////////////
        std::string output_dir =  "FSISimulations/VascularNetwork/GrowingToEqui/ConstantForceArchiving/";
        std::string mesh_file = "/data/vascrem/MeshCollection/VascularNetwork/VascularNetwork.vtu";
        
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(scale, scale, scale);// Scale z by 0.0009 soon 
        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step


        // Create the cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);

        cell_population.SetChasteOutputDirectory(output_dir, 0);
        cell_population.SetInitialAnlgesAcrossMembrane(); // Dont worry about this for now, I think there is something moff
        cell_population.SetRelativePath(output_dir, 0);
        cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
        cell_population.SetBinningIntervals(10, 10, 1);
        cell_population.EdgeLengthVariable(1.2);
        cell_population.SetPrintRemeshedIC(1);
        cell_population.SetTargetRemeshingIterations(10);
        cell_population.SetWriteVtkAsPoints(false);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.SetRemeshingSoftwear("CGAL");
        cell_population.SetOperatingSystem("server");
        // cell_population.ExecuteHistoryDependentRemeshing();
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(SamplingStep);
        simulator.SetDt(dt);
        simulator.SetUpdateCellPopulationRule(false);


        TRACE("First Solve ")
        PRINT_VARIABLE(EndTime)
        cell_population.SetStartTime(EndTime);
        EndTime += 0.5;
        simulator.SetEndTime(EndTime);

        simulator.Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
   

    }



void TestContinuedWithConstantForce() throw(Exception)
   {

        std::string Archieved ="FSISimulations/VascularNetwork/GrowingToEqui/ConstantForceArchiving/";

        double EndTime = 0.5;

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);

        
        // for (int i =0; i<=3; i++)
        //     { 
        
        //         EndTime +=10;
        //         p_simulator->SetEndTime(EndTime);

        //         p_simulator->Solve();
        //         CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        //     }

    }





};

#endif /*TESTRELAXATION_HPP_*/





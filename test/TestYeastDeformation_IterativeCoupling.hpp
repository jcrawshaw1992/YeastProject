#ifndef TESTCYLINDERVALIDATION_HPP_
#define TESTCYLINDERVALIDATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <ctime>
#include <cmath>
#include "OffLatticeSimulation.hpp"
// #include "HoneycombMeshGenerator.hpp"
// #include "CylindricalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "MeshBasedCellPopulation.hpp"

// #include "FixedRegionBoundaryCondition.hpp"

#include "SmartPointers.hpp"
#include "CellLabel.hpp"
#include "VtkMeshWriter.hpp"
#include "Debug.hpp"
#include "CommandLineArguments.hpp"
// #include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "PressureForce.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <boost/filesystem.hpp>

static const double M_TIME_FOR_SIMULATION = 0.010; //50

class YeastDeformation : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }


    public:

    void TestBasicYeastDeformation() //throw (Exception)
    {
    
            std::string mOutputDirectory = "/data/vascrem/testoutput/YeastDeformation/FirstIterationTest/";
           
            /* 1) Set up conditions for FEM */
            std::string mesh_file = "/home/vascrem/Chaste/projects/YeastProject/reaction_diffusion_deforming_membrane/Meshes/YeastMesh.vtu";
            std::string Initial_U = mOutputDirectory+ "Updated_u.xml";
            std::string Initial_V = mOutputDirectory+ "Updated_u.xml";

            TRACE(" Step 1: Attain initial concentration profile ")
            std::string HemeLBCommand =  "cd projects/YeastProject/;./YeastBash ";
            std::string waitFile = mOutputDirectory + "WaitFile.txt";
            HemeLBCommand +=mesh_file+" "+Initial_U+" "+Initial_V+" "+mOutputDirectory + " 0 0.01 0.001 0 "+waitFile;

            // int SystemOutput = std::system(HemeLBCommand.c_str()); 

            // while(! boost::filesystem::exists(waitFile))
            // {
            //     TRACE("waiting within C")
            //     sleep(2); 
            // }
            // // remove(waitFile.c_str())
            // boost::filesystem::remove(waitFile.c_str());

        
            std::string output_directory = "YeastDeformation/FirstIterationTest/";

            // This data file is in mm
            VtkMeshReader<2,3> mesh_reader(mesh_file);
            MutableMesh<2,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            // double scaling = 1e-3;  // so distances are in m

            // mesh.Scale(scaling,scaling,scaling);

            // Create cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
            cell_population.SetWriteVtkAsPoints(true);
            // cell_population.SetOutputMeshInVtk(true);

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellProliferativeTypesWriter>();
            cell_population.AddCellWriter<CellMutationStatesWriter>();

            cell_population.CalculateRestLengths();

            // Set up cell-based simulation
            OffLatticeSimulation<2,3> simulator(cell_population);
            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(100);
            simulator.SetDt(0.01);
            simulator.SetSamplingTimestepMultiple(10);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.

            // Create a force law and pass it to the simulation
                                        // GeneralisedLinearSpringForce
            boost::shared_ptr<GeneralisedLinearSpringForce<2,3> > p_linear_force(new GeneralisedLinearSpringForce<2,3>());
            p_linear_force->SetMeinekeSpringStiffness(0.010);
            simulator.AddForce(p_linear_force);

            

            // Create a force law to apply radial pressure force
            double pressure = 0.01; // to match 80mmhg

            // MAKE_PTR_ARGS(PressureForce, p_radial_force, (pressure));
            // simulator.AddForce(p_radial_force);

            boost::shared_ptr<PressureForce > p_radial_force(new PressureForce());
            // p_radial_force->SetMeinekeSpringStiffness(0.010);
            simulator.AddForce(p_radial_force);
            
            simulator.Solve();
            // To reset before looping: this is usually done by the SetUp and TearDown methods
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        
    }

};

#endif /*TESTCYLINDERVALIDATION_HPP_*/


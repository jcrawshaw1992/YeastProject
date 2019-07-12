#ifndef TESTANGLEFORCE_HPP_
#define TESTANGLEFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <ctime>
#include <cmath>

#include "AppliedForceModifier.hpp"
#include "SpringLengthModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
// #include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "LinearSpringWithRestLengthDependentSpringConstantsForce.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "AppliedForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VtkMeshWriter.hpp"
#include "MembraneStiffnessForce.hpp"
#include "Debug.hpp"


#include "PetscSetupAndFinalize.hpp"



class TestRelaxation : public AbstractCellBasedTestSuite
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

    void TestArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "MembraneStiffnessForce.arch";


        TrianglesMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/open_ended_box");
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        {
        	MembraneStiffnessForce force;
        	force.SetupInitialMembrane(mesh);
            force.SetMembraneStiffness(1.5);

            // Check member variables have been set correctly
            TS_ASSERT_EQUALS(force.GetMembraneStiffness(), 1.5);

            //Check some of the angles are set correctly

            // On 90 degree edge
            std::pair<Node<3>*, Node<3>*> edge = std::pair<Node<3>*, Node<3>*>(mesh.GetNode(0), mesh.GetNode(1));
            TS_ASSERT_DELTA(force.GetOriginalAngle(edge), M_PI/2.0,1e-6);

            // On flat edge
            edge = std::pair<Node<3>*, Node<3>*>(mesh.GetNode(1), mesh.GetNode(3));
            TS_ASSERT_DELTA(force.GetOriginalAngle(edge), 0.0,1e-6);

            // No common edge
            edge = std::pair<Node<3>*, Node<3>*>(mesh.GetNode(0), mesh.GetNode(2));
            TS_ASSERT_THROWS(force.GetOriginalAngle(edge),std::out_of_range);

            // On boundary edge
            edge = std::pair<Node<3>*, Node<3>*>(mesh.GetNode(0), mesh.GetNode(3));
            TS_ASSERT_DELTA(force.GetOriginalAngle(edge), DOUBLE_UNSET,1e-6);


            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer to most abstract class possible
            AbstractForce<2,3>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractForce<2,3>* p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;

            // Check member variables have been correctly archived
            TS_ASSERT_EQUALS(dynamic_cast<MembraneStiffnessForce*>(p_force)->GetMembraneStiffness(), 1.5);

            //Check some of the angles are set correctly

            // On 90 degree edge
            std::pair<Node<3>*, Node<3>*> edge = std::pair<Node<3>*, Node<3>*>(mesh.GetNode(0), mesh.GetNode(1));
            TS_ASSERT_DELTA(dynamic_cast<MembraneStiffnessForce*>(p_force)->GetOriginalAngle(edge), M_PI/2.0,1e-6);

            // On flat edge
            edge = std::pair<Node<3>*, Node<3>*>(mesh.GetNode(1), mesh.GetNode(3));
            TS_ASSERT_DELTA(dynamic_cast<MembraneStiffnessForce*>(p_force)->GetOriginalAngle(edge), 0.0,1e-6);

            // No common edge
            edge = std::pair<Node<3>*, Node<3>*>(mesh.GetNode(0), mesh.GetNode(2));
            TS_ASSERT_THROWS(dynamic_cast<MembraneStiffnessForce*>(p_force)->GetOriginalAngle(edge),std::out_of_range);

            // On boundary edge
            edge = std::pair<Node<3>*, Node<3>*>(mesh.GetNode(0), mesh.GetNode(3));
            TS_ASSERT_DELTA(dynamic_cast<MembraneStiffnessForce*>(p_force)->GetOriginalAngle(edge), DOUBLE_UNSET,1e-6);

            // Tidy up
            delete p_force;
        }

    }

    void  TestOpenEndedBox() throw (Exception)
    {
        double applied_pressure = 5.0;
        double spring_constant = 15.0;

        TrianglesMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/open_ended_box");

        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;    // Jess updated this line because this version of chaste doesnt seem to have the old cell generator, might, i didnt really check ---- CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.CalculateRestLengths();

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory("TestOpenEndedBox");
        simulator.SetEndTime(10.0); //20.0
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetUpdateCellPopulationRule(false);

        // Create a force law and pass it to the simulation
        boost::shared_ptr<GeneralisedLinearSpringForce<2,3> > p_force(new GeneralisedLinearSpringForce<2,3>());
        p_force->SetMeinekeSpringStiffness(spring_constant);
        simulator.AddForce(p_force);

        boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        p_membrane_force->SetupInitialMembrane(mesh);
        simulator.AddForce(p_membrane_force);


        //perturb the box
        simulator.rGetCellPopulation().GetNode(0)->rGetModifiableLocation()[0] +=0.5;
        simulator.rGetCellPopulation().GetNode(0)->rGetModifiableLocation()[1] +=0.5;
        simulator.rGetCellPopulation().GetNode(0)->rGetModifiableLocation()[2] +=0.5;

        simulator.Solve();

        //Check nothings moved
        for (AbstractCellPopulation<2,3>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            c_vector<double,3> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
            cell_location(2)=0.0;
            TS_ASSERT_DELTA(norm_2(cell_location), 1.5, 1e-5);

        }
    }

    void  TestDeformSphere() throw (Exception)
    {
        double spring_constant = 15.0;

        TrianglesMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/sphere");

        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.CalculateRestLengths();

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeformSphere");
        simulator.SetEndTime(20.0); //20.0
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetUpdateCellPopulationRule(false);

        // Create a force law and pass it to the simulation
        boost::shared_ptr<GeneralisedLinearSpringForce<2,3> > p_force(new GeneralisedLinearSpringForce<2,3>());
        p_force->SetMeinekeSpringStiffness(spring_constant);
        simulator.AddForce(p_force);

        boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        p_membrane_force->SetupInitialMembrane(mesh);
        simulator.AddForce(p_membrane_force);

        // Solid base Boundary condition
        boost::shared_ptr<PlaneBoundaryCondition<2,3> > p_bcs(new PlaneBoundaryCondition<2,3>(&cell_population, unit_vector<double>(3,2), -unit_vector<double>(3,2)));
        p_bcs->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bcs);

        simulator.Solve();

//        //Check nothings moved
//        for (AbstractCellPopulation<2,3>::Iterator cell_iter = cell_population.Begin();
//            cell_iter != cell_population.End();
//            ++cell_iter)
//        {
//            c_vector<double,3> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
//            cell_location(2)=0.0;
//            TS_ASSERT_DELTA(norm_2(cell_location), 1.5, 1e-5);
//
//        }
    }
};

#endif /*TESTANGLEFORCE_HPP_*/


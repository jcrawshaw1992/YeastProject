#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <cxxtest/TestSuite.h>
#include "MembraneSurfaceForce.hpp"
#include "UblasCustomFunctions.hpp"

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
#include "MembraneDeformationForce.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"

#include "HemeLBForce.hpp"
#include "RadialForceOnCylinder.hpp"

#include "ElementAnglesWriter.hpp"
#include "RandomNumberGenerator.hpp"

#include "Honeycomb3DMeshGeneratorBentRectangle.hpp"

#include "MembraneBendingForce0TargetAngle.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:
    void TestBending() throw(Exception)
    {
        // TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-BendingParameter"));
        double BendingParameter = 8; //CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-BendingParameter");
        double dt = 0.02; //For most using 0.001, but for apsect ratio 3 and refinemnt 30 need finer

        double EndTime = 300;

        double SamplingTimestepMultiple = 1000; //2000;

        int N_C[3] = { 5, 10, 20 };
        double AspectRatio[3] = { 3, 0.75, 1.5 };

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {

                std::stringstream out;
                out << "AspectRatio" << AspectRatio[j] << "/Refinement" << N_C[i];
                std::string ParameterSet = out.str();
                std::string output_dir = "BendingForceOnBentRectanlge/SweepGrowth/" + ParameterSet;

                int n_z = N_C[i] * AspectRatio[j];
                Honeycomb3DMeshGeneratorBentRectangle generator(N_C[i], N_C[i] * AspectRatio[j], 1.5e-3, 1e-3);
                MutableMesh<2, 3>* p_mesh = generator.GetMesh();

                // Create the cells
                MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
                std::vector<CellPtr> cells;
                CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

                // Create a cell population
                HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*p_mesh, cells);
                cell_population.SetChasteOutputDirectory(output_dir, 0);
                // cell_population.SetUpInitialConfig(0);
                cell_population.SetWriteVtkAsPoints(true);
                cell_population.SetOutputMeshInVtk(true);
                cell_population.AddPopulationWriter<ElementAnglesWriter>();

                // Set population to output all data to results files
                cell_population.AddCellWriter<CellProliferativeTypesWriter>();

                OffLatticeSimulation<2, 3> simulator(cell_population);
                simulator.SetOutputDirectory(output_dir);
                simulator.SetSamplingTimestepMultiple(SamplingTimestepMultiple);
                simulator.SetDt(dt);
                simulator.SetUpdateCellPopulationRule(false);
                simulator.SetEndTime(EndTime);

                boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
                std::map<double, c_vector<long double, 4> > GrowthMaps;
                GrowthMaps[1] = Create_c_vector(pow(10, -8),0,0, 0);
                //Strength , hetro, stepsize, setupsolve
                p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 0, 100, 1);
                simulator.AddSimulationModifier(p_Mesh_modifier);

                /*
                -----------------------------
                Bending forces
                ----------------------------
                */
                boost::shared_ptr<MembraneBendingForce0TargetAngle> p_membrane_force(new MembraneBendingForce0TargetAngle());
                p_membrane_force->SetMembraneStiffness(pow(10, -BendingParameter));
                simulator.AddForce(p_membrane_force);

                /*
                -----------------------------
                Membrane forces
                ----------------------------
                */
                boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
                simulator.AddForce(p_shear_force);

                simulator.Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }
        }
    }

    void offTestBending() throw(Exception)
    {
        // TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-BendingParameter"));
        double BendingParameter = 8; //CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-BendingParameter");

        double dt = 0.0005; //For most using 0.001, but for apsect ratio 3 and refinemnt 30 need finer
        if (CommandLineArguments::Instance()->OptionExists("-dt"))
        {
            dt = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-dt");
        }
        double EndTime = 900;
        if (CommandLineArguments::Instance()->OptionExists("-EndTime"))
        {
            EndTime = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-EndTime");
        }
        double SamplingTimestepMultiple = 1000; //2000;
        if (CommandLineArguments::Instance()->OptionExists("-SamplingTimestepMultiple"))
        {
            SamplingTimestepMultiple = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-SamplingTimestepMultiple");
        }
        bool Archive = 0;
        if (Archive)
        {
            std::string ArchivedDirectory = "/BendingForceOnBentRectanlge/Bend_8";
            double EndTime1 = 700;
            // Load and fix any settings in the simulator
            OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(ArchivedDirectory, EndTime1);
            p_simulator->SetEndTime(EndTime1 + 500);
            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            // SimulationTime::Instance()->Destroy();
            // SimulationTime::Instance()->SetStartTime(0.0);
        }
        else
        {
            int N_C[3] = { 30, 20, 10 };
            double AspectRatio[3] = { 3, 0.75, 1.5 };

            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    if (AspectRatio[j] == 1.5)// & (N_C[i] == 10 | N_C[i] == 20))
                    {
                        continue;
                    }
                    else
                    {

                        std::stringstream out;
                        out << "AspectRatio" << AspectRatio[j] << "/Refinement" << N_C[i];
                        std::string ParameterSet = out.str();
                        std::string output_dir = "BendingForceOnBentRectanlge/SweepGrowth/" + ParameterSet;

                        int n_z = N_C[i] * AspectRatio[j];
                        Honeycomb3DMeshGeneratorBentRectangle generator(N_C[i], N_C[i] * AspectRatio[j], 1e-3, 1e-3);
                        MutableMesh<2, 3>* p_mesh = generator.GetMesh();

                        // Create the cells
                        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
                        std::vector<CellPtr> cells;
                        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
                        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

                        // Create a cell population
                        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*p_mesh, cells);
                        cell_population.SetChasteOutputDirectory(output_dir, 0);
                        // cell_population.SetUpInitialConfig(0);
                        cell_population.SetWriteVtkAsPoints(true);
                        cell_population.SetOutputMeshInVtk(true);
                        cell_population.AddPopulationWriter<ElementAnglesWriter>();

                        // Set population to output all data to results files
                        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

                        OffLatticeSimulation<2, 3> simulator(cell_population);
                        simulator.SetOutputDirectory(output_dir);
                        simulator.SetSamplingTimestepMultiple(SamplingTimestepMultiple);
                        simulator.SetDt(dt);
                        simulator.SetUpdateCellPopulationRule(false);
                        simulator.SetEndTime(EndTime);

                        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
                        std::map<double, c_vector<long double, 4> > GrowthMaps;
                        GrowthMaps[1] = Create_c_vector(pow(10, -7), pow(10, -9), pow(10, -10), 0);
                        //Strength , hetro, stepsize, setupsolve
                        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 0, 100, 1);
                        simulator.AddSimulationModifier(p_Mesh_modifier);

                        /*
                    -----------------------------
                    Bending forces
                    ----------------------------
                    */
                        boost::shared_ptr<MembraneBendingForce0TargetAngle> p_membrane_force(new MembraneBendingForce0TargetAngle());
                        p_membrane_force->SetMembraneStiffness(pow(10, -BendingParameter));
                        simulator.AddForce(p_membrane_force);

                        /*
                    -----------------------------
                    Membrane forces
                    ----------------------------
                    */
                        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
                        simulator.AddForce(p_shear_force);

                        simulator.Solve();
                        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

                        SimulationTime::Instance()->Destroy();
                        SimulationTime::Instance()->SetStartTime(0.0);
                    }
                }
            }
        }
    }
};

#endif /*TESTRELAXATION_HPP_*/

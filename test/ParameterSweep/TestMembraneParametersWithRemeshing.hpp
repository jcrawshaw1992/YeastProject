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
// #include "HemeLBForce.hpp"
#include "MembraneDeformationForce.hpp"
#include "OutwardsPressure.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:
    void TestParametersOverCylinder() throw(Exception)
    {

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-AreaParameter"));
        double AreaParameter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-AreaParameter");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-DilationParameter"));
        double DilationParameter =CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-DilationParameter");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-DeformationParamter"));
        // double DeformationParamter = 10;//CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-DeformationParamter");
        double DeformationParamter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-DeformationParamter");
        PRINT_3_VARIABLES(AreaParameter, DilationParameter, DeformationParamter)
        double dt= 0.01;
         if (CommandLineArguments::Instance()->OptionExists("-dt"))
        {
            dt= CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-dt");
        }
        double startime = 0;
        if (CommandLineArguments::Instance()->OptionExists("-startime"))
        {
            startime = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-startime");
        }
        double EndTime = 30;
        if (CommandLineArguments::Instance()->OptionExists("-EndTime"))
        {
            EndTime = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-EndTime");
        }
        double StartTime = 0;
        if (CommandLineArguments::Instance()->OptionExists("-StartTime"))
        {
            StartTime = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-StartTime");
        }
        double SamplingTimestepMultiple = 100;
        if (CommandLineArguments::Instance()->OptionExists("-SamplingTimestepMultiple"))
        {
            SamplingTimestepMultiple = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-SamplingTimestepMultiple");
        }
        double StartingParameterForSlowIncrease = 1e-8;
        if (CommandLineArguments::Instance()->OptionExists("-StartingParameterForSlowIncrease"))
        {
            StartingParameterForSlowIncrease= CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-StartingParameterForSlowIncrease");
        }

        std::stringstream out;
        out << "Param_" << AreaParameter << "_DilationParam_" << DilationParameter << "_DeformationParam_" << DeformationParamter;
        std::string ParameterSet = out.str();
        std::string output_dir = "ParameterSweepWithRemeshing/Cylinder/"+ParameterSet;


        std::string ArchivedDirectory = "ParameterSweep/Cylinder/";
        if (CommandLineArguments::Instance()->OptionExists("-ArchivedDirectory"))
        {
            ArchivedDirectory  = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-ArchivedDirectory");

            ////////

            // Load and fix any settings in the simulator
            OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(ArchivedDirectory, StartTime);

            /* Update the ouput directory for the population  */
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, StartTime);
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(StartTime);

            p_simulator->SetEndTime(EndTime + StartTime);
            p_simulator->SetOutputDirectory(output_dir + ParameterSet);
            /* 
            -----------------------------
            Update membrane properties
            ----------------------------
            */
            std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
            boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnHeteroMeshModifier<2, 3> >(*iter);
            std::map<double, c_vector<long double, 4> > GrowthMaps;
            GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
            //                                          Strength,hetro,stepsize, setupsolve
            p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 0, 100, 1);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }
        else
        {

            double scale = 1e3;
            double Length = 50e-6 * scale;
            double Radius = 1e-6 * scale; // I want this to grow to 10

            unsigned N_D = 20;//50;
            unsigned N_Z = N_D*2*5;//

            Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
            MutableMesh<2, 3>* p_mesh = generator.GetMesh();
            HistoryDepMutableMesh<2, 3>* mesh = static_cast<HistoryDepMutableMesh<2, 3>*>(p_mesh);


        // Create the cells 
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

            // Create a cell population
            HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
            cell_population.SetChasteOutputDirectory(output_dir, startime);
            cell_population.SetTargetRemeshingEdgeLength(0.3e-6 * scale); //0.2e-6 * scale); 
            cell_population.SetTargetRemeshingIterations(10);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);
            cell_population.SetRemeshingSoftwear("CGAL");
            cell_population.SetOperatingSystem("Server");
            // Set population to output all data to results files
            cell_population.AddCellWriter<CellProliferativeTypesWriter>();

            // Set up cell-based simulation
            OffLatticeSimulation<2,3> simulator(cell_population);
            simulator.SetOutputDirectory(output_dir);
            simulator.SetSamplingTimestepMultiple(SamplingTimestepMultiple);//100);
            simulator.SetDt(dt); // 0.008
            simulator.SetUpdateCellPopulationRule(false);
            simulator.SetEndTime(EndTime);




            /*
            -----------------------------
            RemeshingTriggerOnHeteroMeshModifier
            ----------------------------
            */
            boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
            std::map<double, c_vector<long double, 4> > GrowthMaps;
            GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
            p_Mesh_modifier->SetRemeshingInterval(1000);
            //Strength , hetro, stepsize, setupsolve
            p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 0, 100, 1);
            simulator.AddSimulationModifier(p_Mesh_modifier);


            /*
            -----------------------------
            Constant Compressive tissue pressure
            ----------------------------
            */

            double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
            double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

            boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
            p_ForceOut->SetPressure(P_blood - P_tissue);
            simulator.AddForce(p_ForceOut);

            /*
            -----------------------------
            Membrane forces
            ----------------------------
            */
            boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
            simulator.AddForce(p_shear_force);
            
            /*
            -----------------------------
            Boundary conditions
            ----------------------------
            */

            c_vector<double, 3> PlaneNormal1 = Create_c_vector(0, 0, 1);
            c_vector<double, 3> Point1 = Create_c_vector(0, 0, 1e-6 * scale);

            c_vector<double, 3> PlaneNormal2 = Create_c_vector(0, 0, -1);
            c_vector<double, 3> Point2 = Create_c_vector(0, 0, 49e-6 * scale);


            std::vector<c_vector<double, 3> > boundary_plane_points;
            std::vector<c_vector<double, 3> > boundary_plane_normals;

            boundary_plane_points.push_back(Point1);
            boundary_plane_normals.push_back(PlaneNormal1);

            boundary_plane_points.push_back(Point2);
            boundary_plane_normals.push_back(PlaneNormal2);


            for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
            {
                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&cell_population, boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 0.5));
                simulator.AddCellPopulationBoundaryCondition(p_condition);
            }


            simulator.Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }
    }
};

#endif /*TESTRELAXATION_HPP_*/

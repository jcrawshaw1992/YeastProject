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

// #include "AppliedForceModifier.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "EnclosedRegionBoundaryCondition.hpp"
#include "HemeLBForce.hpp"
#include "MembraneDeformationForce.hpp"
#include "MembraneDeformationForceOnCylinder.hpp"
#include "OutwardsPressureWithBreaks.hpp"
#include "OutwardsPressure.hpp"
#include "RemeshingTriggerOnStepHeteroModifier.hpp"
#include "StepHeteroModifier.hpp"
#include "MembraneBendingForce.hpp"

#include "MembraneBendingForceSensitive.hpp"


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
        cell_population.ExecuteHistoryDependentRemeshing();
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(SamplingStep);
        simulator.SetDt(dt);
        simulator.SetUpdateCellPopulationRule(false);
        /*
        -----------------------------
        StepHeteroModifier
        ----------------------------
        */
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnStepHeteroModifier<2, 3>());
        p_Mesh_modifier->SetMembraneStrength(1);
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        p_Mesh_modifier->SetCollapseType(1);
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); // I have turned this off because I need to know what will happen without remeshing, and then with remeshing
        simulator.AddSimulationModifier(p_Mesh_modifier);

        // /*
        // -----------------------------
        // Constant Compressive tissue pressure
        // ----------------------------
        // */
        // double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        // double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        // boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        // p_ForceOut->SetPressure(P_blood - P_tissue);
        // //simulator.AddForce(p_ForceOut);


        //   /*
        // -----------------------------
        // Membrane forces
        // ----------------------------
        // */
        // boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        // p_shear_force->SetCollapseType(1);
        // //simulator.AddForce(p_shear_force);

        // boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        // p_membrane_force->SetMembraneStiffness(pow(10, -7));
        // p_membrane_force->SetCollapseType(1);
        // //simulator.AddForce(p_membrane_force);


      
        // /*
        // -----------------------------
        // Boundary conditions
        // ----------------------------
        // */

        // std::vector<c_vector<double, 3> > boundary_plane_points;
        // std::vector<c_vector<double, 3> > boundary_plane_normals;


        // boundary_plane_points.push_back(Create_c_vector(0.016595194405737763,0.12835719623300196,-0.002236333573976076   ));
        // boundary_plane_normals.push_back(Create_c_vector( -0.9477460065412346,   0.31902584258272776,-0.000137293563571154   ));

        // boundary_plane_points.push_back(Create_c_vector(0.04573485665015318, 0.2339987318146692,  -0.0015555058516746553 ));
        // boundary_plane_normals.push_back(Create_c_vector(-0.810093202633425, 0.5861764198996681, -0.012091641771441898 ));


        // boundary_plane_points.push_back(Create_c_vector(0.1661163256187882,0.03058976097334303, -0.00045885113080672426 ));
        // boundary_plane_normals.push_back(Create_c_vector(-0.08886496326216727,-0.9954743652165818, -0.03367204331574329 ));


        // for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        // {
        //         boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&cell_population, boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 0.01));
        //         simulator.AddCellPopulationBoundaryCondition(p_condition);
        // }

        TRACE("First Solve ")

    
            // for (int i = 0; i < 5; i++)
            // {
                PRINT_VARIABLE(EndTime)
                cell_population.SetStartTime(EndTime);
                EndTime += 0.5;
                simulator.SetEndTime(EndTime);

                simulator.Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
            // }



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





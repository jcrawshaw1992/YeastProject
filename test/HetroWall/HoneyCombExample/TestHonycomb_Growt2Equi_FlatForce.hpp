#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"

#include "Debug.hpp"
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
#include "MembraneBendingForce.hpp"
#include "MembraneDeformationForce.hpp"
#include "StaticOutwardsPressure.hpp"
#include "OutwardsPressureWithBreaks.hpp"
#include "OutwardsPressure.hpp"
#include "MembraneBendingForce0TargetAngle.hpp"


#include "RemeshingTriggerOnStepHeteroModifier.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:

 void offTestMembraneParameters() throw(Exception)
   {


         TRACE("Jess is good")
        double EndTime = 14;
        double scale = 0.05;
        double SamplingStep = 25;
        double dt = 0.0001;
        double RemeshingTime = 5000000000000;
        double EdgeLength = 0.0004/2;//(2e-6 * scale);

        
        std::string output_dir = "DeformingHoneyComb/TestingMembraneParameters4/";
        std::string mesh_file = "/data/vascrem/MeshCollection/CouseHoneycomb.vtu";
        std::string Archieved =  "DeformingHoneyComb/FlatForce3";

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);

        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
        p_simulator->RemoveAllForces();
        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);



        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -9));
        p_simulator->AddForce(p_membrane_force);
        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(P_blood - P_tissue);
         
        p_simulator->AddForce(p_ForceOut);


        /*
        -----------------------------
        Update membrane properties
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime);

        // double DilationParameter=6;
        // double AreaParameter=7;//4.5;
        // double DeformationParamter=6;//5.3;



        // TestingMembraneParameters
        // std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -8.05), 1e-8) }  };


        // TestingMembraneParameters2
        // std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -7), pow(10, -9), pow(10, -8.05), 1e-8) } };


        // TestingMembraneParameters3
        // std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -7), pow(10, -8), pow(10, -8.05), 1e-8) }};

                // TestingMembraneParameters4
        std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -7), pow(10, -8), pow(10, -7), 1e-8) }};

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        p_Mesh_modifier->SetmSetUpSolve(1);
   


        for (int j = 0; j < 10; j++)
        {
            for (int i = 0; i <= 5; i++)
            {
                // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
                EndTime += 1;
                p_simulator->SetEndTime(EndTime);
                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }
            dt /= 2;
            SamplingStep *= 2;
            RemeshingTime *= 2;
            p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
            p_simulator->SetDt(dt);
        }

    }




   void offTestContinuoing() throw(Exception)
   {


         TRACE("Jess is good")
        double EndTime = 2;
        double scale = 0.05;
        double SamplingStep = 50;
        double dt = 0.001;
        double RemeshingTime = 500;
        double EdgeLength = 0.0004/2;//(2e-6 * scale);

        
        std::string output_dir = "DeformingHoneyComb/FlatForce2";
        std::string mesh_file = "/data/vascrem/MeshCollection/CouseHoneycomb.vtu";
        std::string Archieved =  "DeformingHoneyComb/FlatForce2";

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);

        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);

        // p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        // p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);

        /*
        -----------------------------
        Update membrane properties
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        
        // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 


        for (int j = 0; j < 10; j++)
        {
            for (int i = 0; i <= 50; i++)
            {
                // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
                EndTime += 1;
                p_simulator->SetEndTime(EndTime);
                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }
            dt /= 2;
            SamplingStep *= 2;
            RemeshingTime /= 2;
            p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
            p_simulator->SetDt(dt);
        }

    }


 


    void offTestNoRestraintFlat() throw(Exception)
    {
        TRACE("Jess is good")
        double EndTime = 0;
        double scale = 0.05;
        double SamplingStep = 50;
        double dt = 0.001;
        double RemeshingTime = 500;
        double EdgeLength = 0.0003/2;//(2e-6 * scale);

        
        std::string output_dir = "DeformingHoneyComb/FlatForce5";
        std::string mesh_file = "/data/vascrem/MeshCollection/CouseHoneycomb.vtu";

    
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        mesh.Scale(scale, scale, scale);

        // Create the cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);

        cell_population.SetChasteOutputDirectory(output_dir, 0);
        // cell_population.SetInitialAnlgesAcrossMembrane(); // Dont worry about this for now, I think there is something moff
        cell_population.SetRelativePath(output_dir, 0);
        cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
        cell_population.SetBinningIntervals(15, 5, 1);
        // cell_population.EdgeLengthVariable(1.2);
        cell_population.SetPrintRemeshedIC(1);
        cell_population.SetTargetRemeshingIterations(5);
        cell_population.SetWriteVtkAsPoints(false);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.SetRemeshingSoftwear("CGAL");
        cell_population.SetOperatingSystem("server");
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
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); // I have turned this off because I need to know what will happen without remeshing, and then with remeshing
       

        // std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -8.05), 1e-8) },
        //                                                             {0,  Create_c_vector(pow(10, -6), pow(10, -7), pow(10, -7), 1e-9)}    };

        // p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        simulator.AddSimulationModifier(p_Mesh_modifier);



        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(2 * (P_blood - P_tissue) / 10);
        simulator.AddForce(p_ForceOut);


        boost::shared_ptr<MembraneBendingForce0TargetAngle> p_membrane_force(new MembraneBendingForce0TargetAngle());
        p_membrane_force->SetMembraneStiffness(pow(10, -9));
        simulator.AddForce(p_membrane_force);
        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        // boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        // simulator.AddForce(p_shear_force);

        /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */

        std::vector<c_vector<double, 3> > boundary_plane_points;
        std::vector<c_vector<double, 3> > boundary_plane_normals;

        boundary_plane_points.push_back(Create_c_vector(0.007,0,0 ));
        boundary_plane_normals.push_back(Create_c_vector(-1, 0, 0));

        boundary_plane_points.push_back(Create_c_vector(0.074,0,0));
        boundary_plane_normals.push_back(Create_c_vector(1,0,0));

    
        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
   
                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&cell_population, boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id],2));
                simulator.AddCellPopulationBoundaryCondition(p_condition);

        }

        TRACE("First Solve ")

        // for (int j = 0; j < 10; j++)
        // {
            for (int i = 0; i < 25; i++)
            {
                PRINT_VARIABLE(EndTime)
                cell_population.SetStartTime(EndTime);
                EndTime +=1;
                simulator.SetEndTime(EndTime);

                simulator.Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
            }
            
            // dt /= 5; 
            //  SamplingStep *= 5; RemeshingTime /= 5; EdgeLength*=1.05;
            // simulator.SetSamplingTimestepMultiple(SamplingStep);
            // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime);
            // simulator.SetDt(dt);
            // cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
        // }
    }





    void TestNoRestraintFlat() throw(Exception)
    {
        TRACE("Jess is good")
        double EndTime = 0;
        double scale = 0.05;
        double SamplingStep = 50;
        double dt = 0.005;
        double RemeshingTime = 5000;
        double EdgeLength = 0.0004/2;//(2e-6 * scale);

        
        std::string output_dir = "DeformingHoneyComb/MinimialRemeshing";
        std::string mesh_file = "/data/vascrem/MeshCollection/CouseHoneycomb.vtu";

    
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        mesh.Scale(scale, scale, scale);

        // Create the cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);

        cell_population.SetChasteOutputDirectory(output_dir, 0);
        // cell_population.SetInitialAnlgesAcrossMembrane(); // Dont worry about this for now, I think there is something moff
        cell_population.SetRelativePath(output_dir, 0);
        cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
        cell_population.SetBinningIntervals(15, 5, 1);
        // cell_population.EdgeLengthVariable(1.2);
        cell_population.SetPrintRemeshedIC(1);
        cell_population.SetTargetRemeshingIterations(5);
        cell_population.SetWriteVtkAsPoints(false);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.SetRemeshingSoftwear("CGAL");
        cell_population.SetOperatingSystem("server");
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
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); // I have turned this off because I need to know what will happen without remeshing, and then with remeshing
        simulator.AddSimulationModifier(p_Mesh_modifier);



        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<StaticOutwardsPressure> p_ForceOut(new StaticOutwardsPressure());
        p_ForceOut->SetPressure(2 * (P_blood - P_tissue) / 10);
        p_ForceOut->SetUpNormals(cell_population);
        simulator.AddForce(p_ForceOut);


        boost::shared_ptr<MembraneBendingForce0TargetAngle> p_membrane_force(new MembraneBendingForce0TargetAngle());
        p_membrane_force->SetMembraneStiffness(pow(10, -7));
        simulator.AddForce(p_membrane_force);
        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        // boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        // simulator.AddForce(p_shear_force);

        /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */

        std::vector<c_vector<double, 3> > boundary_plane_points;
        std::vector<c_vector<double, 3> > boundary_plane_normals;

        boundary_plane_points.push_back(Create_c_vector(0.007,0,0 ));
        boundary_plane_normals.push_back(Create_c_vector(-1, 0, 0));

        boundary_plane_points.push_back(Create_c_vector(0.074,0,0));
        boundary_plane_normals.push_back(Create_c_vector(1,0,0));

    
        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
   
                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&cell_population, boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id],2));
                simulator.AddCellPopulationBoundaryCondition(p_condition);

        }

        TRACE("First Solve ")

   
            for (int i = 0; i < 25; i++)
            {
                PRINT_VARIABLE(EndTime)
                cell_population.SetStartTime(EndTime);
                EndTime +=1;
                simulator.SetEndTime(EndTime);

                simulator.Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
            }
            
    
    }


};

#endif /*TESTRELAXATION_HPP_*/

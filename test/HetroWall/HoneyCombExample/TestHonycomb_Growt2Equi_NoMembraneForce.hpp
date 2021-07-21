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
#include "OutwardsPressure.hpp"
#include "OutwardsPressureWithBreaks.hpp"


#include "RemeshingTriggerOnStepHeteroModifier.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:


    void offTestNoRestraint() throw(Exception)
    {
        TRACE("Jess is good")
        double EndTime = 0;
        double scale = 0.05;
        double SamplingStep = 50;
        double dt = 0.001;
        double RemeshingTime = 60;
        double EdgeLength = 0.0003/2;//(2e-6 * scale);

        
        std::string output_dir = "DeformingHoneyComb/NoMembraneForce";
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
        cell_population.SetBinningIntervals(10 ,3, 1);
        // cell_population.EdgeLengthVariable(1.2);
        cell_population.SetPrintRemeshedIC(1);
        cell_population.SetTargetRemeshingIterations(5);
        cell_population.SetWriteVtkAsPoints(true);
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
       

        std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -8.05), 1e-8) },
                                                                    {0,  Create_c_vector(pow(10, -6), pow(10, -7), pow(10, -7), 1e-9)}    };

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        simulator.AddSimulationModifier(p_Mesh_modifier);



        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(2 * (P_blood - P_tissue) / 3);
        simulator.AddForce(p_ForceOut);



        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        // p_membrane_force->SetMembraneStiffness(pow(10, -9));
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

        for (int j = 0; j < 10; j++)
        {
            for (int i = 0; i < 5; i++)
            {
                PRINT_VARIABLE(EndTime)
                cell_population.SetStartTime(EndTime);
                EndTime +=0.1;
                simulator.SetEndTime(EndTime);

                simulator.Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
            }
            
            dt /= 5;  SamplingStep *= 5; RemeshingTime /= 5; EdgeLength*=1.05;
            simulator.SetSamplingTimestepMultiple(SamplingStep);
            simulator.SetDt(dt);
            cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
        }
    }



   void TestContinuoing() throw(Exception)
   {

        TRACE("Jess is good")
        double EndTime = 0.1;
        double scale = 0.05;
        double SamplingStep = 100;
        double dt = 0.0001;
        double RemeshingTime = 600;
        double EdgeLength = 0.0003/2;//(2e-6 * scale);

        
        std::string output_dir = "DeformingHoneyComb/NoMembraneForce";
        std::string Archieved ="DeformingHoneyComb/NoMembraneForce";

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);

        /*
        -----------------------------
        Update membrane properties
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 

        std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -8.05), 1e-10) },
                                                                    {0,  Create_c_vector(pow(10, -6), pow(10, -7), pow(10, -7), 1e-9)}    };

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
      



        for (int j = 0; j < 10; j++)
        {
            for (int i = 0; i <= 5; i++)
            {
                static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
                EndTime += 0.1;
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


 
 
    void offTestSetUpCylinderArchive() throw(Exception)
    {
        TRACE("Jess is good")
        double EndTime = 0;
        double scale = 0.05;
        double SamplingStep = 50;
        double dt = 0.005;
        double RemeshingTime = 100;
        double EdgeLength = 0.0003/2;//(2e-6 * scale);

        
        std::string output_dir = "DeformingHoneyComb/Course";
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
        cell_population.SetBinningIntervals(10 ,3, 1);
        // cell_population.EdgeLengthVariable(1.2);
        cell_population.SetPrintRemeshedIC(1);
        cell_population.SetTargetRemeshingIterations(5);
        cell_population.SetWriteVtkAsPoints(true);
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
        p_Mesh_modifier->SetMembraneStrength(0.5);
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); // I have turned this off because I need to know what will happen without remeshing, and then with remeshing
        simulator.AddSimulationModifier(p_Mesh_modifier);

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(2 * (P_blood - P_tissue) / 3);
        simulator.AddForce(p_ForceOut);

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        // p_membrane_force->SetMembraneStiffness(pow(10, -9));
        simulator.AddForce(p_membrane_force);
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

        std::vector<c_vector<double, 3> > boundary_plane_points;
        std::vector<c_vector<double, 3> > boundary_plane_normals;

        boundary_plane_points.push_back(Create_c_vector(0.011839659156333001,0,0 ));
        boundary_plane_normals.push_back(Create_c_vector(-1, 0, 0));

        boundary_plane_points.push_back(Create_c_vector(0.145298053629932060,0,0));
        boundary_plane_normals.push_back(Create_c_vector(1,0,0));

    
        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
   
                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&cell_population, boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id],2));
                simulator.AddCellPopulationBoundaryCondition(p_condition);

        }

        TRACE("First Solve ")

        for (int j = 0; j < 10; j++)
        {
            for (int i = 0; i < 2; i++)
            {
                PRINT_VARIABLE(EndTime)
                cell_population.SetStartTime(EndTime);
                EndTime += 0.5;
                simulator.SetEndTime(EndTime);

                simulator.Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
            }
            
            dt /= 5;  SamplingStep *= 5; RemeshingTime /= 5; EdgeLength*=1.01;
            simulator.SetSamplingTimestepMultiple(SamplingStep);
            simulator.SetDt(dt);
            cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
        }
    }

    void offTestContinuingHomoArchieve() throw(Exception)
    {

        TRACE("Jess is good")
        double EndTime = 0.5;
        double scale = 0.05;
        double SamplingStep = 100;
        double dt = 0.0001;
        double RemeshingTime = 2000;
        double EdgeLength = 0.00035/2;//(2e-6 * scale);

        
        std::string output_dir = "DeformingHoneyComb/Course";
        std::string Archieved = "DeformingHoneyComb/Course";
  
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);


        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);
        // p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); // I have turned this off because I need to know what will happen without remeshing, and then with remeshing
        


        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);

        for (int j = 0; j < 50; j++)
        {
            for (int i = 0; i <= 10; i++)
            {
                static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
                EndTime += 0.2;
                p_simulator->SetEndTime(EndTime);
                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }
            dt /= 2;
            SamplingStep *= 2;
            RemeshingTime /= 2;
            p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            p_simulator->SetDt(dt);
            p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); // I have turned this off because I need to know what will happen without remeshing, and then with remeshing
        
        }
    }

    //  void offTestIntroduceHetro() throw(Exception)
    //     {
    //         std::string Archieved = "PlexusExample/";
    //         double EndTime = 15;
    //         OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved,EndTime);

    //         c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(-0.8572462913069383, 0.5148872691000456, -0.004460511091465417);
    //         c_vector<double, 3> UpperPlanePoint = Create_c_vector(0.04388350306511342, 0.036447307223086936, 6.2118799792568484e-6);
    //         c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0.8091715380078461,-0.5849424621690752, -0.05553141479916481);
    //         c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.02241255988723196, 0.051968449252480085,0.0014797089022245845);

    //         double dt= 0.0005;
    //         double NewEndTime = EndTime;

    //         double SamplingTimestepMultiple = 1000;

    //         std::string output_dir = "PlexusExample/HetroCollapse/";

    //         /* Update the ouput directory for the population  */
    //         static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);

    //         p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
    //         p_simulator->SetDt(dt);
    //         p_simulator->SetOutputDirectory(output_dir);

    //         /*
    //         -----------------------------
    //         Update membrane properties
    //         ----------------------------
    //         */
    //         std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
    //         boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<StepHeteroModifier<2, 3> >(*iter);
    //         p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );

    //         p_Mesh_modifier->SetUpdateFrequency(0.5/dt);

    //         for (int i =1; i<=5; i++)
    //         {
    //             static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
    //             NewEndTime +=4;
    //             p_simulator->SetEndTime(NewEndTime);

    //             p_simulator->Solve();
    //             CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
    //         }

    //         dt/=10;SamplingTimestepMultiple *= 10;

    //         p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
    //         p_simulator->SetDt(dt);
    //         p_Mesh_modifier->SetUpdateFrequency(0.5/dt);

    //         for (int i =1; i<=40; i++)
    //         {
    //             static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
    //             NewEndTime +=4;
    //             p_simulator->SetEndTime(NewEndTime);

    //             p_simulator->Solve();
    //             CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
    //         }

    //     }

    // void offTestIntroduceHetro() throw(Exception)
    // {
    //     double a =8;
    //     double b =2;
    //     double p =11;
    //     double dt= 0.0001;
    //     std::string output_dir = "HetroPlexus/";
    //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(output_dir, 15);
    //     // Load and fix any settings in the simulator

    //     double scale = 1e3;
    //     c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(-0.8572462913069383, 0.5148872691000456, -0.004460511091465417);
    //     c_vector<double, 3> UpperPlanePoint = Create_c_vector(0.04388350306511342, 0.036447307223086936, 6.2118799792568484e-6);
    //     c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0.8091715380078461,-0.5849424621690752, -0.05553141479916481);
    //     c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.02241255988723196, 0.051968449252480085,0.0014797089022245845);

    //     double DilationParameter=6;
    //     double AreaParameter=7;//4.5;
    //     double DeformationParamter=6;//5.3;

    //     double NewEndTime = 100;
    //     double EndTime = 15;

    //     double SamplingTimestepMultiple = 10000;

    //     /* Update the ouput directory for the population  */
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

    //     p_simulator->SetEndTime(EndTime + NewEndTime);
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
    //     GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), pow(10, -8));
    //     GrowthMaps[0] = Create_c_vector(pow(10, -5), pow(10, -6), pow(10, -4), pow(10, -6));
    //     // Strength,hetro,stepsize, setupsolve
    //     // GrowthMaps, Strength, Hetrogeneous,  StepSize,SetupSolve
    //     p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 1, pow(10,-p), 1);
    //     p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
    //     p_Mesh_modifier->SetBasementMembraneStrength(0);
    //     p_Mesh_modifier->SetPlateauParameters(a, b);
    //     p_Mesh_modifier->SetUpdateFrequency(100);

    //     p_simulator->Solve();
    //     CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);

    // }
};

#endif /*TESTRELAXATION_HPP_*/
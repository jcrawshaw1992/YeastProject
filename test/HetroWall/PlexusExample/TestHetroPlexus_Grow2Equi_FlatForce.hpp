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

#include "RemeshingTriggerOnStepHeteroModifier.hpp"

#include "MembraneBendingForceSensitive.hpp"
class TestRemeshing : public AbstractCellBasedTestSuite
{
public:




    void TestSetUpCylinderArchive2() throw(Exception)
    {
        TRACE("Jess is good")
        double EndTime = 0;
        double scale = 0.00006684491/1.29;

        double SamplingStep = 50;
        double dt = 0.006;
        double RemeshingTime = 250;
        double EdgeLength =0.00045;
        

        std::string output_dir = "DeformingPlexus/FlatForceFINAL8/";
        std::string mesh_file = "/data/vascrem/MeshCollection/Plexus_LongerInlets.vtu";
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
        cell_population.SetBinningIntervals(10, 10, 1);
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

        boost::shared_ptr<MembraneBendingForceSensitive> p_membrane_force(new MembraneBendingForceSensitive());
        p_membrane_force->SetMembraneStiffness(pow(10, -7));
        simulator.AddForce(p_membrane_force);
      
        /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */

        std::vector<c_vector<double, 3> > boundary_plane_points;
        std::vector<c_vector<double, 3> > boundary_plane_normals;




        boundary_plane_points.push_back(Create_c_vector( 0.023667025594574218, 0.03580448541106122, 1.0883377925970793e-6));
        boundary_plane_normals.push_back(-Create_c_vector( 0.9992096159136262, -0.037477299160519,-0.01325124572923273 ));

        boundary_plane_points.push_back(Create_c_vector( 0.02842355040448236, 0.04796528990124343, -0.0012984792809731791 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.6540763785122455,  0.7518206257291297, 0.08336568718943009 ));

        // boundary_plane_points.push_back(Create_c_vector(0.051817210775777105,0.05387570153154082, -0.00019851668457055458 ));
        // boundary_plane_normals.push_back(Create_c_vector(0.8769643127867986, 0.4790312122440197,  0.038245153868211076 ));

        boundary_plane_points.push_back(Create_c_vector(0.06117298902054698,0.03654224598923696,0.0009160147066410285 ));
        boundary_plane_normals.push_back(Create_c_vector(0.7792901394890018, 0.6252088966438896,  -0.04266982601959587  ));

        boundary_plane_points.push_back(Create_c_vector(0.055625630272326206, 0.0200783962665125, -0.000722688871900608 ));
        boundary_plane_normals.push_back(Create_c_vector( 0.7212385318941745, -0.6926132859491577, 0.010090403254884016  ));


        // boundary_plane_points.push_back(Create_c_vector(0.05266362792368705, 0.011971267579691799,-0.0014918860363199145   ));
        // boundary_plane_normals.push_back(Create_c_vector( -0.1352823453802107, -0.9855392137170218, 0.10203501974549702  ));


        boundary_plane_points.push_back(Create_c_vector( 0.026287008430771083, 0.023675441622844417, -0.0007261644860741899 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.6727411476463672, -0.7373576645038901,  0.061016578573512954 ));


        boundary_plane_points.push_back(Create_c_vector( 0.03952667347394293,0.01380981593118016, 0.00035914153313716104  ));
        boundary_plane_normals.push_back(Create_c_vector( -0.14188750875173126, -0.9842863935116474,0.10511056278063648));


        boundary_plane_points.push_back(Create_c_vector( 0.051916367697719554, 0.05396743908086633, 0.0004919967151319857 ));
        boundary_plane_normals.push_back(Create_c_vector(0.8831316656538849, 0.4690231691065414, 0.009784066672481936 ));


        boundary_plane_points.push_back(Create_c_vector(0.036804371668431334, 0.053919748549890005, -0.0007634162035095087 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.6416261436122902, 0.7612012690416773,  0.09427894697418274   ));
  
        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&cell_population, boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 0.01));
                simulator.AddCellPopulationBoundaryCondition(p_condition);
        }

        TRACE("First Solve ")

        // for (int j = 0; j < 10; j++)
        // {
            for (int i = 0; i < 5; i++)
            {
                PRINT_VARIABLE(EndTime)
                // cell_population.SetStartTime(EndTime);
                EndTime += 1;
                simulator.SetEndTime(EndTime);

                simulator.Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
            }
            

             PRINT_VARIABLE(EndTime)
            // cell_population.SetStartTime(EndTime);
            EndTime += 1;
            simulator.SetEndTime(EndTime);

            simulator.Solve();
            p_Mesh_modifier->TurnOffRemeshing(); 
            simulator.RemoveAllForces();
            simulator.RemoveAllCellPopulationBoundaryConditions();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

            // // dt /= 2;  SamplingStep *= 5; RemeshingTime /= 5; 
            // EdgeLength*=1.1; RemeshingTime*=2;
            // // simulator.SetSamplingTimestepMultiple(SamplingStep);
            // // simulator.SetDt(dt);
            // cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
            // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime);
        // }
    }



   void offTestContinuingHomoArchieve333() throw(Exception)
   {

        TRACE("Jess is good")
        double EndTime = 3;
        double scale = 0.00006684491 / 1.29;

        double SamplingStep = 50;
        double dt = 0.001;
        double RemeshingTime = 600;
        double EdgeLength =0.00045;
        

        std::string Archieved = "DeformingPlexus/FlatForceFINAL5/";
        std::string output_dir = "DeformingPlexus/FlatForceFINAL6/";
     
    
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetPrintRemeshedIC(1);
       


        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).ExecuteHistoryDependentRemeshing();

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllCellPopulationBoundaryConditions();
        

          /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */

        std::vector<c_vector<double, 3> > boundary_plane_points;
        std::vector<c_vector<double, 3> > boundary_plane_normals;



        boundary_plane_points.push_back(Create_c_vector( 0.023667025594574218, 0.03580448541106122, 1.0883377925970793e-6));
        boundary_plane_normals.push_back(-Create_c_vector( 0.9992096159136262, -0.037477299160519,-0.01325124572923273 ));

        boundary_plane_points.push_back(Create_c_vector( 0.02842355040448236, 0.04796528990124343, -0.0012984792809731791 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.6540763785122455,  0.7518206257291297, 0.08336568718943009 ));

        boundary_plane_points.push_back(Create_c_vector(0.051817210775777105,0.05387570153154082, -0.00019851668457055458 ));
        boundary_plane_normals.push_back(Create_c_vector(0.8769643127867986, 0.4790312122440197,  0.038245153868211076 ));

        boundary_plane_points.push_back(Create_c_vector(0.06117298902054698,0.03654224598923696,0.0009160147066410285 ));
        boundary_plane_normals.push_back(Create_c_vector(0.7792901394890018, 0.6252088966438896,  -0.04266982601959587  ));

        boundary_plane_points.push_back(Create_c_vector(0.055625630272326206, 0.0200783962665125, -0.000722688871900608 ));
        boundary_plane_normals.push_back(Create_c_vector( 0.7212385318941745, -0.6926132859491577, 0.010090403254884016  ));

        boundary_plane_points.push_back(Create_c_vector(0.05266362792368705, 0.011971267579691799,-0.0014918860363199145   ));
        boundary_plane_normals.push_back(Create_c_vector( -0.1352823453802107, -0.9855392137170218, 0.10203501974549702  ));


        boundary_plane_points.push_back(Create_c_vector( 0.026287008430771083, 0.023675441622844417, -0.0007261644860741899 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.6727411476463672, -0.7373576645038901,  0.061016578573512954 ));


        boundary_plane_points.push_back(Create_c_vector( 0.039700467012298214, 0.015344441098080423, -0.0004565231516265349 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.1706692384045456, -0.9809531883324618,  0.09275156798022458   ));



        boundary_plane_points.push_back(Create_c_vector(0.036804371668431334, 0.053919748549890005, -0.0007634162035095087 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.6416261436122902, 0.7612012690416773,  0.09427894697418274   ));

  
     
        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
      
             boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 0.01));
             p_simulator->AddCellPopulationBoundaryCondition(p_condition);
        }

 
        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime);                                                              
 
            
        for (int i =1; i<=2; i++)
        { 
            EndTime +=1;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

            // dt/=2 ;  SamplingStep*= 2; 
            // // FSIIterations*=2;
            // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
 
            // // p_ForceOut->SetFluidSolidIterations(FSIIterations);
            // p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            // p_simulator->SetDt(dt);
            // p_Mesh_modifier->SetUpdateFrequency(2/dt);
    }






   void offTestContinuingHomoArchieveHetro() throw(Exception)
   {

        TRACE("Jess is good")
        double EndTime =12;
        double SamplingStep = 100;
        double dt = 0.001;
        double RemeshingTime = 1000;
        double EdgeLength =0.0004;

        std::string output_dir = "DeformingPlexus/FlatForce5/";
        std::string Archieved = "DeformingPlexus/FlatForce4/";

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
    
 
        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime);                                                              
 
            
        for (int i =1; i<=100; i++)
        { 
            EndTime +=2;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

            // dt/=2 ;  SamplingStep*= 2; 
            // // FSIIterations*=2;
            // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
 
            // // p_ForceOut->SetFluidSolidIterations(FSIIterations);
            // p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            // p_simulator->SetDt(dt);
            // p_Mesh_modifier->SetUpdateFrequency(2/dt);
    }






    void offTestContinuingHomoArchieveHetro2() throw(Exception)
   {
        std::string Archieved = "DeformingPlexus/BendingForceOnly/MoreCorse3/";
        std::string output_dir = "DeformingPlexus/BendingForceOnly/TestingMembraneForceWanted/";

        double SamplingStep = 250;
        double dt = 0.00001;
        double RemeshingTime = 500000000;
        double EdgeLength =0.0004; // 0.0003;


        double EndTime = 5.3;
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        // p_simulator->RemoveAllForces();

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        // double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        // double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        // boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        // p_ForceOut->SetPressure(2 * (P_blood - P_tissue) / 3);
        // p_simulator->AddForce(p_ForceOut);


        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);


        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_simulator->AddForce(p_membrane_force);
 

        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
         */


        // //1
        std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -7.6), pow(10, -8), pow(10, -6), 1e-10 ) },
                                                                     {0,  Create_c_vector(pow(10, -5), pow(10, -5), pow(10, -5), 1e-10 )}    };


        // 2
        // GrowthMaps =  { {1, Create_c_vector(pow(10, -7), pow(10, -8), pow(10, -6), 1e-10 ) },
        //                                                              {0,  Create_c_vector(pow(10, -5), pow(10, -5), pow(10, -5), 1e-10 )}    };


        // // // // 3
        // GrowthMaps =  { {1, Create_c_vector(pow(10, -7), pow(10, -8), pow(10, -7), 1e-10 ) },
        //                                                              {0,  Create_c_vector(pow(10, -5), pow(10, -5), pow(10, -5), 1e-10 )}    };



        // ---------------------------------------
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime);                                                              
 

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);

        // // p_Mesh_modifier->SetStepSize(pow(10, -8));

        // // Upstream 
        // c_vector<double, 3> UpperPlanePoint = Create_c_vector(0.039304480630177864,0.03322737156824725,0.00018359707230489612);
        // c_vector<double, 3> UpperPlaneNormal = Create_c_vector(-0.8483371554100375,-0.5128640847062302,0.13150856006073028);
        // // Down stream
        // c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.02,0.02,0);
        // c_vector<double, 3> LowerPlaneNormal = Create_c_vector(1,1,0);
        // p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        // p_Mesh_modifier->SetUpdateFrequency(2/dt);
        // p_Mesh_modifier->SetmSetUpSolve(1);

       
        // for (int j =1; j<=40; j++)
        // {

            
            for (int i =1; i<=10; i++)
            { 
                // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
                EndTime +=2;
                p_simulator->SetEndTime(EndTime);

                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }

            dt/=2 ;  SamplingStep*= 2; 
            // FSIIterations*=2;
 
            // p_ForceOut->SetFluidSolidIterations(FSIIterations);
            p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            p_simulator->SetDt(dt);
            p_Mesh_modifier->SetUpdateFrequency(2/dt);


        // }


        // // ---------------------------------------
        // std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        // boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 



        // for (int j = 0; j < 10; j++)
        // {
        //     for (int i = 0; i <= 2; i++)
        //     {
        //         static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
        //         EndTime += 0.5;
        //         p_simulator->SetEndTime(EndTime);
        //         p_simulator->Solve();
        //         CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        //     }
        //     // (p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(0.0003);
        //     // dt /= 2;
        //     // SamplingStep *= 2;
        //     // RemeshingTime /= 2;
        //     // p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        //     // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
        //     // p_simulator->SetDt(dt);
        // }
    }



  void offTestContinuingHomoArchieve() throw(Exception)
   {
        std::string Archieved = "DeformingPlexus/BendingForceOnly/MoreCorse/";
        std::string output_dir = "DeformingPlexus/BendingForceOnly/MoreCorse3/";

        double SamplingStep = 25;
        double dt = 0.0001*5;
        double RemeshingTime = 50;
        double EdgeLength =0.0004; // 0.0003;
 

        double EndTime = 4.3;
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        // p_simulator->RemoveAllForces();

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        // double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        // double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        // boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        // p_ForceOut->SetPressure(2 * (P_blood - P_tissue) / 3);
        // p_simulator->AddForce(p_ForceOut);

        

        /*
        -----------------------------
        Update membrane properties
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 



        for (int j = 0; j < 10; j++)
        {
            for (int i = 0; i <= 2; i++)
            {
                static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
                EndTime += 0.5;
                p_simulator->SetEndTime(EndTime);
                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }
            // (p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(0.0003);
            // dt /= 2;
            // SamplingStep *= 2;
            // RemeshingTime /= 2;
            // p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
            // p_simulator->SetDt(dt);
        }
    }



    void offTestContinuingHomoArchieve1() throw(Exception)
   {
        std::string Archieved = "DeformingPlexus/BendingForceOnly/";
        std::string output_dir = "DeformingPlexus/BendingForceOnly/MoreCorse/";

        double SamplingStep = 5000;
        double dt = 0.00001;
        double RemeshingTime = 15000;
        double EdgeLength =0.0004;


        double EndTime = 3.2;
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllForces();

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(2 * (P_blood - P_tissue) / 3);
        p_simulator->AddForce(p_ForceOut);





        /*
        -----------------------------
        Update membrane properties
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 



        for (int j = 0; j < 10; j++)
        {
            for (int i = 0; i <= 2; i++)
            {
                static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
                EndTime += 0.1;
                p_simulator->SetEndTime(EndTime);
                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }
            // (p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(0.0003);
            dt /= 2;
            SamplingStep *= 2;
            RemeshingTime /= 2;
            p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
            p_simulator->SetDt(dt);
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

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
#include "MeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"
#include "MutableMesh.hpp"


#include "FixedRegionBoundaryCondition.hpp"
#include "MembraneBendingForce.hpp"
#include "MembraneDeformationForce.hpp"
#include "OutwardsPressure.hpp"

#include "RemeshingTriggerOnStepHeteroModifier.hpp"
#include "StepHeteroModifier.hpp"

#include "MembraneBendingForce0TargetAngle.hpp"
// #include "MembraneBendingForceSensitive.hpp"
// #include "NewModifier2.hpp"

// #include "BoundariesModifier.hpp"


class TestRemeshing : public AbstractCellBasedTestSuite
{
public:

    void TestSetUpCylinderArchive2() throw(Exception)
    {
        // TRACE("Jess is good")
        // double EndTime = 0;
        // double scale = 0.00006684491/1.29;

        // double SamplingStep = 50;
        // double dt = 0.05;
        // double RemeshingTime = 250;
        // double EdgeLength =0.00045;
        

        std::string output_dir = "DeformingPlexus/FixingArchieve/";
        std::string mesh_file = "/data/vascrem/MeshCollection/Plexus_LongerInlets.vtu";
        // VtkMeshReader<2, 3> mesh_reader(mesh_file);
        // MutableMesh<2, 3> mesh;
        // mesh.ConstructFromMeshReader(mesh_reader);

        // mesh.Scale(scale, scale, scale);


        double EndTime = 1;
        double scale = 1e3;
        double Length = 40e-6 * scale;
        double Radius = 0.5e-6 * scale; // I want this to grow to 10

        unsigned N_D = 50;
        unsigned N_Z = 150;

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();
        // MutableMesh<2, 3>* mesh = static_cast<MutableMesh<2, 3>*>(p_mesh);

        HistoryDepMutableMesh<2, 3>* mesh = static_cast<HistoryDepMutableMesh<2, 3>*>(p_mesh);


        // Create the cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);


            // MeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);


        // cell_population.SetChasteOutputDirectory(output_dir, 0);
        // // cell_population.SetRelativePath(output_dir, 0);
        // cell_population.SetTargetRemeshingEdgeLength(0.1); // Good
        // cell_population.SetBinningIntervals(10, 10, 1); // Good
        // cell_population.EdgeLengthVariable(1.2);// Good
        // cell_population.SetPrintRemeshedIC(1);// Good
        // cell_population.SetTargetRemeshingIterations(5);// Good
        // cell_population.SetOutputMeshInVtk(true);
        std::string OperatingSystem = "mac";
        // cell_population.SetOperatingSystem(0);
        // Set population to output all data to results files

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetDt(0.05);

        simulator.SetUpdateCellPopulationRule(false);
        /*
        -----------------------------
        StepHeteroModifier
        ----------------------------
        */
        // boost::shared_ptr<NewModifier2<2,3> > p_Mesh_modifier(new NewModifier2<2,3>());
        
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnStepHeteroModifier<2, 3>());
        p_Mesh_modifier->SetMembraneStrength(0.5);
        p_Mesh_modifier->SetRemeshingInterval(10000); // I have turned this off because I need to know what will happen without remeshing, and then with remeshing
        simulator.AddSimulationModifier(p_Mesh_modifier);

        // /*
        // -----------------------------
        // Constant Compressive tissue pressure
        // ----------------------------
        // */
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(2 * (P_blood - P_tissue) / 3);
        simulator.AddForce(p_ForceOut);

        // boost::shared_ptr<MembraneBendingForceSensitive> p_membrane_force(new MembraneBendingForceSensitive());
        // p_membrane_force->SetMembraneStiffness(pow(10, -7));
        // simulator.AddForce(p_membrane_force);
      
        // /*
        // -----------------------------
        // Boundary conditions
        // ----------------------------
        // */

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
            // for (int i = 0; i < 3; i++)
            // {
                PRINT_VARIABLE(EndTime)
                // cell_population.SetStartTime(0);
                // EndTime += 0.1;
                simulator.SetEndTime(EndTime);

                simulator.Solve();
                TRACE("A")
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
// 
                TRACE("B")
            // }

    }




  void TestContinuingHomoArchieve2() throw(Exception)
   {

       TRACE("Jess")
        double DilationParameter = -7;
        double AreaParameter = -6;
        double DeformationParamter = -8;


        TRACE("Jess is good")
        double EndTime = 1;
        double scale = 0.00006684491 / 1.29;

        double SamplingStep = 10;
        double dt = 0.002;
        double RemeshingTime = 1000;
        double EdgeLength =0.00045;
        

        std::string Archieved =   "DeformingPlexus/FixingArchieve";
        // std::string Archieved = "DeformingPlexus/FlatForceFINAL5/";
     
        // std::string output_dir = "DeformingPlexus/Grow2Equi/";
        std::string output_dir = "DeformingPlexus/TestParameters3/";
        PRINT_VARIABLE(Archieved);
     
    
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
 

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
         static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetPrintRemeshedIC(1);

        // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).ExecuteHistoryDependentRemeshing();

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllForces();
        p_simulator->RemoveAllCellPopulationBoundaryConditions();


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
        p_ForceOut->SetPressure((P_blood - P_tissue));
        p_simulator->AddForce(p_ForceOut);


        // /* 
        // -----------------------------
        // Edit  RemeshingTriggerOnStepHeteroModifier
        // ----------------------------
        // */

        


        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
        p_Mesh_modifier->TurnOffRemeshing();   
        //AreaConstant           AreaDilationModulus        ShearModulus
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), 1e-11) } };

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        p_Mesh_modifier->SetmSetUpSolve(1);


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
      
             boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 0.01));
             p_simulator->AddCellPopulationBoundaryCondition(p_condition);
        }

            
        for (int i =1; i<=10; i++)
        { 
            DeformationParamter +=0.1;
            std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), 1e-11) } };

            p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
            p_Mesh_modifier->SetmSetUpSolve(1);

            EndTime +=0.1;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

            // dt/=2 ;  SamplingStep*= 2; 
            // FSIIterations*=2;
            p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
 
            // p_ForceOut->SetFluidSolidIterations(FSIIterations);
            p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            p_simulator->SetDt(dt);
            p_Mesh_modifier->SetUpdateFrequency(2/dt);
    }



};

#endif /*TESTRELAXATION_HPP_*/

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


  void offTestWithConstantForce() throw(Exception)
   {
        double DilationParameter = -5.5;/////
        double AreaParameter = -5;///
        double DeformationParamter = -5;//////
        double BendingParameter = -7;

        // This was the first Idea 
        // double DilationParameter = -6;   
        // double AreaParameter = -6.5;
        // double DeformationParamter = -6;
        // double BendingParameter = -6;


        //AreaConstant           AreaDilationModulus        ShearModulus
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), pow(10, BendingParameter)) },
                                                                    {0,  Create_c_vector(pow(10, -4), pow(10, -4), pow(10, -4),pow(10, BendingParameter))} };
 
 
        TRACE("Jess is good")
        double EndTime = 0;
        double FSI_Iterations = 1000;

        double SamplingStep = 20;
        double dt = 0.001;
        double RemeshingTime = 10000;
        double EdgeLength =0.00045;
        
        std::string output_dir = "SimpleHemeLBPlexus2/GrowingToEqui3/";
     
        std::string mesh_file = "/data/vascrem/testoutput/DeformingPlexus/FlatForceFINAL9/results_from_time_3/mesh_50.vtu";
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);


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
        cell_population.SetTargetRemeshingIterations(10);
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
        // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime);
        p_Mesh_modifier->TurnOffRemeshing();   

        p_Mesh_modifier->SetMembraneStrength(1);
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);

        p_Mesh_modifier->SetmSetUpSolve(1);
        simulator.AddSimulationModifier(p_Mesh_modifier);
       

        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        simulator.AddForce(p_shear_force);

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -7));
        simulator.AddForce(p_membrane_force);


        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood - P_tissue));
        simulator.AddForce(p_ForceOut);

        /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */

        std::vector<c_vector<double, 3> > boundary_plane_points;
        std::vector<c_vector<double, 3> > boundary_plane_normals;


        boundary_plane_points.push_back(Create_c_vector( 0.023667025594574218, 0.03580448541106122, 1.0883377925970793e-6));
        boundary_plane_normals.push_back(-Create_c_vector( 0.9992096159136262, -0.037477299160519,-0.01325124572923273 ));

        boundary_plane_points.push_back(Create_c_vector( 0.036159340056213136,0.0532778119509024,  -0.0011348361020223623 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.6367706683117718, 0.7670909181305777,   0.07806817084681165  ));/////////

        boundary_plane_points.push_back(Create_c_vector(0.06117298902054698,0.03654224598923696,0.0009160147066410285 ));
        boundary_plane_normals.push_back(Create_c_vector(0.7792901394890018, 0.6252088966438896,  -0.04266982601959587  ));

        // boundary_plane_points.push_back(Create_c_vector(0.055625630272326206, 0.0200783962665125, -0.000722688871900608 ));
        // boundary_plane_normals.push_back(Create_c_vector( 0.7212385318941745, -0.6926132859491577, 0.010090403254884016  ));

        boundary_plane_points.push_back(Create_c_vector( 0.026287008430771083, 0.023675441622844417, -0.0007261644860741899 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.6727411476463672, -0.7373576645038901,  0.061016578573512954 ));////

        // 0.03794044772664341, 0.05264219522570732, 7.980308388269382e-5
        // -0.681589647076609, 0.7304037641169716, 0.0441122926375771

        boundary_plane_points.push_back(Create_c_vector( 0.03952667347394293,0.01380981593118016, 0.00035914153313716104  ));
        boundary_plane_normals.push_back(Create_c_vector( -0.14188750875173126, -0.9842863935116474,0.10511056278063648));

        boundary_plane_points.push_back(Create_c_vector( 0.050685301163470184, 0.05368317770751255, -0.000394611429412961 ));
        boundary_plane_normals.push_back(Create_c_vector(0.8831071868527388, 0.4691180632417636, 0.007066772200845511 ));

        // boundary_plane_points.push_back(Create_c_vector(0.036804371668431334, 0.053919748549890005, -0.0007634162035095087 ));
        // boundary_plane_normals.push_back(Create_c_vector(-0.6416261436122902, 0.7612012690416773,  0.09427894697418274   ));////////

        // boundary_plane_points.push_back(Create_c_vector(0.03685664010826702 ,0.05300873980967785, -0.0006933284667226157  ));/////////
        // boundary_plane_normals.push_back(Create_c_vector(-0.6848297687821392, 0.7273063820096684, 0.04509561484011702));////////////////

        boundary_plane_points.push_back(Create_c_vector(0.026864281281697915, 0.037668141898196575, 0.00018974398398148207 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.9998551753153383, -0.015809159893423953, -0.006300703024145976 ));

        boundary_plane_points.push_back(Create_c_vector(0.04013693457760815, 0.015329153153048572, 0.0005807009986400527 )); ///Create_c_vector(0.040547082072423954, 0.018275901698485374, 0.00039994117540888903 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.13760377878446087, -0.9886290666284141, 0.060644609667043325));////

        boundary_plane_points.push_back(Create_c_vector( 0.029862726176558368, 0.022150802525023005,  0.0007074621356822784));
        boundary_plane_normals.push_back(Create_c_vector(-0.6983037235236903,  -0.7144025416789002 ,0.044731623664661664  ));

        boundary_plane_points.push_back(Create_c_vector(  0.05518984915249851, 0.02090387761251939, 0.0007698080819251387  ));
        boundary_plane_normals.push_back(Create_c_vector(0.7069806813154489,  -0.7062411383519247, 0.0374402289806223));

        boundary_plane_points.push_back(Create_c_vector(0.05715768054001756,  0.03813699788011809,  0.0006839495043087638));
        boundary_plane_normals.push_back(Create_c_vector(0.7640527331012904, 0.6448510930472797, -0.01976079037329169 ));


     
        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
             boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&(simulator.rGetCellPopulation()) , boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 0.008));
             simulator.AddCellPopulationBoundaryCondition(p_condition);
        }



         for (int i =1; i<=20; i++)
        { 
            
            EndTime +=0.5;
            simulator.SetEndTime(EndTime);

            simulator.Solve();
            
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

        }

    }




  void TestWithHemeLBForce() throw(Exception)
   {
        double AreaParameter = -5.5;  double DilationParameter = -6; double DeformationParamter = -6; double BendingParameter = -12;
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), pow(10, BendingParameter)) }, {0,  Create_c_vector(pow(10, -4), pow(10, -4), pow(10, -4),pow(10, BendingParameter))} };

        std::string Archieved = "SimpleHemeLBPlexus2/GrowingToEqui3/";//std::string mesh_file = "/data/vascrem/testoutput/DeformingPlexus/FlatForceFINAL9/results_from_time_3/mesh_50.vtu";
        std::string output_dir = "FSISimulations/Plexus/EquiWithHemeLB/";
        double EndTime = 10;
        double SamplingStep = 50;
        double dt = 0.001;
        double RemeshingTime = 10000;
        double EdgeLength =0.00045;
        double FSI_Iterations = 100000000;

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
 

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetPrintRemeshedIC(1);
        p_simulator->RemoveAllForces();

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllForces();
        // p_simulator->RemoveAllCellPopulationBoundaryConditions();


        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        p_Mesh_modifier->SetStepSize(pow(10, -8));

    
        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -7));
        p_simulator->AddForce(p_membrane_force);


        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */
        // Inlet1
       c_vector<double, 3> Point1 = Create_c_vector(0.02297977879048401, 0.037075481876763385,  0.0004635334404919405);
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0.999910929986149, 0.0020912833516396026,-0.013181753607835568  );
        // Inlet2
        c_vector<double, 3> Point2 = Create_c_vector(0.036837793060229,0.055384952301980456,-0.0007597519518940717)  ;
        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0.6732805342808156,-0.7380737547778258,  -0.04405059212657047);

        // Outlet1 This is the bug bar
        c_vector<double, 3> Point3 = Create_c_vector(0.051243051697026636,0.05386889597979771, 0.00016345376729906492 ) ;
        c_vector<double, 3> PlaneNormal3 = Create_c_vector( -0.8734008505817445, -0.4862907639633924, 0.026310588875685135  );

        // Outlet2
        c_vector<double, 3> Point4 = Create_c_vector( 0.05849786867183286,0.04003892834739773, -8.346812953121241e-5) ;
        c_vector<double, 3> PlaneNormal4 = Create_c_vector( -0.7771752850914773, -0.6286228836915505, 0.028841746519592842);

        // Inlet3
        c_vector<double, 3> Point5 = Create_c_vector(0.056296622912376706, 0.020105116941221777, 0.0002243854912119816) ;
        c_vector<double, 3> PlaneNormal5 = Create_c_vector(-0.7208815333497773, 0.6929575165100672, -0.011819271867331964);

        // Outlet3
        c_vector<double, 3> Point6 = Create_c_vector(0.039492636709086544, 0.013468141696478432, -0.00030703284641268327) ;
        c_vector<double, 3> PlaneNormal6 = Create_c_vector(0.14577215937619503, 0.9879552859046741, -0.0519117578571096);

         // Inlet4
        c_vector<double, 3> Point7 = Create_c_vector(0.027954386139121126, 0.02122180981896879, 0.0008357837352478671) ;
        c_vector<double, 3> PlaneNormal7 = Create_c_vector(0.6878807670924608, 0.7247343474980099, -0.03975142539484216);


        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg

        double InletPressure = P_blood; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = P_blood * (0.98);

        boost::shared_ptr<HemeLBForce<2, 3> > p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure*1.1, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, InletPressure*1.05, "Inlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal5, Point5, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal6, Point6, OutletPressure*0.98, "Outlet");
        p_ForceOut->Inlets(PlaneNormal7, Point7, OutletPressure*0.95, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(FSI_Iterations);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"HemeLBForce/", p_simulator->rGetCellPopulation());
        p_simulator->AddForce(p_ForceOut);


      for (int i =1; i<=50; i++)
        { 
    
            EndTime +=1;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

    }




};

#endif /*TESTRELAXATION_HPP_*/





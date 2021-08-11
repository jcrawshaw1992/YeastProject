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


  void TestWithConstantForce() throw(Exception)
   {
        double DilationParameter = -5.5;/////
        double AreaParameter = -5;///
        double DeformationParamter = -5;//////
        double BendingParameter = -8;

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
        double dt = 0.002;
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



};

#endif /*TESTRELAXATION_HPP_*/





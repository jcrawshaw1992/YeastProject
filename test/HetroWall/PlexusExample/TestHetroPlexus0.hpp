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
#include "MembraneDeformationForce.hpp"
#include "OutwardsPressure.hpp"
#include "MembraneBendingForce.hpp"
#include "StepHeteroModifier.hpp"


class TestRemeshing : public AbstractCellBasedTestSuite
{
public:


    void TestSetUpCylinderArchive() throw(Exception)
    {
        double EndTime = 15;
        double scale = 0.00006684491/1.29;

        std::string output_dir = "PlexusExample/";
        // std::string mesh_file = "/Users/jcrawshaw/Documents/Projects/Meshes/Plexus2.vtu";
         std::string mesh_file = "/data/vascrem/MeshCollection/Plexus.vtu";
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        mesh.Scale(scale,scale,scale);
        
        // Create the cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        // cell_population.SetChasteOutputDirectory(output_dir, 0);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(500);
        simulator.SetDt(0.0001);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(EndTime);

        double DilationParameter=6;
        double AreaParameter=7;//4.5;
        double DeformationParamter=6;//5.3;


        /*
        -----------------------------
        StepHeteroModifier
        ----------------------------
        */
        boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier(new StepHeteroModifier<2, 3>());
        p_Mesh_modifier->SetMembraneStrength(1);
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


        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -9));
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

        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;


        boundary_plane_normals.push_back(Create_c_vector( 0.6876729091145614, -0.723181705399979, -0.06414196009395387 ));
        boundary_plane_points.push_back(Create_c_vector( 0.0724675362704713, 0.031054003292400725, -0.0009112922663102466  )/1.29 );

        boundary_plane_points.push_back(Create_c_vector( 0.07302015950969878,  0.04929005203304943 ,  -0.0009321324379871255  )/1.29 );
        boundary_plane_normals.push_back(Create_c_vector(0.8099647858285579, 0.583291387999826, -0.061059007549540904 ));

        boundary_plane_points.push_back(Create_c_vector(0.06070135750737772, 0.06678712447635708 , -0.0012765098114253428  )/1.29 );
        boundary_plane_normals.push_back(Create_c_vector( 0.8198804013531591,0.5725347541133972, -0.0002877660418513367 ));

        boundary_plane_points.push_back(Create_c_vector( 0.0492631187864004, 0.06635229901428957,  -0.0006305583601383584 )/1.29 );
        boundary_plane_normals.push_back(Create_c_vector(-0.6515501955094116, 0.7583557214561671, 0.01946644462514594 ));

        boundary_plane_points.push_back(Create_c_vector(0.03616972612522409,  0.048916381517621316,0.001733557990624898)/1.29 );
        boundary_plane_normals.push_back(Create_c_vector(-0.9994345263050555, 0.033620433310641906, 0.000542303784124859 ));

        boundary_plane_points.push_back(Create_c_vector( 0.03945803742561154, 0.03159871732688262, -0.001004170736673922  )/1.29 );
        boundary_plane_normals.push_back(Create_c_vector(-0.7275651727083541, -0.6859435736276984, -0.011416359346371455 ));

        boundary_plane_points.push_back(Create_c_vector( 0.05271817828553492, 0.023925699842467642,  0.0007781234886797855 )/1.29 );
        boundary_plane_normals.push_back(Create_c_vector(-0.13175654960891212, -0.9844249178554848, 0.11639498589018506 ));

        // -------------------------------------------


        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&cell_population, boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 5));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }

        simulator.Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
}




 void offTestIntroduceHetro() throw(Exception)
    {
        std::string Archieved = "PlexusExample/";
        double EndTime = 15;
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved,EndTime);

        c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(-0.8572462913069383, 0.5148872691000456, -0.004460511091465417);
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0.04388350306511342, 0.036447307223086936, 6.2118799792568484e-6);
        c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0.8091715380078461,-0.5849424621690752, -0.05553141479916481);
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.02241255988723196, 0.051968449252480085,0.0014797089022245845);


        double dt= 0.0005;
        double NewEndTime = EndTime;
    
        double SamplingTimestepMultiple = 1000;
    
        std::string output_dir = "PlexusExample/HetroCollapse/";
        
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
    
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);

        /* 
        -----------------------------
        Update membrane properties
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<StepHeteroModifier<2, 3> >(*iter);
        p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
  
        p_Mesh_modifier->SetUpdateFrequency(0.5/dt);


        for (int i =1; i<=5; i++)
        { 
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
            NewEndTime +=4;
            p_simulator->SetEndTime(NewEndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

        dt/=10;SamplingTimestepMultiple *= 10;
    
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_Mesh_modifier->SetUpdateFrequency(0.5/dt);

        for (int i =1; i<=40; i++)
        { 
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
            NewEndTime +=4;
            p_simulator->SetEndTime(NewEndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }




    }




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

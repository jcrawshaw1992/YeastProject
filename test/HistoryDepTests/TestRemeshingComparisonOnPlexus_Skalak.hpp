#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

//  You might like to fix up the boundaries here

#include <cxxtest/TestSuite.h>
#include <cstdio>
#include <ctime>
#include <cmath>

#include "Debug.hpp"

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VtkMeshWriter.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "XmlTools.hpp"
#include "UblasCustomFunctions.hpp"
#include "VtkMeshReader.hpp"
#include "OutwardsPressure.hpp"
#include "MembraneDeformationForce.hpp"
#include "HistoryDepMutableMesh.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "ElementQualityOutputModifier.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"



using namespace xsd::cxx::tree;

class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:


void TestPlexusDeformationWithRemeshingw() throw (Exception)
  {

        std::string mesh_file ="/Users/jcrawshaw/docker-polnet-master/Plexus_Course_2/Plexus.vtu";
        // std::string mesh_file ="/Users/jcrawshaw/docker-polnet-master/Plexus_Course_2/Plexus.vtu";
        std::string output_dir = "ElementQualityOverPlexus/WithRemeshing14/";//"RemeshingComparison/WithoutRemeshing2";

        double mesh_scale = 1e-3; 
        double startime = 0;
       
        // Read in the plexus 
        VtkMeshReader<2,3> mesh_reader(mesh_file);
        HistoryDepMutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in mm
      // Create the cells 
      
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetChasteOutputDirectory(output_dir, startime);
        cell_population.SetInitialAnlgesAcrossMembrane();
        cell_population.SetRelativePath(output_dir, startime);
        cell_population.SetTargetRemeshingEdgeLength(1.25*1e-2);
        cell_population.SetPrintRemeshedIC(1);
        // cell_population.SetTargetRemeshingIterations(10);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.SetRemeshingSoftwear("CGAL");
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(1000);
        simulator.SetDt(0.2);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(5000);

        /*
          -----------------------------
          RemeshingTriggerOnHeteroMeshModifier
          ----------------------------
        */  
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        p_Mesh_modifier->SetRemeshingInterval(3000);//  p_Mesh_modifier->SetRemeshingInterval(400); //(1000);
        simulator.AddSimulationModifier(p_Mesh_modifier);

        boost::shared_ptr<ElementQualityOutputModifier<2, 3> > p_ElementWriter_modifier(new ElementQualityOutputModifier<2, 3>());
        p_ElementWriter_modifier->SetElementMetricsFileName("/Users/jcrawshaw/Documents/testoutput/"+output_dir);
        p_ElementWriter_modifier->SetWriteoutInterval(500);
        simulator.AddSimulationModifier(p_ElementWriter_modifier);
        /*
        -----------------------------
        Compressive tissue pressure
        ----------------------------
        */        
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood-P_tissue)/50);
        simulator.AddForce(p_ForceOut);
        /*
        -----------------------------
        SMembrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        simulator.AddForce(p_shear_force);

        // Create a plane boundary to represent the inlets/outlets and pass them to the simulation

        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        // TRACE("These BCs need fixing, they no longer cover the plexus -- but this is not the current problem ")
        boundary_plane_points.push_back(Create_c_vector(1.0980424930584722,0.74252,-0.005747));
        boundary_plane_normals.push_back(Create_c_vector(-0.80324620881846220,-0.5905862612435627,0.07748));

        boundary_plane_points.push_back(Create_c_vector(0.9045, 1.00159, 0.004239));
        boundary_plane_normals.push_back(Create_c_vector(-0.79020,-0.61169,-0.03758));

        boundary_plane_points.push_back(Create_c_vector(0.73507,1.0234,0.0068));
        boundary_plane_normals.push_back(Create_c_vector(0.608,-0.78,-0.1215));

        boundary_plane_points.push_back(Create_c_vector(0.5317, 0.74258,0.00973));
        boundary_plane_normals.push_back(Create_c_vector(0.99112,-0.01647,0.1319));

        boundary_plane_points.push_back(Create_c_vector(0.596,0.4537330, 0.000464));
        boundary_plane_normals.push_back(Create_c_vector(0.745, 0.6579,0.1032));

        boundary_plane_points.push_back(Create_c_vector(0.7891,0.3674,-0.0126));
        boundary_plane_normals.push_back(Create_c_vector(0.17634,0.9833, 0.0441));
       // ---
        boundary_plane_points.push_back(Create_c_vector(1.05058,0.42000745,0.42000745));
        boundary_plane_normals.push_back(Create_c_vector(-0.7370,0.665,-0.12232));

        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],-boundary_plane_normals[boundary_id],0.5));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }
     	simulator.Solve();
    }



void TestPlexusDeformationWithOUTRemeshing() throw (Exception)
    {

        
        std::string mesh_file ="/Users/jcrawshaw/docker-polnet-master/Plexus_Course_2/Plexus.vtu";
        std::string output_dir = "ElementQualityOverPlexus/WithOUTRemeshing14";//"RemeshingComparison/WithoutRemeshing2";

        double mesh_scale = 1e-3; 
        // Read in the plexus 
        VtkMeshReader<2,3> mesh_reader(mesh_file);
        HistoryDepMutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in mm
      // Create the cells 
      
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(1000);
        simulator.SetDt(0.2);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(5000);

        /*
          -----------------------------
          RemeshingTriggerOnHeteroMeshModifier
          ----------------------------
        */  
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        simulator.AddSimulationModifier(p_Mesh_modifier);

        boost::shared_ptr<ElementQualityOutputModifier<2, 3> > p_ElementWriter_modifier(new ElementQualityOutputModifier<2, 3>());
        p_ElementWriter_modifier->SetElementMetricsFileName("/Users/jcrawshaw/Documents/testoutput/"+output_dir);
        p_ElementWriter_modifier->SetWriteoutInterval(500);
        simulator.AddSimulationModifier(p_ElementWriter_modifier);
        /*
        -----------------------------
        Compressive tissue pressure
        ----------------------------
        */        
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood-P_tissue)/50);
        simulator.AddForce(p_ForceOut);
        /*
        -----------------------------
        SMembrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        simulator.AddForce(p_shear_force);

        // Create a plane boundary to represent the inlets/outlets and pass them to the simulation

        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        // TRACE("These BCs need fixing, they no longer cover the plexus -- but this is not the current problem ")
        boundary_plane_points.push_back(Create_c_vector(1.0980424930584722,0.74252,-0.005747));
        boundary_plane_normals.push_back(Create_c_vector(-0.80324620881846220,-0.5905862612435627,0.07748));

        boundary_plane_points.push_back(Create_c_vector(0.9045, 1.00159, 0.004239));
        boundary_plane_normals.push_back(Create_c_vector(-0.79020,-0.61169,-0.03758));

        boundary_plane_points.push_back(Create_c_vector(0.73507,1.0234,0.0068));
        boundary_plane_normals.push_back(Create_c_vector(0.608,-0.78,-0.1215));

        boundary_plane_points.push_back(Create_c_vector(0.5317, 0.74258,0.00973));
        boundary_plane_normals.push_back(Create_c_vector(0.99112,-0.01647,0.1319));

        boundary_plane_points.push_back(Create_c_vector(0.596,0.4537330, 0.000464));
        boundary_plane_normals.push_back(Create_c_vector(0.745, 0.6579,0.1032));

        boundary_plane_points.push_back(Create_c_vector(0.7891,0.3674,-0.0126));
        boundary_plane_normals.push_back(Create_c_vector(0.17634,0.9833, 0.0441));
       // ---
        boundary_plane_points.push_back(Create_c_vector(1.05058,0.42000745,0.42000745));
        boundary_plane_normals.push_back(Create_c_vector(-0.7370,0.665,-0.12232));

        boundary_plane_points.push_back(Create_c_vector(0.5480351480035154, 0.40068510040911803,0.004945005753509855));
        boundary_plane_normals.push_back(Create_c_vector(0.7050992912676046,0.7091050070150413,0.0022535483398657854));


        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],-boundary_plane_normals[boundary_id],0.5));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }
     	simulator.Solve();
    }




void offTestPlexusDeformationWithoutRemeshing() throw (Exception)
    {   

        std::string mesh_file ="/Users/jcrawshaw/docker-polnet-master/Plexus/Plexuscourser.vtu";//CGALEdgeLengthControlledPlexus.vtu";
        // std::string mesh_file ="/Users/jcrawshaw/docker-polnet-master/Plexus_Course_2/Plexus.vtu";
        std::string output_dir = "ElementQualityOverPlexus/WithOutRemeshing9/";//"RemeshingComparison/WithoutRemeshing2";

        double mesh_scale = 1e-3; 
        double startime = 0;
       
        // Read in the plexus 
        VtkMeshReader<2,3> mesh_reader(mesh_file);
        HistoryDepMutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in mm
      // Create the cells 
      
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetChasteOutputDirectory(output_dir, startime);
        cell_population.SetInitialAnlgesAcrossMembrane();
        cell_population.SetPrintRemeshedIC(1);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(500);
        simulator.SetDt(0.002);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(1000);
        /*
          -----------------------------
          RemeshingTriggerOnHeteroMeshModifier
          ----------------------------
        */  
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        p_Mesh_modifier->SetMembraneStrength(2.5);
        simulator.AddSimulationModifier(p_Mesh_modifier);

        boost::shared_ptr<ElementQualityOutputModifier<2, 3> > p_ElementWriter_modifier(new ElementQualityOutputModifier<2, 3>());
        p_ElementWriter_modifier->SetElementMetricsFileName("/Users/jcrawshaw/Documents/testoutput/"+output_dir);
        p_ElementWriter_modifier->SetWriteoutInterval(100);
        simulator.AddSimulationModifier(p_ElementWriter_modifier);
        /*
        -----------------------------
        Compressive tissue pressure
        ----------------------------
        */        
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood-P_tissue)/20);
        simulator.AddForce(p_ForceOut);
        // /*
        // -----------------------------
        // SMembrane forces
        // ----------------------------
        // */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        simulator.AddForce(p_shear_force);

        // Create a plane boundary to represent the inlets/outlets and pass them to the simulation

        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        boundary_plane_points.push_back(Create_c_vector(0.04716878340592724, 0.04711,7.20109 ));
        boundary_plane_normals.push_back(-Create_c_vector(0.10911,0.99049,0.083)); 

        boundary_plane_points.push_back(Create_c_vector(0.071635,0.05838,-0.000199));
        boundary_plane_normals.push_back(Create_c_vector(-0.7329,0.67619,0.067)); 

        boundary_plane_points.push_back(Create_c_vector(0.027251,0.05714,0.00068855));
        boundary_plane_normals.push_back(Create_c_vector(0.70831,0.70377,0.0546)); 

        boundary_plane_points.push_back(Create_c_vector(0.017286969824066786, 0.083836,-0.00061));
        boundary_plane_normals.push_back(Create_c_vector(0.9988,0.001155,0.04765)); 

        boundary_plane_points.push_back(Create_c_vector(0.0433,0.1085,-0.0012));
        boundary_plane_normals.push_back(Create_c_vector(0.637,-0.767,0.069)); 

        boundary_plane_points.push_back(Create_c_vector(0.0469,0.045277,0.000412 ));
        boundary_plane_normals.push_back(Create_c_vector(0.18656588,0.9819,0.0321989 )); 

        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],-boundary_plane_normals[boundary_id],1));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }

        boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, Create_c_vector(0.06139278954760502,0.10697694119242532,0.00032382960933312023),-Create_c_vector(-0.777146,-0.62919, -0.0126),1));
        simulator.AddCellPopulationBoundaryCondition(p_condition);

     	simulator.Solve();
    }




void offTestPlexusDeformationWithRemeshing() throw (Exception)
    {

        std::string mesh_file ="/Users/jcrawshaw/docker-polnet-master/Plexus/Plexus.vtu";//CGALEdgeLengthControlledPlexus.vtu";
        // std::string mesh_file ="/Users/jcrawshaw/docker-polnet-master/Plexus_Course_2/Plexus.vtu";
        std::string output_dir = "ElementQualityOverPlexus/WithRemeshing9/";//"RemeshingComparison/WithoutRemeshing2";

        double mesh_scale = 1e-3; 
        double startime = 0;
       
        // Read in the plexus 
        VtkMeshReader<2,3> mesh_reader(mesh_file);
        HistoryDepMutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in mm
      // Create the cells 
      
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetChasteOutputDirectory(output_dir, startime);
        cell_population.SetInitialAnlgesAcrossMembrane();
        cell_population.SetRelativePath(output_dir, startime);
        cell_population.SetTargetRemeshingEdgeLength(1.5*1e-3);
        cell_population.SetPrintRemeshedIC(1);
        // cell_population.SetTargetRemeshingIterations(10);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.SetRemeshingSoftwear("CGAL");
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(1000);
        simulator.SetDt(0.02);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(1000);

        /*
          -----------------------------
          RemeshingTriggerOnHeteroMeshModifier
          ----------------------------
        */  
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        p_Mesh_modifier->SetRemeshingInterval(500);//  p_Mesh_modifier->SetRemeshingInterval(400); //(1000);
        p_Mesh_modifier->SetMembraneStrength(2.5);
        simulator.AddSimulationModifier(p_Mesh_modifier);

        boost::shared_ptr<ElementQualityOutputModifier<2, 3> > p_ElementWriter_modifier(new ElementQualityOutputModifier<2, 3>());
        p_ElementWriter_modifier->SetElementMetricsFileName("/Users/jcrawshaw/Documents/testoutput/"+output_dir);
        p_ElementWriter_modifier->SetWriteoutInterval(100);
        simulator.AddSimulationModifier(p_ElementWriter_modifier);
        /*
        -----------------------------
        Compressive tissue pressure
        ----------------------------
        */        
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood-P_tissue)/10);
        simulator.AddForce(p_ForceOut);
        /*
        -----------------------------
        SMembrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        simulator.AddForce(p_shear_force);

        // Create a plane boundary to represent the inlets/outlets and pass them to the simulation

        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        boundary_plane_points.push_back(Create_c_vector(0.04716878340592724, 0.04711,7.20109 ));
        boundary_plane_normals.push_back(-Create_c_vector(0.10911,0.99049,0.083)); 

        boundary_plane_points.push_back(Create_c_vector(0.071635,0.05838,-0.000199));
        boundary_plane_normals.push_back(Create_c_vector(-0.7329,0.67619,0.067)); 

        boundary_plane_points.push_back(Create_c_vector(0.027251,0.05714,0.00068855));
        boundary_plane_normals.push_back(Create_c_vector(0.70831,0.70377,0.0546)); 

        boundary_plane_points.push_back(Create_c_vector(0.017286969824066786, 0.083836,-0.00061));
        boundary_plane_normals.push_back(Create_c_vector(0.9988,0.001155,0.04765)); 



        boundary_plane_points.push_back(Create_c_vector(0.0433,0.1085,-0.0012));
        boundary_plane_normals.push_back(Create_c_vector(0.637,-0.767,0.069)); 

        boundary_plane_points.push_back(Create_c_vector(0.0469,0.045277,0.000412 ));
        boundary_plane_normals.push_back(Create_c_vector(0.18656588,0.9819,0.0321989 )); 

  
      
        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],-boundary_plane_normals[boundary_id],1));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }

        boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, Create_c_vector(0.06139278954760502,0.10697694119242532,0.00032382960933312023),-Create_c_vector(-0.777146,-0.62919, -0.0126),1));
         simulator.AddCellPopulationBoundaryCondition(p_condition);
        
     	simulator.Solve();
    }






 };



//        boundary_plane_points.push_back(Create_c_vector(0.523,0.734,-2.25));
//         boundary_plane_normals.push_back(Create_c_vector(0.999,-0.03557,0.0024));

//         boundary_plane_points.push_back(Create_c_vector(0.5971,0.45457,0.00166));
//         boundary_plane_normals.push_back(Create_c_vector(0.75169,0.659,-0.0038));

//         boundary_plane_points.push_back(Create_c_vector(0.78357,0.348,-0.00654));
//         boundary_plane_normals.push_back(Create_c_vector(0.1423,0.9849,-0.098));

//         boundary_plane_points.push_back(Create_c_vector(1.0688,0.40478,0.0029 ));
//         boundary_plane_normals.push_back(Create_c_vector( -0.72739,0.6858,-0.0218 ));


//         boundary_plane_points.push_back(Create_c_vector(1.08,0.732,0.004458));
//         boundary_plane_normals.push_back(Create_c_vector(-0.8,-0.59,0.0647));

//         boundary_plane_points.push_back(Create_c_vector(0.90644,1.00186,0.0029));
//         boundary_plane_normals.push_back(Create_c_vector(-0.96445,-0.5025,0.013317));
// // ---
//         boundary_plane_points.push_back(Create_c_vector(0.729,1.016095,0.00905));
//         boundary_plane_normals.push_back(Create_c_vector(0.6199170377,-0.78,0.04815));



#endif /*TESTRELAXATION_HPP_*/



 
// BC for cylinder
        //              //Create a plane boundary to represent the inlet and pass them to the simulation
        // c_vector<long double, 3> Boundary1 = Create_c_vector(0, 0, 0.01);
        // c_vector<long double, 3> Normal1 = -Create_c_vector(0, 0, -1);

        // c_vector<long double, 3> Boundary2 = Create_c_vector(0, 0, 0);
        // c_vector<long double, 3> Normal2 = -Create_c_vector(0, 0, 1);
        
        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 1));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 1));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_2);





//  boundary_plane_points.push_back(Create_c_vector(0.7237464113571987, 1.0242235,0.0014931));
//         boundary_plane_normals.push_back(Create_c_vector(0.62231829,-0.782622,-0.0148935));

//         boundary_plane_points.push_back(Create_c_vector(0.5177698571211321,0.7320682626002549, 0.004367 ));
//         boundary_plane_normals.push_back(Create_c_vector(0.99709303,-0.0030147,-0.076134));


//         boundary_plane_points.push_back(Create_c_vector(0.5935784823895799,0.4543416444213988, -0.0071795883647058825));
//         boundary_plane_normals.push_back(Create_c_vector(0.7443554916120831, 0.6620903684029474, 0.08701290809405006));

//         boundary_plane_points.push_back(Create_c_vector(0.78238,0.3441546,-0.00927  ));
//         boundary_plane_normals.push_back(Create_c_vector(  0.146856101,0.9867,-0.0688));

//         boundary_plane_points.push_back(Create_c_vector( 1.0608065119668328, 0.3975627743768795, 0.023235010360754368 ));
//         boundary_plane_normals.push_back(Create_c_vector(-0.6988607767318602,0.715012888190088, -0.018713216394033155));
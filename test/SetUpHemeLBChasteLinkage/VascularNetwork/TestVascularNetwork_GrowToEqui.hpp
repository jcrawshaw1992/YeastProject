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

        double AreaParameter = -5;  double DilationParameter = -5.5; double DeformationParamter = -5; double BendingParameter = -5;
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), pow(10, BendingParameter)) }, {0,  Create_c_vector(pow(10, -4), pow(10, -4), pow(10, -4),pow(10, BendingParameter))} };

        TRACE("Jess is good")
        double EndTime = 0;
        double scale = 0.0011; 

        double SamplingStep = 200;
        double dt = 0.01;
        double RemeshingTime = 7000;
        double EdgeLength =1.2*scale;
        
        /////////////////////////////////////////////////////////////////////////////////////
        std::string output_dir =  "FSISimulations/VascularNetworkLargerBendingForce/GrowingToEqui/ConstantForceArchiving/";
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
        // cell_population.ExecuteHistoryDependentRemeshing();
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

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(2*(P_blood - P_tissue)/3);
        simulator.AddForce(p_ForceOut);


          /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_shear_force->SetCollapseType(1);
        simulator.AddForce(p_shear_force);

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -7));
        p_membrane_force->SetCollapseType(1);
        simulator.AddForce(p_membrane_force);


      
        /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */

        std::vector<c_vector<double, 3> > boundary_plane_points;
        std::vector<c_vector<double, 3> > boundary_plane_normals;


        boundary_plane_points.push_back(Create_c_vector(0.016595194405737763,0.12835719623300196,-0.002236333573976076   ));
        boundary_plane_normals.push_back(Create_c_vector( -0.9477460065412346,   0.31902584258272776,-0.000137293563571154   ));

        boundary_plane_points.push_back(Create_c_vector(0.04573485665015318, 0.2339987318146692,  -0.0015555058516746553 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.810093202633425, 0.5861764198996681, -0.012091641771441898 ));


        boundary_plane_points.push_back(Create_c_vector(0.1661163256187882,0.03058976097334303, -0.00045885113080672426 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.08886496326216727,-0.9954743652165818, -0.03367204331574329 ));


        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&cell_population, boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 0.01));
                simulator.AddCellPopulationBoundaryCondition(p_condition);
        }

        TRACE("First Solve ")

    
            for (int i = 0; i < 50; i++)
            {
                PRINT_VARIABLE(EndTime)
                //cell_population.SetStartTime(EndTime);
                EndTime += 20;
                simulator.SetEndTime(EndTime);

                simulator.Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
            }



    }

  void TestWithHemeLBForce() throw(Exception)
   {

        double AreaParameter = -5;  double DilationParameter = -5.5; double DeformationParamter = -5; double BendingParameter = -5;
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), pow(10, BendingParameter)) }, {0,  Create_c_vector(pow(10, -4), pow(10, -4), pow(10, -4),pow(10, BendingParameter))} };

        TRACE("Jess is good")
        double EndTime = 20;
        double scale = 0.0011; 

        double SamplingStep = 200;
        double dt = 0.01;
        double RemeshingTime = 10000;
        double FSI_Iterations = 10000;
        double EdgeLength =1.2*scale;
        
        /////////////////////////////////////////////////////////////////////////////////////
        std::string Archieved = "FSISimulations/VascularNetworkLargerBendingForce/GrowingToEqui/ConstantForceArchiving/";;//std::string mesh_file = "/data/vascrem/testoutput/DeformingPlexus/FlatForceFINAL9/results_from_time_3/mesh_50.vtu";
        std::string output_dir = "FSISimulations/VascularNetworkLargerBendingForce/HemeLBEqui/ConstantForceArchiving/";
        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
 

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetPrintRemeshedIC(1);

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
        // p_Mesh_modifier->SetStepSize(pow(10, -8));
        // p_Mesh_modifier->SetmSetUpSolve(1);
    
        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_shear_force->SetCollapseType(1);
        p_simulator->AddForce(p_shear_force);

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -7));
        p_membrane_force->SetCollapseType(1);
        p_simulator->AddForce(p_membrane_force);


        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */
       // Inlet1
       c_vector<double, 3> Point1 = Create_c_vector(0.016595194405737763,0.12835719623300196,-0.002236333573976076);
       c_vector<double, 3> PlaneNormal1 = Create_c_vector(0.9477460065412346, -0.31902584258272776,-0.000137293563571154  );
        

       c_vector<double, 3> Point2 = Create_c_vector(0.04165292612840986, 0.23533064715360683, -0.0007601220860857502);
       c_vector<double, 3> PlaneNormal2 = Create_c_vector(0.7969315092041309, -0.604064075870372, -0.002600361609387508 );
       

       c_vector<double, 3> Point3 = Create_c_vector(0.1661163256187882,0.027296779685866353,0.002303039834163201);
       c_vector<double, 3> PlaneNormal3 = Create_c_vector(0.0807491790098989,0.9963585613635368,  0.027371285808498066 );
       

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg

        double InletPressure = P_blood; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = P_blood * (0.9);

        boost::shared_ptr<HemeLBForce<2, 3> > p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet"); // Issues here 
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet"); //FIne
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");// Issues here 
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetCollapseType(1);
        p_ForceOut->SetFluidSolidIterations(FSI_Iterations);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"HemeLBForce/", p_simulator->rGetCellPopulation());
        p_simulator->AddForce(p_ForceOut);


      for (int i =1; i<=50; i++)
        { 
    
            EndTime +=5;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

    }




};

#endif /*TESTRELAXATION_HPP_*/





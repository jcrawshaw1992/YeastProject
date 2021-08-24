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
        p_Mesh_modifier->SetMembraneStrength(1);
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        p_Mesh_modifier->SetCollapseType(1);
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); // I have turned this off because I need to know what will happen without remeshing, and then with remeshing
        // simulator.AddSimulationModifier(p_Mesh_modifier);

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(P_blood - P_tissue);
        //simulator.AddForce(p_ForceOut);


          /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_shear_force->SetCollapseType(1);
        //simulator.AddForce(p_shear_force);

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -7));
        p_membrane_force->SetCollapseType(1);
        //simulator.AddForce(p_membrane_force);


      
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

    
            // for (int i = 0; i < 5; i++)
            // {
                PRINT_VARIABLE(EndTime)
                // cell_population.SetStartTime(EndTime);
                EndTime += 1;
                simulator.SetEndTime(EndTime);

                simulator.Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
            // }



    }



void TestContinuedWithConstantForce() throw(Exception)
   {

        std::string Archieved ="FSISimulations/VascularNetwork/GrowingToEqui/ConstantForceArchiving/";

        double EndTime = 1;

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);

        
        for (int i =0; i<=3; i++)
            { 
        
                EndTime +=10;
                p_simulator->SetEndTime(EndTime);

                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }

    }






  void offTestWithHemeLBForce() throw(Exception)
   {

        double AreaParameter = -5;  double DilationParameter = -5.5; double DeformationParamter = -5; double BendingParameter = -7;
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), pow(10, BendingParameter)) }, {0,  Create_c_vector(pow(10, -4), pow(10, -4), pow(10, -4),pow(10, BendingParameter))} };

        std::string Archieved = "SimpleHemeLBPlexus2/GrowingToEqui3/";//std::string mesh_file = "/data/vascrem/testoutput/DeformingPlexus/FlatForceFINAL9/results_from_time_3/mesh_50.vtu";
        std::string output_dir = "FSISimulations/Plexus/EquiWithHemeLB3/";
        double EndTime = 10;
        double SamplingStep = 50;
        double dt = 0.001;
        double RemeshingTime = 10000;
        double EdgeLength =0.00045;
        double FSI_Iterations = 100000000;


        

        // 0.00115


        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step



// -0.04

// -0.1993


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
        p_Mesh_modifier->SetStepSize(pow(10, -8));

        p_Mesh_modifier->SetmSetUpSolve(1);
    
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

       c_vector<double, 3> Point1 = Create_c_vector(0.023244704133811703, 0.037075481876763385, 0);
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0 );
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
        double OutletPressure = P_blood * (0.9);

        boost::shared_ptr<HemeLBForce<2, 3> > p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure*1.01, "Inlet"); // Issues here 
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure*1.005, "Inlet"); //FIne
        p_ForceOut->Inlets(PlaneNormal3, Point3, InletPressure, "Inlet");// Issues here 
        p_ForceOut->Inlets(PlaneNormal4, Point4, InletPressure, "Inlet");// Issues here 
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





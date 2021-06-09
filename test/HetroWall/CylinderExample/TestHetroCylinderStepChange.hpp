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
#include "HistoryDepMutableMesh.hpp"


#include "FixedRegionBoundaryCondition.hpp"
#include "MembraneDeformationForceOnCylinder.hpp"
#include "OutwardsPressure.hpp"
#include "OutwardsPressureWithBreaks.hpp"
// #include "StepHeteroModifier.hpp"
#include "StepHeteroModifier.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:
    void offTestSetUpCylinderArchive() throw(Exception)
    {
        double EndTime = 5;
        double scale = 1e3;
        double Length = 40e-6 * scale;
        // double Radius = 1e-6 * scale; // I want this to grow to 10
        double Radius = 0.5e-6 * scale; // I want this to grow to 10

        unsigned N_D = 50;
        unsigned N_Z = 150;

        std::string output_dir = "StepChangeHetroCylinder/";

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();
        HistoryDepMutableMesh<2, 3>* mesh = static_cast<HistoryDepMutableMesh<2, 3>*>(p_mesh);

        // Create the cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
        // cell_population.SetChasteOutputDirectory(output_dir, 0);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(500);
        simulator.SetDt(0.005);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(EndTime);

        /*
        -----------------------------
        StepHeteroModifier
        ----------------------------
        */
        boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier(new StepHeteroModifier<2, 3>());
        p_Mesh_modifier->SetMembraneStrength(1);
        simulator.AddSimulationModifier(p_Mesh_modifier);

        // /*
        // -----------------------------
        // HemeLB Force
        // ----------------------------
        // */
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        boost::shared_ptr<OutwardsPressureWithBreaks> p_ForceOut(new OutwardsPressureWithBreaks());
        p_ForceOut->SetPressure(P_blood - P_tissue);
        p_ForceOut->SetInitialPosition(cell_population, 100);
        p_ForceOut->SetRadiusThreshold(10);
        simulator.AddForce(p_ForceOut);

        // /*
        // -----------------------------
        // Boundary conditions
        // ----------------------------
        // */

        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0, 0, -1);
        c_vector<double, 3> Point1 = Create_c_vector(0, 0, 1e-6 * scale);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0, 0, 1);
        c_vector<double, 3> Point2 = Create_c_vector(0, 0, 39e-6 * scale);


        std::vector<c_vector<double, 3> > boundary_plane_points;
        std::vector<c_vector<double, 3> > boundary_plane_normals;

        boundary_plane_points.push_back(Point1);
        boundary_plane_normals.push_back(PlaneNormal1);

        boundary_plane_points.push_back(Point2);
        boundary_plane_normals.push_back(PlaneNormal2);

        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&cell_population, boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 5));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }

        simulator.Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

}

    void TestRunningArchieve2() throw(Exception)
    {
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load("StepChangeHetroCylinder", 5);
  

        double DilationParameter=8.3;
        double AreaParameter=7.3;
        double DeformationParamter=8.1;

        double dt= 0.01;
        double NewEndTime = 1;
        double EndTime = 5;
        
        double SamplingTimestepMultiple = 100;
        std::string output_dir = "StepChangeHetroCylinder/";
        
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        p_simulator->RemoveAllForces();
        p_simulator->SetEndTime(EndTime + NewEndTime);
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
    

        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForceOnCylinder> p_shear_force(new MembraneDeformationForceOnCylinder());
        p_simulator->AddForce(p_shear_force);


        /*
        -----------------------------
        Outwards force
        ----------------------------
        */


        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_force(new OutwardsPressure());
        p_force->SetPressure(P_blood - P_tissue);
        // p_simulator->AddForce(p_force);


        p_simulator->Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);


    }

    // void TestRunningArchieve3() throw(Exception)
    // {

    //     // std::string Archieved = ;
    //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load("StepChangeHetroCylinder", 7);
    //     // Load and fix any settings in the simulator

   
    //     double DilationParameter=8.3;
    //     double AreaParameter=7.3;
    //     double DeformationParamter=8.1;

     
    //     double NewEndTime = 1;
    //     double EndTime = 7;
        
    //     double SamplingTimestepMultiple = 100;
    //     std::string output_dir = "StepChangeHetroCylinder/";
        
    //     /* Update the ouput directory for the population  */
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

    //     p_simulator->SetEndTime(EndTime + NewEndTime);
    //     p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
    //     p_simulator->SetOutputDirectory(output_dir);
    
    //     /* 
    //     -----------------------------
    //     Update membrane properties
    //     ----------------------------
    //     */
    //     std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
    //     boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<StepHeteroModifier<2, 3> >(*iter);
    //     std::map<double, c_vector<long double, 4> > GrowthMaps;
    //     GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
    //     GrowthMaps[0] = Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -6), 0);
    //     // Strength,hetro,stepsize, setupsolve
    //     // GrowthMaps, Strength, Hetrogeneous,  StepSize,SetupSolve
    //     p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 0, 1e-18, 1);
    //     p_simulator->Solve();
    //     CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);


    // }


    //  void TestIntroduceHetro1() throw(Exception)
    // {
    //     std::string Archieved = "StepChangeHetroCylinder/";
    //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, 8);

    //     c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(0, 0, 1);
    //     // c_vector<double, 3> UpperPlanePoint = Create_c_vector(12e-6 * scale, 0, 0);
    //     c_vector<double, 3> UpperPlanePoint = Create_c_vector(0, 0, 15e-6 * 1e3);
    //     c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0, 0, -1);
    //     c_vector<double, 3> LowerPlanePoint = Create_c_vector(0, 0,25e-6 * 1e3);

    //     double DilationParameter=8.4;
    //     double AreaParameter=7;
    //     double DeformationParamter=8;

     
    //     double dt= 0.001;
    //     double NewEndTime = 100;
    //     double EndTime = 8;
        
    //     double SamplingTimestepMultiple = 100;
    //     std::string output_dir = "StepChangeHetroCylinder/B_2/";
        
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
    //     boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<StepHeteroModifier<2, 3> >(*iter);
    //     std::map<double, c_vector<long double, 4> > GrowthMaps;
    //     GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
    //     GrowthMaps[0] = Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -5), 0);
    //     p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 1, 1e-10, 1);
    //     p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
    //     p_Mesh_modifier->SetBasementMembraneStrength(0);
    //     p_Mesh_modifier->SetPlateauParameters(8, 2);
    //     p_Mesh_modifier->SetUpdateFrequency(100);

    //     p_simulator->Solve();
    //     CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);

    // }


    //  void offTestIntroduceHetro2() throw(Exception)
    // {
    //     std::string Archieved = "StepChangeHetroCylinder/";
    //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, 8);

    //     c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(0, 0, 1);
    //     // c_vector<double, 3> UpperPlanePoint = Create_c_vector(12e-6 * scale, 0, 0);
    //     c_vector<double, 3> UpperPlanePoint = Create_c_vector(0, 0, 15e-6 * 1e3);
    //     c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0, 0, -1);
    //     c_vector<double, 3> LowerPlanePoint = Create_c_vector(0, 0,25e-6 * 1e3);

    //     double DilationParameter=8.4;
    //     double AreaParameter=7;
    //     double DeformationParamter=8;

     
    //     double dt= 0.001;
    //     double NewEndTime = 100;
    //     double EndTime = 8;
        
    //     double SamplingTimestepMultiple = 100;
    //     std::string output_dir = "StepChangeHetroCylinder/B_0.5/";
        
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
    //     boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<StepHeteroModifier<2, 3> >(*iter);
    //     std::map<double, c_vector<long double, 4> > GrowthMaps;
    //     GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
    //     GrowthMaps[0] = Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -5), 0);
    //     p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 1, 1e-10, 1);
    //     p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
    //     p_Mesh_modifier->SetBasementMembraneStrength(0);
    //     p_Mesh_modifier->SetPlateauParameters(8, 0.5);
    //     p_Mesh_modifier->SetUpdateFrequency(100);

    //     p_simulator->Solve();
    //     CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);

    // }



    //  void offTestIntroduceHetro3() throw(Exception)
    // {
    //     std::string Archieved = "StepChangeHetroCylinder/";
    //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, 8);

    //     c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(0, 0, 1);
    //     // c_vector<double, 3> UpperPlanePoint = Create_c_vector(12e-6 * scale, 0, 0);
    //     c_vector<double, 3> UpperPlanePoint = Create_c_vector(0, 0, 15e-6 * 1e3);
    //     c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0, 0, -1);
    //     c_vector<double, 3> LowerPlanePoint = Create_c_vector(0, 0,25e-6 * 1e3);

    //     double DilationParameter=8.4;
    //     double AreaParameter=7;
    //     double DeformationParamter=8;

     
    //     double dt= 0.001;
    //     double NewEndTime = 100;
    //     double EndTime = 8;
        
    //     double SamplingTimestepMultiple = 100;
    //     std::string output_dir = "StepChangeHetroCylinder/B_4/";
        
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
    //     boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<StepHeteroModifier<2, 3> >(*iter);
    //     std::map<double, c_vector<long double, 4> > GrowthMaps;
    //     GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
    //     GrowthMaps[0] = Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -5), 0);
    //     p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 1, 1e-10, 1);
    //     p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
    //     p_Mesh_modifier->SetBasementMembraneStrength(0);
    //     p_Mesh_modifier->SetPlateauParameters(8, 4);
    //     p_Mesh_modifier->SetUpdateFrequency(100);

    //     p_simulator->Solve();
    //     CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);

    // }




};

#endif /*TESTRELAXATION_HPP_*/

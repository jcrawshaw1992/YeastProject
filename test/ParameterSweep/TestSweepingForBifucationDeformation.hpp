#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
// #include <cstdio>
// #include <ctime>
// #include <cmath>
// #include <vector> 

#include "Debug.hpp"
#include "VtkMeshReader.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "SmartPointers.hpp"
 
#include "AbstractCellBasedTestSuite.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"

// #include "AppliedForceModifier.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "MembraneDeformationForce.hpp"
// #include "OutwardsPressureWithBreaks.hpp"
#include "OutwardsPressure.hpp"
#include "MembraneBendingForce.hpp"
#include "HemeLBForce.hpp"


class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:

  void TestParametersOnBifucation() throw (Exception)
    {
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-AreaParameter"));
        double AreaParameter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-AreaParameter");
        double DilationParameter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-DilationParameter");
        double ShearParameter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-ShearParameter");
        double BendingParameter =CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-BendingParameter");

        double dt= 0.000001;
         if (CommandLineArguments::Instance()->OptionExists("-dt"))
        {
            dt= CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-dt");
        }
        PRINT_VARIABLE(dt)
        double EndTime = 50;
        if (CommandLineArguments::Instance()->OptionExists("-EndTime"))
        {
            EndTime = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-EndTime");
        }
        double SamplingTimestepMultiple = 5000;
        if (CommandLineArguments::Instance()->OptionExists("-SamplingTimestepMultiple"))
        {
            SamplingTimestepMultiple = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-SamplingTimestepMultiple");
        }
        std::string mesh_file = "/Users/jcrawshaw/Documents/Projects/MeshCollection/SimpleBifucation/mesh.vtu";
        if (CommandLineArguments::Instance()->OptionExists("-MeshFile"))
        {
            mesh_file  = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-MeshFile");
        }

        //////////
        
        double scale = 1;//1e-3;  
        double scale2 = 1;
        // std::string output_dir = "FSIIdealNetwork/SimpleBifucation_NoRemeshing/";

        std::stringstream out;
        out << "Area_" << AreaParameter<< "_Dil_" << DilationParameter << "_Shear_" << ShearParameter << "_Bend_"<< BendingParameter;
        std::string ParameterSet = out.str();
        std::string output_dir = "BifucationSweep_Finner/"+ParameterSet;

        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(scale*scale2,scale*scale2,scale*scale2);

       // Create the cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetInitialAnlgesAcrossMembrane();
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Binning Functions
        // cell_population.SetBinningIntervals(2,2,1);
        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        simulator.SetDt(dt);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(EndTime);

        /*
        -----------------------------
        RemeshingTriggerOnHeteroMeshModifier
        ----------------------------
        */  

        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        p_Mesh_modifier->SetMembraneParameters(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -ShearParameter), pow(10, -BendingParameter));
        // p_Mesh_modifier->SetRemeshingInterval(200);//
        simulator.AddSimulationModifier(p_Mesh_modifier);

        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */        

        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.0026368426410883503,-0.0707,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(-0.7,-0.7,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.08361380518549465,-0.03684463931177293,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-0.7,0.7,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.08361380518549465,-0.10299085601436264,0);


        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood - P_tissue));// needs to be negative for server ?? 
        simulator.AddForce(p_ForceOut);

        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        simulator.AddForce(p_shear_force);


       /*
        -----------------------------
        Bending forces
        ----------------------------
        */
        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -BendingParameter));
        simulator.AddForce(p_membrane_force);

        /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */
        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        boundary_plane_points.push_back(Point1);
        boundary_plane_normals.push_back(PlaneNormal1);

        boundary_plane_points.push_back(Point2);
        boundary_plane_normals.push_back(PlaneNormal2);
        
        boundary_plane_points.push_back(Point3);
        boundary_plane_normals.push_back(PlaneNormal3);
        

        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
          boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],-boundary_plane_normals[boundary_id],0.1));
          simulator.AddCellPopulationBoundaryCondition(p_condition);
        }
      
     	  simulator.Solve();
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
}


};




#endif /*TESTRELAXATION_HPP_*/

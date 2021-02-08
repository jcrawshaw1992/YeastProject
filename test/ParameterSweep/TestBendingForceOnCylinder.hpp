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
// #include "RadialForceOnCylinderWithBreaks.hpp"
#include "RadialForceOnCylinder.hpp"
#include "MembraneBendingForce.hpp"
#include "HemeLBForce.hpp"


class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:

  void TestParametersOnBifucation() throw (Exception)
    {
        // TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-BendingParameter"));
        // double BendingParameter =CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-BendingParameter");

        double dt= 0.01;
         if (CommandLineArguments::Instance()->OptionExists("-dt"))
        {
            dt= CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-dt");
        }
        double EndTime = 600;
        if (CommandLineArguments::Instance()->OptionExists("-EndTime"))
        {
            EndTime = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-EndTime");
        }
        double SamplingTimestepMultiple = 100;//2000;
        if (CommandLineArguments::Instance()->OptionExists("-SamplingTimestepMultiple"))
        {
            SamplingTimestepMultiple = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-SamplingTimestepMultiple");
        }
        std::string mesh_file = "/Users/jcrawshaw/Documents/Projects/MeshCollection/Blobby.vtu";
        if (CommandLineArguments::Instance()->OptionExists("-MeshFile"))
        {
            mesh_file  = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-MeshFile");
        }
        double BendingParameter = 8;
        std::stringstream out;
        out << "_Bend_"<< BendingParameter;
        std::string ParameterSet = out.str();
        std::string output_dir = "SweepOnCylinder/"+ParameterSet;

        

        Honeycomb3DCylinderMeshGenerator generator(80, 80 *1.5, 1.5e-3, 12e-3);
        MutableMesh<2,3>* p_mesh = generator.GetMesh();
        // p_mesh->Translate(0,0,-6e-3);
        
        // Loop over nodes and translate
        

       // Create the cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*p_mesh, cells);
        cell_population.SetChasteOutputDirectory(output_dir, 0);
        cell_population.SetInitialAnlgesAcrossMembrane();
        cell_population.SetRelativePath(output_dir, 0);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.SetRemeshingSoftwear("CGAL");
        
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

   
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        simulator.SetDt(dt);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(EndTime);

        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */        

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */

        boost::shared_ptr<RadialForceOnCylinder> p_ForceOut(new RadialForceOnCylinder());
        p_ForceOut->SetPressure((P_blood - P_tissue)/10);// needs to be negative for server ?? 
        simulator.AddForce(p_ForceOut);



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

        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0,0,1);
        c_vector<double, 3> Point1 = Create_c_vector(0,0,0.00001);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0,-0,-1);
        c_vector<double, 3> Point2 = Create_c_vector(0,0,0.011);


        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        boundary_plane_points.push_back(Point1);
        boundary_plane_normals.push_back(PlaneNormal1);

        boundary_plane_points.push_back(Point2);
        boundary_plane_normals.push_back(PlaneNormal2);        

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


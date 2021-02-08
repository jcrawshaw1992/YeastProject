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
// #include "RadialForceWithBreaks.hpp"
#include "RadialForce.hpp"
#include "MembraneBendingForce.hpp"
#include "HemeLBForce.hpp"


class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:

  void TestParametersOnBifucation() throw (Exception)
    {
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-BendingParameter"));
        double BendingParameter =CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-BendingParameter");

        double dt= 0.005;
         if (CommandLineArguments::Instance()->OptionExists("-dt"))
        {
            dt= CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-dt");
        }
        double EndTime = 500;
        if (CommandLineArguments::Instance()->OptionExists("-EndTime"))
        {
            EndTime = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-EndTime");
        }
        double SamplingTimestepMultiple = 1000;
        if (CommandLineArguments::Instance()->OptionExists("-SamplingTimestepMultiple"))
        {
            SamplingTimestepMultiple = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-SamplingTimestepMultiple");
        }
        std::string mesh_file = "/Users/jcrawshaw/Documents/Projects/MeshCollection/Blobby.vtu";
        if (CommandLineArguments::Instance()->OptionExists("-MeshFile"))
        {
            mesh_file  = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-MeshFile");
        }

        std::stringstream out;
        out << "_Bend_"<< BendingParameter;
        std::string ParameterSet = out.str();
        std::string output_dir = "SweepOnBlobby/"+ParameterSet;

        // std::string output_dir = "SweepOnBlobby/";

        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        
        // Loop over nodes and translate
        double MeanX = 0;
        double MeanY = 0;
        double MeanZ = 0;
        for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
        {
                c_vector<double, 3> Location = mesh.GetNode(i)->rGetLocation();
                MeanX += Location[0];
                MeanY += Location[1];
                MeanZ += Location[2];
        }
        MeanX/=mesh.GetNumNodes();
        MeanY/=mesh.GetNumNodes();
        MeanZ/=mesh.GetNumNodes();
        mesh.Translate(MeanX,MeanY, MeanZ);

       // Create the cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
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

        boost::shared_ptr<RadialForce> p_ForceOut(new RadialForce());
        p_ForceOut->SetPressure((P_blood - P_tissue));// needs to be negative for server ?? 
        p_ForceOut->SetRadiusThreshold(0.3);
        simulator.AddForce(p_ForceOut);



       /*
        -----------------------------
        Bending forces
        ----------------------------
        */
        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -BendingParameter));
        simulator.AddForce(p_membrane_force);

    
     	  simulator.Solve();
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
}


};




#endif /*TESTRELAXATION_HPP_*/


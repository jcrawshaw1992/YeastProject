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

#include "RemeshingTriggerOnHeteroMeshModifier.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "MembraneDeformationForce.hpp"
#include "RadialForceOnCylinder.hpp"
#include "MembraneBendingForce.hpp"
#include "HemeLBForce.hpp"

#include "RandomNumberGenerator.hpp"
#include "ElementAnglesWriter.hpp"


class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:

  void TestParametersOnBifucation() throw (Exception)
    {
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-BendingParameter"));
        double BendingParameter =CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-BendingParameter");

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
        double SamplingTimestepMultiple = 500;//2000;
        if (CommandLineArguments::Instance()->OptionExists("-SamplingTimestepMultiple"))
        {
            SamplingTimestepMultiple = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-SamplingTimestepMultiple");
        }
        std::string mesh_file = "/Users/jcrawshaw/Documents/Projects/MeshCollection/Sphere.vtu";
        if (CommandLineArguments::Instance()->OptionExists("-MeshFile"))
        {
            mesh_file  = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-MeshFile");
        }
        std::stringstream out;
        out << "_Bend_"<< BendingParameter;
        std::string ParameterSet = out.str();
        std::string output_dir = "SweepOnSphere/"+ParameterSet;

        double Scale = 0.001;
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(Scale,Scale,Scale);

    
       // Create the cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetChasteOutputDirectory(output_dir, 0);
        cell_population.SetUpInitialConfig(0);
        cell_population.SetInitialAnlgesAcrossMembrane();
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.AddPopulationWriter<ElementAnglesWriter>();
        

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
        Bending forces
        ----------------------------
        */
        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -BendingParameter));
        simulator.AddForce(p_membrane_force);


        /*
        -----------------------------
        Set intitial
        ----------------------------
        */

        double Modulo;
        for (int i = 0; i < mesh.GetNumNodes(); i++)
        {

            c_vector<double, 3> InitalLocation = cell_population.GetNode(i)->rGetLocation();
            Modulo = (i*3+2) % 20;


            double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]+ InitalLocation[2] * InitalLocation[2]);
            double NewR = R + (RandomNumberGenerator::Instance()->ranf()-0.5)*Scale/2;//R+ (Modulo-10)*Scale/80;

            double Angle;
            // double Scalled_R;

            double theta = acos(InitalLocation[2]/R);
            double gamma = atan(InitalLocation[1]/InitalLocation[0]);
            if (InitalLocation[0] <0)
            {
                theta = -acos(InitalLocation[2]/R);
            }

            double X = NewR * sin(theta) * cos(gamma); 
            double Y = NewR * sin(theta) * sin(gamma);
            double Z = NewR * cos(theta);

            c_vector<double, 3> DeformedLocation = Create_c_vector(X, Y, Z);
            cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
        }

        VtkMeshWriter<2, 3> mesh_writer1(output_dir, "DeformedOriginalMesh", false);
        mesh_writer1.WriteFilesUsingMesh(mesh);

    
     	  simulator.Solve();
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
}


};




#endif /*TESTRELAXATION_HPP_*/


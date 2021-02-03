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
// #include "CellMutationStatesWriter.hpp"
// #include "CellProliferativeTypesWriter.hpp"

#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "MeshBasedCellPopulation.hpp"

// #include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VtkMeshWriter.hpp"
#include "CommandLineArguments.hpp"
// #include "MembraneStiffnessForce.hpp"
#include "PetscSetupAndFinalize.hpp"

// #include "PressureForce.hpp"
#include "XmlTools.hpp"

#include "UblasCustomFunctions.hpp"


#include "VtkMeshReader.hpp"

#include "RemeshingTriggerModifier.hpp"
#include "HistoryDepMutableMesh.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "Honeycomb3DMeshGenerator.hpp"
//  #include "CGALtypedef.h"


class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:

 void TestMappingOnRectangle() throw (Exception)
    {

       std::string mesh_file ="/Users/jcrawshaw/docker-polnet-master/Plexus_Course_2/Plexus.vtu";
        double mesh_scale = 1e-3; 
        double startime = 0;
       
        // Read in the plexus 
        VtkMeshReader<2,3> mesh_reader(mesh_file);
        HistoryDepMutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in m


        // std::string mesh_file ="/Users/jcrawshaw/docker-polnet-master/Triangle.vtu";
        std::string output_dir = "RemeshingComparison/TrialMeshes/Plexus/";
        VtkMeshWriter<2,3> mesh_writer(output_dir, "OriginalConfig", false);
        mesh_writer.WriteFilesUsingMesh(mesh);
        // double startime = 0;
      // Create the cells 
  
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetRelativePath(output_dir);
        cell_population.SetChasteOutputDirectory(output_dir);

    
            for (int i=0; i<mesh.GetNumNodes(); i++)
            {             
                c_vector<double,3> InitalLocation =  cell_population.GetNode(i)->rGetLocation(); 
                c_vector<double,3>  DeformedLocation =Create_c_vector(InitalLocation[0],InitalLocation[1]*2, InitalLocation[2]);
                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;        
            }


    //     // Deform the nodes so i have a inital condition and deformed condition
        cell_population.SetTargetRemeshingElementArea(100*1e-6);
        cell_population.SetTargetRemeshingIterations(50);

           cell_population.RemeshGeometry();
          cell_population.MappingAdaptedMeshToInitalGeometry();
       cell_population.WriteOutMappedInitalConfig();

    }
};


#endif /*TESTRELAXATION_HPP_*/


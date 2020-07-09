#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

#include <cxxtest/TestSuite.h>
#include <cstdio>
#include <ctime>
#include <cmath>

#include "Debug.hpp"

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellIdWriter.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VtkMeshWriter.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "XmlTools.hpp"
#include "UblasCustomFunctions.hpp"
#include "VtkMeshReader.hpp"
#include "HistoryDepMutableMesh.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"



using namespace xsd::cxx::tree;

class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:


void TestPlexusDeformationWithRemeshingw() throw (Exception)
  {
      std::string output_dir = "TestStrainErrorWithRemeshing/";

      float Refinment[13] = {0.0400,0.0600,0.0800,0.1000,0.1200,0.1400,0.1600,0.1800,0.2000,0.2200,0.2400, 0.2600, 0.2800};

      bool CreateStlsInMatlab = 0;
      // Write these into a file to file to be noted in the remeshing process 
      if (CreateStlsInMatlab ==1)
      {
        ofstream file;
        file.open("/Users/jcrawshaw/Documents/testoutput/" + output_dir +"MeshRefinementlist.txt");
        for (unsigned i = 0; i < 13; i++){ file << std::to_string(Refinment[i])+ " "; }
        file.close();

        std::string CreateAllMeshes = "/Applications/MATLAB_R2018a.app/bin/./matlab  -nodisplay -nosplash -nodesktop -r \"run('~/Documents/Projects/MeshMatlab/Remeshing2D/WriteSquareSTL.m');exit;\"";
        std::system(CreateAllMeshes.c_str());

        for (unsigned i = 0; i < 1; i++)
        {

          std::string FileNumber = std::to_string(Refinment[i]);
          FileNumber.erase ( FileNumber.find_last_not_of('0') + 1, std::string::npos );

          std::string ConvertToVtu1 = "meshio-convert /Users/jcrawshaw/Documents/testoutput/"+output_dir +"Meshes/Square_Deformed_"+FileNumber+".stl  /Users/jcrawshaw/Documents/testoutput/"+output_dir +"Meshes/Square_Deformed_"+FileNumber+".vtu";
          std::system(ConvertToVtu1.c_str());

          std::string ConvertToVtu2 = "meshio-convert /Users/jcrawshaw/Documents/testoutput/"+output_dir +"Meshes/Square_Original_"+FileNumber+".stl  /Users/jcrawshaw/Documents/testoutput/"+output_dir +"Meshes/Square_Original_"+FileNumber+".vtu";
          std::system(ConvertToVtu2.c_str());
        }
      }
        
      for (unsigned i = 3; i < 4; i++)
      {

        std::string FileNumber = std::to_string(Refinment[i]);
        FileNumber.erase ( FileNumber.find_last_not_of('0') + 1, std::string::npos );

        // Read in the original Square 
        std::string mesh_file ="/Users/jcrawshaw/Documents/testoutput/"+output_dir +"Meshes/Square_Original_"+FileNumber+".vtu";
        VtkMeshReader<2,3> mesh_reader(mesh_file);
        HistoryDepMutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
  
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetPrintRemeshedIC(1);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.SetRemeshingSoftwear("PreAllocatedMatlabMesh");
        cell_population.PreAllocatedRemeshedMesh("/Users/jcrawshaw/Documents/testoutput/"+output_dir +"Meshes/Square_Deformed_"+FileNumber+".vtu");
        cell_population.SetPaths(output_dir+"_"+FileNumber+"/");

        for (int j=0; j<mesh.GetNumNodes(); j++)
          {
            c_vector<double,3> InitalLocation =  cell_population.GetNode(j)->rGetLocation();
            InitalLocation[0]*=InitalLocation[0]; InitalLocation[1]*=InitalLocation[1];
            cell_population.GetNode(j)->rGetModifiableLocation() = InitalLocation;        
          }
       

       		VtkMeshWriter<2,3> mesh_writer(output_dir+"_"+FileNumber+"/", "config", false);
	      	mesh_writer.WriteFilesUsingMesh(mesh);

       
        cell_population.ExecuteHistoryDependentRemeshing();

        // -- Need to get matlab remeshing tecnhique in here -- put check in that it is 2d i.e all z values are the same 


        // SO I think everthing Is set up, Here I need to remesh for the first time, make sure the filing here is clear 
        // -- I wont determine strain here because cant be bothered, its nicely set up in matlab to do this -- might even call matlab here 

      }
    }



 };



#endif /*TESTRELAXATION_HPP_*/



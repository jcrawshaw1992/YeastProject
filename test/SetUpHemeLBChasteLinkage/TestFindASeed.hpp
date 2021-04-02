#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
// #include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "Debug.hpp"
#include "VtkMeshReader.hpp"
#include "SmartPointers.hpp"
 
#include "AbstractCellBasedTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"

#include <iostream>
#include <string.h>


class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:

  void TestParametersOnBifucation1() throw (Exception)
    {
        TRACE("PI_3")
        for (double i = 5.5; i< 6.5; i=i+0.1 )
        {      
            PRINT_VARIABLE(i)
            std::stringstream out;
            out << i;
            std::string Mesh = out.str();

            std::string mesh_file = "/Users/jcrawshaw/Downloads/AngleVariation_3X3Network/PI_3/mesh_"+Mesh+ ".vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            MutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            // c_vector<double, 3> Location = mesh.GetNode(100)->rGetLocation();
            c_vector<double, 3> Location = mesh.GetNode(100)->rGetLocation();
            PRINT_VECTOR(Location)
        } 
        
    }

    void TestParametersOnBifucation2() throw (Exception)
    {
        TRACE("PI_6")
        for (double i = 5.5; i< 6.5; i=i+0.1 )
        {      

            std::stringstream out;
            out << i;
            std::string Mesh = out.str();

            std::string mesh_file = "/Users/jcrawshaw/Downloads/AngleVariation_3X3Network/PI_6/mesh_"+Mesh+ ".vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            MutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            c_vector<double, 3> Location = mesh.GetNode(50)->rGetLocation();
            PRINT_VECTOR(Location)
        } 
        
    }

    void TestParametersOnBifucation3() throw (Exception)
    {
        TRACE("PI_2.2 ")
        for (double i = 5.5; i< 6.5; i=i+0.1 )
        {      
            std::stringstream out;
            out << i;
            std::string Mesh = out.str();

            std::string mesh_file = "/Users/jcrawshaw/Downloads/AngleVariation_3X3Network/PI_2.2/mesh_"+Mesh+ ".vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            MutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            c_vector<double, 3> Location = mesh.GetNode(50)->rGetLocation();
            PRINT_VECTOR(Location)
        } 
        
    }

//  void TestParametersOnBifucation2() throw (Exception)
//     {
        
//         std::string mesh_file = "/Users/jcrawshaw/Downloads/ScaledMesh.5.vtu";
//         VtkMeshReader<2, 3> mesh_reader(mesh_file);
//         MutableMesh<2, 3> mesh;
//         mesh.ConstructFromMeshReader(mesh_reader);
 
//        // Create the cells 
//         MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//         std::vector<CellPtr> cells;
//         CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
//         cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

//         // Create a cell population
//         HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);

//         c_vector<double, 3> Location = cell_population.GetNode(150)->rGetLocation();
//         PRINT_VECTOR(Location)
//     }


};




#endif /*TESTRELAXATION_HPP_*/


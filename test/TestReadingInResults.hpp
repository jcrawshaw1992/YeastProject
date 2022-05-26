#ifndef TESTCYLINDERVALIDATION_HPP_
#define TESTCYLINDERVALIDATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <ctime>
#include <cmath>
// #include "OffLatticeSimulation.hpp"
// #include "CellsGenerator.hpp"
// #include "FixedG1GenerationalCellCycleModel.hpp"
// #include "GeneralisedLinearSpringForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
// #include "CellBasedSimulationArchiver.hpp"
#include "MeshBasedCellPopulation.hpp"


#include "SmartPointers.hpp"
#include "CellLabel.hpp"
// #include "VtkMeshWriter.hpp"
#include "Debug.hpp"
#include "CommandLineArguments.hpp"
// #include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "CellIdWriter.hpp"
// #include "CellProliferativeTypesWriter.hpp"
// #include "CellMutationStatesWriter.hpp"
// #include "DifferentiatedCellProliferativeType.hpp"
// #include "MembraneShearForce.hpp"
//#include "BetaCateninOneHitCellMutationState.hpp"
#include "PetscSetupAndFinalize.hpp"



#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>



class YeastDeformation : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }


public:

void TestReadInReactionDiffusionResults() //throw (Exception)
    {
        std::string mesh_file = "projects/YeastProject/reaction_diffusion_deforming_membrane/Output/ReactionDiffusionResults.vtu";

        // This data file is in mm
        VtkMeshReader<2,3> mesh_reader(mesh_file);
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        std::vector<double> ConcentrationVector;
        mesh_reader.GetPointData("Concentration", ConcentrationVector);
        PRINT_VECTOR(ConcentrationVector)

    

    }



};

#endif /*TESTCYLINDERVALIDATION_HPP_*/


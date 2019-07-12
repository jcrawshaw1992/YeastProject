#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

#include <cxxtest/TestSuite.h>
#include <cstdio>
#include <ctime>
#include <cmath>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

// Include to avoid errors loading archives
#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
// #include "MembraneStiffnessForce.hpp"

#include "OffLatticeSimulation.hpp"
#include "AppliedForceModifier.hpp"
#include "SpringLengthModifier.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
// #include "FixedDurationGenerationBasedCellCycleModel.hpp"

#include "FixedG1GenerationalCellCycleModel.hpp"



#include "GeneralisedLinearSpringForce.hpp"
#include "LinearSpringWithRestLengthDependentSpringConstantsForce.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "AppliedForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VtkMeshWriter.hpp"
#include "Debug.hpp"
#include "CommandLineArguments.hpp"
#include "PressureForce.hpp"
#include "PetscSetupAndFinalize.hpp"



 #include "projects/Jess/src/DeformableModel/MembraneShearForce.hpp"
 #include "projects/Jess/src/DeformableModel/MembraneStiffnessForce.hpp"
 #include "projects/Jess/src/DeformableModel/MembraneSurfaceForce.hpp"

// //  #include "HemodynamicPressure.hpp"
//  #include "AppliedForceModifier.hpp"

class TestRunFlowInPipe : public AbstractCellBasedTestSuite
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

	void TestLoadandRunPipe() throw (Exception)
	{

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-start_time"));
		double start_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-start_time");

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-duration"));
		double duration = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-duration");

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-traction_file"));
		std::string traction_file = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-traction_file");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-mesh_scale"));
        double mesh_scale = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-mesh_scale").c_str());

		


		TRACE("Before Load");
		 std::string output_dir = "./"; //"./"
		OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(output_dir, start_time);
		TRACE("After Load");

		 p_simulator->SetEndTime(start_time + duration);
		//  p_simulator->SetEndTime(10);
		// std::string PressureFile = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/DeformingCylindricalMembraneWithFluidFlow/results/Extracted/surface-pressure.xtr"; 
      
		// Update the tractions in the Modifier assumes the AppliedForceModifier is the first in the vector
	    std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
		TRACE("A");
	    //assert(boost::dynamic_pointer_cast<AppliedForceModifier<2,3> >(*iter));
	    boost::shared_ptr<AppliedForceModifier<2,3> > p_force_modifier = boost::static_pointer_cast<AppliedForceModifier<2, 3> >(*iter);
 		TRACE("a");
		p_force_modifier->SetResetTractionsOnCells(true, traction_file);
		TRACE("B");

     	p_simulator->Solve();

		TRACE("Before Save");
		CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(p_simulator);
		TRACE("After Save");

		std::stringstream out;
		out << start_time;
		std::string mesh_size = out.str();


		VtkMeshWriter<2,3> mesh_writer2(output_dir, "config", false);
     	MutableMesh<2,3>* p_mesh2 = &(dynamic_cast<MeshBasedCellPopulation<2,3>*>(&(p_simulator->rGetCellPopulation()))->rGetMesh());
        p_mesh2->Scale(1.0/mesh_scale,1.0/mesh_scale,1.0/mesh_scale); // so distances are back in original scale
		mesh_writer2.WriteFilesUsingMesh(*p_mesh2);

		delete p_simulator;
	}

};

#endif /*TESTRELAXATION_HPP_*/



	// Create an Applied Force modifier to couple to Flow
			// boost::shared_ptr<AppliedForceModifier<2,3> > p_traction_modifier(new AppliedForceModifier<2,3>());
			// p_traction_modifier->SetResetTractionsOnCells(true, traction_file);
			// p_traction_modifier->SetEdgeDivisionThreshold(1e10);
			// p_simulator->AddSimulationModifier(p_traction_modifier);

			// boost::shared_ptr<AppliedForce<2,3>> p_pressure_force(new AppliedForce<2,3>());
			// p_simulator->AddForce(p_pressure_force);


		

		// // Update the traction force assuming the hemeTract is the first in the vector
	    // std::vector<boost::shared_ptr<AbstractForce<2, 3> > > ForceCollection = p_simulator->rGetForceCollection();
		//  double Thesize = p_simulator->rGetForceCollection().size();
		// //  PRINT_VARIABLE(Thesize);

// 		boost::shared_ptr<HemodynamicTraction> p_Traction = boost::static_pointer_cast<HemodynamicTraction>(p_simulator->rGetForceCollection()[0]);

// TRACE("c");

		//  std::string traction_file = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/RegularCylinder/results/Extracted/surface-tractions.xtr";
        //  p_Traction->LoadTractionFromFile(traction_file);

		// boost::shared_ptr<PressureForce> p_pressure_force(new PressureForce());
        // p_simulator->AddForce(p_pressure_force);





     	// VtkMeshWriter<2,3> mesh_writer(output_dir, "config"+mesh_size, false);
     	// MutableMesh<2,3>* p_mesh = &(dynamic_cast<MeshBasedCellPopulation<2,3>*>(&(p_simulator->rGetCellPopulation()))->rGetMesh());
        // p_mesh->Scale(1.0/mesh_scale,1.0/mesh_scale,1.0/mesh_scale); // so distances are back in original scale
		// mesh_writer.WriteFilesUsingMesh(*p_mesh);
	

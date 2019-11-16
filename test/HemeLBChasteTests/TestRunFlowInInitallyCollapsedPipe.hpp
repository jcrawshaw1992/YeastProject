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



// #include "GeneralisedLinearSpringForce.hpp"
// #include "LinearSpringWithRestLengthDependentSpringConstantsForce.hpp"
// #include "LinearSpringWithVariableSpringConstantsForce.hpp"
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


#include "projects/VascularRemodelling/src/MembraneForces/MembraneShearForce.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"
 #include "MembraneSurfaceForce.hpp"

#include "EmptyBasementMatrix.hpp"
#include "LostEndothelialCell.hpp"
#include "HasEndothelialCell.hpp"
 #include "EdgeCorrectionForce.hpp"

// //  #include "HemodynamicPressure.hpp"
//  #include "AppliedForceModifier.hpp"



class RadialForce : public AbstractForce<2, 3>
{
private:
    double mStrength;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractForce<2, 3> >(*this);
        archive& mStrength;
    }

public:
    RadialForce(double strength = 1.0)
            : AbstractForce<2, 3>(),
              mStrength(strength)
    {
        assert(mStrength > 0.0);
    }


    void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation)
    {
        // Helper variables
        MeshBasedCellPopulation<2,3>* p_cell_population = static_cast<MeshBasedCellPopulation<2,3>*>(&rCellPopulation);


        // Calculate midpoint
       for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);
     
            c_vector<double, 3> force = zero_vector<double>(3);

            //Calculate cell normal (average of element normals)
            c_vector<double,3> normal = zero_vector<double>(3);

            std::set<unsigned>&  containing_elements = p_node->rGetContainingElementIndices();
            assert(containing_elements.size()>0);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                    iter != containing_elements.end();
                    ++iter)
            {
                // Negative as normals point inwards for these surface meshes
                normal +=  p_cell_population->rGetMesh().GetElement(*iter)->CalculateNormal();
            }
            normal /= norm_2(normal);

                

            double cell_area = rCellPopulation.GetVolumeOfCell(*cell_iter);
            force = mStrength * normal; // cell_location / norm_2(cell_location);
            cell_iter->GetCellData()->SetItem("area", cell_area);
            
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);

            cell_iter->GetCellData()->SetItem("applied_force_x", force[0]);
            cell_iter->GetCellData()->SetItem("applied_force_y", force[1]);
            cell_iter->GetCellData()->SetItem("applied_force_z", force[2]);
            cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(force));

            cell_iter->GetCellData()->SetItem("norm_x", normal[0]);
            cell_iter->GetCellData()->SetItem("norm_y", normal[1]);
            cell_iter->GetCellData()->SetItem("norm_z", normal[2]);

        }

    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
        AbstractForce<2,3>::OutputForceParameters(rParamsFile);
    }
};


#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(RadialForce)

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RadialForce)





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

		


		// TRACE("Before Load");
		 std::string output_dir = "./"; //"./"
		OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(output_dir, start_time);
		// TRACE("After Load");

		

		 p_simulator->SetEndTime(start_time + duration);
		//  p_simulator->SetEndTime(10);
		// std::string PressureFile = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/DeformingCylindricalMembraneWithFluidFlow/results/Extracted/surface-pressure.xtr"; 
      
		// Update the tractions in the Modifier assumes the AppliedForceModifier is the first in the vector
	    std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
	    //assert(boost::dynamic_pointer_cast<AppliedForceModifier<2,3> >(*iter));
	    boost::shared_ptr<AppliedForceModifier<2,3> > p_force_modifier = boost::static_pointer_cast<AppliedForceModifier<2, 3> >(*iter);
		p_force_modifier->SetResetTractionsOnCells(true, traction_file);
		p_force_modifier->SetupVessel(p_simulator->rGetCellPopulation(), p_simulator->GetOutputDirectory() );
		

		  //  -----------------------------
         //  Radial pressure force  
         //  ----------------------------

            // MAKE_PTR_ARGS(RadialForce, p_radial_force, ( 1.0666e2));
            // p_simulator->AddForce(p_radial_force);



		// Update the traction force assuming the Applied force is the first in the vector
	    // // std::vector<boost::shared_ptr<AbstractForce<2, 3> > > ForceCollection = p_simulator->rGetForceCollection();
		//  double Thesize = p_simulator->rGetForceCollection().size();
		//  PRINT_VARIABLE(Thesize);

		// boost::shared_ptr<MembraneShearForce> p_shear_force = boost::static_pointer_cast<MembraneShearForce>(p_simulator->rGetForceCollection()[1]);
		// p_shear_force->UpdateMembraneProperties(p_simulator->rGetCellPopulation());

		// boost::shared_ptr<MembraneSurfaceForce> p_area_force = boost::static_pointer_cast<MembraneSurfaceForce>(p_simulator->rGetForceCollection()[2]);
		// p_area_force->UpdateMembraneSurfaceForceProperties(p_simulator->rGetCellPopulation());

		// boost::shared_ptr<MembraneStiffnessForce> p_bending_force = boost::static_pointer_cast<MembraneStiffnessForce>(p_simulator->rGetForceCollection()[3]);
		// p_bending_force->UpdateMembraneStiffnessProperties(p_simulator->rGetCellPopulation());

		
		// boost::shared_ptr<MembraneShearForce> p_shear_force(new MembraneShearForce());
	   TRACE("Before Solve");
     	p_simulator->Solve();

		// TRACE("Before Save");
		CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(p_simulator);
		TRACE("After Solve");

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


		


		//  std::string traction_file = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/RegularCylinder/results/Extracted/surface-tractions.xtr";
        //  p_Traction->LoadTractionFromFile(traction_file);

		// boost::shared_ptr<PressureForce> p_pressure_force(new PressureForce());
        // p_simulator->AddForce(p_pressure_force);





     	// VtkMeshWriter<2,3> mesh_writer(output_dir, "config"+mesh_size, false);
     	// MutableMesh<2,3>* p_mesh = &(dynamic_cast<MeshBasedCellPopulation<2,3>*>(&(p_simulator->rGetCellPopulation()))->rGetMesh());
        // p_mesh->Scale(1.0/mesh_scale,1.0/mesh_scale,1.0/mesh_scale); // so distances are back in original scale
		// mesh_writer.WriteFilesUsingMesh(*p_mesh);
	

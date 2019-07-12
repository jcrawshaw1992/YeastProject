#ifndef TESTSPRINGS_HPP_
#define TESTSPRINGS_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include <cstdio>
#include <ctime>
#include <cmath>


#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "AbstractForce.hpp"

#include "OffLatticeSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "LinearSpringWithRestLengthDependentSpringConstantsForce.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "AppliedForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"
#include "Debug.hpp"


class FixedPointBoundaryCondition : public AbstractCellPopulationBoundaryCondition<1>
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<1> >(*this);
    }

public:
    FixedPointBoundaryCondition(AbstractCellPopulation<1>* pCellPopulation)
        : AbstractCellPopulationBoundaryCondition<1>(pCellPopulation)
    {
    }

    void ImposeBoundaryCondition(const std::map<Node<1>*, c_vector<double, 1> >& rOldLocations)
    {
        for (AbstractCellPopulation<1>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
            Node<1>* p_node = this->mpCellPopulation->GetNode(node_index);
            double old_x_coordinate = rOldLocations.find(p_node)->second[0];
            if (old_x_coordinate <1e-5)
            {
                p_node->rGetModifiableLocation()[0] = 0.0;
            }
        }
    }

    bool VerifyBoundaryCondition()
    {
        bool condition_satisfied = true;

        for (AbstractCellPopulation<1>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, 1> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
            double x_coordinate = cell_location(0);

            if (x_coordinate < 0.0)
            {
                //condition_satisfied = false;
                //break;
            }
        }
        return condition_satisfied;
    }

    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
    {
        AbstractCellPopulationBoundaryCondition<1>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
    }
};

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const FixedPointBoundaryCondition * t, const BOOST_PFTO unsigned int file_version)
        {
            const AbstractCellPopulation<1>* const p_cell_population = t->GetCellPopulation();
            ar << p_cell_population;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, FixedPointBoundaryCondition * t, const unsigned int file_version)
        {
            AbstractCellPopulation<1>* p_cell_population;
            ar >> p_cell_population;

            ::new(t)FixedPointBoundaryCondition(p_cell_population);
        }
    }
}

class MotileForce : public AbstractForce<1>
{
private:

    double mStrength;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<1> >(*this);
        archive & mStrength;
    }

public:
    MotileForce(double strength=1.0)
        : AbstractForce<1>(),
          mStrength(strength)
    {
        assert(mStrength > 0.0);
    }

    void AddForceContribution(AbstractCellPopulation<1>& rCellPopulation)
    {
        /* Inside the method, we loop over nodes, and add a constant vector to
         * each node, in the negative ''y''-direction and of magnitude {{{mStrength}}}.
         */
        c_vector<double, 1> force = zero_vector<double>(1);
        force(0) = mStrength;

        for (AbstractCellPopulation<1>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
		{
			if (cell_iter->HasCellProperty<CellLabel>())
			{
				unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
				rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);
			}
		}
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
        AbstractForce<1>::OutputForceParameters(rParamsFile);
    }
};

class RadialForce : public AbstractForce<2>
{
private:

    double mStrength;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<2> >(*this);
        archive & mStrength;
    }

public:
    RadialForce(double strength=1.0)
        : AbstractForce<2>(),
          mStrength(strength)
    {
        assert(mStrength > 0.0);
    }

    void AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
    {
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
		{
        	unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        	c_vector<double,2> cell_location = rCellPopulation.GetNode(node_index)->rGetLocation();

        	c_vector<double, 2> force = zero_vector<double>(2);

        	if (norm_2(cell_location) > 1e-10)
        	{
        		double radius  = norm_2(cell_location);

        		// multiply by radius as Force = P * 2 * R * sin(pi/n)
        		//force = mStrength * radius * cell_location / norm_2(cell_location);
        		// multiply by radius^2 as Force = P * 4 * R^2 * sin(pi/n)^2
        		force = mStrength * radius * radius * cell_location / norm_2(cell_location);
        	//	PRINT_3_VARIABLES(mStrength,radius,norm_2(force));
        	}
//        	PRINT_VARIABLE(rCellPopulation.GetNode(node_index)->rGetAppliedForce());

        	rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);
		}
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
        AbstractForce<2>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(FixedPointBoundaryCondition)
CHASTE_CLASS_EXPORT(MotileForce)
CHASTE_CLASS_EXPORT(RadialForce)

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FixedPointBoundaryCondition)
CHASTE_CLASS_EXPORT(MotileForce)
CHASTE_CLASS_EXPORT(RadialForce)

class TestSprings : public AbstractCellBasedTestSuite
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

//    /*
//     * This test sets up a 1d chain of springs and checks that as we apply a force to the end cell
//     * the chain relaxes to the same rest length. Independent of length.
//     */
//    void noTestChainOfSprings() throw (Exception)
//	{
//    	unsigned num_springs = 10;
//    	double initial_length = 10.0;
//
//    	std::vector<Node<1>*> nodes;
//    	for (unsigned i = 0; i<num_springs+1; i++)
//    	{
//    		nodes.push_back(new Node<1>(i, false, (i)*initial_length/(double)num_springs));
//    	}
//    	MutableMesh<1,1> mesh(nodes);
//
//		// Create cells
//		MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//		std::vector<CellPtr> cells;
//		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);
//
//		// Set up cells so that end cell is labelled
//		boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
//		cells[cells.size()-1]->AddCellProperty(p_label);
//
//		// Create a cell population
//		MeshBasedCellPopulation<1> cell_population(mesh, cells);
//
//		// So using variable rest length
//		cell_population.CalculateRestLengths();
//
//		// So can visualize results //Note doesnt work as vtk doesnt work in 1d!!!
//		cell_population.SetWriteVtkAsPoints(true);
//		cell_population.SetOutputMeshInVtk(true);
//
//		// Set up cell-based simulation
//		OffLatticeSimulation<1> simulator(cell_population);
//		simulator.SetOutputDirectory("TestChainOfCells");
//		simulator.SetEndTime(10.0);
//		simulator.SetDt(0.001);
//		simulator.SetSamplingTimestepMultiple(1);
//		simulator.SetUpdateCellPopulationRule(false); // No remeshing.
//
//		// Create a force law and pass it to the simulation
//		MAKE_PTR(LinearSpringWithRestLengthDependentSpringConstantsForce<1>, p_linear_force);
//		simulator.AddForce(p_linear_force);
//
//		//Fix a cell at x=0 and apply a force to the cell at x=1
//		boost::shared_ptr<FixedPointBoundaryCondition> p_condition(new FixedPointBoundaryCondition(&cell_population));
//		simulator.AddCellPopulationBoundaryCondition(p_condition);
//
//        MAKE_PTR_ARGS(MotileForce, p_motile_force, (1.0));
//        simulator.AddForce(p_motile_force);
//
//		// Run simulation for a short time
//		simulator.Solve();
//
//		//Check all cells are evenly spaced
//		double expected_seperation = simulator.rGetCellPopulation().GetNode(num_springs)->rGetLocation()[0] / (double)num_springs;
//		for (unsigned i = 0; i<num_springs; i++)
//		{
//			double seperation = simulator.rGetCellPopulation().GetNode(i+1)->rGetLocation()[0] -simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0];
//			TS_ASSERT_DELTA(seperation,expected_seperation, 1e-5);
//		}
//		// Check Position of end cell
//		TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(num_springs)->rGetLocation()[0],initial_length*31.0/30.0,1e-5);
//	}


	/*
	 * This test sets up a 1d chain of springs and checks that as we apply a force to the end cell
	 * the chain relaxes to the same rest length. Independent of length.
	 */
	void TestCircleOfSprings() throw (Exception)
	{
		unsigned num_springs = 30;
		double initial_radius = 0.015;

		std::vector<Node<2>*> nodes;
		std::vector<unsigned> location_indices;
		for (unsigned i = 0; i<num_springs; i++)
		{
			location_indices.push_back(i);
			double angle = (double)i*2*M_PI/(double)num_springs;
			nodes.push_back(new Node<2>(i, true, initial_radius * sin(angle), initial_radius * cos(angle)));
		}
		nodes.push_back(new Node<2>(num_springs, false, 0.0,0.0));

		MutableMesh<2,2> mesh(nodes);

		// Create cells
		MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
		std::vector<CellPtr> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_differentiated_type);

		// Create a cell population
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(mesh, cells, location_indices);

		// So using variable rest length
		cell_population.CalculateRestLengths();

		// So can visualize results in VTK
		cell_population.SetWriteVtkAsPoints(true);
		cell_population.SetOutputMeshInVtk(true);

		// Set up cell-based simulation
		OffLatticeSimulation<2> simulator(cell_population);
		simulator.SetOutputDirectory("TestCircleOfCells");
		simulator.SetEndTime(1.0);
		simulator.SetDt(0.0001);
		simulator.SetSamplingTimestepMultiple(100);
		simulator.SetUpdateCellPopulationRule(false); // No remeshing.

		// Create a force law and pass it to the simulation
		//MAKE_PTR(LinearSpringWithRestLengthDependentSpringConstantsForce<2>, p_linear_force);
	    MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
	    p_linear_force->SetMeinekeSpringStiffness(1e3); //15

		simulator.AddForce(p_linear_force);

		double pressure = 1.0666e4; //3
		double scaled_force = pressure*2.0*sin(M_PI/(double)num_springs)*2.0*sin(M_PI/(double)num_springs); // Area, Also change to R^2 in Force
		//double scaled_force = pressure*2.0*sin(M_PI/(double)num_springs); // Length , also change to R in Force
		MAKE_PTR_ARGS(RadialForce, p_radial_force, (scaled_force));
	    simulator.AddForce(p_radial_force);

		// Run simulation for a short time
		simulator.Solve();
//
//		//Check all cells are evenly spaced
//		double expected_seperation = simulator.rGetCellPopulation().GetNode(num_springs)->rGetLocation()[0] / (double)num_springs;
//		for (unsigned i = 0; i<num_springs; i++)
//		{
//			double seperation = simulator.rGetCellPopulation().GetNode(i+1)->rGetLocation()[0] -simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0];
//			TS_ASSERT_DELTA(seperation,expected_seperation, 1e-5);
//		}
//		// Check Position of end cell
//		TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(num_springs)->rGetLocation()[0],initial_length*31.0/30.0,1e-5);
//

	}
};


#endif /*TESTSPRINGS_HPP_*/


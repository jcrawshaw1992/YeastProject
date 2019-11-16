#ifndef TESTCYLINDERVALIDATION_HPP_
#define TESTCYLINDERVALIDATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <ctime>
#include <cmath>
#include "OffLatticeSimulation.hpp"
// #include "CylindricalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"
#include "VtkMeshWriter.hpp"
#include "Debug.hpp"
#include "CommandLineArguments.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "BetaCateninOneHitCellMutationState.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneShearForce.hpp"
#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"
 #include "MembraneSurfaceForce.hpp"

#include "EdgeCorrectionForce.hpp"
// #include "PressureForce.hpp"
#include "XmlTools.hpp"

#include "UblasCustomFunctions.hpp"
#include <boost/serialization/base_object.hpp>

// #include "ConstantPressure.hpp"



using namespace xsd::cxx::tree;

static const double M_TIME_FOR_SIMULATION = 10000; 
static const double M_TIME_STEP = 0.02; 
static const double M_SAMPLING = 5000; 
static const double M_BendingConstant = 1e-9;

static const double ElasticShearModulus = 4.4e-07;
static const double AreaDilationModulus = 0.9e-6;
// static const double M_BendingCOnstant = 0.75e-11 ;
static const double Area_constant = 0.9e-8 ;
static const double pressure = 2.066e1; //1.066e4; // to match 80mmhg



class RadialForce : public AbstractForce<2,3>
{
private:

    double mStrength;
    double mNz;
    double mNc;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<2,3> >(*this);
        archive & mStrength;
    }

public:
    RadialForce(double strength=1.0, double Nc =1.0, double Nz=1.0)
        : AbstractForce<2,3>(),
          mStrength(strength),
          mNz(Nz),
          mNc(Nc)
    {
        assert(mStrength > 0.0);
    }

    void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation)
    {
        // Helper variables
        MeshBasedCellPopulation<2,3>* p_cell_population = static_cast<MeshBasedCellPopulation<2,3>*>(&rCellPopulation);


        // Calculate midpoint
        c_vector<double,3> centroid = zero_vector<double>(3);
        for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            centroid += rCellPopulation.GetNode(node_index)->rGetLocation();
        }
        centroid /= rCellPopulation.GetNumRealCells();

        for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);
            c_vector<double,3> cell_location = p_node->rGetLocation() - centroid;
            cell_location(2) = 0.0;
            cell_location /=norm_2(cell_location);
            c_vector<double, 3> force = zero_vector<double>(3);


            double cell_area = rCellPopulation.GetVolumeOfCell(*cell_iter);

            force =  mStrength * cell_area  * cell_location;

            cell_iter->GetCellData()->SetItem("area", cell_area);
            
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);

        }


   for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
     cell_iter != rCellPopulation.End();
     ++cell_iter)
    {
    
        unsigned ReferenceNode = 0;

        if ((cell_iter)->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())// | (is_boundary_node &&  node_index > (mNc *mNz)- mNc))
        {      
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);

            // PRINT_2_VARIABLES(counter, node_index);
            if (node_index < mNc + 1) // if on lower edge
            {
                ReferenceNode = node_index + (2 * mNc); // select node from two rows up
            }
            else if (node_index > mNc) // if on upper edge
            {
                ReferenceNode = node_index - (2 * mNc); // select node from two rows down
            }
            Node<3>* pReferenceNode = p_cell_population->rGetMesh().GetNode(ReferenceNode);

            // TRACE("clear the force");
            pNode->ClearAppliedForce(); // remove the already present force at this node
            pNode->AddAppliedForceContribution(pReferenceNode->rGetAppliedForce()); // Add the new force

    }

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


class BendingValidation_CorrectedForDrag : public AbstractCellBasedTestSuite
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


     void TestCylinderStretched() throw (Exception)
    {

        unsigned N_D[3] = {40,60 ,80};//,160};
       
        for (unsigned N_D_index = 0; N_D_index < 3; N_D_index++)
        {
            double N_Z =  N_D[N_D_index] *0.75;
                Honeycomb3DCylinderMeshGenerator generator(N_D[N_D_index], N_Z, 1.5e-3, 12e-3);
                MutableMesh<2,3>* p_mesh = generator.GetMesh();
                p_mesh->Translate(-6.0e-3*unit_vector<double>(3,2));
                double MinZ =  -6e-3;
                double MaxZ =  6e-3;
                std::stringstream out;
                out << N_D[N_D_index] << "_" << N_Z;
                std::string mesh_size = out.str();
                std::string output_directory = "BendingValidation/Stretched/" + mesh_size;

                // Create cells
                MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
                MAKE_PTR(BetaCateninOneHitCellMutationState, p_state); //Mutation to mark edge nodes
                std::vector<CellPtr> cells;
                CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

                // Create a cell population
                MeshBasedCellPopulation<2,3> cell_population(*p_mesh, cells);
                cell_population.SetWriteVtkAsPoints(true);
                cell_population.SetOutputMeshInVtk(true);

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellProliferativeTypesWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();

                // Set up cell-based simulation
                OffLatticeSimulation<2,3> simulator(cell_population);
                simulator.SetOutputDirectory(output_directory);
                simulator.SetEndTime(M_TIME_FOR_SIMULATION); //50
                simulator.SetDt(M_TIME_STEP);
                simulator.SetSamplingTimestepMultiple(M_SAMPLING);
                simulator.SetUpdateCellPopulationRule(false); // No remeshing.



                c_vector<long double, 3> Node_location;
                for (AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
                        cell_iter != cell_population.End();
                        ++cell_iter)
                    {
                        unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
                        Node_location = cell_population.GetNode(node_index)->rGetLocation();
                        if (Node_location[2] == MinZ) // These are the nodes along the lower edge and need to be marked as mutated
                        {
                            cell_iter->SetMutationState(p_state);
                        }
                        else if (Node_location[2] == MaxZ) // These are the nodes along the upper edge and need to be marked as mutated
                        {
                            cell_iter->SetMutationState(p_state);
                        }
                    }
                

                // Create a force law to apply radial pressure force
                double TransmuralPressure = pressure;
                MAKE_PTR_ARGS(RadialForce, p_radial_force, (pressure,N_D[N_D_index], N_Z));
                simulator.AddForce(p_radial_force);

                /*
                -----------------------------
                Bending Force
                ----------------------------
                */
                boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
                p_membrane_force->SetMembraneStiffness(M_BendingConstant,N_D[N_D_index], N_Z );
                p_membrane_force->SetupInitialMembrane(*p_mesh, cell_population);
                simulator.AddForce(p_membrane_force);

                // boost::shared_ptr<EdgeCorrectionForce> p_EdgeCorrectionForce(new EdgeCorrectionForce());
                // p_EdgeCorrectionForce->SetMeshType(1, N_D[N_D_index], N_Z );
                // simulator.AddForce(p_EdgeCorrectionForce);    


    
                // Create a plane boundary to represent the inlet and pass them to the simulation
                boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition_1(new FixedRegionBoundaryCondition<2,3>(&cell_population,6.0e-3*unit_vector<double>(3,2),unit_vector<double>(3,2),10));
                simulator.AddCellPopulationBoundaryCondition(p_condition_1);

                boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition_2(new FixedRegionBoundaryCondition<2,3>(&cell_population,6.0e-3*unit_vector<double>(3,2),-unit_vector<double>(3,2),10));
                simulator.AddCellPopulationBoundaryCondition(p_condition_2);

                simulator.Solve();

                // To reset before looping: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
    
        }
    }

      void TestCylinderEqui() throw (Exception)
    {

        unsigned N_D[3] = {40,60 ,80};//,160};
       
        for (unsigned N_D_index = 0; N_D_index < 3; N_D_index++)
        {
            double N_Z =  N_D[N_D_index] *1.5;
                Honeycomb3DCylinderMeshGenerator generator(N_D[N_D_index], N_Z, 1.5e-3, 12e-3);
                MutableMesh<2,3>* p_mesh = generator.GetMesh();
                p_mesh->Translate(-6.0e-3*unit_vector<double>(3,2));
                double MinZ =  -6e-3;
                double MaxZ =  6e-3;
                std::stringstream out;
                out << N_D[N_D_index] << "_" << N_Z;
                std::string mesh_size = out.str();
                std::string output_directory = "BendingValidation/Equilaterial/" + mesh_size;

                // Create cells
                MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
                MAKE_PTR(BetaCateninOneHitCellMutationState, p_state); //Mutation to mark edge nodes
                std::vector<CellPtr> cells;
                CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

                // Create a cell population
                MeshBasedCellPopulation<2,3> cell_population(*p_mesh, cells);
                cell_population.SetWriteVtkAsPoints(true);
                cell_population.SetOutputMeshInVtk(true);

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellProliferativeTypesWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();

                // Set up cell-based simulation
                OffLatticeSimulation<2,3> simulator(cell_population);
                simulator.SetOutputDirectory(output_directory);
                simulator.SetEndTime(M_TIME_FOR_SIMULATION); //50
                simulator.SetDt(M_TIME_STEP);
                simulator.SetSamplingTimestepMultiple(M_SAMPLING);
                simulator.SetUpdateCellPopulationRule(false); // No remeshing.



                c_vector<long double, 3> Node_location;
                for (AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
                        cell_iter != cell_population.End();
                        ++cell_iter)
                    {
                        unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
                        Node_location = cell_population.GetNode(node_index)->rGetLocation();
                        if (Node_location[2] == MinZ) // These are the nodes along the lower edge and need to be marked as mutated
                        {
                            cell_iter->SetMutationState(p_state);
                        }
                        else if (Node_location[2] == MaxZ) // These are the nodes along the upper edge and need to be marked as mutated
                        {
                            cell_iter->SetMutationState(p_state);
                        }
                    }
                

                // Create a force law to apply radial pressure force
                double TransmuralPressure = pressure;
                MAKE_PTR_ARGS(RadialForce, p_radial_force, (pressure,N_D[N_D_index], N_Z));
                simulator.AddForce(p_radial_force);

                /*
                -----------------------------
                Bending Force
                ----------------------------
                */
                boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
                p_membrane_force->SetMembraneStiffness(M_BendingConstant,N_D[N_D_index], N_Z );
                p_membrane_force->SetupInitialMembrane(*p_mesh, cell_population);
                simulator.AddForce(p_membrane_force);

    
                // Create a plane boundary to represent the inlet and pass them to the simulation
                boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition_1(new FixedRegionBoundaryCondition<2,3>(&cell_population,6.0e-3*unit_vector<double>(3,2),unit_vector<double>(3,2),10));
                simulator.AddCellPopulationBoundaryCondition(p_condition_1);

                boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition_2(new FixedRegionBoundaryCondition<2,3>(&cell_population,6.0e-3*unit_vector<double>(3,2),-unit_vector<double>(3,2),10));
                simulator.AddCellPopulationBoundaryCondition(p_condition_2);

                simulator.Solve();

                // To reset before looping: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
    
        }
    }



};

#endif /*TESTCYLINDERVALIDATION_HPP_*/


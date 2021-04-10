#ifndef TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_
#define TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_
//TestCylindricalGrowthDeformableMembrane
#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>
#include <ctime>
#include "CellsGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
//#include "GeneralisedLinearSpringForce_CorrectedForArea.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "CellIdWriter.hpp"
#include "CellLabel.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CommandLineArguments.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VtkMeshWriter.hpp"

#include "MembraneSurfaceForceTest.hpp"
//#include "MembraneSurfaceForce.hpp"
#include "MembraneStiffnessForce.hpp"
// #include "MembraneShearForce.hpp"
#include "NewForce.hpp"

//#include "TriangleVertexMesh.hpp"

#include "Debug.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "VertexMeshReader.hpp"

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

    c_vector<long double, 3> MatrixMultiplication(c_vector<c_vector<long double, 3>, 3> Matrix, c_vector<long double, 3> Vector)
    {

        c_vector<long double, 3> Answer;

        for (int i = 0; i < 3; i++)
        {
            Answer[i] = inner_prod(Matrix[i], Vector);
        }

        return Answer;
    }

    void AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
    {
        // Helper variables
        MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

        for (typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
             elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
             ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();

            if (elem_index == 1)
            {
                double theta = mStrength;
                c_vector<c_vector<long double, 3>, 3> Rotation;
                
                Rotation[0] = Create_c_vector(cos(theta), 0, sin(theta));
                Rotation[1] = Create_c_vector(0, 1, 0);
                Rotation[2] = Create_c_vector(-sin(theta), 0, cos(theta));

                Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));

                c_vector<long double, 3> x3 = pNode2->rGetLocation();
                c_vector<long double, 3> X3 = MatrixMultiplication(Rotation, x3);
                c_vector<long double, 3> Force3 = (-x3 + X3)/0.2;
                PRINT_VECTOR(x3);
                PRINT_VECTOR(X3);
                PRINT_VECTOR(Force3)

                pNode2->AddAppliedForceContribution(Force3);
            }
        }
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
        AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(RadialForce)

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RadialForce)

class CylinderValidation_CorrectedForDrag : public AbstractCellBasedTestSuite
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
        double elapsed_time = (time - mLastStartTime) / (CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:
    void TestMembraneForce() throw(Exception)
    {

        static const unsigned P[12] = {  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
        double P_max = 13;
        for (int P_index = 0; P_index < P_max - 1; P_index++)
        {
            TrianglesMeshReader<2, 3> mesh_reader("mesh/test/data/Two2dTriangles_in_3d");
            MutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            
            std::stringstream out;
            out << P[P_index];
            std::string Force = out.str();
            std::string output_directory = "BendingTest/" + Force;

            PRINT_VARIABLE(P[P_index]);
            // Create cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, 4, p_differentiated_type);

            // Create a cell population
            MeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellProliferativeTypesWriter>();
            cell_population.AddCellWriter<CellMutationStatesWriter>();

            cell_population.CalculateRestLengths();

            // Set up cell-based simulation
            OffLatticeSimulation<2, 3> simulator(cell_population);
            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(10); //50
            simulator.SetDt(0.002);
            simulator.SetSamplingTimestepMultiple(20);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.


            /*
            -----------------------------
            Bending Force
            ----------------------------
*/

            double membrane_constant = 15;//1e-12;
            boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
            p_membrane_force->SetupInitialMembrane(mesh);
            p_membrane_force->SetMembraneStiffness(membrane_constant);
            simulator.AddForce(p_membrane_force);

            /*
        
            -----------------------------
            Radial pressure force 
            ----------------------------
*/
            MAKE_PTR_ARGS(RadialForce, p_radial_force, (P[P_index]*M_PI/6));
            simulator.AddForce(p_radial_force);

            simulator.Solve();

            // To reset before looping: this is usually done by the SetUp and TearDown methods
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }
    }
};

#endif


#ifndef ConstantPressure_HPP_
#define ConstantPressure_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>
#include <ctime>
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "Debug.hpp"

#include "OffLatticeSimulation.hpp"

#include "SmartPointers.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

class ConstantPressure : public AbstractForce<2, 3>
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
    ConstantPressure(double strength = 1.0)
            : AbstractForce<2, 3>(),
              mStrength(strength)
    {
        assert(mStrength > 0.0);
    }

    void AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
    {
        MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
          
        for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {

            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);


            c_vector<long double, 3> Normal = zero_vector<long double>(3);
            // long double Area = 0;
            std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
            assert(containing_elements.size() > 0);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                 iter != containing_elements.end();
                 ++iter)
            {
                Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
                Node<3>* pNode1= p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
                Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

                c_vector<long double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
                c_vector<long double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

                c_vector<long double, 3> normalVector = VectorProduct(vector_12, vector_13);
                // Area+= norm_2(normalVector)/6;
                 Normal += (normalVector/norm_2(normalVector));
            }
            Normal /= norm_2(Normal );
             c_vector<long double, 3> ConstantPressure = mStrength * Normal;  // / norm_2(cell_location);
         //   PRINT_VECTOR(force);

            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(ConstantPressure);

            cell_iter->GetCellData()->SetItem("Compressive_force_x", ConstantPressure[0]);
            cell_iter->GetCellData()->SetItem("Compressive_force_y", ConstantPressure[1]);
            cell_iter->GetCellData()->SetItem("Compressive_force_z", ConstantPressure[2]);
            cell_iter->GetCellData()->SetItem("CompressivePressure", norm_2(ConstantPressure));
        }
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
        AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(ConstantPressure)

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ConstantPressure)
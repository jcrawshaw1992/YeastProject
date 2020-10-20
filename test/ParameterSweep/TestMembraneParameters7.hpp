#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
// #include <cstdio>
// #include <ctime>
// #include <cmath>
// #include <vector>

#include "Debug.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "SmartPointers.hpp"
#include "VtkMeshReader.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"


#include "FixedRegionBoundaryCondition.hpp"
#include "HemeLBForce.hpp"
#include "MembraneDeformationForce.hpp"
#include "OutwardsPressure.hpp"
#include "OutwardsPressureWithBreaks.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:

void TestParametersOverCylinder3B() throw(Exception)
    {

        double AreaParameter[1] =  {8};//{6, 6.5, 7, 8};
        double DilationParameter[3] = {7,7.5, 8};
        double DeformationParamter[3] = { 7,7.5, 8};

        double NewEndTime = 15;
        double EndTime = 30;
        std::string output_dir = "ParameterSweep/Cylinder/";

        double P_blood = 0.002133152; double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        for (int A = 0; A < 1; A++)
        {
            for (int Di = 0; Di < 3; Di++)
            {
                for (int Def = 0; Def < 3; Def++)
                {
                    std::stringstream out;
                    out << "Param_" << AreaParameter[A] << "_DilationParam_" << DilationParameter[Di] << "_DeformationParam_" << DeformationParamter[Def];
                    std::string ParameterSet = out.str();

                    // Load and fix any settings in the simulator
                    OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(output_dir, EndTime);

                    /* Update the ouput directory for the population  */
                    static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
                    static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

                    p_simulator->RemoveAllForces();
                    p_simulator->SetEndTime(EndTime + NewEndTime);
                    p_simulator->SetSamplingTimestepMultiple(1000);
                    p_simulator->SetDt(0.001);
                    p_simulator->SetOutputDirectory(output_dir + "Parameteres2/"+ParameterSet );

                    /*
                    -----------------------------
                    Constant Compressive tissue pressure
                    ----------------------------
                    */

                    double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
                    double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

                    boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
                    p_ForceOut->SetPressure(P_blood - P_tissue);
                    p_simulator->AddForce(p_ForceOut);

                    /*
                    -----------------------------
                    Membrane forces
                    ----------------------------
                    */
                    boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
                    p_simulator->AddForce(p_shear_force);

                    /* 
                    -----------------------------
                    Update membrane properties
                    ----------------------------
                    */
                    std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
                    boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnHeteroMeshModifier<2, 3> >(*iter);
                    std::map<double, c_vector<long double, 4> > GrowthMaps;
                    GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter[A]), pow(10, -DilationParameter[Di]), pow(10, -DeformationParamter[Def]), 0);
                    //                                          Strength,hetro,stepsize, setupsolve
                    p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 0, 100, 1);

                    if (AreaParameter[A] < (double)7 || DilationParameter[Di] < (double)7 || DeformationParamter[Def] < (double)7)
                    {
                        p_Mesh_modifier->SetSlowIncreaseInMembraneStrength(1, 1);
                    }
                    p_simulator->Solve();
                    CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
                }
            }
        }
    }

};

#endif /*TESTRELAXATION_HPP_*/

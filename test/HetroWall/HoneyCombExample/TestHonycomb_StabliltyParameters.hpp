#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"

#include "Debug.hpp"
#include "SmartPointers.hpp"
// #include "VtkMeshReader.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"

#include "FixedRegionBoundaryCondition.hpp"
#include "MembraneBendingForce.hpp"
#include "MembraneBendingForce0TargetAngle.hpp"
#include "MembraneDeformationForce.hpp"
#include "OutwardsPressure.hpp"
#include "OutwardsPressureWithBreaks.hpp"

#include "RemeshingTriggerOnStepHeteroModifier.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:

   void TestMembraneParameters() throw(Exception)
    {

        double N_D[2] = { -8,-7};


                    double DilationParameter = -8;
                    double AreaParameter = -8;
                    double DeformationParamter = -7.5;

                    std::stringstream out;
                    out << "DilationParameter_" << DilationParameter << "AreaParameter" << AreaParameter << "DeformationParamter" << DeformationParamter;
                    std::string ParameterSet = out.str();
                    std::string output_dir = "DeformingHoneyComb/TestSetParameter2/" ;//+ ParameterSet;

                    TRACE("Jess is good")
                    double EndTime = 12;
                    double SamplingStep = 100;
                    double dt = 0.0001;
                    double NewEndTime = EndTime+10;

                    std::string Archieved = "DeformingHoneyComb/FlatForce3";

                    OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
                    /* Update the ouput directory for the population  */
                    static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
                    double EdgeLength = 0.0004;//(2e-6 * scale); 0.00045/2;//(2e-6 * scale);
                    static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
                    static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).ExecuteHistoryDependentRemeshing();

                    p_simulator->RemoveAllForces();
                    p_simulator->SetSamplingTimestepMultiple(SamplingStep);
                    p_simulator->SetDt(dt);
                    p_simulator->SetOutputDirectory(output_dir);

                    boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
                    p_membrane_force->SetMembraneStiffness(pow(10, -9));
                    p_simulator->AddForce(p_membrane_force);
                    /*
                    -----------------------------
                    Membrane forces
                    ----------------------------
                    */
                    boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
                    p_simulator->AddForce(p_shear_force);

                    /*
                    -----------------------------
                    Constant Compressive tissue pressure
                    ----------------------------
                    */
                    double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
                    double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

                    boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
                    p_ForceOut->SetPressure((P_blood - P_tissue));
                    p_simulator->AddForce(p_ForceOut);

                    /*
                    -----------------------------
                    Update membrane properties
                    ----------------------------
                    */
                    std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
                    boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);
                    p_Mesh_modifier->TurnOffRemeshing();
                    // p_Mesh_modifier->SetRemeshingInterval(10);
               

                    //AreaConstant           AreaDilationModulus        ShearModulus
                    std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), 1e-10) } };

                    p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
                    p_Mesh_modifier->SetmSetUpSolve(1);


                    for (int i = 0; i < 5; i++)
                    {
                        PRINT_VARIABLE(EndTime)
                        EndTime +=1;
                        p_simulator->SetEndTime(EndTime);
                        p_simulator->Solve();
                        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);

                    }
            
    }



    void offTestMembraneParameters() throw(Exception)
    {

        double N_D[2] = { -8,-7};

        for (unsigned i = 0; i < 2; i++)
        {
            for (unsigned j = 0; j < 2; j++)
            {
                for (unsigned k = 0; k < 2; k++)
                {

                    double DilationParameter = N_D[i];
                    double AreaParameter = N_D[j];
                    double DeformationParamter = N_D[k];

                    std::stringstream out;
                    out << "DilationParameter_" << DilationParameter << "AreaParameter" << AreaParameter << "DeformationParamter" << DeformationParamter;
                    std::string ParameterSet = out.str();
                    std::string output_dir = "DeformingHoneyComb/ParameterSweep2/" + ParameterSet;

                    TRACE("Jess is good")
                    double EndTime = 12;
                    double SamplingStep = 100;
                    double dt = 0.0001;
                    double NewEndTime = EndTime+0.01;

                    std::string Archieved = "DeformingHoneyComb/FlatForce3";

                    OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
                    /* Update the ouput directory for the population  */
                    static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);

                    p_simulator->RemoveAllForces();
                    p_simulator->SetSamplingTimestepMultiple(SamplingStep);
                    p_simulator->SetDt(dt);
                    p_simulator->SetOutputDirectory(output_dir);

                    boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
                    p_membrane_force->SetMembraneStiffness(pow(10, -9));
                    p_simulator->AddForce(p_membrane_force);
                    /*
                    -----------------------------
                    Membrane forces
                    ----------------------------
                    */
                    boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
                    p_simulator->AddForce(p_shear_force);

                    /*
                    -----------------------------
                    Constant Compressive tissue pressure
                    ----------------------------
                    */
                    double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
                    double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

                    boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
                    p_ForceOut->SetPressure((P_blood - P_tissue));
                    p_simulator->AddForce(p_ForceOut);

                    /*
                    -----------------------------
                    Update membrane properties
                    ----------------------------
                    */
                    std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
                    boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);
                    p_Mesh_modifier->TurnOffRemeshing();
                    // p_Mesh_modifier->SetRemeshingInterval(10);
               

                    //AreaConstant           AreaDilationModulus        ShearModulus
                    std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), 1e-10) } };

                    p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
                    p_Mesh_modifier->SetmSetUpSolve(1);

                    p_simulator->SetEndTime(NewEndTime);
                    p_simulator->Solve();
                    CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
                }
            }
        }
    }
};

#endif /*TESTRELAXATION_HPP_*/

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

#include "EnclosedRegionBoundaryCondition.hpp"
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

            double DilationParameter = -7.4;
            double AreaParameter = -7;
            double DeformationParamter = -8;

            std::stringstream out;
            out << "DilationParameter_" << DilationParameter << "AreaParameter" << AreaParameter << "DeformationParamter" << DeformationParamter;
            std::string ParameterSet = out.str();
            std::string output_dir = "DeformingHoneyComb/Bending/" ;//+ ParameterSet;

            TRACE("Jess is good")
            double EndTime = 11;
            double SamplingStep = 1;
            double dt = 0.0005;
            double NewEndTime = EndTime+10;

            std::string Archieved = "DeformingHoneyComb/FlatForce3";

            OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
            /* Update the ouput directory for the population  */
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()). SetWriteVtkAsPoints(false);
            double EdgeLength = 0.0003;//(2e-6 * scale); 0.00045/2;//(2e-6 * scale);
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
            // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).ExecuteHistoryDependentRemeshing();

            p_simulator->RemoveAllForces();
            p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            p_simulator->SetDt(dt);
            p_simulator->SetOutputDirectory(output_dir);

            boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
            p_membrane_force->SetMembraneStiffness(pow(10, -8));
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
            std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), 1e-9) } };

            p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
            p_Mesh_modifier->SetmSetUpSolve(1);

             ///////////////////////////////////////////////////////////////////////////////////////

            std::vector<c_vector<double, 3> > boundary_plane_points1;
            std::vector<c_vector<double, 3> > boundary_plane_normals1;

            std::vector<c_vector<double, 3> > boundary_plane_points2;
            std::vector<c_vector<double, 3> > boundary_plane_normals2;


            /* */
            boundary_plane_points1.push_back(Create_c_vector(0.017985684246494982, 0,0));
            boundary_plane_normals1.push_back(Create_c_vector(1,0,0));
            /// ----------------------------------------
            boundary_plane_points2.push_back(Create_c_vector(0.026139036127872058, 0,0 ));
            boundary_plane_normals2.push_back(Create_c_vector(1,0,0));

            /* */
            boundary_plane_points1.push_back(Create_c_vector(0.034963365591332625, 0,0));
            boundary_plane_normals1.push_back(Create_c_vector(1,0,0));
            /// ----------------------------------------
            boundary_plane_points2.push_back(Create_c_vector(0.04307533991933138, 0,0 ));
            boundary_plane_normals2.push_back(Create_c_vector(1,0,0));


            /* */
            boundary_plane_points1.push_back(Create_c_vector(0.05203772956749085, 0,0));
            boundary_plane_normals1.push_back(Create_c_vector(1,0,0));
            /// ----------------------------------------
            boundary_plane_points2.push_back(Create_c_vector(0.06008748970025353, 0,0 ));
            boundary_plane_normals2.push_back(Create_c_vector(1,0,0));

    
            unsigned counter =0;
            for (unsigned boundary_id = 0; boundary_id < boundary_plane_points1.size()-1; boundary_id++)
            {
                counter +=1;
                boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition(new EnclosedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , boundary_plane_points1[boundary_id], boundary_plane_normals1[boundary_id], 1)); //0.01));

                p_condition->SetPointOnPlane2( boundary_plane_points2[boundary_id]);
                p_condition->SetNormalToPlane2(boundary_plane_normals2[boundary_id]);


                p_simulator->AddCellPopulationBoundaryCondition(p_condition);
            }



            std::vector<c_vector<double, 3> > boundary_plane_points1D1;
            std::vector<c_vector<double, 3> > boundary_plane_normals1D1;

            std::vector<c_vector<double, 3> > boundary_plane_points2D1;
            std::vector<c_vector<double, 3> > boundary_plane_normals2D1;


            
            /* Diagonal 1 */
            boundary_plane_points1D1.push_back(Create_c_vector(0.012427020567408235, -0.014536134484497733 , -3.0128174549181702e-5 ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector( 0.017158353759139884,  -0.01926746767622938  , -3.0128174549181702e-5  ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));

            /* Diagonal 2 */
            boundary_plane_points1D1.push_back(Create_c_vector(0.011316836157625475, -0.0018692261147451787, 0.001297028002651899 ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.01511884615100136, -0.005671236108121062,0.001297028002651899   ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));

            /* Diagonal 3 */
            boundary_plane_points1D1.push_back(Create_c_vector(0.027451585483880808,  -0.009569147737977221, 0.0024318538097381783   ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.03145119971027805,-0.013568761964374502,  0.0024318538097381783  ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


            /* Diagonal 4 */
            boundary_plane_points1D1.push_back(Create_c_vector( 0.02708341835052748,0.005187797143975152,  0.00235468249113942   ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.031625346773351694, 0.000645868721150932, 0.00235468249113942  ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


            /* Diagonal 5 */
            boundary_plane_points1D1.push_back(Create_c_vector( 0.04497342599254667, -0.0024329349981689874, 0.003183050459649167   ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.04858572474811474,  -0.006045233753737065, 0.003183050459649167 ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


             /* Diagonal 5 */
            boundary_plane_points1D1.push_back(Create_c_vector( 0.04645169728971372, -0.014822167982229684, 0.0035170578062890282  ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector( 0.050977788945146305,-0.019348259637662277,  0.0035170578062890282 ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


             /* Diagonal 5 */
            // boundary_plane_points1D1.push_back(Create_c_vector( 0.04537342524530815, -0.0013781875058484624, 0.004500704378491893   ));
            // boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            // /// ----------------------------------------
            // boundary_plane_points2D1.push_back(Create_c_vector(0.04945107627070632, -0.005455838531246647, 0.004500704378491893  ));
            // boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


             /* Diagonal 6 */
            boundary_plane_points1D1.push_back(Create_c_vector(0.06296600825580362, -0.00829049002434537, 0.005023984536790793    ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.06665564222165832, -0.011980123990200122, 0.005023984536790793  ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


            /* Diagonal 7 */
            boundary_plane_points1D1.push_back(Create_c_vector(0.06173093960628514, 0.005548447313304878, 0.0016269652296710439     ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.06665564222165832, -0.011980123990200122, 0.005023984536790793  ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));

            /* Diagonal 8 */
            boundary_plane_points1D1.push_back(Create_c_vector(0.06291175102417185, -0.008252548152637524, -0.0003746202233619752 ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.06669267744109122, -0.012397007197412713,-0.00033079634065997506  ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


    
            // unsigned counter =0;
            // for (unsigned boundary_id = 0; boundary_id < boundary_plane_points1D1.size()-1; boundary_id++)
            // {
            //     // counter +=1;
            //     boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition(new EnclosedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , boundary_plane_points1D1[boundary_id], boundary_plane_normals1D1[boundary_id], 0.05)); //0.01));

            //     p_condition->SetPointOnPlane2( boundary_plane_points2D1[boundary_id]);
            //     p_condition->SetNormalToPlane2(boundary_plane_normals2D1[boundary_id]);


            //     p_simulator->AddCellPopulationBoundaryCondition(p_condition);
            // }


           ///////////////////////////////////////////////////////////////////////////////////////


            for (int i = 0; i < 5; i++)
            {
                PRINT_VARIABLE(EndTime)
                EndTime +=1;
                p_simulator->SetEndTime(EndTime);
                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }
            
    }



};

#endif /*TESTRELAXATION_HPP_*/

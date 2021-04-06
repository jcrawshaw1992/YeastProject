
#ifndef TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
#define TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>
#include <ctime>
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "CellsGenerator.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"

#include "OffLatticeSimulation.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

// #include "projects/VascularRemodelling/src/MembraneForces/MembraneShearForce.hpp"
// #include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"
// #include "projects/VascularRemodelling/src/MembraneForces/MembraneSurfaceForce.hpp"

#include <vector>
#include "AppliedForce.hpp"
#include "AppliedForceModifier.hpp"
#include "VtkMeshWriter.hpp"
// #include <iostream>
#include <fstream>

#include "CSVReader.hpp"
//  #include "EmptyBasementMatrix.hpp"
//  #include "LostEndothelialCell.hpp"
//  #include "HasEndothelialCell.hpp"

static const double M_TIME_FOR_SIMULATION = 0.01; //40; //50
static const double M_SAMPLING_TIME_STEP = 100; //50
static const double M_TIME_STEP = 0.002;

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

    void AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
    {
        MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

        // Calculate midpoint

        for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {

            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);

            c_vector<long double, 3> normal = zero_vector<long double>(3);

            std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
            assert(containing_elements.size() > 0);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                 iter != containing_elements.end();
                 ++iter)
            {
                // Negative as normals point inwards for these surface meshes
                c_vector<long double, 3> NewNormal = -p_cell_population->rGetMesh().GetElement(*iter)->CalculateNormal();

                NewNormal /= norm_2(NewNormal);
                normal += NewNormal;
            }
            normal /= norm_2(normal);

            c_vector<long double, 3> force = mStrength * normal; // / norm_2(cell_location);
            //   PRINT_VECTOR(force);

            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);

            cell_iter->GetCellData()->SetItem("applied_force_x", force[0]);
            cell_iter->GetCellData()->SetItem("applied_force_y", force[1]);
            cell_iter->GetCellData()->SetItem("applied_force_z", force[2]);
            cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(force));

            // cell_iter->GetCellData()->SetItem("norm_x", normal[0]);
            // cell_iter->GetCellData()->SetItem("norm_y", normal[1]);
            // cell_iter->GetCellData()->SetItem("norm_z", normal[2]);
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

class SurfaceAreaForceOnCylinder : public AbstractCellBasedTestSuite
{

private:
    double mLastStartTime;
    //double mEndTime;
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
    bool IsnumberInVector(std::vector<unsigned> Vector, unsigned Number)
    {
        bool IsInVector = 0;

        for (std::vector<unsigned>::iterator it = Vector.begin(); it != Vector.end(); ++it)
        {
            if (*it == Number)
            {
                IsInVector = 1;
                return IsInVector;
            }
        }
        return IsInVector;
    }

    std::vector<unsigned> RemoveInternalEdges(std::vector<unsigned> Neighbourhood)
    {

        sort(Neighbourhood.begin(), Neighbourhood.end());
        std::vector<unsigned> DuplicatesRemoved = Neighbourhood;
        DuplicatesRemoved.erase(unique(DuplicatesRemoved.begin(), DuplicatesRemoved.end()), DuplicatesRemoved.end());

        std::vector<unsigned> InternalPairs;
        std::set_difference(Neighbourhood.begin(), Neighbourhood.end(), DuplicatesRemoved.begin(), DuplicatesRemoved.end(),
                            std::inserter(InternalPairs, InternalPairs.begin()));

        // I now have the original vector, a vector with the duplicates removed, and a vector with only the repeating elements
        // I can use the set_difference to find the difference between the InternalPairs and the DuplicatesRemoved vectors
        std::vector<unsigned> ExternalPairs;
        std::set_difference(DuplicatesRemoved.begin(), DuplicatesRemoved.end(), InternalPairs.begin(), InternalPairs.end(),
                            std::inserter(ExternalPairs, ExternalPairs.begin()));

        return ExternalPairs;
    }

    void TestAreaFroceDragCorrectedEqui() throw(Exception)
    {

        double mesh_scale = 1;//1e-3;
    ''
        // std::string mesh_file = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/Bifurcation/SetUpData/config.vtu";
        std::string mesh_file = "/Users/jcrawshaw/docker-polnet-master/ScalledMesh/Plexus.vtu";

        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(mesh_scale, mesh_scale, mesh_scale); // so distances are in m

        std::stringstream out;
        std::string output_directory = "BifucationShrinking"; //Shrinking"; // + Parameters + "/";

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(BetaCateninOneHitCellMutationState, p_state); //Mutation to mark nodes

        // MAKE_PTR(WildTypeCellMutationState, p_WildTypeState); //Mutation to mark nodes

        //(0.023 +0.0225)/133.3223874 average pressure
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(0.002); //(M_TIME_FOR_SIMULATION);
        simulator.SetDt(0.002);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

        // Creating an object of CSVWriter
        CSVReader reader("/Users/jcrawshaw/Documents/ChasteWorkingDirectory/Bifurcation/Center.csv");
        // Get the data from CSV File
        std::vector<std::vector<double> > CenterLine = reader.getData();

        // Need to identify the bifucation region, and change the radius 

        //  Iterate over the nodes, and find the closest centerline points
        std::map<unsigned, int> NearestPoint;
        std::map<unsigned, double> RadiusFromCenter;

        std::map<unsigned, c_vector<double, 3> > NewPositions;
        std::map<unsigned, c_vector<double, 3> > LatticeNormal;
        std::map<unsigned, c_vector<double, 3> > AverageNormal;

        for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {

            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = cell_population.GetNode(node_index);

            std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
            assert(containing_elements.size() > 0);
            c_vector<long double, 3> normal = zero_vector<long double>(3);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                 iter != containing_elements.end();
                 ++iter)
            {
                // Negative as normals point inwards for these surface meshes
                normal += -cell_population.rGetMesh().GetElement(*iter)->CalculateNormal();
            }
            normal /= containing_elements.size();
            normal /= norm_2(normal);
            LatticeNormal[node_index] = normal;
        }

        for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = cell_population.GetNode(node_index);

            std::set<unsigned> NeighbouringNodeIndices = cell_population.GetNeighbouringNodeIndices(node_index);
            std::vector<unsigned> Neighbourhood;

            for (std::set<unsigned>::iterator iter = NeighbouringNodeIndices.begin();
                 iter != NeighbouringNodeIndices.end();
                 ++iter)
            {
                Neighbourhood.push_back(*iter);
                std::set<unsigned> ExtendedFamily = cell_population.GetNeighbouringNodeIndices(*iter);
                for (std::set<unsigned>::iterator iter2 = ExtendedFamily.begin();
                     iter2 != ExtendedFamily.end();
                     ++iter2)
                {

                    Neighbourhood.push_back(*iter2);
                    std::set<unsigned> Extended2Family = cell_population.GetNeighbouringNodeIndices(*iter2);
                    for (std::set<unsigned>::iterator iter3 = Extended2Family.begin();
                         iter3 != Extended2Family.end();
                         ++iter3)
                    {
                        Neighbourhood.push_back(*iter3);
                        std::set<unsigned> Extended3Family = cell_population.GetNeighbouringNodeIndices(*iter3);
                        for (std::set<unsigned>::iterator iter4 = Extended3Family.begin();
                            iter4 != Extended3Family.end();
                            ++iter4)
                        {
                            Neighbourhood.push_back(*iter4);
                        }
                    }
                }
            }
            sort(Neighbourhood.begin(), Neighbourhood.end());
            Neighbourhood.erase(unique(Neighbourhood.begin(), Neighbourhood.end()), Neighbourhood.end());

            c_vector<double, 3> AveragedLatticeNormal = Create_c_vector(0, 0, 0);

            for (std::vector<unsigned>::iterator iter = Neighbourhood.begin();
                 iter != Neighbourhood.end();
                 ++iter)
            {

                AveragedLatticeNormal += LatticeNormal[*iter];
            }

            AveragedLatticeNormal /= Neighbourhood.size();
            AveragedLatticeNormal /= norm_2(AveragedLatticeNormal);
            AverageNormal[node_index] = AveragedLatticeNormal;

            c_vector<double, 3> Location = p_node->rGetLocation();
            double MinLength = 100;
            int ClosestPoint = 0;
            double Radius;
            double counter = 0;
            for (std::vector<std::vector<double> >::iterator iter = CenterLine.begin();
                 iter != CenterLine.end();
                 ++iter)
            {
                c_vector<double, 3> Point;

                Point[0] = (*iter)[1] * mesh_scale;
                Point[1] = (*iter)[2] * mesh_scale;
                Point[2] = (*iter)[3] * mesh_scale;
                // PRINT_VARIABLE(CenterLine[iter][1]);

                c_vector<double, 3> Distance = Location - Point;
                // double value = std::strtod(iter->c_str(), NULL);
                // PRINT_VARIABLE( value );
                double Angle = acos(inner_prod(Distance, AveragedLatticeNormal) / norm_2(Distance));
                double Length = norm_2(Distance);
                // PRINT_VARIABLE((*iter)[0] );
                if (Length < MinLength && (*iter)[0] > 0)
                {
                    MinLength = Length;
                    ClosestPoint = counter;
                    Radius = (*iter)[0] * mesh_scale;
                }
                counter += 1;
            }

            // PRINT_VARIABLE(ClosestPoint);
            NearestPoint[node_index] = ClosestPoint;
            RadiusFromCenter[node_index] = Radius;

            std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
            assert(containing_elements.size() > 0);
            c_vector<long double, 3> normal = zero_vector<long double>(3);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                 iter != containing_elements.end();
                 ++iter)
            {
                // Negative as normals point inwards for these surface meshes
                normal += -cell_population.rGetMesh().GetElement(*iter)->CalculateNormal();
            }
            normal /= containing_elements.size();
            normal /= norm_2(normal);
            LatticeNormal[node_index] = normal;
            // PRINT_VARIABLE(Radius);
            // if (Radius>  0.8 * mesh_scale)
            // {

            NewPositions[node_index] = Location - Radius * 0.7 * AveragedLatticeNormal;
            // }
            // else{

            // NewPositions[node_index] = Location - Radius *0.2* AveragedLatticeNormal ;

            // }
            c_vector<double, 3> Difference = NewPositions[node_index] - Location;
            // PRINT_VECTOR( Difference);
        }

        for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
            // Node<3>* p_node = cell_population.GetNode(node_index);
            ChastePoint<3> NewLocation; // Needs to be a chaste point im setting
            NewLocation.rGetLocation()[0] = NewPositions[node_index][0];
            NewLocation.rGetLocation()[1] = NewPositions[node_index][1];
            NewLocation.rGetLocation()[2] = NewPositions[node_index][2];

            mesh.SetNode(node_index, NewLocation, false);

            (cell_iter)->GetCellData()->SetItem("Initial_Location_X", NewPositions[node_index][0]);
            (cell_iter)->GetCellData()->SetItem("Initial_Location_Y", NewPositions[node_index][1]);
            (cell_iter)->GetCellData()->SetItem("Initial_Location_Z", NewPositions[node_index][2]);
        }

        // //  -----------------------------
        // //  Save original mesh before I shrink it
        // //----------------------------

        // VtkMeshWriter<2, 3> mesh_writer(output_directory, "Original", false);
        // MutableMesh<2, 3>* p_mesh = &(dynamic_cast<MeshBasedCellPopulation<2, 3>*>(&(simulator.rGetCellPopulation()))->rGetMesh());
        // p_mesh->Scale(1.0 / mesh_scale, 1.0 / mesh_scale, 1.0 / mesh_scale); // so distances are back in original scal
        // mesh_writer.WriteFilesUsingMesh(*p_mesh);

        // //  -----------------------------
        // //  Parameters
        // //----------------------------

        // double ElasticShearModulus = 4.4e-9;
        // double AreaDilationModulus = 0.9e-11;
        // double membrane_constant = 5;
        // double Area_constant = 0.9e-15;
        // double Pressure = 0.001;

        // //  -----------------------------
        // //  Constant pressure force
        // //----------------------------

        // std::map<unsigned, c_vector<double, 3> > NewPositions;
        // std::map<unsigned, c_vector<double, 3> > LatticeNormal;
        // std::map<unsigned, c_vector<double, 3> > NewLatticeNormal;
        // std::map<unsigned, double> InitalAverageRadius;

        // std::map<unsigned, double> AverageRadius;
        // std::map<unsigned, double> LocalAverageRadius;

        // for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
        //      cell_iter != cell_population.End();
        //      ++cell_iter)
        // {
        //     unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
        //     Node<3>* p_node = cell_population.GetNode(node_index);

        //     c_vector<double, 3> NodeLocation = p_node->rGetLocation();

        //     c_vector<long double, 3> normal = zero_vector<long double>(3);
        //     double NumberOfNormals = 0;

        //     std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
        //     assert(containing_elements.size() > 0);
        //     for (std::set<unsigned>::iterator iter = containing_elements.begin();
        //          iter != containing_elements.end();
        //          ++iter)
        //     {
        //         // Negative as normals point inwards for these surface meshes
        //         normal += -cell_population.rGetMesh().GetElement(*iter)->CalculateNormal();
        //         NumberOfNormals += 1;
        //     }
        //     normal /= NumberOfNormals;
        //     normal /= norm_2(normal);

        //     LatticeNormal[node_index] = normal;
        // }

        // // Need to find R, if it is reasonable. Use the same stuff as what I used for the curviture, but rather than the average curvature
        // double RunningAverage = 0;
        // double CurrentNumberOfNodes = 0;
        // // Want average r for each node

        // for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
        //      cell_iter != cell_population.End();
        //      ++cell_iter)
        // {
        //     unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
        //     Node<3>* p_node = cell_population.GetNode(node_index);

        //     unsigned NumberOfNeighbours = 0;
        //     double AverageR = 0;
        //     double InitalAverageR = 0;

        //     CurrentNumberOfNodes += 1;

        //     std::set<unsigned> NeighbouringNodeIndices = cell_population.GetNeighbouringNodeIndices(node_index);

        //     std::vector<unsigned> Neighbourhood;

        //     for (std::set<unsigned>::iterator iter = NeighbouringNodeIndices.begin();
        //          iter != NeighbouringNodeIndices.end();
        //          ++iter)
        //     {
        //         Neighbourhood.push_back(*iter);
        //         std::set<unsigned> ExtendedFamily = cell_population.GetNeighbouringNodeIndices(*iter);
        //         for (std::set<unsigned>::iterator iter2 = ExtendedFamily.begin();
        //              iter2 != ExtendedFamily.end();
        //              ++iter2)
        //         {
        //             if (*iter2 != node_index)
        //             {
        //                 Neighbourhood.push_back(*iter2);
        //                 std::set<unsigned> Extended2Family = cell_population.GetNeighbouringNodeIndices(*iter2);
        //                 for (std::set<unsigned>::iterator iter3 = Extended2Family.begin();
        //                      iter3 != Extended2Family.end();
        //                      ++iter3)
        //                 {
        //                     if (*iter3 != node_index)
        //                     {
        //                         Neighbourhood.push_back(*iter3);
        //                     }
        //                 }
        //             }
        //         }
        //     }

        //     sort(Neighbourhood.begin(), Neighbourhood.end());
        //     Neighbourhood.erase(unique(Neighbourhood.begin(), Neighbourhood.end()), Neighbourhood.end());

        //     // PRINT_VECTOR(Neighbourhood);
        //     NumberOfNeighbours = 0;
        //     // std::set<unsigned> NeighbouringNodeIndices  = cell_population.GetNeighbouringNsnode_index);
        //     // for (std::set<unsigned>::iterator iter = NeighbouringNodeIndices.begin();
        //     //      iter != NeighbouringNodeIndices.end();
        //     //      ++iter)
        //     // {
        //     for (std::vector<unsigned>::iterator iter = Neighbourhood.begin();
        //          iter != Neighbourhood.end();
        //          ++iter)
        //     {
        //         Node<3>* pNeighbourNode = cell_population.GetNode(*iter);
        //         //  Node<3>* pNeighbourNode = * iter;
        //         c_vector<double, 3> Vector01 = pNeighbourNode->rGetLocation() - p_node->rGetLocation();

        //         double angle = acos(inner_prod(LatticeNormal[node_index], LatticeNormal[*iter]));

        //         double R = norm_2(Vector01) / (2 * sin(angle / 2));
        //         //    PRINT_VARIABLE(R);
        //         if (std::abs(angle) > 0.02) // edge is too small and the angle between the normals of the two nodes are nearly parallell
        //         {
        //             InitalAverageR += R;
        //             NumberOfNeighbours += 1;
        //         }
        //     }
        //     InitalAverageR /= NumberOfNeighbours;
        //     NumberOfNeighbours = 0;

        //     // for (std::set<unsigned>::iterator iter = NeighbouringNodeIndices.begin();
        //     //         iter != NeighbouringNodeIndices.end();
        //     //         ++iter)
        //     //     {
        //     for (std::vector<unsigned>::iterator iter = Neighbourhood.begin();
        //          iter != Neighbourhood.end();
        //          ++iter)
        //     {

        //         //  unsigned Neighbour_node_index = cell_population.GetLocationIndexUsingCell(* iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
        //         Node<3>* pNeighbourNode = cell_population.GetNode(*iter);
        //         c_vector<double, 3> Vector01 = pNeighbourNode->rGetLocation() - p_node->rGetLocation();
        //         double angle = acos(inner_prod(LatticeNormal[node_index], LatticeNormal[*iter]));
        //         double R = norm_2(Vector01) / (2 * sin(angle / 2));
        //         //  PRINT_3_VARIABLES(angle, R, norm_2(Vector01));

        //         if (angle > 0.02 && R < 1.5 * InitalAverageR) // edge is too small and the angle between the normals of the two nodes are nearly parallell
        //         {

        //             NumberOfNeighbours += 1;
        //             AverageR += R;
        //         }
        //         else
        //         {
        //             TRACE("Catch 1");
        //         }
        //     }
        //     AverageRadius[node_index] = AverageR / NumberOfNeighbours;
        //     PRINT_VARIABLE(AverageRadius[node_index]);
        // }

        // double TotalAverageRadius = 0;
        // double MaxRadius = 0;
        // double MinRadius = 0;
        // unsigned MinNode;
        // unsigned MaxNode;
        // std::vector<unsigned> MinNodeNeighbourhood ;

        // std::map<unsigned, std::vector<unsigned> > ExtendedNeighbood;
        // std::map<unsigned, std::vector<unsigned> > DistantNeighbood;
        // for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
        //      cell_iter != cell_population.End();
        //      ++cell_iter)
        // {
        //     std::vector<unsigned> ExtendedNeighbours;
        //     std::vector<unsigned> DistantNeighbours;

        //     unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
        //     Node<3>* p_node = cell_population.GetNode(node_index);

        //     unsigned NumberOfNeighbours = 0;
        //     double LocalAverage = 0;
        //     //  LocalAverage += AverageRadius[node_index ];
        //     std::vector<unsigned> LocalNeighbours;
        //     // std::set<unsigned> NeighbouringNodeIndices  = GetNeighbouringNodeIndices(node_index);
        //     std::set<unsigned> NeighbouringNodeIndices = cell_population.GetNeighbouringNodeIndices(node_index);
        //     for (std::set<unsigned>::iterator iter = NeighbouringNodeIndices.begin();
        //          iter != NeighbouringNodeIndices.end();
        //          ++iter)
        //     {
        //         LocalNeighbours.push_back(*iter);
        //         ExtendedNeighbours.push_back(*iter);
        //         std::set<unsigned> ExtendedFamily = cell_population.GetNeighbouringNodeIndices(*iter);
        //         for (std::set<unsigned>::iterator iter2 = ExtendedFamily.begin();
        //              iter2 != ExtendedFamily.end();
        //              ++iter2)
        //         {
        //             // if (IsnumberInVector(ExtendedNeighbours, * iter2) ==0)
        //             // {

        //             ExtendedNeighbours.push_back(*iter2);
        //             std::set<unsigned> ExtendedExtendedFamily = cell_population.GetNeighbouringNodeIndices(*iter2);
        //             for (std::set<unsigned>::iterator iter3 = ExtendedExtendedFamily.begin();
        //                  iter3 != ExtendedExtendedFamily.end();
        //                  ++iter3)
        //             {
        //                 // if (IsnumberInVector(ExtendedNeighbours, * iter3) ==0)
        //                 // {
        //                 ExtendedNeighbours.push_back(*iter3);
        //                 // DistantNeighbours.push_back(* iter3);
        //                 std::set<unsigned> ExtendedExtendedExtendedFamily = cell_population.GetNeighbouringNodeIndices(*iter3);
        //                 for (std::set<unsigned>::iterator iter4 = ExtendedExtendedExtendedFamily.begin();
        //                      iter4 != ExtendedExtendedExtendedFamily.end();
        //                      ++iter4)
        //                 {
        //                     //   if (IsnumberInVector(ExtendedNeighbours, * iter4) ==0)
        //                     //   {
        //                     ExtendedNeighbours.push_back(*iter4);
        //                     // DistantNeighbours.push_back(* iter4);
        //                     std::set<unsigned> Extended4Family = cell_population.GetNeighbouringNodeIndices(*iter4);
        //                     for (std::set<unsigned>::iterator iter5 = Extended4Family.begin();
        //                          iter5 != Extended4Family.end();
        //                          ++iter5)
        //                     {
        //                         ExtendedNeighbours.push_back(*iter5);
        //                         DistantNeighbours.push_back(* iter5);
        //                         std::set<unsigned> Extended5Family  = cell_population.GetNeighbouringNodeIndices(* iter5);
        //                         // for (std::set<unsigned>::iterator iter6 = Extended5Family.begin();
        //                         //         iter6 != Extended5Family.end();
        //                         //         ++iter6)
        //                         //             {
        //                         //                 ExtendedNeighbours.push_back(* iter6);
        //                         //                 // DistantNeighbours.push_back(* iter6);
        //                         //                 std::set<unsigned> Extended6Family  = cell_population.GetNeighbouringNodeIndices(* iter5);
        //                         //                 // for (std::set<unsigned>::iterator iter7 = Extended6Family.begin();
        //                         //                 //         iter7 != Extended6Family.end();
        //                         //                 //         ++iter7)
        //                         //                 //             {
        //                         //                 //             ExtendedNeighbours.push_back(* iter7);
        //                         //                 //             }
        //                         //              }

        //                     }
        //                 }
        //             }
        //         }
        //     }

        //     sort(ExtendedNeighbours.begin(), ExtendedNeighbours.end());
        //     ExtendedNeighbours.erase(unique(ExtendedNeighbours.begin(), ExtendedNeighbours.end()), ExtendedNeighbours.end());

        //     std::vector<unsigned> Neighbourhood = ExtendedNeighbours;

        //     for (std::set<unsigned>::iterator iter = NeighbouringNodeIndices.begin();
        //          iter != NeighbouringNodeIndices.end();
        //          ++iter)
        //     {
        //         Neighbourhood.push_back(*iter);
        //         std::set<unsigned> ExtendedFamily = cell_population.GetNeighbouringNodeIndices(*iter);
        //         for (std::set<unsigned>::iterator iter2 = ExtendedFamily.begin();
        //              iter2 != ExtendedFamily.end();
        //              ++iter2)
        //         {
        //             Neighbourhood.push_back(*iter2);
        //             std::set<unsigned> Extended2Family = cell_population.GetNeighbouringNodeIndices(*iter2);
        //             for (std::set<unsigned>::iterator iter3 = Extended2Family.begin();
        //                  iter3 != Extended2Family.end();
        //                  ++iter3)
        //             {
        //                 Neighbourhood.push_back(*iter3);
        //             }
        //         }
        //     }

        //     DistantNeighbours = RemoveInternalEdges(Neighbourhood);

        //     ExtendedNeighbood[node_index] = ExtendedNeighbours;
        //     DistantNeighbood[node_index] = DistantNeighbours;
        //     double InitalLocalAverage = 0;
        //     c_vector<double, 3> AveragedLatticeNormal =Create_c_vector(0,0,0);

        //     for (std::vector<unsigned>::iterator iter = ExtendedNeighbours.begin();
        //          iter != ExtendedNeighbours.end();
        //          ++iter)
        //     {
        //         CellPtr p_neighbour_cell = cell_population.GetCellUsingLocationIndex(*iter);
        //         //   double edge = (p_neighbour_cell)->GetCellData()->GetItem("Edge");

        //         NumberOfNeighbours += 1;
        //         InitalLocalAverage += AverageRadius[*iter];
        //     }
        //      NumberOfNeighbours = 0;
        //     for (std::vector<unsigned>::iterator iter = ExtendedNeighbours.begin();
        //          iter != ExtendedNeighbours.end();
        //          ++iter)
        //     {
        //         CellPtr p_neighbour_cell = cell_population.GetCellUsingLocationIndex(*iter);
        //         //   double edge = (p_neighbour_cell)->GetCellData()->GetItem("Edge");
        //         if (AverageRadius[*iter] < 3 * InitalLocalAverage)
        //         {
        //             NumberOfNeighbours += 1;
        //             LocalAverage += AverageRadius[*iter];

        //         }
        //          AveragedLatticeNormal += LatticeNormal[*iter];

        //     }
        //     LocalAverageRadius[node_index] = LocalAverage / NumberOfNeighbours;
        //     AveragedLatticeNormal /=ExtendedNeighbours.size() ;
        //     NewLatticeNormal[node_index] = AveragedLatticeNormal/ norm_2(AveragedLatticeNormal);

        //     if (LocalAverageRadius[node_index] > MaxRadius)
        //     {
        //         MaxRadius = LocalAverageRadius[node_index];
        //         MaxNode = node_index;
        //     }
        //     else if (LocalAverageRadius[node_index] < MinRadius)
        //     {
        //         MaxRadius = LocalAverageRadius[node_index];
        //         MinNode = node_index;
        //         MinNodeNeighbourhood = LocalNeighbours;
        //     }

        //     TotalAverageRadius += LocalAverageRadius[node_index];
        // }

        // TotalAverageRadius = TotalAverageRadius - (MaxRadius + MinRadius);
        // TotalAverageRadius = TotalAverageRadius / (mesh.GetNumNodes() - 2);

        // double Scalling = 2;

        // for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
        //      cell_iter != cell_population.End();
        //      ++cell_iter)
        // {
        //     unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
        //     Node<3>* p_node = cell_population.GetNode(node_index);

        //     c_vector<double, 3> NodeLocation = p_node->rGetLocation();

        //     // Need to walk backward into the mesh by the scalling factor

        //     // if (std::abs(LocalAverageRadius[node_index] - MaxRadius) < 0.001 || std::abs(LocalAverageRadius[node_index] - MinRadius) < 0.001)
        //     // {
        //     //     LocalAverageRadius[node_index] = TotalAverageRadius;
        //     // }

        // // if (IsnumberInVector(MinNodeNeighbourhood , node_index ) ==1)
        // // {
        // //     LatticeNormal[node_index] = LatticeNormal[MinNode];
        // // }

        //     c_vector<double, 3> PositionVector = NodeLocation - (LocalAverageRadius[node_index] * 0.3) * NewLatticeNormal[node_index];

        //     NewPositions[node_index] = PositionVector;
        //     // LatticeNormal[node_index] = normal ;

        //     ChastePoint<3> NewLocation; // Needs to be a chaste point im setting
        //     NewLocation.rGetLocation()[0] = NewPositions[node_index][0];
        //     NewLocation.rGetLocation()[1] = NewPositions[node_index][1];
        //     NewLocation.rGetLocation()[2] = NewPositions[node_index][2];

        //     mesh.SetNode(node_index, NewLocation, false);

        //     (cell_iter)->GetCellData()->SetItem("Initial_Location_X", NewPositions[node_index][0]);
        //     (cell_iter)->GetCellData()->SetItem("Initial_Location_Y", NewPositions[node_index][1]);
        //     (cell_iter)->GetCellData()->SetItem("Initial_Location_Z", NewPositions[node_index][2]);

        //     (cell_iter)->GetCellData()->SetItem("Original_Location_X", NodeLocation[0]);
        //     (cell_iter)->GetCellData()->SetItem("Original_Location_Y", NodeLocation[1]);
        //     (cell_iter)->GetCellData()->SetItem("Original_Location_Z", NodeLocation[2]);
        //     TRACE("UPDATED");
        // }

        //-----------------------

        // If node is higher than all of its second neighbours

        // for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
        //      cell_iter != cell_population.End();
        //      ++cell_iter)
        // {
        //     unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
        //     Node<3>* p_node = cell_population.GetNode(node_index);

        // }

        // // //  -----------------------------
        // // //  Shearing Force
        // //----------------------------

        // // boost::shared_ptr<MembraneShearForce> p_shear_force(new MembraneShearForce());
        // // p_shear_force->SetScallingShear(Scalling);
        // // p_shear_force->SetAreaDilationModulus(AreaDilationModulus);
        // // p_shear_force->SetElasticShearModulus(ElasticShearModulus);
        // // p_shear_force->SetupMembraneConfiguration(cell_population);

        // // simulator.AddForce(p_shear_force);

        // //  -----------------------------
        // //  Area Force
        // //  ----------------------------

        // // double Area_constant = 5 * 1e-5;

        // // boost::shared_ptr<MembraneSurfaceForce> p_surface_force(new MembraneSurfaceForce());
        // // p_surface_force->SetScallingArea(Scalling);
        // // p_surface_force->SetupInitialAreas(cell_population);
        // // p_surface_force->SetMembraneStiffness(Area_constant);

        // // simulator.AddForce(p_surface_force);

        // // // //   -----------------------------
        // // // //   Bending Force
        // // // //  ----------------------------

        // // boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        // // p_membrane_force->SetScallingBending(Scalling);
        // // p_membrane_force->SetupInitialMembrane(mesh, cell_population);
        // // p_membrane_force->SetMembraneStiffness(membrane_constant);

        // // simulator.AddForce(p_membrane_force);

        // // //         /*
        // // //         -----------------------------
        // // //         Boundaries
        // // //         ----------------------------
        // //         */

        // //Create a plane boundary to represent the inlet and pass them to the simulation

        // c_vector<long double, 3> Boundary1 = Create_c_vector(12.48824184355049, 34.75302061558864, 41.78949821113195);
        // c_vector<long double, 3> Normal1 = Create_c_vector(-0.05607774749413225, 0.762765339938692, 0.6442393362906335);
        // double Radius1 = 1.7;

        // c_vector<long double, 3> Boundary2 = Create_c_vector(12.597373380655702, 48.382440094438316, 42.984851419357064);
        // c_vector<long double, 3> Normal2 = Create_c_vector(-0.04847413662454751, -0.989768366942236, -0.13419701143842153);
        // double Radius2 = 1.2;

        // c_vector<long double, 3> Boundary3 = Create_c_vector(-0.04710924565104161, -0.9898343148598956, -0.13419667693363677);
        // c_vector<long double, 3> Normal3 = Create_c_vector(-0.04710924565104161, -0.9898343148598956, -0.13419667693363677);
        // double Radius3 = 1.3;

        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, Radius1));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, Radius2));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_2);
        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_3(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary3, Normal3, Radius3));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_3);

        simulator.Solve();

        // To reset before looping: this is usually done by the SetUp and TearDown methods
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/

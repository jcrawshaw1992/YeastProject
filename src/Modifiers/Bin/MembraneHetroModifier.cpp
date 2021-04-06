
/*

This code controles the membrame properties. 

It sets the regions on one branch to have a greater stiffness than the rest of a birfucation

*/

#include "MembraneHetroModifier.hpp"
#include <algorithm>
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

#include <math.h>
#include "Debug.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MembraneHetroModifier<ELEMENT_DIM, SPACE_DIM>::MembraneHetroModifier()
        : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MembraneHetroModifier<ELEMENT_DIM, SPACE_DIM>::~MembraneHetroModifier()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembraneHetroModifier<ELEMENT_DIM, SPACE_DIM>::SetMembranePropeties(std::map<double, c_vector<long double, 4> > GrowthMaps, double Strength, bool Hetrogeneous)
{
    mGrowthMaps = GrowthMaps;
    TRACE("New membrane properties")
    PRINT_2_VARIABLES(mGrowthMaps[10][3], GrowthMaps[10][3]);

    mOn = 1;
    mAchievedTargetK = 0;
    mStrength = Strength;
    mHetro = Hetrogeneous;
    mCounter =50;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembraneHetroModifier<ELEMENT_DIM, SPACE_DIM>::SetBoundaryNode(unsigned Node_Index)
{
    // I need location and neighbours here
    mStoredBoundaryDefiningNodes.push_back(Node_Index);
    TRACE("node stored")
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembraneHetroModifier<ELEMENT_DIM, SPACE_DIM>::SetBoundary( c_vector<long double, 3>  Point, c_vector<long double, 3>  normal)
{
    
    mBoundaryPoint = Point;
    mBoundaryNormal = normal;
}





template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembraneHetroModifier<ELEMENT_DIM, SPACE_DIM>::SetAConstantHetrogenaity(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
 assert(SPACE_DIM == 3); assert(ELEMENT_DIM == 2);
 MAKE_PTR(EmptyBasementMatrix, p_Basement); MAKE_PTR(HasEndothelialCell, p_EC); MAKE_PTR(LostEndothelialCell, p_NextToBasement); //Mutations to mark 


    if(mHetro)
    {
        std::set<unsigned> MutantNodeIndices; std::set<unsigned> EdgeMutantNodeIndices;
    
        // if above plane with a normal 
        c_vector<double, 3> BasementNormal  = Create_c_vector(0.7747460075692358,0.6253214339470684,0.09349720852470317);
        // with an orgin or 
        c_vector<double, 3> BasementPoint =  Create_c_vector(0.10565108560001966,0.15455588216280267,0.017040208357096718);
        // then it is a basement matrix

        // if above plane with a normal and a point of ..., then its an edge
        c_vector<double, 3> EdgeNormal  = Create_c_vector(-0.7752270779441546, -0.6178870924813352,-0.131295539022199 );
        c_vector<double, 3> EdgePoint =  Create_c_vector(0.10660370304382119,0.18311312173972846, -0.0019649382624181404);  

        mDistanceBetweenPlanes = std::abs(norm_2(BasementPoint - EdgePoint ));
        PRINT_VARIABLE(norm_2(BasementPoint - EdgePoint ));
        double Ks = mGrowthMaps[mStrength](2);  double Ka = mGrowthMaps[mStrength](1); double KA = mGrowthMaps[mStrength](0);  double d = mDistanceBetweenPlanes;
        double Kbs =  (double)mGrowthMaps[1.1](2); double Kba =  (double)mGrowthMaps[1.1](1); double KbA =  (double)mGrowthMaps[1.1](0);   
        

        c_vector<double, 3> X;
        for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("Boundary", 0);
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            X = rCellPopulation.GetNode(node_index)->rGetLocation();
        
            if (BasementNormal[0]*(X[0]-BasementPoint[0]) + BasementNormal[1]*(X[1]-BasementPoint[1])  + BasementNormal[2]*(X[2]-BasementPoint[2]) < 0)
            {
                cell_iter->SetMutationState(p_Basement);
                MutantNodeIndices.insert(node_index);
                cell_iter->GetCellData()->SetItem("NodeIndex",  node_index);
                cell_iter->GetCellData()->SetItem("ShearModulus", Kbs);
                cell_iter->GetCellData()->SetItem("AreaDilationModulus", Kba);
                cell_iter->GetCellData()->SetItem("AreaConstant", KbA);
                cell_iter->GetCellData()->SetItem("BendingConstant", mGrowthMaps[1.2](3));
            }
            else
            {
                cell_iter->SetMutationState(p_EC);
                cell_iter->GetCellData()->SetItem("NodeIndex",  node_index  );
                cell_iter->GetCellData()->SetItem("ShearModulus", mGrowthMaps[mStrength](2));
                cell_iter->GetCellData()->SetItem("AreaDilationModulus", mGrowthMaps[mStrength](1));
                cell_iter->GetCellData()->SetItem("AreaConstant", mGrowthMaps[mStrength](0));
                cell_iter->GetCellData()->SetItem("BendingConstant", mGrowthMaps[mStrength](3));
            }
        
        }
                 
            // The endothelial Z mMinZ;
            //Shear K
            double c1S = Kbs;  double a1S = (Ks - c1S) /(d*d);
            //Alpha 
            double c1a = Kba;  double a1a = (Ka - c1a) / (d *d);//   
            //Area  
            double c1A = KbA; double a1A = (KA - c1A) / (d* d);//     
            double ShearMod; double AlphaMod; double AreaMod;

            double A = BasementNormal[0]*BasementPoint[0]; double B = BasementNormal[1]*BasementPoint[1]; double C = BasementNormal[2]*BasementPoint[2];

        for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {   
            //  if cell type cell_iter->SetMutationState(p_EC);  and above the plane, then set at edge mutation
            if (cell_iter->GetMutationState()->template IsType<HasEndothelialCell>())
            {
                unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                X = rCellPopulation.GetNode(node_index)->rGetLocation();
                if (EdgeNormal[0]*(X[0]-EdgePoint[0]) + EdgeNormal[1]*(X[1]-EdgePoint[1])  + EdgeNormal[2]*(X[2]-EdgePoint[2]) > 0)
                {
                    cell_iter->SetMutationState(p_NextToBasement);
                    MutantNodeIndices.insert(node_index);
                    cell_iter->GetCellData()->SetItem("NodeIndex",  node_index);
                    //DistanceFromBasement  
                    // double D =  std::abs(A*X[0]+B*X[1]+C*X[2])/sqrt(A*A+B*B+C*C);

                    c_vector<double, 3> DistanceFromCentroid =  X -  BasementPoint;
                    // DitanceToPlane
                    double D = inner_prod(DistanceFromCentroid, BasementNormal);
                    // PRINT_VARIABLE(D);
                    //K =a(Z)^2 +c -- distance to plane
                    ShearMod = a1S * D*D + c1S;  AlphaMod = a1a * D*D + c1a ;  AreaMod = a1A * D*D + c1A;
        
                    cell_iter->GetCellData()->SetItem("ShearModulus", ShearMod);  cell_iter->GetCellData()->SetItem("AreaDilationModulus", AlphaMod); cell_iter->GetCellData()->SetItem("AreaConstant", AreaMod);
                }
            }  
        }
    }
    else
    {
       
        double Ks = mGrowthMaps[mStrength](2);  double Ka = mGrowthMaps[mStrength](1); double KA = mGrowthMaps[mStrength](0);  double d = mDistanceBetweenPlanes;
        

        c_vector<double, 3> X;
        for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            cell_iter->SetMutationState(p_EC);
            cell_iter->GetCellData()->SetItem("NodeIndex",  node_index  );
            cell_iter->GetCellData()->SetItem("ShearModulus", mGrowthMaps[mStrength](2));
            cell_iter->GetCellData()->SetItem("AreaDilationModulus", mGrowthMaps[mStrength](1));
            cell_iter->GetCellData()->SetItem("AreaConstant", mGrowthMaps[mStrength](0));
            cell_iter->GetCellData()->SetItem("BendingConstant", mGrowthMaps[mStrength](3));
        }
    }
}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembraneHetroModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{

 // PRINT_2_VARIABLES(ELEMENT_DIM, SPACE_DIM);
    assert(SPACE_DIM == 3);
    assert(ELEMENT_DIM == 2);
    MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
    MAKE_PTR(EmptyBasementMatrix, p_Basement); //Mutation to mark nodes one basement
    MAKE_PTR(HasEndothelialCell, p_EC); //Mutation to mark nodes with ECs
    MAKE_PTR(LostEndothelialCell, p_NextToBasement); //Mutation to mark nodes with ECs next to the basement
    c_vector<long double, 3> Node_location;
    c_vector<long double, 3> P;
 
    std::set<unsigned> MutantNodeIndices;
    std::set<unsigned> EdgeMutantNodeIndices;

    double Distance = 0.08;
    double A = mBoundaryNormal[0]*mBoundaryPoint[0]; double B = mBoundaryNormal[1]*mBoundaryPoint[1];  double C = mBoundaryNormal[2]*mBoundaryPoint[2];
    double DistanceFromPlane;

    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        P = rCellPopulation.GetNode(node_index)->rGetLocation();
        DistanceFromPlane = std::abs(A*P[0] + B*P[1] + C*P[2])/sqrt(A*A+B*B+C*C);
        
        double DistanceToPoint = std::abs(norm_2(P - mBoundaryPoint));
        if (DistanceFromPlane < Distance)   
        {
            cell_iter->SetMutationState(p_Basement);
            MutantNodeIndices.insert(node_index);
            mSamplebasementNode = node_index;
          
            // TRACE("Mutant")
        }
        else
        {
            cell_iter->SetMutationState(p_EC);
            mSampleECNode = node_index;
            
        }
        cell_iter->GetCellData()->SetItem("NodeIndex",  node_index  );
        cell_iter->GetCellData()->SetItem("ShearModulus", mGrowthMaps[mStrength](2));
        cell_iter->GetCellData()->SetItem("AreaDilationModulus", mGrowthMaps[mStrength](1));
        cell_iter->GetCellData()->SetItem("AreaConstant", mGrowthMaps[mStrength](0));
        cell_iter->GetCellData()->SetItem("BendingConstant", mGrowthMaps[mStrength](3));

    }

    
    if (mHetro)
    {
    
        // Loop over the mutants and get the nodes along the boundary
        for (std::set<unsigned>::iterator iter = MutantNodeIndices.begin();
             iter != MutantNodeIndices.end();
             ++iter)
        {


            // Need to check if any of the neighbours are an EC, if so
            // Save the EC as a mutant on the boundary cell_iter->SetMutationState(p_NextToBasement);
            mBasementNodes.push_back(*iter);
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
            std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringLocationIndices(p_cell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
            for (std::set<unsigned>::iterator N_iter = neighbouring_node_indices.begin();
                 N_iter != neighbouring_node_indices.end();
                 ++N_iter)
            {

                CellPtr N_cell = rCellPopulation.GetCellUsingLocationIndex(*N_iter);
                if (N_cell->GetMutationState()->template IsType<HasEndothelialCell>())
                {
                    // This cell is type EC and next to the basement matrix, which means the
                    // this basement mutant is on the edge :)
                    EdgeMutantNodeIndices.insert(*iter);
                    break;
                }
            }
        }



        // SO now we have the membrane region and the lining edge of the membrane region -- can clean it up a bit


        // Find the centroid of the edges
        c_vector<double, 3> EdgeCentroid= Create_c_vector(0,0,0);

        for (std::set<unsigned>::iterator iter = EdgeMutantNodeIndices.begin();
             iter != EdgeMutantNodeIndices.end();
             ++iter)
        {
                EdgeCentroid += rCellPopulation.GetNode(*iter)->rGetLocation();
        }   
        EdgeCentroid /= EdgeMutantNodeIndices.size();

        // Now loop over and find the two minimum points

        double Min1 = 1000;
        double Min2 = 1000000;
        unsigned NodeA;
        unsigned NodeB;

         for (std::set<unsigned>::iterator iter = EdgeMutantNodeIndices.begin();
             iter != EdgeMutantNodeIndices.end();
             ++iter)
        {
               double distance = norm_2(EdgeCentroid - rCellPopulation.GetNode(*iter)->rGetLocation() );
                   if (distance < Min1)
                   {    
                        Min1 = distance;
                        NodeA = *iter;
                   } 
        }  
        c_vector<double, 3> VectorA = EdgeCentroid  -rCellPopulation.GetNode(NodeA)->rGetLocation(); 

          for (std::set<unsigned>::iterator iter = EdgeMutantNodeIndices.begin();
             iter != EdgeMutantNodeIndices.end();
             ++iter)
        {   
                c_vector<double, 3> VectorB = EdgeCentroid  - rCellPopulation.GetNode(*iter)->rGetLocation(); 
                double Angle = std::abs(acos(inner_prod(VectorA,VectorB )/(norm_2(VectorA)*norm_2(VectorB)) )) ;
                if (Angle > 3*M_PI/4)
                {
                    double distance = norm_2(VectorB);
                   if (distance < Min2)
                   {    
                        Min2 = distance;
                        NodeB = *iter;
                   } 

                }
               
        }  

        c_vector<double, 3> AX = EdgeCentroid - rCellPopulation.GetNode(NodeA)->rGetLocation();
        c_vector<double, 3> BX = EdgeCentroid - rCellPopulation.GetNode(NodeB)->rGetLocation();
        c_vector<double, 3> n = VectorProduct(AX,BX);
        // Equation of plane 
        // Plane = n[0](x-EdgeCentroid[0]) + n[1](x-EdgeCentroid[1])  +n[2](x-EdgeCentroid[2])  =0

        // mEdgeOfBasementPlane
        mEdgeOfBasementPlane[0] = n[0];  mEdgeOfBasementPlane[1] =EdgeCentroid[0];  mEdgeOfBasementPlane[2] =n[1];
        mEdgeOfBasementPlane[3] = EdgeCentroid[1]; mEdgeOfBasementPlane[4] =n[2];  mEdgeOfBasementPlane[5] = EdgeCentroid[2];
        mNormalToPlane = n/norm_2(n);

        // mEdgeOfBasementPlane.push_back(n[0]);  mEdgeOfBasementPlane.push_back(EdgeCentroid[0]); mEdgeOfBasementPlane.push_back(n[1]);  mEdgeOfBasementPlane.push_back(EdgeCentroid[1]);  mEdgeOfBasementPlane.push_back(n[2]); mEdgeOfBasementPlane.push_back(EdgeCentroid[2]); 



        // Check if above or below the plane

        std::set<unsigned> MutantNodeIndicesToRemove;
        for (std::set<unsigned>::iterator iter = MutantNodeIndices.begin();
            iter != MutantNodeIndices.end();
            ++iter)
        {
             c_vector<double, 3> X = rCellPopulation.GetNode(*iter)->rGetLocation();
            if (n[0]*(X[0]-EdgeCentroid[0]) + n[1]*(X[1]-EdgeCentroid[1])  +n[2]*(X[2]-EdgeCentroid[2]) > 0)
            {
                // Want to remove???
                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                p_cell->SetMutationState(p_EC);
                MutantNodeIndicesToRemove.insert(*iter);
            }
        }
        TRACE("DONE");

        // Loop over nodes to remove and remove them from the MutantNodeIndices
        for (std::set<unsigned>::iterator iter = MutantNodeIndicesToRemove.begin();
            iter != MutantNodeIndicesToRemove.end();
            ++iter)
        {
            MutantNodeIndices.erase(*iter);
        }
        // Now again need to iterate around and refind the edge nodes -- is there an easier way? 
        EdgeMutantNodeIndices.clear();

          // Loop over the mutants and get the nodes along the boundary
        for (std::set<unsigned>::iterator iter = MutantNodeIndices.begin();
             iter != MutantNodeIndices.end();
             ++iter)
        {


            // Need to check if any of the neighbours are an EC, if so
            // Save the EC as a mutant on the boundary cell_iter->SetMutationState(p_NextToBasement);
            mBasementNodes.push_back(*iter);
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
            std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringLocationIndices(p_cell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
            for (std::set<unsigned>::iterator N_iter = neighbouring_node_indices.begin();
                 N_iter != neighbouring_node_indices.end();
                 ++N_iter)
            {

                CellPtr N_cell = rCellPopulation.GetCellUsingLocationIndex(*N_iter);
                if (N_cell->GetMutationState()->template IsType<HasEndothelialCell>())
                {
                    // This cell is type EC and next to the basement matrix, which means the
                    // this basement mutant is on the edge :)
                    EdgeMutantNodeIndices.insert(*iter);
                    // Node_location = rCellPopulation.GetNode(*iter)->rGetLocation();

                    Node_location = rCellPopulation.GetNode(*N_iter)->rGetLocation();
                }
            }
        }

        // Now find a region that is neigbhoring the membrane region

        // How about we say everything 

        // std::map<unsigned, std::vector<unsigned> > NeighboursBoundary;

        TRACE("Now about to find the edge region")
        //---------------------------------------
        //---------------------------------------
        // Now here we get the neighbourhood
        //---------------------------------------
        //---------------------------------------

        // std::vector<unsigned> NeighboursBoundary;
        std::set<unsigned> NeighboursBoundary;
        double MaxZ = 0;
        double Minz = 100;
        CellPtr EcCell;
        // vector<unsigned> NeighbourCollection;
        for (std::set<unsigned>::iterator iter = EdgeMutantNodeIndices.begin();
             iter != EdgeMutantNodeIndices.end();
             ++iter)
        {

            // Here we have an edge with a EC region, so we should grab all the ECs in the
            // between the nodies location and  Z away in either direction that is an EC
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
    
            std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringLocationIndices(p_cell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);

            for (std::set<unsigned>::iterator N_iter = neighbouring_node_indices.begin();
                 N_iter != neighbouring_node_indices.end();
                 ++N_iter)
            {

                CellPtr N_cell = rCellPopulation.GetCellUsingLocationIndex(*N_iter);
                if (N_cell->GetMutationState()->template IsType<HasEndothelialCell>())
                {
                    // Set this as a boundary
                    // N_cell->SetMutationState(p_NextToBasement);
                    // Save in the collection vector
                    NeighboursBoundary.insert(*N_iter);
                    // EcCell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                    // EcCell->SetMutationState(p_NextToBasement);
                    std::set<unsigned> neighbouring_node_indices2 = rCellPopulation.GetNeighbouringLocationIndices(N_cell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                    for (std::set<unsigned>::iterator N_iter_2 = neighbouring_node_indices2.begin();
                         N_iter_2 != neighbouring_node_indices2.end();
                         ++N_iter_2)
                    {
                        CellPtr N_cell_2 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_2);
                        if (N_cell_2->GetMutationState()->template IsType<HasEndothelialCell>())
                        {

                            NeighboursBoundary.insert(*N_iter_2);

                            std::set<unsigned> neighbouring_node_indices3 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_2); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);

                            for (std::set<unsigned>::iterator N_iter_3 = neighbouring_node_indices3.begin();
                                 N_iter_3 != neighbouring_node_indices3.end();
                                 ++N_iter_3)
                            {

                                CellPtr N_cell_3 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_3);

                                if (N_cell_3->GetMutationState()->template IsType<HasEndothelialCell>())
                                {
                                    //  TRACE("if endothilial state");
                                    NeighboursBoundary.insert(*N_iter_3);

                                    std::set<unsigned> neighbouring_node_indices4 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_3); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                    for (std::set<unsigned>::iterator N_iter_4 = neighbouring_node_indices4.begin();
                                         N_iter_4 != neighbouring_node_indices4.end();
                                         ++N_iter_4)
                                    {

                                        CellPtr N_cell_4 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_4);
                                        if (N_cell_4->GetMutationState()->template IsType<HasEndothelialCell>())
                                        {
                                            NeighboursBoundary.insert(*N_iter_4);

                                            // --------------------
                                            std::set<unsigned> neighbouring_node_indices5 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_4); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                            for (std::set<unsigned>::iterator N_iter_5 = neighbouring_node_indices5.begin();
                                                 N_iter_5 != neighbouring_node_indices5.end();
                                                 ++N_iter_5)
                                            {

                                                CellPtr N_cell_5 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_5);
                                                if (N_cell_5->GetMutationState()->template IsType<HasEndothelialCell>())
                                                {
                                                    NeighboursBoundary.insert(*N_iter_5);

                                                    // --------------------
                                                    std::set<unsigned> neighbouring_node_indices6 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_5); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                                    for (std::set<unsigned>::iterator N_iter_6 = neighbouring_node_indices6.begin();
                                                         N_iter_6 != neighbouring_node_indices6.end();
                                                         ++N_iter_6)
                                                    {

                                                        CellPtr N_cell_6 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_6);                 
                                                    }
                                                    // -------------------
                                                }
                                            }
                                            // -------------------
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Need the edges of each !!!
        
        c_vector<double, 3> BEdgeCentroid = Create_c_vector(0,0,0);
        std::set<unsigned> EdgeNeighboursBoundary;

        // Loop over the neighbours and mark as next to martix.
        for (std::set<unsigned>::iterator it = NeighboursBoundary.begin(); it != NeighboursBoundary.end(); ++it)
        {
            CellPtr BoundaryCell = rCellPopulation.GetCellUsingLocationIndex(*it);
            BoundaryCell->SetMutationState(p_NextToBasement);
            BoundaryCell->GetCellData()->SetItem("ShearModulus", mGrowthMaps[2](2));
            BoundaryCell->GetCellData()->SetItem("AreaDilationModulus", mGrowthMaps[2](1));
            BoundaryCell->GetCellData()->SetItem("AreaConstant", mGrowthMaps[2](0));
        }

        for (std::set<unsigned>::iterator it = NeighboursBoundary.begin(); it != NeighboursBoundary.end(); ++it)
        {
             CellPtr BoundaryCell = rCellPopulation.GetCellUsingLocationIndex(*it);
             std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringLocationIndices(BoundaryCell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);

             // loop over neighbours and check which ones are next to the EC 
 
            for (std::set<unsigned>::iterator N_iter = neighbouring_node_indices.begin();
                N_iter != neighbouring_node_indices.end();
                ++N_iter)
                {

                    CellPtr N_cell = rCellPopulation.GetCellUsingLocationIndex(*N_iter);
                    if (N_cell->GetMutationState()->template IsType<HasEndothelialCell>())
                    {
                        EdgeNeighboursBoundary.insert(*it);
                        // BoundaryCell->SetMutationState(p_Basement);
                        BEdgeCentroid += rCellPopulation.GetLocationOfCellCentre(BoundaryCell);
                        break ;
                    }

                }

        }
        // THis is theretically the centeral point of the out edge
        BEdgeCentroid /= EdgeNeighboursBoundary.size();


        double Min1a = 1000;
        Min2 = 1000000;
        double Max =0;
        double AngleMax =0;

        // Get the points in the edge of the edge that are the closest and furtherst from the center point
        for (std::set<unsigned>::iterator iter = EdgeNeighboursBoundary.begin(); iter != EdgeNeighboursBoundary.end(); ++iter)
        {
             double distance = norm_2(BEdgeCentroid  - rCellPopulation.GetNode(*iter)->rGetLocation() );
                   if (distance < Min1a)
                   {    
                        Min1a = distance;
                        NodeA = *iter;
                   } else if ( distance > Max)
                   {
                       Max = distance;
                       NodeB = *iter;
                   }

        }

        // Have the first point, now want a second that is a certain angle away from the first 
        c_vector<long double, 3> Vector1 = BEdgeCentroid  - rCellPopulation.GetNode(NodeA)->rGetLocation();


        for (std::set<unsigned>::iterator iter = EdgeNeighboursBoundary.begin(); iter != EdgeNeighboursBoundary.end(); ++iter)
        {
                c_vector<long double, 3> Vector2 = BEdgeCentroid  - rCellPopulation.GetNode(*iter)->rGetLocation();
                double Angle = std::abs(acos(inner_prod( Vector1, Vector2)/(norm_2(Vector1)*norm_2(Vector2)) ));
                // How about finding the largest angle??
                // if (Angle > AngleMax)
                // {
                //      AngleMax = Angle;
                //      NodeB = *iter;
                // }

                if ( Angle > 3.5*M_PI/4  )
                {
                    double distance = norm_2(Vector2 );
                   if (distance < Min2)
                   {    
                        Min2 = distance;
                        NodeB = *iter;
                   } 
                }

        }




        AX = BEdgeCentroid - rCellPopulation.GetNode(NodeA)->rGetLocation();
        BX = BEdgeCentroid - rCellPopulation.GetNode(NodeB)->rGetLocation();
        c_vector<double, 3> Bn = VectorProduct(AX,BX);
        // Equation of plane 
        // BoundaryPlane = Bn[0](x-EdgeCentroid[0]) + Bn[1](x-EdgeCentroid[1])  + Bn[2](x-EdgeCentroid[2])  =0

        // mEdgeOfECPlane
        mEdgeOfECPlane[0] = Bn[0];  mEdgeOfECPlane[1] =BEdgeCentroid[0];  mEdgeOfECPlane[2] =Bn[1];
        mEdgeOfECPlane[3] = BEdgeCentroid[1]; mEdgeOfECPlane[4] =Bn[2];  mEdgeOfECPlane[5] = BEdgeCentroid[2];


            // Check if above or below the plane

            MutantNodeIndicesToRemove.clear();
            std::vector<unsigned> NeighboursBoundary2;
            for (std::set<unsigned>::iterator iter = NeighboursBoundary.begin();
                iter != NeighboursBoundary.end();
                ++iter)
            {
                c_vector<double, 3> X = rCellPopulation.GetNode(*iter)->rGetLocation();
                if (Bn[0]*(X[0]-BEdgeCentroid[0]) + Bn[1]*(X[1]-BEdgeCentroid[1])  + Bn[2]*(X[2]-BEdgeCentroid[2]) > 0)
                {
                    // Want to remove???
                    CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                    p_cell->SetMutationState(p_EC);
                    MutantNodeIndicesToRemove.insert(*iter);

                }
                else
                {
                        NeighboursBoundary2.push_back(*iter);  
                }
                
            }
        // I need The distance between the planes --  Set this as the distance between the centroids.
        mDistanceBetweenPlanes = norm_2(BEdgeCentroid - EdgeCentroid);
        mCentroidEC = BEdgeCentroid;
        mCentroidBasement = EdgeCentroid;      
        ////-----------------------------------------
         TRACE("DONE");
        mNodesNextToBasement = NeighboursBoundary2;
    }
    else
    {
        TRACE("NOw in the else part")
        // mMinZ = 0;
        // mMaxZ = 0;
    }




    // Now impose the conditions 
    bool ImposeConditionsHetrogenatityInstantly = true;
    if (ImposeConditionsHetrogenatityInstantly)
    {
    //-----------------------------------------------------

        double Kbs =  (double)mGrowthMaps[1.2](2);
        double Kba =  (double)mGrowthMaps[1.2](1);
        double KbA =  (double)mGrowthMaps[1.2](0);
        // double Kbb = std::min((double)Step_Kb, (double)mGrowthMaps[1.2](3));
        // PRINT_4_VARIABLES(Kbs, Kba, KbA, Kbb)
        
        double Ks = mGrowthMaps[mStrength](2);
        double Ka = mGrowthMaps[mStrength](1);
        double KA = mGrowthMaps[mStrength](0);
        double d = mDistanceBetweenPlanes;

        // The endothelial Z mMinZ;
        //Shear K
        double c1S = Kbs;  double a1S = (Ks - c1S) /(d*d);

        //Alpha 
        double c1a = Kba;  double a1a = (Ka - c1a) / (d *d);//   

        //Area  
        double c1A = KbA; double a1A = (KA - c1A) / (d* d);//     
    
        double ShearMod; double AlphaMod; double AreaMod;


        for (std::vector<unsigned>::iterator it = mBasementNodes.begin(); it != mBasementNodes.end(); ++it)
        {
            CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*it);
            cell_iter->GetCellData()->SetItem("ShearModulus", Kbs);
            cell_iter->GetCellData()->SetItem("AreaDilationModulus", Kba);
            cell_iter->GetCellData()->SetItem("AreaConstant", KbA );
        }

        for (std::vector<unsigned>::iterator it = mNodesNextToBasement.begin(); it != mNodesNextToBasement.end(); ++it)
        {
            CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*it);

            // Change the membrane parameters of the nodes neighbouring the basement mutants -- do this linearly in time and quadratically in space
            c_vector<double, SPACE_DIM> location = rCellPopulation.GetLocationOfCellCentre(cell_iter);
                //K =a*Z^2 +c
                c_vector<double, 3> DistanceFromCentroid =  location -  mCentroidBasement;
                // DitanceToPlane
                double D = inner_prod(DistanceFromCentroid, mNormalToPlane);
                ShearMod = a1S * D*D + c1S;
                AlphaMod = a1a * D*D + c1a;
                AreaMod = a1A * D*D + c1A;
      
            cell_iter->GetCellData()->SetItem("ShearModulus", ShearMod);
            cell_iter->GetCellData()->SetItem("AreaDilationModulus", AlphaMod);
            cell_iter->GetCellData()->SetItem("AreaConstant", AreaMod);
        }
        //-------------------------------------------------------
       
    }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembraneHetroModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolveBasedOnGivenRegionOrNodes(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    // TRACE("In the SetUpSolve");
    // PRINT_2_VARIABLES(ELEMENT_DIM, SPACE_DIM);
    assert(SPACE_DIM == 3);
    assert(ELEMENT_DIM == 2);
    MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
    MAKE_PTR(EmptyBasementMatrix, p_Basement); //Mutation to mark nodes one basement
    MAKE_PTR(HasEndothelialCell, p_EC); //Mutation to mark nodes with ECs
    MAKE_PTR(LostEndothelialCell, p_NextToBasement); //Mutation to mark nodes with ECs next to the basement
    c_vector<long double, 3> Node_location;
 
    std::set<unsigned> MutantNodeIndices;
    std::set<unsigned> EdgeMutantNodeIndices;

    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node_location = rCellPopulation.GetNode(node_index)->rGetLocation();
        if (Node_location[1] > 0.2 && node_index < 2367 )   
        {
            cell_iter->SetMutationState(p_Basement);
            MutantNodeIndices.insert(node_index);
            mSamplebasementNode = node_index;
        }
        else
        {
            cell_iter->SetMutationState(p_EC);
            mSampleECNode = node_index;
        }
        cell_iter->GetCellData()->SetItem("NodeIndex",  node_index  );

        cell_iter->GetCellData()->SetItem("ShearModulus", mGrowthMaps[mStrength](2));
        cell_iter->GetCellData()->SetItem("AreaDilationModulus", mGrowthMaps[mStrength](1));
        cell_iter->GetCellData()->SetItem("AreaConstant", mGrowthMaps[mStrength](0));
        cell_iter->GetCellData()->SetItem("BendingConstant", mGrowthMaps[mStrength](3));

    }

    
    if (mHetro)
    {
    
        // Loop over the mutants and get the nodes along the boundary
        for (std::set<unsigned>::iterator iter = MutantNodeIndices.begin();
             iter != MutantNodeIndices.end();
             ++iter)
        {


            // Need to check if any of the neighbours are an EC, if so
            // Save the EC as a mutant on the boundary cell_iter->SetMutationState(p_NextToBasement);
            mBasementNodes.push_back(*iter);
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
            std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringLocationIndices(p_cell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
            for (std::set<unsigned>::iterator N_iter = neighbouring_node_indices.begin();
                 N_iter != neighbouring_node_indices.end();
                 ++N_iter)
            {

                CellPtr N_cell = rCellPopulation.GetCellUsingLocationIndex(*N_iter);
                if (N_cell->GetMutationState()->template IsType<HasEndothelialCell>())
                {
                    // This cell is type EC and next to the basement matrix, which means the
                    // this basement mutant is on the edge :)
                    EdgeMutantNodeIndices.insert(*iter);
                    break;
                }
            }
        }



        // SO now we have the membrane region and the lining edge of the membrane region -- can clean it up a bit


        // Find the centroid of the edges
        c_vector<double, 3> EdgeCentroid= Create_c_vector(0,0,0);

        for (std::set<unsigned>::iterator iter = EdgeMutantNodeIndices.begin();
             iter != EdgeMutantNodeIndices.end();
             ++iter)
        {
                EdgeCentroid += rCellPopulation.GetNode(*iter)->rGetLocation();
        }   
        EdgeCentroid /= EdgeMutantNodeIndices.size();

        // Now loop over and find the two minimum points

        double Min1 = 1000;
        double Min2 = 1000000;
        unsigned NodeA;
        unsigned NodeB;

         for (std::set<unsigned>::iterator iter = EdgeMutantNodeIndices.begin();
             iter != EdgeMutantNodeIndices.end();
             ++iter)
        {
               double distance = norm_2(EdgeCentroid - rCellPopulation.GetNode(*iter)->rGetLocation() );
                   if (distance < Min1)
                   {    
                        Min1 = distance;
                        NodeA = *iter;
                   } 
        }  
        c_vector<double, 3> VectorA = EdgeCentroid  -rCellPopulation.GetNode(NodeA)->rGetLocation(); 

          for (std::set<unsigned>::iterator iter = EdgeMutantNodeIndices.begin();
             iter != EdgeMutantNodeIndices.end();
             ++iter)
        {   
                c_vector<double, 3> VectorB = EdgeCentroid  - rCellPopulation.GetNode(*iter)->rGetLocation(); 
                double Angle = std::abs(acos(inner_prod(VectorA,VectorB )/(norm_2(VectorA)*norm_2(VectorB)) )) ;
                if (Angle > 3*M_PI/4)
                {
                    double distance = norm_2(VectorB);
                   if (distance < Min2)
                   {    
                        Min2 = distance;
                        NodeB = *iter;
                   } 

                }
               
        }  

        c_vector<double, 3> AX = EdgeCentroid - rCellPopulation.GetNode(NodeA)->rGetLocation();
        c_vector<double, 3> BX = EdgeCentroid - rCellPopulation.GetNode(NodeB)->rGetLocation();
        c_vector<double, 3> n = VectorProduct(AX,BX);
        // Equation of plane 
        // Plane = n[0](x-EdgeCentroid[0]) + n[1](x-EdgeCentroid[1])  +n[2](x-EdgeCentroid[2])  =0

        // mEdgeOfBasementPlane
        mEdgeOfBasementPlane[0] = n[0];  mEdgeOfBasementPlane[1] =EdgeCentroid[0];  mEdgeOfBasementPlane[2] =n[1];
        mEdgeOfBasementPlane[3] = EdgeCentroid[1]; mEdgeOfBasementPlane[4] =n[2];  mEdgeOfBasementPlane[5] = EdgeCentroid[2];
        mNormalToPlane = n/norm_2(n);

        // mEdgeOfBasementPlane.push_back(n[0]);  mEdgeOfBasementPlane.push_back(EdgeCentroid[0]); mEdgeOfBasementPlane.push_back(n[1]);  mEdgeOfBasementPlane.push_back(EdgeCentroid[1]);  mEdgeOfBasementPlane.push_back(n[2]); mEdgeOfBasementPlane.push_back(EdgeCentroid[2]); 



        // Check if above or below the plane

        std::set<unsigned> MutantNodeIndicesToRemove;
        for (std::set<unsigned>::iterator iter = MutantNodeIndices.begin();
            iter != MutantNodeIndices.end();
            ++iter)
        {
             c_vector<double, 3> X = rCellPopulation.GetNode(*iter)->rGetLocation();
            if (n[0]*(X[0]-EdgeCentroid[0]) + n[1]*(X[1]-EdgeCentroid[1])  +n[2]*(X[2]-EdgeCentroid[2]) > 0)
            {
                // Want to remove???
                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                p_cell->SetMutationState(p_EC);
                MutantNodeIndicesToRemove.insert(*iter);
            }
        }
        TRACE("DONE");

        // Loop over nodes to remove and remove them from the MutantNodeIndices
        for (std::set<unsigned>::iterator iter = MutantNodeIndicesToRemove.begin();
            iter != MutantNodeIndicesToRemove.end();
            ++iter)
        {
            MutantNodeIndices.erase(*iter);
        }
        // Now again need to iterate around and refind the edge nodes -- is there an easier way? 
        EdgeMutantNodeIndices.clear();

          // Loop over the mutants and get the nodes along the boundary
        for (std::set<unsigned>::iterator iter = MutantNodeIndices.begin();
             iter != MutantNodeIndices.end();
             ++iter)
        {


            // Need to check if any of the neighbours are an EC, if so
            // Save the EC as a mutant on the boundary cell_iter->SetMutationState(p_NextToBasement);
            mBasementNodes.push_back(*iter);
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
            std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringLocationIndices(p_cell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
            for (std::set<unsigned>::iterator N_iter = neighbouring_node_indices.begin();
                 N_iter != neighbouring_node_indices.end();
                 ++N_iter)
            {

                CellPtr N_cell = rCellPopulation.GetCellUsingLocationIndex(*N_iter);
                if (N_cell->GetMutationState()->template IsType<HasEndothelialCell>())
                {
                    // This cell is type EC and next to the basement matrix, which means the
                    // this basement mutant is on the edge :)
                    EdgeMutantNodeIndices.insert(*iter);
                    // Node_location = rCellPopulation.GetNode(*iter)->rGetLocation();

                    Node_location = rCellPopulation.GetNode(*N_iter)->rGetLocation();
                }
            }
        }

        // Now find a region that is neigbhoring the membrane region

        // How about we say everything 

        // std::map<unsigned, std::vector<unsigned> > NeighboursBoundary;


        //---------------------------------------
        //---------------------------------------
        //---------------------------------------
        //---------------------------------------






        // std::vector<unsigned> NeighboursBoundary;
        std::set<unsigned> NeighboursBoundary;
        double MaxZ = 0;
        double Minz = 100;
        CellPtr EcCell;
        // vector<unsigned> NeighbourCollection;
        for (std::set<unsigned>::iterator iter = EdgeMutantNodeIndices.begin();
             iter != EdgeMutantNodeIndices.end();
             ++iter)
        {
            // Just collect the boundaries for maxiumum ease of interpretation later

            // Here we have an edge with a EC region, so we should grab all the ECs in the
            // between the nodies location and  Z away in either direction that is an EC
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
    
            std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringLocationIndices(p_cell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);

            for (std::set<unsigned>::iterator N_iter = neighbouring_node_indices.begin();
                 N_iter != neighbouring_node_indices.end();
                 ++N_iter)
            {

                CellPtr N_cell = rCellPopulation.GetCellUsingLocationIndex(*N_iter);
                if (N_cell->GetMutationState()->template IsType<HasEndothelialCell>())
                {
                    // Set this as a boundary
                    // N_cell->SetMutationState(p_NextToBasement);
                    // Save in the collection vector
                    NeighboursBoundary.insert(*N_iter);
                    // EcCell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                    // EcCell->SetMutationState(p_NextToBasement);
                    std::set<unsigned> neighbouring_node_indices2 = rCellPopulation.GetNeighbouringLocationIndices(N_cell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                    for (std::set<unsigned>::iterator N_iter_2 = neighbouring_node_indices2.begin();
                         N_iter_2 != neighbouring_node_indices2.end();
                         ++N_iter_2)
                    {
                        CellPtr N_cell_2 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_2);
                        if (N_cell_2->GetMutationState()->template IsType<HasEndothelialCell>())
                        {

                            NeighboursBoundary.insert(*N_iter_2);

                            std::set<unsigned> neighbouring_node_indices3 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_2); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);

                            for (std::set<unsigned>::iterator N_iter_3 = neighbouring_node_indices3.begin();
                                 N_iter_3 != neighbouring_node_indices3.end();
                                 ++N_iter_3)
                            {

                                CellPtr N_cell_3 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_3);

                                if (N_cell_3->GetMutationState()->template IsType<HasEndothelialCell>())
                                {
                                    //  TRACE("if endothilial state");
                                    NeighboursBoundary.insert(*N_iter_3);

                                    std::set<unsigned> neighbouring_node_indices4 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_3); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                    for (std::set<unsigned>::iterator N_iter_4 = neighbouring_node_indices4.begin();
                                         N_iter_4 != neighbouring_node_indices4.end();
                                         ++N_iter_4)
                                    {

                                        CellPtr N_cell_4 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_4);
                                        if (N_cell_4->GetMutationState()->template IsType<HasEndothelialCell>())
                                        {
                                            NeighboursBoundary.insert(*N_iter_4);

                                            // --------------------
                                            std::set<unsigned> neighbouring_node_indices5 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_4); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                            for (std::set<unsigned>::iterator N_iter_5 = neighbouring_node_indices5.begin();
                                                 N_iter_5 != neighbouring_node_indices5.end();
                                                 ++N_iter_5)
                                            {

                                                CellPtr N_cell_5 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_5);
                                                if (N_cell_5->GetMutationState()->template IsType<HasEndothelialCell>())
                                                {
                                                    NeighboursBoundary.insert(*N_iter_5);

                                                    // --------------------
                                                    std::set<unsigned> neighbouring_node_indices6 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_5); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                                    for (std::set<unsigned>::iterator N_iter_6 = neighbouring_node_indices6.begin();
                                                         N_iter_6 != neighbouring_node_indices6.end();
                                                         ++N_iter_6)
                                                    {

                                                        CellPtr N_cell_6 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_6);
                                                        if (N_cell_6->GetMutationState()->template IsType<HasEndothelialCell>())
                                                        {
                                                            NeighboursBoundary.insert(*N_iter_6);

                                                            // -------------------
                                                            // -------------------
                                                            std::set<unsigned> neighbouring_node_indices7 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_6); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                                            for (std::set<unsigned>::iterator N_iter_7 = neighbouring_node_indices7.begin();
                                                                N_iter_7 != neighbouring_node_indices7.end();
                                                                ++N_iter_7)
                                                            {

                                                                CellPtr N_cell_7 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_7);
                                                                if (N_cell_7->GetMutationState()->template IsType<HasEndothelialCell>())
                                                                {
                                                                    NeighboursBoundary.insert(*N_iter_7);
                                                                }
                                                            }


                                                        }
                                                    }
                                                    // -------------------
                                                }
                                            }
                                            // -------------------
                                        }
                                    }
                                }
                            }
                        }
                        // Loop over the neighbours for several levels
                    }
                }
            }
        }

        // sort(NeighboursBoundary.begin(), NeighboursBoundary.end());
        // NeighboursBoundary.erase(unique(NeighboursBoundary.begin(), NeighboursBoundary.end()), NeighboursBoundary.end());
       


        // want to get the edge



        // Need the edges of each !!!
        
        c_vector<double, 3> BEdgeCentroid = Create_c_vector(0,0,0);
        std::set<unsigned> EdgeNeighboursBoundary;

        // Loop over the neighbours and mark as next to martix.
        for (std::set<unsigned>::iterator it = NeighboursBoundary.begin(); it != NeighboursBoundary.end(); ++it)
        {
            CellPtr BoundaryCell = rCellPopulation.GetCellUsingLocationIndex(*it);
            BoundaryCell->SetMutationState(p_NextToBasement);
            // BoundaryCell->GetCellData()->SetItem("ShearModulus", mGrowthMaps[2](2));
            // BoundaryCell->GetCellData()->SetItem("AreaDilationModulus", mGrowthMaps[2](1));
            // BoundaryCell->GetCellData()->SetItem("AreaConstant", mGrowthMaps[2](0));
        }

        for (std::set<unsigned>::iterator it = NeighboursBoundary.begin(); it != NeighboursBoundary.end(); ++it)
        {
             CellPtr BoundaryCell = rCellPopulation.GetCellUsingLocationIndex(*it);
             std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringLocationIndices(BoundaryCell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);

             // loop over neighbours and check which ones are next to the EC 
 
            for (std::set<unsigned>::iterator N_iter = neighbouring_node_indices.begin();
                N_iter != neighbouring_node_indices.end();
                ++N_iter)
                {

                    CellPtr N_cell = rCellPopulation.GetCellUsingLocationIndex(*N_iter);
                    if (N_cell->GetMutationState()->template IsType<HasEndothelialCell>())
                    {
                        EdgeNeighboursBoundary.insert(*it);
                        // BoundaryCell->SetMutationState(p_Basement);
                        BEdgeCentroid += rCellPopulation.GetLocationOfCellCentre(BoundaryCell);
                        break ;
                    }

                }

        }
        BEdgeCentroid /= EdgeNeighboursBoundary.size();


        double Min1a = 1000;
        Min2 = 1000000;
        double Max =0;
        double AngleMax =0;
        for (std::set<unsigned>::iterator iter = EdgeNeighboursBoundary.begin(); iter != EdgeNeighboursBoundary.end(); ++iter)
        {
             double distance = norm_2(BEdgeCentroid  - rCellPopulation.GetNode(*iter)->rGetLocation() );
                   if (distance < Min1a)
                   {    
                        Min1a = distance;
                        NodeA = *iter;
                   } else if ( distance > Max)
                   {
                       Max = distance;
                       NodeB = *iter;
                   }

        }
        // Have the first point, now want a second that is a certain angle away from the first 
        // c_vector<long double, 3> Vector1 = BEdgeCentroid  - rCellPopulation.GetNode(NodeA)->rGetLocation();


        // for (std::set<unsigned>::iterator iter = EdgeNeighboursBoundary.begin(); iter != EdgeNeighboursBoundary.end(); ++iter)
        // {
        //         c_vector<long double, 3> Vector2 = BEdgeCentroid  - rCellPopulation.GetNode(*iter)->rGetLocation();
        //         double Angle = std::abs(acos(inner_prod( Vector1, Vector2)/(norm_2(Vector1)*norm_2(Vector2)) ));
        //         // How about finding the largest angle??
        //         if (Angle > AngleMax)
        //         {
        //              AngleMax = Angle;
        //              NodeB = *iter;
        //         }

        //         // if ( Angle > 3*M_PI/4  )
        //         // {
        //         //     double distance = norm_2(Vector2 );
        //         //    if (distance < Min2)
        //         //    {    
        //         //         Min2 = distance;
        //         //         NodeB = *iter;
        //         //    } 
        //         // }

        // }




        AX = BEdgeCentroid - rCellPopulation.GetNode(NodeA)->rGetLocation();
        BX = BEdgeCentroid - rCellPopulation.GetNode(NodeB)->rGetLocation();
        c_vector<double, 3> Bn = VectorProduct(AX,BX);
        // Equation of plane 
        // BoundaryPlane = Bn[0](x-EdgeCentroid[0]) + Bn[1](x-EdgeCentroid[1])  + Bn[2](x-EdgeCentroid[2])  =0

        // mEdgeOfECPlane
        mEdgeOfECPlane[0] = Bn[0];  mEdgeOfECPlane[1] =BEdgeCentroid[0];  mEdgeOfECPlane[2] =Bn[1];
        mEdgeOfECPlane[3] = BEdgeCentroid[1]; mEdgeOfECPlane[4] =Bn[2];  mEdgeOfECPlane[5] = BEdgeCentroid[2];


            // Check if above or below the plane

            MutantNodeIndicesToRemove.clear();
            std::vector<unsigned> NeighboursBoundary2;
            for (std::set<unsigned>::iterator iter = NeighboursBoundary.begin();
                iter != NeighboursBoundary.end();
                ++iter)
            {
                c_vector<double, 3> X = rCellPopulation.GetNode(*iter)->rGetLocation();
                if (Bn[0]*(X[0]-BEdgeCentroid[0]) + Bn[1]*(X[1]-BEdgeCentroid[1])  + Bn[2]*(X[2]-BEdgeCentroid[2]) > 0)
                {
                    // Want to remove???
                    CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                    p_cell->SetMutationState(p_EC);
                    MutantNodeIndicesToRemove.insert(*iter);

                }
                else
                {
                        NeighboursBoundary2.push_back(*iter);  
                }
                
            }


            // SO now I have both planes, in the next step I can set the properties gradient and let it go 
            // NeighboursBoundary.clear();
            // std::vector<unsigned> NeighboursBoundary = NeighboursBoundary2;
            TRACE("DONE");

            // I need The distance between the planes --  Set this as the distance between the centroids. 

            mDistanceBetweenPlanes = norm_2(BEdgeCentroid - EdgeCentroid);

            mCentroidEC = BEdgeCentroid;
            mCentroidBasement = EdgeCentroid;
            
        ////-----------------------------------------



        // Now loop over and find the two minimum points

        // Now need to figure out the interplation 

        TRACE("Out of loop");
   
        mNodesNextToBasement = NeighboursBoundary2;
        // mBasementNodes = MutantNodeIndices;
    }
    else
    {
        TRACE("NOw in the else part")
        mMinZ = 0;
        mMaxZ = 0;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembraneHetroModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{

    // // TRACE("Should be empty");

    //   assert(SPACE_DIM == 3);

    // // See if any edges are too long and if so divide them
    // double num_cells = rCellPopulation.GetNumRealCells();

    // // MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_rCellPopulation = dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(rCellPopulation));
    // if (mOn == 1 && mAchievedTargetK == 0)
    // {
    //     // if (mCounter ==100)
    //     // {
    //             UpdateCellData(rCellPopulation);
    //             // mCounter =0;
    //     // } 
    //     // else 
    //     // {
    //     //     mCounter +=1;
    //     // }
        
    // }

}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembraneHetroModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{

    // TRACE("Should be empty");
   
   
    /*
        This needs to be properly implimented 
        -- Temp --
        Initally the membrane properties are set for 20X Growth, 
        When the radius gets to 20X inital, the membrane properties in the center will be changed to decreased

        -- To Impliment -- 
        Need to adapt this to sweep over all nodes and adapt the membrane properties of the ones that have lost or gained a Potts lattice
	*/


    assert(SPACE_DIM == 3); // Currently assumes that SPACE_DIM = 3

     if (mAchievedTargetK == 0)
    {

        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(mSamplebasementNode);

        double Step_Kbs = p_cell->GetCellData()->GetItem("ShearModulus") + 1e-12;
        double Step_Kba = p_cell->GetCellData()->GetItem("AreaDilationModulus") + 1e-12;
        double Step_KbA = p_cell->GetCellData()->GetItem("AreaConstant") + 1e-12;

        double Kbs = std::min((double)Step_Kbs, (double)mGrowthMaps[1.2](2));
        double Kba = std::min((double)Step_Kba, (double)mGrowthMaps[1.2](1));
        double KbA = std::min((double)Step_KbA, (double)mGrowthMaps[1.2](0));
        // double Kbb = std::min((double)Step_Kb, (double)mGrowthMaps[1.2](3));
        // PRINT_4_VARIABLES(Kbs, Kba, KbA, Kbb)
        double Ks = mGrowthMaps[5](2);
        double Ka = mGrowthMaps[5](1);
        double KA = mGrowthMaps[5](0);
        // double KB = mGrowthMaps[5](3);

            double d = mDistanceBetweenPlanes;

            // The endothelial Z mMinZ;

            //Shear K
            double c1S = Kbs;      
            double a1S = (Ks - c1S) /(d*d);

            //Alpha 
            double c1a = Kba;    
            double a1a = (Ka - c1a) / (d *d);//   
            

            //Area  
            double c1A = KbA;       
            double a1A = (KA - c1A) / (d* d);//     
            

            double ShearMod;
            double AlphaMod;
            double AreaMod;
            // double BendingMod;

    
            // MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);

        for (std::vector<unsigned>::iterator it = mBasementNodes.begin(); it != mBasementNodes.end(); ++it)
        {

            CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*it);
            cell_iter->GetCellData()->SetItem("ShearModulus", Kbs);
            cell_iter->GetCellData()->SetItem("AreaDilationModulus", Kba);
            cell_iter->GetCellData()->SetItem("AreaConstant", KbA );
  
        }

        for (std::vector<unsigned>::iterator it = mNodesNextToBasement.begin(); it != mNodesNextToBasement.end(); ++it)
        {

            CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*it);

            // Change the membrane parameters of the nodes neighbouring the basement mutants -- do this linearly in time and quadratically in space
            c_vector<double, SPACE_DIM> location = rCellPopulation.GetLocationOfCellCentre(cell_iter);
           
                //K =a(Z -b)^2 +c
                c_vector<double, 3> DistanceFromCentroid =  location -  mCentroidBasement;
                // DitanceToPlane
                double D = inner_prod(DistanceFromCentroid, mNormalToPlane);
                ShearMod = a1S * D*D + c1S;
                AlphaMod = a1a * D*D + c1a;
                AreaMod = a1A * D*D + c1A;
      
            cell_iter->GetCellData()->SetItem("ShearModulus", ShearMod);
            cell_iter->GetCellData()->SetItem("AreaDilationModulus", AlphaMod);
            cell_iter->GetCellData()->SetItem("AreaConstant", AreaMod);
        }
        if (Step_KbA > mGrowthMaps[1.2](0) && Step_Kba > mGrowthMaps[1.2](1) && Step_Kbs > mGrowthMaps[1.2](2))
            {
                mAchievedTargetK = 1;
                TRACE("Area, shear and dilation hit limit");
                mOn = 0;
            }
    }
    else if (mAchievedTargetK == 1)
    {
        TRACE("Target achieved, shouldnt get here ")
    }
}

// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// double  MembraneHetroModifier<ELEMENT_DIM, SPACE_DIM>::FindMin(double N1, double N2)
// {
//         if  (N1 > N2)
//         {
//             return N2;
//         } else if  (N2 > N1)
//         {
//             return N1;
//         }

// }

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembraneHetroModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class MembraneHetroModifier<1, 1>;
template class MembraneHetroModifier<1, 2>;
template class MembraneHetroModifier<2, 2>;
template class MembraneHetroModifier<1, 3>;
template class MembraneHetroModifier<2, 3>;
template class MembraneHetroModifier<3, 3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MembraneHetroModifier)

// Left side
// double Zb = mMinBasementZ;
// double m1S = (Ks - Kbs) / (mMinZ - Zb);
// double c1S = Ks - mMinZ * m1S;

// double m1a = (Ka - Kba) / (mMinZ - Zb);
// double c1a = Ka - mMinZ * m1a;

// double m1A = (KA - KbA) / (mMinZ - Zb);
// double c1A = KA - mMinZ * m1A;

// // Right side
// Zb = mMaxBasementZ;//1e-3;

// double m2S = (Ks - Kbs) / (mMaxZ - Zb);
// double c2S = Ks - m2S * mMaxZ;

// double m2a = (Ka - Kba) / (mMaxZ - Zb);
// double c2a = Ka - m2a * mMaxZ;

// double m2A = (KA - KbA) / (mMaxZ - Zb);
// double c2A = KA - m2A * mMaxZ;

// // Quadradtic

// // Left side
// // Zb = -1e-3;
// Zb =  mMinBasementZ;
// m1S = (Ks - Kbs) / (mMinZ*mMinZ - Zb*Zb);
// c1S =  Ks - m1S * mMinZ*mMinZ;

// m1a = (Ka - Kba) / (mMinZ*mMinZ - Zb*Zb);
// c1a =  Ka - m1a * mMinZ*mMinZ;

// m1A = (KA - KbA) / (mMinZ*mMinZ - Zb*Zb);
// c1A =  KA - m1A * mMinZ*mMinZ;

// // Right side
// // Zb = 1e-3;
// Zb = mMaxBasementZ;
// m2S = (Ks - Kbs) / (mMaxZ*mMaxZ - Zb*Zb);
// c2S =  Ks - m2S * mMaxZ*mMaxZ;

// m2a = (Ka - Kba) / (mMaxZ*mMaxZ - Zb*Zb);
// c2a =  Ka - m2a * mMaxZ*mMaxZ;

// m2A = (KA - KbA) / (mMaxZ*mMaxZ - Zb*Zb);
// c2A =  KA - m2A * mMaxZ*mMaxZ;

// Zb = -1e-3;
// Kbs = mGrowthMaps[1.2](2);
// Kba = mGrowthMaps[1.2](2);
// KbA = mGrowthMaps[1.2](2);

// / // Quadradtic y =a(x -b)^2 +c

//     // Left side
//     // Zb = -1e-3;
//     double Zb
//     = mMinBasementZ;
// double dz = Zb - mMinZ;

// double a1S = (Kbs - Ks) / (mMinZ * mMinZ - Zb * Zb); // (-Ks + Kbs)/( pow(mMinZ- Zb,2) -1 );
// double b1 = mMinZ + Zb; //Zb;//
// double c1S = Ks - a1S * Zb * Zb;

// double a1a = (Kba - Ka) / (mMinZ * mMinZ - Zb * Zb); //  (Ka - Kba)/( pow(mMinZ- Zb,2) -1 );
// double c1a = Ka - a1a * Zb * Zb; //Kba - a1a;

// double a1A = (KbA - KA) / (mMinZ * mMinZ - Zb * Zb); //   (KA - KbA)/( pow(mMinZ- Zb,2) -1 );
// double c1A = KA - a1A * Zb * Zb; //KbA - a1A;

// // RIght side
// // Zb = -1e-3;
// Zb = mMaxBasementZ;
// dz = Zb - mMaxZ;

// double a2S = (Kbs - Ks) / (mMaxZ * mMaxZ - Zb * Zb); // (-Ks + Kbs)/( pow(mMaxZ- Zb,2) -1 );
// double b2 = mMaxZ + Zb; //Zb;//
// // c1S = Kbs - a1S;
// double c2S = Ks - a1S * Zb * Zb;

// double a2a = (Kba - Ka) / (mMaxZ * mMaxZ - Zb * Zb); //  (Ka - Kba)/( pow(mMaxZ- Zb,2) -1 );
// double c2a = Ka - a1a * Zb * Zb; //Kba - a1a;

// double a2A = (KbA - KA) / (mMaxZ * mMaxZ - Zb * Zb); //   (KA - KbA)/( pow(mMinZ- Zb,2) -1 );
// double c2A = KA - a1A * Zb * Zb; //KbA - a1A;

// // Quadradtic

// MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
// for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
//      cell_iter != rCellPopulation.End();
//      ++cell_iter)
// {
//     if (cell_iter->GetMutationState()->template IsType<EmptyBasementMatrix>())
//     {
//         cell_iter->GetCellData()->SetItem("ShearModulus", mGrowthMaps[2](2));
//         cell_iter->GetCellData()->SetItem("AreaDilationModulus", mGrowthMaps[2](1));
//         cell_iter->GetCellData()->SetItem("AreaConstant", mGrowthMaps[2](0));
//     }
//     else if (cell_iter->GetMutationState()->template IsType<LostEndothelialCell>())
//     {
//         // SO now easy K = m1 Z + c1
//         c_vector<double, SPACE_DIM> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
//         if (location[2] < 0)
//         {
//             ShearMod = m1S * location[2] + c1S;
//             AlphaMod = m1a * location[2] + c1a;
//             AreaMod = m1A * location[2] + c1A;

//             // ShearMod = m1S *pow(location[2]  - c1S,2);
//             // AlphaMod = m1a *pow(location[2]  - c1a,2);
//             // AreaMod  = m1A *pow(location[2]  - c1A,2);
//         }
//         else
//         {
//             ShearMod = m2S * location[2] + c2S;
//             AlphaMod = m2a * location[2] + c2a;
//             AreaMod = m2A * location[2] + c2A;

//             // ShearMod = m2S *pow(location[2]  - c2S,2);
//             // AlphaMod = m2a *pow(location[2]  - c2a,2);
//             // AreaMod  = m2A *pow(location[2]  - c2A,2);
//         }

//         cell_iter->GetCellData()->SetItem("ShearModulus", ShearMod);
//         cell_iter->GetCellData()->SetItem("AreaDilationModulus", AlphaMod);
//         cell_iter->GetCellData()->SetItem("AreaConstant", AreaMod);
//     }
// }

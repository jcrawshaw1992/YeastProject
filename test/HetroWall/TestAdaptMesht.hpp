
#ifndef TestBendingForceParameterSweep_HPP_
#define TestBendingForceParameterSweep_HPP_

#include <cmath>
#include <cstdio>
#include <ctime>
#include <cxxtest/TestSuite.h>
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


#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "VtkMeshReader.hpp"
#include "VtkMeshWriter.hpp"
#include "MembraneForcesBasic.hpp"
#include "MembraneHetroModifier.hpp"

// #include "ConstantPressure.hpp"
#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"
#include "CellMutationStatesWriter.hpp"

#include "OutwardsPressure.hpp"

#include "projects/VascularRemodelling/src/mechanics/MembraneStiffnessForce.hpp"

static const double M_TIME_FOR_SIMULATION = 100; //40; //50
static const double M_SAMPLING_TIME_STEP = 100; //50
static const double M_TIME_STEP = 0.002;


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
    void TestCollapseCode() throw(Exception)
    {
        
            // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
            // in um will be too large and break chaste without carefull playing with or a tiny time step

        // Birfucation
        std::string mesh_file = "projects/VascularRemodelling/test/data/bifurcation_cut/Scalled/configChaste.vtu";

        // This data file is in mm??
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        // scale = 1e-3; // so distances are in m
    // for (typename AbstractMesh<2,3>::NodeIterator node_iter = p_mesh.GetNodeIteratorBegin();
    //         node_iter != p_mesh.GetNodeIteratorEnd();
    //         ++node_iter)
    // {
    //     // Want to find all of the node areas, save them 
    // }
    double AverageArea = 0;
    double MinArea = 1000;
    double MaxArea = 0;
    double N = p_mesh.GetNumElements();
    double NodeN = p_mesh.GetNumNodes();
    std::map<unsigned , double > AreaMap;
     for (typename MutableMesh<2,3>::ElementIterator elem_iter = p_mesh.GetElementIteratorBegin();
         elem_iter != p_mesh.GetElementIteratorEnd();
         ++elem_iter)
        {
            // Need to get the three nodes associated with this element 
                // unsigned element_index = *(elem_iter);
                unsigned element_index =elem_iter->GetIndex();
                Node<3>* pNode0 = p_mesh.GetNode(elem_iter->GetNodeGlobalIndex(0));
                Node<3>* pNode1 = p_mesh.GetNode(elem_iter->GetNodeGlobalIndex(1));
                Node<3>* pNode2 = p_mesh.GetNode(elem_iter->GetNodeGlobalIndex(2));

                c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
                c_vector<double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3
                c_vector<double, 3> normalVector = VectorProduct(vector_12, vector_13);
                double Area = norm_2(normalVector)/2;
                AreaMap[element_index] = Area;
                if (Area < MinArea)
                {
                    MinArea = Area;
                } else if (Area > MaxArea)
                {
                    MaxArea = Area;
                }
                AverageArea +=Area;
        }

        AverageArea/=N;


    PRINT_4_VARIABLES(MinArea, MaxArea, AverageArea, NodeN);


    // Loop over everything and find what you need to get rid of 
    std::set<unsigned> TaggedElements;
    std::map<unsigned, std::vector<unsigned> > NeighbourElementMap;
    for (typename MutableMesh<2,3>::ElementIterator elem_iter = p_mesh.GetElementIteratorBegin();
            elem_iter != p_mesh.GetElementIteratorEnd();
            ++elem_iter)
            {
                // Need to get the three nodes associated with this element 
                    // unsigned element_index = *(elem_iter);
                    unsigned element_index =elem_iter->GetIndex();
                    if ( AreaMap[element_index] <3*MinArea)
                    {
                        // Element needs to be saved 
                        TaggedElements.insert(element_index);
                    }
            }

    for (std::set<unsigned>::iterator Tagged_iter = TaggedElements.begin();
                            Tagged_iter != TaggedElements.end();
                            ++Tagged_iter)
                    {
                     TRACE("I want to get rid of this element");

                     Element<2, 3>* elem_iter = p_mesh.GetElement(*Tagged_iter);

                     // Save it in a new map, and save the assocciated nodes 

                     // get the nodes in this element 

                        Node<3>* pNode0 = p_mesh.GetNode(elem_iter->GetNodeGlobalIndex(0));
                        Node<3>* pNode1 = p_mesh.GetNode(elem_iter->GetNodeGlobalIndex(1));
                        Node<3>* pNode2 = p_mesh.GetNode(elem_iter->GetNodeGlobalIndex(2));

                        c_vector<double, 3> CenterPoint = (pNode2->rGetLocation() + pNode1->rGetLocation() + pNode0->rGetLocation())/3;
                        
                        Node<3>* p_new_node = new Node<3>(0u, CenterPoint);
                        p_mesh.AddNode(p_new_node);

                        // I want to loop over the containing elements for each of these nodes and save them in a new set
                        // in a map. Then from this map I can remove things??? 
                        for (int i= 0; i<3; i++)
                        {
                            std::set<unsigned> containing_elem_indices = p_mesh.GetNode(elem_iter->GetNodeGlobalIndex(i))->rGetContainingElementIndices();
                            unsigned ReplaceableNode = elem_iter->GetNodeGlobalIndex(i);
                            
                            for (std::set<unsigned>::iterator neighbour_iter = containing_elem_indices.begin();
                            neighbour_iter != containing_elem_indices.end();
                            ++neighbour_iter)
                            {
                                if ( *neighbour_iter != *Tagged_iter)
                                {
                                  NeighbourElementMap[*Tagged_iter].push_back(*neighbour_iter);
                                }
                            }
                        }

                        std::vector<unsigned> OrginalElements = NeighbourElementMap[*Tagged_iter];
                        sort(NeighbourElementMap[*Tagged_iter].begin(), NeighbourElementMap[*Tagged_iter].end());
                        NeighbourElementMap[*Tagged_iter].erase(unique(NeighbourElementMap[*Tagged_iter].begin(), NeighbourElementMap[*Tagged_iter].end()), NeighbourElementMap[*Tagged_iter].end());

                        std::vector<unsigned> diff;
                        std::set_difference(NeighbourElementMap[*Tagged_iter].begin(), NeighbourElementMap[*Tagged_iter].end(), OrginalElements.begin(), OrginalElements.end(),
                                std::inserter(diff, diff.begin()));

                        for (std::vector<unsigned>::iterator Terminal_Elements =  NeighbourElementMap[*Tagged_iter].begin();
                                    Terminal_Elements != NeighbourElementMap[*Tagged_iter].end();
                                    ++Terminal_Elements)
                            {
                              // Get the element, and the nodes attached to the element 
                                Element<2, 3>* p_element = p_mesh.GetElement(*Terminal_Elements);
                                unsigned elem_index = p_element->GetIndex();
                                std::set<unsigned> NodeIndexSet;
                                

                                NodeIndexSet.insert(p_element->GetNodeGlobalIndex(0));
                                NodeIndexSet.insert(p_element->GetNodeGlobalIndex(1));
                                NodeIndexSet.insert(p_element->GetNodeGlobalIndex(2));
                                PRINT_VARIABLE(NodeIndexSet.size());
                                // Remove which ever of the original nodes is the one of problems
                                NodeIndexSet.erase(elem_iter->GetNodeGlobalIndex(0));
                                NodeIndexSet.erase(elem_iter->GetNodeGlobalIndex(1));
                                NodeIndexSet.erase(elem_iter->GetNodeGlobalIndex(2));
                                if (NodeIndexSet.size() ==2)
                                {   
                                    std::vector<Node<3>* > NodesForNewElement;                     
                                    for (std::set<unsigned>::iterator iter = NodeIndexSet.begin();
                                        iter != NodeIndexSet.end();
                                        ++iter)
                                    {
                                        NodesForNewElement.push_back( p_mesh.GetNode(*iter ) );
                                    }
                                    NodesForNewElement.push_back(p_new_node);
                                    Element<2,3>* NewElem = new Element<2,3>(0, NodesForNewElement,1);
                                    // Add the new element and remove the old
                                    p_mesh.AddElement(NewElem);
                                    p_mesh.DeleteElement(elem_index);
                                } else if (NodeIndexSet.size() ==1)
                                {
                                    Element<2, 3>* p_element = p_mesh.GetElement(*Terminal_Elements);
                                    unsigned elem_index = p_element->GetIndex();
                                    p_mesh.DeleteElement(elem_index);
                               
                                }
                    
                            }
                        // for (std::vector<unsigned>::iterator Terminal_Elements =  diff.begin();
                        //             Terminal_Elements != diff.end();
                        //             ++Terminal_Elements)
                        //     {
                        //         Element<2, 3>* p_element = p_mesh.GetElement(*Terminal_Elements);
                        //             unsigned elem_index = p_element->GetIndex();
                        //             p_mesh.DeleteElement(elem_index);
                        //     }
                        

                        //  now delete the old nodes 
                        // p_mesh.DeleteNode(elem_iter->GetNodeGlobalIndex(0));
                        // p_mesh.DeleteNode(elem_iter->GetNodeGlobalIndex(1));
                        // p_mesh.DeleteNode(elem_iter->GetNodeGlobalIndex(2));



                      }

                      std::string output_dir= "TryingToMakeANewMesh/";
                      VtkMeshWriter<2,3> mesh_writer(output_dir, "JessNewMeshconfig", false);
                      mesh_writer.WriteFilesUsingMesh(p_mesh);

                    

                        // -- If there is repetition, then this element has appeared twice 



                        // // Get the neighbouring elements 
                        // std::set<unsigned> containing_elem_indices = pNode0->rGetContainingElementIndices();
                        // unsigned ReplaceableNode = elem_iter->GetNodeGlobalIndex(0);
                        // // Loop over set and get the nodes in this set carefully -- skip element index

                        // for (std::set<unsigned>::iterator neighbour_iter = containing_elem_indices.begin();
                        //     neighbour_iter != containing_elem_indices.end();
                        //     ++neighbour_iter)
                        // {
                    
                        //     if ( *neighbour_iter != element_index)
                        //     {
                        //       // Get the element, and the nodes attached to the element 
                        //         Element<2, 3>* p_element = p_mesh.GetElement(*neighbour_iter);
                        //         unsigned elem_index = p_element->GetIndex();
                        //         std::set<unsigned> NodeIndexSet;

                        //         NodeIndexSet.insert(p_element->GetNodeGlobalIndex(0));
                        //         NodeIndexSet.insert(p_element->GetNodeGlobalIndex(1));
                        //         NodeIndexSet.insert(p_element->GetNodeGlobalIndex(2));
                        //         NodeIndexSet.erase(ReplaceableNode );
                        //         std::vector<Node<3>* > NodesForNewElement;                     
                        //         for (std::set<unsigned>::iterator iter = NodeIndexSet.begin();
                        //             iter != NodeIndexSet.end();
                        //             ++iter)
                        //         {
                        //             NodesForNewElement.push_back( p_mesh.GetNode(*iter ) );
                        //         }
                        //         NodesForNewElement.push_back(p_new_node);

                        //         Element<2,3>* NewElem = new Element<2,3>(0, NodesForNewElement,1);
                        //         // Add the new element and remove the old
                        //         p_mesh.AddElement(NewElem);
                        //         p_mesh.DeleteElement(elem_index);
                        //         TRACE("ADDed an elemenet???");
                        //     }
                        // //     TRACE("Out of IF");


                        // }
                        TRACE("Out of neighbours loop ");
                   
                    TRACE("Out of IF area test ");
                    

            TRACE("Out of element loop ");


    //  NodeN = p_mesh.GetNumNodes();
    // PRINT_VARIABLE(NodeN );
    }
};

// Add the new element
                        //  mrMesh.AddElement(p_element);
                        // mpMutableMesh->AddNode(pNewNode)

                        // I need to remove the three nodes above from the nodes list and add the Center point. 

                        // First add center point to the nodes list 


                          //     THis didnt work, dont know why 


                            // Am up to the bit where I have created  a new node and the new elements for only one of the nodes 
                            // in the orriginal elemet that I am getting rid of 

                            // Now I neew to remove the elements, remove the nodes, but first fixing up the edges
                                // Get the nodes 
                                // for (unsigned i =1; i<3 ; i++)
                                // {
                                //     if (ReplaceableNode != p_element->GetNodeGlobalIndex(i) )
                                //     {
                                //         if (HaveANode ==0)
                                //         {
                                //          pNeighbourNode0 = p_mesh.GetNode(p_element->GetNodeGlobalIndex(0));
                                //          HaveANode =1;
                                //         }
                                //         else
                                //         {
                                //         pNeighbourNode1 = p_mesh.GetNode(p_element->GetNodeGlobalIndex(0));
                                //         }
                                //     }
                                // }
                                
                                // Node<3>* pNode0 = p_mesh.GetNode(elem_iter->GetNodeGlobalIndex(0));
                                // Node<3>* pNode1 = p_mesh.GetNode(elem_iter->GetNodeGlobalIndex(1));
                                // Node<3>* pNode2 = p_mesh.GetNode(elem_iter->GetNodeGlobalIndex(2));



#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/


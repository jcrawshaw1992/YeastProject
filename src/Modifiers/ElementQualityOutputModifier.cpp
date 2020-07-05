
/*

This code controles the membrame properties. 

At the start of a simulation, the intial conditions delievered at 20(?)X less than original starting configuration of the vessel 
The system then develops to an equlibirum (that will be the realisitc size of a network). After this point the memnbrane properties will 
be variable.

For the inital impletmentation of this code, the membrane properties in the center will be changed to decrease the vessel size to 2X the intial 
condition, with some smoothing on the vessel properties just outside the region. 

After this, the membrane constants will be dependent on the number of cells in the region and the denoted mutation state. 

*/

#include "ElementQualityOutputModifier.hpp"
#include <algorithm>
#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

#include "OutsideFLuidSimulationMutation.hpp"

#include <math.h>
#include "Debug.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementQualityOutputModifier<ELEMENT_DIM, SPACE_DIM>::ElementQualityOutputModifier()
        : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementQualityOutputModifier<ELEMENT_DIM, SPACE_DIM>::~ElementQualityOutputModifier()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementQualityOutputModifier<ELEMENT_DIM, SPACE_DIM>::SetWriteoutInterval(double WriteoutInterval)
{
    // TRACE("SetWriteoutInterval")
    mWriteoutInterval = WriteoutInterval;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementQualityOutputModifier<ELEMENT_DIM, SPACE_DIM>::SetElementMetricsFileName(std::string FileName)
{
    // TRACE("SetElementMetricsFileName")
    // mFileName = FileName + "/ElementMetricsFile.txt";
    mAreaFileName = FileName + "/ElementAreas.txt";
    mInternalAnglesFileName = FileName + "/ElementInternalAngles.txt";
    mAspectRatioFileName = FileName + "/ElementAspectRatio.txt";

    mRadiusRatioFileName = FileName + "/ElementRadiusRatio.txt";
    mEdgeLengthsFileName = FileName + "/ElementEdgeLengths.txt";
    mMeanChangeInLocalCurvatureFileName = FileName + "/MeanChangeInLocalCurvature.txt";
    mElementRadiFileName = FileName + "/ElementElementRadi.txt";
    mReadMe = FileName + "/readme.txt";

    mEdgeAnglefile = FileName + "/EdgeAnglefile.txt";
    mLocalNodesForEdege = FileName + "/LocalNodesForEdege.txt";
    mContainingElementsForEdges = FileName + "/ContainingElementsForEdge.txt";
    mRegionofContainingEdgesForEdges = FileName + "/RegionofContainingEdgesForEdges.txt";
    mTimeFiles = FileName + "/OutPutTimePoints.txt";
    mTimeCounter =0;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementQualityOutputModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    // Write a discriptive title at the top of the output files 

    // Record information for each of the output files in a readme
    mOutputFile.open (mReadMe);
    mOutputFile << mAreaFileName << "\n Area of each element at each time step. Elements are separated by a space, and the next time step is denoted by a line break. The number of elements will change in the event of remeshing.\n\n";
    mOutputFile << mInternalAnglesFileName << "\n The three internal angles in each triangular element at each time step. Elements are separated by a commer, the angles for each elements are separated by a space,  and the next time step is denoted by a line break. The number of elements will change in the event of remeshing.\n\n";
    mOutputFile << mAspectRatioFileName << "\n Aspect ratio of each element at each time step. Aspect ratio is the ratio betwee the max and min axis of the element. AR = 4/sqrt(3) * ElementArea/(MaxEdgeLength * MaxEdgeLength); Elements are separated by a space, and the next time step is denoted by a line break. The number of elements will change in the event of remeshing.\n\n";
    mOutputFile << mRadiusRatioFileName << " \n Radius ratio of each element at each time step.  RR = 16* ElementArea*ElementArea /(l1*l2*l3*( l1+l2+l3 ) ); Elements are separated by a space, and the next time step is denoted by a line break. The number of elements will change in the event of remeshing.\n\n";
    mOutputFile << mEdgeLengthsFileName << "\n  Edge Lengths of each element at each time step. Elements are separated by a space, edges are separated by a commer, and the next time step is denoted by a line break. The number of elements will change in the event of remeshing.\n\n";
    mOutputFile << mMeanChangeInLocalCurvatureFileName << "\n The mean angle between the normal of each element from that of the surrounding elements (any element sharing a node with the element in question). Elements are separated by a space, and the next time step is denoted by a line break. The number of elements will change in the event of remeshing.\n\n";
    mOutputFile<< mElementRadiFileName << "\n The inner and outer radi of each element. Elements are separated by a commer, and the next time step is denoted by a line break. The number of elements will change in the event of remeshing.\n\n";
    mOutputFile << " Output interval = "<<mWriteoutInterval <<"\n\n" ;
    mOutputFile.close();
    mOutPutCounter =0;

    WriteOutDataToFiles(rCellPopulation);
    
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementQualityOutputModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    mTimeCounter +=1;
    if  (mWriteoutIntervalCounter == mWriteoutInterval-1)  
    {
        WriteOutDataToFiles(rCellPopulation);
        mWriteoutIntervalCounter =0;
    }
    else
    {
        mWriteoutIntervalCounter +=1;
    }
        
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementQualityOutputModifier<ELEMENT_DIM, SPACE_DIM>::WriteOutDataToFiles(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    
    assert(SPACE_DIM ==3);
    HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation = static_cast<HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
    // Here we need to write and append 
  
    // Open the file 
    ofstream oAreafile(mAreaFileName, ios::out | ios::app);
    ofstream oInternalAnglesfile(mInternalAnglesFileName, ios::out | ios::app);
    ofstream oAspectRatiofile(mAspectRatioFileName, ios::out | ios::app);

    ofstream oRadiusRatiofile(mRadiusRatioFileName, ios::out | ios::app);
    ofstream oEdgeLengthsfile(mEdgeLengthsFileName, ios::out | ios::app);
    ofstream oMeanChangeInLocalCurvaturefile(mMeanChangeInLocalCurvatureFileName, ios::out | ios::app);
    ofstream oElementRadifile(mElementRadiFileName, ios::out | ios::app);

    ofstream oEdgeAnglefile(mEdgeAnglefile, ios::out | ios::app);
    ofstream oLocalNodesForEdege(mLocalNodesForEdege, ios::out | ios::app);
    ofstream oContainingElementsForEdges(mContainingElementsForEdges, ios::out | ios::app);
    ofstream oRegionofContainingEdgesForEdges(mRegionofContainingEdgesForEdges, ios::out | ios::app);

    ofstream oTimeFile(mTimeFiles, ios::out | ios::app);
    oContainingElementsForEdges << mTimeCounter << " ";
    oTimeFile << mTimeCounter << "\n ";
    

    // Loop over the elements :) 
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = pCellPopulation->rGetMesh().GetElementIteratorBegin();
        elem_iter != pCellPopulation->rGetMesh().GetElementIteratorEnd();
        ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        Node<SPACE_DIM>* pNode0 = pCellPopulation->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0));
        Node<SPACE_DIM>* pNode1 = pCellPopulation->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(1));
        Node<SPACE_DIM>* pNode2 = pCellPopulation->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));

        c_vector<double, SPACE_DIM> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
        c_vector<double, SPACE_DIM> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3
        c_vector<double, SPACE_DIM> vector_23 = pNode2->rGetLocation() - pNode1->rGetLocation(); // Vector 1 to 3

        double l1 = norm_2(vector_12); // Edge Length 1
        double l2 = norm_2(vector_13); // l2
        double l3 = norm_2(vector_23); // EdgeLength3

        double ElementArea = 0.5*norm_2(VectorProduct(vector_12, vector_13));
        oAreafile << elem_index<< " " << ElementArea << " "; 

        double Angle1 = acos(inner_prod(vector_12, vector_13) / (l1 * l2));
        double Angle2 = acos(inner_prod(vector_23, vector_13) / (l3 * l2));
        double Angle3 = acos(inner_prod(vector_12, vector_23) / (l1 * l3));

        double MaxEdgeLength = std::max(std::max(l1,  l2) , std::max(l2,l3)); 

        double AspectRatio = 4/sqrt(3) * ElementArea/(MaxEdgeLength * MaxEdgeLength);
        double RadiusRatio = 16* ElementArea*ElementArea /(l1*l2*l3*( l1+l2+l3 ) );
        double RadiusInner = 2* ElementArea/(l1+l2+l3);
        double RadiusOuter = (l1*l2*l3)/(4*ElementArea);
 
        oAreafile << ElementArea << ", "; 
        oInternalAnglesfile << Angle1 << " " <<Angle2 << " " <<Angle3 << ", "; 
        oAspectRatiofile << AspectRatio << ", ";
    
        oRadiusRatiofile << RadiusRatio << ", ";
        oEdgeLengthsfile << l1 << " " <<l2 << " " <<l3 << ", "; 
        // oMeanChangeInLocalCurvaturefile << elem_index << " "<< AverageAngleBetweenSurroundingNormals << " ";
        oElementRadifile << RadiusInner << " " <<RadiusOuter << ", "; 
          

    }

    // Need to loop over the edges 
    for (typename MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SpringIterator spring_iterator = pCellPopulation->SpringsBegin();
        spring_iterator != pCellPopulation->SpringsEnd();
        ++spring_iterator)
        {

        Node<SPACE_DIM>* pNode1 = spring_iterator.GetNodeA();
        Node<SPACE_DIM>* pNode3 = spring_iterator.GetNodeB();

        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edge = std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>(pNode1, pNode3);
        std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > nonUnitNormals;
        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> otherNodes;

        bool boundary_edge_found = pCellPopulation->CalculateElementNormals( edge, nonUnitNormals, otherNodes);

        Node<SPACE_DIM>* pNode2 = otherNodes.first;
        Node<SPACE_DIM>* pNode4 = otherNodes.second;

        if (boundary_edge_found)
        {
            continue;
        }

        c_vector<double, SPACE_DIM> normal_1 = (nonUnitNormals.first)/norm_2(nonUnitNormals.first);
        c_vector<double, SPACE_DIM> normal_2 = (nonUnitNormals.second)/norm_2(nonUnitNormals.second);

        double LocalAngle = acos(inner_prod(normal_1, normal_2));

        // Here are all the nodes connected to the two elements attached to this edge
        unsigned node_index1 = pNode1->GetIndex();
        unsigned node_index2 = pNode2->GetIndex();
        unsigned node_index3 = pNode3->GetIndex();
        unsigned node_index4 = pNode4->GetIndex();
        /*
        *  Find common Elements      
        */
        // Find the indices of the elements owned by each node
        std::set<unsigned> elements_containing_node1 = pNode1->rGetContainingElementIndices();
        std::set<unsigned> elements_containing_node3 = pNode3->rGetContainingElementIndices();

        //  Here are the two elements attached to this edge
        std::set<unsigned> shared_elements;
        std::set_intersection(elements_containing_node1.begin(),
                            elements_containing_node1.end(),
                            elements_containing_node3.begin(),
                            elements_containing_node3.end(),
                            std::inserter(shared_elements, shared_elements.begin()));

        //  Here are the all the elements contained to the two nodes making to this edges
        std::set<unsigned> AllNeighbouring_elements;
        std::set_union(elements_containing_node1.begin(), elements_containing_node1.end(),
        elements_containing_node3.begin(), elements_containing_node3.end(),
        std::inserter(AllNeighbouring_elements, AllNeighbouring_elements.begin()));
    
        oEdgeAnglefile << LocalAngle << ", "; 
        oLocalNodesForEdege << node_index1 << " " <<node_index2 << " " <<node_index3 << " " <<node_index4 << " ";  

        for (typename std::set<unsigned>::iterator i = shared_elements.begin();  i != shared_elements.end();  ++i)
            {

                oContainingElementsForEdges << *i << " ";
            }

        for (typename std::set<unsigned>::iterator i = AllNeighbouring_elements.begin();  i != AllNeighbouring_elements.end();  ++i)
            {
                    oRegionofContainingEdgesForEdges << *i << " ";
        
            }
            oRegionofContainingEdgesForEdges <<", ";
    }
    oAreafile << "\n";
    oInternalAnglesfile << "\n";
    oAspectRatiofile << "\n";

    oRadiusRatiofile<< "\n";
    oEdgeLengthsfile<< "\n";
    oMeanChangeInLocalCurvaturefile<< "\n";
    oElementRadifile<< "\n";

    oEdgeAnglefile<< "\n";
    oLocalNodesForEdege<< "\n";
    oContainingElementsForEdges<< "\n";
    oRegionofContainingEdgesForEdges<< "\n";

    // TRACE("Finished Write data to file")
    mOutPutCounter+=1;

}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementQualityOutputModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ElementQualityOutputModifier<1, 1>;
template class ElementQualityOutputModifier<1, 2>;
template class ElementQualityOutputModifier<2, 2>;
template class ElementQualityOutputModifier<1, 3>;
template class ElementQualityOutputModifier<2, 3>;
template class ElementQualityOutputModifier<3, 3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ElementQualityOutputModifier)



// Old ways to get boundaries 

    // Want to see what happens if i identify the boundaries by the size of the local elements 

//     std::map<unsigned, double> AreaOfCells;
//     double AverageCellArea = 0;
//     double NumberOfCells = 0;
//     MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
//     for ( typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
//             cell_iter != rCellPopulation.End();
//             ++cell_iter)
//     {

//         unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
//         Node<SPACE_DIM>* p_node = rCellPopulation.GetNode(node_index);

//         std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();

//         double Area =0;
//         assert(containing_elements.size() > 0);
//         for (std::set<unsigned>::iterator iter = containing_elements.begin();
//                 iter != containing_elements.end();
//                 ++iter)
//         {
//             Node<SPACE_DIM>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
//             Node<SPACE_DIM>* pNode1= p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
//             Node<SPACE_DIM>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

//             c_vector<double, SPACE_DIM> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
//             c_vector<double, SPACE_DIM> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

//             c_vector<double, SPACE_DIM> normalVector = VectorProduct(vector_12, vector_13);
//             Area+= 0.5*norm_2(normalVector)/3;
//         }
//         NumberOfCells +=1;
//         AreaOfCells[node_index] = Area;
//         AverageCellArea +=Area;
//     }

//     AverageCellArea/= NumberOfCells;

//    for ( typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
//             cell_iter != rCellPopulation.End();
//             ++cell_iter)
//         {
            
            
//             unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
//             if (AreaOfCells[node_index] < 1e-1*AverageCellArea)
//             {
//                 TRACE("Trip")
//                 cell_iter->GetCellData()->SetItem("SecondBoundary", 1);

//             }
//             else{
//                 cell_iter->GetCellData()->SetItem("SecondBoundary", 0);

//             }

//         }
    

       // The ratio of the circumradius to twice the inradius is the correct definition. There is a slightly simpler formula: 
        // AR = abc/((b+c-a)(c+a-b)(a+b-c))

        // Get the variation of the normals in the local vorrinoi region

        // std::set<unsigned> Containing_elements1 = pCellPopulation->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0))->rGetContainingElementIndices();
        // std::set<unsigned> Containing_elements2 = pCellPopulation->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0))->rGetContainingElementIndices();
        // std::set<unsigned> Containing_elements3 = pCellPopulation->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0))->rGetContainingElementIndices();

        // std::set_union(Containing_elements1.begin(), Containing_elements1.end(),
        //         Containing_elements2.begin(), Containing_elements2.end(),
        //         std::inserter(Containing_elements1, Containing_elements1.begin()));

        // std::set<unsigned> Neighbouring_elements;
        // std::set_union(Containing_elements1.begin(), Containing_elements1.end(),
        //         Containing_elements3.begin(), Containing_elements3.end(),
        //         std::inserter(Neighbouring_elements, Neighbouring_elements.begin()));


        // Want to get the average angle between the normal of the current element, and the surrounding elements 

        // c_vector<double, SPACE_DIM> CurrentNormal = - pCellPopulation->rGetMesh().GetElement(elem_index)->CalculateNormal();
        // CurrentNormal /=norm_2(CurrentNormal);
        // double AverageAngleBetweenSurroundingNormals =0;

        // for (typename std::set<unsigned>::iterator iter = Neighbouring_elements.begin();
        //     iter != Neighbouring_elements.end();
        //     ++iter)
        //     {
        //         if (*iter!=elem_index)
        //         {
        //             // need to get the normals for this element .... 
        //             // Negative as normals point inwards for these surface meshes
        //             c_vector<double, SPACE_DIM> NeighbourNormal = - pCellPopulation->rGetMesh().GetElement(*iter)->CalculateNormal();
        //             NeighbourNormal /=norm_2(NeighbourNormal);
        //             double AngleBetween = acos(inner_prod(CurrentNormal,NeighbourNormal));
        //             AverageAngleBetweenSurroundingNormals +=AngleBetween;
        //         }
                
        //     }
        //     AverageAngleBetweenSurroundingNormals/=(Neighbouring_elements.size()-1);

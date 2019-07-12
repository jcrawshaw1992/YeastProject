#include "FlowBasedCDC42StimulusProtocol.hpp"
#include "UblasCustomFunctions.hpp"
#include <numeric>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FlowBasedCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::FlowBasedCDC42StimulusProtocol(const TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, double stimulusDuration, const c_vector<double,ELEMENT_DIM>& flowDirection)
    : AbstractCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>(rMesh, stimulusDuration),
      mFlowDirection(flowDirection)
{
    // Flow direction must be a normal vector
    assert(fabs(norm_2(mFlowDirection) - 1.0) < 1e-9);

    PopulateEdgeNormalForNodesOnBoundary(rMesh);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FlowBasedCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::~FlowBasedCDC42StimulusProtocol()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double FlowBasedCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::GetCDC42ActivationRate(const Node<SPACE_DIM>* pNode, double t) const
{
    if (t < this->mStimulusDuration)
    {
        typename EdgeMapType::const_iterator edge_normal = mEdgeNormalForNodesOnBoundary.find(pNode);

        if (edge_normal != mEdgeNormalForNodesOnBoundary.end())
        {
            double dot_product = inner_prod(mFlowDirection, mEdgeNormalForNodesOnBoundary.at(pNode));
            assert(fabs(dot_product) <= 1.0);

            // The stimulus is zero in nodes with a normal defining a component in the direction of the flow, [0, 1] for those against the flow
            return -std::min(0.0, dot_product);
        }
    }

    return 0.0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FlowBasedCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::PopulateEdgeNormalForNodesOnBoundary(const TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{
    typedef std::map<unsigned, std::vector<c_vector<double,ELEMENT_DIM> > > NormalsList;
    NormalsList all_normals_found_for_node_on_boundary;

    for (TetrahedralMesh<2,2>::BoundaryElementIterator it = rMesh.GetBoundaryElementIteratorBegin();
         it != rMesh.GetBoundaryElementIteratorEnd();
         ++it)
    {
        assert((*it)->GetNumNodes() == 2);
        for (unsigned i=0; i<(*it)->GetNumNodes(); i++)
        {
            all_normals_found_for_node_on_boundary[(*it)->GetNode(i)->GetIndex()].push_back((*it)->CalculateNormal());
        }
    }

    for (typename NormalsList::const_iterator iter = all_normals_found_for_node_on_boundary.begin();
         iter != all_normals_found_for_node_on_boundary.end(); ++iter)
    {
        c_vector<double,ELEMENT_DIM> mean_vector = std::accumulate(iter->second.begin(), iter->second.end(), Create_c_vector(0,0));
        if(norm_2(mean_vector) == 0.0)
        {
            // This happens in some degenerate (yet still correct) cases, specially on regular grids.
            std::cout << "WARNING: zero normal computed for node " << iter->first << std::endl;
        }
        else
        {
            mean_vector /= norm_2(mean_vector);
        }
        mEdgeNormalForNodesOnBoundary[rMesh.GetNode(iter->first)] = mean_vector;
    }
}

// Explicit instantiation
template class FlowBasedCDC42StimulusProtocol<2, 2>;

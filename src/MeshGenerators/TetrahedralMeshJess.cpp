/*

THis version of the tet mesh doent assign boundaries 

*/

#include "TetrahedralMeshJess.hpp"
#include "TetrahedralMesh.hpp"
#include <cassert>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>

#include "BoundaryElement.hpp"
#include "Element.hpp"
#include "Exception.hpp"
#include "Node.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "RandomNumberGenerator.hpp"
#include "Warnings.hpp"

// Jonathan Shewchuk's triangle and Hang Si's tetgen
#define REAL double
#define VOID void
#include "tetgen.h"
#include "triangle.h"
#undef REAL
#undef VOID

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TetrahedralMeshJess<ELEMENT_DIM, SPACE_DIM>::TetrahedralMeshJess()
{
    Clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMeshJess<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader)
{
    assert(rMeshReader.HasNodePermutation() == false);
    this->mMeshFileBaseName = rMeshReader.GetMeshFileBaseName();

    // Record number of corner nodes
    unsigned num_nodes = rMeshReader.GetNumNodes();

    /*
     * Reserve memory for nodes, so we don't have problems with
     * pointers stored in elements becoming invalid.
     */
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    //typename std::map<std::pair<unsigned,unsigned>,unsigned>::const_iterator iterator;
    //std::map<std::pair<unsigned,unsigned>,unsigned> internal_nodes_map;

    // Add nodes
    std::vector<double> coords;
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        coords = rMeshReader.GetNextNode();
        Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index, coords, false);

        for (unsigned i = 0; i < rMeshReader.GetNodeAttributes().size(); i++)
        {
            double attribute = rMeshReader.GetNodeAttributes()[i];
            p_node->AddNodeAttribute(attribute);
        }
        this->mNodes.push_back(p_node);
    }

    //unsigned new_node_index = mNumCornerNodes;

    rMeshReader.Reset();
    // Add elements
    //new_node_index = mNumCornerNodes;
    this->mElements.reserve(rMeshReader.GetNumElements());

    for (unsigned element_index = 0; element_index < (unsigned)rMeshReader.GetNumElements(); element_index++)
    {
        ElementData element_data = rMeshReader.GetNextElementData();
        std::vector<Node<SPACE_DIM>*> nodes;

        /*
         * NOTE: currently just reading element vertices from mesh reader - even if it
         * does contain information about internal nodes (ie for quadratics) this is
         * ignored here and used elsewhere: ie don't do this:
         *   unsigned nodes_size = node_indices.size();
         */
        for (unsigned j = 0; j < ELEMENT_DIM + 1; j++) // num vertices=ELEMENT_DIM+1, may not be equal to nodes_size.
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        Element<ELEMENT_DIM, SPACE_DIM>* p_element = new Element<ELEMENT_DIM, SPACE_DIM>(element_index, nodes);

        this->mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            double attribute_value = element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }

    // Add boundary elements and nodes
    for (unsigned face_index = 0; face_index < (unsigned)rMeshReader.GetNumFaces(); face_index++)
    {
        ElementData face_data = rMeshReader.GetNextFaceData();
        std::vector<unsigned> node_indices = face_data.NodeIndices;

        /*
         * NOTE: unlike the above where we just read element *vertices* from mesh reader, here we are
         * going to read a quadratic mesh with internal elements.
         * (There are only a few meshes with internals in the face file that we might as well use them.)
         *
         */
        std::vector<Node<SPACE_DIM>*> nodes;
        for (unsigned node_index = 0; node_index < node_indices.size(); node_index++)
        {
            assert(node_indices[node_index] < this->mNodes.size());
            // Add Node pointer to list for creating an element
            nodes.push_back(this->mNodes[node_indices[node_index]]);
        }

        // This is a boundary face, so ensure all its nodes are marked as boundary nodes
        for (unsigned j = 0; j < nodes.size(); j++)
        {
            if (!nodes[j]->IsBoundaryNode())
            {
                nodes[j]->SetAsBoundaryNode();
                this->mBoundaryNodes.push_back(nodes[j]);
            }

            // Register the index that this bounday element will have with the node
            nodes[j]->AddBoundaryElement(face_index);
        }

        // The added elements will be deleted in our destructor
        BoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>* p_boundary_element = new BoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>(face_index, nodes);
        this->mBoundaryElements.push_back(p_boundary_element);

        if (rMeshReader.GetNumFaceAttributes() > 0)
        {
            assert(rMeshReader.GetNumFaceAttributes() == 1);
            double attribute_value = face_data.AttributeValue;
            p_boundary_element->SetAttribute(attribute_value);
        }
    }

    RefreshJacobianCachedData();

    rMeshReader.Reset();
}


// Explicit instantiation
template class TetrahedralMeshJess<1, 1>;
template class TetrahedralMeshJess<1, 2>;
template class TetrahedralMeshJess<1, 3>;
template class TetrahedralMeshJess<2, 2>;
template class TetrahedralMeshJess<2, 3>;
template class TetrahedralMeshJess<3, 3>;

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */


#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(TetrahedralMeshJess)



#include "ElementQualityWriter.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementQualityWriter<ELEMENT_DIM, SPACE_DIM>::ElementQualityWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("ElementQuality.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementQualityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    assert(SPACE_DIM==2 || SPACE_DIM==3); // LCOV_EXCL_LINE
    VertexMesh<ELEMENT_DIM,SPACE_DIM>* voronoi_tesselation = pCellPopulation->GetVoronoiTessellation();

    // Loop over elements of voronoi_tesselation
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem_iter = voronoi_tesselation->GetElementIteratorBegin();
         elem_iter != voronoi_tesselation->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in voronoi_tesselation
        unsigned elem_index = elem_iter->GetIndex();

        // Get the index of the corresponding node in mrMesh
        unsigned node_index = voronoi_tesselation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

        // Write node index and location to file
        *this->mpOutStream << node_index << " ";
        c_vector<double, SPACE_DIM> node_location = pCellPopulation->GetNode(node_index)->rGetLocation();
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *this->mpOutStream << node_location[i] << " ";
        }

        double cell_volume = voronoi_tesselation->GetVolumeOfElement(elem_index);
        double cell_surface_area = voronoi_tesselation->GetSurfaceAreaOfElement(elem_index);
        *this->mpOutStream << cell_volume << " " << cell_surface_area << " ";
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementQualityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("ElementQualityWriter cannot be used with a CaBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementQualityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("ElementQualityWriter cannot be used with a NodeBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementQualityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("ElementQualityWriter cannot be used with a PottsBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementQualityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("ElementQualityWriter cannot be used with a VertexBasedCellPopulation");
}

// Explicit instantiation
template class ElementQualityWriter<1,1>;
template class ElementQualityWriter<1,2>;
template class ElementQualityWriter<2,2>;
template class ElementQualityWriter<1,3>;
template class ElementQualityWriter<2,3>;
template class ElementQualityWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ElementQualityWriter)

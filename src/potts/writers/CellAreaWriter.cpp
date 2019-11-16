
#include "CellAreaWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellAreaWriter<ELEMENT_DIM, SPACE_DIM>::CellAreaWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("CellArea.dat")
{
    this->mVtkCellDataName = "Area";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellAreaWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
     double Area =  pCell->GetCellData()->GetItem("Area");
    return Area;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellAreaWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
                                                                                        // AbstractCellPopulation<DIM,DIM>& rCellPopulation
{

    unsigned cell_id = pCell->GetCellId();    
    double Area  =  pCell->GetCellData()->GetItem("Area");
    *this->mpOutStream << " " << Area;

}

// Explicit instantiation
template class CellAreaWriter<1,1>;
template class CellAreaWriter<1,2>;
template class CellAreaWriter<2,2>;
template class CellAreaWriter<1,3>;
template class CellAreaWriter<2,3>;
template class CellAreaWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellAreaWriter)

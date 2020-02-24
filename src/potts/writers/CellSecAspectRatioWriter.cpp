
#include "CellSecAspectRatioWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellSecAspectRatioWriter<ELEMENT_DIM, SPACE_DIM>::CellSecAspectRatioWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("CellAspectRatio.dat")
{                          
    this->mVtkCellDataName = "AspectRatio";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellSecAspectRatioWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
     double AspectRatio =  pCell->GetCellData()->GetItem("AspectRatio");
    return AspectRatio;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellSecAspectRatioWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
                                                                                        // AbstractCellPopulation<DIM,DIM>& rCellPopulation
{

    unsigned cell_id = pCell->GetCellId();    
    double AspectRatio  =  pCell->GetCellData()->GetItem("AspectRatio");
    *this->mpOutStream << " " << AspectRatio;

}

// Explicit instantiation
template class CellSecAspectRatioWriter<1,1>;
template class CellSecAspectRatioWriter<1,2>;
template class CellSecAspectRatioWriter<2,2>;
template class CellSecAspectRatioWriter<1,3>;
template class CellSecAspectRatioWriter<2,3>;
template class CellSecAspectRatioWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellSecAspectRatioWriter)

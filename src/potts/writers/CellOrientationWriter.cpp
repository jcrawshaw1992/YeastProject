
#include "CellOrientationWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellOrientationWriter<ELEMENT_DIM, SPACE_DIM>::CellOrientationWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("CellArea.dat")
{
    this->mVtkCellDataName = "Area";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellOrientationWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
     double Orientation =  pCell->GetCellData()->GetItem("Orientation");
    return Orientation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellOrientationWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
                                                                                        // AbstractCellPopulation<DIM,DIM>& rCellPopulation
{

    unsigned cell_id = pCell->GetCellId();    
    double Orientation  =  pCell->GetCellData()->GetItem("Orientation");
    *this->mpOutStream << " " << Orientation;

}

// Explicit instantiation
template class CellOrientationWriter<1,1>;
template class CellOrientationWriter<1,2>;
template class CellOrientationWriter<2,2>;
template class CellOrientationWriter<1,3>;
template class CellOrientationWriter<2,3>;
template class CellOrientationWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellOrientationWriter)


#include "CellPerimeterWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPerimeterWriter<ELEMENT_DIM, SPACE_DIM>::CellPerimeterWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("CellPerimeter.dat")
{
    this->mVtkCellDataName = "Perimeter";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellPerimeterWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double Perimeter =  pCell->GetCellData()->GetItem("Perimeter");
    return Perimeter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPerimeterWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
                                                                                        // AbstractCellPopulation<DIM,DIM>& rCellPopulation
{

    unsigned cell_id = pCell->GetCellId();    
    double Perimeter  =  pCell->GetCellData()->GetItem("Perimeter");
    *this->mpOutStream << " " << Perimeter;

}

// Explicit instantiation
template class CellPerimeterWriter<1,1>;
template class CellPerimeterWriter<1,2>;
template class CellPerimeterWriter<2,2>;
template class CellPerimeterWriter<1,3>;
template class CellPerimeterWriter<2,3>;
template class CellPerimeterWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPerimeterWriter)

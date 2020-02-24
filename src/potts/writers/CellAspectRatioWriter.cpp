
#include "CellAspectRatioWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellAspectRatioWriter<ELEMENT_DIM, SPACE_DIM>::CellAspectRatioWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("CellAspectRatio.dat")
{
    this->mVtkCellDataName = "Aspect Ratio";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellAspectRatioWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
     double ApsectRatio  =  pCell->GetCellData()->GetItem("AspectRatio");
    return ApsectRatio;
}


/* I have several parameters writen into the CellParameteres.dat file 
At each time step I have 
time cell i.d AspectRatio MajorAxisAngle id AspectRatio MajorAxisAngle
*/

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellAspectRatioWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
                                                                                        // AbstractCellPopulation<DIM,DIM>& rCellPopulation
{
    unsigned cell_id = pCell->GetCellId();
    

    double ApsectRatio  =  pCell->GetCellData()->GetItem("AspectRatio");
    *this->mpOutStream << " " << cell_id << " " << ApsectRatio;
}

// Explicit instantiation
template class CellAspectRatioWriter<1,1>;
template class CellAspectRatioWriter<1,2>;
template class CellAspectRatioWriter<2,2>;
template class CellAspectRatioWriter<1,3>;
template class CellAspectRatioWriter<2,3>;
template class CellAspectRatioWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellAspectRatioWriter)

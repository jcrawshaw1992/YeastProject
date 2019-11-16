
#include "CellCenterWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellCenterWriter<ELEMENT_DIM, SPACE_DIM>::CellCenterWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("CellCenter.dat")
{
    this->mVtkCellDataName = "Center";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellCenterWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double Center =  pCell->GetCellData()->GetItem("CenterX");
    return Center;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> CellCenterWriter<ELEMENT_DIM, SPACE_DIM>::GetVectorCellDataForVtkOutput( CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    c_vector<double, SPACE_DIM> CenterPoint;//= Create_c_vector(0,0,0);
    CenterPoint[0] =  pCell->GetCellData()->GetItem("CenterX");
    CenterPoint[1] =  pCell->GetCellData()->GetItem("CenterX");
    CenterPoint[2] =  pCell->GetCellData()->GetItem("CenterX");
    return CenterPoint;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellCenterWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
                                                                                        // AbstractCellPopulation<DIM,DIM>& rCellPopulation
{
    c_vector<double, SPACE_DIM> CenterPoint;
    CenterPoint[0] =  pCell->GetCellData()->GetItem("CenterX");
    CenterPoint[1] =  pCell->GetCellData()->GetItem("CenterY");
    CenterPoint[2] =  pCell->GetCellData()->GetItem("CenterZ");

    // unsigned cell_id = pCell->GetCellId();    
    *this->mpOutStream << " " << CenterPoint[0] << " " << CenterPoint[1] << " "  << CenterPoint[2] << " "  ;

}

// Explicit instantiation
template class CellCenterWriter<1,1>;
template class CellCenterWriter<1,2>;
template class CellCenterWriter<2,2>;
template class CellCenterWriter<1,3>;
template class CellCenterWriter<2,3>;
template class CellCenterWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellCenterWriter)

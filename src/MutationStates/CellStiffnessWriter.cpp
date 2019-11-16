/*
Writes out the membrane constants :) --- Jess
*/


#include "CellStiffnessWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellLabel.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellStiffnessWriter<ELEMENT_DIM, SPACE_DIM>::CellStiffnessWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.MembraneConstants")
{
    this->mVtkCellDataName = "MembraneProperties";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellStiffnessWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double KA = pCell->GetCellData()->GetItem("AreaConstant");

    return KA;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellStiffnessWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    
    double KS =     pCell->GetCellData()->GetItem("ShearModulus"); 
    double Kalpha = pCell->GetCellData()->GetItem("AreaDilationModulus"); 
    double KA =     pCell->GetCellData()->GetItem("AreaConstant");
    
    *this->mpOutStream <<  KS << " " << Kalpha << " " << KA << " ";
}

// Explicit instantiation
template class CellStiffnessWriter<1,1>;
template class CellStiffnessWriter<1,2>;
template class CellStiffnessWriter<2,2>;
template class CellStiffnessWriter<1,3>;
template class CellStiffnessWriter<2,3>;
template class CellStiffnessWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellStiffnessWriter)

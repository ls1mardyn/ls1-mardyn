#include "AdResSCellProcessor.h"

AdResSCellProcessor::AdResSCellProcessor(const double cfRad, const double ljcRad ):CellProcessor{cfRad,ljcRad}{

}

CellProcessor* AdResSCellProcessor::GetCellProcessor(){
    return this->cell_processor;
}

void AdResSCellProcessor::initTraversal(){
    this->cell_processor->initTraversal();
}

void AdResSCellProcessor::preprocessCell(ParticleCell& cell){
    cell_processor->preprocessCell(cell);
}

void AdResSCellProcessor::processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll){
    //if is pure interaction
    cell_processor->processCellPair(cell1,cell2,sumAll);
    //else it is adressing time 
}
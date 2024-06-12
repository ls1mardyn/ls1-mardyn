#include"InteractionCellProcessor.h"
#include"PMF.h"

InteractionCellProcessor::InteractionCellProcessor(const double rc1,const double rc2):CellProcessor{rc1,rc2},internal_cell_processor{nullptr}{

}

void InteractionCellProcessor::initTraversal(){
    this->internal_cell_processor->initTraversal();
}

void InteractionCellProcessor::preprocessCell(ParticleCell& cell){
    this->internal_cell_processor->preprocessCell(cell);
}

void InteractionCellProcessor::processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll){
    this->internal_cell_processor->processCellPair(cell1,cell2,sumAll);
}

void InteractionCellProcessor::processCell(ParticleCell& cell){
    this->internal_cell_processor->processCell(cell);
}

double InteractionCellProcessor::processSingleMolecule(Molecule* m1, ParticleCell& cell2){
    return this->internal_cell_processor->processSingleMolecule(m1,cell2);
}

void InteractionCellProcessor::postprocessCell(ParticleCell& cell){
    this->internal_cell_processor->postprocessCell(cell);
}

void InteractionCellProcessor::endTraversal(){
    this->internal_cell_processor->endTraversal();
}
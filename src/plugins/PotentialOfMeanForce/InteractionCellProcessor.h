/*
 * Created on Tue Jun 11 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

class PMF;

class InteractionCellProcessor:public CellProcessor{


    private:
    CellProcessor* internal_cell_processor;//For equal component interactions
    PMF* adres;
    public:
    InteractionCellProcessor(const double cfRad, const double ljcRad);
    CellProcessor* GetCellProcessor(){
        return internal_cell_processor;
    }

    void SetCellProcessor(CellProcessor* cp){
        delete internal_cell_processor;
        internal_cell_processor = cp;
    }

    virtual void initTraversal() override;
    virtual void preprocessCell(ParticleCell& cell) override;
    virtual void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll = false) override;
    virtual void processCell(ParticleCell& cell) override;
    virtual double processSingleMolecule(Molecule* m1, ParticleCell& cell2) override;
    virtual void postprocessCell(ParticleCell& cell) override;
    virtual void endTraversal() override;

};
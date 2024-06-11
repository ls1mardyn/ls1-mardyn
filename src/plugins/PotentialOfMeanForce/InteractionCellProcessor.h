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


class InteractionCellProcessor:public CellProcessor{


    private:
    CellProcessor* cell_processor;//For equal component interactions

    public:
    InteractionCellProcessor(const double cfRad, const double ljcRad);
    CellProcessor* GetCellProcessor(){
        return cell_processor;
    }

    //functions from the interface
    //Must something particular be done here?
    virtual void initTraversal() override{}
    //Ideally does nothing?
    virtual void preprocessCell(ParticleCell& cell) override{}
    //Call original processor if pure interaction
    virtual void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll = false) override{}
    virtual void processCell(ParticleCell& cell) override{}
    virtual double processSingleMolecule(Molecule* m1, ParticleCell& cell2) override{}
    virtual void postprocessCell(ParticleCell& cell) override{}
    virtual void endTraversal() override{}

};
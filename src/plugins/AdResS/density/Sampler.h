/*
 * Created on Wed May 08 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */


//Abstract class?
#pragma once
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

#include "plugins/AdResS/util/grid3d/Grid.h"
#include "plugins/AdResS/util/grid3d/GridHandler.h"
using namespace Grid3D;

//TODO: this has to be the base class for all the samplers we work with, not an enum based profiler (since we are not profiling)
class Sampler{

    private:
    std::vector<double> sampled_data;

    public:
    /**
     * 
     * Move data member inits to constructor, call other methods within
    */
    virtual void init()=0;

    /**
     * lets settle for a tab separator, alright?
    */
    virtual void writeDensity(const std::string& filename)=0;
    std::vector<double>& GetSampledData(){
        return sampled_data;
    }
    virtual void SampleData(ParticleContainer* pc)=0;
};

class GridSampler2:public Sampler{

    private:

    Grid* grid;
    GridHandler* handler;
    double measure_radius;
    public:

    GridSampler2(Grid* grid, double rad);

    /**
     * 
     * Do resizing and other stuff
    */
    virtual void init() override;
    virtual void SampleData(ParticleContainer* pc) override;
    virtual void writeDensity(const std::string& filename) override;

    private:
    //TODO: I insist on this implementation. The grid at some point is distributed, so to get the correct nodal values, it should be based on grid parallelization. Not on globalizing the particle positions on every step. 
    void SampleAtNodes(ParticleContainer* pc);

    /**
     * \brief Check if particle within distance using Euclidean norm
    */
    bool ParticleInsideMeasuringSpace(std::array<double,3> n_pos, std::array<double,3> l_pos);
};



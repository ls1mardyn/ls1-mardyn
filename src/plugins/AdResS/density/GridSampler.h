/*
 * Created on Wed Apr 24 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

#include "plugins/AdResS/util/grid3d/Grid.h"
#include "plugins/AdResS/util/grid3d/GridHandler.h"

namespace Grid3D {
	class GridSampler{
	private:
		//TODO:split into node and element sample types
		struct Samples{
			std::vector<int> particles_per_cell;
			std::vector<double> material_density;
			std::vector<int> particles_per_node;
			double target_number_density;
		} samples;

		Grid* grid;
		double measure_radius;
		GridHandler handler;

	public:
		GridSampler() = default;
		void init(Grid* grid);
		GridHandler& GetGridHandler();
		void SetMeasureRadius(double r);
		void SampleAtNodes(const std::array<std::vector<double>,3>& positions);
		void SetSubsetMaterialDensityValues();
		void SetTargetValue();
		Samples& GetSamples(){ return samples; }
		template<class T>
		std::ostream& WritePlaneSamples(std::ostream& out, T data);
		//TODO: really has to be templatized or so, code repetition
		std::ostream& WriteSample(std::ostream& out, std::vector<double>& smpl);
		std::ostream& WriteSample(std::ostream& out, std::vector<int>& smpl);
		std::ostream& WriteInfo(std::ostream& out);
	private:
		IdxArray GetParticleLocalCellIndices(ParticleIterator it);
		bool ParticleInsideMeasuringSpace(std::array<double, 3>, std::array<double, 3> par_pos);
	};
};

template<class T>
std::ostream& Grid3D::GridSampler::WritePlaneSamples(std::ostream& out, T data){
    std::string prefix = "//[Sampler]: ";
    std::vector<int> plane_nodes = handler.GetNodesOnPlane(grid,10);
    out<<prefix+" target value: "<<this->samples.target_number_density<<"\n";
    for(int i=0;i<plane_nodes.size();i++){
        double rx = grid->GetNodes()[plane_nodes[i]].GetPosition()[0];
        double ry = grid->GetNodes()[plane_nodes[i]].GetPosition()[1];
        double val = data[plane_nodes[i]];
        out<<plane_nodes[i]<<"\t"<<rx<<"\t"<<ry<<"\t"<<val<<"\n";
    }
    return out;
}
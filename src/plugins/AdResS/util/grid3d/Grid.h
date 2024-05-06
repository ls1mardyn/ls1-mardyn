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
#include "Helper.h"
#include <iostream>
#include "plugins/AdResS/util/Subset.h"

namespace Grid3D {
	//Should handle both other classes entirely
	// Is blind wrt the Sampler or its actual use
	class Grid{
	private:
		struct NodeInfo{
			IdxArray nodes_per_dimension{1,1,1};
			int total_nodes = 1;
		} nodal_info;
		struct ElementInfo{
			IdxArray elements_per_dimension{1,1,1};
			DataArray element_width_per_dimension{0,0,0};
			int total_elements;
			int largest_index;
			double volume;
		} element_info;

		std::vector<Node> nodes;
		std::vector<Element> elements;
		DataArray lower_corner;
		DataArray upper_corner;
		ParticleContainer* particle_container;
		std::vector<Subset> subsets;

	public:

		Grid() = default;
		Grid(DataArray lower, DataArray upper, int x, int y, int z);
		~Grid() = default;
		void SetMeshLimits(DataArray lower, DataArray upper);
		void MeshAllDomain();
		void StartGrid(int els_x, int els_y, int els_z);
		void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
		std::array<double, 3>& GetUpperCorner();
		std::array<double, 3>& GetLowerCorner();
		ElementInfo& GetElementInfo(){
			return element_info;
		}
		NodeInfo& GetNodeInfo(){
			return nodal_info;
		}
		std::vector<Node>& GetNodes(){
			return this->nodes;
		}

		std::array<int, 8> GetElementGlobalNodeIndices(int el);
		void AddSubset(const std::vector<int>& list, const std::string& Name);
		Subset& GetSubset(std::string& Name);
		std::vector<Subset>& GetSubsets();
		friend std::ostream& operator<<(std::ostream& out, Grid& grid);

	private:
		void InitNodePositions();
	};
};

/*
 * Created on Mon Apr 29 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */


#pragma once
#include<vector>
#include "Grid.h"

namespace Grid3D {
	class GridHandler{
	public:
		//TODO:This method must go
		std::vector<int> GetNodesOnPlane(Grid* grid, double pz);
		std::vector<int> GetNodesOnPlanebyComponent(Grid* grid, double p, int component);
		//TODO:Have to make nodes exclusive to one subset
		void SetGridBoundarySubsets(Grid* grid);
	};
};
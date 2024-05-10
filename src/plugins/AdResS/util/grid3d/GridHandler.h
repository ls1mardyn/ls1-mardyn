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
	/**
	 * \brief Helpful methods to manipulate a grid
	*/
	class GridHandler{
	public:

		/**
		 * \brief Get a list of global node indeces based on a normal direction 
		 * 
		 * \param component must be between 0 and 2 depending on which normal direction is desired
		 * \param p is the point coordinate defining the plane
		*/
		std::vector<int> GetNodesOnPlanebyComponent(Grid* grid, double p, int component);
		//TODO: nodes must be exclusive to one subset (do they?)
		void SetGridBoundarySubsets(Grid* grid);
	};
};
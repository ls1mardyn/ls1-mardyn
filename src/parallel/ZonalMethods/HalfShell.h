/*
 * CommunicationScheme.h
 *
 *  Created on: 14.07.2017
 *      Author: sauermann
 */

#pragma once

#include "ZonalMethod.h"

/**
 * This class implements the HalfShell Method. Halo is only imported, if it lies within the half shell.
 * The half shell is defined using the same traversal scheme as the VectorizedCellProcessor uses for determining,
 * whether or whether not to calculate the macroscopic values, i.e: cellIndex1 < cellIndex2.
 * As LinkedCells defines the cellIndex using z as the major array index (then y, then x), we will do the same here:
 * if z<z_domain, we don't import the region. if z>z_domain we import it.
 * if z within z_domain: y determines the procedure.
 * The HalfShell method requires forces to be exchanged.
 */
class HalfShell: public ZonalMethod {
public:
	HalfShell(){}
	virtual ~HalfShell(){}

	virtual std::vector<HaloRegion> getHaloImportForceExportRegions(HaloRegion& initialRegion, double cutoffRadius,
			bool coversWholeDomain[3]) override {
		const std::function<bool(const int[3])> condition = [](const int d[3])->bool {
			int pseudoCellIndex = ((d[2] * 2) + d[1]) * 2 + d[2];
			return pseudoCellIndex > 0;
		};
		return getHaloRegionsConditional(initialRegion, cutoffRadius, coversWholeDomain, condition);
	}
};


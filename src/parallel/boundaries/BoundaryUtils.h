
#include <vector>
#include <string>
#include <tuple>
#include <array>

#include "BoundaryType.h"
#include "DimensionType.h"

#include "molecules/Molecule.h"

namespace BoundaryUtils
{
	const std::array<std::string, 6> permissibleDimensionsString = {"+x", "+y", "+z", "-x", "-y", "-z"};
	const std::array<int, 6> permissibleDimensionsInt = {-1, -2, -3, 1, 2, 3};

	bool checkIfDimensionStringPermissible(std::string dimension);
	bool checkIfDimensionNumericPermissible(int dim);

	DimensionType convertStringToDimension(std::string dimension);
	DimensionType convertNumericToDimension(int dim);
	DimensionType convertLS1DimsToDimensionPos(int dim);
	std::vector<DimensionType> convertHaloOffsetToDimensionVector(int* offset);

	std::string convertDimensionToString(DimensionType dimension);
	std::string convertDimensionToStringAbs(DimensionType dimension);
	int convertDimensionToNumeric(DimensionType dimension);
	int convertDimensionToNumericAbs(DimensionType dimension);
	int convertDimensionToLS1Dims(DimensionType dimension);

	BoundaryType convertStringToBoundary(std::string boundary);

	std::tuple<std::array<double,3>, std::array<double,3>> getInnerBuffer(const std::array<double,3> givenRegionBegin, const std::array<double,3> givenRegionEnd, DimensionType dimension, double regionWidth);
	bool isMoleculeLeaving(const Molecule molecule, const std::array<double,3> regionBegin, const std::array<double,3> regionEnd, DimensionType dimension, double timestepLength, double nextStepVelAdjustment);
	std::tuple<std::array<double,3>, std::array<double,3>> getOuterBuffer(const std::array<double,3> givenRegionBegin, const std::array<double,3> givenRegionEnd, DimensionType dimension, double* regionWidth);

	inline int findSign(int n) { return n < 0 ? -1 : 1; }
	inline int findSign(DimensionType dimension) { return findSign(convertDimensionToNumeric(dimension)); }
}
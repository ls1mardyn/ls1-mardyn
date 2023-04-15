
#include <vector>
#include <string>

#include "BoundaryType.h"
#include "DimensionType.h"

#include "molecules/Molecule.h"

namespace BoundaryUtils
{
	const std::vector<std::string> permissibleDimensionsString = {"+x", "+y", "+z", "-x", "-y", "-z"};
	const std::vector<int> permissibleDimensionsInt = {-1, -2, -3, 1, 2, 3};

	bool checkIfDimensionStringPermissible(std::string dimension);
	bool checkIfDimensionNumericPermissible(int dim);

	DimensionType convertStringToDimension(std::string dimension);
	DimensionType convertNumericToDimension(int dim);

	std::string convertDimensionToString(DimensionType dimension);
	std::string convertDimensionToStringAbs(DimensionType dimension);
	int convertDimensionToNumeric(DimensionType dimension);
	int convertDimensionToNumericAbs(DimensionType dimension);
	int convertDimensionToLS1Dims(DimensionType dimension);

	BoundaryType convertStringToBoundary(std::string boundary);

	bool getInnerBuffer(double* givenRegionBegin, double* givenRegionEnd, DimensionType dimension, double regionWidth, double* returnRegionBegin, double* returnRegionEnd);
	bool isMoleculeLeaving(Molecule molecule, double* regionBegin, double* regionEnd, DimensionType dimension, double timestepLength);
	bool getOuterBuffer(double* givenRegionBegin, double* givenRegionEnd, DimensionType dimension, double regionWidth, double* returnRegionBegin, double* returnRegionEnd);

	inline int findSign(int n) { return n < 0 ? -1 : 1; }
	inline int findSign(DimensionType dimension) { return findSign(convertDimensionToNumeric(dimension)); }
}
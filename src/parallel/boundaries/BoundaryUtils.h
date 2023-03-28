
#include <vector>
#include <string>

#include "BoundaryType.h"
#include "DimensionType.h"

#include "molecules/Molecule.h"

namespace BoundaryUtils
{
	static const std::vector<std::string> permissibleDimensionsString = {"+x", "+y", "+z", "-x", "-y", "-z"};
	static const std::vector<int> permissibleDimensionsInt = {-1, -2, -3, 1, 2, 3};

	static bool checkIfDimensionStringPermissible(std::string dimension);
	static bool checkIfDimensionNumericPermissible(int dim);

	static DimensionType convertStringToDimension(std::string dimension);
	static DimensionType convertNumericToDimension(int dim);

	static std::string convertDimensionToString(DimensionType dimension);
	static std::string convertDimensionToStringAbs(DimensionType dimension);
	static int convertDimensionToNumeric(DimensionType dimension);
	static int convertDimensionToNumericAbs(DimensionType dimension);

	

	static int findSign(int n) { return n < 0 ? -1 : 1; }
}
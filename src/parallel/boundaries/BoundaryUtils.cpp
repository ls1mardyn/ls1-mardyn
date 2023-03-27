#include <algorithm>
#include <string>


#include "BoundaryUtils.h"
#include "utils/Logger.h"

bool BoundaryUtils::checkIfDimensionStringPermissible(std::string dimension)
{
	return std::find(permissibleDimensionsString.begin(),permissibleDimensionsString.end(),dimension) == permissibleDimensionsString.end();
}

bool BoundaryUtils::checkIfDimensionNumericPermissible(int dim)
{ 
	return (dim >= -3 && dim <= 3 && dim != 0); 
}

DimensionType BoundaryUtils::convertStringToDimension(std::string dimension) {
	if(checkIfDimensionStringPermissible(dimension))
	{
		Log::global_log->error() << "Invalid dimension passed for enum conversion" << std::endl;
		return DimensionType::ERROR;
	}
	if(dimension == "+x")
		return DimensionType::POSX;
	if(dimension == "+y")
		return DimensionType::POSY;
	if(dimension == "+z")
		return DimensionType::POSZ;
	if(dimension == "-x")
		return DimensionType::NEGX;
	if(dimension == "-y")
		return DimensionType::NEGY;
	//if(dimension == "-z")
	return DimensionType::NEGZ;
}

DimensionType BoundaryUtils::convertNumericToDimension(int dim)
{
	if(checkIfDimensionNumericPermissible(dim))
	{
		Log::global_log->error() << "Invalid dimension passed for enum conversion" << std::endl;
		return DimensionType::ERROR;
	}
	switch(findSign(dim))
	{
		case -1:
			switch(dim)
			{
				case 1:	return DimensionType::NEGX;
				case 2: return DimensionType::NEGY;
				default: return DimensionType::NEGZ; //case 3
				
			}
		case 1:
			switch(dim)
			{
				case 1: return DimensionType::POSX;
				case 2: return DimensionType::POSY;
				default: return DimensionType::POSZ; //case 3
			}
	}
}

std::string BoundaryUtils::convertDimensionToString(DimensionType dimension) 
{ 
	switch (dimension)
	{
	case DimensionType::POSX:
		return "+x";
	
	case DimensionType::POSY:
		return "+y";

	case DimensionType::POSZ:
		return "+z";

	case DimensionType::NEGX:
		return "-x";
	
	case DimensionType::NEGY:
		return "-y";
	
	case DimensionType::NEGZ:
		return "-z";

	default:
		Log::global_log->error() << "Invalid dimension passed for enum conversion" << std::endl;
		return "error";
	}
}

int BoundaryUtils::convertDimensionToNumeric(DimensionType dimension)
{
	switch (dimension)
	{
	case DimensionType::POSX:
		return 1;
	
	case DimensionType::POSY:
		return 2;

	case DimensionType::POSZ:
		return 3;

	case DimensionType::NEGX:
		return -1;
	
	case DimensionType::NEGY:
		return -2;
	
	case DimensionType::NEGZ:
		return -3;

	default:
		Log::global_log->error() << "Invalid dimension passed for enum conversion" << std::endl;
		return 0;
	}
}
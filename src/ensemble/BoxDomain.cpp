#include "BoxDomain.h"

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

#include "Simulation.h"

using namespace std;
using Log::global_log;

enum BoxConfig
{
	BLANK = 0,
	L = 1, //001b
	BOXMIN = 2, //010b
	L_BOXMIN = 3, //011b
	BOXMAX = 4, //100b
	L_BOXMAX = 5, //101b
	BOXMIN_BOXMAX = 6, //110b
	L_BOXMIN_BOXMAX = 7 //111b
};

BoxDomain::BoxDomain() {
	for( int d = 0; d < 3; d++)
		_rmin[d] = 0;
}

void BoxDomain::readXML(XMLfileUnits& xmlconfig) {
	std::array<double, 3> l;
	int tempConfig = 0;

	bool l1, l2, l3, lExists;
	l1 = xmlconfig.getNodeValueReduced("lx",l[0]);
	l2 = xmlconfig.getNodeValueReduced("ly",l[1]);
	l3 = xmlconfig.getNodeValueReduced("lz",l[2]);
	lExists = l1 & l2 & l3;
	if(lExists)
		tempConfig |= 1;

	if(xmlconfig.changecurrentnode("boxMin"))
	{
		xmlconfig.getNodeValueReduced("x", _rmin[0]);
		xmlconfig.getNodeValueReduced("y", _rmin[1]);
		xmlconfig.getNodeValueReduced("z", _rmin[2]);
		xmlconfig.changecurrentnode("..");
		tempConfig |= 2;
	}
	if(xmlconfig.changecurrentnode("boxMax"))
	{
		xmlconfig.getNodeValueReduced("x", _rmax[0]);
		xmlconfig.getNodeValueReduced("y", _rmax[1]);
		xmlconfig.getNodeValueReduced("z", _rmax[2]);
		xmlconfig.changecurrentnode("..");
		tempConfig |= 4;
	}

	BoxConfig boxConfig = static_cast<BoxConfig>(tempConfig); //value will be within range, so not undefined
	switch(boxConfig)
	{
		case BLANK:
			global_log->error() << "Domain specification does not have enough info to define box dimensions. Please make sure the domain specification is correct." << endl;
			Simulation::exit(1);
		
		case L:					// old behaviour
			for(int d = 0; d < 3; d++)
			{
				_rmin[d] = 0;
				_rmax[d] = l[d];
			}
			break;
		
		case L_BOXMIN:
		case L_BOXMIN_BOXMAX: 	// overwriting _rmax
			for(int d = 0; d < 3; d++)
				_rmax[d] = _rmin[d] + l[d];
			break;

		case L_BOXMAX:
			for(int d = 0; d < 3; d++)
				_rmin[d] = _rmax[d] - l[d];
			break;

		case BOXMAX:			// constructor already set _rmin to zero, so nothing to do here
		case BOXMIN_BOXMAX: 	// both values have been read and set, nothing to calculate
			break;

		default:				// not enough data
			global_log->error() << "Domain specification does not have enough info to define box dimensions. Please make sure the domain specification is correct." << endl;
			Simulation::exit(1);
	}

	global_log->info() << "Box lower corner (x,y,z): " << _rmin[0] << "," << _rmin[1] << "," << _rmin[2] << endl;
	global_log->info() << "Box upper corner (x,y,z): " << _rmax[0] << "," << _rmax[1] << "," << _rmax[2] << endl;
}

double BoxDomain::V() {
	return (_rmax[0] - _rmin[0]) * (_rmax[1] - _rmin[1]) * (_rmax[2] - _rmin[2]);
}

void BoxDomain::setLength(int d, double l) {
	_rmax[d] = l + _rmin[d];
}

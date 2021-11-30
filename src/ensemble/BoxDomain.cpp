#include "BoxDomain.h"

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

#include "Simulation.h"

using namespace std;
using Log::global_log;

BoxDomain::BoxDomain() {
	for( int d = 0; d < 3; d++)
		_rmin[d] = 0;
}

void BoxDomain::readXML(XMLfileUnits& xmlconfig) {
	double l[3];
	bool l1, l2, l3, lExists;
	l1 = xmlconfig.getNodeValueReduced("lx",l[0]);
	l2 = xmlconfig.getNodeValueReduced("ly",l[1]);
	l3 = xmlconfig.getNodeValueReduced("lz",l[2]);
	lExists = l1 & l2 & l3;
	bool boxMinExists, boxMaxExists;
	if(boxMinExists = xmlconfig.changecurrentnode("boxMin"))
	{
		xmlconfig.getNodeValueReduced("x", _rmin[0]);
		xmlconfig.getNodeValueReduced("y", _rmin[1]);
		xmlconfig.getNodeValueReduced("z", _rmin[2]);
		xmlconfig.changecurrentnode("..");
	}
	if(boxMaxExists = xmlconfig.changecurrentnode("boxMax"))
	{
		xmlconfig.getNodeValueReduced("x", _rmax[0]);
		xmlconfig.getNodeValueReduced("y", _rmax[1]);
		xmlconfig.getNodeValueReduced("z", _rmax[2]);
		xmlconfig.changecurrentnode("..");
	}

	//calculations for missing elements
	//earlier implementation did not check for existence of l values, but this implementation is checking for it now
	if(!lExists && !boxMinExists && !boxMaxExists)
	{
		global_log->error() << "Domain specification defines neither bounds not lengths" << endl;
		Simulation::exit(1);
	}
	else if(lExists && !boxMinExists && !boxMaxExists) //old behavior
	{
		for(int d = 0; d < 3; d++)
		{
			_rmin[d] = 0;
			_rmax[d] = l[d];
		}
	}
	else if(boxMinExists && lExists) //does not matter if boxMax exists, overwrite anyway
	{
		for(int d = 0; d < 3; d++)
			_rmax[d] = _rmin[d] + l[d];
	}
	else if(boxMaxExists && lExists) //'all existing' case covered by earlier case
	{
		for(int d = 0; d < 3; d++)
			_rmin[d] = _rmax[d] - l[d];
	}
	else if(boxMaxExists && boxMinExists) //do nothing
	{}
	else //only boxmin or boxmax exists, overwrite
	{
		global_log->error() << "Domain specification does not have enough info to define box dimensions" << endl;
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

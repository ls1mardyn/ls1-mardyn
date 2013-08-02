#include "molecules/Site.h"
#include "Simulation.h"

void LJcenter::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValueReduced("coords/x", _r[0]);
	xmlconfig.getNodeValueReduced("coords/y", _r[1]);
	xmlconfig.getNodeValueReduced("coords/z", _r[2]);
	xmlconfig.getNodeValueReduced("mass", _m);
	xmlconfig.getNodeValueReduced("epsilon", _eps);
	xmlconfig.getNodeValueReduced("sigma", _sigma);
	_rc = _simulation.getLJCutoff();
	xmlconfig.getNodeValueReduced("cutoff", _rc);
	_uLJshift6 = 0.0;
	xmlconfig.getNodeValueReduced("shifted", _uLJshift6);
}

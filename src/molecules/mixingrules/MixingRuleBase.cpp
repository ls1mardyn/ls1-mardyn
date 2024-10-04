#include "MixingRuleBase.h"

#include <algorithm>
#include <ostream>

#include "utils/Logger.h"
#include "utils/mardyn_assert.h"
#include "utils/xmlfileUnits.h"

void MixingRuleBase::readXML(const XMLfileUnits& xmlconfig) {
	int cid1;
	int cid2;
	xmlconfig.getNodeValue("@cid1", cid1);
	xmlconfig.getNodeValue("@cid2", cid2);

    // catch invalid inputs
	if (cid1 == cid2) {
		std::ostringstream error_message;
		error_message << "Mixing rules: cid1 and cid2 must not be the same but are both " << cid1 << std::endl;
		MARDYN_EXIT(error_message);
	} else if (std::min(cid1, cid2) < 0) {
		std::ostringstream error_message;
		error_message << "Mixing rules: cid1 and cid2 must be greater than zero" << std::endl;
		MARDYN_EXIT(error_message);
	} 
	
	// Symmetry for mixing rules is assumed
	// _cid1 must be less than _cid2
	std::tie(_cid1, _cid2) = std::minmax(cid1, cid2);
	// component id - 1 to convert to internal format starting with 0
	--_cid1;
	--_cid2;
}

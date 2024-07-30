#include "MixingRuleBase.h"

#include <ostream>

#include "utils/Logger.h"
#include "utils/mardyn_assert.h"
#include "utils/xmlfileUnits.h"

void MixingRuleBase::readXML(const XMLfileUnits& xmlconfig) {
	int cid1;
	int cid2;
	xmlconfig.getNodeValue("@cid1", cid1);
	xmlconfig.getNodeValue("@cid2", cid2);

	if (cid1 == cid2) {
		Log::global_log->error() << "Mixing rules: cid1 and cid2 must not be the same but are both " << cid1 << std::endl;
		mardyn_exit(1);
	} else if ((cid1 <= 0) or (cid2 <= 0)) {
		Log::global_log->error() << "Mixing rules: cid1 and cid2 must be greater than zero" << std::endl;
		mardyn_exit(1);
	} else if (cid1 > cid2) {
	// Symmetry for mixing rules is assumed
	// cid1 must be less than cid2
		_cid1 = cid2 - 1;  // component id - 1 to convert to internal format starting with 0
		_cid2 = cid1 - 1;
	} else {
		_cid1 = cid1 - 1;
		_cid2 = cid2 - 1;
	}
}

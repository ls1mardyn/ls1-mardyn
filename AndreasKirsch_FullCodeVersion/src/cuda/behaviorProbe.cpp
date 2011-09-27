#include "behaviorProbe.h"

int BehaviorProbe::counter;
BehaviorProbe::Type BehaviorProbe::lastType;
const char *BehaviorProbe::msg[BPT_COUNT] = {
		"", "FM cleared", "FM used", "FM set"
};

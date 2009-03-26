#include "utils/Watch.h"

#include <sstream>

utils::Watch::Watch(const std::string& className, const std::string& operationName):
  _log( className ),
  _operationName( operationName ),
  _startClockTicks(std::clock()) {
}

utils::Watch::~Watch() {
 clock_t totalClockTicks = clock() - _startClockTicks;
 
 std::ostringstream message;
 message << "total number of clock ticks within block: " << totalClockTicks
         << " (" << double(totalClockTicks)/CLOCKS_PER_SEC << " sec)";
 _log.info( _operationName, message.str() );
}

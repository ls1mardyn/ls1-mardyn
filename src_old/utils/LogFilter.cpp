#include "utils/LogFilter.h"

#ifdef Parallel
#include "parallel/Node.h"
#endif

//configuration::LogFilterConfiguration::BlackList utils::LogFilter::_blacklist;

// void utils::LogFilter::configureLogFilter( const configuration::LogFilterConfiguration& configuration ) {
//   _blacklist = configuration.getBlacklist();
// }

std::string utils::LogFilter::extractNamespace(const std::string& className) {
  int separator = className.find( "::" );
  return className.substr(0,separator);
}

utils::LogFilter::LogFilter() {
}

utils::LogFilter::~LogFilter() {
}

bool utils::LogFilter::writeDebug(const std::string& className) {
//   for (configuration::LogFilterConfiguration::BlackList::const_iterator p = _blacklist.begin(); p!=_blacklist.end(); p++ ) {
//   	#ifdef Parallel
//   	if ( (*p)._namespaceName == extractNamespace(className) ) {
//   	  if (
//   	   ( (*p)._rank == -1 ) ||
//   	   ( (*p)._rank == parallel::Node::getInstance().getRank() )
//   	  ) {
//   	    return false;
//   	  }
//   	}
//   	#else
//   	if ( (*p)._namespaceName == extractNamespace(className) ) {
//   	  return false;
//   	}
//   	#endif
//   }
  return true;
}

bool utils::LogFilter::writeInfo(const std::string& className) {
  return true;
}

bool utils::LogFilter::writeWarning(const std::string& className) {
  return true;
}

bool utils::LogFilter::writeError(const std::string& className) {
  return true;
}

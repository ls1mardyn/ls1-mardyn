#ifndef _UTILS_LOGFILTER_H_
#define _UTILS_LOGFILTER_H_

//#include "configuration/LogFilterConfiguration.h"

#include <set>
#include <string>

namespace utils {
  class LogFilter;
}

/**
 * Filters the debug messages according to a black list.
 *
 * @version $Revision: 1.2 $
 * @author  Tobias Weinzierl
 */
class utils::LogFilter {
  private:
  //static configuration::LogFilterConfiguration::BlackList _blacklist;
    
    static std::string extractNamespace(const std::string& className);

    LogFilter();
    virtual ~LogFilter();
  public:
    //static void configureLogFilter( const configuration::LogFilterConfiguration& configuration );
    static bool writeDebug(const std::string& className);
    static bool writeInfo(const std::string& className);
    static bool writeWarning(const std::string& className);
    static bool writeError(const std::string& className);
};

#endif

#ifndef _UTILS_WATCH_H_
#define _UTILS_WATCH_H_

#include <ctime>
#include "utils/Log.h"

namespace utils {
  class Watch;
}


/**
 * A simple class that has to be included to measure the clock ticks required
 * for an operation. To use it you've to create an instance of this class at
 * the beginning of the operation block you want to measure the time spent
 * within it:
 *
 * <pre>
 *   void anyOperation() {
 *     utils::Watch watch("MyClass","anyOperation()");
 *     ...
 *   }
 * </pre>
 *
 * The result the of the measurement is written to log.info level. Note that
 * this operation works within operation blocks (e.g. a for-loop), too.
 *
 * @version $Revision: 1.3 $
 * @author  Tobias Weinzierl
 */
class utils::Watch {
  private:
    /**
     * Log device the result is written to.
     */
    utils::Log     _log;

    /**
     * Stores the name of the operation the watch is used within.
     */
    std::string    _operationName;

    /**
     * Holds the clock ticks at the beginning of the time measurement.
     */
    std::clock_t   _startClockTicks;
  public:
    /**
     * Construct a watch and start measuring time immediatly.
     *
     * @param className     Name of the class the watch is used within.
     * @param operationName Name of the operation the watch is used within.
     */
    Watch(const std::string& className, const std::string& operationName);

    /**
     * Stops the watch and plots the time spent.
     */
    virtual ~Watch();
};

#endif

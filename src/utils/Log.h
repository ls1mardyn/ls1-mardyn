#ifndef _UTILS_LOG_H_
#define _UTILS_LOG_H_

#include <iostream>
#include "utils/CommandLineLogger.h"

namespace utils {
  class Log;
}

/**
 * Log is the class all logging classes should use. To use the logging api they
 * have to create an instance by their own. It is suggested to hold this
 * instance static for the constructor of the Log class has to be given the
 * class name of the logging class. The Log class itself is stateless.
 * The log requests on this instance are processed here and forwarded to the
 * assigned logger (an internal attribute).
 *
 * Which concrete implementation has to be used for logging is switched using
 * a compiler attribute. Since the logging is used extremly often, this is
 * better than dynamic binding.
 *
 * There are four different log levels, the user may write any output:
 *
 * - error: Here only errors have to be written to. Error messages may never
 *          be oppressed.
 * - warning:
 * - debug: Debug information that is switched off normally.
 * - info:  Statistical information, copyright and similar information. Should
 *          be used rather seldom.
 *
 * @version $Revision: 1.3 $
 * @author  Tobias Weinzierl
 */
class utils::Log {
  private:
    /**
     * Name of the class that is using the interface.
     */
    std::string _className;

    static utils::CommandLineLogger _logger;
  public:
    /**
     * Constructor.
     *
     * @param className Name of the class that is using the logging component.
     *                  Please specify both class name and namespace using the
     *                  format namespace::classname
     */
    Log(const std::string& className);
    virtual ~Log();

   /**
     * Log a debug information. The method has to be given the method name. Often it
     * is useful to pass the whole method signature, that means e.g. if the
     * method info itself would log, the method invocation would look like this:
     * <pre>
     *   info("info(std::string,std::string)", "just an example text");
     * </pre>
     * So you are able to identify methods in your log although you use
     * overloading.
     *
     * @param methodName method name
     * @param message    log message
     */
    void debug(const std::string& methodName, const std::string& message) const;

    /**
     * Log an information. The method has to be given the method name. Often it
     * is useful to pass the whole method signature, that means e.g. if the
     * method info itself would log, the method invocation would look like this:
     * <pre>
     *   info("info(std::string,std::string)", "just an example text");
     * </pre>
     * So you are able to identify methods in your log although you use
     * overloading.
     *
     * @param methodName method name
     * @param message    log message
     */
    void info(const std::string& methodName, const std::string& message) const;

    /**
     * Log a warning. The method has to be given the method name. Often it
     * is useful to pass the whole method signature, that means e.g. if the
     * method info itself would log, the method invocation would look like this:
     * <pre>
     *   info("info(std::string,std::string)", "just an example text");
     * </pre>
     * So you are able to identify methods in your log although you use
     * overloading.
     *
     * @param methodName method name
     * @param message    log message
     */
    void warning(const std::string& methodName, const std::string& message) const;

    /**
     * Log an error. The method has to be given the method name. Often it
     * is useful to pass the whole method signature, that means e.g. if the
     * method info itself would log, the method invocation would look like this:
     * <pre>
     *   info("info(std::string,std::string)", "just an example text");
     * </pre>
     * So you are able to identify methods in your log although you use
     * overloading.
     *
     * @param methodName method name
     * @param message    log message
     */
    void error(const std::string& methodName, const std::string& message) const;
};


#endif

#ifndef _COMMANDLINELOGGER_H_
#define _COMMANDLINELOGGER_H_

#include <iostream>

namespace utils {
  class CommandLineLogger;
}

/**
 * Standard log output device. This device is configured using some compiler
 * switches that define the type of information written. Usually there's only
 * one instance of this class for only the Log type aggregates it in a static
 * way.
 *
 * @version $Revision: 1.6 $
 * @author  Tobias Weinzierl
 */
class utils::CommandLineLogger {
  private:
    /**
     * output stream
     */
    std::ostream& _out;

    /**
     * Flag indicating wheather the device should write information
     * about the computer the information is written by.
     *
     * Defines weather the machine information should be written. This flag is
     * useful for cluster applications.
     */
    #ifdef LogMachine
    static const bool WRITE_MACHINE_INFORMATION = true;
    #else
    static const bool WRITE_MACHINE_INFORMATION = false;
    #endif

    /**
     * Flag indicating wheather the device should write information
     * about the computer the information is written by.
     *
     * Augment every log entry with a time stamp.
     */
    #ifdef LogTime
    static const bool WRITE_TIME_STAMP = true;
    #else
    static const bool WRITE_TIME_STAMP = false;
    #endif

    /**
     * Flag indicating wheather the device should write information
     * about the computer the information is written by.
     *
     * Write a trace information giving you the class and line the log is
     * written from.
     */
    #ifdef LogTrace
    static const bool WRITE_TRACE_INFORMATION = true;
    #else
    static const bool WRITE_TRACE_INFORMATION = false;
    #endif

    /**
     * Flag indicating wheather the device should write information
     * about the computer the information is written by.
     *
     * Add a information on the log level if one writes a log entry.
     */
    #ifdef LogMessageType
    static const bool WRITE_LOG_MESSAGE_TYPE = true;
    #else
    static const bool WRITE_LOG_MESSAGE_TYPE = false;
    #endif

    /**
     * If this flag is set, the complete logging is switched off.
     */
    #ifdef LogOff
    static const bool WRITE_LOG_INFORMATION = false;
    #else
    static const bool WRITE_LOG_INFORMATION = true;
    #endif

    /**
     * Flag indicating wheather the device should write information
     * about the computer the information is written by.
     *
     * If the Debug flag is set, all debug information is written. Otherwise
     * the log interface won't plot your debug outputs.
     */
    #ifdef Debug
    static const bool LOG_DEBUG_INFO = true;
    #else
    static const bool LOG_DEBUG_INFO = false;
    #endif

    /**
     * Writes information about the computer the output is written from.
     * The information string contains the operation system name, the computer
     * name and the cpu name.
     */
    void writeMachineInformation();

    /**
     * Writes a timestamp to the standard output.
     */
    void writeTimeStamp();

    void writeTraceInformation( const std::string& className, const std::string& methodName );

    void writeLogMessageType( const std::string& type );

    /**
     * Declared private since assignment does not make sense for an output
     * class (output information mismatch).
     */
    CommandLineLogger& operator=(const CommandLineLogger& rhs);

    /**
     * Declared private since copying does not make sense for an output
     * class (output information mismatch).
     */
    CommandLineLogger(const CommandLineLogger& param);
  public:
   /**
     * Create a LineLogger.
     *
     * @param out output stream logger should write to
     */
    CommandLineLogger( std::ostream& out = std::cout );

    ~CommandLineLogger();

    void debug(const std::string& className, const std::string& methodName, const std::string& message);
    void info(const std::string& className, const std::string& methodName, const std::string& message);
    void warning(const std::string& className, const std::string& methodName, const std::string& message);
    void error(const std::string& className, const std::string& methodName, const std::string& message);
};

#endif

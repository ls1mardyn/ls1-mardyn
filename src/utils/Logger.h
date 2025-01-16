#ifndef LOGGER_H_
#define LOGGER_H_

#define USE_GETTIMEOFDAY

#include <chrono>
#include <ctime>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>

#ifdef ENABLE_MPI
#include <mpi.h>

#ifndef NDEBUG
/** When NDEBUG macro is undefined check the expression return value to be MPI_SUCCESS.
 *
 * Check expression to return MPI_SUCCESS. Can only be used after MPI_Init is called
 * because the Logger cannot be initialized before with the current implementation.
 */
#define MPI_CHECK(x) do {                   \
    int __ret;                              \
    if (MPI_SUCCESS != (__ret = (x)))       \
    Log::global_log->error() << "MPI returned with error code " << __ret << std::endl;  \
} while (0)

#else
#define MPI_CHECK(x) x
#endif

#endif  /* ENABLE_MPI */


/* we use a separate namespace because we have some global definitions for
 * the log level */
namespace Log {

class Logger;

/**
 * Global logger variable for use in the entire program.
 * Must be initialized with constructor
 * Namespace visibility:
 */
#ifndef LOGGER_SRC
extern std::unique_ptr<Log::Logger> global_log;
#endif

/**
 * list of available log levels
 * For each level a name has to be specified in the constructor of the Logger() class.
 * This name will be prepended later on to the log message. */
typedef enum {
	None    = 0,  /* supress output */
	Fatal   = 1,  /* program exit */
	Error   = 2,  /* program corrected */
	Warning = 4,  /* perhaps wrong */
	Info    = 8,  /* user info */
	Debug   = 16, /* detailed info for debugging */
	All     = 32,
} logLevel;

/** @brief The Logger class provides a simple interface to handle log messages.
 *
 * The logger provides an easy way to output logging messages with various log levels.
 * For writing log messages use fatal(), error(), warning(), info() or debug() as
 * with normal streams, e.g.
 * > log.error() << "Wrong parameter." << std::endl;
 * Please include std::endl statements at the end of output as they will flush
 * the stream buffers.
 * If ENABLE_MPI is enabled logger initialization has to take place after the
 * MPI_Init call.
 */
class Logger {
private:
	logLevel _log_level;
	logLevel _msg_log_level;
	bool _do_output;
	std::shared_ptr<std::ostream> _log_stream;
	std::map<logLevel, std::string> logLevelNames;

	std::chrono::system_clock::time_point _starttime;
	int _rank;

	/// initialize the list of log levels with the corresponding short names
	void init_log_levels() {
		logLevelNames.insert(std::pair<logLevel, std::string>(None,    "NONE"));
		logLevelNames.insert(std::pair<logLevel, std::string>(Fatal,   "FATAL ERROR"));
		logLevelNames.insert(std::pair<logLevel, std::string>(Error,   "ERROR"      ));
		logLevelNames.insert(std::pair<logLevel, std::string>(Warning, "WARNING"    ));
		logLevelNames.insert(std::pair<logLevel, std::string>(Info,    "INFO"       ));
		logLevelNames.insert(std::pair<logLevel, std::string>(Debug,   "DEBUG"      ));
		logLevelNames.insert(std::pair<logLevel, std::string>(All,     "ALL"      ));
	}

public:
	/**
	 * Constructor for a logger to a stream.
	 *
	 * Initializes the log level, log stream and the list of log level names.
	 * If ENABLE_MPI is enabled by default, all process perform logging output.
	 * Note: The default stream used (std::cout) cannot be deleted. Therefore the
	 * passed shared pointer to it uses a no-op deleter function.
	 */
	Logger(logLevel level = Log::Error, std::shared_ptr<std::ostream> os = std::shared_ptr<std::ostream>(&std::cout, [](void*){ /* no-op */}));

	~Logger() = default;

	/// General output template for variables, strings, etc.
	template<typename T>
	Logger &operator<<(T const& t) {
		if (_msg_log_level <= _log_level && _do_output)
			*_log_stream << t;
		return *this;
	}

	/* Specialized versions for manipulators.  */
	// e.g. endl
	Logger &operator<<(std::ostream& (*f)(std::ostream&)) {
		if (_msg_log_level <= _log_level && _do_output)
			*_log_stream << f;
		return *this;
	}
	// e.g. hex.
	Logger &operator<<(std::ios_base& (*f)(std::ios_base&)) {
		if (_msg_log_level <= _log_level && _do_output)
			f(*_log_stream);
		return *this;
	}
	template<class Ch, class Tr>
	Logger &operator<<(std::basic_ios<Ch, Tr>& (*f)(std::basic_ios<Ch, Tr>&)) {
		if (_msg_log_level <= _log_level && _do_output)
			f(*_log_stream);
		return *this;
	}

	/// Add log info in front of messages
	Logger& msg_level(logLevel level) {
				_msg_log_level = level;
		if (_msg_log_level <= _log_level && _do_output) {
			// Include timestamp
			const auto now = std::chrono::system_clock::now();
			const auto time_since_start = now - _starttime;

			const auto now_time_t = std::chrono::system_clock::to_time_t(now);
			std::tm now_local{};
			localtime_r(&now_time_t, &now_local);
			*_log_stream << logLevelNames[level] << ":\t" << std::put_time(&now_local, "%Y-%m-%dT%H:%M:%S") << " ";
			*_log_stream << std::setw(8) << std::chrono::duration<double>(time_since_start).count() << " ";

			*_log_stream << "[" << _rank << "]\t";
		}
		return *this;
	}

	/// shorthand versions for easy usage of the different log levels
	Logger& fatal() {
		return msg_level(Fatal);
	}
	Logger& error() {
		return msg_level(Error);
	}
	Logger& error_always_output() {
		_do_output = true;
		return msg_level(Error);
	}
	Logger& warning() {
		return msg_level(Warning);
	}
	Logger& info() {
		return msg_level(Info);
	}
	Logger& debug() {
		return msg_level(Debug);
	}

	/// set log level
	logLevel set_log_level(logLevel l) {
		_log_level = l;
		return _log_level;
	}
	/// set log level from string
	logLevel set_log_level(std::string l) {
		for (const auto& [lvl, name] : logLevelNames) {
			if (name == l) {
				_log_level = lvl;
				return _log_level;
			}
		}
		error_always_output() << "Logger::set_log_level() unknown argument '" << l
							  << "'. Falling back to output everything!" << std::endl;
		return Log::All;
	}

	/// return log level
	logLevel get_log_level() {
		return _log_level;
	}

	void set_log_stream(std::shared_ptr<std::ostream> os) {
		_log_stream = os;
	}

	/// switch on / off output
	bool set_do_output(bool val) {
		return _do_output = val;
	}

	bool get_do_output() {
		return _do_output;
	}

	/// initialize starting time
	void init_starting_time() {
		_starttime = std::chrono::system_clock::now();
	}


}; /* end of class Logger */
} /* end of namespace */

#endif /*LOGGER_H_*/

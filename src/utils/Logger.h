#ifndef LOGGER_H_
#define LOGGER_H_

#define USE_GETTIMEOFDAY

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <map>
#include <ctime>
#include <string>
#include <sstream>
#include <chrono>

#ifdef USE_GETTIMEOFDAY
#include <sys/time.h>
#endif

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
 * Gobal logger variable for use in the entire program.
 * Must be initialized with constructor e.g. new Log::Logger().
 * Namespace visibility:
 *    */
#ifndef LOGGER_SRC
extern Log::Logger *global_log;
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
 * Provides easy interface to handle log messages. Initialize either with
 * output level and stream or output level and filename or use default constructor
 * values (Error, &(std::cout)). With a given file basename and MPI Support each rank will
 * create and write to his own file.
 * For writing log messages use fatal(), error(), warning(), info() or debug() as
 * with normal streams, e.g.
 * > log.error() << "Wrong parameter." << std::endl;
 * For easy handling of output within MPI applications there are the following methods:
 * set_mpi_output_root(int root)
 * set_mpi_output_rall()
 * set_mpi_output_ranks(int num_ranks, int * ranks)
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
	std::string _filename;
	std::ostream *_log_stream;
	std::map<logLevel, std::string> logLevelNames;
#ifdef USE_GETTIMEOFDAY
	timeval _starttime;
#else
	time_t _starttime;
#endif

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

	// don't allow copy-construction
	Logger(const Logger&) : _log_level(Log::Error), _msg_log_level(Log::Error), _do_output(true),
			_filename(""), _log_stream(0), logLevelNames(), _starttime(), _rank(0)
	{ }

	// don't allow assignment
	Logger& operator=(const Logger&) { return *this; }

public:
	/** Initializes the log level, log stream and the list of log level names.
	 * If ENABLE_MPI is enabled by default all process perform logging output. */
	Logger(logLevel level = Log::Error, std::ostream *os = &(std::cout));

	Logger(logLevel level, std::string prefix);

	/// Destructor flushes stream
	~Logger();

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
			const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
			tm unused{};
			const auto* lt = localtime_r(&now, &unused);
			//*_log_stream << ctime(&t) << " ";
			std::stringstream timestampstream;
			// maybe sprintf is easier here...
			timestampstream << std::setfill('0') << std::setw(4) << (1900 + lt->tm_year) << std::setw(2) << (1 + lt->tm_mon) << std::setw(2) << lt->tm_mday << "T" << std::setw(2) << lt->tm_hour << std::setw(2) << lt->tm_min << std::setw(2) << lt->tm_sec;
			*_log_stream << logLevelNames[level] << ":\t" << timestampstream.str() << " ";
			//timestampstream.str(""); timestampstream.clear();
#ifdef USE_GETTIMEOFDAY
			timeval tod;
			gettimeofday(&tod, 0);
			*_log_stream << std::setw(8) << tod.tv_sec - _starttime.tv_sec + (tod.tv_usec - _starttime.tv_usec) / 1.E6 << " ";
#else
			*_log_stream << t-_starttime << "\t";
#endif

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
		// identify the loglevel by comparing the first char of the names
		for (const auto& [lvl, name] : logLevelNames) {
			// case-insensitive matching
			if (std::toupper(name[0]) == std::toupper(l[0])) {
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

	/// switch on / off output
	bool set_do_output(bool val) {
		return _do_output = val;
	}

	bool get_do_output() {
		return _do_output;
	}

	/// initialize starting time
	void init_starting_time() {
#ifdef USE_GETTIMEOFDAY
		gettimeofday(&_starttime, 0);
#else
		_starttime = time(NULL);
#endif
	}

	/* methods for easy handling of output processes */

	/// allow logging only for a single process
	void set_mpi_output_root(int root = 0);

	/// all processes shall perform logging
	void set_mpi_output_all();

	/// allow a set of processes for logging
	bool set_mpi_output_ranks(int num_nums, int* nums);



}; /* end of class Logger */
} /* end of namespace */

#endif /*LOGGER_H_*/

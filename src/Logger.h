#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>

/* some applications define PARALLEL to enable MPI */
#ifdef PARALLEL
#define MPI_SUPPORT
#endif

#ifdef MPI_SUPPORT
#include <mpi.h>
#endif

/* we use a seperate namespace because we have some global definitions for
 * the log level */
namespace Log {

	class Logger;

/* Gobal logger variable for use in the entire program.
 * Must be initialized with constructor e.g. new Log::Logger().
 * Namespace visibility:
 *	 using Log::global_log;
 */
#ifndef LOGGER_SRC
	extern Logger *global_log;
#endif

/* list of available log levels
 * For each level a name has to be specified in the Logger::init_levels() method.
 * This name will be prepended to the log message. */
typedef enum {
	LLNULL	= 0,	/* (LogLevel)Null to supress output */
	FATAL	= 1,	/* program exit */
	ERROR	= 2,	/* program corrected */
	WARNING = 4,	/* perhaps wrong */
	INFO	= 8,	/* user info */
	DEBUG	= 16,	/* detailed info for debugging */
	ALL
} logLevel;

/* Logging class 
 * Provides easy interface to handle log messages. Initialize either with 
 * output level and stream or output level and filename or use default constructor
 * values (WARNING, cout). With a given file basename and MPI Support each rank will
 * write to his own file.
 * For writing log messages use stream reference returned from fatal(), error(),
 * warning(), info() or debug() as normal, e.g. 
 * > log.error() << "Wrong parameter." << endl; 
 * For easy handling of output within MPI applications there  are
 * log.set_mpi_output_root(int root)
 * log.set_mpi_output_ranks(int num_ranks, int * ranks)
 * Please include std::endl statements at the end of output as they will flush
 * the stream buffers.
 * If MPI_SUPPORT is enabled logger initialization has to take place after the 
 * MPI_Init call.
 */
class Logger {
private:
	std::ostream* _logStream;	   /* output stream */
	std::ostream* _logNullStream;  /* dummy stream */
	logLevel _logLevel;			   /* log level up to which we output */
	std::map <logLevel, std::string> logLevelNames; /* list with log level names */
	bool _no_output;  /* true if no log messages shall be printed */
	int _comm_rank;   /* rank in MPI_COMM_WORLD; 0 in non MPI case */

	void init(logLevel maxlevel, std::ostream* logStream);
	std::ostream& msg_level(logLevel level);

public:
	Logger (logLevel maxlevel = WARNING, std::ostream* logStream = &std::cout);

	Logger (logLevel maxlevel, std::string filename);

	logLevel set_log_level(logLevel maxlevel) {
		return _logLevel  = maxlevel;
	}

	logLevel get_log_level() {
		return _logLevel;
	}

	~Logger();

	std::ostream& fatal()	{ return msg_level(FATAL);	 }
	std::ostream& error()	{ return msg_level(ERROR);	 }
	std::ostream& warning() { return msg_level(WARNING); }
	std::ostream& info()	{ return msg_level(INFO);	 }
	std::ostream& debug()	{ return msg_level(DEBUG);	 }

	/* switch on / off output */
	bool set_no_output(bool val) {
		return _no_output = val;
	}

	/* methods for easy handling of output ranks */
	bool set_mpi_output_root(int root = 0);
	bool set_mpi_output_ranks(int num_nums, int* nums);
	/* parameters:	nums[num_nums]	array containing the numbers of ranks which
	 * shall print messages 
	 */

}; /* end of class Logger */

} /* end of namespace */
#endif

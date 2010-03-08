#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
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
  extern Log::Logger *global_log;
#endif

  /* list of available log levels
   * For each level a name has to be specified in the constructor of the Logger() class.
   * This name will be prepended later on to the log message. */
  typedef enum {
    NONE        = 0,    /* supress output */
    FATAL       = 1,    /* program exit */
    ERROR       = 2,    /* program corrected */
    WARNING     = 4,    /* perhaps wrong */
    INFO        = 8,    /* user info */
    MARDYN_DEBUG = 16,   /* detailed info for debugging */
    ALL
  } logLevel;

  /* Logging class 
   * Provides easy interface to handle log messages. Initialize either with 
   * output level and stream or output level and filename or use default constructor
   * values (WARNING, cout). With a given file basename and MPI Support each rank will
   * write to his own file.
   * For writing log messages use fatal(), error(), warning(), info() or debug() as 
   * with normal streams, e.g. 
   * > log.error() << "Wrong parameter." << endl; 
   * For easy handling of output within MPI applications there  are
   * log.set_mpi_output_root(int root)
   * log.set_mpi_output_rall()
   * log.set_mpi_output_ranks(int num_ranks, int * ranks)
   * Please include std::endl statements at the end of output as they will flush
   * the stream buffers.
   * If MPI_SUPPORT is enabled logger initialization has to take place after the 
   * MPI_Init call.
   */
  class Logger{
    public:

      /* Initializes the log level, log stream and the list of log level names.
       * If MPI_SUPPORT is enabled by default all process perform logging output. */
      Logger( logLevel level = Log::ERROR, std::ostream& os = std::cout ): _log_level(level), _log_stream(os) { 
	_do_output  = true;
	logLevelNames.insert(std::pair<logLevel, std::string>(FATAL,   "FATAL ERROR"  ));
	logLevelNames.insert(std::pair<logLevel, std::string>(ERROR,   "ERROR"        ));
	logLevelNames.insert(std::pair<logLevel, std::string>(WARNING, "WARNING"      ));
	logLevelNames.insert(std::pair<logLevel, std::string>(INFO,    "INFO"         ));
	logLevelNames.insert(std::pair<logLevel, std::string>(MARDYN_DEBUG, "DEBUG"        ));
#ifdef MPI_SUPPORT
	MPI_Comm_rank( MPI_COMM_WORLD, &_rank );
#endif
      }

      /* Destructor flushes stream */
      ~Logger(){ 
	_log_stream << std::flush;
      }

      /* General output template for variables, strings, etc. */
      template<typename T> 
	Logger &operator<<(T const& t) { 
	  if ( _msg_log_level <= _log_level && _do_output )
	    _log_stream << t; 
	  return *this;
	}

      /* Specialized versions for manipulators.  */
      // e.g. endl 
      Logger &operator<<(std::ostream& (*f)(std::ostream&)) {
	if ( _msg_log_level <= _log_level && _do_output )
	  _log_stream <<  f;
	return *this; 
      }
      // e.g. hex.
      Logger &operator<<(std::ios_base& (*f)(std::ios_base&)) {
	if ( _msg_log_level <= _log_level && _do_output )
	  f(_log_stream);
	return *this;
      }
      template <class Ch, class Tr >
	Logger &operator<<(std::basic_ios<Ch,Tr>& (*f)(std::basic_ios<Ch,Tr>&)) {
	  if ( _msg_log_level <= _log_level && _do_output )
	    f(_log_stream);
	  return *this;
	}

      /* Add log info in front of messages */
      Logger& msg_level( logLevel level ){
	_msg_log_level = level;
	if ( _msg_log_level <= _log_level && _do_output ) {
#ifdef MPI_SUPPORT
	  _log_stream << "[rank " << _rank << "] ";
#endif
	  _log_stream << logLevelNames[level] << ": ";
	}
	return *this;
      }

      /* shorthand versions for easy usage of the different log levels */
      Logger& fatal()  { return msg_level(FATAL);   }
      Logger& error()  { return msg_level(ERROR);   }
      Logger& warning(){ return msg_level(WARNING); }
      Logger& info()   { return msg_level(INFO);    }
      Logger& debug()  { return msg_level(MARDYN_DEBUG);   }


      // set log level
      logLevel set_log_level( logLevel l ){
	_log_level = l;
	return _log_level;
      }
      // return log level
      logLevel get_log_level(){
	return _log_level;
      }

      /* switch on / off output */
      bool set_do_output(bool val) {
	return _do_output = val;
      }

#ifdef MPI_SUPPORT
      /* methods for easy handling of output processes */

      /* allow logging only for a single process */
      bool set_mpi_output_root(int root = 0) {
	if (_rank != root)
	  _do_output = false;
	return _do_output;
      }

      /* all processes shall perform logging */
      bool set_mpi_output_all() {
	_do_output = true;
	return _do_output;
      }

      /* allow a set of processes for logging */
      bool set_mpi_output_ranks(int num_nums, int* nums) {
	int i;
	for(i = 0; i < num_nums; i++)
	  if (nums[i] == _rank)
	    _do_output = true;
	return _do_output;
      }
#endif


    private:        
      std::ostream& _log_stream;
      logLevel _log_level;
      logLevel _msg_log_level;
      bool _do_output;
#ifdef MPI_SUPPORT
      int _rank;
#endif

      std::map<logLevel, std::string> logLevelNames;
  }; /* end of class Logger */
} /* end of namespace */

#endif

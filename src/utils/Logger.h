#ifndef LOGGER_H
#define LOGGER_H

#define USE_GETTIMEOFDAY

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <map>
#include <ctime>
#include <string>
#include <sstream>

#ifdef USE_GETTIMEOFDAY
#include <sys/time.h>
#endif

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

	/** Gobal logger variable for use in the entire program.
	 * Must be initialized with constructor e.g. new Log::Logger().
	 * Namespace visibility:
	 *	 using Log::global_log;
	 */
#ifndef LOGGER_SRC
	extern Log::Logger *global_log;
#endif

	/** list of available log levels
	 * For each level a name has to be specified in the constructor of the Logger() class.
	 * This name will be prepended later on to the log message. */
	typedef enum {
		None        = 0,    /* supress output */
		Fatal       = 1,    /* program exit */
		Error       = 2,    /* program corrected */
		Warning     = 4,    /* perhaps wrong */
		Info        = 8,    /* user info */
		Debug       = 16,   /* detailed info for debugging */
		All
	} logLevel;

	/** Logging class 
	 * Provides easy interface to handle log messages. Initialize either with 
	 * output level and stream or output level and filename or use default constructor
	 * values (Error, &(std::cout)). With a given file basename and MPI Support each rank will
	 * create and write to his own file.
	 * For writing log messages use fatal(), error(), warning(), info() or debug() as 
	 * with normal streams, e.g. 
	 * > log.error() << "Wrong parameter." << endl; 
	 * For easy handling of output within MPI applications there are the following methods:
	 * set_mpi_output_root(int root)
	 * set_mpi_output_rall()
	 * set_mpi_output_ranks(int num_ranks, int * ranks)
	 * Please include std::endl statements at the end of output as they will flush
	 * the stream buffers.
	 * If MPI_SUPPORT is enabled logger initialization has to take place after the 
	 * MPI_Init call.
	 */
	class Logger{
		private:
			std::string _filename;
			std::ostream *_log_stream;
			logLevel _log_level;
			logLevel _msg_log_level;
			bool _do_output;
			std::map<logLevel, std::string> logLevelNames;
#ifdef USE_GETTIMEOFDAY
			timeval _starttime;
#else
			time_t _starttime;
#endif

#ifdef MPI_SUPPORT
			int _rank;
#endif

			/// initilaize the list of log levels with the corresponding short names
			void init_log_levels() {
				logLevelNames.insert(std::pair<logLevel, std::string>(Fatal,   "FATAL ERROR"  ));
				logLevelNames.insert(std::pair<logLevel, std::string>(Error,   "ERROR"        ));
				logLevelNames.insert(std::pair<logLevel, std::string>(Warning, "WARNING"      ));
				logLevelNames.insert(std::pair<logLevel, std::string>(Info,    "INFO"         ));
				logLevelNames.insert(std::pair<logLevel, std::string>(Debug,   "DEBUG"        ));
			}


		public:
			/** Initializes the log level, log stream and the list of log level names.
			 * If MPI_SUPPORT is enabled by default all process perform logging output. */
			Logger( logLevel level = Log::Error, std::ostream *os = &(std::cout) ): _log_level(level), _log_stream(os) {
				init_starting_time();
				_filename = "";
				_do_output  = true;
				this->init_log_levels();
#ifdef MPI_SUPPORT
				MPI_Comm_rank( MPI_COMM_WORLD, &_rank );
#endif
			}

			Logger( logLevel level, std::string prefix ): _log_level(level) {
				init_starting_time();
				std::stringstream filenamestream;
				filenamestream << prefix;
				_do_output  = true;
#ifdef MPI_SUPPORT
				MPI_Comm_rank( MPI_COMM_WORLD, &_rank );
				filenamestream << "_R" << _rank;
#endif
				filenamestream << ".log";
				_filename = filenamestream.str();
				_log_stream = new std::ofstream( _filename.c_str() );
			}

			/// Destructor flushes stream
			~Logger(){ 
				*_log_stream << std::flush;
				if ( _filename != "" )
					(static_cast<std::ofstream*>(_log_stream))->close();
			}

			/// General output template for variables, strings, etc.
			template<typename T> 
				Logger &operator<<(T const& t) { 
					if ( _msg_log_level <= _log_level && _do_output )
						*_log_stream << t; 
					return *this;
				}

			/* Specialized versions for manipulators.  */
			// e.g. endl 
			Logger &operator<<(std::ostream& (*f)(std::ostream&)) {
				if ( _msg_log_level <= _log_level && _do_output )
					*_log_stream <<  f;
				return *this; 
			}
			// e.g. hex.
			Logger &operator<<(std::ios_base& (*f)(std::ios_base&)) {
				if ( _msg_log_level <= _log_level && _do_output )
					f(*_log_stream);
				return *this;
			}
			template <class Ch, class Tr >
				Logger &operator<<(std::basic_ios<Ch,Tr>& (*f)(std::basic_ios<Ch,Tr>&)) {
					if ( _msg_log_level <= _log_level && _do_output )
						f(*_log_stream);
					return *this;
				}

			/// Add log info in front of messages
			Logger& msg_level( logLevel level ){
				using namespace std;
				_msg_log_level = level;
				if ( _msg_log_level <= _log_level && _do_output ) {
					// Include timestamp
					time_t t;
					t = time( NULL );
					tm* lt = localtime(&t);
					//*_log_stream << ctime(&t) << " ";
					stringstream timestampstream;
					// maybe sprintf is easier here...
					timestampstream << setfill('0') << setw(4) << (1900 + lt->tm_year) 
						<< setw(2) << (1 + lt->tm_mon) << setw(2) << lt->tm_mday << "T" 
						<< setw(2) << lt->tm_hour << setw(2) << lt->tm_min << setw(2) << lt->tm_sec;
					*_log_stream << timestampstream.str() << " ";
					//timestampstream.str(""); timestampstream.clear();
#ifdef USE_GETTIMEOFDAY
					timeval tod;
					gettimeofday(&tod, 0);
					*_log_stream << setw(8) << tod.tv_sec-_starttime.tv_sec+(tod.tv_usec-_starttime.tv_usec)/1.E6 << " ";
#else
					*_log_stream << t-_starttime << "\t";
#endif
#ifdef MPI_SUPPORT
					*_log_stream << "[" << _rank << "]\t";
#endif
					*_log_stream << logLevelNames[level] << ": ";
				}
				return *this;
			}

			/// shorthand versions for easy usage of the different log levels
			Logger& fatal()  { return msg_level(Fatal);   }
			Logger& error()  { return msg_level(Error);   }
			Logger& warning(){ return msg_level(Warning); }
			Logger& info()   { return msg_level(Info);    }
			Logger& debug()  { return msg_level(Debug);   }


			/// set log level
			logLevel set_log_level( logLevel l ){
				_log_level = l;
				return _log_level;
			}
			/// return log level
			logLevel get_log_level(){
				return _log_level;
			}

			/// switch on / off output
			bool set_do_output(bool val) {
				return _do_output = val;
			}

			/// initialize starting time
			void init_starting_time(){
#ifdef USE_GETTIMEOFDAY
				gettimeofday(&_starttime, 0);
#else
				_starttime = time( NULL );
#endif
			}

#ifdef MPI_SUPPORT
			/* methods for easy handling of output processes */

			/// allow logging only for a single process
			bool set_mpi_output_root(int root = 0) {
				if (_rank != root)
					_do_output = false;
				return _do_output;
			}

			/// all processes shall perform logging
			bool set_mpi_output_all() {
				_do_output = true;
				return _do_output;
			}

			/// allow a set of processes for logging
			bool set_mpi_output_ranks(int num_nums, int* nums) {
				int i;
				for(i = 0; i < num_nums; i++)
					if (nums[i] == _rank)
						_do_output = true;
				return _do_output;
			}
#endif
	}; /* end of class Logger */
} /* end of namespace */

#endif

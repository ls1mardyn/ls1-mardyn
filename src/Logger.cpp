#define LOGGER_SRC

#include "Logger.h"

namespace Log {

Logger *global_log;

/* Nullstream to send supressed output
 * We build a buffer returning no badbit() and use this for the
 * stream initialization */
typedef struct nullstream: std::ostream {
	struct nullbuf: std::streambuf {
		int overflow(int c) { return traits_type::not_eof(c); }
	} m_sbuf;
	nullstream(): std::ios(&m_sbuf), std::ostream(&m_sbuf) {}
} Nullstream; 


/* initialization called by constructor
 * Initializes the log level, log stream, nullstream and the 
 * list of log level names.
 * If MPI_SUPPORT is enabled all processes print log messages by default. */
void Logger::init(logLevel maxlevel, std::ostream* logStream) {
	_logLevel  = maxlevel;
	_logStream = logStream;
	_logNullStream = new Nullstream();
	_no_output = false;
	logLevelNames.insert(std::pair<logLevel, std::string>(FATAL,   "FATAL ERROR"));
	logLevelNames.insert(std::pair<logLevel, std::string>(ERROR,   "ERROR"		));
	logLevelNames.insert(std::pair<logLevel, std::string>(WARNING, "WARNING"	));
	logLevelNames.insert(std::pair<logLevel, std::string>(INFO,    "INFO"		));
	logLevelNames.insert(std::pair<logLevel, std::string>(DEBUG,   "DEBUG"		));
#ifdef MPI_SUPPORT
	MPI_Comm_rank(MPI_COMM_WORLD, &_comm_rank);
#else
	_comm_rank = 0;
#endif
}

/* helper function to simplify the following fatal, error, etc. methods */
std::ostream& Logger::msg_level(logLevel level){
	if(_no_output) return *_logNullStream;
	if(level <= _logLevel) {
#ifdef MPI_SUPPORT
		*_logStream << "[rank " << _comm_rank << "] ";
#endif
		*_logStream << logLevelNames[level] << ": ";
		return *_logStream;
	}
	else 
		return *_logNullStream;
}

Logger::Logger (Log::logLevel maxlevel, std::ostream* logStream) {
	init(maxlevel, logStream);
}

Logger::Logger (logLevel maxlevel, std::string filename) {
#ifdef MPI_SUPPORT
	char rank_str[16];
	/* _comm_rank set in init, but we need it earlier in this case. */
	MPI_Comm_rank(MPI_COMM_WORLD, &_comm_rank);
	sprintf(rank_str, "%d", _comm_rank);
	filename = filename + "_" + rank_str;
#endif
	filename = filename + ".log";
	std::ofstream *logStream = new std::ofstream(filename.c_str(), std::ios::out);
	this->init(maxlevel, logStream);
}

	bool Logger::set_mpi_output_root(int root) {
		if (_comm_rank != root) 
			_no_output = true;
		return _no_output;
	}

bool Logger::set_mpi_output_ranks(int num_nums, int* nums) {
	int i;
	for(i = 0; i < num_nums; i++)
		if (nums[i] == _comm_rank)
			_no_output = false;
	return _no_output;
}

Logger::~Logger(){
	_logStream->flush();
	/* workaround: endl forces buffer flush in cases flushing does not work */
	*_logStream << std::endl; 
}

} /* end namespace Log */

#define LOGGER_SRC

#include "Logger.h"

namespace Log {

std::unique_ptr<Logger> global_log;

// Write to stream
Logger::Logger(logLevel level, std::shared_ptr<std::ostream> os) :
	_log_level(level), _msg_log_level(Log::Error),
	_do_output(true),
	_log_stream(os),
	logLevelNames(), _starttime(), _rank(0)
{
	init_starting_time();
	this->init_log_levels();
#ifdef ENABLE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif
	*_log_stream << std::boolalpha;  // Print boolean as true/false
}

} /* end namespace Log */

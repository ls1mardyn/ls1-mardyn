#define LOGGER_SRC

#include "Logger.h"

namespace Log {

std::unique_ptr<Logger> global_log;

// Write to stream
Logger::Logger(logLevel level, std::shared_ptr<std::ostream> os) :
	_log_level(level), _msg_log_level(Log::Error),
	_do_output(true),
	_log_stream(os),
	logLevelNames(), _starttime(), _msg_prefix()
{
	this->init_starting_time();
	this->init_log_levels();
	*_log_stream << std::boolalpha;  // Print boolean as true/false
}

} /* end namespace Log */

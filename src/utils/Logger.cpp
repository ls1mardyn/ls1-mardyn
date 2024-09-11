#define LOGGER_SRC

#include "Logger.h"

#include <memory>

namespace Log {

std::unique_ptr<Logger> global_log;

// Write to stream
Logger::Logger(logLevel level, std::ostream *os) :
	_log_level(level), _msg_log_level(Log::Error),
	_do_output(true), _filename(""),
	// std::cout is managed globally,
	// so do nothing when _log_stream goes out of scope
	_log_stream(os, [](std::ostream*){/* no-op deleter */}),
	logLevelNames(), _starttime(), _rank(0)
{
	init_starting_time();
	this->init_log_levels();
#ifdef ENABLE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif
	*_log_stream << std::boolalpha;  // Print boolean as true/false
}

// Write to file
Logger::Logger(logLevel level, std::string prefix) :
	_log_level(level), _msg_log_level(Log::Error),
	_do_output(true), _filename(""), _log_stream(nullptr),
	logLevelNames(), _starttime(), _rank(0)
{
	init_starting_time();
	this->init_log_levels();
	std::stringstream filenamestream;
	filenamestream << prefix;
#ifdef ENABLE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
	filenamestream << "_R" << _rank;
#endif
	filenamestream << ".log";
	_filename = filenamestream.str();
	
	_log_stream = std::make_shared<std::ofstream>(_filename.c_str());
	*_log_stream << std::boolalpha;  // Print boolean as true/false
}

/// allow logging only for a single process
void Logger::set_mpi_output_root(int root) {
	if (_rank != root)
		_do_output = false;
	else
		_do_output = true;
}

/// all processes shall perform logging
void Logger::set_mpi_output_all() {
	_do_output = true;
}

/// allow a set of processes for logging
bool Logger::set_mpi_output_ranks(int num_nums, int* nums) {
	int i;
	for(i = 0; i < num_nums; i++)
		if (nums[i] == _rank)
			_do_output = true;
	return _do_output;
}

} /* end namespace Log */

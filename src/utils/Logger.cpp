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

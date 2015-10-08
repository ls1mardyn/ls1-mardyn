#ifndef TIMER_H_
#define TIMER_H_

#include <iostream>

/* We use MPIs Wtime in parallel application, else clock */
#ifdef ENABLE_MPI 
#include <mpi.h>
#else
#include <sys/time.h>
#endif

/* PAPI hardware performance counter support */
/*
 * Take care when using multiple timers because calls to the constructor 
 * reset the hw counters to zero! 
 */
#if WITH_PAPI
#include <papi.h>

#ifndef NDEBUG
#define PAPI_CHECK(BOOL,MSG) do { \
			if ((BOOL)) { \
				std::cerr << "PAPI_ERROR: " << MSG << std::endl; \
			} \
		} while (0);
#else /* NDEBUG */
#define PAPI_CHECK(BOOL,MSG) (BOOL);
#endif /* NDEBUG */

#endif /* WITH_PAPI */

typedef enum {
	TIMER_HALTED  = 0,
	TIMER_RUNNING = 1
} timer_state;

//! @brief This class is used to messure times in sequential and parallel versions
//! @author Christoph Niethammer
class Timer {
	double _start;      // stop time
	double _stop;       // start time
	double _etime;      // elapsed time
	timer_state _state; // timer state
	bool _synced;       // timer should be synced at start and end accross processes/threads

#if WITH_PAPI
	long long *_papi_start;
	long long *_papi_stop;
	long long *_papi_counter;
	int _papi_num_counters;
	int _papi_num_avail_counters;
	int _papi_EventSet;
	bool _collect_papi;
#endif /* WITH_PAPI */

public:
	Timer() : 
		_start(0), _stop(0), _etime(0), _state(TIMER_HALTED), _synced(false)
#if WITH_PAPI
		, _papi_start(0), _papi_stop(0), _papi_counter(0), _papi_num_counters(0), _papi_num_avail_counters(0), _papi_EventSet(0), _collect_papi(false)
#endif /* WITH_PAPI */
	{
		reset();
	}

	~Timer() {
#if WITH_PAPI
		delete[] _papi_start;
		delete[] _papi_stop;
		delete[] _papi_counter;
#endif /* WITH_PAPI */
	}

	void start() {
		_stop = 0.;
		if (_state == TIMER_HALTED) {
			_start = timer();
			_state = TIMER_RUNNING;
#if WITH_PAPI
			if(_collect_papi) {
				PAPI_CHECK( (PAPI_read( _papi_EventSet, _papi_start)), "Failed reading counters.");
			}
#endif /* WITH_PAPI */
		}
		else
			std::cerr << "WARNING: Timer already running" << std::endl;
	}

	void stop() {
		if (_state == TIMER_RUNNING) {
			_stop = timer();
			_state = TIMER_HALTED;
			_etime += _stop - _start;
#if WITH_PAPI
			if(_collect_papi) {
				PAPI_CHECK((PAPI_read( _papi_EventSet, _papi_stop)), "Failed reading counters.");
				for (int i = 0; i < _papi_num_counters; i++) {
					_papi_counter[i] += _papi_stop[i] - _papi_start[i];
				}
			}
#endif /* WITH_PAPI */
		}
		else
			std::cerr << "WARNING: Timer not running" << std::endl;
	}

	void reset() {
		_state = TIMER_HALTED;
		_start = _stop = _etime = 0.;
	}

	double get_start() {
		return _start;
	}
	double get_end() {
		return _stop;
	}
	double get_etime() {
		return _etime;
	}
	timer_state get_state() {
		return _state;
	}

/** @brief Synchronize counter (accross MPI processes) */
	void set_sync(bool sync) {
		_synced = sync;
	}

#if WITH_PAPI
	int add_papi_counters(int n, char *papi_event_list[]) {
		if(_collect_papi) {
			std::cerr << "PAPI ERROR: PAPI counters already started." << std::endl;
		}
		_papi_num_avail_counters = PAPI_num_counters();
		if(_papi_num_avail_counters < 0) {
			std::cerr << "PAPI ERROR: This machine does not provide hardware counters.";
		}
		else {
			_papi_start = new long long[_papi_num_avail_counters];
			_papi_stop = new long long[_papi_num_avail_counters];
			_papi_counter = new long long[_papi_num_avail_counters];
			_papi_num_counters = n;
			_papi_EventSet = PAPI_NULL;
			PAPI_CHECK( (PAPI_create_eventset(&_papi_EventSet) != PAPI_OK), "Failed creating event set.");
			if (_papi_num_avail_counters < _papi_num_counters) {
				std::cerr << "PAPI WARNING: Not enough hw counter available. Skipping counters " << _papi_num_avail_counters << " - " << _papi_num_counters << std::endl;
				_papi_num_counters = _papi_num_avail_counters;
			}
			for (int i = 0; i < _papi_num_counters; i++) {
#ifndef NDEBUG
				std::cerr << "PAPI INFO: adding HW Counter [" << i << "]  " << papi_event_list[i] << std::endl;
#endif
				PAPI_CHECK( (PAPI_add_named_event(_papi_EventSet, papi_event_list[i]) != PAPI_OK), "Could not add counter to event set.");
				_papi_counter[i] = 0;
			}
			PAPI_CHECK( (PAPI_start(_papi_EventSet) != PAPI_OK), "Could not start PAPI counters.");
		}
		_collect_papi = true;
		return _papi_num_counters;
	}
	/* get number of used papi_counters */
	int get_papi_num_counters() {
		return _papi_num_counters;
	}
	/* get counter value between stop and last start */
	long long get_papi_counter(int index) {
		return (_papi_counter[index]);
	}
	long long get_global_papi_counter(int index) {
		long long counter = _papi_counter[index];
#if ENABLE_MPI
		MPI_Allreduce(MPI_IN_PLACE, &counter, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif /* ENABLE_MPI */
		return counter;
	}
#endif /* WITH_PAPI */

private:
	double timer() {
		double time;
#ifdef ENABLE_MPI
		if (_synced)
			MPI_Barrier(MPI_COMM_WORLD);
		time = MPI_Wtime();
#else
		struct timeval tmp_time;
		gettimeofday(&tmp_time, NULL);
		time = (1.0e6 * (double) tmp_time.tv_sec + (double) tmp_time.tv_usec) / 1.0e6;
#endif
		return time;
	}
};
#endif /*TIMER_H_*/

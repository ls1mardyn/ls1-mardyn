#ifndef TIMER_H_
#define TIMER_H_

#include <assert.h>
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
#ifdef WITH_PAPI
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


class PAPI_Initializer{
public:
	PAPI_Initializer() {
		static bool initialized = false;
		if (!initialized) {
			initialized = true;
			PAPI_CHECK( PAPI_library_init(PAPI_VER_CURRENT) != PAPI_OK, "Failed initializing PAPI library." );
			PAPI_CHECK( PAPI_thread_init(pthread_self) != PAPI_OK, "Failed initializing PAPI thread support." );
			this->initialized = true;
		}
	}

	bool papi_initialized() { return this->initialized; }

private:
	bool initialized = false;
};
#endif /* WITH_PAPI */

typedef enum {
	TIMER_HALTED  = 0,
	TIMER_RUNNING = 1
} timer_state;

//! @brief This class is used to measure times in sequential and parallel versions
//! @author Christoph Niethammer
class Timer {
	double _start;      // stop time
	double _stop;       // start time
	double _etime;      // elapsed time
	timer_state _state; // timer state
	bool _synced;       // timer should be synced at start and end across processes/threads

#ifdef WITH_PAPI
	long long *_papi_start;
	long long *_papi_stop;
	long long *_papi_counter;
	int _papi_num_counters;
	int _papi_num_avail_counters;
	int _papi_EventSet;
	bool _collect_papi;
#endif /* WITH_PAPI */

private:
	bool _active;		// timer can be active or not; if not active, then all function calls will have no effect on the timer

public:
	Timer() : 
		_start(0), _stop(0), _etime(0), _state(TIMER_HALTED), _synced(false), _active(true)
#ifdef WITH_PAPI
		, _papi_start(0), _papi_stop(0), _papi_counter(0), _papi_num_counters(0), _papi_num_avail_counters(0), _papi_EventSet(0), _collect_papi(false)
#endif /* WITH_PAPI */
	{
#ifdef WITH_PAPI
		// static dummy used to make PAPI initialization thread-safe (according C++11)
		static PAPI_Initializer pi;
		if (!pi.papi_initialized()) {
			std::cerr << "PAPI not initialized!!!1"  << std::endl;
		}
#endif
		reset();
	}

	~Timer() {
#ifdef WITH_PAPI
		delete[] _papi_start;
		delete[] _papi_stop;
		delete[] _papi_counter;
#endif /* WITH_PAPI */
	}
	void start();
	void stop();
	void reset();
	
	double get_start() {
		return _start;
	}
	double get_end() {
		return _stop;
	}
	double get_etime() {
		return _etime;
	}
	double get_etime_running() {
		return timer() - _start;
	}
	timer_state get_state() {
		return _state;
	}

	void activateTimer(){ _active = true; }
	void deactivateTimer(){ _active = false; }
	bool isActive(){ return _active; }

/** @brief Synchronize counter (across MPI processes) */
	void set_sync(bool sync) {
		_synced = sync;
	}

#ifdef WITH_PAPI
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

/*
 * Christoph Niethammer
 */
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
#ifdef PAPI
#include <papi.h>

struct Papi_Event {
	int event_code;
	char event_name[256];
};


/* Defines the list of used PAPI_counters */
/******************************************/
struct Papi_Event papi_event_list[] = {
//{ PAPI_TOT_CYC, "PAPI_TOT_CYC" },
	{ PAPI_TOT_INS, "PAPI_TOT_INS" },
//{ PAPI_VEC_DP,  "PAPI_VEC_DP"  },
	{ PAPI_L2_DCM,  "PAPI_L2_DCM"  },
	{ PAPI_L2_ICM,  "PAPI_L2_ICM"  },
	{ PAPI_L1_ICM,  "PAPI_L1_ICM"  }
//{ PAPI_DP_OPS,  "PAPI_DP_OPS"  }
};
/******************************************/

#endif /* PAPI */

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
#ifdef ENABLE_MPI
	bool _synced;       // timer syncs processes
#endif /* ENABLE_MPI */
#ifdef PAPI
	long long *_papi_start;
	long long *_papi_stop;
	int _papi_num_counters;
	int _papi_num_avail_counters;
	int _papi_EventSet;

#ifdef DEBUG
#define PAPI_CHECK(BOOL,MSG) do { \
			if ((BOOL)) { \
				std::cerr << "PAPI_ERROR: " << MSG << std::endl; \
			} \
		} while (0);
#else /* DEBUG */
#define PAPI_CHECK(BOOL,MSG) do { (BOOL); } while (0);
#endif /* DEBUG */
#endif /* PAPI */

public:
	Timer() {
		reset();
#ifdef ENABLE_MPI
		set_sync(true);
#endif /* ENABLE_MPI */
#ifdef PAPI
		if ( (_papi_num_avail_counters = PAPI_num_counters()) < 0 ) {
#ifdef DEBUG
			std::cerr << "PAPI ERROR: This machine does not provide hardware counters.";
#endif /* DEBUG */
		}
		else {
#ifdef DEBUG
			std::cout << "PAPI INFO: " << _papi_num_avail_counters << " hw counters available" << std::endl;
#endif /* DEBUG */
			_papi_start = new long long[_papi_num_avail_counters];
			_papi_stop = new long long[_papi_num_avail_counters];
			_papi_EventSet = PAPI_NULL;
			PAPI_CHECK((PAPI_create_eventset(&_papi_EventSet) != PAPI_OK), "Failed creating event set.");
			_papi_num_counters = sizeof(papi_event_list) / sizeof(papi_event_list[0]);
			if (_papi_num_avail_counters < _papi_num_counters) {
#ifdef DEBUG
				std::cerr << "PAPI WARNING: Not enough hw counter available. Skipping counters " << _papi_num_avail_counters << " - " << _papi_num_counters << std::endl;
#endif /* DEBUG */
				_papi_num_counters = _papi_num_avail_counters;
			}
			for (int i = 0; i < _papi_num_counters; i++) {
#ifdef DEBUG
				std::cerr << "PAPI INFO: adding HW Counter [" << i << "]  " << papi_event_list[i].event_name << " (" << papi_event_list[i].event_code << ")" << std::endl;
#endif /* DEBUG */
				if (PAPI_add_event(_papi_EventSet, papi_event_list[i].event_code) != PAPI_OK) {
#ifdef DEBUG
					std::cerr << "Could not add event " << papi_event_list[i].event_name << std::endl;
#endif /* DEBUG */
				}
			}
			PAPI_start(_papi_EventSet);
		}
#endif
	}

	~Timer() {
#ifdef PAPI
		delete[] _papi_start;
		delete[] _papi_stop;
#endif
	}

	void start() {
		_stop = 0.;
		if (_state == TIMER_HALTED) {
			_start = timer();
			_state = TIMER_RUNNING;
#ifdef PAPI
			PAPI_CHECK( (PAPI_read( _papi_EventSet, _papi_start)), "Failed reading counters.");
#endif
		}
		else
			std::cerr << "WARNING: Timer already running" << std::endl;
	}

	void stop() {
		if (_state == TIMER_RUNNING) {
			_stop = timer();
			_state = TIMER_HALTED;
			_etime += _stop - _start;
#ifdef PAPI
			PAPI_CHECK((PAPI_read( _papi_EventSet, _papi_stop)), "Failed reading counters.");
#endif
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

#ifdef ENABLE_MPI
	void set_sync(bool sync) {
		_synced = sync;
	}
#endif

#ifdef PAPI
	/* get number of used papi_counters */
	int get_papi_num_counters() {
		return _papi_num_counters;
	}
	/* get counter value between stop and last start */
	long long get_papi_counter(int index) {
		return (_papi_stop[index] - _papi_start[index]);
	}
	/* get counter name */
	char* get_papi_counter_name(int index) {
		return (papi_event_list[index].event_name);
	}
#endif /* PAPI */

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

/*
 * Christoph Niethammer
 */
#ifndef TIMER_H
#define TIMER_H

/* We use MPIs Wtime in parallel application, else clock */
#ifdef PARALLEL 
	#include <mpi.h>
#else
	#include <sys/time.h>
#endif

typedef enum {
	TIMER_HALTED = 0,
	TIMER_RUNNING = 1
} timer_state;

//! @brief This class is used to messure times in sequential and parallel versions
//! @author Christoph Niethammer
class Timer {
	double _start;		// stop time
	double _stop;		// start time
	double _etime;		// elapsed time
	timer_state _state;	// timer state
#ifdef PARALLEL
	bool _synced;		// timer syncs processes
#endif

public:
	Timer() {
		reset();
		#ifdef PARALLEL
		set_sync(true);
		#endif
	}

	void start() { 
		_stop = 0.; 
		if (_state == TIMER_HALTED) {
			_start = timer(); 
			_state = TIMER_RUNNING;
		} else std::cerr << "WARNING: Timer already running" << std::endl;
	}

	void stop()  { 
		if (_state == TIMER_RUNNING) {
			_stop  = timer(); 
			_state = TIMER_HALTED;
			_etime += _stop - _start;
		} else std::cerr << "WARNING: Timer not running" << std::endl;
	}

	void reset() { 
		_state = TIMER_HALTED;
		_start = _stop = _etime =  0.;
	}

	double get_start() { return _start; }
	double get_end()   { return _stop; }
	double get_etime() { return _etime; }
	timer_state get_state() { return _state; }

#ifdef PARALLEL
	void set_sync (bool sync) { _synced = sync; }
#endif

private:
	double timer() { 
		double time;
		#ifdef PARALLEL
			if (_synced)
				MPI_Barrier(MPI_COMM_WORLD);
			time = MPI_Wtime();
		#else
			struct timeval tmp_time;
			gettimeofday(&tmp_time, NULL);
			time = (1.0e6 * (double)tmp_time.tv_sec + (double)tmp_time.tv_usec) / 1.0e6;
		#endif
		return time;
	}
};
#endif

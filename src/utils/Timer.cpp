#include "utils/Timer.h"
#include "mardyn_assert.h"
	void Timer::start() {
		if (!isActive()) {
			mardyn_assert(_state == TIMER_HALTED);
			return;
		}

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

	void Timer::stop() {
		if (!isActive()) {
			mardyn_assert(_state == TIMER_HALTED);
			return;
		}

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

	void Timer::reset() {
		if (!isActive()) {
			mardyn_assert(_state == TIMER_HALTED);
			return ;
		}

		_state = TIMER_HALTED;
		_start = _stop = _etime = 0.;
#if WITH_PAPI
		if (_collect_papi) {
			PAPI_reset(_papi_EventSet);
            for (int i = 0; i < _papi_num_counters; i++) {
				_papi_counter[i] = 0;
			}
		}
#endif /* WITH_PAPI */
	}

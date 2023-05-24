#include "io/LoadBalanceWriter.h"

#include "Simulation.h"
#include "TimerProfiler.h"
#include "parallel/DomainDecompBase.h"

LoadbalanceWriter::LoadbalanceWriter::LoadbalanceWriter() :
	_writeFrequency(1),
	_outputFilename("mardyn-LB.dat"),
	_defaultTimer(nullptr),
	_timerNames(0),
	_times(0),
	_global_times(0),
	_simsteps(0),
	_warninglevels()
{}

void LoadbalanceWriter::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << std::endl;
	xmlconfig.getNodeValue("averageLength", _averageLength);
	global_log->info() << "Average length: " << _averageLength << std::endl;
	xmlconfig.getNodeValue("outputfilename", _outputFilename);
	global_log->info() << "Output filename: " << _outputFilename << std::endl;

	XMLfile::Query query = xmlconfig.query("timers/timer");
	std::string oldpath = xmlconfig.getcurrentnodepath();
	for(auto timerIter = query.begin(); timerIter; ++timerIter ) {
		xmlconfig.changecurrentnode(timerIter);
		std::string timername;
		xmlconfig.getNodeValue("name", timername);
		if(timername != LB_WRITER_DEFAULT_TIMER_NAME) {
			_timerNames.push_back(timername);
		}
		double warninglevel = 0;
		xmlconfig.getNodeValue("warninglevel", warninglevel);
		_warninglevels[timername] = warninglevel;

		bool incrementalTimer = false;  // if the timer is increasing in every time step, this should be true
		xmlconfig.getNodeValue("incremental", incrementalTimer);
		_incremental[timername] = incrementalTimer;
		if (incrementalTimer) {
			_incremental_previous_times[timername] = 0.;
		}

		global_log->info() << "Added timer for LB monitoring: " << timername << ", warninglevel: " << warninglevel
				<< ", incremental: " << incrementalTimer << std::endl;

	}
	xmlconfig.changecurrentnode(oldpath);
}

void LoadbalanceWriter::init(ParticleContainer */*particleContainer*/,
                             DomainDecompBase *domainDecomp, Domain */*domain*/) {
	std::string default_timer_name(LB_WRITER_DEFAULT_TIMER_NAME);

	// memory of this timer is managed by TimerProfiler.
	_defaultTimer = new Timer();

	global_simulation->timers()->registerTimer(default_timer_name, std::vector<std::string>{"SIMULATION"}, _defaultTimer);
	_timerNames.push_back(default_timer_name);

	size_t timestep_entry_offset = 2 * _timerNames.size();
	_simsteps.reserve(_writeFrequency);
	_times.reserve(timestep_entry_offset * _writeFrequency);
	_sum_times.reserve(_timerNames.size() * _writeFrequency);

	// global_times is needed at every process!
	_global_times.reserve(timestep_entry_offset * _writeFrequency);

	if (domainDecomp->getRank() == 0) {
		_global_sum_times.reserve(_timerNames.size() * _writeFrequency);
		_global_max_average_times.reserve(_timerNames.size() * (_writeFrequency+_averageLength));
		_global_sum_average_times.reserve(_timerNames.size() * (_writeFrequency+_averageLength));
	}

	_defaultTimer->reset();
	_defaultTimer->start();

	if(domainDecomp->getRank() == 0) {
		writeOutputFileHeader();
	}
}

void LoadbalanceWriter::endStep(
        ParticleContainer */*particleContainer*/,
        DomainDecompBase *domainDecomp, Domain */*domain*/,
        unsigned long simstep
)  {
	_defaultTimer->stop();
	recordTimes(simstep);
	_defaultTimer->reset();

	if((simstep % _writeFrequency) == 0) {
		  LoadbalanceWriter::flush(domainDecomp);
	}
	_defaultTimer->start();
}

void LoadbalanceWriter::writeOutputFileHeader() {
	std::ofstream outputfile(_outputFilename, std::ofstream::out);
	outputfile << "#simstep";
	for(const auto& timername : _timerNames) {
		outputfile << "\t#" << timername <<"#\t\t";
	}
	outputfile << "\nstep";
	for(const auto& timername : _timerNames) {
		outputfile << "\tmin\tmax\tf_LB\timbalance\timbalance_average";
	}
	outputfile << std::endl;
	outputfile.close();
}

void LoadbalanceWriter::recordTimes(unsigned long simstep) {
	_simsteps.push_back(simstep);
	for(const auto& timername : _timerNames) {
		Timer *timer = global_simulation->timers()->getTimer(timername);
		double time = timer->get_etime();
		if(_incremental[timername]){
			// timer is incremental, so we should subtract the previous time of the timer from the current time.
			time -= _incremental_previous_times[timername];
			_incremental_previous_times[timername] += time;  // add the time from the current
		}
		_times.push_back(time);  // needed for maximum
		_times.push_back(-time); // needed for minimum
		_timesForAverageBuffer.emplace_back(time);  // needed for average
	}
}

void LoadbalanceWriter::flush(DomainDecompBase* domainDecomp) {
#ifdef ENABLE_MPI
	//! @todo If this shall become a general LB monitor/manager a MPI_Allreduce will be needed here
	MPI_CHECK(
			MPI_Reduce( _times.data(), _global_times.data(), _times.size(), MPI_DOUBLE, MPI_MAX,0, domainDecomp->getCommunicator()));

	// only sum over the positive times, so "+2"
	for (size_t ind = 0; ind < _times.size(); ind += 2) {
		// sum of +time (not the negative ones!)
		_sum_times.push_back(_times[ind]);
	}
	MPI_CHECK(MPI_Reduce(_sum_times.data(), _global_sum_times.data(), _sum_times.size(), MPI_DOUBLE, MPI_SUM, 0,
						 domainDecomp->getCommunicator()));

	auto averagedTimes = getAveragedTimes();

	_global_max_average_times.resize(averagedTimes.size(), 0);
	MPI_CHECK(MPI_Reduce(averagedTimes.data(), _global_max_average_times.data(), averagedTimes.size(), MPI_DOUBLE, MPI_MAX, 0,
						 domainDecomp->getCommunicator()));
	_global_sum_average_times.resize(averagedTimes.size(), 0);
	MPI_CHECK(MPI_Reduce(averagedTimes.data(), _global_sum_average_times.data(), averagedTimes.size(), MPI_DOUBLE, MPI_SUM, 0,
						 domainDecomp->getCommunicator()));
#else
	//! @todo in case we do not reuse _times later on, we may use here  _global_times = std::move(_times);
	_global_times = _times;
	for(size_t ind = 0; ind < _times.size(); ind += 2) {
		_sum_times.push_back(_times[ind]);
		_global_sum_times.push_back(_times[ind]);
	}

	auto averagedTimes = getAveragedTimes();

	// sum and max are actually the same if only one process is used.
	_global_max_average_times = averagedTimes;
	_global_sum_average_times = averagedTimes;
#endif
	if(0 == domainDecomp->getRank()) {
		std::ofstream outputfile(_outputFilename,  std::ofstream::app);
		for(size_t i = 0; i < _simsteps.size(); ++i) {
			writeLBEntry(i, outputfile, domainDecomp->getNumProcs());
		}
		outputfile.close();
	}
	resetTimes();
}

std::vector<double> LoadbalanceWriter::getAveragedTimes() {
	std::vector<double> averagedTimes(_timesForAverageBuffer.size(), 0.);
	const auto numTimers = _timerNames.size();
	mardyn_assert(_timesForAverageBuffer.size() % numTimers == 0);
	std::vector<double> currentAverageSum(numTimers, 0.);
	for (size_t ind = 0, averageIndex = 0; ind < _timesForAverageBuffer.size(); ind += numTimers, ++averageIndex) {
		for (size_t timer_id = 0; timer_id < numTimers; ++timer_id) {
			currentAverageSum[timer_id] += _timesForAverageBuffer[ind + timer_id];
			if (averageIndex >= _averageLength) {
				currentAverageSum[timer_id] -= _timesForAverageBuffer[ind + timer_id - _averageLength * numTimers];
				averagedTimes[ind + timer_id] = currentAverageSum[timer_id] / static_cast<double>(_averageLength);
			} else {
				averagedTimes[ind + timer_id] = currentAverageSum[timer_id] / static_cast<double>(averageIndex + 1);
			}
		}
	}
	return averagedTimes;
}

void LoadbalanceWriter::writeLBEntry(size_t id, std::ofstream &outputfile, int numRanks) {
	outputfile << _simsteps[id];
	size_t timestep_entry_offset = 2* _timerNames.size();
	size_t base_offset = timestep_entry_offset * id;
	size_t single_offset = _timerNames.size() * id;
	size_t single_offset_average = _timerNames.size() * id + _lastTimesOldValueCount;
	for(size_t timer_id = 0; timer_id < _timerNames.size(); ++timer_id) {
		size_t timer_offset = 2 * timer_id;
		size_t offset = base_offset + timer_offset;
		double max = _global_times[offset];
		double min = -_global_times[offset + 1];
		// imbalance = 1 - (average of time) / maximal time -- fraction: time spent useless / really needed time
		double optimaltime = (_global_sum_times[timer_id + single_offset] / numRanks);
		double optimaltime_average = (_global_sum_average_times[timer_id + single_offset_average] / numRanks);
		double imbalance = 1 - optimaltime / max;
		double imbalance_average = 1 - optimaltime_average / _global_max_average_times[timer_id + single_offset_average];
		double f_LB = max / min;
		std::string timername = _timerNames[timer_id];
		if(_warninglevels.count(timername) > 0 && _warninglevels[timername] < f_LB) {
			displayWarning(_simsteps[id], timername, f_LB);
		}
		outputfile << "\t" << min << "\t" << max << "\t" << f_LB << "\t" << imbalance << "\t" << imbalance_average;
	}
	outputfile << std::endl;
}

void LoadbalanceWriter::displayWarning(unsigned long simstep, const std::string& timername, double f_LB) {
	global_log->warning() << "Load balance limit exceeded in simstep " << simstep
		<< " for timer " << timername
		<< ", limit: " << _warninglevels[timername]
		<< "  value: " << f_LB << std::endl;
}

void LoadbalanceWriter::resetTimes() {
	_times.clear();
	_global_times.clear();
	_simsteps.clear();
	_sum_times.clear();
	_global_sum_times.clear();
	_global_max_average_times.clear();
	_global_sum_average_times.clear();

	// Removes first few elements, s.t., only the last _averageLength - 1 elements remain.
	auto numTimers = _timerNames.size();
	if (_timesForAverageBuffer.size() > (_averageLength - 1) * numTimers) {
		_timesForAverageBuffer.erase(_timesForAverageBuffer.begin(),
									 _timesForAverageBuffer.end() - numTimers * (_averageLength - 1));
		mardyn_assert(_timesForAverageBuffer.size() == (_averageLength - 1) * numTimers);
	}
	_lastTimesOldValueCount = _timesForAverageBuffer.size();
}

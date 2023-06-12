#ifdef TASKTIMINGPROFILE

#include <bhfmm/FastMultipoleMethod.h>
#include <Simulation.h>
#include <climits>
#include "TaskTimingProfiler.h"

TaskTimingProfiler::TaskTimingProfiler() {

}

void TaskTimingProfiler::init(unsigned long estimatedNumberOfTasks) {
	if (mardyn_get_max_threads() == 1){
		timings.push_back(new std::vector<TaskTimingProfiler::timing>);
		timings[0]->reserve(estimatedNumberOfTasks * global_simulation->getNumberOfTimesteps());
		return;
	}

	for (int i = 0; i < mardyn_get_max_threads(); ++i) {
		timings.push_back(new std::vector<TaskTimingProfiler::timing>);
		// assume no thread will bear more than 70% of the total workload
		// it would be preferable if the vectors are not rescaled since this leads to measurement overhead
		// TODO: better estimation to save memory...
		timings[i]->reserve(std::ceil(estimatedNumberOfTasks * .7) * global_simulation->getNumberOfTimesteps());
	}
}

TaskTimingProfiler::~TaskTimingProfiler() {
	for (unsigned int i = 0; i < timings.size(); ++i) {
		delete timings[i];
	}
}

unsigned long TaskTimingProfiler::start() {
	return _rdtsc();
}

void TaskTimingProfiler::stop(unsigned long start, int taskId) {
	unsigned long stop = _rdtsc();
	timing t {
			start,
			stop,
			taskId
	};
	timings[mardyn_get_thread_num()]->push_back(t);
}

void TaskTimingProfiler::dump(std::string filename) {
    Log::global_log->info() << "Writing task timings to " << filename << " ... " << std::flush;

    // find first tick to subtract from all others
    unsigned long startTick = ULLONG_MAX;
    for (unsigned int i = 0; i < timings.size(); ++i) {
        startTick = std::min(startTick, timings[i]->at(0).start);
    }

	std::ofstream outputFile(filename);
	outputFile << std::setw(14) << "Start"
			   << std::setw(18) << "Stop"
			   << std::setw(7) << "Type"
			   << std::setw(7) << "Thread"
			   << std::endl;

	for (unsigned int threadId = 0; threadId < timings.size(); ++threadId) {
		for (timing t : *timings[threadId]) {
			outputFile << std::setw(14) << t.start - startTick
					   << std::setw(18) << t.stop - startTick
					   << std::setw(7) << t.taskId
					   << std::setw(7) << threadId
					   << std::endl;
		}
	}
	outputFile.close();

    std::cout << "done!" << std::endl;

#ifdef QUICKSCHED
	Log::global_log->info() << "Task ids:" << std::endl
					   << std::setw(42) << "" << "P2PPreprocessSingleCell  = "
					   << bhfmm::FastMultipoleMethod::taskType::P2PPreprocessSingleCell << std::endl
					   << std::setw(42) << "" << "P2Pc08StepBlock          = "
					   << bhfmm::FastMultipoleMethod::taskType::P2Pc08StepBlock << std::endl
					   << std::setw(42) << "" << "P2PPostprocessSingleCell = "
					   << bhfmm::FastMultipoleMethod::taskType::P2PPostprocessSingleCell << std::endl
					   << std::setw(42) << "" << "P2MCompleteCell          = "
					   << bhfmm::FastMultipoleMethod::taskType::P2MCompleteCell << std::endl
					   << std::setw(42) << "" << "M2MCompleteCell          = "
					   << bhfmm::FastMultipoleMethod::taskType::M2MCompleteCell << std::endl
					   << std::setw(42) << "" << "M2LInitializeCell        = "
					   << bhfmm::FastMultipoleMethod::taskType::M2LInitializeCell << std::endl
					   << std::setw(42) << "" << "M2LInitializeSource      = "
					   << bhfmm::FastMultipoleMethod::taskType::M2LInitializeSource << std::endl
					   << std::setw(42) << "" << "M2LTranslation           = "
					   << bhfmm::FastMultipoleMethod::taskType::M2LTranslation << std::endl
					   << std::setw(42) << "" << "M2LPair2Way              = "
					   << bhfmm::FastMultipoleMethod::taskType::M2LPair2Way << std::endl
					   << std::setw(42) << "" << "M2LFinalizeCell          = "
					   << bhfmm::FastMultipoleMethod::taskType::M2LFinalizeCell << std::endl
					   << std::setw(42) << "" << "L2LCompleteCell          = "
					   << bhfmm::FastMultipoleMethod::taskType::L2LCompleteCell << std::endl
					   << std::setw(42) << "" << "L2PCompleteCell          = "
					   << bhfmm::FastMultipoleMethod::taskType::L2PCompleteCell << std::endl;
#endif /* QUICKSCHED */
}

#endif /* TASKTIMINGPROFILE */

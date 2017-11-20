/*
 * TaskTimingProfiler.h
 *
 *  Created on: August 13, 2017
 *      Author: gratlf
 */

#ifndef TASKTIMINGPROFILER_H
#define TASKTIMINGPROFILER_H

#ifdef TASKTIMINGPROFILE

#include <vector>
#include <string>
#include <cmath>
#include <x86intrin.h>

#include "WrapOpenMP.h"

class TaskTimingProfiler {
public:
	TaskTimingProfiler();
	~TaskTimingProfiler();

	void init(unsigned long estimatedNumberOfTasks);
	unsigned long start();
	void stop(unsigned long start, int taskId);
	void dump(std::string filename);

private:
	struct timing {
		unsigned long start,
			          stop;
		int taskId;
	};
	// one timing vector per thread
	std::vector<std::vector<TaskTimingProfiler::timing>*> timings;
};

#endif /* TASKTIMINGPROFILE */
#endif //TASKTIMINGPROFILER_H

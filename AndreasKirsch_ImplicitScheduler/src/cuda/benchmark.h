/*
 * benchmark.h
 *
 *  Created on: Aug 15, 2011
 *      Author: andreas
 */

#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include <vector>
#include <numeric>
#include <string>

#include "config.h"
#include "particleContainer/Cell.h"

class Measure {
protected:
	std::vector< double > values;
public:
	void addDataPoint(const double value) {
		values.push_back( value );
	}

	double getAverage() const {
		return std::accumulate( values.begin(), values.end(), 0.0 ) / values.size();
	}

	bool isUsed() const {
		return !values.empty();
	}

	int getCount() const {
		return values.size();
	}

	double operator [] (int index) const {
		return values[index];
	}

	operator double() const {
		return getAverage();
	}
};

struct SimulationStats {
	// run measurements
	Measure totalTime;
	Measure CUDA_frameTime, CUDA_preTime, CUDA_postTime, CUDA_singleTime, CUDA_pairTime, CUDA_processingTime;

	// frame measurements
	Measure potentials, virials;
	Measure forceRMSError, torqueRMSError;

	// run/build info
	int timeSteps;
	int moleculeCount;
	double cutOffRadius;

	std::string name;

	void writeFrameStats( const std::string &frameFile );
	void writeRunStats( const std::string &buildFile );

	static void writeCellStats( const std::string &domainInfoFile, const std::vector<Cell> &cells, const std::string &domainInfo );
};

extern SimulationStats simulationStats;

#endif /* BENCHMARK_H_ */

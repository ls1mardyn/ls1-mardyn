/*
 * LoadCalc.cpp
 *
 *  Created on: 10.06.2017
 *      Author: griebel_s
 */

#include <armadillo>
#include "mpi.h"

#include "LoadCalc.h"
#include "DomainDecompBase.h"



std::vector<double> TunerLoad::readVec(std::istream& in, int& count1, int& count2) {
	std::vector<double> vec;
	int tempCount1 = 0;
	int tempCount2 = -1;
	while (std::isdigit(in.peek())) {
		++tempCount1;
		std::string line;
		std::getline(in, line);
		auto pos = line.begin();
		int i = 0;

		while (true) {
			++i;
			auto new_pos = std::find(pos, line.end(), ';');
			auto subStr = std::string(pos, new_pos);
			vec.push_back(std::stod(subStr));
			if (new_pos == line.end()) {
				break;
			}
			pos = new_pos + 1;
		}

		if (tempCount2 == -1) {
			tempCount2 = i;
		} else {
			if (i != tempCount2) {
				Log::global_log->error_always_output()
						<< "The file contains data of 2D-vectors with a different amounts of elements in the second dimension!"
						<< std::endl;
				Log::global_log->error_always_output()
						<< "This means the files is corrupted. Please remove it (or disallow the tuner to read from inputfiles) before restarting!"
						<< std::endl;
				Simulation::exit(1);
			}
		}
	}
	count1 = tempCount1;
	count2 = tempCount2;
	return vec;
}

std::array<double, 3> TunerLoad::calcConsts(const std::vector<double>& timeVec, bool inner) {
	//take the median of the last 10% of the times
	std::array<double,3> consts {};
	//the blocks are here to avoid at least some copy paste mistakes
	{
		std::vector<double> temp1 {};
		for (int i = timeVec.size() / _count2 * 9 / 10;i < _count1;++i) {
			double quot = inner ? i * (i - 1)/2 : i * i;
			//if quot is zero (which can happen if i == 0 || (i == 1 && own)), than Infinity would be pushed into the vector
			//so this is caught before it happens
			if (quot == 0) {
				quot = 1;
			}
			temp1.push_back(accessVec(timeVec, i, 0) / quot);
		}
		//calculates the median
		std::sort(temp1.begin(), temp1.end());
		consts[0] = temp1.at(temp1.size() / 2);
	}
	{
		std::vector<double> temp2 {};
		for (int i = timeVec.size() / _count1 * 9 / 10;i < _count2;++i) {
			double quot = inner ? i * (i - 1)/2 : i * i;
			if (quot == 0) {
				quot = 1;
			}
			temp2.push_back(accessVec(timeVec, 0, i) / quot);
		}
		std::sort(temp2.begin(), temp2.end());
		consts[1] = temp2.at(temp2.size() / 2);
	}
	{
		std::vector<double> temp3 {};
		for (int index1 = timeVec.size() / _count2 * 9 / 10;index1 < _count1;++index1) {
			for (int index2 = timeVec.size() / _count1 * 9 / 10;index2 < _count2;++index2) {
				double fac1;
				double fac2;
				double quot;
				if (inner) {
					fac1 = index1 * (index1 - 1);
					fac2 = index2 * (index2 - 1);
					quot = index1 * index2;
				} else {
					fac1 = index1 * index1;
					fac2 = index2 * index2;
					quot = index1 * index2 * 2; // my own type 1 with neighbor type 2, my own type 2 with neighbor type 1
				}
				if (quot == 0) {
					quot = 1;
				}
				temp3.push_back((accessVec(timeVec, index1, index2) - fac1 * consts[0] - fac2 * consts[1]) / quot);
			}
		}
		std::sort(temp3.begin(), temp3.end());
		consts[2] = temp3.at(temp3.size() / 2);
	}
	return consts;
}

TunerLoad::TunerLoad(int count1, int count2,
		std::vector<double>&& ownTime, std::vector<double>&& faceTime,
		std::vector<double>&& edgeTime, std::vector<double>&& cornerTime) :
		_count1 { count1 }, _count2 { count2 },
		_ownTime { std::move(ownTime) }, _faceTime {std::move(faceTime) },
		_edgeTime { std::move(edgeTime) }, _cornerTime { std::move(cornerTime) },
		_ownConst ( calcConsts(_ownTime, true) ), _faceConst ( calcConsts(_faceTime, false) ),
		_edgeConst ( calcConsts(_edgeTime, false) ), _cornerConst ( calcConsts(_cornerTime, false) ) {

	if (_ownTime.size() != size_t(_count1 * _count2)) {
		global_log->error_always_output()
				<< "_edgeTime was initialized with the wrong size of "
				<< _ownTime.size() << " expected: " << _count1 * _count2;
	}

	if (_faceTime.size() != size_t(count1 * _count2)) {
		global_log->error_always_output()
				<< "_edgeTime was initialized with the wrong size of "
				<< _faceTime.size() << " expected: " << _count1 * _count2;
	}

	if (_edgeTime.size() != size_t(_count1 * _count2)) {
		global_log->error_always_output()
				<< "_edgeTime was initialized with the wrong size of "
				<< _edgeTime.size() << " expected: " << _count1 * _count2;
	}

	if (_cornerTime.size() != size_t(_count1 * _count2)) {
		global_log->error_always_output()
				<< "_edgeTime was initialized with the wrong size of "
				<< _cornerTime.size() << " expected: " << _count1 * _count2;
	}
}

TunerLoad TunerLoad::read(std::istream& stream){
	std::string inStr {};
	std::getline(stream, inStr);
	if(inStr != "Vectorization Tuner File"){
		Log::global_log->error() << "The tunerfile is corrupted! Missing header \"Vectorization Tuner File\"";
		Log::global_log->error() << "Please remove it or fix it before restarting!";
		Simulation::exit(1);
	}

	int count1;
	int count2;

	std::getline(stream, inStr);
	if(inStr != "own"){
		Log::global_log->error() << "The tunerfile is corrupted! Missing Section \"own\"";
		Log::global_log->error() << "Please remove it or fix it before restarting!";
		Simulation::exit(1);
	}
	auto ownTime = readVec(stream, count1, count2);
	std::getline(stream, inStr);

	if(inStr != "face"){
		Log::global_log->error() << "The tunerfile is corrupted! Missing Section \"face\"";
		Log::global_log->error() << "Please remove it or fix it before restarting!";
		Simulation::exit(1);
	}
	auto faceTime = readVec(stream, count1, count2);
	std::getline(stream, inStr);

	if(inStr != "edge"){
		Log::global_log->error() << "The tunerfile is corrupted! Missing Section \"edge\"";
		Log::global_log->error() << "Please remove it or fix it before restarting!";
		Simulation::exit(1);
	}
	auto edgeTime = readVec(stream, count1, count2);
	std::getline(stream, inStr);

	if(inStr != "corner"){
		Log::global_log->error() << "The tunerfile is corrupted! Missing Section \"corner\"";
		Log::global_log->error() << "Please remove it or fix it before restarting!";
		Simulation::exit(1);
	}
	auto cornerTime = readVec(stream, count1, count2);
	return TunerLoad {count1, count2, std::move(ownTime), std::move(faceTime), std::move(edgeTime), std::move(cornerTime)};
}


// MEASURELOAD
std::string MeasureLoad::TIMER_NAME = "SIMULATION_FORCE_CALCULATION";

void MeasureLoad::prepareLoads(DomainDecompBase* decomp, MPI_Comm& comm) {
	int numRanks = decomp->getNumProcs();

	// owntime = time since last check
	double ownTime = global_simulation->timers()->getTime(TIMER_NAME) - _previousTime;
	_previousTime += ownTime;

	std::vector<double> neededTimes(numRanks);
	MPI_Gather(&ownTime, 1, MPI_DOUBLE, neededTimes.data(), numRanks, MPI_DOUBLE, 0, comm);


	std::vector<unsigned long> statistics = global_simulation->getMoleculeContainer()->getParticleCellStatistics();
	int maxParticleCount = statistics.size();
	if(numRanks < maxParticleCount){
		Log::global_log->fatal() << "Not enough processes to sample from. Aborting!" << std::endl;
		global_simulation->exit(559);
	}

	int global_maxParticleCount = 0;
	MPI_Allreduce(&maxParticleCount, &global_maxParticleCount, 1, MPI_INT, MPI_MAX, comm);
	statistics.resize(global_maxParticleCount);
	if (decomp->getRank() == 0) {
		std::vector<unsigned long> global_statistics(global_maxParticleCount * numRanks);
		MPI_Gather(statistics.data(), statistics.size(), MPI_UINT64_T, global_statistics.data(),
				global_statistics.size(), MPI_UINT64_T, 0, comm);

		// right hand side = global_statistics ^ T \cdot neededTimes
		std::vector<double> right_hand_side(global_maxParticleCount, 0.);
		for (int particleCount = 0; particleCount < global_maxParticleCount; particleCount++) {
			for (int rank = 0; rank < numRanks; rank++) {
				right_hand_side[particleCount] += global_statistics[global_maxParticleCount * rank + particleCount]
						* neededTimes[rank];
			}
		}

		// system matrix is: global_statistics ^ T \cdot global_statistics
		std::vector<unsigned long> system_matrix(global_maxParticleCount * global_maxParticleCount, 0ul);
		for (int particleCount1 = 0; particleCount1 < global_maxParticleCount; particleCount1++) {
			for (int rank = 0; rank < numRanks; rank++) {
				for (int particleCount2 = 0; particleCount2 < global_maxParticleCount; particleCount2++) {
					system_matrix[particleCount1 * global_maxParticleCount + particleCount2] +=
							global_statistics[global_maxParticleCount * rank + particleCount1]
									* global_statistics[global_maxParticleCount * rank + particleCount2];
				}
			}
		}

		// now we have to solve: system_matrix \cdot cell_time_vector = right_hand_side
		arma::mat arma_system_matrix(global_maxParticleCount, global_maxParticleCount);
		for(int row = 0; row < global_maxParticleCount; row ++){
			for(int column = 0; column < global_maxParticleCount; column ++){
				arma_system_matrix[row * global_maxParticleCount + column] = system_matrix[row * global_maxParticleCount
						+ column];
			}
		}
		arma::vec arma_rhs(right_hand_side);
		arma::vec cell_time_vec = arma::solve(arma_system_matrix, arma_rhs);

	} else {
		MPI_Gather(statistics.data(), statistics.size(), MPI_UINT64_T, nullptr,
				0 /*here insignificant*/, MPI_UINT64_T, 0, comm);
	}

}

/*
 * LoadCalc.cpp
 *
 *  Created on: 10.06.2017
 *      Author: griebel_s
 */
#ifdef MARDYN_ARMADILLO
#include <armadillo>
#endif

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
	// take the median of the last 10% of the times
	std::array<double, 3> consts { };
	// the blocks are here to avoid at least some copy paste mistakes
	{
		std::vector<double> temp1 { };
		for (int i = timeVec.size() / _count2 * 9 / 10; i < _count1; ++i) {
			double quot = inner ? i * (i - 1) / 2 : i * i;
			//if quot is zero (which can happen if i == 0 || (i == 1 && own)), than Infinity would be pushed into the vector
			//so this is caught before it happens
			if (quot == 0) {
				quot = 1;
			}
			temp1.push_back(accessVec(timeVec, i, 0) / quot);
		}
		// calculates the median
		std::sort(temp1.begin(), temp1.end());
		consts[0] = temp1.at(temp1.size() / 2);
	}
	{
		std::vector<double> temp2 { };
		for (int i = timeVec.size() / _count1 * 9 / 10; i < _count2; ++i) {
			double quot = inner ? i * (i - 1) / 2 : i * i;
			if (quot == 0) {
				quot = 1;
			}
			temp2.push_back(accessVec(timeVec, 0, i) / quot);
		}
		std::sort(temp2.begin(), temp2.end());
		consts[1] = temp2.at(temp2.size() / 2);
	}
	{
		std::vector<double> temp3 { };
		for (int index1 = timeVec.size() / _count2 * 9 / 10; index1 < _count1; ++index1) {
			for (int index2 = timeVec.size() / _count1 * 9 / 10; index2 < _count2; ++index2) {
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

TunerLoad::TunerLoad(int count1, int count2, std::vector<double>&& ownTime, std::vector<double>&& faceTime,
		std::vector<double>&& edgeTime, std::vector<double>&& cornerTime) :
		_count1 { count1 }, _count2 { count2 }, _ownTime { std::move(ownTime) }, _faceTime { std::move(faceTime) }, _edgeTime {
				std::move(edgeTime) }, _cornerTime { std::move(cornerTime) }, _ownConst(calcConsts(_ownTime, true)), _faceConst(
				calcConsts(_faceTime, false)), _edgeConst(calcConsts(_edgeTime, false)), _cornerConst(
				calcConsts(_cornerTime, false)) {

	if (_ownTime.size() != size_t(_count1 * _count2)) {
		global_log->error_always_output() << "_edgeTime was initialized with the wrong size of " << _ownTime.size()
				<< " expected: " << _count1 * _count2;
	}

	if (_faceTime.size() != size_t(count1 * _count2)) {
		global_log->error_always_output() << "_edgeTime was initialized with the wrong size of " << _faceTime.size()
				<< " expected: " << _count1 * _count2;
	}

	if (_edgeTime.size() != size_t(_count1 * _count2)) {
		global_log->error_always_output() << "_edgeTime was initialized with the wrong size of " << _edgeTime.size()
				<< " expected: " << _count1 * _count2;
	}

	if (_cornerTime.size() != size_t(_count1 * _count2)) {
		global_log->error_always_output() << "_edgeTime was initialized with the wrong size of " << _cornerTime.size()
				<< " expected: " << _count1 * _count2;
	}
}

TunerLoad TunerLoad::read(std::istream& stream) {
	std::string inStr { };
	std::getline(stream, inStr);
	if (inStr != "Vectorization Tuner File") {
		Log::global_log->error() << "The tunerfile is corrupted! Missing header \"Vectorization Tuner File\"";
		Log::global_log->error() << "Please remove it or fix it before restarting!";
		Simulation::exit(1);
	}

	int count1;
	int count2;

	std::getline(stream, inStr);
	if (inStr != "own") {
		Log::global_log->error() << "The tunerfile is corrupted! Missing Section \"own\"";
		Log::global_log->error() << "Please remove it or fix it before restarting!";
		Simulation::exit(1);
	}
	auto ownTime = readVec(stream, count1, count2);
	std::getline(stream, inStr);

	if (inStr != "face") {
		Log::global_log->error() << "The tunerfile is corrupted! Missing Section \"face\"";
		Log::global_log->error() << "Please remove it or fix it before restarting!";
		Simulation::exit(1);
	}
	auto faceTime = readVec(stream, count1, count2);
	std::getline(stream, inStr);

	if (inStr != "edge") {
		Log::global_log->error() << "The tunerfile is corrupted! Missing Section \"edge\"";
		Log::global_log->error() << "Please remove it or fix it before restarting!";
		Simulation::exit(1);
	}
	auto edgeTime = readVec(stream, count1, count2);
	std::getline(stream, inStr);

	if (inStr != "corner") {
		Log::global_log->error() << "The tunerfile is corrupted! Missing Section \"corner\"";
		Log::global_log->error() << "Please remove it or fix it before restarting!";
		Simulation::exit(1);
	}
	auto cornerTime = readVec(stream, count1, count2);
	return TunerLoad { count1, count2, std::move(ownTime), std::move(faceTime), std::move(edgeTime), std::move(
			cornerTime) };
}

// MEASURELOAD
std::string MeasureLoad::TIMER_NAME = "SIMULATION_FORCE_CALCULATION";

MeasureLoad::MeasureLoad() :
		_times() {
	_previousTime = global_simulation->timers()->getTime(TIMER_NAME);
	_interpolationConstants = {0., 0., 0.};
	_preparedLoad = false;
}

#ifdef MARDYN_ARMADILLO
// this non-negative least-squares (nnls) algorithm is taken from:
// https://github.com/linxihui/Misc/blob/master/Practice/NMF/nnls.cpp
arma::vec nnls(const arma::mat &A, const arma::vec &b, int max_iter = 500, double tol = 1e-8) {
	/*
	 * Description: sequential Coordinate-wise algorithm for non-negative least square regression A x = b, s.t. x >= 0
	 * Reference: http://cmp.felk.cvut.cz/ftp/articles/franc/Franc-TR-2005-06.pdf
	 */

	arma::vec mu = -A.t() * b;
	arma::mat H = A.t() * A;
	arma::vec x(A.n_cols), x0(A.n_cols);
	x.fill(0);
	x0.fill(-9999);

	int i = 0;
	double tmp;
	while (i < max_iter && max(abs(x - x0)) > tol) {
		x0 = x;
		for (unsigned int k = 0; k < A.n_cols; k++) {
			tmp = x[k] - mu[k] / H.at(k, k);
			if (tmp < 0)
				tmp = 0;
			if (tmp != x[k])
				mu += (tmp - x[k]) * H.col(k);
			x[k] = tmp;
		}
		++i;
	}
	return x;
}
#endif

#ifdef ENABLE_MPI
int MeasureLoad::prepareLoads(DomainDecompBase* decomp, MPI_Comm& comm) {

#ifndef MARDYN_ARMADILLO
	global_log->info() << "not compiled with armadillo. MeasureLoad not usable." << std::endl;
	return 1;
#else
	int numRanks = decomp->getNumProcs();

	// owntime = time since last check
	double ownTime = global_simulation->timers()->getTime(TIMER_NAME) - _previousTime;
	_previousTime += ownTime;

	std::vector<double> neededTimes(numRanks);
	MPI_Gather(&ownTime, 1, MPI_DOUBLE, neededTimes.data(), 1, MPI_DOUBLE, 0, comm);

	std::vector<unsigned long> statistics = global_simulation->getMoleculeContainer()->getParticleCellStatistics();
	int maxParticlesP1 = statistics.size();

	int global_maxParticlesP1 = 0;  // maxParticle Count + 1 = degrees of freedom
	MPI_Allreduce(&maxParticlesP1, &global_maxParticlesP1, 1, MPI_INT, MPI_MAX, comm);

	if (_alwaysUseInterpolation) {
		if (numRanks < _interpolationConstants.size()) {
			Log::global_log->warning()
				<< "MeasureLoad: Not enough processes to sample from (needed _interpolationConstants.size(): "
				<< _interpolationConstants.size() << ", numRanks: " << numRanks << ")." << std::endl;
			return 1;
		}
	} else {
		// Here, we require as many ranks as particles!
		if (numRanks < global_maxParticlesP1) {
			Log::global_log->warning() << "MeasureLoad: Not enough processes to sample from (needed=maxParticlesP1: "
									   << global_maxParticlesP1 << ", numRanks: " << numRanks << ")." << std::endl;
			return 1;
		}
	}

	statistics.resize(global_maxParticlesP1);
	if (decomp->getRank() == 0) {
		std::vector<unsigned long> global_statistics(global_maxParticlesP1 * numRanks);
		MPI_Gather(statistics.data(), statistics.size(), MPI_UINT64_T, global_statistics.data(), statistics.size(),
				MPI_UINT64_T, 0, comm);

		// right hand side = global_statistics ^ T \cdot neededTimes
		std::vector<double> right_hand_side(global_maxParticlesP1, 0.);
		for (int particleCount = 0; particleCount < global_maxParticlesP1; particleCount++) {
			for (int rank = 0; rank < numRanks; rank++) {
				right_hand_side[particleCount] += global_statistics[global_maxParticlesP1 * rank + particleCount]
						* neededTimes[rank];
			}
		}

		// system matrix is: global_statistics ^ T \cdot global_statistics
		std::vector<unsigned long> system_matrix(global_maxParticlesP1 * global_maxParticlesP1, 0ul);
		for (int particleCount1 = 0; particleCount1 < global_maxParticlesP1; particleCount1++) {
			for (int rank = 0; rank < numRanks; rank++) {
				for (int particleCount2 = 0; particleCount2 < global_maxParticlesP1; particleCount2++) {
					system_matrix[particleCount1 * global_maxParticlesP1 + particleCount2] +=
							global_statistics[global_maxParticlesP1 * rank + particleCount1]
									* global_statistics[global_maxParticlesP1 * rank + particleCount2];
				}
			}
		}

		// now we have to solve: system_matrix \cdot cell_time_vector = right_hand_side
		arma::mat arma_system_matrix(global_maxParticlesP1, global_maxParticlesP1);
		for (int row = 0; row < global_maxParticlesP1; row++) {
			for (int column = 0; column < global_maxParticlesP1; column++) {
				arma_system_matrix[row * global_maxParticlesP1 + column] = system_matrix[row * global_maxParticlesP1
						+ column];
			}
		}
		if(_alwaysUseInterpolation) {
			arma::mat quadratic_equation_fit_matrix(global_maxParticlesP1, 3);
			for (int row = 0; row < global_maxParticlesP1; row++) {
				quadratic_equation_fit_matrix.at(row, 0) = row * row;
				quadratic_equation_fit_matrix.at(row, 1) = row;
				quadratic_equation_fit_matrix.at(row, 2) = 1;
			}
			arma_system_matrix = quadratic_equation_fit_matrix.t() * arma_system_matrix * quadratic_equation_fit_matrix;
			arma::vec arma_rhs(right_hand_side);
			arma_rhs = quadratic_equation_fit_matrix.t() * arma_rhs;
			arma::vec coefficient_vec = nnls(arma_system_matrix, arma_rhs);
			mardyn_assert(coefficient_vec.size() == 3);
			global_log->info() << "coefficient_vec: " << std::endl << coefficient_vec << std::endl;
			for (int i = 0; i < 3; i++) {
				_interpolationConstants[i] = coefficient_vec[i];
			}
			MPI_Bcast(_interpolationConstants.data(), 3, MPI_DOUBLE, 0, comm);
		} else {
			// old version
			arma::vec arma_rhs(right_hand_side);
			arma::vec cell_time_vec = arma::solve(arma_system_matrix, arma_rhs);

			global_log->info() << "cell_time_vec: " << cell_time_vec << std::endl;
			_times = arma::conv_to<std::vector<double> >::from(cell_time_vec);
			mardyn_assert(_times.size() == global_maxParticlesP1);
			MPI_Bcast(_times.data(), global_maxParticlesP1, MPI_DOUBLE, 0, comm);
		}
	} else {
		MPI_Gather(statistics.data(), statistics.size(), MPI_UINT64_T, nullptr, 0 /*here insignificant*/, MPI_UINT64_T,
				0, comm);
		if (_alwaysUseInterpolation) {
			MPI_Bcast(_interpolationConstants.data(), 3, MPI_DOUBLE, 0, comm);
		} else {
			// old version
			_times.resize(global_maxParticlesP1);
			MPI_Bcast(_times.data(), global_maxParticlesP1, MPI_DOUBLE, 0, comm);
		}
	}

	// interpolation constants:
	if(not _alwaysUseInterpolation) {
		calcConstants();
	}
	_preparedLoad = true;
	return 0;
#endif
}
#endif  // ENABLE_MPI
void MeasureLoad::calcConstants() {
	// we do a least squares fit of: y = a x^2 + b x + c

	// we need at least three entries for that!
	mardyn_assert(_times.size() >= 3);

	// we do a least square fit of the last 50% of the data:
	size_t start, numElements;
	{
		size_t one = _times.size() - 3;
		size_t two = _times.size() * 0.5;
		start = std::min(one, two);
		numElements = _times.size() - start;
	}

	std::array<double, 5> momentsX{};  // stores the moments of x: sum{t^i}
	std::array<double, 3> momentsYX{};  // stores the following: sum{d* t^i}

	momentsX[0] = numElements;

	for (size_t i = start; i < _times.size(); i++) {
		double x = i;
		double x2 = x * x;
		double x3 = x2 * x;
		double x4 = x2 * x2;
		double y = _times[i];
		double yx = y * x;
		double yx2 = y * x2;

		momentsX[1] += x;
		momentsX[2] += x2;
		momentsX[3] += x3;
		momentsX[4] += x4;
		momentsYX[0] += y;
		momentsYX[1] += yx;
		momentsYX[2] += yx2;
	}
#ifdef MARDYN_ARMADILLO
	// 3x3 matrix:
	arma::mat system_matrix(3, 3);

	// 3x1 vector from moments:
	arma::vec rhs(3);

	for (size_t row = 0; row < 3ul; ++row) {
		rhs[row] = momentsYX[row];
		for (size_t col = 0; col < 3ul; ++col) {
			system_matrix.at(row, col) = momentsX[row + col];
		}
	}

	arma::vec solution = arma::solve(system_matrix, rhs);

	for (size_t row = 0; row < 3ul; row++) {
		_interpolationConstants[row] = solution[2 - row];
	}
	global_log->info() << "extrapolationconst: " << solution << std::endl;
#endif

}

double MeasureLoad::getValue(int numParticles) const {
	mardyn_assert(numParticles >= 0);
	mardyn_assert(_preparedLoad);

	if(_alwaysUseInterpolation) {
		return _interpolationConstants[0] * numParticles * numParticles + _interpolationConstants[1] * numParticles +
			   _interpolationConstants[2];
	} else {
		size_t numPart = numParticles;
		if (numPart < _times.size()) {
			// if we are within the known (i.e. measured) particle count, just use the known values
			return _times[numPart];
		} else {
			// otherwise we interpolate
			return _interpolationConstants[0] * numPart * numPart + _interpolationConstants[1] * numPart +
				_interpolationConstants[2];
		}
	}
}

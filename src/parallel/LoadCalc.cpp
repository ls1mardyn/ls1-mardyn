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

MeasureLoad::MeasureLoad(bool timeValuesShouldBeIncreasing, int interpolationStartsAt)
	: _timeValuesShouldBeIncreasing{timeValuesShouldBeIncreasing},
	  _interpolationStartsAt{interpolationStartsAt} {
	_previousTime = global_simulation->timers()->getTime(TIMER_NAME);
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

template<typename iterator1, typename iterator2>
bool isFinite(iterator1 start, iterator2 end) {
	bool isProper = true;
	for(auto iter = start; iter != end; ++iter){
		isProper &= std::isfinite(*iter);
	}
	return isProper;
}

/**
 * Generate an approximate solution vector (using least-squares) based on the following equation:
 * arma_system_matrix * x = arma_rhs
 * The following restrictions are applied:
 * - x has to be component-wise non-negative.
 * - The first elements (up to notIncreasingStartAt) are increasing. The elements after that are not restricted.
 * @param arma_system_matrix The system matrix.
 * @param arma_rhs The right hand side.
 * @param quadraticInterpolationStartAt The index starting at which the quadratic interpolation starts. It is not
 * allowed to be bigger than arma_system_matrix.n_cols.
 * @return The solution as described above.
 */
inline arma::vec getIncreasingSolutionVec(arma::mat arma_system_matrix, const arma::vec& arma_rhs,
								   const int quadraticInterpolationStartAt) {
	// We multiply the matrix increasing_time_values_matrix (a lower-triangular matrix with only ones) to the
	// system matrix, s.t., the solution vector is: [t_0, t_1-t_0, t_2-t_1, ..., t_n-t_n-1]
	// (t_i is the time needed for a cell with i particles.)
	// If we then solve this equation using a non-negative-least squares, we can guarantee that t_i+1 > t-i.
	arma::mat increasing_time_values_matrix(arma_system_matrix.n_cols, arma_system_matrix.n_cols);

	for (int row = 0; row < quadraticInterpolationStartAt; row++) {
		for (int column = 0; column <= row; ++column) {
			increasing_time_values_matrix(row, column) = 1.;
		}
		for (int column = row + 1; column < arma_system_matrix.n_cols; ++column) {
			increasing_time_values_matrix(row, column) = 0.;
		}
	}

	int num_extra = arma_system_matrix.n_cols - quadraticInterpolationStartAt;
	mardyn_assert(num_extra == 0 or num_extra == 3);
	if(num_extra == 3){
		// Starting with notIncreasingStartAt, we want, that the quadratic fit is bigger than the last tn. So we modify
		// the solution vector further. The last three lines will look like this:
		//
		// 0 ... 0   1   0  0
		// 0 ... 0   0   1  0
		// 1 ... 1 -N^2 -N  1
		//
		// N = quadraticInterpolationStartAt.
		// If the last three values of the solution are positive, it can be guaranteed that the time values will further
		// increase and are also bigger than t_{N-1}.
		for (int last_rows_index = 0; last_rows_index < num_extra; last_rows_index++) {
			for (int column = 0; column < quadraticInterpolationStartAt; ++column) {
				increasing_time_values_matrix(quadraticInterpolationStartAt + last_rows_index, column) =
					(last_rows_index == 2 ? 1. : 0.);
			}
		}
		int N = quadraticInterpolationStartAt;

		increasing_time_values_matrix(quadraticInterpolationStartAt + 0, quadraticInterpolationStartAt + 0) = 1;
		increasing_time_values_matrix(quadraticInterpolationStartAt + 0, quadraticInterpolationStartAt + 1) = 0;
		increasing_time_values_matrix(quadraticInterpolationStartAt + 0, quadraticInterpolationStartAt + 2) = 0;

		increasing_time_values_matrix(quadraticInterpolationStartAt + 1, quadraticInterpolationStartAt + 0) = 0;
		increasing_time_values_matrix(quadraticInterpolationStartAt + 1, quadraticInterpolationStartAt + 1) = 1;
		increasing_time_values_matrix(quadraticInterpolationStartAt + 1, quadraticInterpolationStartAt + 2) = 0;

		increasing_time_values_matrix(quadraticInterpolationStartAt + 2, quadraticInterpolationStartAt + 0) = -N * N;
		increasing_time_values_matrix(quadraticInterpolationStartAt + 2, quadraticInterpolationStartAt + 1) = -N;
		increasing_time_values_matrix(quadraticInterpolationStartAt + 2, quadraticInterpolationStartAt + 2) = 1;
	}

	arma_system_matrix = arma_system_matrix * increasing_time_values_matrix;

	arma::vec increasing_cell_time_vec = nnls(arma_system_matrix, arma_rhs);
	// multiply increasing_time_values_matrix to increasing_cell_time_vec to obtain the original cell time
	// values.
	return increasing_time_values_matrix * increasing_cell_time_vec;
}


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

	int interpolationStartsAt = std::min(_interpolationStartsAt, global_maxParticlesP1);

	if (_interpolationStartsAt >= 0) {
		if (numRanks < _interpolationConstants.size() + interpolationStartsAt) {
			Log::global_log->warning()
				<< "MeasureLoad: Not enough processes to sample from (needed _interpolationConstants.size() + interpolationStartsAt: "
				<< _interpolationConstants.size() + interpolationStartsAt << ", numRanks: " << numRanks << ")." << std::endl;
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

		// right hand side = neededTimes
		arma::vec arma_rhs(neededTimes);

		// system matrix is = global_statistics
		arma::mat arma_system_matrix(numRanks, global_maxParticlesP1);
		for (int rank /*row*/ = 0; rank < numRanks; rank++) {
			for (int particleCount /*column*/ = 0; particleCount < global_maxParticlesP1; particleCount++) {
				arma_system_matrix(rank, particleCount) = global_statistics[global_maxParticlesP1 * rank + particleCount];
			}
		}
		if (interpolationStartsAt >= 0) {
			int num_dof = interpolationStartsAt + 3;
			arma::mat quadratic_equation_fit_matrix(global_maxParticlesP1, num_dof);
			/*
			 * The matrix looks like this (example for interpolationStartsAt = 3):
			 * 1 0 0 0   0   0
			 * 0 1 0 0   0   0
			 * 0 0 1 0   0   0
			 * 0 0 0 3^2 3^1 3^0
			 * 0 0 0 4^2 4^1 4^0
			 * 0 0 0 5^2 5^1 5^0
			 * 0 0 0 6^2 6^1 6^0
			 * ...
			 *
			 * The solution vector thus first contains the time values for the cells with < interpolationStartsAt
			 * particles and, afterwards, the interpolation constants a, b and c for which a^2 * N + b * N + c holds.
			 */
			for (int row = 0; row < interpolationStartsAt; ++row) {
				for (int col = 0; col < num_dof; ++col) {
					quadratic_equation_fit_matrix.at(row, col) = (row == col ? 1. : 0.);
				}
			}
			for (int row = interpolationStartsAt; row < global_maxParticlesP1; ++row) {
				for (int col = 0; col < interpolationStartsAt; ++col) {
					quadratic_equation_fit_matrix.at(row, col) = 0.;
				}
				quadratic_equation_fit_matrix.at(row, interpolationStartsAt + 0) = row * row;
				quadratic_equation_fit_matrix.at(row, interpolationStartsAt + 1) = row;
				quadratic_equation_fit_matrix.at(row, interpolationStartsAt + 2) = 1.;
			}
			// We multiply the system matrix with quadratic_equation_fit_matrix to be able to get the coefficients of a
			// quadratic least squares fit (non-negative)
			arma_system_matrix = arma_system_matrix * quadratic_equation_fit_matrix;

			arma::vec coefficient_vec;
			if (_timeValuesShouldBeIncreasing) {
				coefficient_vec = getIncreasingSolutionVec(arma_system_matrix, arma_rhs, interpolationStartsAt);
			} else {
				coefficient_vec = nnls(arma_system_matrix, arma_rhs);
			}
			mardyn_assert(coefficient_vec.size() == num_dof);
			global_log->info() << "coefficient_vec: " << std::endl;
			coefficient_vec.raw_print(std::cout);
			std::cout << std::endl;
			_times.resize(interpolationStartsAt);
			for (int i = 0; i < interpolationStartsAt; ++i) {
				_times[i] = coefficient_vec[i];
			}
			MPI_Bcast(_times.data(), interpolationStartsAt, MPI_DOUBLE, 0, comm);
			for (int i = 0; i < 3; i++) {
				_interpolationConstants[i] = coefficient_vec[interpolationStartsAt + i];
			}
			MPI_Bcast(_interpolationConstants.data(), 3, MPI_DOUBLE, 0, comm);
		} else if (_timeValuesShouldBeIncreasing) {
			arma::vec cell_time_vec = getIncreasingSolutionVec(arma_system_matrix, arma_rhs, global_maxParticlesP1);

			global_log->info() << "cell_time_vec:" << std::endl;
			cell_time_vec.raw_print(std::cout);
			_times = arma::conv_to<std::vector<double> >::from(cell_time_vec);
			mardyn_assert(_times.size() == global_maxParticlesP1);
			MPI_Bcast(_times.data(), global_maxParticlesP1, MPI_DOUBLE, 0, comm);
		} else {
			arma::vec cell_time_vec = nnls(arma_system_matrix, arma_rhs);

			global_log->info() << "cell_time_vec:\n" << cell_time_vec << std::endl;
			_times = arma::conv_to<std::vector<double> >::from(cell_time_vec);
			mardyn_assert(_times.size() == global_maxParticlesP1);
			MPI_Bcast(_times.data(), global_maxParticlesP1, MPI_DOUBLE, 0, comm);
		}
	} else {
		MPI_Gather(statistics.data(), statistics.size(), MPI_UINT64_T, nullptr, 0 /*here insignificant*/, MPI_UINT64_T,
				0, comm);
		if (interpolationStartsAt >= 0) {
			_times.resize(interpolationStartsAt);
			MPI_Bcast(_times.data(), interpolationStartsAt, MPI_DOUBLE, 0, comm);
			MPI_Bcast(_interpolationConstants.data(), 3, MPI_DOUBLE, 0, comm);
		} else {
			_times.resize(global_maxParticlesP1);
			MPI_Bcast(_times.data(), global_maxParticlesP1, MPI_DOUBLE, 0, comm);
		}
	}

	// interpolation constants:
	if (interpolationStartsAt < 0) {
		calcConstants();
	}

	if (not isFinite(_times.begin(), _times.end())) {
		global_log->warning() << "Detected non-finite number in MeasureLoad" << std::endl;
		return 1;
	}

	if (not isFinite(_interpolationConstants.begin(), _interpolationConstants.end())) {
		global_log->warning() << "Detected non-finite number in MeasureLoad" << std::endl;
		return 1;
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

	// We do a least square fit of the last 50% of the data or at least three elements.
	size_t start = std::min(_times.size() - 3, _times.size() / 2);
	size_t numElements = _times.size() - start;

#ifdef MARDYN_ARMADILLO
	// Nx3 matrix:
	arma::mat system_matrix(numElements, 3);

	// Nx1 vector from moments:
	arma::vec rhs(numElements);

	for (size_t row = 0; row < numElements; ++row) {
		rhs[row] = _times[start + row];
		system_matrix.at(row, 0) = static_cast<double>(row * row);
		system_matrix.at(row, 1) = static_cast<double>(row);
		system_matrix.at(row, 2) = 1;
	}

	arma::vec solution = nnls(system_matrix, rhs);

	for (size_t row = 0; row < 3ul; row++) {
		_interpolationConstants[row] = solution[2 - row];
	}
	global_log->info() << "_interpolationConstants: " << std::endl << solution << std::endl;
#endif

}

double MeasureLoad::getValue(int numParticles) const {
	mardyn_assert(numParticles >= 0);
	mardyn_assert(_preparedLoad);

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

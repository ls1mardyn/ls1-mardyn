/*
 * LoadCalc.h
 *
 *  Created on: 10.06.2017
 *      Author: griebel_s
 */

#pragma once

#include <array>
#include <vector>
#include <algorithm>
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <utils/Logger.h>

class DomainDecompBase;

#include "Simulation.h"


class LoadCalc {
public:

	virtual ~LoadCalc() {
	}

	virtual double getOwn(int index1, int index2) const = 0;

	virtual double getFace(int index1, int index2) const = 0;

	virtual double getEdge(int index1, int index2) const = 0;

	virtual double getCorner(int index1, int index2) const = 0;
};

/**
 * Stores the measured time for the vectorization tuner.
 *
 * Extrapolates times when needed
 */
class TunerLoad: public LoadCalc {
public:

	/**
	 * The whole idea of this class is to take the ownership of the measured values of the tuner, so they are given by rvalue reference.
	 * This constructor is used when tuning values are read from a file.
	 */
	TunerLoad(int count1, int count2, std::vector<double>&& ownTime, std::vector<double>&& faceTime,
			std::vector<double>&& edgeTime, std::vector<double>&& cornerTime);

	/**
	 * This constructor is normally used when creating new files.
	 */
	TunerLoad() :
			_count1 { 0 }, _count2 { 0 } {
	}

	/**
	 * These function get the load for inner cell interactions and for the different neighbor types,
	 * using extrapolation if necessary
	 */
	double getOwn(int index1, int index2) const override {
		if (index2 < _count2 && index1 < _count1) {
			return accessVec(_ownTime, index1, index2);
		} else {
			// extrapolation
			return index1 * (index1 - 1) / 2 * _ownConst[0] + index2 * (index2 - 1) / 2 * _ownConst[1]
					+ index2 * index1 * _ownConst[2];
		}
	}

	double getFace(int index1, int index2) const override {
		if (index2 < _count2 && index1 < _count1) {
			return accessVec(_faceTime, index1, index2);
		} else {
			return costsNeighbour(index1, index2, _faceConst);
		}
	}

	double getEdge(int index1, int index2) const override {
		if (index2 < _count2 && index1 < _count1) {
			return accessVec(_edgeTime, index1, index2);
		} else {
			return costsNeighbour(index1, index2, _edgeConst);
		}
	}

	double getCorner(int index1, int index2) const override {
		if (index2 < _count2 && index1 < _count1) {
			return accessVec(_cornerTime, index1, index2);
		} else {
			return costsNeighbour(index1, index2, _cornerConst);
		}
	}

	int getCount1() const noexcept {
		return _count1;
	}

	int getCount2() const noexcept {
		return _count2;
	}

	/**
	 * Writes the given TunerTime object to the given Stream
	 */
	static void write(std::ostream& stream, const TunerLoad& time) {
		stream << "Vectorization Tuner File" << std::endl;
		stream << "own" << std::endl;
		time.writeVec(stream, time._ownTime);
		stream << "face" << std::endl;
		time.writeVec(stream, time._faceTime);
		stream << "edge" << std::endl;
		time.writeVec(stream, time._edgeTime);
		stream << "corner" << std::endl;
		time.writeVec(stream, time._cornerTime);
	}

	/**
	 * Reads the tuner values from a given file and returns the corresponding TunerTimes object
	 */
	static TunerLoad read(std::istream& stream);

private:
	/**
	 * Reads a single two dimensional vector (double values separated by ; )
	 * whose dimensions are stored in the given integers
	 */
	static std::vector<double> readVec(std::istream& in, int& count1, int& count2);

	/**
	 * Extrapolation for vectors storing the neighbor times
	 */
	static double costsNeighbour(int index1, int index2, std::array<double, 3> consts) noexcept {
		return index1 * index1 * consts[0] + index2 * index2 * consts[1] + index1 * index2 * 2 * consts[2];
	}

	/**
	 * Writes the given two dimensional array to the given stream
	 *
	 * The data of each dimension is stored in a line, where the single elements are separated by a ;
	 */
	void writeVec(std::ostream& out, const std::vector<double>& vec) const {
		for (int index1 = 0; index1 < _count1; ++index1) {
			for (int index2 = 0; index2 < _count2 - 1; ++index2) {
				out << accessVec(vec, index1, index2) << ";";
			}
			//no terminating ; for the last entry
			if (vec.size() != 0) {
				out << accessVec(vec, index1, _count2 - 1);
			}
			out << std::endl;
		}
	}

	/**
	 * calculates the time needed for a single interaction for the neighbourType represented by the given vector
	 *
	 * Neighbor and inner cell interactions have to be treated differently since, if they have the same amount of particles
	 * there are different amounts of interactions
	 */
	std::array<double, 3> calcConsts(const std::vector<double>& timeVec, bool inner);

	double accessVec(const std::vector<double>& vec, int index1, int index2) const {
		return vec[index2 + _count2 * index1];
	}

	/*
	 * Number of particles up to which the times where still measured and are stored in the vectors (per type)
	 * This means each vector has the size _count1*_count2.
	 */
	int _count1;
	int _count2;
	/*
	 * Storing the measured times for the different neighbor types.
	 *
	 * These are 2D-arrays
	 */
	std::vector<double> _ownTime;
	std::vector<double> _faceTime;
	std::vector<double> _edgeTime;
	std::vector<double> _cornerTime;

	/*
	 * Store the constants needed for the extrapolation. These represent the time needed for a single transaction,
	 * so to extrapolate then simply multiplies these constants times the number of
	 *
	 * The first value is the constant for interactions between particles of type 1, the second between particles of type 2
	 * and the the third for an interaction between particles of type 1 and 2
	 */
	std::array<double, 3> _ownConst{};
	std::array<double, 3> _faceConst{};
	std::array<double, 3> _edgeConst{};
	std::array<double, 3> _cornerConst{};
};

/**
 * The load estimation function previously used in MarDyn
 */
class TradLoad: public LoadCalc {
public:
	double getOwn(int index1, int index2) const override {
		return (index1 + index2) * (index1 + index2);
	}

	double getFace(int index1, int index2) const override {
		return 0.5 * (index1 + index2) * (index1 + index2);
	}

	double getEdge(int index1, int index2) const override {
		return 0.5 * (index1 + index2) * (index1 + index2);
	}

	double getCorner(int index1, int index2) const override {
		return 0.5 * (index1 + index2) * (index1 + index2);
	}
};


/**
 * This class provides loads by time-measurements of the real simulation.
 * The time needed by each process is measured and based on that the time needed for each cell is reconstructed.
 */
class MeasureLoad: public LoadCalc {
	static std::string TIMER_NAME;

public:
	MeasureLoad(bool alwaysUseInterpolation);

	double getOwn(int index1, int index2) const override {
		return getValue(index1);
	}

	double getFace(int index1, int index2) const override {
		return getValue(index1);
	}

	double getEdge(int index1, int index2) const override {
		return getValue(index1);
	}

	double getCorner(int index1, int index2) const override {
		return getValue(index1);
	}

#ifdef ENABLE_MPI
	int prepareLoads(DomainDecompBase* decomp, MPI_Comm& comm);
#endif
private:
	double getValue(int numParticles) const;
	void calcConstants();

	/// stores the times needed for the cells
	std::vector<double> _times;

	/**
	 * Stores values of the extrapolation of the form y = a x^2 + b x + c.
	 * Hereby, x is the number of particles and y the expected runtime of the simulation.
	 */
	std::array<double, 3> _interpolationConstants{};

	bool _preparedLoad{false};
	double _previousTime{0.};
	bool _alwaysUseInterpolation{true};
};

/**
 * \file
 *
 * TODO: documentation severely outdated
 *
 * \brief A CellProcessor that produces Flop information.
 * \details The LJFlopCounter calculates the number of flops that are required to
 * calculate the LJ force calculation for the given Molecules. Then it calls another
 * CellProcessor and measures its runtime.<br>
 * The output shows the flop counts for each of the parts of the LJ flop calculation,
 * for the current iteration and accumulated over all iterations performed so far.
 * It also prints the runtime of the child CellProcessor and calculates the
 * flops per second for the ChildCellProcessor.<br>
 * <br>
 * If you wish to exclude some parts of the calculation or change the flop counts for
 * an individual part of the calculation, you have to sum up the specific flop counts
 * manually and calculate the flops per second using the provided time.<br>
 * If you wish to calculate the effective flops per second based on the total runtime of
 * the simulation, including the traversal of the CellContainer and other things,
 * you can run the simulation once with the LJFlopCounter, and once without, using the same
 * CellProcessor and the same input for the actual calculation. Then you can take the
 * flop counts and the runtime of the simulation without the LJFlopCounter and calculate
 * the flops per second based on these.<br>
 * <br>
 * The LJFlopCounter can also be used to find out how many of the checked MoleculePairs
 * actually are inside the cutoff radius. Note that this only works for single centered
 * molecules. For this you divide the flop counts for
 * molecule distance calculation and center distance calculation by the weights for a
 * single of these operations to acquire the number of such calculations done. The
 * molecule distances give the number of Molecule pairs checked, the center distances
 * provide the number of Molecule pairs that are close enough to each other.<br>
 * <br>
 * The flop values for each operation are based on the unvectorized calculation
 * in VectorizedLJCellProcessor. Details can be found for each value individually in the
 * class documentation.
 *
 * \author Johannes Heckl
 */

#ifndef FLOPCOUNTER_H_
#define FLOPCOUNTER_H_

#include "CellProcessor.h"
#include <memory>
#include <vector>
#include <stddef.h>
#include <string>
#include <sstream>

#include "molecules/MoleculeForwardDeclaration.h"

class CellDataSoA;
class CellDataSoARMM;

/**
 * \brief
 * \author Johannes Heckl
 * \author Nikola Tchipev
 */
class FlopCounter: public CellProcessor {
public:
	/**
	 * \brief Set up the counter.
	 * \details The CellProcessor is called after counting has been finished and<br>
	 * the runtime of the CellProcessor is measured. The CellProcessor is deleted<br>
	 * when LJFlopCounter is destructed.<br>
	 * \param cutoffRadius Must be the same as the cutoff radius for the CellProcessor.
	 * \param LJcutoffRadius Must be the same as the LJ cutoff radius for the CellProcessor.
	 */
	FlopCounter(double cutoffRadius, double LJcutoffRadius);

	/**
	 * free threadData
	 */
	~FlopCounter();

	/**
	 * \brief Initializes the internal counters.
	 */
	void initTraversal();

	/**
	 * \brief Only pass through to child.
	 */
	void preprocessCell(ParticleCell& cell) {}


	double processSingleMolecule(Molecule* /*m1*/, ParticleCell& /*cell2*/) { return 0.0; }  // why 0.0 flops???

	/**
	 * \brief Count flops for this cell.
	 */
	void processCell(ParticleCell& cell);

        void processCellPair(ParticleCell& c1, ParticleCell& c2, bool sumAll = false);

	/**
	 * \brief Only pass through to child.
	 */
	void postprocessCell(ParticleCell& cell) {}

	/**
	 * \brief Print results.
	 */
	void endTraversal();

	double getTotalFlopCount() const {
		return _totalFlopCount;
	}

	double getTotalMoleculeDistanceFlopCount() const {
		return _currentCounts.getMoleculeDistanceFlops();
	}

	void resetCounters() {
		_currentCounts.clear();
//		_totalCounts.clear();
		_totalFlopCount = 0.;
		_myFlopCount = 0.;
	}

	double getMyFlopCount() const {
		return _myFlopCount;
	}

	void printStats() const;

private:
	template<class ForcePolicy, bool CalculateMacroscopic>
	void _calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2);
	template<class ForcePolicy, bool CalculateMacroscopic>
	void _calculatePairs(const CellDataSoARMM & soa1, const CellDataSoARMM & soa2);

	void handlePair(const Molecule& Mi, const Molecule& Mj,
			bool addMacro = true);

	// used for indices within an array!
	enum PotentialIndices {
		I_LJ = 0,
		I_CHARGE,
		I_CHARGE_DIPOLE,
		I_DIPOLE,
		I_CHARGE_QUADRUPOLE,
		I_DIPOLE_QUADRUPOLE,
		I_QUADRUPOLE,
		NUM_POTENTIALS
	};

	class _PotentialCounts {
	public:
		void init(const std::string& n, int kM, int mM, int sFTM, int sMM) {
			_name = n;
			clear();
			_kernelMultiplier = kM;
			_macroMultiplier = mM;
			_sumForceTorqueMultiplier = sFTM;
			_sumMacroMultiplier = sMM;
		}
		void clear() {
			_numKernelCalls = 0.0;
			_numMacroCalls = 0.0;
		}
		void addPotentialCounts(const _PotentialCounts& pc) {
			_numKernelCalls += pc._numKernelCalls;
			_numMacroCalls += pc._numMacroCalls;
		}
		void addKernelAndMacro(double valueBoth, bool addMacro) {
			_numKernelCalls += valueBoth;
			if (addMacro)
				_numMacroCalls += valueBoth;
		}
		void collCommAppend();
		void collCommGet();
		double getKernelAndMacroFlops() const {
			return _numKernelCalls * _kernelMultiplier + _numMacroCalls * _macroMultiplier;
		}
		double getForceTorqueSums() const {
			return _numKernelCalls * _sumForceTorqueMultiplier;
		}
		double getMacroValueSums() const {
			return _numMacroCalls * _sumMacroMultiplier;
		}
		std::string printNameKernelAndMacroCalls() const {
			std::ostringstream ostr;

			if (_numKernelCalls == 0) { return ostr.str(); } // potential is very likely not present

			ostr << " " << _name
				<< ": kernel calls: " << _numKernelCalls
				<< " macro calls: " << _numMacroCalls
				<< std::endl;
			return ostr.str();
		}

		/* Fields */
		double _numKernelCalls;
		double _numMacroCalls;

		/* Multipliers */
		int _kernelMultiplier;
		int _macroMultiplier;
		int _sumForceTorqueMultiplier;
		int _sumMacroMultiplier;

		/* name */
		std::string _name;
	};

	class _Counts {
	public:
		_Counts();
		void clear() {
			_moleculeDistances = 0;

			for (int i = 0; i < NUM_POTENTIALS; ++i) {
				_potCounts[i].clear();
			}
		}
		void addCounts(const _Counts& c) {
			_moleculeDistances += c._moleculeDistances;

			for (int i = 0; i < NUM_POTENTIALS; ++i) {
				_potCounts[i].addPotentialCounts(c._potCounts[i]);
			}
		}
		void allReduce();
		double sumKernelCalls() const {
			double ret = 0.;
			for (int i = 0; i < NUM_POTENTIALS; ++i) {
				ret += _potCounts[i]._numKernelCalls;
			}
			return ret;
		}
		double sumMacros() const {
			double ret = 0.;
			for (int i = 0; i < NUM_POTENTIALS; ++i) {
				ret += _potCounts[i]._numMacroCalls;
			}
			return ret;
		}
		void print() const;

		double getMoleculeDistanceFlops() const {
			return _moleculeDistances * _distanceMultiplier;
		}
		double getCenterDistanceFlops() const {
			return sumKernelCalls() * _distanceMultiplier;
		}
		double getForceTorqueSumFlops() const {
			double ret = 0.;
			for (int i = 0; i < NUM_POTENTIALS; ++i) {
				ret += _potCounts[i].getForceTorqueSums();
			}
			return ret;
		}
		double getMacroValueSumFlops() const {
			double ret = 0.;
			for (int i = 0; i < NUM_POTENTIALS; ++i) {
				ret += _potCounts[i].getMacroValueSums();
			}
			return ret;
		}
		double getTotalFlops() const {
			double ret = getMoleculeDistanceFlops() + getCenterDistanceFlops();
			for (int i = 0; i < NUM_POTENTIALS; ++i) {
				ret += _potCounts[i].getKernelAndMacroFlops();
			}
			ret += getForceTorqueSumFlops() + getMacroValueSumFlops();
			return ret;
		}
		void addKernelAndMacro(PotentialIndices i, double valueBoth, bool addMacro) {
			_potCounts[i].addKernelAndMacro(valueBoth, addMacro);
		}
		void initPotCounter(PotentialIndices i, const std::string& n, int kM, int mM, int sFTM, int sMM) {
			_potCounts[i].init(n, kM, mM, sFTM, sMM);
		}

		/* Fields which are summed explicitly */
		double _moleculeDistances;
		int _distanceMultiplier;

		_PotentialCounts _potCounts[NUM_POTENTIALS];
	};

	std::vector<_Counts *> _threadData;

	_Counts _currentCounts;
//	_Counts _totalCounts; TODO: is this needed?
	double _totalFlopCount;
	double _myFlopCount;
	bool _synchronized;
};

#endif /* FLOPCOUNTER_H_ */

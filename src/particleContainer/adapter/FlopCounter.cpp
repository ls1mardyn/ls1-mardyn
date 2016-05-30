/**
 * \file
 * \brief A CellProcessor that produces Flop information.
 * \author Johannes Heckl
 */

#include "FlopCounter.h"

#include "particleContainer/ParticleCell.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"

FlopCounter::_Counts::_Counts():
	_moleculeDistances(0.), _distanceMultiplier(8){
	// 3 sub + 3 square + 2 add


//inverse R squared is one, because only 1/(R^2) has to be calculated, while R^2 already is calculated.

	// Kernel: 15 = 1 (inverse R squared) + 8 (compute scale) + 3 (apply scale) + 3 (virial tensor)
	// Macro: 4 = 2 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 0 (torques)
	// sum Macro: 2 (upot + virial) + 0 (RF)
	initPotCounter(I_LJ, "Lennard-Jones", 15, 4, 12, 2);

	// Kernel: 10 = 1 (inverse R squared) + 1 (square root) + 2 (compute scale) + 3 (apply scale) + 3 (virial tensor)
	// Macro: 2 = 0 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 0 (torques)
	// sum Macro: 2 (upot + virial) + 0 (RF)
	initPotCounter(I_CHARGE, "Charge", 10, 2, 12, 2);

	// Kernel: 34 = 1 (inverse R squared) + 1 (square root) + 29 + 3 (virial tensor)
	// Macro: 3 = 1 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 3 (torques)
	// sum Macro: 2 (upot + virial) + 0 (RF)
	initPotCounter(I_CHARGE_DIPOLE, "Charge-Dipole", 34, 3, 15, 2);

	// Kernel: 101 = 1 (inverse R squared) + 1 (square root) + 96 + 3 (virial tensor)
	// Macro: 5 = 3 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 6 (torques)
	// sum Macro: 2 (upot + virial) + 1 (RF)
	initPotCounter(I_DIPOLE, "Dipole", 101, 5, 18, 3);

	// Kernel: 52 = 1 (inverse R squared) + 1 (square root) + 47 + 3 (virial tensor)
	// Macro: 2 = 0 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 3 (torques)
	// sum Macro: 2 (upot + virial) + 0 (RF)
	initPotCounter(I_CHARGE_QUADRUPOLE, "Charge-Quadrupole", 52, 2, 15, 2);

	// Kernel: 121 = 1 (inverse R squared) + 1 (square root) + 116 + 3 (virial tensor)
	// Macro: 2 = 0 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 6 (torques)
	// sum Macro: 2 (upot + virial) + 0 (RF)
	initPotCounter(I_DIPOLE_QUADRUPOLE, "Dipole-Quadrupole", 121, 2, 18, 2);

	// Kernel: 131 = 1 (inverse R squared) + 1 (square root) + 126 + 3 (virial tensor)
	// Macro: 2 = 0 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 6 (torques)
	// sum Macro: 2 (upot + virial) + 0 (RF)
	initPotCounter(I_QUADRUPOLE, "Quadrupole", 131, 2, 18, 2);
}

void FlopCounter::_PotentialCounts::collCommAppend() {
	DomainDecompBase& domainDecomp =  global_simulation->domainDecomposition();
	domainDecomp.collCommAppendDouble(_numKernelCalls);
	domainDecomp.collCommAppendDouble(_numMacroCalls);
}

void FlopCounter::_PotentialCounts::collCommGet() {
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	_numKernelCalls = domainDecomp.collCommGetDouble();
	_numMacroCalls = domainDecomp.collCommGetDouble();
}

void FlopCounter::_Counts::allReduce() {
	DomainDecompBase& domainDecomp =  global_simulation->domainDecomposition();
	domainDecomp.collCommInit(15);

	domainDecomp.collCommAppendDouble(_moleculeDistances);

	for (int i = 0; i < NUM_POTENTIALS; ++i) {
		_potCounts[i].collCommAppend();//adds 2 values each
	}

	domainDecomp.collCommAllreduceSum();
	_moleculeDistances = domainDecomp.collCommGetDouble();

	for (int i = 0; i < NUM_POTENTIALS; ++i) {
		_potCounts[i].collCommGet();
	}

	domainDecomp.collCommFinalize();
}


double FlopCounter::_Counts::process() const {
	using std::endl;
	using Log::global_log;

	global_log->info() << " Molecule distance: " << getMoleculeDistanceFlops() << endl;
	global_log->info() << " Center distance: " << getCenterDistanceFlops() << endl;
	for (int i = 0; i < NUM_POTENTIALS; ++i) {
		std::string str = _potCounts[i].printNameKernelAndMacroFlops();
		if(str.length() > 0) global_log->info() << str;
	}
	global_log->info() << " Force/Torque Sum: " << getForceTorqueSumFlops() << endl;
	global_log->info() << " Macroscopic value sum: " << getMacroValueSumFlops() << endl;
	global_log->info() << "Total: " << getTotalFlops() << endl;

	return getTotalFlops();
}

FlopCounter::FlopCounter(double cutoffRadius, double LJCutoffRadius) : CellProcessor(cutoffRadius, LJCutoffRadius),
		_currentCounts(), _totalFlopCount(0.), _myFlopCount(0.), _synchronized(true){
}

void FlopCounter::initTraversal(const size_t /*numCells*/) {
	_currentCounts.clear();
}


void FlopCounter::endTraversal() {
	_myFlopCount = _currentCounts.getTotalFlops();
	if(_synchronized){
		_currentCounts.allReduce();
	}

	Log::global_log->info()
				<< "FLOP counts in force calculation for this iteration:" << std::endl;
	_totalFlopCount = _currentCounts.process();

//	_totalCounts.addCounts(_currentCounts);
//	Log::global_log->info()
//			<< "Accumulated FLOP counts in force calculation:" << std::endl;
//	_totalFlopCount = _totalCounts.process();
}

void FlopCounter::preprocessCell(ParticleCell & /*c*/) {
}

void FlopCounter::postprocessCell(ParticleCell & /*c*/) {
}

void FlopCounter::handlePair(const Molecule& Mi, const Molecule& Mj, bool addMacro) {
	// Have to compare the distance between 2 molecules.
	_currentCounts._moleculeDistances += 1.;

	double distSquared = 0.0;
	for (int dim = 0; dim < 3; ++dim) {
		const double dist = Mi.r(dim) - Mj.r(dim);
		distSquared += dist * dist;
	}

	const size_t numLJcenters_i 	= Mi.numLJcenters();
	const size_t numCharges_i 		= Mi.numCharges();
	const size_t numDipoles_i 		= Mi.numDipoles();
	const size_t numQuadrupoles_i 	= Mi.numQuadrupoles();

	const size_t numLJcenters_j 	= Mj.numLJcenters();
	const size_t numCharges_j 		= Mj.numCharges();
	const size_t numDipoles_j 		= Mj.numDipoles();
	const size_t numQuadrupoles_j 	= Mj.numQuadrupoles();


	if (distSquared < _LJCutoffRadiusSquare) {
		_currentCounts.addKernelAndMacro(I_LJ, numLJcenters_i * numLJcenters_j, addMacro);
	}

	if (distSquared < _cutoffRadiusSquare) {
		_currentCounts.addKernelAndMacro(I_CHARGE, numCharges_i * numCharges_j, addMacro);

		_currentCounts.addKernelAndMacro(I_CHARGE_DIPOLE, numCharges_i * numDipoles_j + numDipoles_i * numCharges_j, addMacro);
		_currentCounts.addKernelAndMacro(I_DIPOLE, numDipoles_i * numDipoles_j, addMacro);

		_currentCounts.addKernelAndMacro(I_CHARGE_QUADRUPOLE, numCharges_i * numQuadrupoles_j + numQuadrupoles_i * numCharges_j, addMacro);
		_currentCounts.addKernelAndMacro(I_DIPOLE_QUADRUPOLE, numDipoles_i * numQuadrupoles_j + numQuadrupoles_i * numDipoles_j, addMacro);
		_currentCounts.addKernelAndMacro(I_QUADRUPOLE, numQuadrupoles_i * numQuadrupoles_j, addMacro);
	}

}

void FlopCounter::processCell(ParticleCell & c) {
	using std::vector;

	// we don't execute any flops if cell is a halo cell
	if (c.isHaloCell())
		return;

	const vector<Molecule *> & molecules = c.getParticlePointers();
	const vector<Molecule *>::size_type numMolecules = molecules.size();

	if (numMolecules == 0) {
		return;
	}

	// Macroscopic values are always added for processCell
	const bool addMacro = true;

	for (vector<Molecule *>::size_type i = 0; i < numMolecules; ++i) {
		const Molecule& Mi = *(molecules[i]);
		for (vector<Molecule *>::size_type j = i + 1; j < numMolecules; ++j) {
			const Molecule& Mj = *(molecules[j]);
			handlePair(Mi, Mj, addMacro);
		}
	}
}

void FlopCounter::processCellPair(ParticleCell & c1, ParticleCell & c2) {
	using std::vector;

	// we don't execute any flops if both cells are halo cells
	if (c1.isHaloCell() and c2.isHaloCell())
		return;

	const vector<Molecule *> & molecules1 = c1.getParticlePointers();
	const vector<Molecule *> & molecules2 = c2.getParticlePointers();

	const vector<Molecule *>::size_type numMolecules1 = molecules1.size();
	const vector<Molecule *>::size_type numMolecules2 = molecules2.size();

	if ((numMolecules1 == 0) or (numMolecules2 == 0)) {
		return;
	}

	// if both cells are non-halo cells, all macros are added,
	// otherwise an "isLessThan" check is needed
	const bool allMacro = (!c1.isHaloCell()) and (!c2.isHaloCell());

	for (vector<Molecule *>::size_type i = 0; i < numMolecules1; ++i) {
		const Molecule & Mi = *(molecules1[i]);
		for (vector<Molecule *>::size_type j = 0; j < numMolecules2; ++j) {
			const Molecule & Mj = *(molecules2[j]);
			const bool addMacro = allMacro or c1.getCellIndex() < c2.getCellIndex();
			handlePair(Mi, Mj, addMacro);
		}
	}
}

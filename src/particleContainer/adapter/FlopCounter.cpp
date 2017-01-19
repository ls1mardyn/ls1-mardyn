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

FlopCounter::FlopCounter(double cutoffRadius, double LJcutoffRadius) : _cutoffRadiusSquare(cutoffRadius * cutoffRadius), _LJcutoffRadiusSquare(LJcutoffRadius * LJcutoffRadius) {
	_totalCounts.clear();
}

void FlopCounter::initTraversal(const size_t numCells) {
	_currentCounts.clear();
}

void FlopCounter::endTraversal() {

	DomainDecompBase& domainDecomp =  global_simulation->domainDecomposition();
	domainDecomp.collCommInit(15);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_molDist);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_LJ);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_Charges);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_Dipoles);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_ChargesDipoles);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_LJMacro);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_ChargesMacro);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_DipolesMacro);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_ChargesDipolesMacro);

	domainDecomp.collCommAppendDouble(_currentCounts.calc_Quadrupoles);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_ChargesQuadrupoles);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_DipolesQuadrupoles);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_QuadrupolesMacro);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_ChargesQuadrupolesMacro);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_DipolesQuadrupolesMacro);
	domainDecomp.collCommAllreduceSum();
	_currentCounts.calc_molDist = domainDecomp.collCommGetDouble();
	_currentCounts.calc_LJ = domainDecomp.collCommGetDouble();
	_currentCounts.calc_Charges = domainDecomp.collCommGetDouble();
	_currentCounts.calc_Dipoles = domainDecomp.collCommGetDouble();
	_currentCounts.calc_ChargesDipoles = domainDecomp.collCommGetDouble();
	_currentCounts.calc_LJMacro = domainDecomp.collCommGetDouble();
	_currentCounts.calc_ChargesMacro = domainDecomp.collCommGetDouble();
	_currentCounts.calc_DipolesMacro = domainDecomp.collCommGetDouble();
	_currentCounts.calc_ChargesDipolesMacro = domainDecomp.collCommGetDouble();

	_currentCounts.calc_Quadrupoles = domainDecomp.collCommGetDouble();
	_currentCounts.calc_ChargesQuadrupoles = domainDecomp.collCommGetDouble();
	_currentCounts.calc_DipolesQuadrupoles = domainDecomp.collCommGetDouble();
	_currentCounts.calc_QuadrupolesMacro = domainDecomp.collCommGetDouble();
	_currentCounts.calc_ChargesQuadrupolesMacro = domainDecomp.collCommGetDouble();
	_currentCounts.calc_DipolesQuadrupolesMacro = domainDecomp.collCommGetDouble();
	domainDecomp.collCommFinalize();

	_totalCounts.addCounts(_currentCounts);

	const double cflMolDist = _currentCounts.calc_molDist * _flops_MolDist;
	const double cflCenterDist = (_currentCounts.calc_LJ + _currentCounts.calc_Charges + _currentCounts.calc_Dipoles + _currentCounts.calc_ChargesDipoles +
			_currentCounts.calc_Quadrupoles + _currentCounts.calc_ChargesQuadrupoles + _currentCounts.calc_DipolesQuadrupoles) * _flops_CenterDist;
	const double cflLJKernel = _currentCounts.calc_LJ * _flops_LJKernel;
	const double cflChargesKernel = _currentCounts.calc_Charges * _flops_ChargesKernel;
	const double cflDipolesKernel = _currentCounts.calc_Dipoles * _flops_DipolesKernel;
	const double cflChargesDipolesKernel = _currentCounts.calc_ChargesDipoles * _flops_ChargesDipolesKernel;
	const double cflQuadrupolesKernel = _currentCounts.calc_Quadrupoles * _flops_QuadrupolesKernel;
	const double cflChargesQuadrupolesKernel = _currentCounts.calc_ChargesQuadrupoles * _flops_ChargesQuadrupolesKernel;
	const double cflDipolesQuadrupolesKernel = _currentCounts.calc_DipolesQuadrupoles * _flops_DipolesQuadrupolesKernel;
	const double cflSum = (_currentCounts.calc_LJ*2. + _currentCounts.calc_Charges*2. + _currentCounts.calc_Dipoles*4. + _currentCounts.calc_ChargesDipoles*3. +
			_currentCounts.calc_Quadrupoles*4. + _currentCounts.calc_ChargesQuadrupoles*3. + _currentCounts.calc_DipolesQuadrupoles*4. ) * _flops_ForcesSum;
	const double cflMacro = _currentCounts.calc_LJMacro * _flops_LJMacroValues + _currentCounts.calc_ChargesMacro * _flops_ChargesMacroValues + _currentCounts.calc_DipolesMacro * _flops_DipolesMacroValues + _currentCounts.calc_ChargesDipolesMacro * _flops_ChargesDipolesMacroValues +
			_currentCounts.calc_QuadrupolesMacro * _flops_QuadrupolesMacroValues + _currentCounts.calc_ChargesQuadrupolesMacro * _flops_ChargesQuadrupolesMacroValues + _currentCounts.calc_DipolesQuadrupolesMacro * _flops_DipolesQuadrupolesMacroValues;
	const double cflMacroSum = (_currentCounts.calc_LJMacro + _currentCounts.calc_ChargesMacro + _currentCounts.calc_DipolesMacro + _currentCounts.calc_ChargesDipolesMacro +
			_currentCounts.calc_QuadrupolesMacro + _currentCounts.calc_ChargesQuadrupolesMacro + _currentCounts.calc_DipolesQuadrupolesMacro) * _flops_MacroSum + _currentCounts.calc_DipolesMacro * _flops_MacroSumRF;

	const double cflTotal = cflMolDist + cflCenterDist + cflLJKernel + cflChargesKernel + cflDipolesKernel + cflChargesDipolesKernel +
			cflQuadrupolesKernel + cflChargesQuadrupolesKernel + cflDipolesQuadrupolesKernel + cflSum + cflMacro + cflMacroSum;

	Log::global_log->info()
			<< "FLOP counts in force calculation for this iteration:" << std::endl
			<< " Molecule distance: " << cflMolDist
			<< " Center distance: " << cflCenterDist
			<< " LJ Kernel: " << cflLJKernel
			<< " Charges Kernel: " << cflChargesKernel
			<< " Dipoles Kernel: " << cflDipolesKernel
			<< " Charges-Dipoles Kernel: " << cflChargesDipolesKernel
			<< " Quadrupoles Kernel: " << cflQuadrupolesKernel
			<< " Charges-Quadrupoles Kernel: " << cflChargesQuadrupolesKernel
			<< " Dipoles-Quadrupoles Kernel: " << cflDipolesQuadrupolesKernel
			<< " Force/Torque Sum: "  << cflSum
			<< " Macroscopic values: " << cflMacro
			<< " Macroscopic value sum: " << cflMacroSum << std::endl
			<< "Current total FLOPS: " << cflTotal << std::endl;

	const double flMolDist = _totalCounts.calc_molDist * _flops_MolDist;
	const double flCenterDist = (_totalCounts.calc_LJ + _totalCounts.calc_Charges + _totalCounts.calc_Dipoles + _totalCounts.calc_ChargesDipoles +
			_totalCounts.calc_Quadrupoles + _totalCounts.calc_ChargesQuadrupoles + _totalCounts.calc_DipolesQuadrupoles) * _flops_CenterDist;
	const double flLJKernel = _totalCounts.calc_LJ * _flops_LJKernel;
	const double flChargesKernel = _totalCounts.calc_Charges * _flops_ChargesKernel;
	const double flDipolesKernel = _totalCounts.calc_Dipoles * _flops_DipolesKernel;
	const double flChargesDipolesKernel = _totalCounts.calc_ChargesDipoles * _flops_ChargesDipolesKernel;
	const double flQuadrupolesKernel = _totalCounts.calc_Quadrupoles * _flops_QuadrupolesKernel;
	const double flChargesQuadrupolesKernel = _totalCounts.calc_ChargesQuadrupoles * _flops_ChargesQuadrupolesKernel;
	const double flDipolesQuadrupolesKernel = _totalCounts.calc_DipolesQuadrupoles * _flops_DipolesQuadrupolesKernel;
	const double flSum = (_totalCounts.calc_LJ*2. + _totalCounts.calc_Charges*2. + _totalCounts.calc_Dipoles*4. + _totalCounts.calc_ChargesDipoles*3. +
			_totalCounts.calc_Quadrupoles*4. + _totalCounts.calc_ChargesQuadrupoles*3. + _totalCounts.calc_DipolesQuadrupoles*4. ) * _flops_ForcesSum;
	const double flMacro = _totalCounts.calc_LJMacro * _flops_LJMacroValues + _totalCounts.calc_ChargesMacro * _flops_ChargesMacroValues + _totalCounts.calc_DipolesMacro * _flops_DipolesMacroValues + _totalCounts.calc_ChargesDipolesMacro * _flops_ChargesDipolesMacroValues +
			_totalCounts.calc_QuadrupolesMacro * _flops_QuadrupolesMacroValues + _totalCounts.calc_ChargesQuadrupolesMacro * _flops_ChargesQuadrupolesMacroValues + _totalCounts.calc_DipolesQuadrupolesMacro * _flops_DipolesQuadrupolesMacroValues;
	const double flMacroSum = (_totalCounts.calc_LJMacro + _totalCounts.calc_ChargesMacro + _totalCounts.calc_DipolesMacro + _totalCounts.calc_ChargesDipolesMacro +
			_totalCounts.calc_QuadrupolesMacro + _totalCounts.calc_ChargesQuadrupolesMacro + _totalCounts.calc_DipolesQuadrupolesMacro) * _flops_MacroSum + _totalCounts.calc_DipolesMacro * _flops_MacroSumRF;

	_totalFlopCount = flMolDist + flCenterDist + flLJKernel + flChargesKernel + flDipolesKernel + flChargesDipolesKernel +
			flQuadrupolesKernel + flChargesQuadrupolesKernel + flDipolesQuadrupolesKernel + flSum + flMacro + flMacroSum;

	Log::global_log->info()
			<< "Accumulated FLOP counts in force calculation for this iteration:" << std::endl
			<< " Molecule distance: " << flMolDist
			<< " Center distance: " << flCenterDist
			<< " LJ Kernel: " << flLJKernel
			<< " Charges Kernel: " << flChargesKernel
			<< " Dipoles Kernel: " << flDipolesKernel
			<< " Charges-Dipoles Kernel: " << flChargesDipolesKernel
			<< " Quadrupoles Kernel: " << flQuadrupolesKernel
			<< " Charges-Quadrupoles Kernel: " << flChargesQuadrupolesKernel
			<< " Dipoles-Quadrupoles Kernel: " << flDipolesQuadrupolesKernel
			<< " Force/Torque Sum: "  << flSum
			<< " Macroscopic values: " << flMacro
			<< " Macroscopic value sum: " << flMacroSum << std::endl
			<< "Accumulated total FLOPS: " << _totalFlopCount << std::endl;
}

void FlopCounter::preprocessCell(ParticleCell & c) {
}

void FlopCounter::postprocessCell(ParticleCell & c) {
}

void FlopCounter::processCell(ParticleCell & c) {
	// don't count halo cells
	if (c.isHaloCell())
		return;

	const MoleculeList & molecules = c.getParticlePointers();
	if (molecules.size() > 1) {
		const MoleculeList::const_iterator end_i = --(molecules.end());
		const MoleculeList::const_iterator end_j = molecules.end();

		for (MoleculeList::const_iterator i = molecules.begin(); i != end_i; ++i) {
			for (MoleculeList::const_iterator j = molecules.begin(); j != end_j; ++j) {

				// Have to compare the distance between 2 molecules.
				_currentCounts.calc_molDist += 1;

				const double d_x = (*i)->r(0) - (*j)->r(0);
				const double d_y = (*i)->r(1) - (*j)->r(1);
				const double d_z = (*i)->r(2) - (*j)->r(2);
				const double d2 = d_x * d_x + d_y * d_y + d_z * d_z;

				const size_t numLJcenters_i = (*i)->numLJcenters();
				const size_t numLJcenters_j = (*j)->numLJcenters();
				const size_t numCharges_i = (*i)->numCharges();
				const size_t numCharges_j = (*j)->numCharges();
				const size_t numDipoles_i = (*i)->numDipoles();
				const size_t numDipoles_j = (*j)->numDipoles();
				const size_t numQuadrupoles_i = (*i)->numQuadrupoles();
				const size_t numQuadrupoles_j = (*j)->numQuadrupoles();


				if (j > i) {
					if (d2 < _LJcutoffRadiusSquare) {
						const size_t numLJcenterPairs = numLJcenters_i * numLJcenters_j;
						// Have to calculate the LJ force for each pair of centers.
						_currentCounts.calc_LJ += numLJcenterPairs;
						// Have to calculate macroscopic values for each pair of centers.
						_currentCounts.calc_LJMacro += numLJcenterPairs;
					}

					if (d2 < _cutoffRadiusSquare) {
						const size_t numChargesPairs = numCharges_i * numCharges_j;
						const size_t numDipolesPairs = numDipoles_i * numDipoles_j;
						const size_t numQuadrupolesPairs = numQuadrupoles_i * numQuadrupoles_j;

						// Have to calculate the charge force for each pair of centers.
						_currentCounts.calc_Charges += numChargesPairs;
						_currentCounts.calc_Dipoles += numDipolesPairs;
						_currentCounts.calc_Quadrupoles += numQuadrupolesPairs;
						// Have to calculate macroscopic values for each pair of centers.
						_currentCounts.calc_ChargesMacro += numChargesPairs;
						_currentCounts.calc_DipolesMacro += numDipolesPairs;
						_currentCounts.calc_QuadrupolesMacro += numQuadrupolesPairs;
					}
				}
				if (i != j)
				{
					if (d2 < _cutoffRadiusSquare) {
							const size_t numChargesDipolesPairs = numCharges_i * numDipoles_j;
							const size_t numChargesQuadrupolesPairs = numCharges_i * numQuadrupoles_j;
							const size_t numDipolesQuadrupolesPairs = numDipoles_i + numQuadrupoles_j;

							// Have to calculate the charge force for each pair of centers.
							_currentCounts.calc_ChargesDipoles += numChargesDipolesPairs;
							_currentCounts.calc_ChargesQuadrupoles += numChargesQuadrupolesPairs;
							_currentCounts.calc_DipolesQuadrupoles += numDipolesQuadrupolesPairs;
							// Have to calculate macroscopic values for each pair of centers.
							_currentCounts.calc_ChargesDipolesMacro += numChargesDipolesPairs;
							_currentCounts.calc_ChargesQuadrupolesMacro += numChargesQuadrupolesPairs;
							_currentCounts.calc_DipolesQuadrupolesMacro += numDipolesQuadrupolesPairs;
					}
				}
			}
		}
	}
}

void FlopCounter::processCellPair(ParticleCell & c1, ParticleCell & c2) {
	// don't count halo pairs
	if (c1.isHaloCell() and c2.isHaloCell())
		return;

	const MoleculeList & molecules1 = c1.getParticlePointers();
	const MoleculeList & molecules2 = c2.getParticlePointers();
	if ((molecules1.size() > 0) && (molecules2.size() > 0)) {
		const MoleculeList::const_iterator end_i = molecules1.end();
		const MoleculeList::const_iterator end_j = molecules2.end();

		for (MoleculeList::const_iterator i = molecules1.begin(); i != end_i;
				++i) {
			for (MoleculeList::const_iterator j = molecules2.begin();
					j != end_j; ++j) {

				// Have to compare the distance between 2 molecules.
				_currentCounts.calc_molDist += 1;

				const double d_x = (*i)->r(0) - (*j)->r(0);
				const double d_y = (*i)->r(1) - (*j)->r(1);
				const double d_z = (*i)->r(2) - (*j)->r(2);
				const double d2 = d_x * d_x + d_y * d_y + d_z * d_z;

				if (d2 < _LJcutoffRadiusSquare) {
					const size_t numLJcenters_i = (*i)->numLJcenters();
					const size_t numLJcenters_j = (*j)->numLJcenters();
					const size_t numLJcenterPairs = numLJcenters_i * numLJcenters_j;

					// Have to calculate the LJ force for each pair of centers.
					_currentCounts.calc_LJ += numLJcenterPairs;

					if ((c1.isHaloCell() == (!c2.isHaloCell()))
							&& ((*i)->isLessThan(**j))) {
						// Have to calculate macroscopic values for each pair of centers.
						_currentCounts.calc_LJMacro += numLJcenterPairs;
					}
				}

				if (d2 < _cutoffRadiusSquare) {
					const size_t numCharges_i = (*i)->numCharges();
					const size_t numCharges_j = (*j)->numCharges();
					const size_t numDipoles_i = (*i)->numDipoles();
					const size_t numDipoles_j = (*j)->numDipoles();
					const size_t numQuadrupoles_i = (*i)->numQuadrupoles();
					const size_t numQuadrupoles_j = (*j)->numQuadrupoles();
					const size_t numChargesPairs = numCharges_i * numCharges_j;
					const size_t numDipolesPairs = numDipoles_i * numDipoles_j;
					const size_t numQuadrupolesPairs = numQuadrupoles_i * numQuadrupoles_j;
					const size_t numChargesDipolesPairs = numCharges_i * numDipoles_j + numDipoles_i * numCharges_j;
					const size_t numChargesQuadrupolesPairs = numCharges_i * numQuadrupoles_j + numQuadrupoles_i * numCharges_j;
					const size_t numDipolesQuadrupolesPairs = numDipoles_i + numQuadrupoles_j + numQuadrupoles_i + numDipoles_j;

					// Have to calculate the charge force for each pair of centers.
					_currentCounts.calc_Charges += numChargesPairs;
					_currentCounts.calc_Dipoles += numDipolesPairs;
					_currentCounts.calc_Quadrupoles += numQuadrupolesPairs;
					_currentCounts.calc_ChargesDipoles += numChargesDipolesPairs;
					_currentCounts.calc_ChargesQuadrupoles += numChargesQuadrupolesPairs;
					_currentCounts.calc_DipolesQuadrupoles += numDipolesQuadrupolesPairs;


					if ((c1.isHaloCell() == (!c2.isHaloCell()))
							&& ((*i)->isLessThan(**j))) {
						// Have to calculate macroscopic values for each pair of centers.
						_currentCounts.calc_ChargesMacro += numChargesPairs;
						_currentCounts.calc_DipolesMacro += numDipolesPairs;
						_currentCounts.calc_QuadrupolesMacro += numQuadrupolesPairs;
						_currentCounts.calc_ChargesDipolesMacro += numChargesDipolesPairs;
						_currentCounts.calc_ChargesQuadrupolesMacro += numChargesQuadrupolesPairs;
						_currentCounts.calc_DipolesQuadrupolesMacro += numDipolesQuadrupolesPairs;
					}
				}
			}
		}
	}
}

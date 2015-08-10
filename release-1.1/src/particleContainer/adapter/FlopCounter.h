/**
 * \file
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

#ifndef LJFLOPCOUNTER_H_
#define LJFLOPCOUNTER_H_

#include "CellProcessor.h"
#include <memory>
#include <vector>
#include <stddef.h>
#include "../../utils/Timer.h"

class Molecule;

/**
 * \brief
 * \author Johannes Heckl
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
	 * \brief Initializes the internal counters.
	 */
	void initTraversal(const size_t);

	/**
	 * \brief Only pass through to child.
	 */
	void preprocessCell(ParticleCell& cell);

	/**
	 * \brief Count flops for this pair.
	 */
	void processCellPair(ParticleCell& cell1, ParticleCell& cell2);

	double processSingleMolecule(Molecule* m1, ParticleCell& cell2) { return 0.0; }

	/**
	 * \brief Count flops for this cell.
	 */
	void processCell(ParticleCell& cell);

	/**
	 * \brief Only pass through to child.
	 */
	void postprocessCell(ParticleCell& cell);

	/**
	 * \brief Print results.
	 */
	void endTraversal();

	double getTotalFlopCount() const {
		return _totalFlopCount;
	}

private:
	// The cutoff radius squared.
	const double _cutoffRadiusSquare;
	const double _LJcutoffRadiusSquare;

	typedef std::vector<Molecule *> MoleculeList;

	// Collection of counts that are tracked.
	class _Counts {
	public:
		void clear() {
			calc_molDist = 0;
			calc_LJ = 0;
			calc_Charges = 0;
			calc_Dipoles = 0;
			calc_ChargesDipoles = 0;
			calc_LJMacro = 0;
			calc_ChargesMacro = 0;
			calc_DipolesMacro = 0;
			calc_ChargesDipolesMacro = 0;

			calc_Quadrupoles = 0;
			calc_ChargesQuadrupoles = 0;
			calc_DipolesQuadrupoles = 0;
			calc_QuadrupolesMacro = 0;
			calc_ChargesQuadrupolesMacro = 0;
			calc_DipolesQuadrupolesMacro = 0;


		}
		void addCounts(const _Counts c) {
			calc_molDist += c.calc_molDist;
			calc_LJ += c.calc_LJ;
			calc_Charges += c.calc_Charges;
			calc_Dipoles += c.calc_Dipoles;
			calc_ChargesDipoles += c.calc_ChargesDipoles;
			calc_LJMacro += c.calc_LJMacro;
			calc_ChargesMacro += c.calc_ChargesMacro;
			calc_DipolesMacro += c.calc_ChargesMacro;
			calc_ChargesDipolesMacro += c.calc_ChargesDipolesMacro;

			calc_Quadrupoles += c.calc_Quadrupoles;
			calc_ChargesQuadrupoles += c.calc_ChargesQuadrupoles;
			calc_DipolesQuadrupoles += c.calc_DipolesQuadrupoles;
			calc_QuadrupolesMacro += c.calc_QuadrupolesMacro;
			calc_ChargesQuadrupolesMacro += c.calc_ChargesQuadrupolesMacro;
			calc_DipolesQuadrupolesMacro += c.calc_DipolesQuadrupolesMacro;
		}
		double calc_molDist;
		double calc_LJ;
		double calc_Charges;
		double calc_Dipoles;
		double calc_ChargesDipoles;
		double calc_LJMacro;
		double calc_ChargesMacro;
		double calc_DipolesMacro;
		double calc_ChargesDipolesMacro;

		double calc_Quadrupoles;
		double calc_ChargesQuadrupoles;
		double calc_DipolesQuadrupoles;
		double calc_QuadrupolesMacro;
		double calc_ChargesQuadrupolesMacro;
		double calc_DipolesQuadrupolesMacro;
	};

	_Counts _totalCounts;
	_Counts _currentCounts;
	double _totalFlopCount;

	/**
	 * \brief The calculation of the distance between 2 molecules.
	 * \details 3 + 3 + 2 for 3 differences between molecule coordinates,
	 * squaring these differences, and adding the 3 values together with 2 additions.
	 */
	static const size_t _flops_MolDist = 8;
	/**
	 * \brief The calculation of the distance between 2 centers.
	 * \details Same as _flops_MolDist for distances between centers.
	 */
	static const size_t _flops_CenterDist = 8;
	/**
	 * \brief The calculation of the force between 2 LJ centers.
	 * \details 1 for inverse of r squared, 8 for calculation of scale, 3 for applying the
	 * scale.
	 */
	static const size_t _flops_LJKernel = 12;
	/**
	 * \brief The calculation of the force between 2 charges.
	 * \details 1 for inverse of r squared, 1 for square root, 2 for calculation of scale, 3 for applying the
	 * scale.
	 */
	static const size_t _flops_ChargesKernel = 7;
	/**
	 * \brief The calculation of the force between 2 dipoles.
	 * \details 1 for inverse of r squared, 1 for square root, 96 for calculation of forces and torques
	 * scale.
	 */
	static const size_t _flops_DipolesKernel = 98;
	/**
	 * \brief The calculation of the force between 1 charge and 1 dipole.
	 * \details 1 for inverse of r squared, 1 for square root, 29 for calculation forces and torque
	 * scale.
	 */
	static const size_t _flops_ChargesDipolesKernel = 31;
	/**
	 * \brief The calculation of the force between 2 quadrupoles.
	 * \details 1 for inverse of r squared, 1 for square root, 126 for calculation of forces and torques
	 * scale.
	 */
	static const size_t _flops_QuadrupolesKernel = 128;
	/**
	 * \brief The calculation of the force between 1 charge and 1 quadrupole.
	 * \details 1 for inverse of r squared, 1 for square root, 47 for calculation forces and torque
	 * scale.
	 */
	static const size_t _flops_ChargesQuadrupolesKernel = 49;
	/**
	 * \brief The calculation of the force between 1 charge and 1 quadrupole.
	 * \details 1 for inverse of r squared, 1 for square root, 116 for calculation forces and torque
	 * scale.
	 */
	static const size_t _flops_DipolesQuadrupolesKernel = 118;
	/**
	 * \brief Summing up the forces.
	 * \details 3 for adding or subtracting the resulting force or torque to one molecule
	 * the other.
	 */
	static const size_t _flops_ForcesSum = 3;
	/**
	 * \brief The calculation of Virial and UPot for LJ centers.
	 * \details 2 for upot, 5 for virial.
	 */
	static const size_t _flops_LJMacroValues = 7;
	/**
	 * \brief The calculation of Virial and UPot for charges.
	 * \details 0 for upot, 5 for virial.
	 */
	static const size_t _flops_ChargesMacroValues = 5;
	/**
	 * \brief The calculation of Virial and UPot and reaction field for dipoles.
	 * \details 3 for upot, 5 for virial, 1 for rf.
	 */
	static const size_t _flops_DipolesMacroValues = 9;
	/**
	 * \brief The calculation of Virial and UPot
	 * \details 1 for upot, 5 for virial
	 */
	static const size_t _flops_ChargesDipolesMacroValues = 6;
	/**
	 * \brief The calculation of Virial
	 * \details 5 for virial.
	 */
	static const size_t _flops_QuadrupolesMacroValues = 5;
	/**
	 * \brief The calculation of Virial
	 * \details 5 for virial.
	 */
	static const size_t _flops_ChargesQuadrupolesMacroValues = 5;
	/**
	 * \brief The calculation of Virial
	 * \details 5 for virial.
	 */
	static const size_t _flops_DipolesQuadrupolesMacroValues = 5;
	/**
	 * \brief Summing up Virial and UPot.
	 * \details 2 additions.
	 */
	static const size_t _flops_MacroSum = 2;
	/**
	 * \brief Summing up rf.
	 * \details 1 addition.
	 */
	static const size_t _flops_MacroSumRF = 1;
};

#endif

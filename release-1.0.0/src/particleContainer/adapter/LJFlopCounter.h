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
class LJFlopCounter: public CellProcessor {
public:
	/**
	 * \brief Set up the counter.
	 * \details The CellProcessor is called after counting has been finished and<br>
	 * the runtime of the CellProcessor is measured. The CellProcessor is deleted<br>
	 * when LJFlopCounter is destructed.<br>
	 * \param cutoffRadius Must be the same as the LJ cutoff radius for the CellProcessor.
	 */
	LJFlopCounter(double cutoffRadius);

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

	double processSingleMolecule(Molecule* m1, ParticleCell& cell2);

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
	const double _rc2;

	typedef std::vector<Molecule *> MoleculeList;

	// Collection of counts that are tracked.
	class _Counts {
	public:
		void clear() {
			calc_LJ = 0;
			calc_Macro = 0;
			calc_molDist = 0;
		}
		void addCounts(const _Counts c) {
			calc_LJ += c.calc_LJ;
			calc_Macro += c.calc_Macro;
			calc_molDist += c.calc_molDist;
		}
		double calc_molDist;
		double calc_LJ;
		double calc_Macro;
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
	 * \brief The calculation of the distance between 2 LJ centers.
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
	 * \brief Summing up the LJ forces.
	 * \details 3 for adding the resulting force to one molecule, 3 for subtracting from
	 * the other.
	 */
	static const size_t _flops_LJSum = 6;
	/**
	 * \brief The calculation of Virial and UPotLJ.
	 * \details 2 for upot, 5 for virial.
	 */
	static const size_t _flops_MacroValues = 7;
	/**
	 * \brief Summing up Virial and UPotLJ.
	 * \details 2 additions.
	 */
	static const size_t _flops_MacroSum = 2;
};

#endif

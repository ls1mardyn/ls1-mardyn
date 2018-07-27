#ifndef SRC_MOLECULES_MOLECULEIDPOOL_H_
#define SRC_MOLECULES_MOLECULEIDPOOL_H_

#include "utils/mardyn_assert.h"

/** @brief The MoleculeIdPool manages molecule ID handling.
*
* The MoleculeIdPool divides the available range of molecule IDs across processes.
* It keeps track of how many molecule IDs were obtained by the process.
*
* @note this is not a singleton and you have to ensure to use the same pool object yourself.
* @todo this scheme can be used to improve the check for duplicate IDs after the run.
*/
class MoleculeIdPool {
	/* Local pool range is implemented by subdividing the total range of IDs into
	 * numProcesses blocks of poolSize / numProcesses size assigning them consecutively
	 * to the processes. */
public:
	MoleculeIdPool(unsigned long poolsize, int numProcs, int myProcRank) :
		_poolSize(poolsize), _numProcesses(numProcs), _myProcRank(myProcRank), _moleculesFromThisProcess(0) {}

	/** get a new, jet unused molecule ID */
	unsigned long getNewMoleculeId() {
		mardyn_assert(_moleculesFromThisProcess < localIdRangeSize());
		return myIDoffset() + _moleculesFromThisProcess++;
	}
	/** get the rank of the process belonging to molecule ID */
	int getOwnerRank(unsigned long id) {
		return id / localIdRangeSize();
	}
	/** get the number of molecule IDs this process has been assigned */
	unsigned long localIdRangeSize(){
		static unsigned long localIdRangeSize = _poolSize / _numProcesses;
		return localIdRangeSize;
	}

private:
	unsigned long myIDoffset() {
		static unsigned long myIDoffset = localIdRangeSize() * _myProcRank;
		return myIDoffset;
	}

	unsigned long _poolSize;  //!< overall pool size
	int _numProcesses;  //!< number of process sharing this pool
	int _myProcRank;  //!< the rank of the local process
	unsigned long _moleculesFromThisProcess;  //!< number of molecules IDs already provided
};

#endif  // SRC_MOLECULES_MOLECULEIDPOOL_H_

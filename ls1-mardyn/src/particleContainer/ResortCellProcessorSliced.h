/*
 * ResortCellProcessorSliced.h
 *
 *  Created on: 10 Oct 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_RESORTCELLPROCESSORSLICED_H_
#define SRC_PARTICLECONTAINER_RESORTCELLPROCESSORSLICED_H_

#include "adapter/CellProcessor.h"
#include "LinkedCells.h"

// special update cell processor for the REDUCED_MEMORY_MODE
// to be used only when the SLICED traversal is being used.

class ResortCellProcessorSliced: public CellProcessor {
public:
	ResortCellProcessorSliced(LinkedCells * container) :
			CellProcessor(0.0, 0.0), _container(container) {
		// allocate threadData
		_threadData.resize(mardyn_get_max_threads());
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			ThreadData * myown = new ThreadData();
			const int myid = mardyn_get_thread_num();
			_threadData[myid] = myown;
		} // end pragma omp parallel

	}
	~ResortCellProcessorSliced() {
		// free threadData
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			const int myid = mardyn_get_thread_num();
			delete _threadData[myid];
		}
	};
	void initTraversal() {}
	void preprocessCell(ParticleCell& ) {}
	void processCellPair(ParticleCell&, ParticleCell&) {}
	void processCell(ParticleCell& cell) {
		// get leaving molecules
		std::vector<Molecule>& b = _threadData[mardyn_get_thread_num()]->_buffer;

		cell.getLeavingMolecules(b);

		// add leaving molecules
		for (auto & mol : b) {
			_container->addParticle(mol);
		}
		b.clear();
	}
	double processSingleMolecule(Molecule*, ParticleCell& ) { return 0.0;}
	void postprocessCell(ParticleCell& ) {}
	void endTraversal() {}
private:
	class ThreadData {
	public:
		ThreadData() {
			_buffer.reserve(2);
		}
		std::vector<Molecule> _buffer;
	};

	std::vector<ThreadData * > _threadData;

	LinkedCells * _container;

};

#endif /* SRC_PARTICLECONTAINER_RESORTCELLPROCESSORSLICED_H_ */

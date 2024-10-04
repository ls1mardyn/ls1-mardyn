#ifndef TRAVERSALTUNER_H_
#define TRAVERSALTUNER_H_

#include <algorithm>
#include <utility>
#include <vector>

#include <utils/Logger.h>
#include <Simulation.h>
#include "LinkedCellTraversals/CellPairTraversals.h"
#include "LinkedCellTraversals/QuickschedTraversal.h"
#include "LinkedCellTraversals/C08CellPairTraversal.h"
#include "LinkedCellTraversals/C04CellPairTraversal.h"
#include "LinkedCellTraversals/OriginalCellPairTraversal.h"
#include "LinkedCellTraversals/HalfShellTraversal.h"
#include "LinkedCellTraversals/MidpointTraversal.h"
#include "LinkedCellTraversals/NeutralTerritoryTraversal.h"
#include "LinkedCellTraversals/SlicedCellPairTraversal.h"


template<class CellTemplate>
class TraversalTuner {
	friend class LinkedCellsTest;

public:
	// Probably remove this once autotuning is implemented
	enum traversalNames {
		ORIGINAL = 0,
		C08      = 1,
		C04      = 2,
		SLICED   = 3,
		HS       = 4,
		MP       = 5,
		C08ES    = 6,
		NT       = 7,
		// quicksched has to be the last traversal!
		QSCHED   = 8,
	};

	TraversalTuner();

	~TraversalTuner();

	void findOptimalTraversal();

	void readXML(XMLfileUnits &xmlconfig);

	/**
	 * Rebuild the traversals.
	 * @param cells The vector of cells.
	 * @param dims The dimensions (cells per dimension, including halo!)
	 * @param cellLength The length of the cells.
	 * @param cutoff The cutoff radius.
	 */
	void rebuild(std::vector<CellTemplate> &cells,
				 const std::array<unsigned long, 3> &dims, double cellLength[3], double cutoff);

	void traverseCellPairs(CellProcessor &cellProcessor);

	void traverseCellPairs(traversalNames name, CellProcessor &cellProcessor);

	void traverseCellPairsOuter(CellProcessor &cellProcessor);

	void traverseCellPairsInner(CellProcessor &cellProcessor, unsigned stage, unsigned stageCount);


	bool isTraversalApplicable(traversalNames name, const std::array<unsigned long, 3> &dims) const; // new


	traversalNames getSelectedTraversal() const {
		return selectedTraversal;
	}

	CellPairTraversals<ParticleCell> *getCurrentOptimalTraversal() { return _optimalTraversal; }

private:
	std::vector<CellTemplate>* _cells;
	std::array<unsigned long, 3> _dims;

	traversalNames selectedTraversal;

	std::vector<std::pair<CellPairTraversals<CellTemplate> *, CellPairTraversalData *> > _traversals;

	CellPairTraversals<CellTemplate> *_optimalTraversal;

	unsigned _cellsInCutoff = 1;
};

template<class CellTemplate>
TraversalTuner<CellTemplate>::TraversalTuner() : _cells(nullptr), _dims(), _optimalTraversal(nullptr) {
	// defaults:
	selectedTraversal = {
			mardyn_get_max_threads() > 1 ? C08 : SLICED
	};
	auto *c08Data = new C08CellPairTraversalData;
	auto *c04Data = new C04CellPairTraversalData;
	auto *origData = new OriginalCellPairTraversalData;
	auto *slicedData = new SlicedCellPairTraversalData;
	auto *hsData = new HalfShellTraversalData;
	auto *mpData = new MidpointTraversalData;
	auto *ntData = new NeutralTerritoryTraversalData;
	auto *c08esData = new C08CellPairTraversalData;

	_traversals = {
			std::make_pair(nullptr, origData),
			std::make_pair(nullptr, c08Data),
			std::make_pair(nullptr, c04Data),
			std::make_pair(nullptr, slicedData),
			std::make_pair(nullptr, hsData),
			std::make_pair(nullptr, mpData),
			std::make_pair(nullptr, ntData),
			std::make_pair(nullptr, c08esData)
	};
#ifdef QUICKSCHED
	struct QuickschedTraversalData *quiData = new QuickschedTraversalData;
	quiData->taskBlockSize = {{2, 2, 2}};
	if (std::is_base_of<ParticleCellBase, CellTemplate>::value) {
		_traversals.push_back(std::make_pair(nullptr, quiData));
	}
#endif
}

template<class CellTemplate>
TraversalTuner<CellTemplate>::~TraversalTuner() {
	for (auto t : _traversals) {
		if (t.first != nullptr)
			delete (t.first);
		if (t.second != nullptr)
			delete (t.second);
	}
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::findOptimalTraversal() {
	// TODO implement autotuning here! At the moment the traversal is chosen via readXML!

	_optimalTraversal = _traversals[selectedTraversal].first;

	// log traversal
	if (dynamic_cast<HalfShellTraversal<CellTemplate> *>(_optimalTraversal))
		Log::global_log->info() << "Using HalfShellTraversal." << std::endl;
	else if (dynamic_cast<OriginalCellPairTraversal<CellTemplate> *>(_optimalTraversal))
		Log::global_log->info() << "Using OriginalCellPairTraversal." << std::endl;
	else if (dynamic_cast<C08CellPairTraversal<CellTemplate> *>(_optimalTraversal))
		Log::global_log->info() << "Using C08CellPairTraversal without eighthShell." << std::endl;
	else if (dynamic_cast<C08CellPairTraversal<CellTemplate, true> *>(_optimalTraversal))
		Log::global_log->info() << "Using C08CellPairTraversal with eighthShell." << std::endl;
	else if (dynamic_cast<C04CellPairTraversal<CellTemplate> *>(_optimalTraversal))
		Log::global_log->info() << "Using C04CellPairTraversal." << std::endl;
	else if (dynamic_cast<MidpointTraversal<CellTemplate> *>(_optimalTraversal))
		Log::global_log->info() << "Using MidpointTraversal." << std::endl;
	else if (dynamic_cast<NeutralTerritoryTraversal<CellTemplate> *>(_optimalTraversal))
		Log::global_log->info() << "Using NeutralTerritoryTraversal." << std::endl;
	else if (dynamic_cast<QuickschedTraversal<CellTemplate> *>(_optimalTraversal)) {
		Log::global_log->info() << "Using QuickschedTraversal." << std::endl;
#ifndef QUICKSCHED
		std::ostringstream error_message;
		error_message << "MarDyn was compiled without Quicksched Support. Aborting!" << std::endl;
		MARDYN_EXIT(error_message);
#endif
	} else if (dynamic_cast<SlicedCellPairTraversal<CellTemplate> *>(_optimalTraversal))
		Log::global_log->info() << "Using SlicedCellPairTraversal." << std::endl;
	else
		Log::global_log->warning() << "Using unknown traversal." << std::endl;

	if (_cellsInCutoff > _optimalTraversal->maxCellsInCutoff()) {
		std::ostringstream error_message;		error_message << "Traversal supports up to " << _optimalTraversal->maxCellsInCutoff()
							<< " cells in cutoff, but value is chosen as " << _cellsInCutoff << std::endl;		MARDYN_EXIT(error_message);
	}
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::readXML(XMLfileUnits &xmlconfig) {
	std::string oldPath(xmlconfig.getcurrentnodepath());
	// read traversal type default values
	std::string traversalType;

	xmlconfig.getNodeValue("traversalSelector", traversalType);
	transform(traversalType.begin(), traversalType.end(), traversalType.begin(), ::tolower);

	if (traversalType.find("c08es") != std::string::npos)
		selectedTraversal = C08ES;
	else if (traversalType.find("c08") != std::string::npos)
		selectedTraversal = C08;
	else if (traversalType.find("c04") != std::string::npos)
		selectedTraversal = C04;
	else if (traversalType.find("qui") != std::string::npos)
		selectedTraversal = QSCHED;
	else if (traversalType.find("slice") != std::string::npos)
		selectedTraversal = SLICED;
	else if (traversalType.find("ori") != std::string::npos)
		selectedTraversal = ORIGINAL;
	else if (traversalType.find("hs") != std::string::npos)
		selectedTraversal = HS;
	else if (traversalType.find("mp") != std::string::npos)
		selectedTraversal = MP;
	else if (traversalType.find("nt") != std::string::npos) {
		selectedTraversal = NT;
	} else {
		// selector already set in constructor, just print a warning here
		if (mardyn_get_max_threads() > 1) {
			Log::global_log->warning() << "No traversal type selected. Defaulting to c08 traversal." << std::endl;
		} else {
			Log::global_log->warning() << "No traversal type selected. Defaulting to sliced traversal." << std::endl;
		}
	}

	_cellsInCutoff = xmlconfig.getNodeValue_int("cellsInCutoffRadius", 1); // This is currently only used for an assert

	// workaround for stupid iterator:
	// since
	// xmlconfig.changecurrentnode(traversalIterator);
	// does not work resolve paths to traversals manually
	// use iterator only to resolve number of traversals (==iterations)
	std::string basePath(xmlconfig.getcurrentnodepath());

	int i = 1;
	XMLfile::Query qry = xmlconfig.query("traversalData");
	for (XMLfile::Query::const_iterator traversalIterator = qry.begin(); traversalIterator; ++traversalIterator) {
		std::string path(basePath + "/traversalData[" + std::to_string(i) + "]");
		xmlconfig.changecurrentnode(path);

		traversalType = xmlconfig.getNodeValue_string("@type", "NOTHING FOUND");
		transform(traversalType.begin(), traversalType.end(), traversalType.begin(), ::tolower);
		if (traversalType == "c08") {
			// nothing to do
		} else if (traversalType.find("qui") != std::string::npos) {
#ifdef QUICKSCHED
			if (not std::is_base_of<ParticleCellBase, CellTemplate>::value) {
				Log::global_log->warning() << "Attempting to use Quicksched with cell type that does not store task data!"
									  << std::endl;
			}
			for (auto p : _traversals) {
				if (struct QuickschedTraversalData *quiData = dynamic_cast<QuickschedTraversalData *>(p.second)) {
					// read task block size
					std::string tag       = "taskBlockSize/l";
					char   dimension = 'x';

					for (int j = 0; j < 3; ++j) {
						tag += (dimension + j);
						xmlconfig.getNodeValue(tag, quiData->taskBlockSize[j]);
						if (quiData->taskBlockSize[j] < 2) {
							std::ostringstream error_message;							error_message << "Task block size in "
												<< (char) (dimension + j)
												<< " direction is <2 and thereby invalid! ("
												<< quiData->taskBlockSize[j] << ")"
												<< std::endl;							MARDYN_EXIT(error_message);
						}
					}
					break;
				}
			}
#else
			Log::global_log->warning() << "Found quicksched traversal data in config "
								  << "but mardyn was compiled without quicksched support! "
								  << "(make ENABLE_QUICKSCHED=1)" << std::endl;
#endif
		} else {
			Log::global_log->warning() << "Unknown traversal type: " << traversalType << std::endl;
		}
		++i;
	}
	xmlconfig.changecurrentnode(oldPath);
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::rebuild(std::vector<CellTemplate> &cells, const std::array<unsigned long, 3> &dims,
										   double cellLength[3], double cutoff) {
	_cells = &cells; // new - what for?
	_dims = dims; // new - what for?

	for (size_t i = 0ul; i < _traversals.size(); ++i) {
		auto& [traversalPointerReference, traversalData] = _traversals[i];
		// decide whether to initialize or rebuild
		if (traversalPointerReference == nullptr) {
			switch (i) {
				case traversalNames ::ORIGINAL:
					traversalPointerReference = new OriginalCellPairTraversal<CellTemplate>(cells, dims);
					break;
				case traversalNames::C08:
					traversalPointerReference = new C08CellPairTraversal<CellTemplate>(cells, dims);
					break;
				case traversalNames::C04:
					traversalPointerReference = new C04CellPairTraversal<CellTemplate>(cells, dims);
					break;
				case traversalNames::SLICED:
					traversalPointerReference = new SlicedCellPairTraversal<CellTemplate>(cells, dims);
					break;
				case traversalNames::HS:
					traversalPointerReference = new HalfShellTraversal<CellTemplate>(cells, dims);
					break;
				case traversalNames::MP:
					traversalPointerReference = new MidpointTraversal<CellTemplate>(cells, dims);
					break;
				case traversalNames::NT:
					traversalPointerReference =
						new NeutralTerritoryTraversal<CellTemplate>(cells, dims, cellLength, cutoff);
					break;
				case traversalNames::C08ES:
					traversalPointerReference = new C08CellPairTraversal<CellTemplate, true>(cells, dims);
					break;
				case traversalNames::QSCHED: {
					mardyn_assert((std::is_base_of<ParticleCellBase, CellTemplate>::value));
					auto *quiData = dynamic_cast<QuickschedTraversalData *>(traversalData);
					traversalPointerReference = new QuickschedTraversal<CellTemplate>(cells, dims, quiData->taskBlockSize);
				} break;
				default:
					std::ostringstream error_message;
					error_message << "Unknown traversal data found in TraversalTuner._traversals!" << std::endl;
					MARDYN_EXIT(error_message);
			}
		}
		traversalPointerReference->rebuild(cells, dims, cellLength, cutoff, traversalData);
	}
	_optimalTraversal = nullptr;
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::traverseCellPairs(CellProcessor &cellProcessor) {
	if (not _optimalTraversal) {
		findOptimalTraversal();
	}
	_optimalTraversal->traverseCellPairs(cellProcessor);
}

template<class CellTemplate>
inline void TraversalTuner<CellTemplate>::traverseCellPairs(traversalNames name,
		CellProcessor& cellProcessor) {
	if (name == getSelectedTraversal()) {
		traverseCellPairs(cellProcessor);
	} else {
		SlicedCellPairTraversal<CellTemplate> slicedTraversal(*_cells, _dims);
		switch(name) {
		case SLICED:
			slicedTraversal.traverseCellPairs(cellProcessor);
			break;
		default:
			std::ostringstream error_message;
			error_message<< "Calling traverseCellPairs(traversalName, CellProcessor&) for something else than the Sliced Traversal is disabled for now. Aborting." << std::endl;
			MARDYN_EXIT(error_message);
			break;
		}
	}
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::traverseCellPairsOuter(CellProcessor &cellProcessor) {
	if (not _optimalTraversal) {
		findOptimalTraversal();
	}
	_optimalTraversal->traverseCellPairsOuter(cellProcessor);
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::traverseCellPairsInner(CellProcessor &cellProcessor, unsigned stage,
														  unsigned stageCount) {
	if (not _optimalTraversal) {
		findOptimalTraversal();
	}
	_optimalTraversal->traverseCellPairsInner(cellProcessor, stage, stageCount);
}

template<class CellTemplate>
inline bool TraversalTuner<CellTemplate>::isTraversalApplicable(
		traversalNames name, const std::array<unsigned long, 3> &dims) const {
	bool ret = true;
	switch(name) {
	case SLICED:
		ret = SlicedCellPairTraversal<CellTemplate>::isApplicable(dims);
		break;
	case QSCHED:
#ifdef QUICKSCHED
		ret = true;
#else
		ret = false;
#endif
		break;
	case C08:
		ret = true;
		break;
	case C04:
		ret = true;
		break;
	case ORIGINAL:
		ret = true;
		break;
	default:
		Log::global_log->warning() << "unknown traversal given in TraversalTuner::isTraversalApplicable, assuming that is applicable" << std::endl;
	}
	return ret;
}

#endif //TRAVERSALTUNER_H_

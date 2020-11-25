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
#include "LinkedCellTraversals/SlicedCellPairTraversal.h"

using Log::global_log;

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
		QSCHED   = 7,
	};

	TraversalTuner();

	~TraversalTuner();

	void findOptimalTraversal();

	void readXML(XMLfileUnits &xmlconfig);

	void rebuild(std::vector<CellTemplate> &cells,
				 const std::array<unsigned long, 3> &dims);

	void traverseCellPairs(CellProcessor &cellProcessor);

	void traverseCellPairs(traversalNames name, CellProcessor &cellProcessor);

	void traverseCellPairsOuter(CellProcessor &cellProcessor);

	void traverseCellPairsInner(CellProcessor &cellProcessor, unsigned stage, unsigned stageCount);


	bool isTraversalApplicable(traversalNames name, const std::array<unsigned long, 3> &dims) const; // new


	traversalNames getSelectedTraversal() const {
		return selectedTraversal;
	}
        
    CellPairTraversals<ParticleCell>* getCurrentOptimalTraversal() {
            return _optimalTraversal;
    }
        
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
	struct C08CellPairTraversalData      *c08Data    = new C08CellPairTraversalData;
	struct C04CellPairTraversalData      *c04Data    = new C04CellPairTraversalData;
	struct OriginalCellPairTraversalData *origData   = new OriginalCellPairTraversalData;
	struct SlicedCellPairTraversalData   *slicedData = new SlicedCellPairTraversalData;
    struct HalfShellTraversalData    	 *hsData   	 = new HalfShellTraversalData;
    struct MidpointTraversalData    	 *mpData   	 = new MidpointTraversalData;
	struct C08CellPairTraversalData    	 *c08esData  = new C08CellPairTraversalData;


	_traversals = {
			make_pair(nullptr, origData),
			make_pair(nullptr, c08Data),
			make_pair(nullptr, c04Data),
			make_pair(nullptr, slicedData),
            make_pair(nullptr, hsData),
            make_pair(nullptr, mpData),
			make_pair(nullptr, c08esData)
	};
#ifdef QUICKSCHED
	struct QuickschedTraversalData *quiData = new QuickschedTraversalData;
	quiData->taskBlockSize = {{2, 2, 2}};
	if (is_base_of<ParticleCellBase, CellTemplate>::value) {
		_traversals.push_back(make_pair(nullptr, quiData));
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
		global_log->info() << "Using HalfShellTraversal." << endl;
	else if (dynamic_cast<OriginalCellPairTraversal<CellTemplate> *>(_optimalTraversal))
		global_log->info() << "Using OriginalCellPairTraversal." << endl;
	else if (dynamic_cast<C08CellPairTraversal<CellTemplate> *>(_optimalTraversal))
		global_log->info() << "Using C08CellPairTraversal without eighthShell." << endl;
	else if (dynamic_cast<C08CellPairTraversal<CellTemplate, true> *>(_optimalTraversal))
		global_log->info() << "Using C08CellPairTraversal with eighthShell." << endl;
	else if (dynamic_cast<C04CellPairTraversal<CellTemplate> *>(_optimalTraversal))
		global_log->info() << "Using C04CellPairTraversal." << endl;
	else if (dynamic_cast<MidpointTraversal<CellTemplate> *>(_optimalTraversal))
		global_log->info() << "Using MidpointTraversal." << endl;
	else if (dynamic_cast<QuickschedTraversal<CellTemplate> *>(_optimalTraversal)) {
		global_log->info() << "Using QuickschedTraversal." << endl;
#ifndef QUICKSCHED
		global_log->error() << "MarDyn was compiled without Quicksched Support. Aborting!" << endl;
		Simulation::exit(1);
#endif
	} else if (dynamic_cast<SlicedCellPairTraversal<CellTemplate> *>(_optimalTraversal))
		global_log->info() << "Using SlicedCellPairTraversal." << endl;
	else
		global_log->warning() << "Using unknown traversal." << endl;

	if (_optimalTraversal->maxCellsInCutoff() < _cellsInCutoff) {
		global_log->error() << "Traversal supports up to " << _optimalTraversal->maxCellsInCutoff()
							<< " cells in cutoff, but value is chosen as " << _cellsInCutoff << std::endl;
		Simulation::exit(45);
	}
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::readXML(XMLfileUnits &xmlconfig) {
	string oldPath(xmlconfig.getcurrentnodepath());
	// read traversal type default values
	string traversalType;

	xmlconfig.getNodeValue("traversalSelector", traversalType);
	transform(traversalType.begin(),
			  traversalType.end(),
			  traversalType.begin(),
			  ::tolower);

	if (traversalType.find("c08es") != string::npos)
		selectedTraversal = C08ES;
	else if (traversalType.find("c08") != string::npos)
		selectedTraversal = C08;
	else if (traversalType.find("c04") != string::npos)
		selectedTraversal = C04;
	else if (traversalType.find("qui") != string::npos)
		selectedTraversal = QSCHED;
	else if (traversalType.find("slice") != string::npos)
		selectedTraversal = SLICED;
	else if (traversalType.find("ori") != string::npos)
		selectedTraversal = ORIGINAL;
	else if (traversalType.find("hs") != string::npos)
		selectedTraversal = HS;
	else if (traversalType.find("mp") != string::npos)
		selectedTraversal = MP;
	else if (traversalType.find("nt") != string::npos) {
		global_log->error() << "nt method not yet properly implemented. please select a different method." << std::endl;
		Simulation::exit(1);
	} else {
		// selector already set in constructor, just print a warning here
		if (mardyn_get_max_threads() > 1) {
			global_log->warning() << "No traversal type selected. Defaulting to c08 traversal." << endl;
		} else {
			global_log->warning() << "No traversal type selected. Defaulting to sliced traversal." << endl;
		}
	}

	_cellsInCutoff = xmlconfig.getNodeValue_int("cellsInCutoffRadius", 1); // This is currently only used for an assert

	// workaround for stupid iterator:
	// since
	// xmlconfig.changecurrentnode(traversalIterator);
	// does not work resolve paths to traversals manually
	// use iterator only to resolve number of traversals (==iterations)
	string basePath(xmlconfig.getcurrentnodepath());

	int i = 1;
	XMLfile::Query qry = xmlconfig.query("traversalData");
	for (XMLfile::Query::const_iterator traversalIterator = qry.begin(); traversalIterator; ++traversalIterator) {
		string path(basePath + "/traversalData[" + to_string(i) + "]");
		xmlconfig.changecurrentnode(path);

		traversalType = xmlconfig.getNodeValue_string("@type", "NOTHING FOUND");
		transform(traversalType.begin(),
				  traversalType.end(),
				  traversalType.begin(),
				  ::tolower);
		if (traversalType == "c08") {
			// nothing to do
		} else if (traversalType.find("qui") != string::npos) {
#ifdef QUICKSCHED
			if (not is_base_of<ParticleCellBase, CellTemplate>::value) {
				global_log->warning() << "Attempting to use Quicksched with cell type that does not store task data!"
									  << endl;
			}
			for (auto p : _traversals) {
				if (struct QuickschedTraversalData *quiData = dynamic_cast<QuickschedTraversalData *>(p.second)) {
					// read task block size
					string tag       = "taskBlockSize/l";
					char   dimension = 'x';

					for (int j = 0; j < 3; ++j) {
						tag += (dimension + j);
						xmlconfig.getNodeValue(tag, quiData->taskBlockSize[j]);
						if (quiData->taskBlockSize[j] < 2) {
							global_log->error() << "Task block size in "
												<< (char) (dimension + j)
												<< " direction is <2 and thereby invalid! ("
												<< quiData->taskBlockSize[j] << ")"
												<< endl;
							Simulation::exit(1);
						}
					}
					break;
				}
			}
#else
			global_log->warning() << "Found quicksched traversal data in config "
								  << "but mardyn was compiled without quicksched support! "
								  << "(make ENABLE_QUICKSCHED=1)"
								  << endl;
#endif
		} else {
			global_log->warning() << "Unknown traversal type: " << traversalType << endl;
		}
		++i;
	}
	xmlconfig.changecurrentnode(oldPath);
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::rebuild(std::vector<CellTemplate> &cells, 
        const std::array<unsigned long, 3> &dims) {
    
	_cells = &cells; // new - what for?
	_dims = dims; // new - what for?

	for (size_t i = 0ul; i < _traversals.size(); ++i) {
		auto& tPair = _traversals[i];
		// decide whether to initialize or rebuild
		if (tPair.first == nullptr) {
			switch (i) {
				case traversalNames ::ORIGINAL:
					tPair.first = new OriginalCellPairTraversal<CellTemplate>(cells, dims);
					break;
				case traversalNames::C08:
					tPair.first = new C08CellPairTraversal<CellTemplate>(cells, dims);
					break;
				case traversalNames::C04:
					tPair.first = new C04CellPairTraversal<CellTemplate>(cells, dims);
					break;
				case traversalNames::SLICED:
					tPair.first = new SlicedCellPairTraversal<CellTemplate>(cells, dims);
					break;
				case traversalNames::HS:
					tPair.first = new HalfShellTraversal<CellTemplate>(cells, dims);
					break;
				case traversalNames::MP:
					tPair.first = new MidpointTraversal<CellTemplate>(cells, dims);
					break;
				case traversalNames::C08ES:
					tPair.first = new C08CellPairTraversal<CellTemplate, true>(cells, dims);
					break;
				case traversalNames::QSCHED: {
					mardyn_assert((is_base_of<ParticleCellBase, CellTemplate>::value));
					QuickschedTraversalData *quiData = dynamic_cast<QuickschedTraversalData *>(tPair.second);
					tPair.first = new QuickschedTraversal<CellTemplate>(cells, dims, quiData->taskBlockSize);
				} break;
				default:
					global_log->error() << "Unknown traversal data found in TraversalTuner._traversals!" << endl;
					Simulation::exit(1);
			}
		}
		tPair.first->rebuild(cells, dims, tPair.second);
	}
	_optimalTraversal = nullptr;
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::traverseCellPairs(CellProcessor &cellProcessor) {
	if (_optimalTraversal == nullptr)
		findOptimalTraversal();
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
			Log::global_log->error()<< "Calling traverseCellPairs(traversalName, CellProcessor&) for something else than the Sliced Traversal is disabled for now. Aborting." << std::endl;
			mardyn_exit(1);
			break;
		}
	}
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::traverseCellPairsOuter(CellProcessor &cellProcessor) {
	if (_optimalTraversal == nullptr)
		findOptimalTraversal();
	_optimalTraversal->traverseCellPairsOuter(cellProcessor);
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::traverseCellPairsInner(CellProcessor &cellProcessor, unsigned stage,
														  unsigned stageCount) {
	if (_optimalTraversal == nullptr)
		findOptimalTraversal();
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
		global_log->warning() << "unknown traversal given in TraversalTuner::isTraversalApplicable, assuming that is applicable" << std::endl;
	}
	return ret;
}

#endif //TRAVERSALTUNER_H_

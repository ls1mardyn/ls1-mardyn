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
#include "LinkedCellTraversals/OriginalCellPairTraversal.h"
#include "LinkedCellTraversals/HalfShellTraversal.h"
#include "LinkedCellTraversals/MidpointTraversal.h"
#include "LinkedCellTraversals/SlicedCellPairTraversal.h"

using Log::global_log;

template<class CellTemplate>
class TraversalTuner {
	friend class LinkedCellsTest;

public:
	TraversalTuner();

	~TraversalTuner();

	void findOptimalTraversal();

	void readXML(XMLfileUnits &xmlconfig);

	void rebuild(std::vector<CellTemplate> &cells,
				 const std::array<unsigned long, 3> &dims);

	void traverseCellPairs(CellProcessor &cellProcessor);

	void traverseCellPairsOuter(CellProcessor &cellProcessor);

	void traverseCellPairsInner(CellProcessor &cellProcessor, unsigned stage, unsigned stageCount);

    CellPairTraversals<ParticleCell>* getSelectedTraversal(){
    	return _optimalTravesal;
    }


private:

    // Probably remove this once autotuning is implemented
    enum traversalNames {
        ORIGINAL = 0,
        C08      = 1,
        SLICED   = 2,
		HS	 	 = 3,
		MP	 	 = 4,
        QSCHED   = 5
    };

	traversalNames selectedTraversal;

	std::vector<std::pair<CellPairTraversals<CellTemplate> *, CellPairTraversalData *> > _traversals;

	CellPairTraversals<CellTemplate> *_optimalTravesal;
};

template<class CellTemplate>
TraversalTuner<CellTemplate>::TraversalTuner() : _optimalTravesal(nullptr) {
	// defaults:
	selectedTraversal = {
			mardyn_get_max_threads() > 1 ? C08 : SLICED
	};
	struct C08CellPairTraversalData      *c08Data    = new C08CellPairTraversalData;
	struct OriginalCellPairTraversalData *origData   = new OriginalCellPairTraversalData;
	struct SlicedCellPairTraversalData   *slicedData = new SlicedCellPairTraversalData;
    struct HalfShellTraversalData    	 *hsData   	 = new HalfShellTraversalData;
    struct MidpointTraversalData    	 *mpData   	 = new MidpointTraversalData;


	_traversals = {
			make_pair(nullptr, origData),
			make_pair(nullptr, c08Data),
			make_pair(nullptr, slicedData),
            make_pair(nullptr, hsData),
            make_pair(nullptr, mpData)
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

	_optimalTravesal = _traversals[selectedTraversal].first;

	// log traversal
    if (dynamic_cast<HalfShellTraversal<CellTemplate> *>(_optimalTravesal))
           global_log->info() << "Using HalfShellTraversal." << endl;
    else if (dynamic_cast<OriginalCellPairTraversal<CellTemplate> *>(_optimalTravesal))
        global_log->info() << "Using OriginalCellPairTraversal." << endl;
    else if (dynamic_cast<C08CellPairTraversal<CellTemplate> *>(_optimalTravesal))
        global_log->info() << "Using C08CellPairTraversal." << endl;
    else if (dynamic_cast<MidpointTraversal<CellTemplate> *>(_optimalTravesal))
        global_log->info() << "Using MidpointTraversal." << endl;

	else if (dynamic_cast<QuickschedTraversal<CellTemplate> *>(_optimalTravesal)) {
		global_log->info() << "Using QuickschedTraversal." << endl;
#ifndef QUICKSCHED
		global_log->error() << "MarDyn was compiled without Quicksched Support. Aborting!" << endl;
		Simulation::exit(1);
#endif
	} else if (dynamic_cast<SlicedCellPairTraversal<CellTemplate> *>(_optimalTravesal))
		global_log->info() << "Using SlicedCellPairTraversal." << endl;
	else
		global_log->warning() << "Using unknown traversal." << endl;

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

	if (traversalType.find("c08") != string::npos)
		selectedTraversal = C08;
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
	else {
		// selector already set in constructor, just print a warning here
		if (mardyn_get_max_threads() > 1) {
			global_log->warning() << "No traversal type selected. Defaulting to c08 traversal." << endl;
		} else {
			global_log->warning() << "No traversal type selected. Defaulting to sliced traversal." << endl;
		}
	}

	// workaround for stupid iterator:
	// since
	// xmlconfig.changecurrentnode(traversalIterator);
	// does not work resolve paths to traversals manually
	// use iterator only to resolve number of traversals (==iterations)
	string basePath(xmlconfig.getcurrentnodepath());

	int i = 1;

	for (XMLfile::Query::const_iterator traversalIterator = xmlconfig.query("traversalData").begin();
		 traversalIterator;
		 ++traversalIterator) {
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
	for (auto &tPair : _traversals) {
		// decide whether to initialize or rebuild
		if (tPair.first == nullptr) {
			if (dynamic_cast<C08CellPairTraversalData *>(tPair.second)) {
				tPair.first = new C08CellPairTraversal<CellTemplate>(cells,
																	 dims);
			} else if (QuickschedTraversalData *quiData = dynamic_cast<QuickschedTraversalData *>(tPair.second)) {
				mardyn_assert((is_base_of<ParticleCellBase, CellTemplate>::value));
				tPair.first = new QuickschedTraversal<CellTemplate>(cells,
																	dims,
																	quiData->taskBlockSize);
			} else if (dynamic_cast<HalfShellTraversalData *>(tPair.second)) {
			    tPair.first = new HalfShellTraversal<CellTemplate>(cells, dims);
		    } else if (dynamic_cast<MidpointTraversalData *>(tPair.second)) {
			    tPair.first = new MidpointTraversal<CellTemplate>(cells, dims);
			} else if (dynamic_cast<OriginalCellPairTraversalData *>(tPair.second)) {
				tPair.first = new OriginalCellPairTraversal<CellTemplate>(cells,
																		  dims);
			} else if (dynamic_cast<SlicedCellPairTraversalData *>(tPair.second)) {
				tPair.first = new SlicedCellPairTraversal<CellTemplate>(cells,
																		dims);
			} else {
				global_log->error() << "Unknown traversal data found in TraversalTuner._traversals!" << endl;
				Simulation::exit(1);
			}
		}
		tPair.first->rebuild(cells, dims, tPair.second);
	}
	_optimalTravesal = nullptr;
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::traverseCellPairs(CellProcessor &cellProcessor) {
	if (_optimalTravesal == nullptr)
		findOptimalTraversal();
	_optimalTravesal->traverseCellPairs(cellProcessor);
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::traverseCellPairsOuter(CellProcessor &cellProcessor) {
	if (_optimalTravesal == nullptr)
		findOptimalTraversal();
	_optimalTravesal->traverseCellPairsOuter(cellProcessor);
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::traverseCellPairsInner(CellProcessor &cellProcessor, unsigned stage,
														  unsigned stageCount) {
	if (_optimalTravesal == nullptr)
		findOptimalTraversal();
	_optimalTravesal->traverseCellPairsInner(cellProcessor, stage, stageCount);
}


#endif //TRAVERSALTUNER_H_

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

using Log::global_log;

class TraversalTuner {
	friend class LinkedCellsTest;

public:
    TraversalTuner();

    ~TraversalTuner();

    void findOptimalTraversal();

    void readXML(XMLfileUnits &xmlconfig);

    void rebuild(std::vector<ParticleCell> &cells,
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
        QSCHED   = 3,
		HS		 = 4
    };

    traversalNames selectedTraversal;

    std::vector<std::pair<CellPairTraversals<ParticleCell> *, CellPairTraversalData *> > _traversals;

    CellPairTraversals<ParticleCell> *_optimalTravesal;
};


#endif //TRAVERSALTUNER_H_

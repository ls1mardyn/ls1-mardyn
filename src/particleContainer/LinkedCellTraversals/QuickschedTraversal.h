/*
 * QuickschedTraversal.h
 *
 *  Created on: 20 May 2017
 *      Author: gratlf
 */

#ifndef SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_QUICKSCHEDTRAVERSAL_H
#define SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_QUICKSCHEDTRAVERSAL_H

#include "C08BasedTraversals.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h"
#include "particleContainer/LinkedCells.h"

#ifdef QUICKSCHED
#include "quicksched.h"
#endif


struct QuickschedTraversalData : CellPairTraversalData {
    std::array<unsigned long, 3> taskBlockSize;
};

template<class CellTemplate>
class QuickschedTraversal : public C08BasedTraversals<CellTemplate> {
public:
    QuickschedTraversal(std::vector<CellTemplate> &cells,
                        const std::array<unsigned long, 3> &dims,
                        const std::array<unsigned long, 3> &taskBlockSize);

    ~QuickschedTraversal();


    virtual void rebuild(std::vector<CellTemplate> &cells,
                         const std::array<unsigned long, 3> &dims, double cellLength[3], double cutoff,
                         CellPairTraversalData *data);

    void traverseCellPairs(CellProcessor &cellProcessor);

    void traverseCellPairsOuter(CellProcessor &cellProcessor);

    void traverseCellPairsInner(CellProcessor &cellProcessor, unsigned stage, unsigned stageCount);

private:
    enum taskType {
        PackedAdjustable
    };

    void init();

    static void runner(int type, void *data);

    std::array<unsigned long, 3> _taskBlocksize;
    CellProcessor           *_contextCellProcessor;

    struct qsched *_scheduler;
    taskType      _taskTypeSelector;
};

template<class CellTemplate>
QuickschedTraversal<CellTemplate>::QuickschedTraversal(std::vector<CellTemplate> &cells,
                                                       const std::array<unsigned long, 3> &dims,
                                                       const std::array<unsigned long, 3> &taskBlockSize)
#ifdef QUICKSCHED
        : C08BasedTraversals<CellTemplate>(cells, dims),
          _taskBlocksize(taskBlockSize),
          _contextCellProcessor(nullptr),
          _scheduler(new struct qsched),
          _taskTypeSelector(PackedAdjustable) {
    mardyn_assert((std::is_base_of<ParticleCellBase, CellTemplate>::value));
    qsched_init(_scheduler, mardyn_get_max_threads(), qsched_flag_none);
}

#else
        : C08BasedTraversals<CellTemplate>(cells, dims) {}
#endif /* QUICKSCHED */

template<class CellTemplate>
void QuickschedTraversal<CellTemplate>::rebuild(std::vector<CellTemplate> &cells,
                                                const std::array<unsigned long, 3> &dims, double cellLength[3], double cutoff,
                                                CellPairTraversalData *data) {
#ifdef QUICKSCHED
    if (QuickschedTraversalData *qui_data = dynamic_cast<QuickschedTraversalData *>(data)) {
        CellPairTraversals<CellTemplate>::rebuild(cells, dims, cellLength, cutoff, data);
        qsched_reset(_scheduler);
        _taskBlocksize = qui_data->taskBlockSize;
        init();
    } else {
        Log::global_log->error() << "QuickschedTraversal::rebuild was called with incompatible Traversal data!" << std::endl;
    }
#endif /* QUICKSCHED */
}

template<class CellTemplate>
void QuickschedTraversal<CellTemplate>::init() {
#ifdef QUICKSCHED
    qsched_res_t  resourceId;
    qsched_task_t taskId;
    unsigned long cellIndex;
    // macro for easier access and to avoid aliasing
//#define m_cells (*((std::vector<ParticleCellBase> *)(this->_cells)))
//    std::vector<ParticleCellBase> m_cells = *(dynamic_cast<std::vector<ParticleCellBase> *>(this->_cells));
    std::vector<ParticleCell> m_cells = *((std::vector<ParticleCell> *)(this->_cells));

    switch (_taskTypeSelector) {
        case PackedAdjustable: {
            // check that blocksize is within domain size
            for (int i = 0; i < 3; ++i) {
                if (_taskBlocksize[i] > this->_dims[i]) {
                    std::ostringstream error_message;
                    error_message << "Blocksize is bigger than number of cells in dimension "
                                        << (char) ('x' + i) << ". (" << _taskBlocksize[i] << " > "
                                        << this->_dims[i] << ")" << std::endl;
                    MARDYN_EXIT(error_message.str());
                }
            }

            Log::global_log->info() << "Generating resource and task ids" << std::endl;
            for (unsigned long z = 0; z < this->_dims[2]; ++z) {
                for (unsigned long y = 0; y < this->_dims[1]; ++y) {
                    for (unsigned long x = 0; x < this->_dims[0]; ++x) {
                        cellIndex  = threeDimensionalMapping::threeToOneD(x, y, z, this->_dims);
                        resourceId = qsched_addres(_scheduler, qsched_owner_none, qsched_res_none);
                        m_cells[cellIndex].setResourceId(resourceId);
                        // only create tasks with offset blocksize-1.
                        // -1 because they need to overlap
                        // skip tasks for rear halo layer as they would only contain halo cells
                        if ((z % (_taskBlocksize[2] - 1) == 0
                             && y % (_taskBlocksize[1] - 1) == 0
                             && x % (_taskBlocksize[0] - 1) == 0)
                            &&
                            (x < this->_dims[0] - 1
                             && y < this->_dims[1] - 1
                             && z < this->_dims[2] - 1)) {
                            // also save the pointers as long
                            unsigned long payload[]{x, y, z, (unsigned long) this};
                            taskId = qsched_addtask(_scheduler,
                                                    PackedAdjustable,
                                                    task_flag_none,
                                                    payload,
                                                    sizeof(payload),
                                                    1);
                            m_cells[cellIndex].setTaskId(taskId);
                        }
                    } /* end for-x */
                } /* end for-y*/
            } /* end for-z */

            // set dependencies
            Log::global_log->info() << "Setting task dependencies" << std::endl;
            for (unsigned long z = 0; z < this->_dims[2] - 1; z += _taskBlocksize[2] - 1) {
                for (unsigned long y = 0; y < this->_dims[1] - 1; y += _taskBlocksize[1] - 1) {
                    for (unsigned long x = 0; x < this->_dims[0] - 1; x += _taskBlocksize[0] - 1) {
                        cellIndex = threeDimensionalMapping::threeToOneD(x, y, z, this->_dims);

                        // create locks for all 8 corners of the block
                        for (unsigned long i = 0; i < _taskBlocksize[0]
                                                  && x + i < this->_dims[0]; i += _taskBlocksize[0] - 1) {
                            for (unsigned long j = 0; j < _taskBlocksize[1]
                                                      && y + j < this->_dims[1]; j += _taskBlocksize[1] - 1) {
                                for (unsigned long k = 0; k < _taskBlocksize[2]
                                                          && z + k < this->_dims[2]; k += _taskBlocksize[2] - 1) {
                                    qsched_addlock(_scheduler,
                                                   m_cells[cellIndex].getTaskId(),
                                                   m_cells[threeDimensionalMapping::threeToOneD(x + i,
                                                                                                y + j,
                                                                                                z + k,
                                                                                                this->_dims)].getRescourceId());
                                }
                            }
                        }
                    } /* end for-x */
                } /* end for-y*/
            } /* end for-z */
            break;
        } /* end case PackedAdjustable */
        default:
            Log::global_log->error() << "QuickschedHandler::init() received non existing task type!" << std::endl;
    }
#endif // QUICKSCHED
}

template<class CellTemplate>
QuickschedTraversal<CellTemplate>::~QuickschedTraversal() {
#ifdef QUICKSCHED
    qsched_free(_scheduler);
    delete (_scheduler);
#endif // QUICKSCHED
}

template<class CellTemplate>
void QuickschedTraversal<CellTemplate>::runner(int type, void *data) {
#ifdef QUICKSCHED

    QuickschedTraversal *context = ((QuickschedTraversal **) data)[3];
#ifdef PRINT_SCHEDULING_TIMINGS
    struct VectorizedCellProcessor::Timings timing ;
            if(_simulation.getSimStep() > 10){
                timing.start = _rdtsc();
            }
#endif
    switch (type) {
        case PackedAdjustable: {
            // TODO optimize calculation order (1. corners 2. edges 3. rest) and gradually release resources
            unsigned long x = ((unsigned long *) data)[0];
            unsigned long y = ((unsigned long *) data)[1];
            unsigned long z = ((unsigned long *) data)[2];

            // traverse over block
            for (unsigned long i = 0; i < context->_taskBlocksize[0] - 1
                                      && i < context->_dims[0] - 1; ++i) {
                for (unsigned long j = 0; j < context->_taskBlocksize[1] - 1
                                          && j < context->_dims[1] - 1; ++j) {
                    for (unsigned long k = 0; k < context->_taskBlocksize[2] - 1
                                              && k < context->_dims[2] - 1; ++k) {

                        //process cell
                        unsigned long baseIndex = threeDimensionalMapping::threeToOneD(x + i,
                                                                                       y + j,
                                                                                       z + k,
                                                                                       context->_dims);
                        context->processBaseCell(
                                *(context->_contextCellProcessor),
                                baseIndex);
                    }
                }
            }
            break;
        } /* end case PackedAdjustable */
        default:
            Log::global_log->error() << "Undefined Quicksched task type: " << type << std::endl;
    }
#ifdef PRINT_SCHEDULING_TIMINGS
    if(_simulation.getSimStep() > 10){
        timing.end = _rdtsc();
        (dynamic_cast<VectorizedCellProcessor*>(_contextCellProcessor))->getThreadData()[omp_get_thread_num()]->_timings.push_back(timing);
    }
#endif
#endif /* QUICKSCHED */
}

template<class CellTemplate>
void QuickschedTraversal<CellTemplate>::traverseCellPairs(CellProcessor &cellProcessor) {
    _contextCellProcessor = &cellProcessor;
#ifdef QUICKSCHED
    qsched_run(_scheduler, mardyn_get_max_threads(), runner);
#endif /* QUICKSCHED */
}

template<class CellTemplate>
void QuickschedTraversal<CellTemplate>::traverseCellPairsInner(CellProcessor &cellProcessor, unsigned stage,
                                                               unsigned stageCount) {
    Log::global_log->error() << "QuickschedTraversal::traverseCellPairsInner is not implemented!" << std::endl;
}

template<class CellTemplate>
void QuickschedTraversal<CellTemplate>::traverseCellPairsOuter(CellProcessor &cellProcessor) {
    Log::global_log->error() << "QuickschedTraversal::traverseCellPairsOuter is not implemented!" << std::endl;
}

#endif //SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_QUICKSCHEDTRAVERSAL_H

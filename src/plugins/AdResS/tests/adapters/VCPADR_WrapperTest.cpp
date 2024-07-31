//
// Created by alex on 7/30/24.
//

#include "VCPADR_WrapperTest.h"
#include "../../adapters/VCPADR_Wrapper.h"
#include "plugins/AdResS/AdResS.h"

TEST_SUITE_REGISTRATION(VCPADR_WrapperTest);

void VCPADR_WrapperTest::processorSelectionTest() {
    // Init required objects
    const char* filename = "AdResS-empty-10x10x10.inp";
    ParticleContainer* container(initializeFromFile(ParticleContainerFactory::LinkedCell, filename, 2.5));
    _simulation.getEnsemble()->getComponent(0)->setName("FP_M");
    _simulation.getEnsemble()->getComponent(1)->setName("H_M");
    _simulation.getEnsemble()->getComponent(2)->setName("CG_M");

    std::unique_ptr<AdResS> plugin = std::make_unique<AdResS>();
    Resolution::Config resConf {
        _simulation.getEnsemble()->getComponents(),
        _domain
    };
    resConf.comp_to_res.resize(resConf.components->size(), Resolution::FullParticle);
    Resolution::FPRegion fpRegion { {4,4,4}, {6,6,6}, {1,1,1}};
    fpRegion.init();
    resConf.fpRegions.push_back(fpRegion);
    plugin->_resolutionHandler.init(resConf);

    // cells with lower coord >= (3,3,3) and upper coord <= (7,7,7) contain hybrid particles
    // => must be handled via AdResS Cell Processor
    VCPADR_Wrapper cell_processor {*_domain, 2.5, 2.5, plugin->_resolutionHandler};
    cell_processor.init();
    cell_processor._adr_processor = std::make_unique<CellProcMock>(2.5, 2.5);
    cell_processor._reference_processor = std::make_unique<CellProcMock>(2.5, 2.5);

    // perform actual testing...
    container->traverseCells(cell_processor);

    for (auto* cell_ptr : dynamic_cast<CellProcMock*>(cell_processor._adr_processor.get())->get_cells()) {
        for (int dim = 0; dim < 3; dim++) {
            ASSERT_TRUE(cell_ptr->getBoxMin(dim) >= 2.5);
            ASSERT_TRUE(cell_ptr->getBoxMax(dim) <= 7.5);
        }
    }

    for (auto* cell_ptr : dynamic_cast<CellProcMock*>(cell_processor._reference_processor.get())->get_cells()) {
        bool not_H = false;
        for (int dim = 0; dim < 3; dim++) {
            not_H |= !(cell_ptr->getBoxMin(dim) >= 2.5 && cell_ptr->getBoxMax(dim) <= 7.5);
        }
        ASSERT_TRUE(not_H);
    }

    for (auto cell_ptr_arr : dynamic_cast<CellProcMock*>(cell_processor._adr_processor.get())->get_cellPairs()) {
        // at least one of these cells must be in H region
        bool in_H_0 = false;
        bool in_H_1 = false;
        for (int dim = 0; dim < 3; dim++) {
            in_H_0 |= (cell_ptr_arr[0]->getBoxMin(dim) >= 2.5) && (cell_ptr_arr[0]->getBoxMax(dim) <= 7.5);
            in_H_1 |= (cell_ptr_arr[1]->getBoxMin(dim) >= 2.5) && (cell_ptr_arr[1]->getBoxMax(dim) <= 7.5);
        }

        ASSERT_TRUE(in_H_0 | in_H_1);
    }

    for (auto cell_ptr_arr : dynamic_cast<CellProcMock*>(cell_processor._reference_processor.get())->get_cellPairs()) {
        // at least one of these cells must be in H region
        bool in_H_0 = false;
        bool in_H_1 = false;
        for (int dim = 0; dim < 3; dim++) {
            in_H_0 |= (cell_ptr_arr[0]->getBoxMin(dim) >= 2.5) && (cell_ptr_arr[0]->getBoxMax(dim) <= 7.5);
            in_H_1 |= (cell_ptr_arr[1]->getBoxMin(dim) >= 2.5) && (cell_ptr_arr[1]->getBoxMax(dim) <= 7.5);
        }

        ASSERT_TRUE(!(in_H_0 | in_H_1));
    }

    delete container;
}

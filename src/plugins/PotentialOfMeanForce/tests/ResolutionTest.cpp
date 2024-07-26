#include "ResolutionTest.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"

TEST_SUITE_REGISTRATION(ResolutionTest);

void ResolutionTest::testCheckResolution(){
    
    const std::string file_name = "pmf.inp";
    ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell,file_name,5);

    std::cout<<container->getCellLength()<<std::endl;
}

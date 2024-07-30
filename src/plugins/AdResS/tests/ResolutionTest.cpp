#include "ResolutionTest.h"
#include "particleContainer/ParticleContainer.h"

TEST_SUITE_REGISTRATION(ResolutionTest);

ResolutionTest::ResolutionTest(){
    //Config config;
    // config.domain=_simulation.getDomain();
    FPRegion region{{4,0,0},{6,10,10},{2,2,2}};
    region.init();
    config.fpRegions.emplace_back(region);
    Component fp{0},hy{1},cg{2};
    fp.setName("FP_Dummy");
    fp.addLJcenter(0,0,0,1,1,1,0);
    fp.addLJcenter(0,0,0,1,1,1,0);
    
    hy.setName("H_Dummy");
    hy.addLJcenter(0,0,0,1,1,1,0);
    hy.addLJcenter(0,0,0,1,1,1,0);
    hy.addLJcenter(0,0,0,1,1,1,0);

    cg.setName("CG_Dummy");
    cg.addLJcenter(0,0,0,1,1,1,0);

    config.components = new (std::vector<Component>);
    config.components->emplace_back(fp);
    config.components->emplace_back(hy);
    config.components->emplace_back(cg);

    handler.init(config);
}

void ResolutionTest::testCheckResolution(){
    Component fp = (config.components)->at(0);
    Component cg = (config.components)->at(2);

    Molecule m1{1,&fp,5,5,5,0,0,0,0,0,0,0,0,0,0};
    Molecule m2{2,&fp,1,1,1,0,0,0,0,0,0,0,0,0,0};
    Molecule m3{3,&fp,7,7,7,0,0,0,0,0,0,0,0,0,0};
    std::array<double, 3> boxMin = {0.0, 0.0, 0.0};
	std::array<double, 3> boxMax = {10.0, 10.0, 10.0};
    LinkedCells lc{boxMin.data(),boxMax.data(),5};
    lc.addParticle(m1);
    lc.addParticle(m2);
    lc.addParticle(m3);
    handler.checkResolution(lc);

    std::array<double,3> pos = {5,5,5};
    std::variant<ParticleIterator, SingleCellIterator<ParticleCell>> v = lc.getMoleculeAtPosition(pos.data());
    SingleCellIterator<ParticleCell> m_temp = std::get<SingleCellIterator<ParticleCell>>(v);

    ASSERT_EQUAL(m_temp->componentid(),0u);

    pos = {1,1,1};
    v = lc.getMoleculeAtPosition(pos.data());
    m_temp = std::get<SingleCellIterator<ParticleCell>>(v);

    ASSERT_EQUAL(m_temp->componentid(),2u);


    pos = {7,7,7};
    v = lc.getMoleculeAtPosition(pos.data());
    m_temp = std::get<SingleCellIterator<ParticleCell>>(v);

    ASSERT_EQUAL(m_temp->componentid(),1u);


}
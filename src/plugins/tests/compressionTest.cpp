//
// Created by fernanor on 2018-12-13
//

#include "compressionTest.h"


TEST_SUITE_REGISTRATION(compressionTest);

compressionTest::compressionTest() {}

compressionTest::~compressionTest() {}

void compressionTest::testNone() {
    std::unique_ptr<Compression> compInstance = Compression::create("None");
    ASSERT_TRUE(compInstance.get());

    constexpr size_t datasize = 1024*1024;
    std::vector<char> someData(datasize);
    std::mt19937 gen(0);
    std::normal_distribution<> charsonly(0, 255);
    for (auto& elem : someData) {
        elem = charsonly(gen);
    }    
    std::vector<char> compressedData;
    std::vector<char> decompressedData;

    //compression
    compInstance->compress(someData.begin(), someData.end(), compressedData);
    size_t correctUncompressedSize = datasize;
    size_t correctCompressedSize = datasize;
    ASSERT_EQUAL(someData.size(), compInstance->getUncompressedSize());
    ASSERT_EQUAL(correctUncompressedSize, compInstance->getUncompressedSize());
    ASSERT_EQUAL(compressedData.size(), compInstance->getCompressedSize());
    ASSERT_EQUAL(correctCompressedSize, compInstance->getCompressedSize());

    //decompression
    compInstance->decompress(compressedData.begin(), compressedData.end(), decompressedData);
    ASSERT_EQUAL(correctUncompressedSize, decompressedData.size());

    //assert symmetry
    for (auto i=0; i<decompressedData.size(); ++i) {
        ASSERT_EQUAL(decompressedData[i], someData[i]);
    }
}

#ifdef ENABLE_LZ4
void compressionTest::testLz4() {
    std::unique_ptr<Compression> compInstance = Compression::create("LZ4");
    ASSERT_TRUE(compInstance.get());

    constexpr size_t datasize = 1024*1024;
    std::vector<char> someData(datasize);
    std::mt19937 gen(0);
    std::normal_distribution<> charsonly(0, 255);
    for (auto& elem : someData) {
        elem = charsonly(gen);
    }    
    std::vector<char> compressedData;
    std::vector<char> decompressedData;

    //compression
    compInstance->compress(someData.begin(), someData.end(), compressedData);
    size_t correctUncompressedSize = datasize;
    std::cout << "acs:" << compInstance->getCompressedSize() << std::endl;
    // ASSERT_EQUAL(someData.size(), compInstance->getUncompressedSize());
    // ASSERT_EQUAL(correctUncompressedSize, compInstance->getUncompressedSize());
    // ASSERT_EQUAL(compressedData.size(), compInstance->getCompressedSize());
    // ASSERT_EQUAL(correctCompressedSize, compInstance->getCompressedSize());

    //decompression
    // compInstance->decompress(compressedData.begin(), compressedData.end(), decompressedData);
    // ASSERT_EQUAL(correctUncompressedSize, decompressedData.size());

    // assert symmetry
    // for (auto i=0; i<decompressedData.size(); ++i) {
    //     ASSERT_EQUAL(decompressedData[i], someData[i]);
    // }
}
#endif

void compressionTest::testThrowOnFailedCreate() {
    try {
        std::unique_ptr<Compression> compInstance = Compression::create("FailMe");
    }
    catch (std::invalid_argument ia) {
        std::string const correctMsg("CompressionWrapper error > Invalid encoding: FailMe");
        std::string const exceptionMsg(ia.what());
        ASSERT_EQUAL(correctMsg, exceptionMsg);
    }
}

/*
    const char* filename = "1clj-regular-2x2x2-offset.inp";
    double cutoff = .5;
    ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff);


    COMaligner* plugin = new COMaligner();

    plugin->init(container, _domainDecomposition, _domain);
    plugin->beforeForces(container, _domainDecomposition, 1);

    double m = plugin->_mass;
    if (_domainDecomposition->getNumProcs() != 1) {
        test_log->info() << "compressionTest::testCOMalign: Mass Check SKIPPED (required exactly 1 process but was run with " <<  _domainDecomposition->getNumProcs() << " processes)" << std::endl;
    }
    else{
		double expectedMass;
		expectedMass = 8.0;
		ASSERT_EQUAL_MSG("Mass does not match number of particles", expectedMass, m);
    }

    // TEST MOTION
    ASSERT_EQUAL_MSG("x motion is wrong", -.25, plugin->_motion[0]);
    ASSERT_EQUAL_MSG("y motion is wrong", -.25, plugin->_motion[1]);
    ASSERT_EQUAL_MSG("z motion is wrong", -.25, plugin->_motion[2]);

    // initialize oldContainer only now, to prevent it from interfering with anything relevant!
    ParticleContainer* oldContainer = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff);
    // TEST IF MOTION WAS APPLIED
    auto newPos = container->iterator();
    auto oldPos = oldContainer->iterator();
    while(newPos.isValid()){
        for(int d = 0; d < 3; d++){
            ASSERT_EQUAL_MSG("Motion has not been properly applied" ,oldPos->r(d) - .25, newPos->r(d));
        }
        ++newPos;
        ++oldPos;
    }*/

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
    std::uniform_int_distribution<> charsonly(0, 255);
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
void compressionTest::testLz4Random() {
    std::unique_ptr<Compression> compInstance = Compression::create("LZ4");
    ASSERT_TRUE(compInstance.get());

    constexpr size_t datasize = 1024*1024;
    std::vector<char> someData(datasize);
    //fill input with random data
    std::mt19937 gen(0);
    std::uniform_int_distribution<> charsonly(0, 255);
    for (auto& elem : someData) {
        elem = charsonly(gen);
    }    
    std::vector<char> compressedData;
    std::vector<char> decompressedData;

    //compression
    compInstance->compress(someData.begin(), someData.end(), compressedData);
    size_t correctUncompressedSize = datasize;
    ASSERT_EQUAL(someData.size(), compInstance->getUncompressedSize());
    ASSERT_EQUAL(correctUncompressedSize, compInstance->getUncompressedSize());
    ASSERT_EQUAL(compressedData.size(), compInstance->getCompressedSize()); //this is actually larger than uncompressed, Kudos to mt19937 entropy :P

    //decompression
    compInstance->decompress(compressedData.begin(), compressedData.end(), decompressedData);
    ASSERT_EQUAL(correctUncompressedSize, decompressedData.size());

    // assert symmetry
    for (auto i=0; i<decompressedData.size(); ++i) {
        ASSERT_EQUAL(decompressedData[i], someData[i]);
    }
}

void compressionTest::testLz4Sine() {
    std::unique_ptr<Compression> compInstance = Compression::create("LZ4");
    ASSERT_TRUE(compInstance.get());

    constexpr size_t datasize = 1024*1024;
    std::vector<char> someData(datasize);
    // fill with 4 full sine periods
    constexpr double dphi = 8.*M_PI/datasize;
    int t = 0;
    for (auto& elem : someData) {
        elem = static_cast<char>(sin(dphi*t++/datasize)*127.);
    }    
    std::vector<char> compressedData;
    std::vector<char> decompressedData;

    //compression
    compInstance->compress(someData.begin(), someData.end(), compressedData);
    size_t correctUncompressedSize = datasize;
    ASSERT_EQUAL(someData.size(), compInstance->getUncompressedSize());
    ASSERT_EQUAL(correctUncompressedSize, compInstance->getUncompressedSize());
    ASSERT_EQUAL(compressedData.size(), compInstance->getCompressedSize());

    //decompression
    compInstance->decompress(compressedData.begin(), compressedData.end(), decompressedData);
    ASSERT_EQUAL(correctUncompressedSize, decompressedData.size());

    // assert symmetry
    for (auto i=0; i<decompressedData.size(); ++i) {
        ASSERT_EQUAL(decompressedData[i], someData[i]);
    }
}
#endif

void compressionTest::testThrowOnFailedCreate() {
    try {
        std::unique_ptr<Compression> compInstance = Compression::create("FailMe");
    }
    catch (const std::invalid_argument& ia) {
        std::string const correctMsg("CompressionWrapper error > Invalid encoding: FailMe");
        std::string const exceptionMsg(ia.what());
        ASSERT_EQUAL(correctMsg, exceptionMsg);
    }
}
#include "compression.h"

int Lz4Compression::compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) {
    std::cout << "Compressing using LZ4" << std::endl;
    size_t uncompressedSize = uncompressedEnd-uncompressedStart;
    compressed.resize(uncompressedSize);
    LZ4_compress_default(&(*uncompressedStart), compressed.data(), uncompressedSize, uncompressedSize);
    return 0;
}

int Lz4Compression::decompress(ByteIterator compressedStart, ByteIterator compressedEnd, std::vector<char>& decompressed) {
    std::cout << "Decompressing LZ4 data" << std::endl;
    size_t compressedSize = compressedEnd-compressedStart;
    size_t decompressedSize = 10000000;
    decompressed.resize(decompressedSize);
    LZ4_decompress_safe (&(*compressedStart), decompressed.data(), compressedSize, decompressedSize);
    return 0;
}

int NoCompression::compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) {
    std::cout << "Compressing" << std::endl;
    compressed.resize(uncompressedEnd-uncompressedStart);
    auto curPosCompressed = compressed.begin();
    while (uncompressedStart != uncompressedEnd) {
        *curPosCompressed = *uncompressedStart;
        ++curPosCompressed;
        ++uncompressedStart;
    }
    return 0;
}

int NoCompression::decompress(ByteIterator compressedStart, ByteIterator compressedEnd, std::vector<char>& decompressed) {
    std::cout << "Decompressing" << std::endl;
    decompressed.resize(compressedEnd-compressedStart);
    auto curPosDecompressed = decompressed.begin();
    while (compressedStart != compressedEnd) {
        *curPosDecompressed = *compressedStart;
        ++curPosDecompressed;
        ++compressedStart;
    }
    return 0;
}
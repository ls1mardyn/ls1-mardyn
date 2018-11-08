#include "compression.h"

int Lz4Compression::compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) {
    std::cout << "Compressing using LZ4" << std::endl;
    size_t uncompressedSize = uncompressedEnd-uncompressedStart;
    compressed.resize(uncompressedSize+sizeof(size_t));
    // compress the uncompressed source
    auto compressedSize = LZ4_compress_default(
            &(*uncompressedStart),
            compressed.data()+sizeof(size_t),
            uncompressedSize,
            uncompressedSize);
    // add the uncompressed size to the beginning of the result array
    std::cout << "Size after compressing: " << compressedSize << std::endl;
    std::copy(reinterpret_cast<char*>(&uncompressedSize),
            reinterpret_cast<char*>(&uncompressedSize)+sizeof(size_t),
            compressed.data());
    compressed.resize(sizeof(size_t)+compressedSize);
    return 0;
}

int Lz4Compression::decompress(ByteIterator compressedStart, ByteIterator compressedEnd, std::vector<char>& decompressed) {
    std::cout << "Decompressing LZ4 data" << std::endl;
    size_t compressedSize = compressedEnd-compressedStart-sizeof(size_t);
    // get the decompressed size
    size_t decompressedSize;
    std::copy(&(*compressedStart),
            &(*compressedStart)+sizeof(size_t),
            reinterpret_cast<char*>(&decompressedSize));
    std::cout << "Decompressed size read from data chunk: " << decompressedSize << std::endl; 
    decompressed.resize(decompressedSize);
    LZ4_decompress_safe(&(*compressedStart)+sizeof(size_t), decompressed.data(), compressedSize, decompressedSize);
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
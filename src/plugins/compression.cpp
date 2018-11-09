#include "compression.h"

std::unique_ptr<Compression> Compression::create(std::string encoding) {
    if (encoding.compare("LZ4") == 0) {
        return std::make_unique<Lz4Compression>();
    }
    else if (encoding.compare("None") == 0) {
        return std::make_unique<NoCompression>();
    }
    else {
        return std::unique_ptr<Compression>(nullptr);
    }
}

int Lz4Compression::compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) {
    uncompressedSize = uncompressedEnd-uncompressedStart;
    compressed.resize(uncompressedSize+sizeof(size_t));
    // compress the uncompressed source
    compressedSize = LZ4_compress_default(
            &(*uncompressedStart),
            compressed.data()+sizeof(size_t),
            uncompressedSize,
            uncompressedSize);
    if (compressedSize < 0) {
        error = CompressionError::COMP_ERR_COMPRESSION_FAILED;
        return static_cast<int>(compressedSize);
    }
    // add the uncompressed size to the beginning of the result array
    std::copy(reinterpret_cast<char*>(&uncompressedSize),
            reinterpret_cast<char*>(&uncompressedSize)+sizeof(size_t),
            compressed.data());
    compressedSize += sizeof(size_t);
    compressed.resize(compressedSize);
    error = CompressionError::COMP_SUCCESS;
    return 0;
}

int Lz4Compression::decompress(ByteIterator compressedStart, ByteIterator compressedEnd, std::vector<char>& decompressed) {
    compressedSize = compressedEnd-compressedStart;
    // get the decompressed size
    std::copy(&(*compressedStart),
            &(*compressedStart)+sizeof(size_t),
            reinterpret_cast<char*>(&uncompressedSize));
    decompressed.resize(uncompressedSize);
    auto actualDecompressedSize = LZ4_decompress_safe(&(*compressedStart)+sizeof(size_t),
            decompressed.data(),
            compressedSize-sizeof(size_t),
            uncompressedSize);
    if (actualDecompressedSize != uncompressedSize) {
        error = CompressionError::COMP_ERR_DECOMPRESSION_FAILED;
        return static_cast<int>(actualDecompressedSize);
    }
    error = CompressionError::COMP_SUCCESS;
    return 0;
}

int NoCompression::compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) {
    uncompressedSize = uncompressedEnd-uncompressedStart;
    compressedSize = uncompressedSize;
    compressed.resize(compressedSize);
    auto curPosCompressed = compressed.begin();
    while (uncompressedStart != uncompressedEnd) {
        *curPosCompressed = *uncompressedStart;
        ++curPosCompressed;
        ++uncompressedStart;
    }
    error = CompressionError::COMP_SUCCESS;
    return 0;
}

int NoCompression::decompress(ByteIterator compressedStart, ByteIterator compressedEnd, std::vector<char>& decompressed) {
    compressedSize = compressedEnd-compressedStart;
    uncompressedSize = compressedSize;
    decompressed.resize(uncompressedSize);
    auto curPosDecompressed = decompressed.begin();
    while (compressedStart != compressedEnd) {
        *curPosDecompressed = *compressedStart;
        ++curPosDecompressed;
        ++compressedStart;
    }
    error = CompressionError::COMP_SUCCESS;
    return 0;
}
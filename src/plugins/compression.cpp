#include "compression.h"

std::unique_ptr<Compression> Compression::create(std::string encoding) {
#ifdef ENABLE_LZ4
    if (encoding.compare("LZ4") == 0) {
        return std::unique_ptr<Lz4Compression>(new Lz4Compression());
    }
    else
#endif /* ENABLE_LZ4 */
    if (encoding.compare("None") == 0) {
        return std::unique_ptr<NoCompression>(new NoCompression());
    }
    else {
        std::string const errormsg("CompressionWrapper error > Invalid encoding: "+encoding);
        throw std::invalid_argument(errormsg);
    }
}

#ifdef ENABLE_LZ4
void Lz4Compression::compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) {
    decltype(_uncompressedSize) uncompressedSize = uncompressedEnd-uncompressedStart;
    decltype(_compressedSize) compressedSize;
    // create a temp vector to store the compressed data
    std::vector<char> temp;
    temp.resize(uncompressedSize+sizeof(decltype(_compressedSize)));
    // compress the uncompressed source
    auto actualCompressedSize = LZ4_compress_default(
            &(*uncompressedStart),
            temp.data()+sizeof(decltype(_compressedSize)),
            uncompressedSize,
            uncompressedSize);
    //following may throw (duh)
    if (actualCompressedSize < 0) {
        std::string const errormsg(
            "CompressionWrapper error > LZ4_compress_default failed with error code: "
            +std::to_string(actualCompressedSize));
        throw std::runtime_error(errormsg);
    }
    else {
        compressedSize = actualCompressedSize;
    }
    // add the uncompressed size to the beginning of the result array as meta data(may throw)
    std::copy(reinterpret_cast<char*>(&uncompressedSize),
            reinterpret_cast<char*>(&uncompressedSize)+sizeof(decltype(_uncompressedSize)),
            temp.data());
    // add these extra bytes for the meta data to the tracked size
    compressedSize += sizeof(decltype(_uncompressedSize));
    temp.resize(compressedSize);
    // If vector allocators match (and they have to as method is not templated),
    // swap can't throw -> method has strong guarantee
    _uncompressedSize = uncompressedSize;
    _compressedSize = compressedSize;
    compressed.swap(temp);
}

void Lz4Compression::decompress(ByteIterator compressedStart, ByteIterator compressedEnd, std::vector<char>& decompressed) {
    decltype(_compressedSize) compressedSize = compressedEnd-compressedStart;
    decltype(_uncompressedSize) uncompressedSize;
    // get the decompressed size, may throw
    std::copy(&(*compressedStart),
            &(*compressedStart)+sizeof(decltype(_uncompressedSize)),
            reinterpret_cast<char*>(&uncompressedSize));
    // resize may throw
    std::vector<char> temp;
    temp.resize(uncompressedSize);
    auto actualDecompressedSize = LZ4_decompress_safe(&(*compressedStart)+sizeof(decltype(_uncompressedSize)),
            temp.data(),
            compressedSize-sizeof(decltype(_uncompressedSize)),
            uncompressedSize);
    // following may throw (duh)
    if (actualDecompressedSize != uncompressedSize) {
        std::string const errormsg(
            "CompressionWrapper error > LZ4_decompress_default failed with error code: "
            +std::to_string(actualDecompressedSize));
        throw std::runtime_error(errormsg);
    }
    // If vector allocators match (and they have to as method is not templated),
    // swap can't throw -> method has strong guarantee
    _compressedSize = compressedSize;
    _uncompressedSize = uncompressedSize;
    decompressed.swap(temp);
}
#endif /* ENABLE_LZ4 */

void NoCompression::compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) {
    // Everything is either const on operating on temporaries, until the swap -> strong guarantee
    std::vector<char> temp;
    auto uncompressedSize = uncompressedEnd-uncompressedStart;
    auto compressedSize = uncompressedSize;
    temp.resize(compressedSize);
    auto curPosCompressed = temp.begin();
    while (uncompressedStart != uncompressedEnd) {
        *curPosCompressed = *uncompressedStart;
        ++curPosCompressed;
        ++uncompressedStart;
    }
    _compressedSize = compressedSize;
    _uncompressedSize = uncompressedSize;
    compressed.swap(temp);
}

void NoCompression::decompress(ByteIterator compressedStart, ByteIterator compressedEnd, std::vector<char>& decompressed) {
    // Everything is either const on operating on temporaries, until the swap -> strong guarantee
    std::vector<char> temp;
    auto compressedSize = compressedEnd-compressedStart;
    auto uncompressedSize = compressedSize;
    temp.resize(uncompressedSize);
    auto curPosDecompressed = temp.begin();
    while (compressedStart != compressedEnd) {
        *curPosDecompressed = *compressedStart;
        ++curPosDecompressed;
        ++compressedStart;
    }
    _uncompressedSize = uncompressedSize;
    _compressedSize = compressedSize;
    decompressed.swap(temp);
}

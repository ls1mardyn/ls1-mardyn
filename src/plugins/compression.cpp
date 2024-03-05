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
    // create a temp vector to store the compressed data (may throw)
    std::vector<char> temp;
    temp.resize(LZ4_compressBound(uncompressedSize));
    // compress the uncompressed source (may throw)
    auto actualCompressedSize = LZ4_compress_default(
            &(*uncompressedStart),                           //input data begin of type char*
            temp.data()+sizeof(_compressedSize),             //output data pointer of type char*
            uncompressedSize,                                //input data size (decltype of _uncompressedSize, usually size_t), imp. cast to int
            temp.size());                                    //output data pointer of decltype(temp.size()), most probably size_t, imp. cast to int
    //following may throw (duh)
    if (actualCompressedSize == 0) {
        std::string const errormsg(
            "CompressionWrapper error > LZ4_compress_default failed with error code: "
            +std::to_string(actualCompressedSize));
        throw std::runtime_error(errormsg);
    }
    else {
        compressedSize = actualCompressedSize;
    }
    // add the uncompressed size to the beginning of the result array as meta data (may throw)
    std::copy(reinterpret_cast<char*>(&uncompressedSize),
            reinterpret_cast<char*>(&uncompressedSize)+sizeof(_uncompressedSize),
            temp.data());
    // add these extra bytes for the meta data to the tracked size
    compressedSize += sizeof(_uncompressedSize);
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
            &(*compressedStart)+sizeof(_uncompressedSize),
            reinterpret_cast<char*>(&uncompressedSize));
    // resize may throw
    std::vector<char> temp;
    temp.resize(uncompressedSize);
    auto actualDecompressedSize = LZ4_decompress_safe(&(*compressedStart)+sizeof(_uncompressedSize),
            temp.data(),
            compressedSize-sizeof(_uncompressedSize),
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

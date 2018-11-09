/*
 * Compression Library Wrapper
 *
 *  Created on: 08 Nov 2018
 *      Author: Oliver Fernandes
 */

#ifndef SRC_PLUGINS_COMPRESSION_H_
#define SRC_PLUGINS_COMPRESSION_H_

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>


#include "lz4.h"

using ByteIterator = std::vector<char>::iterator;

enum CompressionError {
	COMP_SUCCESS = 0,
	COMP_ERR_COMPRESSION_FAILED = 200,
	COMP_ERR_DECOMPRESSION_FAILED = 300
};

class Compression {
public:
	/**
	 * Compresses a series of bytes.
     * Call to compress data.
	 * @param[in] uncompressedStart Iterator pointing to the start of the array to be compressed
	 * @param[in] uncompressedEnd Iterator pointing to the element past the last element of the array to be compressed
	 * @param[in] compressed A std::vector<char> holding the compressed result. The vector will be appropriately resized.
     *                       Any previous contents will be destroyed.
	 */
    virtual int compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) = 0;
	/**
	 * Decompresses a series of bytes.
     * Call to decompress data.
	 * @param[in] compressedStart Iterator pointing to the start of the array to be decompressed
	 * @param[in] compressedEnd Iterator pointing to the element past the last element of the array to be decompressed
	 * @param[in] decompressed A std::vector<char> holding the decompressed result. The vector will be appropriately resized.
     *                         Any previous contents will be destroyed. The first sizeof(size_t) bytes hold the uncompressed
	 *                         size.
	 */
    virtual int decompress(ByteIterator compressedStart, ByteIterator compressedEnd, std::vector<char>& decompressed) = 0;
	/**
	 * Returns the uncompressed size of the data.
     * Returns the uncompressed size of the data. This will only give sane results if either compression or decompression has
	 * been attempted. Negative values mean that there is no data, or something went wrong during the processing calls.
	 * @return Uncompressed size of data. Depending on if the instance was used for compression or decompression, this is the
	 *         input or output size, respectively.
	 */
	long long getUncompressedSize(void) const {
		return uncompressedSize;
	};
	/**
	 * Returns the compressed size of the data.
     * Returns the compressed size of the data. This will only give sane results if either compression or decompression has
	 * been attempted. Negative values mean that there is no data, or something went wrong during the processing calls.
	 * @return Uncompressed size of data. Depending on if the instance was used for compression or decompression, this is the
	 *         output or input size, respectively.
	 */
	long long getCompressedSize(void) const {
		return compressedSize;
	};
	/**
	 * Returns the error state of the instance.
     * Returns the last error encountered. Should be checked after every call to compress/decompress.
	 * If this returns something else than CompressionError::COMP_SUCCESS, the return value of the compress/decompress call may
	 * help identifying the issue.
	 * @return Uncompressed size of data. Depending on if the instance was used for compression or decompression, this is the
	 *         output or input size, respectively.
	 */
	CompressionError getError(void) const {
		return error;
	}
	/**
	 * Create an instance of the Compression class.
     * Use this to instantiate an object which is able to do compression on arrays. The argument passed is a std::string
	 * containing the tag of the compression algorithm. At the moment, this is just "LZ4" and "None".
	 * @param[in] encoding A tag from a list of tags associated with various compression algorithms.
	 * @return A std::unique_ptr to the Compression object.
	 */
	static std::unique_ptr<Compression> create(std::string encoding);
protected:
	long long uncompressedSize = -1;
	long long compressedSize = -1;
	CompressionError error;
};

class NoCompression : public Compression {
public:
    int compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) override;
    int decompress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) override;
};

class Lz4Compression : public Compression {
public:
    int compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) override;
    int decompress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) override;
};

#endif /* SRC_PLUGINS_COMPRESSION_H_ */

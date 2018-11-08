/*
 * Compression Library Wrapper
 *
 *  Created on: 08 Nov 2018
 *      Author: Oliver Fernandes
 */

#ifndef SRC_PLUGINS_COMPRESSION_H_
#define SRC_PLUGINS_COMPRESSION_H_

#include <vector>
#include <iostream>
#include <fstream>

#include "lz4.h"

using ByteIterator = std::vector<char>::iterator;

class CompressionInterface {
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
     *                         Any previous contents will be destroyed.
	 */
    virtual int decompress(ByteIterator compressedStart, ByteIterator compressedEnd, std::vector<char>& decompressed) = 0;
};

class NoCompression : public CompressionInterface {
public:
    int compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) override;
    int decompress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) override;
};

class Lz4Compression : public CompressionInterface {
public:
    int compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) override;
    int decompress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) override;
};

#endif /* SRC_PLUGINS_COMPRESSION_H_ */

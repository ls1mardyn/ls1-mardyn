/*
 * Compression Library Wrapper
 *
 *  Created on: 08 Nov 2018
 *      Author: Oliver Fernandes
 */

#ifndef SRC_PLUGINS_COMPRESSION_H_
#define SRC_PLUGINS_COMPRESSION_H_

#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#ifdef ENABLE_LZ4
#include "lz4.h"
#endif

/** @brief The Compression class provides easy to use compression methods.
 *
 * The idea of this class is to provide a "drop-in" possibility to (de-)compress a std::vector<char>, without worrying
 * about to much boiler plate code. The usage is simple, consider the following example:
 *
 * \code{.cpp}
 *     #include "compression.h"
 *     std::unique_ptr<Compression> comp = Compression::create("LZ4"); // Possible tags are "LZ4" and "None" atm
 *     std::vector<char> myData = {'1','2','3','a','b','c'};
 *     std::vector<char> compressedData;
 *     std::vector<char> decompressedData;
 *     comp.compress(someData.begin(), someData.end(), compressedData);
 *     comp.decompress(compressedData.begin(), compressedData.end(), decompressedData);
 *     for (size_t i = 0; i<decompressedData.size(); ++i) {
 *         assert(decompressedData[i] == myData[i]);
 *     }
 * \endcode
 */

class Compression {
public:
	virtual ~Compression() = default;
    /**
     * Just a type alias to avoid some typing
     */
    using ByteIterator = std::vector<char>::iterator;
	/**
	 * Compresses a series of bytes.
     *
     * Call to compress data.
	 * @param[in] uncompressedStart Iterator pointing to the start of the array to be compressed
	 * @param[in] uncompressedEnd Iterator pointing to the element past the last element of the array to be compressed
	 * @param[in] compressed A std::vector<char> holding the compressed result. The vector will be appropriately resized.
     *                       Any previous contents will be destroyed.
     * @return Error codes. Returns 0 on success. Error codes are specific to the employed algorithm.
	 */
    virtual void compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) = 0;
	/**
	 * Decompresses a series of bytes.
     *
     * Call to decompress data.
	 * @param[in] compressedStart Iterator pointing to the start of the array to be decompressed
	 * @param[in] compressedEnd Iterator pointing to the element past the last element of the array to be decompressed
	 * @param[in] decompressed A std::vector<char> holding the decompressed result. The vector will be appropriately resized.
     *                         Any previous contents will be destroyed. The first sizeof(size_t) bytes hold the uncompressed
	 *                         size.
     * @return Error codes. Returns 0 on success. Error codes are specific to the employed algorithm.
	 */
    virtual void decompress(ByteIterator compressedStart, ByteIterator compressedEnd, std::vector<char>& decompressed) = 0;
	/**
	 * Returns the uncompressed size of the data.
     *
     * Returns the uncompressed size of the data. This will only give sane results if either compression or decompression has
	 * been attempted. Zero means that there is no data, or something went wrong during the processing calls.
	 * @return Uncompressed size of data. Depending on if the instance was used for compression or decompression, this is the
	 *         input or output size, respectively.
	 */
	size_t getUncompressedSize(void) const {
		return _uncompressedSize;
	};
	/**
	 * Returns the compressed size of the data.
     *
     * Returns the compressed size of the data. This will only give sane results if either compression or decompression has
	 * been attempted. Zero means that there is no data, or something went wrong during the processing calls.
	 * @return Uncompressed size of data. Depending on if the instance was used for compression or decompression, this is the
	 *         output or input size, respectively.
	 */
	size_t getCompressedSize(void) const {
		return _compressedSize;
	};
	/**
	 * Create an instance of the Compression class.
     *
     * Use this to instantiate an object which is able to do compression on arrays. The argument passed is a std::string
	 * containing the tag of the compression algorithm. At the moment, this is just "LZ4" and "None".
     * If the tag is not recognized a nullptr will be returned.
	 * @param[in] encoding A tag from a list of tags associated with various compression algorithms.
	 * @return A std::unique_ptr to the Compression object.
	 */
	static std::unique_ptr<Compression> create(std::string encoding);
protected:
	size_t _uncompressedSize = 0;
	size_t _compressedSize = 0;
};

class NoCompression : public Compression {
public:
	/**
	 * Copies a series of bytes.
     * As the class name says, it just does a copy.
	 * @param[in] compressedStart Iterator pointing to the start of the array to be copied.
	 * @param[in] compressedEnd Iterator pointing to the element past the last element of the array to be compressed
	 * @param[in] decompressed A std::vector<char> holding the copied result. The vector will be appropriately resized.
     *                         Any previous contents will be destroyed. The first sizeof(size_t) bytes hold the uncompressed
	 *                         size.
     * @return Error codes. Returns 0 on success. May provide additional information when getError() is not COMP_SUCCESS.
	 */
    void compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) override;
	/**
	 * Copies a series of bytes back.
     * Call this to decompress data.
	 * @param[in] compressedStart Iterator pointing to the start of the array to be decompressed
	 * @param[in] compressedEnd Iterator pointing to the element past the last element of the array to be decompressed
	 * @param[in] decompressed A std::vector<char> holding the copied result. The vector will be appropriately resized.
     *                         Any previous contents will be destroyed. The first sizeof(size_t) bytes hold the uncompressed
	 *                         size.
     * @return Error codes. Returns 0 on success. May provide additional information when getError() is not COMP_SUCCESS.
	 */
    void decompress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) override;
};

#ifdef ENABLE_LZ4
class Lz4Compression : public Compression {
public:
	/**
	 * Copies a series of bytes.
     * As the class name says, it just does a copy.
	 * @param[in] compressedStart Iterator pointing to the start of the array to be copied.
	 * @param[in] compressedEnd Iterator pointing to the element past the last element of the array to be compressed
	 * @param[in] decompressed A std::vector<char> holding the compressed result. The vector will be appropriately resized.
     *                         Any previous contents will be destroyed. The first sizeof(size_t) bytes hold the uncompressed
	 *                         size.
     * @return Error codes. Returns 0 on success. Error codes match the LZ4 documentation. May provide additional information
     *         when getError() is not COMP_SUCCESS.
	 */
    void compress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) override;
	/**
	 * Decompresses a series of bytes.
     * Call this to decompress data.
	 * @param[in] compressedStart Iterator pointing to the start of the array to be decompressed
	 * @param[in] compressedEnd Iterator pointing to the element past the last element of the array to be decompressed
	 * @param[in] decompressed A std::vector<char> holding the decompressed result. The vector will be appropriately resized.
     *                         Any previous contents will be destroyed. The first sizeof(size_t) bytes hold the uncompressed
	 *                         size.
     * @return Error codes. Returns 0 on success. Error codes match the LZ4 documentation. May provide additional information
     *         when getError() is not COMP_SUCCESS.
	 */
    void decompress(ByteIterator uncompressedStart, ByteIterator uncompressedEnd, std::vector<char>& compressed) override;
};
#endif /* ENABLE_LZ4 */

#endif /* SRC_PLUGINS_COMPRESSION_H_ */

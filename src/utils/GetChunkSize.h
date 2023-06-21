#pragma once

namespace chunk_size {

/**
 * Formula to calculate a chunk size for OpenMP loop parallelization based on the given magic numbers.
 * The chunk size is the number of loop iterations scheduled together as one task.
 * @param loop_size The number of iterations.
 * @param max_num_chunks There should be enough chunks so that load can be balanced, but not too many since scheduling creates overhead.
 * @param max_chunk_size When iterating the specified number of work items, the overhead should be negligible. 
 *  Bigger numbers should not produce any better results.
 * @return The chunk size based on the given numbers.
 * @note max_chunk_size This is the more dominant factor, i.e., if max_num_chunks would indicate bigger chunks,
 * max_chunk_size is restricting the chunk size.
 */
constexpr int getChunkSize(size_t loop_size, size_t max_num_chunks, size_t max_chunk_size) {
	return static_cast<int>(std::max(std::min(loop_size / max_num_chunks, max_chunk_size), 1ul));
}

}  // namespace chunk_size

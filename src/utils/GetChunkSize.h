#pragma once

namespace chunk_size {

/**
 * Get the chunk size based on the given magic numbers.
 * @param loop_size The number of iterations.
 * @param max_num_chunks should be big enough to provide sufficient scheduling, but not too big, as too much overhead
 * might exist.
 * @param max_chunk_size should not be bigger than needed. When iterating the specified number of work items, the
 * overhead should be negligible. Bigger numbers should not produce any better results.
 * @return The chunk size based on the given numbers.
 * @note max_chunk_size is the more dominant factor, i.e., if max_num_chunks would indicate bigger chunks,
 * max_chunk_size is restricting the chunk size.
 */
constexpr int getChunkSize(size_t loop_size, size_t max_num_chunks, size_t max_chunk_size) {
	return static_cast<int>(std::max(std::min(loop_size / max_num_chunks, max_chunk_size), 1ul));
}

}  // namespace chunk_size
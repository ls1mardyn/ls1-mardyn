#pragma once

#include <deque>

/**
 * Class that handles a history with a maximal number of elements of type T.
 * The oldest entries are overwritten if new entries are added and the maximal capacity has been reached.
 * @tparam T The type of elements to store.
 */
template <typename T>
class RotatingHistory {
public:
	explicit RotatingHistory(size_t capacity = 1ul) : _capacity{capacity} {};

	void insert(const T& t) {
		if (_current_size == _capacity) {
			_storage.push_back(t);
			_storage.pop_front();
		} else {
			_storage.push_back(t);
			++_current_size;
		}
	}

	auto begin() { return _storage.begin(); }

	auto begin() const { return _storage.cbegin(); }

	auto end() { return _storage.end(); }

	auto end() const { return _storage.cend(); }

private:
	/**
	 * The capacity of the history.
	 */
	size_t _capacity{0ul};

	/**
	 * The number of stored elements.
	 */
	size_t _current_size{0ul};

	/**
	 * The actual storage.
	 */
	std::deque<T> _storage;
};
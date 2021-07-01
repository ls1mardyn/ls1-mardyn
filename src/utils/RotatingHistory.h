#pragma once

#include <deque>

/**
 * Class that handles a history with a maximal number of elements of type T.
 * The oldest entries are overwritten if new entries are added and the maximal capacity has been reached.
 * @tparam T The type of elements to store.
 * @note This is not the most efficient implementation! For a more efficient implementation, a std::vector should be
 * used. In that case, either an own iterator class is needed or the order of the elements is lost through the
 * iterators.
 */
template <typename T>
class RotatingHistory {
public:
	/**
	 * Constructor.
	 * @param capacity The number of history entries that can be stored.
	 */
	explicit RotatingHistory(size_t capacity = 1ul) : _capacity{capacity} {};

	/**
	 * Sets the capacity to new_capacity.
	 * @param new_capacity
	 */
	void setCapacity(size_t new_capacity){
		_capacity = new_capacity;
		// shrink to fit:
		while(_current_size > _capacity){
			_storage.pop_front();
			--_current_size;
		}
	}

	/**
	 * Inserts an element.
	 * Potentially deletes old elements if the capacity was reached.
	 * @param t
	 */
	void insert(const T& t) {
		if (_current_size == _capacity) {
			_storage.push_back(t);
			_storage.pop_front();
		} else {
			_storage.push_back(t);
			++_current_size;
		}
	}

	/**
	 * Returns the begin iterator.
	 * This initially points to the oldest entry.
	 * @return The iterator.
	 */
	auto begin() { return _storage.begin(); }

	/**
	 * @copydoc RotatingHistory::begin()
	 * @note const version
	 */
	auto begin() const { return _storage.cbegin(); }

	/**
	 * Returns the end iterator.
	 * @return
	 */
	auto end() { return _storage.end(); }

	/**
	 * @copydoc RotatingHistory::end()
	 * @note const version
	 */
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
#pragma once

#include <any>
#include <functional>

/**
 * Class that wraps a function.
 */
class FunctionWrapper {
public:
	/**
	 * Default constructor.
	 */
	FunctionWrapper() = default;

	/**
	 * Move constructor.
	 * @param other
	 */
	FunctionWrapper(FunctionWrapper&& other) = default;

	/**
	 * Copy constructor.
	 * @param other
	 */
	FunctionWrapper(const FunctionWrapper& other) = default;

	/**
	 * Copy assignment operator.
	 * @param other
	 * @return
	 */
	FunctionWrapper& operator=(const FunctionWrapper& other) = default;

	/**
	 * Move assignment operator.
	 * @param other
	 * @return
	 */
	FunctionWrapper& operator=(FunctionWrapper&& other) = default;

	/**
	 * Constructor for use with a lambda or function pointer.
	 * Example syntax:
	 * - As assignment operator, which will store a std::function<int(int)> :
	 * \code
	 *   FunctionWrapper wrapper;
	 *   wrapper = [] (int i) { return 4 + i; };
	 * \endcode
	 * - As constructor, which will store a std::function<int()> :
	 * \code
	 *   FunctionWrapper wrapper([]{return 4;})
	 * \endcode
	 * @tparam T The type of the lambda or function pointer.
	 * @param myFunctionObject The function object, can, e.g., be a lambda.
	 */
	template <typename T>
	FunctionWrapper(T&& myFunctionObject) : _myFunction{std::function(myFunctionObject)} {}

	/**
	 * Get the stored function object.
	 * Example syntax:
	 * - to get a function without input that returns an unsigned long:
	 * \code
	 *   wrapper.get<unsigned long>();
	 * \endcode
	 * - to get a function with unsigned long input that returns an unsigned long:
	 * \code
	 *   wrapper.get<unsigned long, unsigned long>();
	 * \endcode
	 * - to get a function with unsigned long as input, but without output:
	 * \code
	 *   wrapper.get<void, unsigned long>();
	 * \endcode
	 * - to get a function with two int as input that returns another int:
	 * \code
	 *   wrapper.get<int, int, int>();
	 * \endcode
	 * @tparam R The return type.
	 * @tparam Args The arguments.
	 * @return A std::function<R(Args...)> object that wraps the stored function.
	 */
	template <class R, class... Args>
	auto get() const {
		return std::any_cast<std::function<R(Args...)>>(_myFunction);
	}

private:
	std::any _myFunction;
};
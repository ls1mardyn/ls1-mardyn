#pragma once

#include <any>
#include <functional>

class FunctionWrapper {
public:
	FunctionWrapper() = default;

	FunctionWrapper(FunctionWrapper&& other) = default;

	FunctionWrapper(const FunctionWrapper& other) = default;

	FunctionWrapper& operator=(const FunctionWrapper& other) = default;

	FunctionWrapper& operator=(FunctionWrapper&& other) = default;

	template <typename T>
	FunctionWrapper(T&& myLambda) : _myFunction{std::function(myLambda)} {}

	template <class R, class... Args>
	auto get() const {
		return std::any_cast<std::function<R(Args...)>>(_myFunction);
	}

private:
	std::any _myFunction;
};
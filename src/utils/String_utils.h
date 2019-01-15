#pragma once

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace string_utils {

template <typename T>
static std::string join(const std::vector<T>& vec, const std::string delimiter) {
	std::stringstream ss;
	for(auto elem : vec) {
		ss << elem;
		if (elem != vec.back()) {
			ss << delimiter;
		}
	}
	return ss.str();
}

/**
 * Converts a string to lower case.
 * @param s input string.
 * @return Copy of s in lower case.
 */
static std::string toLowercase(const std::string& s) {
	std::string ret;
	ret.resize(s.length());
	std::transform(s.begin(), s.end(), ret.begin(), ::tolower);
	return ret;
}

}  // namespace string_utils

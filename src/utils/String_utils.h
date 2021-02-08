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

/**
 * Removes all leading and trailing whitespaces from a std::string.
 * @param str
 * @param whitespace characters to be removed. Default is whitespaces and tabs
 * @return trimmed string
 */
static std::string trim(const std::string& str, const std::string &whitespace = " \t") {
	const auto strBegin = str.find_first_not_of(whitespace);
	if (strBegin == std::string::npos)
		return "";

	const auto strEnd = str.find_last_not_of(whitespace);
	const auto strRange = strEnd - strBegin + 1;

	return str.substr(strBegin, strRange);
}

/**
 * Split a given input string using a delimiter into a vector of substrings.
 * @param input The input string.
 * @param delimiter The delimiter.
 * @return Vector of substrings.
 */
static std::vector<std::string> split(const std::string& input, const char delimiter ){
	size_t last = 0;
	size_t next = 0;
	std::vector<std::string> output;
	while ((next = input.find(delimiter, last)) != std::string::npos) {
		output.emplace_back(input.substr(last, next - last));
		last = next + 1;
	}
	output.emplace_back(input.substr(last));
	return output;
}

}  // namespace string_utils

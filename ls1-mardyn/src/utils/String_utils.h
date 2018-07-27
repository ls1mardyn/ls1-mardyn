#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <string>
#include <sstream>
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

}

#endif

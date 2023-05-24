#pragma once

#include <string>

const std::string MARDYN_VERSION_MAJOR = "1";
const std::string MARDYN_VERSION_MINOR = "2";
const std::string MARDYN_VERSION_PATCH = "0";
const std::string MARDYN_VERSION_BRANCH = "_master";
// this is the current git commit id
const std::string MARDYN_VERSION_HASH = "_3a847be73";
// append "dirty" if uncommited code was used
const std::string MARDYN_VERSION_IS_DIRTY = "_dirty";
// final combined version string
const std::string MARDYN_VERSION =
	MARDYN_VERSION_MAJOR + "."
	+ MARDYN_VERSION_MINOR + "."
	+ MARDYN_VERSION_PATCH
	+ MARDYN_VERSION_BRANCH
	+ MARDYN_VERSION_HASH
	+ MARDYN_VERSION_IS_DIRTY;

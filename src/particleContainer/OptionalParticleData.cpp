/**
 * @file OptionalParticleData.cpp
 * @author seckler
 * @date 13.08.20
 */

#include <tuple>

#include "OptionalParticleData.h"
void OptionalParticleData::createMap(const std::string& dataName,
									 OptionalParticleDataTypes::SupportedTypesVariant defaultValue) {
	// Add default value and check whether it doesn't exist, yet.
	auto [iterator, inserted] = _defaultValues.emplace(dataName, defaultValue);
	if (not inserted) {
		throw std::runtime_error("OptionalParticleData::createMap: map with key " + dataName + " already exists.");
	}
	std::visit(
		[&](auto value) {
			// Creates map with correct type.
			_particleData.emplace(dataName, std::unordered_map<size_t, decltype(value)>(value));
		},
		defaultValue);
}

void OptionalParticleData::createParticleDataWithDefaultValues(size_t id) {
	auto defaultValueIter = _defaultValues.begin();
	auto mapIter = _particleData.begin();
	for (; mapIter != _particleData.end(); ++mapIter, ++defaultValueIter) {
		std::visit(
			[=](auto& dataMap) {
				using map_type = std::remove_reference_t<decltype(dataMap)>;
				// Emplace default value with key id.
				dataMap.emplace(id, std::get<typename map_type::mapped_type>(defaultValueIter->second));
			},
			mapIter->second);
	}
}

void OptionalParticleData::removeParticleData(size_t id) {
	auto defaultValueIter = _defaultValues.begin();
	auto mapIter = _particleData.begin();
	for (; mapIter != _particleData.end(); ++mapIter, ++defaultValueIter) {
		std::visit(
			[=](auto& dataMap) {
				// Remove value with key id.
				dataMap.erase(id);
			},
			mapIter->second);
	}
}

OptionalParticleDataTypes::SupportedTypesVariant OptionalParticleData::getData(const std::string& dataName, size_t id) {
	return std::visit(
		[=](auto& dataMap) {
			// Get value of map with key id.
			return dataMap[id];
		},
		_particleData[dataName]);
}

void OptionalParticleData::setData(const std::string& dataName, size_t id,
								   OptionalParticleDataTypes::SupportedTypesVariant value) {
	std::visit(
		[=](auto& dataMap) {
			using map_type = std::remove_reference_t<decltype(dataMap)>;
			// Set value of map with key id.
			dataMap[id] = std::get<typename map_type::mapped_type>(value);
		},
		_particleData[dataName]);
}

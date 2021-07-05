/**
 * @file OptionalParticleData.h
 * @author seckler
 * @date 13.08.20
 */

#pragma once

#include <any>
#include <unordered_map>
#include <variant>

/**
 * Helper class to generate valid types.
 * @tparam SupportedTypes
 */
template <typename... SupportedTypes>
struct OptionalParticleDataTypesWithTemplate {
	/**
	 * Internal type of the map from a particle ID to the data.
	 * @note Here we choose an unordered_map for fast access.
	 */
	template <typename Type>
	using ParticleIdToDataMap = std::unordered_map<size_t, Type>;

	/**
	 * Map From name of the data to ParticleIdToDataMap.
	 */
	using NameToParticleDataMap = std::unordered_map<std::string, std::variant<ParticleIdToDataMap<SupportedTypes>...>>;

	/**
	 * Variant for the basic types.
	 */
	using SupportedTypesVariant = std::variant<SupportedTypes...>;

	/**
	 * Map From name of the data to default value.
	 */
	using NameToDefaultValueMap = std::unordered_map<std::string, SupportedTypesVariant>;
};

// ADD a type here if you want more supported types:
using OptionalParticleDataTypes = OptionalParticleDataTypesWithTemplate<double>;

/**
 * Class to store optional particle data.
 *
 * This class assumes that before any particles are added, all maps are created.
 */
class OptionalParticleData {
public:
	/**
	 * Creates a map and associates a default value to it.
	 * @param dataName
	 * @param defaultValue
	 */
	void createMap(const std::string& dataName, OptionalParticleDataTypes::SupportedTypesVariant defaultValue);

	/**
	 * Creates default values for all present maps.
	 * @param id
	 * @note: do not use from plugin!
	 */
	void createParticleDataWithDefaultValues(size_t id);

	/**
	 * Removes the data of one particle.
	 * @param id
	 * @note: do not use from plugin!
	 */
	void removeParticleData(size_t id);

	/**
	 * Get the data of the data associated with dataName for one particle.
	 * @param dataName The name of the data.
	 * @param id The id of the particle.
	 * @return Returns a variant of all supported values.
	 */
	OptionalParticleDataTypes::SupportedTypesVariant getData(const std::string& dataName, size_t id);

	/**
	 * Set the data associated with dataName for one particle.
	 * @param dataName The name of the data.
	 * @param id The id of the particle.
	 * @param value The value to which the data should be set.
	 */
	void setData(const std::string& dataName, size_t id, OptionalParticleDataTypes::SupportedTypesVariant value);

private:
	OptionalParticleDataTypes::NameToParticleDataMap _particleData;
	OptionalParticleDataTypes::NameToDefaultValueMap _defaultValues;
};

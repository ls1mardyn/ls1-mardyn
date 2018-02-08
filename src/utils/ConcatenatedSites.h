/**
 * ConcatenatedSites.h
 * Micha Mueller
 */

#ifndef CONCATENATEDSITES_H
#define CONCATENATEDSITES_H

#include "AlignedArrayTriplet.h"
#include "utils/mardyn_assert.h"
#include "../particleContainer/adapter/vectorization/SIMD_TYPES.h"
#include <array>

namespace ConcSites {
	/**
	 * \brief What coordinate would you like to have?
	 */
	enum class CoordinateType {
		X, Y, Z
	};

	/**
	 * \brief Specify which of the 4 data-categories is needed
	 */
	enum class SiteType {
		LJC, CHARGE, DIPOLE, QUADRUPOLE
	};
} /*namespace ConcSites*/

/**
 * \brief	Class to manage the storage of ljc-, charge-, dipole- and quadrupole-data in one single AlignedArrayTriplet
 * \details
 * \tparam 	T The type which is to be stored in the Arrays
 * \author	Micha Mueller
 */
template<typename T>
class ConcatenatedSites {
public:

	/*
	 * \brief Constructor
	 */
	ConcatenatedSites(size_t ljc_num = 0, size_t charges_num = 0,
			size_t dipoles_num = 0, size_t quadrupoles_num = 0) {
		resize(ljc_num, charges_num, dipoles_num, quadrupoles_num);
	}

	/**
	 * \brief 	Get a Pointer to the beginning of the specified data
	 * \details
	 * \tparam st		Specify which of the four site-types you would like to access. See enum ConcSites::SiteType for possible values
	 * \tparam coord	Indicates which of the 3 coordinates one needs. Has to be a value as specified in enum ConcSites::CoordinateType
	 * \return	Pointer to the first value of the data, as indicated by the parameters
	 */
	vcp_inline T* getBeginPointer (ConcSites::SiteType st, ConcSites::CoordinateType coord) {
		T* returnPointer = nullptr;

		switch (coord) {
		case ConcSites::CoordinateType::X:
			returnPointer = _data.xBegin();
			break;
		case ConcSites::CoordinateType::Y:
			returnPointer = _data.yBegin();
			break;
		case ConcSites::CoordinateType::Z:
			returnPointer = _data.zBegin();
			break;
		}

		mardyn_assert(returnPointer != nullptr);

		switch (st) {
		case ConcSites::SiteType::QUADRUPOLE:
			returnPointer += AlignedArray<T>::_round_up(_dipoles_num); // fallthrough
			/* no break */
		case ConcSites::SiteType::DIPOLE:
			returnPointer += AlignedArray<T>::_round_up(_charges_num); // fallthrough
			/* no break */
		case ConcSites::SiteType::CHARGE:
			returnPointer += AlignedArray<T>::_round_up(_ljc_num);
			/* no break */
		case ConcSites::SiteType::LJC:
			/* no break */ ; /* ; needed to compile here */
		}

		return returnPointer;
	}

	vcp_inline const T* getBeginPointer (ConcSites::SiteType st, ConcSites::CoordinateType coord) const {
		const T* returnPointer = nullptr;
		size_t offset = 0;

		switch (st) {
		case ConcSites::SiteType::QUADRUPOLE:
			offset += AlignedArray<T>::_round_up(_dipoles_num); // fallthrough
			/* no break */
		case ConcSites::SiteType::DIPOLE:
			offset += AlignedArray<T>::_round_up(_charges_num); // fallthrough
			/* no break */
		case ConcSites::SiteType::CHARGE:
			offset += AlignedArray<T>::_round_up(_ljc_num);
			/* no break */
		case ConcSites::SiteType::LJC:
			/* no break */ ; /* ; needed to compile here */
		}

		switch (coord) {
		case ConcSites::CoordinateType::X:
			returnPointer = _data.xBegin() + offset;
			break;
		case ConcSites::CoordinateType::Y:
			returnPointer = _data.yBegin() + offset;
			break;
		case ConcSites::CoordinateType::Z:
			returnPointer = _data.zBegin() + offset;
			break;
		}

		mardyn_assert(returnPointer != nullptr);

		return returnPointer;
	}

	/**
	 * \brief	Get the value triplet X,Y,Z of ConcSites::SiteType st at position index
	 */
	vcp_inline std::array<T, 3> getTriplet(ConcSites::SiteType st, size_t index) const {
		std::array<T, 3> retArray;
		retArray[0] = getBeginPointer(st, ConcSites::CoordinateType::X)[index];
		retArray[1] = getBeginPointer(st, ConcSites::CoordinateType::Y)[index];
		retArray[2] = getBeginPointer(st, ConcSites::CoordinateType::Z)[index];
		return retArray;
	}

	/**
	 * \brief	Set the value triplet X,Y,Z of ConcSites::SiteType st at position index to given values
	 */
	vcp_inline void setTriplet(std::array<T, 3> values, ConcSites::SiteType st, size_t index) {
		getBeginPointer(st, ConcSites::CoordinateType::X)[index] = values[0];
		getBeginPointer(st, ConcSites::CoordinateType::Y)[index] = values[1];
		getBeginPointer(st, ConcSites::CoordinateType::Z)[index] = values[2];
	}

	/**
	 * \brief 	Resize the ConcatenatedSites to have enough space for the given number of elements
	 * \tparam ljc_num
	 * \tparam charges
	 * \tparam dipoles_num
	 * \tparam quadrupoles_num
	 */
	void resize(size_t ljc_num, size_t charges_num, size_t dipoles_num, size_t quadrupoles_num) {
		_ljc_num = ljc_num;
		_charges_num = charges_num;
		_dipoles_num = dipoles_num;
		_quadrupoles_num = quadrupoles_num;

		size_t num_centers =
				  AlignedArray<T>::_round_up(_ljc_num)
				+ AlignedArray<T>::_round_up(_charges_num)
				+ AlignedArray<T>::_round_up(_dipoles_num)
				+ AlignedArray<T>::_round_up(_quadrupoles_num);

		_data.resize_zero_shrink(num_centers);
		setPaddingToZero(_data);
	}

	/**
	 * \brief	Get the size of currently occupied memory
	 * \return	Number of allocated bytes
	 */
	size_t get_dynamic_memory() const {	return _data.get_dynamic_memory(); }

private:

	AlignedArrayTriplet<T> _data;

// Is there a better solution than keeping them duplicated here and in CellDataSoA.h?
// But need them unless getBeginPointer shall be moved to CellDataSoA.h
// Sizes are duplicated in CellDataSoA, but there's no pretty way around this.
	size_t _ljc_num;
	size_t _charges_num;
	size_t _dipoles_num;
	size_t _quadrupoles_num;

	/*
	 * \brief	Set the unused memory to zero
	 */
	void setPaddingToZero(AlignedArray<T>& t) const {
		size_t ljc_size = t._round_up(_ljc_num);
		size_t charges_size = t._round_up(_charges_num);
		size_t dipoles_size = t._round_up(_dipoles_num);

		t.zero(_ljc_num);
		t.zero(ljc_size + _charges_num);
		t.zero(ljc_size + charges_size + _dipoles_num);
		t.zero(ljc_size + charges_size + dipoles_size + _quadrupoles_num);
	}

};

#endif /* CONCATENATEDSITES_H */

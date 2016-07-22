/*
 * Vector3.h
 *
 *  Created on: Nov 21, 2014
 *      Author: tchipevn
 */

#ifndef VECTOR3_H_
#define VECTOR3_H_

#include <cmath>
#include <iostream>

namespace bhfmm {

template<typename type>
class Vector3;

/** Global operators first */

template<typename type>
std::ostream& operator<<(std::ostream& stream, const Vector3<type>& v) {

	stream << "[" << v._content[0] << "; " << v._content[1] << "; " << v._content[2] << "]";
	return stream;
}

template<typename type>
Vector3<type> operator*(double scalar, const Vector3<type>& v) {
	return v * scalar;
}

template<typename type>
class Vector3 {
public:
	Vector3() {
	}

	Vector3(type arg) {
		_content[0] = arg;
		_content[1] = arg;
		_content[2] = arg;
	}

	Vector3(type arg0, type arg1, type arg2) {
		_content[0] = arg0;
		_content[1] = arg1;
		_content[2] = arg2;
	}

	Vector3(type args[3]) {
		_content[0] = args[0];
		_content[1] = args[1];
		_content[2] = args[2];
	}

	Vector3(const Vector3& other) {
		_content[0] = other._content[0];
		_content[1] = other._content[1];
		_content[2] = other._content[2];
	}

	Vector3 operator+(const Vector3& rhs) const {
		type result[3];
		result[0] = _content[0] + rhs._content[0];
		result[1] = _content[1] + rhs._content[1];
		result[2] = _content[2] + rhs._content[2];
		return Vector3(result);
	}

	void operator+=(const Vector3& rhs) {
		_content[0] += rhs._content[0];
		_content[1] += rhs._content[1];
		_content[2] += rhs._content[2];
	}

	Vector3 operator-(const Vector3& rhs) const {
		type result[3];
		result[0] = _content[0] - rhs._content[0];
		result[1] = _content[1] - rhs._content[1];
		result[2] = _content[2] - rhs._content[2];
		return Vector3(result);
	}

	Vector3 operator*(double scalar) const {
		type result[3];
		result[0] = _content[0] * scalar;
		result[1] = _content[1] * scalar;
		result[2] = _content[2] * scalar;
		return Vector3(result);
	}

	void operator*=(double scalar) {
			type result[3];
			_content[0] = _content[0] * scalar;
			_content[1] = _content[1] * scalar;
			_content[2] = _content[2] * scalar;
	}
	type MaxNorm() const {
		type norm;
		norm = std::max(abs(_content[0]), abs(_content[1]));
		norm = std::max(norm, abs(_content[2]));
		return norm;
	}

	type L2NormSquare() const {
		return _content[0] * _content[0] + _content[1] * _content[1] + _content[2] * _content[2];
	}

	double L2Norm() const {
		return sqrt(L2NormSquare());
	}

	Vector3& operator=(const Vector3& rhs) {
		if (this != &rhs) {
			_content[0] = rhs._content[0];
			_content[1] = rhs._content[1];
			_content[2] = rhs._content[2];
		}
		return *this;
	}

	Vector3& operator=(type rhs) {
		_content[0] = rhs;
		_content[1] = rhs;
		_content[2] = rhs;
		return *this;
	}

	type& operator[](int i) {
		return _content[i];
	}

	type operator[](int i) const {
		return _content[i];
	}

	bool operator==(const Vector3& rhs) const {
		if (_content[0] != rhs._content[0])
			return false;
		if (_content[1] != rhs._content[1])
			return false;
		if (_content[2] != rhs._content[2])
			return false;
		return true;
	}

	~Vector3() {
	}

	const type * data() const {
		return _content;
	}

	friend std::ostream& operator<<<type>(std::ostream& stream, const Vector3& v);
//	friend Vector3 operator*<type>(double scalar, const Vector3& v);

	/* used in std::map, returns true if lhs < rhs */
	class compare {
	public:
		bool operator()(const Vector3& lhs, const Vector3& rhs) const {
			return lhs[0] < rhs[0]
                || ( lhs[0] == rhs[0] && ( lhs[1] < rhs[1]
                || ( lhs[1] == rhs[1] && lhs[2] < rhs[2])));
		}
	};

private:
	type _content[3];
};

} /* namespace bhfmm */

#endif /* VECTOR3_H_ */

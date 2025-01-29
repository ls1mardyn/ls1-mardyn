/*
 * HaloBufferNoOverlap.h
 *
 *  Created on: Apr 23, 2016
 *      Author: obi
 */

#ifndef HALOBUFFERNOOVERHEAD_H_
#define HALOBUFFERNOOVERHEAD_H_

template <class T> class HaloBufferNoOverlap {
public:
	HaloBufferNoOverlap(int xHaloSize, int yHaloSize, int zHaloSize);
	virtual ~HaloBufferNoOverlap() = default;

	auto & getFrontBuffer(){
		return _frontBuffer;
	}
	auto & getBackBuffer(){
		return _backBuffer;
	}
	auto & getTopBuffer(){
		return _topBuffer;
	}
	auto & getBottomBuffer(){
		return _bottomBuffer;
	}
	auto & getLeftBuffer(){
		return _leftBuffer;
	}
	auto & getRightBuffer(){
		return _rightBuffer;
	}
	int  getXSize(){
		return _xHaloSize;
	}
	int  getYSize(){
		return _yHaloSize;
	}
	int  getZSize(){
		return _zHaloSize;
	}

	void clear();
private:
	// arrays for MPI halo transfer (send)
	std::vector<T> _leftBuffer, _rightBuffer, _topBuffer, _bottomBuffer, _frontBuffer, _backBuffer;
	int _xHaloSize,_yHaloSize,_zHaloSize;
};

#include <algorithm>


template <class T>
HaloBufferNoOverlap<T>::HaloBufferNoOverlap(int xHaloSize, int yHaloSize, int zHaloSize) {
	_xHaloSize = xHaloSize;
	_yHaloSize = yHaloSize;
	_zHaloSize = zHaloSize;

	_leftBuffer = std::vector<T>(_xHaloSize);
	_rightBuffer = std::vector<T>(_xHaloSize);

	_bottomBuffer = std::vector<T>(_yHaloSize);
	_topBuffer = std::vector<T>(_yHaloSize);

	_backBuffer = std::vector<T>(_zHaloSize);
	_frontBuffer = std::vector<T>(_zHaloSize);
}

template <class T>
void HaloBufferNoOverlap<T>::clear(){
	std::fill(_leftBuffer.begin(), _leftBuffer.end(), 0.0);
	std::fill(_rightBuffer.begin(), _rightBuffer.end(), 0.0);
	std::fill(_frontBuffer.begin(), _frontBuffer.end(), 0.0);
	std::fill(_backBuffer.begin(), _backBuffer.end(), 0.0);
	std::fill(_topBuffer.begin(), _topBuffer.end(), 0.0);
	std::fill(_bottomBuffer.begin(), _bottomBuffer.end(), 0.0);
}
#endif /* HALOBUFFER_H_ */

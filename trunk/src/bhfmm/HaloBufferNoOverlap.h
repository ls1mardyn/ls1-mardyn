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
	virtual ~HaloBufferNoOverlap();
	void initOverlap();
	//void initNoOverlap(int xHaloSize, int yHaloSize, int zHaloSize);
	T * getFrontBuffer(){
		return _frontBuffer;
	}
	T * getBackBuffer(){
		return _backBuffer;
	}
	T * getTopBuffer(){
		return _topBuffer;
	}
	T * getBottomBuffer(){
		return _bottomBuffer;
	}
	T * getLeftBuffer(){
		return _leftBuffer;
	}
	T * getRightBuffer(){
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
	T* _leftBuffer, * _rightBuffer, * _topBuffer, * _bottomBuffer, * _frontBuffer, * _backBuffer; //arrays for MPI halo transfer (send)
	int _xHaloSize,_yHaloSize,_zHaloSize;


};

#include <algorithm>


template <class T>
HaloBufferNoOverlap<T>::HaloBufferNoOverlap(int xHaloSize, int yHaloSize, int zHaloSize) {
	_xHaloSize = xHaloSize;
	_yHaloSize = yHaloSize;
	_zHaloSize = zHaloSize;

	_leftBuffer = new T[_xHaloSize];
	_rightBuffer = new T[_xHaloSize];

	_bottomBuffer = new T[_yHaloSize];
	_topBuffer = new T[_yHaloSize];

	_backBuffer = new T[_zHaloSize];
	_frontBuffer = new T[_zHaloSize];
}

template <class T>
HaloBufferNoOverlap<T>::~HaloBufferNoOverlap() {
	delete[] _leftBuffer;
	delete[] _rightBuffer;
	delete[] _bottomBuffer;
	delete[] _topBuffer;
	delete[] _backBuffer;
	delete[] _frontBuffer;
}

template <class T>
void HaloBufferNoOverlap<T>::clear(){
	std::fill(_leftBuffer, _leftBuffer + _xHaloSize , 0.0);
	std::fill(_rightBuffer, _rightBuffer + _xHaloSize , 0.0);
	std::fill(_frontBuffer, _frontBuffer + _zHaloSize, 0.0);
	std::fill(_backBuffer, _backBuffer + _zHaloSize, 0.0);
	std::fill(_topBuffer, _topBuffer + _yHaloSize, 0.0);
	std::fill(_bottomBuffer, _bottomBuffer + _yHaloSize, 0.0);
}
#endif /* HALOBUFFER_H_ */

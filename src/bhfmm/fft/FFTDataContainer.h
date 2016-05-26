/*
 * FFTDataContainer.h
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */
#ifndef FFTDATA_H_
#define FFTDATA_H_

/**
 * Abstract class of the storage of FFT expansion's related data
 * (simple case: two matrices, real and imag part of the expansion in Fourier space)
 */
class FFTDataContainer {
public:
	//force the call to the correct destructor
	virtual ~FFTDataContainer() {
	}

	/**
	 * Should return a deep copy of the container
	 *
	 * @return FFTDataContainer*, deep copy of this object
	 */
	virtual FFTDataContainer* copyContainer() = 0;

};

#endif

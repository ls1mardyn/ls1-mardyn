//
// Created by alex on 06.05.24.
//
#ifndef MARDYN_AVERAGER_H
#define MARDYN_AVERAGER_H

#include <algorithm>
#include <ostream>
#include <iostream>
#include "utils/Logger.h"

/**
 * Averages a container type T over every provided sample of T.
 * */
template<typename T>
class Averager {
public:
	Averager() : _sum_data(), _step_count(0) {	}

	/**
	 * Sets the capture size to the size of the container data.
	 * */
	void setDataSize(const T& data){
		_sum_data.resize(data.size(), 0.0);
	}

	void reset() {
		std::fill(_sum_data.begin(), _sum_data.end(), 0);
		_step_count = 0;
	}

	/**
	 * Captures data and increases capture count by one.
	 * */
	void averageData(const T& data) {
		if(data.size() != _sum_data.size()){
			Log::global_log->warning() << "[TimeAveraging] data structure mismatch" << std::endl;
			setDataSize(data);
			_step_count = 0;
		}

		_step_count++;
		for(std::size_t i = 0; i < data.size(); i++){
			_sum_data[i] = _sum_data[i] + data[i];
		}

	}

	/**
	 * Gets an immutable reference to the sum of all captured instances.
	 * */
	const T& getSumData() const {
		return _sum_data;
	}

	/**
	 * Gets the current averaged quantity over all captured instances.
	 * */
	T getAveragedDataCopy() const {
		T copy = _sum_data;
		for(int i = 0; i < copy.size(); i++) {
			copy[i] = copy[i] / static_cast<double>(_step_count);
		}
		return copy;
	}

	/**
	 * Writes the current averaged quantity into the provided buffer.\n
	 * Must have correct size
	 * */
	 void getAveragedData(T& dst) {
		for(int i = 0; i < _sum_data.size(); i++) {
			dst[i] = _sum_data[i] / static_cast<double>(_step_count);
		}
	 }

	/**
	 * Gets the total amount of captured frames
	 * */
	int getStepCount(){
		return _step_count;
	}

private:
	//! @brief capture buffer of data
	T _sum_data;
	//! @brief count of captured frames
	int _step_count;
};
#endif //MARDYN_AVERAGER_H

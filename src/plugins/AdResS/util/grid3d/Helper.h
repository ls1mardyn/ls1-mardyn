/*
 * Created on Wed Apr 24 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#ifndef MARDYN_HELPER_H
#define MARDYN_HELPER_H

#include <array>
#include<cmath>

namespace Grid3D {
	//Set array type
	using IdxArray=std::array<int, 3>;
	using DataArray = std::array<double, 3>;

	//TODO: all these need iterators and that kind of stuff
	//Use vector only where necessary, stop using tuple its stupid
	class Element{
	private:
		struct Properties{
			int index;
			double x_width;//Not set
			double y_width;//Not set
			double z_width;//Not set
			double volume;//Not set
			IdxArray local_indeces;
		} properties;

	public:
		Properties& GetProperties(){
			return properties;
		}
	};

	class Node{
	public:
		DataArray& GetPosition(){
			return this->position;
		}
		void SetPosition(double x, double y, double z){
			this->position={x,y,z};
		}
	private:
		DataArray position;
	};
};

template<class T>
std::array<int, 3> MapGlobalToLocal(int idx, T vec){
	std::array<int, 3> local_indeces{ 0 };
	local_indeces[0]= std::fmod(idx, vec[0]);
	local_indeces[1]= std::fmod(idx/vec[0],vec[1]);
	local_indeces[2]= idx/(vec[0]*vec[1]);

	return local_indeces;
}

template<class T>
int MapLocalToGlobal(std::array<int, 3>& local_indeces, T vec){
	int global_index = local_indeces[0] + local_indeces[1]*vec[0] + local_indeces[2]*vec[0]*vec[1];

	return global_index;
}

template<class T>
int MapLocalToGlobal(int i, int j, int k, T vec){
	T idcs = {i,j,k};
	return MapLocalToGlobal(idcs, vec);
}

#endif //MARDYN_HELPER_H

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
	using i3 =std::array<int, 3>;
	using d3 = std::array<double, 3>;

	//TODO: all these need iterators and that kind of stuff
	//Use vector only where necessary, stop using tuple its stupid
	//TODO: do these belong here? If you are the user of grid class, where would you expect to see this?
	class Element{
	private:
		struct Properties{
			int index;
			i3 local_indeces;
		} properties;

	public:
		Properties& GetProperties(){
			return properties;
		}
	};

	class Node{
	public:
		d3& GetPosition(){
			return this->position;
		}
		void SetPosition(double x, double y, double z){
			this->position={x,y,z};
		}
	private:
		d3 position;
	};
};

template<class T>
std::array<int, 3> MapGlobalToLocal(int idx, T dims){
	std::array<int, 3> local_indeces{ 0 };
	local_indeces[0]= std::fmod(idx, dims[0]);
	local_indeces[1]= std::fmod(idx / dims[0], dims[1]);
	local_indeces[2]= idx/(dims[0] * dims[1]);

	return local_indeces;
}

template<class T>
int MapLocalToGlobal(std::array<int, 3>& local_indeces, T dims){
	int global_index = local_indeces[0] + local_indeces[1] * dims[0] + local_indeces[2] * dims[0] * dims[1];

	return global_index;
}

template<class T>
int MapLocalToGlobal(int i, int j, int k, T dims){
	T idcs = {i,j,k};
	return MapLocalToGlobal(idcs, dims);
}

#endif //MARDYN_HELPER_H

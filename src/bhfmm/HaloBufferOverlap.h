/*
 * HaloBufferOverlap.h
 *
 *  Created on: Apr 26, 2016
 *      Author: obi
 */

#ifndef HALOBUFFEROVERLAP_H_
#define HALOBUFFEROVERLAP_H_
#include "bhfmm/utils/Vector3.h"
#include <math.h>
#include <vector>
#include <memory>
#ifdef ENABLE_MPI

#include "mpi.h"

namespace bhfmm{
template <class T>
class HaloBufferOverlap {
public:
	HaloBufferOverlap(Vector3<int> areaHaloSize, Vector3<int> edgeHaloSize,
		int cornerHaloSize, MPI_Comm comm, std::vector<int>& areaNeighbours,std::vector<int>& edgeNeighbours,std::vector<int>& cornerNeighbours, bool isSend, bool doNT,
		int areaNumber = 6, int edgeNumber = 12, int cornerNumber = 8 , std::vector<std::vector<std::vector<int>>> allRanks = std::vector<std::vector<std::vector<int>>>(0), Vector3<int> numCellsOnGlobalLevel = Vector3<int>(1), bool fuseGlobalCommunication = false);
	virtual ~HaloBufferOverlap();
	void startCommunication();
	//communicate without persistent sends and receives
	void communicate(bool postProcessing);
	void communicateGlobalLevels(int globalLevel, int stopLevel = 1, bool backCommunication = false);

	void wait();
	int testIfFinished();
	std::vector<T *>& getAreaBuffers(){
		return _areaBuffers;
	}
	std::vector<T *>& getEdgeBuffers(){
		return _edgeBuffers;
	}
	std::vector<T *>& getCornerBuffers(){
		return _cornerBuffers;
	}

	void setNumberOfGlobalLevelsInBuffer(int number){
		_globalLevelsInBuffer = number;
	}

	void clear();
private:
	void communicateLevelGlobal(int level, int globalLevel, int offset, bool backCommunication);
	void initCommunicationDouble();
	void fillArraySizes(Vector3<int> areaSizes, Vector3<int> edgeSizes);
	std::vector<T *>  _areaBuffers,  _edgeBuffers,  _cornerBuffers; //arrays for MPI halo transfer (send)
	Vector3<int> _areaHaloSize,_edgeHaloSize;
	int _cornerHaloSize;
	//these arrays save for every halo element the specific size -> if neighbour rank order is changed these arrays have to be changed too but nothing else
	std::vector<int> _areaHaloSizes, _edgeHaloSizes; //corner Arrays are always of the same size!
	std::vector<int> _areaNeighbours, _edgeNeighbours,_cornerNeighbours;
	MPI_Request * _areaRequests, *_edgeRequests, *_cornerRequests;
	MPI_Comm _comm;
	bool _isSend;
	bool _doNT;
	bool _isGlobal;
	bool _importWholeGlobalRegion;
	std::vector<std::vector<std::vector<int>>> _allRanks;
	Vector3<int> _numCellsOnGlobalLevel;
	int _globalLevelsInBuffer;
	int _offsetFactor;
	//if this flag is on then instead of 216 only 26 communications are performed and every processor sends 8 cell values instead of 1
	bool _fuseGlobalCommunication;
};

#include <algorithm>


template <class T>
HaloBufferOverlap<T>::HaloBufferOverlap(Vector3<int> areaHaloSize, Vector3<int> edgeHaloSize,
		int cornerHaloSize, MPI_Comm comm, std::vector<int>& areaNeighbours,std::vector<int>& edgeNeighbours,std::vector<int>& cornerNeighbours, bool isSend, bool doNT, int areaNumber, int edgeNumber, int cornerNumber, std::vector<std::vector<std::vector<int>>> allRanks, Vector3<int> numCellsOnGlobalLevel, bool fuseGlobalCommunication):
_areaBuffers(areaNumber), _edgeBuffers(edgeNumber), _cornerBuffers(cornerNumber),  _areaHaloSize(areaHaloSize), _edgeHaloSize(edgeHaloSize), _areaNeighbours(areaNeighbours), _edgeNeighbours(edgeNeighbours), _cornerNeighbours(cornerNeighbours), _doNT(doNT), _allRanks(allRanks), _numCellsOnGlobalLevel(numCellsOnGlobalLevel), _fuseGlobalCommunication(fuseGlobalCommunication) {

	_cornerHaloSize = cornerHaloSize;
	if(edgeNumber == 0){
		_isGlobal = true;
		if(numCellsOnGlobalLevel[0] > 1 or numCellsOnGlobalLevel[1] > 1 or numCellsOnGlobalLevel[2] > 1){
			_importWholeGlobalRegion = 1;
		}
		else{
			_importWholeGlobalRegion = 0;

		}
	}
	else{
		_isGlobal = false;
	}

	if(areaNumber != 0){
		std::vector<MPI_Request> _areaRequests(_areaBuffers.size());
	}
	if(edgeNumber != 0){
		std::vector<MPI_Request> _edgeRequests(_edgeBuffers.size());
	}
	if(cornerNumber != 0){
		std::vector<MPI_Request> _cornerRequests(_cornerBuffers.size());
	}

	fillArraySizes(areaHaloSize,edgeHaloSize);

//	_cornerHaloSizes = new int[_cornerBuffers.size()];
	_comm = comm;
	if(areaNumber != 0){
		for(unsigned int i=0; i<_areaBuffers.size();i++){
			if(!_isGlobal){
				_areaBuffers[i] = new T[_areaHaloSizes[i]];
			}
			else{
				_areaBuffers[i] = new T[_areaHaloSizes[0]];
			}
		}
	}

	if(edgeNumber != 0){
		for(unsigned int i=0; i<_edgeBuffers.size();i++){
			_edgeBuffers[i] = new T[_edgeHaloSizes[i]];
		}
	}
	if(cornerNumber != 0){
		for(unsigned int i=0; i<_cornerBuffers.size();i++){
			_cornerBuffers[i] = new T[cornerHaloSize];
		}
	}
	_isSend = isSend;
	clear();
	//initCommunicationDouble();
	if(_isGlobal){
		if(_importWholeGlobalRegion and not _fuseGlobalCommunication){
			if(_doNT){
				_offsetFactor = 25;
				if(_numCellsOnGlobalLevel[1] == 2){
					_offsetFactor += 19 ;
				}
				if(_numCellsOnGlobalLevel[0] == 2){
					_offsetFactor += 5 - _numCellsOnGlobalLevel[1] + 1;
				}
				if(_numCellsOnGlobalLevel[2] == 2){
					_offsetFactor += 5 - _numCellsOnGlobalLevel[1] + 1 ;
				}
			}
			else{
				_offsetFactor = (_fuseGlobalCommunication)? 26 : 216;
			}
		}
		else{
			if(_doNT){
				_offsetFactor = (_fuseGlobalCommunication)? 7 : 25;
			}
			else{
				_offsetFactor = (_fuseGlobalCommunication)? 26 : 189;
			}
		}
	}
	else{
		_offsetFactor = 0;
	}
}

template <class T>
void HaloBufferOverlap<T>::fillArraySizes(Vector3<int> areaSizes, Vector3<int> edgeSizes){
	if(!_isGlobal){
		_areaHaloSizes.resize(_areaBuffers.size());
		for(unsigned int i = 0; i < _areaBuffers.size(); i++){
			_areaHaloSizes[i] = areaSizes[i/2];
		}
	}
	else{
		_areaHaloSizes.resize(1);
		_areaHaloSizes[0] = areaSizes[0];

	}
	if(!_isGlobal){
		_edgeHaloSizes.resize(_edgeBuffers.size());
		for(unsigned int i = 0; i < _edgeBuffers.size(); i++){
			_edgeHaloSizes[i] = edgeSizes[2-i/4];
		}
	}
	else{
//		_edgeHaloSizes = new int[_edgeBuffers.size()];
//		for(unsigned int i = 0; i < _edgeBuffers.size(); i++){
//			_edgeHaloSizes[i] = edgeSizes[0];
//		}
	}
}
template <class T>
HaloBufferOverlap<T>::~HaloBufferOverlap() {

	for(unsigned int i=0; i<_areaBuffers.size();i++){
//		MPI_Request_free(&_areaRequests[i]);
		delete[] (_areaBuffers[i]);
	}
	if(!_isGlobal){

		for(unsigned int i=0; i<_edgeBuffers.size();i++){
	//		MPI_Request_free(&_edgeRequests[i]);
			delete[] _edgeBuffers[i];
		}
		for(unsigned int i=0; i<_cornerBuffers.size();i++){
	//		MPI_Request_free(&_cornerRequests[i]);
			delete[] _cornerBuffers[i];
		}
	}
}

template <class T>
void HaloBufferOverlap<T>::clear(){
	for(unsigned int i = 0; i < _areaBuffers.size(); i++){
		if(!_isGlobal){
			std::fill(_areaBuffers[i], _areaBuffers[i] + _areaHaloSizes[i] , 0.0);
		}
		else{
			std::fill(_areaBuffers[i], _areaBuffers[i] + _areaHaloSizes[0] , 0.0);
		}
	}
	for(unsigned int i = 0; i < _edgeBuffers.size(); i++){
		std::fill(_edgeBuffers[i], _edgeBuffers[i] + _edgeHaloSizes[i] , 0.0);
	}
	for(unsigned int i = 0; i < _cornerBuffers.size(); i++){
		std::fill(_cornerBuffers[i], _cornerBuffers[i] + _cornerHaloSize , 0.0);
	}
}

//we assume here that the neighbour arrays are sorted in the way that sites and their opposite sites are always alternating in the array
template <class T>
void HaloBufferOverlap<T>::initCommunicationDouble(){
	for (unsigned int i = 0; i < _areaBuffers.size(); i++){
		if(_isSend){
			MPI_Rsend_init(_areaBuffers[i], _areaHaloSizes[i], MPI_DOUBLE, _areaNeighbours[i], i + 42, _comm, &_areaRequests[i]);
			//MPI_Rsend_init(_areaBuffers[i], _areaHaloSize, MPI_DOUBLE, _areaNeighbours[i], i + 42, _comm, &_areaRequests[i]);

		}
		else{
			//adjusts that the tag of receive corresponds to send
			int indexShift = (i%2 == 0)? +1: -1;
			MPI_Recv_init(_areaBuffers[i], _areaHaloSizes[i], MPI_DOUBLE, _areaNeighbours[i], i + 42 + indexShift, _comm, &_areaRequests[i]);
		}
	}
	for (unsigned int i = 0; i < _edgeBuffers.size(); i++){
		if(_isSend){
			MPI_Rsend_init(_edgeBuffers[i], _edgeHaloSizes[i], MPI_DOUBLE, _edgeNeighbours[i], i + 42, _comm, &_edgeRequests[i]);
			//MPI_Rsend_init(_edgeBuffers[i], _edgeHaloSize, MPI_DOUBLE, _edgeNeighbours[i], i + 42, _comm, &_edgeRequests[i]);

		}
		else{
			int indexShift = (i%2 == 0)? +1: -1;
			MPI_Recv_init(_edgeBuffers[i], _edgeHaloSizes[i], MPI_DOUBLE, _edgeNeighbours[i], i + 42 + indexShift, _comm, &_edgeRequests[i]);
		}
	}
	for (unsigned int i = 0; i < _cornerBuffers.size(); i++){
		if(_isSend){
			MPI_Rsend_init(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42, _comm, &_cornerRequests[i]);
		//	MPI_Rsend_init(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42, _comm, &_cornerRequests[i]);
		}
		else{
			int indexShift = (i%2 == 0)? +1: -1;
			MPI_Recv_init(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42 + indexShift, _comm, &_cornerRequests[i]);
		}
	}
}

template <class T>
void HaloBufferOverlap<T>::communicate(bool postProcessing){
	int requestIndex = 0;
	for (unsigned int i = 0; i < _areaBuffers.size(); i++){
		if(_isSend){
			if(_doNT){
				if(postProcessing){
					if(i == 0 or i == 4){
						continue;
					}
				}
				else{
					if(i == 1 or i == 5){
						continue;
					}
				}
			}
			if(postProcessing)
				MPI_Isend(_areaBuffers[i], _areaHaloSizes[i], MPI_DOUBLE, _areaNeighbours[i], i + 42, _comm, &_areaRequests[requestIndex]);
			else
				MPI_Irsend(_areaBuffers[i], _areaHaloSizes[i], MPI_DOUBLE, _areaNeighbours[i], i + 42, _comm, &_areaRequests[requestIndex]);

			requestIndex++;
			//MPI_Rsend_init(_areaBuffers[i], _areaHaloSize, MPI_DOUBLE, _areaNeighbours[i], i + 42, _comm, &_areaRequests[i]);

		}
		else{
			if(_doNT){
				if(postProcessing){
					if(i == 1 or i == 5){
						continue;
					}
				}
				else{
					if(i == 0 or i == 4){
						continue;
					}
				}
			}
			//adjusts that the tag of receive corresponds to send
			int indexShift = (i%2 == 0)? +1: -1;
			MPI_Irecv(_areaBuffers[i], _areaHaloSizes[i], MPI_DOUBLE, _areaNeighbours[i], i + 42 + indexShift, _comm, &_areaRequests[requestIndex]);
			requestIndex++;
		}
	}
	requestIndex = 0;
	for (unsigned int i = 0; i < _edgeBuffers.size(); i++){
		if(_isSend){
			if(_doNT){
				if(not(postProcessing)){
					if(not(i == 4 or i == 6)){
						continue;
					}
				}
				else{
					if(not(i == 5 or i == 7)){
						continue;
					}
				}
			}
			if(postProcessing)
				MPI_Isend(_edgeBuffers[i], _edgeHaloSizes[i], MPI_DOUBLE, _edgeNeighbours[i], i + 42, _comm, &_edgeRequests[requestIndex]);
			else
				MPI_Irsend(_edgeBuffers[i], _edgeHaloSizes[i], MPI_DOUBLE, _edgeNeighbours[i], i + 42, _comm, &_edgeRequests[requestIndex]);
			requestIndex++;
			//MPI_Rsend_init(_edgeBuffers[i], _edgeHaloSize, MPI_DOUBLE, _edgeNeighbours[i], i + 42, _comm, &_edgeRequests[i]);

		}
		else{
			if(_doNT){
				if(not(postProcessing)){
					if(not(i == 5 or i == 7)){
						continue;
					}
				}
				else{
					if(not(i == 4 or i == 6)){
						continue;
					}
				}
			}
			int indexShift = (i%2 == 0)? +1: -1;
			MPI_Irecv(_edgeBuffers[i], _edgeHaloSizes[i], MPI_DOUBLE, _edgeNeighbours[i], i + 42 + indexShift, _comm, &_edgeRequests[requestIndex]);
			requestIndex++;
		}
	}
	requestIndex = 0;
	if(not(_doNT)){
		for (unsigned int i = 0; i < _cornerBuffers.size(); i++){
			if(_isSend){
				MPI_Irsend(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42, _comm, &_cornerRequests[requestIndex]);
				requestIndex++;
			//	MPI_Rsend_init(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42, _comm, &_cornerRequests[i]);
			}
			else{
				int indexShift = (i%2 == 0)? +1: -1;
				MPI_Irecv(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42 + indexShift, _comm, &_cornerRequests[requestIndex]);
				requestIndex++;
			}
		}
	}

}

template <class T>
void HaloBufferOverlap<T>::communicateGlobalLevels(int globalLevel, int stopLevel, bool backCommunication){
	int minimumLevel = (_fuseGlobalCommunication)? 2:1;
	stopLevel = (stopLevel < minimumLevel)? minimumLevel: stopLevel;
	for(int l = globalLevel; l >= stopLevel ; l--){
		int offset;
		if(_doNT and not _fuseGlobalCommunication){
			offset = (globalLevel == l) ? 0 : _offsetFactor + (globalLevel - l - 1) * 25;
		}
		else{
			offset = _offsetFactor * (globalLevel - l);
		}
		communicateLevelGlobal(l,globalLevel,offset, backCommunication);
	}
}

template <class T>
void HaloBufferOverlap<T>::communicateLevelGlobal(int level, int globalLevel, int offset, bool backCommunication){

	int stride = pow(2,globalLevel - level);
	int myRank;
	int coords[3];
	int coordsFloored[3];
	int coordsLevel[3];
	int coordsRemainder[3];
	MPI_Comm_rank(_comm,&myRank);
	int indexPosition = 0;
	int rank;
	MPI_Cart_coords(_comm, myRank, 3, coords);
	for(int d = 0; d < 3; d++){
		coordsFloored[d] = ((coords[d] * _numCellsOnGlobalLevel[d]) / (2 * stride)) * 2 * stride;
		coordsRemainder[d] = (coords[d] * _numCellsOnGlobalLevel[d]) % (stride);
		coordsLevel[d] = ((coords[d] * _numCellsOnGlobalLevel[d]) / (stride));

	}
	int start, end;
	if(_fuseGlobalCommunication){
		start = -1;
		end = 1;
	}
	else{
		start = -2;
		end = 3;
	}
	for(int x = start; x <= end; x++ ){
		for(int y = start; y <= end; y++){
			for(int z = start; z <= end; z++){
				bool condition;
				if(_doNT){ //allow only cells in plate or tower
					int cellsX = _numCellsOnGlobalLevel[0];
					int cellsY = _numCellsOnGlobalLevel[1];
					int cellsZ = _numCellsOnGlobalLevel[2];

					if((_isSend and !backCommunication) or (backCommunication and !_isSend)){ //reversed plate in this case
						if(!_fuseGlobalCommunication){
							condition = ( floor((x * stride + coordsFloored[0] + coordsRemainder[0])/(1.0 * cellsX)) == coords[0] and floor((z * stride + coordsFloored[2] + coordsRemainder[2])/(1.0 * cellsZ)) == coords[2]) //tower
										or (((y < 2 and y >= 0 and cellsY == 2 and level == globalLevel) or floor((y * stride + coordsFloored[1] + coordsRemainder[1])/(1.0 * cellsY)) == coords[1]) and x < 2 and not (x >= 0 and z >= 2)); //plate (reversed for send)
						}
						else{
							condition = (x == 0 and z == 0) // tower
										or (y == 0 and x < 1 and not ( x==0 and z == 1)); //plate (reversed for send)
						}
					}
					else{
						if(!_fuseGlobalCommunication){

						condition = ( floor((x * stride + coordsFloored[0] + coordsRemainder[0])/(1.0 * cellsX)) == coords[0] and floor((z * stride + coordsFloored[2] + coordsRemainder[2])/(1.0 * cellsZ)) == coords[2]) //tower
										or (((y < 2 and y >= 0 and cellsY == 2 and level == globalLevel) or floor((y * stride + coordsFloored[1] + coordsRemainder[1])/(1.0 * cellsY)) == coords[1]) and x >= 0 and not (x < 2 and z < 0)); //plate
						}
						else{
							condition = (x == 0 and z == 0) // tower
										or (y == 0 and x >= 0 and not ( x < 1 and z == -1)); //plate
						}
					}
//					condition = condition and (x >= 2 or x < 0 or z >= 2 or z < 0 or y >= 2 or y < 0); // do not send within parent cell

				}
				else{
					if(!_fuseGlobalCommunication){//import all 189 (if ranks are a power of 8) or 216 (if multiple cells on global level) values
						condition = _importWholeGlobalRegion or abs(x +(coordsFloored[0] - coords[0])/stride) >= 2 or abs(y +(coordsFloored[1] - coords[1])/stride) >= 2 or abs(z + (coordsFloored[2] - coords[2])/stride) >= 2;
					}
					else{ //imort all cells from 26 neighbours (each one sends 8 cells)
						condition = x != 0 or y != 0 or z != 0;
					}
				}
				if(condition){
					int xIndex, yIndex, zIndex;
					if(!_fuseGlobalCommunication){
						xIndex = (int) floor(((coordsFloored[0] + (x * stride + coordsRemainder[0]) * 1.0) / _numCellsOnGlobalLevel[0]) + _allRanks.size()) % _allRanks.size();
						yIndex = (int) floor(((coordsFloored[1] + (y * stride + coordsRemainder[1]) * 1.0) / _numCellsOnGlobalLevel[1]) + _allRanks[0].size()) % _allRanks[0].size();
						zIndex = (int) floor(((coordsFloored[2] + (z * stride + coordsRemainder[2]) * 1.0) / _numCellsOnGlobalLevel[2]) + _allRanks[0][0].size()) % _allRanks[0][0].size();
					}
					else{
						int xLocal = coordsLevel[0] % 2;
						int yLocal = coordsLevel[1] % 2;
						int zLocal = coordsLevel[2] % 2;
						xIndex = (int) floor(((coordsFloored[0] + ((2 * x + xLocal) * stride + coordsRemainder[0]) * 1.0) / _numCellsOnGlobalLevel[0]) + _allRanks.size()) % _allRanks.size();
						yIndex = (int) floor(((coordsFloored[1] + ((2 * y + yLocal) * stride + coordsRemainder[1]) * 1.0) / _numCellsOnGlobalLevel[1]) + _allRanks[0].size()) % _allRanks[0].size();
						zIndex = (int) floor(((coordsFloored[2] + ((2 * z + zLocal) * stride + coordsRemainder[2]) * 1.0) / _numCellsOnGlobalLevel[2]) + _allRanks[0][0].size()) % _allRanks[0][0].size();
					}
					rank = _allRanks[xIndex][yIndex][zIndex];

					const int xOffset = (_fuseGlobalCommunication)? 0 : abs(x * stride + coordsRemainder[0]) % _numCellsOnGlobalLevel[0];
					const int yOffset = (_fuseGlobalCommunication)? 0 : abs(y * stride + coordsRemainder[1]) % _numCellsOnGlobalLevel[1];
					const int zOffset = (_fuseGlobalCommunication)? 0 : abs(z * stride + coordsRemainder[2]) % _numCellsOnGlobalLevel[2];
//					if(myRank == 1){
////						std::cout << myRank << " "<<rank << " " << xIndex << " " << yIndex << " " << zIndex << " " << x << " " << y << " "<< z << " " << coordsFloored[0] << " " <<  coordsFloored[1] << " " << coordsFloored[2] << " " << xOffset << " " << yOffset << " " << zOffset << "\n";
//					}
					//number of cells the every processor is responsible for on current level
					int numCellsOnLevel = (level == globalLevel)? _numCellsOnGlobalLevel[0] * _numCellsOnGlobalLevel[1] * _numCellsOnGlobalLevel[2] : 1;
					if(_isSend){
						if(!backCommunication){
							MPI_Irsend(_areaBuffers[8*(globalLevel-level) + 4*zOffset + 2*yOffset + xOffset], _areaHaloSizes[0], MPI_DOUBLE, rank, 1000 + zOffset * 4 + yOffset * 2 + xOffset + 8 * (globalLevel - level), _comm, &_areaRequests[indexPosition + offset]);
						}
						else{
							if(_fuseGlobalCommunication){ //in the back communication only 1 (or equal to the number of cells one owns on global level) cell is send instead of 8
								MPI_Isend(_areaBuffers[indexPosition + offset], _areaHaloSizes[0] / 8 * numCellsOnLevel , MPI_DOUBLE, rank, 1000 + zOffset * 4 + yOffset * 2 + xOffset + 8 * (globalLevel - level), _comm, &_areaRequests[indexPosition + offset]);
							}
							else{
								MPI_Isend(_areaBuffers[indexPosition + offset], _areaHaloSizes[0], MPI_DOUBLE, rank, 1000 + zOffset * 4 + yOffset * 2 + xOffset + 8 * (globalLevel - level), _comm, &_areaRequests[indexPosition + offset]);
							}
						}
						indexPosition++;
					//	MPI_Rsend_init(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42, _comm, &_cornerRequests[i]);
					}
					else{
//						std::cout << indexPosition << "\n";
//						MPI_Barrier(_comm);
						if(!backCommunication){
							MPI_Irecv(_areaBuffers[indexPosition + offset], _areaHaloSizes[0], MPI_DOUBLE, rank, 1000 + zOffset * 4 + yOffset * 2 + xOffset + 8 * (globalLevel - level), _comm, &_areaRequests[indexPosition + offset]);
						}
						else{
							if(_fuseGlobalCommunication){ //in the back communication only 1 (or equal to the number of cells one owns on global level) cell is send instead if 8
								MPI_Irecv(_areaBuffers[indexPosition + offset], _areaHaloSizes[0] / 8 * numCellsOnLevel, MPI_DOUBLE, rank, 1000 + zOffset * 4 + yOffset * 2 + xOffset + 8 * (globalLevel - level), _comm, &_areaRequests[indexPosition + offset]);
							}
							else{
								MPI_Irecv(_areaBuffers[indexPosition + offset], _areaHaloSizes[0], MPI_DOUBLE, rank, 1000 + zOffset * 4 + yOffset * 2 + xOffset + 8 * (globalLevel - level), _comm, &_areaRequests[indexPosition + offset]);
							}
						}
						indexPosition++;
					}
				}
			}
		}
	}
//	if(_doNT)
//		std::cout << " indexposition: "<< indexPosition <<" " <<_isSend<< " "<< backCommunication << " offset "<< offset<< " offsetFactor: "<< _offsetFactor  << " rank: "<< myRank <<"\n";
	if((indexPosition != _offsetFactor and (not(_doNT) or _fuseGlobalCommunication or level == globalLevel)) or(indexPosition != 25 and _doNT and level != globalLevel and not _fuseGlobalCommunication)){
		std::cout << "Error offsetFactor is calculated wrong or too few sends!!! -> synchronization possibly broken!!! \n";
		std::cout << indexPosition << " " << _offsetFactor << " \n";
	}
//	MPI_Barrier(_comm);
//	std::cout << indexPosition << "\n";
}

template <class T>
void HaloBufferOverlap<T>::startCommunication(){
	//outdated!!!
	if(not(_doNT)){
		 MPI_Startall(_areaBuffers.size(), _areaRequests);
		 MPI_Startall(_edgeBuffers.size(), _edgeRequests);
		 MPI_Startall(_cornerBuffers.size(), _cornerRequests);
	}
	else{
		MPI_Startall(4, _areaRequests);
		MPI_Startall(2, _edgeRequests);
	}
//	 std::cout << _areaBuffers.size() << _edgeBuffers.size() << _cornerBuffers.size() <<"\n";
}

template <class T>
void HaloBufferOverlap<T>::wait(){
	//outdated!!!
	if(not(_doNT)){
		MPI_Status * areaStatusArray = new MPI_Status[_areaBuffers.size()];
		MPI_Waitall(_areaBuffers.size(),_areaRequests, areaStatusArray);

		MPI_Status * edgeStatusArray = new MPI_Status[_edgeBuffers.size()];
		MPI_Waitall(_edgeBuffers.size(),_edgeRequests, edgeStatusArray);

		MPI_Status * cornerStatusArray = new MPI_Status[_cornerBuffers.size()];
		MPI_Waitall(_cornerBuffers.size(),_cornerRequests, cornerStatusArray);
	}
	else{
		MPI_Status * areaStatusArray = new MPI_Status[4];
		MPI_Waitall(4,_areaRequests, areaStatusArray);

		MPI_Status * edgeStatusArray = new MPI_Status[2];
		MPI_Waitall(2,_edgeRequests, edgeStatusArray);
	}
}

template <class T>
int HaloBufferOverlap<T>::testIfFinished(){
	int areaFlag, edgeFlag, cornerFlag;
	if(not(_doNT)){
		if(!_isGlobal){

			std::vector<MPI_Status> areaStatusArray(_areaBuffers.size());
			MPI_Testall(_areaBuffers.size(),_areaRequests, &areaFlag, areaStatusArray.data());

			std::vector<MPI_Status> edgeStatusArray(_edgeBuffers.size());
			MPI_Testall(_edgeBuffers.size(),_edgeRequests, &edgeFlag, edgeStatusArray.data());

			std::vector<MPI_Status> cornerStatusArray(_cornerBuffers.size());
			MPI_Testall(_cornerBuffers.size(),_cornerRequests, &cornerFlag, cornerStatusArray.data());

			return areaFlag * edgeFlag * cornerFlag;
		}
		else{
//			std::cout << _areaBuffers.size() << "\n";
			if(_areaBuffers.size() == 0) return true;
			std::vector<MPI_Status> areaStatusArray(_areaBuffers.size());
			MPI_Testall(_areaBuffers.size(),_areaRequests, &areaFlag, areaStatusArray.data());
			return areaFlag;
		}
	}
	else{
		if(!_isGlobal){
			MPI_Status areaStatusArray[4];
			MPI_Testall(4,_areaRequests, &areaFlag, areaStatusArray);

			MPI_Status edgeStatusArray[2];
			MPI_Testall(2,_edgeRequests, &edgeFlag, edgeStatusArray);
			return areaFlag * edgeFlag;
		}
		else{
			int numRequests;
//			std::cout << "Testing! \n" ;
			if(!_fuseGlobalCommunication){
				numRequests = (_globalLevelsInBuffer == 1) ? _offsetFactor : _offsetFactor + (_globalLevelsInBuffer - 1) * 25;
			}
			else{
				numRequests = _globalLevelsInBuffer * _offsetFactor ;
			}
			if(numRequests == 0) return true;
			std::vector<MPI_Status> areaStatusArray(numRequests);
			MPI_Testall(numRequests,_areaRequests, &areaFlag, areaStatusArray.data());
//			MPI_Status status;
//
//			for(int i= 0; i<numberOfSends; i++){
//				MPI_Test(&_areaRequests[i],&areaFlag,&status);
//				if(!areaFlag){
////					std::cout << i <<"\n" ;
//					return areaFlag;
//				}
//			}

			return areaFlag;
		}
	}
}
}
#endif
#endif /* HALOBUFFEROVERLAP_H_ */

// Andreas Kirsch <kirschan@tum.de>

#include "LinkedCellsCUDAold.h"
#include "LinkedCells.h"
#include "molecules/potforce.h"
#include "handlerInterfaces/ParticlePairsHandler.h"
#include "Cell.h"
#include "molecules/Molecule.h"
#include "ensemble/GrandCanonical.h"
#include "io/OneCLJGenerator.h"
#include <sys/time.h>

#include <string>

#include <cuda_runtime.h>
#include "cutil_math.h"

struct CUDAException : public std::exception {
	const cudaError_t errorCode;
	const std::string errorSource;

	CUDAException( cudaError_t errorCode, const std::string &errorSource = "" )
	: errorCode( errorCode ), errorSource( errorSource ) {}

	~CUDAException() throw() {}

    /** Returns a C-style character string describing the general cause
     *  of the current error.  */
    virtual const char* what() const throw() {
    	return errorSource.c_str();
    }
};

#define CUDA_THROW_ON_ERROR( expr ) \
	do { \
		cudaError_t errorCode = (expr); \
		if( errorCode != cudaSuccess ) { \
			throw CUDAException( errorCode, #expr ); \
		} \
	} while( 0 )


extern __shared__ float4 pblock[];
__global__ void LJ( float4* pos,			//positions of the particles
					int* pairs,			//the cell pairs
					float4* force,
					int numberOfPairs, int maxCellSize,float cutOffRadius,float eps24,float sig2
					)
{
	int nt = 64;						//nt threads are responsible for each cell pair
	int n = gridDim.x * blockDim.x; //get_global_size(0); 		//size of the index space
  	int gti = (blockIdx.x * blockDim.x + threadIdx.x) / nt;		//index of the cell pair of this thread (global id divided by number of threads per pair)
  	int ti = threadIdx.x;			//index of this thread inside its workgroup
  	int mi= threadIdx.x % nt;		//index of the particle inside the first cell for which this thread is responsible
  										//if the cell contains more than nt particles, then mi is used as an offset

	if(gti < numberOfPairs) {
		int cell1=pairs[gti*5];				//startindex of first cell
		int cell2=pairs[gti*5+1];			//startindex of second cell
		int cellSize1=pairs[gti*5+2];		//size of first cell
		int cellSize2=pairs[gti*5+3];		//size of second cell
		//int offset=pairs[gti*5+4];
		float4 p;							//stores the current particle position
		float4 f;							//stores the current force of the particle
		float virial=0.;

		//if the cell contains more than nt particles, each thread has to be responsible for several particles
		//this loop iterates over them
		for(int i1=0; i1 < ceil((float)cellSize1/(float)nt)+1; i1++) {
			int index1=i1*nt+mi;		//index of current particle inside the first cell
			if(index1 < cellSize1) {	//check if index1 is inside the first cell
				p=pos[cell1+index1];
			}
			f.x=f.y=f.z=f.w=0.0f;
			for(int i2=0; i2 < ceil((float)cellSize2/(float)nt)+1; i2++) {
				int index2=i2*nt;
				if(index2+mi < cellSize2) {
					pblock[ti] = pos[cell2+index2+mi];
				}
				__syncthreads();
					if(index1 < cellSize1) {
						for(int m=0; m < nt; m++) {		 //iterate over all cached particles
							if(index2+m < cellSize2) {
								float4 p2=pblock[ti-mi+m];
								//calculate distance
								float4 d = p2 - p;
								d.w = d.x*d.x +d.y*d.y + d.z*d.z;
								if(d.w < cutOffRadius && d.w != 0) {	// pairs are inside the cut of radius and not the same particle
									//calculate the LJ Potential and sum it up:
									float invdr2=1.f/d.w;
  									float lj6=sig2*invdr2; lj6=lj6*lj6*lj6;
  									float lj12=lj6*lj6;
  									float lj12m6=lj12-lj6;
  									float u6=eps24*lj12m6+f.w;
  									float fac=eps24*(lj12+lj12m6)*invdr2;
  									//sum up the forces
  									f=f+fac*d;
  									f.w=u6;
  									virial+=f.x*d.x+f.y*d.y+f.z*d.z;  //needed for virial in
								}
							}
						}
					}
				__syncthreads();
			}
			if(index1 < cellSize1) {
			    int i=(cell1+index1)*2;
				force[i]=f;
				force[i+1].x=virial;
			}
		}
	}
}

LinkedCellsOpenCL::LinkedCellsOpenCL(double bBoxMin[3], double bBoxMax[3], double cutoffRadius, double LJCutoffRadius,
	     double tersoffCutoffRadius, double cellsInCutoffRadius,
	     ParticlePairsHandler* partPairsHandler)
			: LinkedCells(bBoxMin,bBoxMax,cutoffRadius,LJCutoffRadius,tersoffCutoffRadius,cellsInCutoffRadius,partPairsHandler) {
	LinkedCellsOpenCL::numberOfParticles=numberOfParticles;

	totalTime=0;
	numberOfParticles = 0;
}

LinkedCellsOpenCL::~LinkedCellsOpenCL() {
	cout << "LinkedCellsOpenCL::deconstruct" << endl;
}

void LinkedCellsOpenCL::initCUDA() {
	cout << "Initializing CUDA" << endl;

	try {
		int	deviceCount;
		CUDA_THROW_ON_ERROR( cudaGetDevice( &deviceCount ) );

		cout << deviceCount << " devices found" << endl;
		cout << "Using device 0" << endl;
		CUDA_THROW_ON_ERROR( cudaSetDevice( 0 ) );

		// TODO: output more information using cudaGetDeviceProps
	}
	catch( const CUDAException &cudaError ) {
		std::cerr << "ERROR: " << cudaError.errorSource << "(" << cudaError.errorCode << ")"
				<< std::endl;

		exit(EXIT_FAILURE);
	}
}

void LinkedCellsOpenCL::invertPairs() {
	for(int i=0; i < numberOfPairs;++i) {
		int swapId=pairs[i*5];
		int swapSize=pairs[i*5+2];
		pairs[i*5]=pairs[i*5+1];
		pairs[i*5+2]=pairs[i*5+3];
		pairs[i*5+1]=swapId;
		pairs[i*5+3]=swapSize;
	}
}

//traverse pairs an calculate LJ on GPU
void LinkedCellsOpenCL::traversePairs() {
	timeval start,end;
	int numberOfCells=_cells.size();
	cout << "LinkedCellsOpenCL::traversePairs()" << endl;
	if(numberOfParticles == 0) {
		initCUDA();
		cout << "LinkedCellsOpenCL::countPairs()" << endl;
		countPairs();
		cout << "LinkedCellsOpenCL::countPairs():" << numberOfPairs << "  numberOfCells " << numberOfCells<< endl;
	}
	numberOfParticles=countParticles();
	cout << "LinkedCellsOpenCL::countParticles in Cells " << numberOfParticles<< endl;
	m_positions = (float*) memalign(16, numberOfParticles * sizeof(float4));
	cout << "LinkedCellsOpenCL::traversePairs()" << endl;
	cout << "LinkedCellsOpenCL::traversePairsInit()" << endl;
	traversePairsInit();




	 int maxCellSize=getMaxCellSizeAndPositions();


	pairs = (int*) memalign(16, 2*numberOfPairs*5* sizeof(int));

	cout << "createPairs-" << endl;
	createPairs();
	cout << "createPairs+" << endl;
	for(int m=0; m < numberOfPairs; ++m) {
		//cout << "pair: cell1 " <<  pairs[m*5] <<" cell2 " <<  pairs[m*5+1]<<" cell1 size " <<  pairs[m*5+2]<<" cell2 size "  <<  pairs[m*5+3]<<" offset " <<  pairs[m*5+4] << endl;
	}

	cout << "numberOfPairs " << numberOfPairs << " resultSize " << resultSize << endl;

	forces = (float*)  memalign(16, 2*numberOfParticles  * sizeof(float4));


	cout << "LinkedCellsOpenCL::calculate distances" << endl;

	double dif;
	gettimeofday(&start,NULL);
	try {
		void *memPositions, *memPairs, *memForces;

		// TODO: move the device memory allocation and release into the init/deinit functions
		CUDA_THROW_ON_ERROR( cudaMalloc( &memPositions, numberOfParticles * sizeof(float4) ) );
		CUDA_THROW_ON_ERROR( cudaMalloc( &memPairs, numberOfPairs *10* sizeof(int) ) );
		CUDA_THROW_ON_ERROR( cudaMalloc( &memForces, 2*numberOfParticles * sizeof(float4) ) );

		CUDA_THROW_ON_ERROR( cudaMemcpy( memPositions, m_positions, numberOfParticles * sizeof(float4), cudaMemcpyHostToDevice ) );
		CUDA_THROW_ON_ERROR( cudaMemcpy( memPairs, pairs, numberOfPairs *10* sizeof(int), cudaMemcpyHostToDevice ) );
		CUDA_THROW_ON_ERROR( cudaMemcpy( memForces, forces, 2*numberOfParticles * sizeof(float4), cudaMemcpyHostToDevice ) );

		cout << "LinkedCellsOpenCL maxCellSize " << maxCellSize<< endl;

		int gridSize = (numberOfPairs-1) / max_threads + 1;

		cout << "LinkedCellsOpenCL::gridSize:" << gridSize << endl;

		size_t sharedMemorySize = max_threads * sizeof(float4);
		// TODO: refactor *64..
		LJ<<<gridSize*64, max_threads, sharedMemorySize>>>( (float4*) memPositions, (int*) memPairs, (float4*) memForces, numberOfPairs, maxCellSize,
				(float) cutoffRadiusSquare, (float) 24, (float) 1);

		CUDA_THROW_ON_ERROR( cudaMemcpy( forces, memForces, 2*numberOfParticles	* sizeof(float4), cudaMemcpyDeviceToHost ) );

		CUDA_THROW_ON_ERROR( cudaFree( memPairs ) );
		CUDA_THROW_ON_ERROR( cudaFree( memForces ) );

		invertPairs();
		float* tempForces=forces;

		forces = (float*)  memalign(16, 2*numberOfParticles  * sizeof(float4));

		void *memPairs2, *memForces2;
		CUDA_THROW_ON_ERROR( cudaMalloc( &memPairs2, numberOfPairs *10* sizeof(int) ) );
		CUDA_THROW_ON_ERROR( cudaMalloc( &memForces2, 2*numberOfParticles * sizeof(float4) ) );

		CUDA_THROW_ON_ERROR( cudaMemcpy( memPairs2, pairs, numberOfPairs *10* sizeof(int), cudaMemcpyHostToDevice ) );
		CUDA_THROW_ON_ERROR( cudaMemcpy( memForces2, forces, 2*numberOfParticles * sizeof(float4), cudaMemcpyHostToDevice ) );

		// TODO: refactor *64..
		LJ<<<gridSize*64, max_threads, sharedMemorySize>>>( (float4*) memPositions, (int*) memPairs2, (float4*) memForces2, numberOfPairs, maxCellSize,
				(float) cutoffRadiusSquare, (float) 24, (float) 1);

		CUDA_THROW_ON_ERROR( cudaMemcpy( forces, memForces2, 2*numberOfParticles * sizeof(float4), cudaMemcpyDeviceToHost ) );

		CUDA_THROW_ON_ERROR( cudaFree( memPairs2 ) );
		CUDA_THROW_ON_ERROR( cudaFree( memForces2 ) );

		//sum up the forces of the two steps:
		for(int m=0; m < numberOfParticles*8; ++m) {
				forces[m]+=tempForces[m];
		}
		free(tempForces);
		CUDA_THROW_ON_ERROR( cudaFree( memPositions ) );
	}
	catch( const CUDAException &cudaError ) {
		std::cerr << "ERROR in LinkedCellsOpenCL::traversePairs(): " << cudaError.errorSource << "(" << cudaError.errorCode << ")"
				<< std::endl;

		exit(EXIT_FAILURE);
	}

	gettimeofday(&end,NULL);
	dif=((end.tv_sec - start.tv_sec) * 1000000) +(end.tv_usec - start.tv_usec);
	dif= dif / 1000000;
	totalTime+=dif;
	cout << "calculated distances in  " << dif << " seconds total:  " << totalTime << endl;

	cout << "LinkedCellsOpenCL::calculate pairs" << endl;

	//copy the forces to the molecules:
	int index=0;
	std::list<Molecule*>::iterator molIter1;
	double uPot=0;
	float virial=0;
	for(unsigned i = 0; i < _cells.size(); i++ ){
					Cell& currentCell = _cells[i];
					for( molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++ ) {
						Molecule& molecule1 = **molIter1;
						// TODOAK: hack hack hack.. constness is cast away
						double* Fsite=(double*) molecule1.ljcenter_F(0);
						Fsite[0]=(double)forces[index*8];
						Fsite[1]=(double)forces[index*8+1];
						Fsite[2]=(double)forces[index*8+2];
						uPot=uPot+(double)forces[index*8+3];
						virial+=forces[index*8+4];
						index++;
					}
	}

	cout << "LinkedCellsOpenCL::finish" << endl;
	traversePairsFinish();

	free(pairs);
	free(distances);
	free(m_positions);
}

void LinkedCellsOpenCL::traversePairsInit() {

	LinkedCells::_particlePairsHandler->init();

	// loop over all cells
	vector<Cell>::iterator cellIter;
	std::list<Molecule*>::iterator molIter1;
	std::list<Molecule*>::iterator molIter2;
	for( cellIter = _cells.begin(); cellIter !=  _cells.end(); cellIter++ ) {
		for( molIter1 = cellIter->getParticlePointers().begin(); molIter1 != cellIter->getParticlePointers().end(); molIter1++ ) {
			double zero[] = {0,0,0};
			(*molIter1)->setF( zero );
			(*molIter1)->setM( zero );
		}
	}

	vector<unsigned long>::iterator cellIndexIter;
	vector<unsigned long>::iterator neighbourOffsetsIter;

	// sqare of the cutoff radius
	 cutoffRadiusSquare = _cutoffRadius * _cutoffRadius;
	 LJCutoffRadiusSquare = _LJCutoffRadius * _LJCutoffRadius;
	tersoffCutoffRadiusSquare = _tersoffCutoffRadius * _tersoffCutoffRadius;


	for( unsigned i = 0; i < _cells.size(); i++ )
	{
		Cell& currentCell = _cells[i];
		for( molIter1 = currentCell.getParticlePointers().begin();
				molIter1!=currentCell.getParticlePointers().end();
				molIter1++ )
		{
			Molecule& molecule1 = **molIter1;
			molecule1.clearTersoffNeighbourList();
		}
	}
}

void LinkedCellsOpenCL::countPairs() {

	numberOfPairs=0;
	numberOfSelfs=0;




	vector<unsigned long>::iterator cellIndexIter;
	vector<unsigned long>::iterator neighbourOffsetsIter;
	std::list<Molecule*>::iterator molIter1;
		// loop over all inner cells and calculate forces to forward neighbours
		for( cellIndexIter = _innerCellIndices.begin(); cellIndexIter != _innerCellIndices.end(); cellIndexIter++ ) {
			Cell& currentCell = _cells[*cellIndexIter];
			numberOfSelfs++;
			numberOfPairs++;
			// loop over all neighbours
			for( neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++ )
			{
				Cell& neighbourCell = _cells[*cellIndexIter + *neighbourOffsetsIter];
				if(neighbourCell.getParticlePointers().size() != 0)numberOfPairs++;
			}
		}


		// loop over all boundary cells and calculate forces to forward and backward neighbours
		for( cellIndexIter = _boundaryCellIndices.begin(); cellIndexIter != _boundaryCellIndices.end(); cellIndexIter++ )
		{
			Cell& currentCell = _cells[*cellIndexIter];
			numberOfSelfs++;
			numberOfPairs++;
			// loop over all forward neighbours
			for( neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++ )
			{
				Cell& neighbourCell = _cells[*cellIndexIter + *neighbourOffsetsIter];
				numberOfPairs++;
			}

			// loop over all backward neighbours. calculate only forces
			// to neighbour cells in the halo region, all others already have been calculated
			for( neighbourOffsetsIter = _backwardNeighbourOffsets.begin(); neighbourOffsetsIter != _backwardNeighbourOffsets.end(); neighbourOffsetsIter++ )
			{
				Cell& neighbourCell = _cells[*cellIndexIter+*neighbourOffsetsIter];
				if(neighbourCell.isHaloCell())
				{
					numberOfPairs++;
				}
			}
		}

}
void LinkedCellsOpenCL::createPairs() {
	resultSize=0;
	seflResultSize=0;
	int pi=0;
	vector<unsigned long>::iterator cellIndexIter;
	vector<unsigned long>::iterator neighbourOffsetsIter;
	std::list<Molecule*>::iterator molIter1;
		// loop over all inner cells and calculate forces to forward neighbours
		for( cellIndexIter = _innerCellIndices.begin(); cellIndexIter != _innerCellIndices.end(); cellIndexIter++ ) {
			Cell& currentCell = _cells[*cellIndexIter];
			pairs[pi++]=currentCell.id();
			pairs[pi++]=currentCell.id();
			pairs[pi++]=currentCell.getParticlePointers().size();
			pairs[pi++]=currentCell.getParticlePointers().size();
			pairs[pi++]=resultSize;
			resultSize+=currentCell.getParticlePointers().size()*currentCell.getParticlePointers().size();

			seflResultSize+=currentCell.getParticlePointers().size();//*currentCell.getParticlePointers().size();

			// loop over all neighbours
			for( neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++ )
			{
				Cell& neighbourCell = _cells[*cellIndexIter + *neighbourOffsetsIter];
				if(neighbourCell.getParticlePointers().size() != 0) {
				pairs[pi++]=currentCell.id();
				pairs[pi++]=neighbourCell.id();
				pairs[pi++]=currentCell.getParticlePointers().size();
				pairs[pi++]=neighbourCell.getParticlePointers().size();
				pairs[pi++]=resultSize;
				resultSize+=currentCell.getParticlePointers().size()*neighbourCell.getParticlePointers().size();
				}
			}
		}


		// loop over all boundary cells and calculate forces to forward and backward neighbours
		for( cellIndexIter = _boundaryCellIndices.begin(); cellIndexIter != _boundaryCellIndices.end(); cellIndexIter++ )
		{
			Cell& currentCell = _cells[*cellIndexIter];
			pairs[pi++]=currentCell.id();
			pairs[pi++]=currentCell.id();
			pairs[pi++]=currentCell.getParticlePointers().size();
			pairs[pi++]=currentCell.getParticlePointers().size();
			pairs[pi++]=resultSize;
			resultSize+=currentCell.getParticlePointers().size()*currentCell.getParticlePointers().size();

			seflResultSize+=currentCell.getParticlePointers().size();//*currentCell.getParticlePointers().size();
			// loop over all forward neighbours
			for( neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++ )
			{
				Cell& neighbourCell = _cells[*cellIndexIter + *neighbourOffsetsIter];
				pairs[pi++]=currentCell.id();
				pairs[pi++]=neighbourCell.id();
				pairs[pi++]=currentCell.getParticlePointers().size();
				pairs[pi++]=neighbourCell.getParticlePointers().size();
				pairs[pi++]=resultSize;
				resultSize+=currentCell.getParticlePointers().size()*neighbourCell.getParticlePointers().size();
			}
			// loop over all backward neighbours. calculate only forces
			// to neighbour cells in the halo region, all others already have been calculated
			for( neighbourOffsetsIter = _backwardNeighbourOffsets.begin(); neighbourOffsetsIter != _backwardNeighbourOffsets.end(); neighbourOffsetsIter++ )
			{
				Cell& neighbourCell = _cells[*cellIndexIter+*neighbourOffsetsIter];
				if(neighbourCell.isHaloCell())
				{
				 pairs[pi++]=currentCell.id();
				 pairs[pi++]=neighbourCell.id();
				 pairs[pi++]=currentCell.getParticlePointers().size();
				 pairs[pi++]=neighbourCell.getParticlePointers().size();
				 pairs[pi++]=resultSize;
				 resultSize+=currentCell.getParticlePointers().size()*neighbourCell.getParticlePointers().size();
				}
			}
		}
}

void LinkedCellsOpenCL::traversePairsFinish() {
	double params[15];
		double delta_r;
		bool knowparams = false;
		std::list<Molecule*>::iterator molIter1;
		vector<unsigned long>::iterator cellIndexIter;
		for( cellIndexIter = _innerCellIndices.begin(); cellIndexIter != _boundaryCellIndices.end(); cellIndexIter++ )
		{
			if( cellIndexIter == _innerCellIndices.end() )
				cellIndexIter = _boundaryCellIndices.begin();
			Cell& currentCell = _cells[*cellIndexIter];
			for( molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++ )
			{
				Molecule& molecule1 = **molIter1;
				if( molecule1.numTersoff() == 0 ) continue;
				if( !knowparams )
				{
					delta_r = molecule1.tersoffParameters(params);
					knowparams = true;
				}
				_particlePairsHandler->processTersoffAtom(molecule1, params, delta_r);
			}
		}

		_particlePairsHandler->finish();
}

int LinkedCellsOpenCL::countParticles() {
	int result=0;
	std::list<Molecule*>::iterator molIter1;
	for(unsigned i = 0; i < _cells.size(); i++ ){
		result+=_cells[i].getParticlePointers().size();
	}
	return result;
}

int LinkedCellsOpenCL::getMaxCellSizeAndPositions() {
	int maxSize=0;
	int index=0;
	int offset=0;
	std::list<Molecule*>::iterator molIter1;
	for(unsigned i = 0; i < _cells.size(); i++ ){
			Cell& currentCell = _cells[i];
			if (_cellsValid == false) {
					cout << "Cell structure in LinkedCells (traversePairs) invalid, call update first" << endl;
					exit(1);
				}

			currentCell.setId(offset);
			offset+=currentCell.getParticlePointers().size();
			//if(currentCell.getParticlePointers().size() != 0)cout << "cell size: " << currentCell.getParticlePointers().size()<< endl;
			maxSize=max(maxSize,(int)currentCell.getParticlePointers().size());
			if( currentCell.getParticlePointers().size() == 0 ) {
				continue;
			}
			for( molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++ ) {
				Molecule& molecule1 = **molIter1;
				if(index >= numberOfParticles) {
					cout << "found more particles then declared  "<< index<< endl;
					//exit(1);
				} else {
					m_positions[index*4]=molecule1.r(0);
					m_positions[index*4+1]=molecule1.r(1);
					m_positions[index*4+2]=molecule1.r(2);
					m_positions[index*4+3]=index;
				}
				index++;

			}
	}
	return maxSize;
}

#ifndef CUDA_HELPERS_H_
#define CUDA_HELPERS_H_

#include <cuda.h>
#include <vector_types.h>

#include <malloc.h>
#include <vector>

#include <assert.h>
#include <string>
#include <stddef.h>

// TODO: use a namespace instead of the slightly stupid CUDA prefix

// TODO: move Unpacked/PackerVector into their own file (and move this project-specific include with them)
#include "config.h"

namespace CUDADetail {
	template<typename DataType> struct TypeInfo {
		const static size_t size = sizeof( DataType );
		const static size_t alignof = __alignof( DataType );
	};

	template<> struct TypeInfo<double3> {
		const static size_t size = sizeof( double3 );
		const static size_t alignof = 8;
	};

	template<typename DataType> struct TypeInfo<DataType *> {
		const static size_t size = sizeof( CUdeviceptr );
		const static size_t alignof = __alignof( CUdeviceptr );
	};
}

class CUDA {
public:
	struct Exception : public std::exception {
		const CUresult errorCode;
		const std::string errorSource;

		Exception( CUresult errorCode, const std::string &errorSource = "" )
		: errorCode( errorCode ), errorSource( errorSource ) {}

		~Exception() throw() {}

		/** Returns a C-style character string describing the general cause
		 *  of the current error.  */
		virtual const char* what() const throw() {
			return errorSource.c_str();
		}
	};

	#define CUDA_THROW_ON_ERROR_EX( expr, msgformat, ... ) do {\
			CUresult result = (expr);\
			if( result != CUDA_SUCCESS ) {\
				printf( "CUDA error %i in '" #expr "'!\n" msgformat, result, ##__VA_ARGS__ );\
				throw CUDA::Exception( result, #expr ); \
			}\
		} while(false)

	#define CUDA_THROW_ON_ERROR( expr ) CUDA_THROW_ON_ERROR_EX( expr, "" )

	template<typename DataType>
	class DeviceBuffer {
	protected:
		CUdeviceptr _deviceBuffer;

		int _byteSize;

	public:
		DeviceBuffer() : _deviceBuffer( 0 ), _byteSize( 0 ) {
		}

		~DeviceBuffer() {
			if( _byteSize )
				CUDA_THROW_ON_ERROR( cuMemFree( _deviceBuffer ) );
		}

		void resize(int count) {
			if( _byteSize )
				CUDA_THROW_ON_ERROR( cuMemFree( _deviceBuffer ) );

			_byteSize = count * sizeof( DataType );
			if( _byteSize > 0 ) {
				CUDA_THROW_ON_ERROR_EX( cuMemAlloc( &_deviceBuffer, _byteSize ), "\twhile allocating %i MB\n", _byteSize / (1<<20) );
			}
			else {
				_deviceBuffer = 0;
				// if _byteSize < 0
				_byteSize = 0;
			}
		}

		void zeroDevice() const {
			CUDA_THROW_ON_ERROR( cuMemsetD8( _deviceBuffer, 0, _byteSize ) );
		}

		// returns true if the buffer has been reallocated
		bool copyToDevice(const DataType *data, int count) {
			const int byteSize = count * sizeof( DataType );
			const bool needsResize = _byteSize < byteSize;
			if( needsResize )
				resize( count );

			CUDA_THROW_ON_ERROR( cuMemcpyHtoD( _deviceBuffer, data, byteSize ) );
			return needsResize;
		}

		bool copyToDevice(const std::vector<DataType> &collection) {
			return copyToDevice( &collection.front(), collection.size() );
		}

		void copyToHost(std::vector<DataType> &collection) const {
			if( !_byteSize ) {
				return;
			}

			int count = _byteSize / sizeof( DataType );
			collection.resize( count );

			CUDA_THROW_ON_ERROR( cuMemcpyDtoH( &collection.front(), _deviceBuffer, count * sizeof( DataType ) ) );
		}

		CUdeviceptr devicePtr() const {
			return _deviceBuffer;
		}
	};

	class Timer {
	private:
		CUevent _startEvent, _endEvent;

		Timer( const Timer & );
		Timer & operator =( const Timer & );
	public:
		Timer() {
			CUDA_THROW_ON_ERROR( cuEventCreate( &_startEvent, ::CU_EVENT_DEFAULT ) );
			CUDA_THROW_ON_ERROR( cuEventCreate( &_endEvent, ::CU_EVENT_DEFAULT ) );
		}

		~Timer() {
			CUDA_THROW_ON_ERROR( cuEventDestroy( _startEvent ) );
			CUDA_THROW_ON_ERROR( cuEventDestroy( _endEvent ) );
		}

		void begin() {
			CUDA_THROW_ON_ERROR( cuEventRecord( _startEvent, 0 ) );
		}

		void end() {
			CUDA_THROW_ON_ERROR( cuEventRecord( _endEvent, 0 ) );
		}

		float getElapsedTime() {
			CUDA_THROW_ON_ERROR( cuEventSynchronize( _endEvent ) );

			float elapsedTime;
			CUDA_THROW_ON_ERROR( cuEventElapsedTime( &elapsedTime, _startEvent, _endEvent ) );

			return elapsedTime;
		}
	};

	class NullTimer {
	public:
		void begin() {
		}

		void end() {
		}

		float getElapsedTime() {
			return 0.0f;
		}
	};

	typedef Timer EventTimer;

	class Function;
	class Module;

	template<typename DataType>
	class Global {
	protected:
		CUdeviceptr _dataPointer;

		friend class Module;

		Global( CUdeviceptr dataPointer ) : _dataPointer( dataPointer ) {}

	public:
		void set( const DataType &data ) const {
			CUDA_THROW_ON_ERROR( cuMemcpyHtoD( _dataPointer, &data, sizeof( DataType ) ) );
		}
	};

	template<typename DataType>
	class Global<DataType *> {
	protected:
		CUdeviceptr _dataPointer;

		friend class Module;

		Global( CUdeviceptr dataPointer ) : _dataPointer( dataPointer ) {}

	public:
		void set( const CUdeviceptr data ) const {
			CUDA_THROW_ON_ERROR( cuMemcpyHtoD( _dataPointer, &data, sizeof( CUdeviceptr ) ) );
		}

		void set( const DeviceBuffer<DataType> &buffer ) const {
			set( buffer.devicePtr() );
		}
	};

	template<typename DataType, uint size>
	class GlobalArray {
	protected:
		CUdeviceptr _dataPointer;

		friend class Module;

		GlobalArray( CUdeviceptr dataPointer ) : _dataPointer( dataPointer ) {}

	public:
		void set( int index, const DataType &data ) const {
			CUDA_THROW_ON_ERROR( cuMemcpyHtoD( _dataPointer + index * sizeof( DataType ), &data, sizeof( DataType ) ) );
		}
	};

	template<typename DataType, uint size>
	class GlobalArray< DataType *, size > {
	protected:
		CUdeviceptr _dataPointer;

		friend class Module;

		GlobalArray( CUdeviceptr dataPointer ) : _dataPointer( dataPointer ) {}

	public:
		static const uint ArraySize = size;

		void set( int index, const CUdeviceptr data ) const {
			CUDA_THROW_ON_ERROR( cuMemcpyHtoD( _dataPointer + index * sizeof( CUdeviceptr ), &data, sizeof( CUdeviceptr ) ) );
		}

		void set( int index, const DeviceBuffer<DataType> &buffer ) const {
			set( index, buffer.devicePtr() );
		}
	};

	template<typename DataType, uint size>
	class _UnpackedGlobalVectorBase {
	protected:
		GlobalArray< DataType *, size > _devicePointers;
		DeviceBuffer< DataType > _deviceBuffers[size];

		friend class Module;

		std::vector<DataType> _hostBuffers[size];

	public:
		_UnpackedGlobalVectorBase(const GlobalArray< DataType *, size > &devicePointers)
			: _devicePointers( devicePointers )
		{
		}

		void clear() {
			for( int i = 0; i < size ; i++ ) {
				_hostBuffers[ i ].clear();
			}
		}

		void push( const DataType entry[size] ) {
			for( int i = 0 ; i < size ; i++ ) {
				_hostBuffers[ i ].push_back( entry[ i ] );
			}
		}

		void push( int startIndex, const DataType *entry, int count = 1 ) {
			assert( size >= startIndex + count );
			for( int i = 0 ; i < count ; i++ ) {
				const int j = startIndex + i;
				_hostBuffers[ j ].push_back( entry[ j ] );
			}
		}

		void updateDevice() {
			for( int i = 0 ; i < size ; i++ ) {
				_deviceBuffers[ i ].copyToDevice( _hostBuffers[ i ] );
				_devicePointers.set( i, _deviceBuffers[ i ] );
			}
		}

		void zeroDevice() const {
			for( int i = 0 ; i < size ; i++ ) {
				_deviceBuffers[ i ].zeroDevice();
			}
		}

		void resize( int numElements ) {
			for( int i = 0 ; i < size ; i++ ) {
				_hostBuffers[ i ].resize( numElements );

				_deviceBuffers[ i ].resize( numElements );
				_devicePointers.set( i, _deviceBuffers[ i ] );
			}
		}

		void copyToHost() {
			for( int i = 0 ; i < size ; i++ ) {
				_deviceBuffers[ i ].copyToHost( _hostBuffers[ i ] );
			}
		}

		void readBack( int index, DataType *entry ) {
			for( int i = 0 ; i < size ; i++ ) {
				entry[ i ] = _hostBuffers[ i ][ index ];
			}
		}
	};

	template<typename DataType, uint size>
	class UnpackedGlobalVector : public _UnpackedGlobalVectorBase<DataType, size> {
	public:
		UnpackedGlobalVector(const GlobalArray< DataType *, size > &devicePointers)
			: _UnpackedGlobalVectorBase<DataType, size>( devicePointers )
		{
		}
	};

	template<uint size>
	class UnpackedGlobalVector<floatType, size> : public _UnpackedGlobalVectorBase<floatType, size> {
		typedef _UnpackedGlobalVectorBase<floatType, size> base;
	public:
		UnpackedGlobalVector(const GlobalArray<floatType*, size > &devicePointers)
			: _UnpackedGlobalVectorBase<floatType, size>( devicePointers )
		{
		}

		void push( int startIndex, const floatType3 &entry ) {
			assert( size >= startIndex + 3 );
			this->_hostBuffers[ startIndex ].push_back( entry.x );
			this->_hostBuffers[ startIndex + 1 ].push_back( entry.y );
			this->_hostBuffers[ startIndex + 2 ].push_back( entry.z );
		}

		void push( const floatType3 &entry ) {
			assert( size == 3 );
			push( 0, entry );
		}

		void push( const QuaternionStorage &quaternion ) {
			assert( size == 4 );
			base::push( 0, (const floatType*) &quaternion, 4 );
		}

		void push( const Matrix3x3Storage &mat ) {
			assert( size == 9 );
			push( 0, mat.rows[0] );
			push( 3, mat.rows[0] );
			push( 6, mat.rows[0] );
		}

		void copyToHost( std::vector<floatType3> &collection ) {
			assert( size == 3 );

			base::copyToHost();

			collection.clear();

			const int numElements = this->_hostBuffers[0].size();
			collection.resize( numElements );
			for( int index = 0 ; index < numElements ; index++ ) {
				collection[index] = make_floatType3(
						this->_hostBuffers[0][index],
						this->_hostBuffers[1][index],
						this->_hostBuffers[2][index]
					);
			}
		}
	};

	template<typename DataType>
	class PackedGlobalVector {
	protected:
		Global<DataType *> _devicePointer;
		DeviceBuffer<DataType> _deviceBuffer;

		friend class Module;

		std::vector<DataType> _hostBuffer;

	public:
		PackedGlobalVector(const Global< DataType * > &devicePointer)
			: _devicePointer( devicePointer )
		{
		}

		void clear() {
			_hostBuffer.clear();
		}

		void push( const DataType &entry) {
			_hostBuffer.push_back( entry );
		}

		void updateDevice() {
			_deviceBuffer.copyToDevice( _hostBuffer );
			_devicePointer.set( _deviceBuffer );
		}

		void zeroDevice() const {
			_deviceBuffer.zeroDevice();
		}

		void resize( int numElements ) {
			_hostBuffer.resize( numElements );

			_deviceBuffer.resize( numElements );
			_devicePointer.set( _deviceBuffer );
		}

		std::vector<DataType> & copyToHost() {
			_deviceBuffer.copyToHost( _hostBuffer );
			return getHostBuffer();
		}

		void copyToHost( std::vector< DataType > &collection ) {
			_deviceBuffer.copyToHost( _hostBuffer );
			collection.assign( _hostBuffer.begin(), _hostBuffer.end() );
		}

		void readBack( int index, DataType &entry ) {
			entry = _hostBuffer[index];
		}

		std::vector<DataType> & getHostBuffer() {
			return _hostBuffer;
		}
	};

	template<typename DataType>
	class GlobalDeviceBuffer {
	protected:
		Global<DataType *> _devicePointer;
		DeviceBuffer<DataType> _deviceBuffer;

		friend class Module;
	public:
		GlobalDeviceBuffer(const Global< DataType * > &devicePointer)
			: _devicePointer( devicePointer )
		{
		}

		void zeroDevice() const {
			_deviceBuffer.zeroDevice();
		}

		void resize( int numElements ) {
			_deviceBuffer.resize( numElements );
			_devicePointer.set( _deviceBuffer );
		}
	};

	class FunctionCall {
	protected:
		CUfunction _function;
		int _offset;
		bool _blockShapeSet;

		friend class Function;

		FunctionCall( const CUfunction &function ) : _function( function ), _offset( 0 ), _blockShapeSet( false ) {}

		void defaultBlockShape() const {
			if( !_blockShapeSet ) {
				CUDA_THROW_ON_ERROR( cuFuncSetBlockShape( _function, 1, 1, 1 ) );
			}
		}
	public:
		template<typename T>
		FunctionCall & parameter( const T &param ) {
			using namespace CUDADetail;

			// align with parameter size
			_offset = (_offset + (TypeInfo<T>::alignof - 1)) & ~(TypeInfo<T>::alignof - 1);
			CUDA_THROW_ON_ERROR( cuParamSetv( _function, _offset, (void *) &param, TypeInfo<T>::size ) );
			_offset += TypeInfo<T>::size;

			return *this;
		}

		template<typename T>
		FunctionCall & parameter( const DeviceBuffer<T> &param ) {
			using namespace CUDADetail;

			// align with parameter size
			_offset = (_offset + (TypeInfo<CUdeviceptr>::alignof - 1)) & ~(TypeInfo<CUdeviceptr>::alignof - 1);

			CUdeviceptr devicePtr = param.devicePtr();
			CUDA_THROW_ON_ERROR( cuParamSetv( _function, _offset, (void *) &devicePtr, TypeInfo<CUdeviceptr>::size ) );

			_offset += TypeInfo<CUdeviceptr>::size;

			return *this;
		}

		FunctionCall & setBlockShape( int x, int y, int z ) {
			CUDA_THROW_ON_ERROR( cuFuncSetBlockShape( _function, x, y, z ) );

			_blockShapeSet = true;
			return *this;
		}

		FunctionCall & setBlockShape( const dim3 &shape) {
			CUDA_THROW_ON_ERROR( cuFuncSetBlockShape( _function, shape.x, shape.y, shape.z ) );

			_blockShapeSet = true;
			return *this;
		}

		void execute(int gridWidth, int gridHeight) const {
			CUDA_THROW_ON_ERROR( cuParamSetSize( _function, _offset ) );

			defaultBlockShape();

			CUDA_THROW_ON_ERROR_EX( cuLaunchGrid( _function, gridWidth, gridHeight ), "\twhile launching a grid of size %i x %i\n", gridWidth, gridHeight );
		}

		void execute(int numJobs) const {
			execute( numJobs, 1 );
		}

		// execute with at least numJobs thread blocks
		// at most (1<<(ceil(log2(numJobs)-15)) - 1 thread blocks too many
		void executeAtLeast( int numJobs ) const {
			for( int log2_y = 0 ; log2_y < 15 ; log2_y++ ) {
				int y = 1 << log2_y;
				int x = ((numJobs + y - 1) >> log2_y);
				if( x < 65536 ) {
					execute( x, y );
					return;
				}
			}

			CUDA_THROW_ON_ERROR_EX( CUDA_ERROR_INVALID_VALUE, "\tnumJobs = %i > 65535^2\n", numJobs );
		}

		void execute() const {
			CUDA_THROW_ON_ERROR( cuParamSetSize( _function, _offset ) );

			defaultBlockShape();

			CUDA_THROW_ON_ERROR( cuLaunch( _function ) );
		}
	};

	class Function {
	protected:
		CUfunction _function;

		friend class Module;

		Function( CUfunction function ) : _function( function ) {}

	public:
		FunctionCall call() const {
			return FunctionCall( _function );
		}
	};

	class Module {
	protected:
		CUmodule _module;

		friend class CUDA;

		Module(CUmodule module) : _module( module ) {}

	public:

		template<typename DataType>
		Global<DataType> getGlobal(const char *name) const {
			CUdeviceptr dptr;
			size_t dataSize;

			CUDA_THROW_ON_ERROR_EX( cuModuleGetGlobal( &dptr, &dataSize, _module, name ), "\twhen trying to access global '%s'\n", name );

			assert( dataSize == CUDADetail::TypeInfo<DataType>::size );

			return Global<DataType>(dptr);
		}

		template<typename DataType, uint size>
		GlobalArray<DataType, size> getGlobalArray(const char *name) const {
			CUdeviceptr dptr;
			size_t dataSize;

			CUDA_THROW_ON_ERROR_EX( cuModuleGetGlobal( &dptr, &dataSize, _module, name ), "\twhen trying to access array global '%s'\n", name );

			assert( dataSize == CUDADetail::TypeInfo<DataType>::size * size );

			return GlobalArray<DataType, size>(dptr);
		}

		Function getFunction(const char *name) const {
			CUfunction function;

			CUDA_THROW_ON_ERROR( cuModuleGetFunction( &function, _module, name ) );

			return Function( function );
		}
	};

protected:
	CUdevice device;
	CUcontext context;

	CUDA(int deviceIndex, bool preferL1Cache, size_t printfBufferSize) {
		CUDA_THROW_ON_ERROR( cuInit( 0 ) );

		CUDA_THROW_ON_ERROR( cuDeviceGet( &device, deviceIndex ) );

		printMemoryInfo();

		CUDA_THROW_ON_ERROR( cuCtxCreate( &context, 0, device ) );

		CUDA_THROW_ON_ERROR( cuCtxSetCacheConfig( preferL1Cache ? CU_FUNC_CACHE_PREFER_L1 : CU_FUNC_CACHE_PREFER_SHARED ) );

		CUDA_THROW_ON_ERROR( cuCtxSetLimit( CU_LIMIT_PRINTF_FIFO_SIZE, printfBufferSize ) );

		size_t actualLimit;
		CUDA_THROW_ON_ERROR( cuCtxGetLimit( &actualLimit, CU_LIMIT_PRINTF_FIFO_SIZE ) );
		printf( "actual CUDA printf buffer size: %i KB (wanted: %i KB)\n", actualLimit >> 10, printfBufferSize >> 10 );
	}

	~CUDA() {
		CUDA_THROW_ON_ERROR( cuCtxDestroy( context ) );

		printMemoryInfo();
	}

	void printMemoryInfo() {
		CUcontext infoContext;
		CUDA_THROW_ON_ERROR( cuCtxCreate( &infoContext, 0, device ) );

		size_t free, total;
		CUDA_THROW_ON_ERROR( cuMemGetInfo( &free, &total ) );

		CUDA_THROW_ON_ERROR( cuCtxDestroy( infoContext ) );

		printf( "CUDA memory info: %i MB free (%i MB total); %f%% free\n", free / (1<<20), total / (1<<20), (float) free / total );
	}

	static CUDA *singleton;

public:

	static CUDA &get() {
		assert( singleton );
		return *singleton;
	}

	static void create(int deviceIndex, bool preferL1Cache, size_t printfBufferSize) {
		assert( !singleton );
		singleton = new CUDA(deviceIndex, preferL1Cache, printfBufferSize);
	}

	static void destruct() {
		delete singleton;
		singleton = NULL;
	}

	Module loadModule(const char *name) {
		CUmodule module;
		CUDA_THROW_ON_ERROR( cuModuleLoad( &module, name ) );

		return Module( module );
	}

	Module loadModuleData(const char *imageStart, const char *imageEnd, int maxRegisterCount) {
		std::string image(imageStart, imageEnd);

		CUmodule module;
		CUjit_option options = CU_JIT_MAX_REGISTERS;
		CUDA_THROW_ON_ERROR( cuModuleLoadDataEx( &module, image.c_str(), 1, &options, (void**) &maxRegisterCount ) );

		return Module( module );
	}
};

inline CUDA &cuda() {
	return CUDA::get();
}

#endif

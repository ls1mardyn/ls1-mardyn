#ifndef MOLECULEPAIRHANDLER_H_
#define MOLECULEPAIRHANDLER_H_

#include "cudaComponent.h"

class MoleculePairHandler : public CUDAStaticDataComponent {
public:
	MoleculePairHandler( const CUDAComponent &component ) :
		CUDAStaticDataComponent( component ),
		_cutOffRadiusSquared( _module.getGlobal<floatType>("cutOffRadiusSquared") ),
		_epsRFInvrc3( _module.getGlobal<floatType>("epsRFInvrc3") ) {
	}

	virtual void upload() {
		const floatType cutOffRadius = _linkedCells.getCutoff();
		_cutOffRadiusSquared.set( cutOffRadius * cutOffRadius );

		const floatType epsRF = _domain.getepsilonRF();
		const floatType epsRFInvrc3 = 2. * (epsRF - 1.) / ((cutOffRadius * cutOffRadius * cutOffRadius) * (2. * epsRF + 1.));
		_epsRFInvrc3.set( epsRFInvrc3 );
	}

protected:
	CUDA::Global<floatType> _cutOffRadiusSquared;
	CUDA::Global<floatType> _epsRFInvrc3;
};

#endif /* GLOBALSTATS_H_ */

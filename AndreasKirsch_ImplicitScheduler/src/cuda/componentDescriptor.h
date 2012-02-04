/*
 * componentDescriptor.h
 *
 *  Created on: Jun 21, 2011
 *      Author: andreas
 */

#ifndef COMPONENTDESCRIPTOR_H_
#define COMPONENTDESCRIPTOR_H_

#include <vector>

#include "utils/Logger.h"

#include "cutil_math.h"

#include "sharedDecls.h"

class ComponentDescriptorStorage : public CUDAStaticDataComponent {
public:
	ComponentDescriptorStorage( const CUDAComponent &component ) :
		CUDAStaticDataComponent( component ),

		_componentDescriptors( _module.getGlobal<ComponentDescriptors>("componentDescriptors") ),
		_componentMixXis( _module.getGlobal<ComponentMixCoefficients>("componentMixXis") ),
		_componentMixEtas( _module.getGlobal<ComponentMixCoefficients>("componentMixEtas") ),

		_debugComponentDescriptors( _module.getFunction("debugComponentDescriptors") ) {
	}

	virtual void upload() {
		ComponentMixCoefficients xis;
		ComponentMixCoefficients etas;
		ComponentDescriptors componentDescriptors;

		const std::vector<Component> &components = _domain.getComponents();

		const int numComponents = components.size();
		if( numComponents > MAX_NUM_COMPONENTS ) {
			// error
			Log::global_log->fatal() << "Domain has " << numComponents <<
					" components, but MAX_NUM_COMPONENTS = " << MAX_NUM_COMPONENTS << std::endl;
			exit(-1);
		}

		for( int i = 0 ; i < numComponents ; i++ ) {
			const Component &component = components[i];
			ComponentDescriptor &componentDescriptor = componentDescriptors[i];

			componentDescriptor.numLJCenters = component.numLJcenters();
			componentDescriptor.numCharges = component.numCharges();
			componentDescriptor.numDipoles = component.numDipoles();
			componentDescriptor.numQuadrupoles = component.numQuadrupoles();

			if( componentDescriptor.numLJCenters > MAX_NUM_LJCENTERS ) {
				// error
				Log::global_log->fatal() << "A component has " << componentDescriptor.numLJCenters <<
						" lj centers, but MAX_NUM_LJCENTERS = " << MAX_NUM_LJCENTERS << std::endl;
				exit(-1);
			}
			if( componentDescriptor.numCharges > MAX_NUM_CHARGES ) {
				// error
				Log::global_log->fatal() << "A component has " << componentDescriptor.numCharges <<
						" charges, but MAX_NUM_CHARGES = " << MAX_NUM_CHARGES << std::endl;
				exit(-1);
			}
			if( componentDescriptor.numDipoles > MAX_NUM_DIPOLES ) {
				// error
				Log::global_log->fatal() << "A component has " << componentDescriptor.numDipoles <<
						" dipoles, but MAX_NUM_DIPOLES = " << MAX_NUM_DIPOLES << std::endl;
				exit(-1);
			}
			if( componentDescriptor.numQuadrupoles > MAX_NUM_QUADRUPOLES ) {
				// error
				Log::global_log->fatal() << "A component has " << componentDescriptor.numQuadrupoles <<
						" quadrupoles, but MAX_NUM_QUADRUPOLES = " << MAX_NUM_QUADRUPOLES << std::endl;
				exit(-1);
			}

			// TODO: use inheritance for relativePosition?
#if MAX_NUM_LJCENTERS > 0
			for( int ljCenterIndex = 0 ; ljCenterIndex < componentDescriptor.numLJCenters ; ljCenterIndex++ ) {
				ComponentDescriptor::LJCenter &ljCenter = componentDescriptor.ljCenters[ljCenterIndex];

				const LJcenter &cLjCenter = component.ljcenter(ljCenterIndex);
				ljCenter.ljParameters.epsilon = cLjCenter.eps();
				ljCenter.ljParameters.sigma = cLjCenter.sigma();
				ljCenter.relativePosition = make_floatType3( cLjCenter.rx(), cLjCenter.ry(), cLjCenter.rz() );
			}
#endif

#if MAX_NUM_CHARGES > 0
			for( int chargeIndex = 0 ; chargeIndex < componentDescriptor.numCharges ; chargeIndex++ ) {
				ComponentDescriptor::Charge &charge = componentDescriptor.charges[chargeIndex];

				const Charge &cCharge = component.charge(chargeIndex);
				charge.relativePosition = make_floatType3( cCharge.rx(), cCharge.ry(), cCharge.rz() );

				charge.q = cCharge.q();
			}
#endif

#if MAX_NUM_DIPOLES > 0
			for( int dipoleIndex = 0 ; dipoleIndex < componentDescriptor.numDipoles ; dipoleIndex++ ) {
				ComponentDescriptor::Dipole &dipole = componentDescriptor.dipoles[dipoleIndex];

				const Dipole &cDipole = component.dipole(dipoleIndex);
				dipole.relativePosition = make_floatType3( cDipole.rx(), cDipole.ry(), cDipole.rz() );
				dipole.relativeE = make_floatType3( cDipole.ex(), cDipole.ey(), cDipole.ez() );

				dipole.absMu = cDipole.absMy();
			}
#endif

#if MAX_NUM_QUADRUPOLES > 0
			for( int quadrupoleIndex = 0 ; quadrupoleIndex < componentDescriptor.numQuadrupoles ; quadrupoleIndex++ ) {
				ComponentDescriptor::Quadrupole &quadrupole = componentDescriptor.quadrupoles[quadrupoleIndex];

				const Quadrupole &cQuadrupole = component.quadrupole(quadrupoleIndex);
				quadrupole.relativePosition = make_floatType3( cQuadrupole.rx(), cQuadrupole.ry(), cQuadrupole.rz() );
				quadrupole.relativeE = make_floatType3( cQuadrupole.ex(), cQuadrupole.ey(), cQuadrupole.ez() );

				quadrupole.absQ = cQuadrupole.absQ();
			}
#endif
		}

		// set mix coeff tables
		std::vector<double> &dmixcoeff = _domain.getmixcoeff();
		int index = 0;
		for( int i = 0 ; i < components.size() ; i++ ) {
			xis[i][i] = etas[i][i] = 1.0;
			for( int j = i + 1 ; j < components.size() ; j++ ) {
				floatType xi = dmixcoeff[index++];
				floatType eta = dmixcoeff[index++];

				xis[i][j] = xis[j][i] = xi;
				etas[i][j] = etas[j][i] = eta;
			}
		}

		_componentDescriptors.set( componentDescriptors );

		_componentMixXis.set( xis );

		_componentMixEtas.set( etas );

#ifdef DEBUG_COMPONENT_DESCRIPTORS
		_debugComponentDescriptors.call().parameter( components.size() ).execute();
#endif
	}

protected:
	CUDA::Global<ComponentDescriptors> _componentDescriptors;
	CUDA::Global<ComponentMixCoefficients> _componentMixXis, _componentMixEtas;

	CUDA::Function _debugComponentDescriptors;
};


#endif /* COMPONENTDESCRIPTOR_H_ */

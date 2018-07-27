#include "Wall.h"

#include <vector>
#include <array>
#include <cmath>
#include <cstdint>

using namespace std;
using Log::global_log;

Wall::Wall()
	: _rhoW(0.),
	_yc(0.),
	_yOff(0.),
	Delta(0.),
	_eps_wi(NULL),
	_sig3_wi(NULL),
	_sig2_wi(NULL),
	_sig_wi(NULL),
	_uShift_9_3(NULL),
	_uPot_9_3(NULL),
	_uShift_10_4(NULL),
	_uPot_10_4(NULL),
	_nc(0),
	_dWidth(0.),
	_dWidthHalf(0.)
{
}

Wall::~Wall()
{
	// free memory
	delete [] _eps_wi;
	delete [] _sig3_wi;
	delete [] _sig2_wi;
	delete [] _sig_wi;
	delete [] _uShift_9_3;
	delete [] _uPot_9_3;
	delete [] _uShift_10_4;
	delete [] _uPot_10_4;
}

void Wall::readXML(XMLfileUnits& xmlconfig)
{
	double density, sigma, epsilon, yoff, ycut;
	xmlconfig.getNodeValue("density", density);
	xmlconfig.getNodeValue("sigma", sigma);
	xmlconfig.getNodeValue("epsilon", epsilon);
	xmlconfig.getNodeValue("yoff", yoff);
	xmlconfig.getNodeValue("ycut", ycut);
	xmlconfig.getNodeValue("width", _dWidth);
	_dWidthHalf = _dWidth * 0.5;
	global_log->info() << "Using feature 'Wallfun' with parameters: density=" << density << ", "
			"sigma=" << sigma << ", epsilon=" << epsilon << ", yoff=" << yoff << ", ycut=" << ycut << ", "
			"width=" << _dWidth << endl;

	XMLfile::Query query = xmlconfig.query("component");
	unsigned int numComponentsConsidered = query.card();
	global_log->info() << "Wallfun: Setting parameters 'xi', 'eta' for " << numComponentsConsidered << " components." << endl;

	std::vector<Component>* components = global_simulation->getEnsemble()->getComponents();
	unsigned int numComponents = components->size();
	std::vector<double> xi_sf(numComponents);
	std::vector<double> eta_sf(numComponents);

	_bConsiderComponent.resize(numComponents);
	for(auto&& bi : _bConsiderComponent)
		bi = false;

	string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator componentIter;
	for( componentIter = query.begin(); componentIter; componentIter++ )
	{
		xmlconfig.changecurrentnode( componentIter );
		unsigned int cid;
		xmlconfig.getNodeValue("@id", cid);
		xmlconfig.getNodeValue("xi",   xi_sf.at(cid-1) );
		xmlconfig.getNodeValue("eta", eta_sf.at(cid-1) );
		_bConsiderComponent.at(cid-1) = true;
	}
	xmlconfig.changecurrentnode(oldpath);

	this->initializeLJ93(components, density, sigma, epsilon, xi_sf, eta_sf, yoff, ycut);
	global_log->info() << "Wallfun initialized." << endl;
}

void Wall::initializeLJ93(const std::vector<Component>* components, 
		 double in_rhoWall, double in_sigWall, double in_epsWall, std::vector<double> in_xi, std::vector<double> in_eta,
		 double in_yOffWall, double in_yWallCut) {
	global_log->info() << "Initializing the wall function LJ-9-3.\n";
	this->_rhoW = in_rhoWall;
	this->_yc = in_yWallCut;
	this->_yOff = in_yOffWall;

	/*!*** So far: only 1CLJ components allowed ****/
	_nc = components->size();
	_eps_wi = new double[_nc];
	_sig3_wi = new double[_nc];
	_uShift_9_3 = new double[_nc];
	_uPot_9_3 = new double[_nc];

	for (unsigned i = 0; i < _nc; i++) {
		_eps_wi[i] = in_xi.at(i) * sqrt(in_epsWall * (components->at(i)).ljcenter(0).eps());

		double sig_wi = 0.5 * in_eta[i] * (in_sigWall + (components->at(i)).ljcenter(0).sigma());
		_sig3_wi[i] = sig_wi * sig_wi * sig_wi;
		double sig9_sf = _sig3_wi[i] * _sig3_wi[i] * _sig3_wi[i];

		double y = _yc - _yOff;
		double y3 = y * y * y;
		double y9 = y3 * y3 * y3;

		_uShift_9_3[i] = 4.0 / 3.0 * M_PI * _rhoW * _eps_wi[i] *_sig3_wi[i] * (sig9_sf / 15.0 / y9 - _sig3_wi[i] / 2.0 / y3);
		_uPot_9_3[i] = 0.0;
	}
}

void Wall::initializeLJ104(const std::vector<Component>* components, 
		 double in_rhoWall, double in_sigWall, double in_epsWall, std::vector<double> in_xi, std::vector<double> in_eta,
		 double in_yOffWall, double in_yWallCut) {
	global_log->info() << "Initializing the wall function LJ-10-4.\n";
	this->_rhoW = in_rhoWall;
	this->_yc = in_yWallCut;
	this->_yOff = in_yOffWall;

	/*!*** So far: only 1CLJ components allowed ****/
	_nc = components->size();
	_eps_wi = new double[_nc];
	_sig3_wi = new double[_nc];
	_sig2_wi = new double[_nc];
	_sig_wi = new double[_nc];
	_uShift_10_4 = new double[_nc];
	_uPot_10_4 = new double[_nc];

	for (unsigned i = 0; i < _nc; i++) {
		_eps_wi[i] = in_xi.at(i) * sqrt(in_epsWall * (components->at(i)).ljcenter(0).eps());

		double sig_wi = 0.5 * in_eta[i] * (in_sigWall + (components->at(i)).ljcenter(0).sigma());
		_sig_wi[i] = sig_wi;
		_sig2_wi[i] = sig_wi * sig_wi;
		_sig3_wi[i] = sig_wi * sig_wi * sig_wi;
		double sig10_sf = _sig3_wi[i] * _sig3_wi[i] * _sig3_wi[i] * _sig_wi[i];
		double sig4_sf = _sig2_wi[i] * _sig2_wi[i];

		double y = _yc - _yOff;
		double y2 = y * y;
		double y4 = y2 * y2;
		double y10 = y4 * y4 * y2;

		double bracket = y + 0.61 * Delta;
		double bracket3 = bracket * bracket * bracket;

		_uShift_10_4[i] = 2 * M_PI * _eps_wi[i] * _rhoW * _sig2_wi[i] * Delta * (2 / 5 * sig10_sf / y10 - sig4_sf / y4 - sig4_sf / (3 * Delta * bracket3));
		_uPot_10_4[i] = 0.0;
	}
}

void Wall::calcTSLJ_9_3( ParticleContainer* partContainer, Domain* domain)
{
	double regionLowCorner[3], regionHighCorner[3];

	/*! LJ-9-3 potential applied in y-direction */
	for(uint8_t d=0; d<3; ++d){
		regionLowCorner[d]  = partContainer->getBoundingBoxMin(d);
		regionHighCorner[d] = partContainer->getBoundingBoxMax(d);
	}
	double dWallRangeLeft  = _yOff-_dWidthHalf - _yc;
	double dWallRangeRight = _yOff+_dWidthHalf + _yc;
	regionLowCorner[1]  = (dWallRangeLeft  > regionLowCorner[1]  ) ? dWallRangeLeft  : regionLowCorner[1];
	regionHighCorner[1] = (dWallRangeRight < regionHighCorner[1] ) ? dWallRangeRight : regionHighCorner[1];

	//perform a check if the region is contained by the particleContainer???
	if ( false == partContainer->isRegionInBoundingBox(regionLowCorner, regionHighCorner) )
		return;

	#if defined (_OPENMP)
	#pragma omp parallel shared(regionLowCorner, regionHighCorner)
	#endif
	{
		RegionParticleIterator begin = partContainer->iterateRegionBegin(regionLowCorner, regionHighCorner);
		RegionParticleIterator end = partContainer->iterateRegionEnd();

		double f[3];
		for(unsigned d=0; d<3; d++)
			f[d] = 0.;

		for(RegionParticleIterator mi = begin; mi != end; ++mi)
		{
			unsigned int cid = mi->componentid();
			if(false == _bConsiderComponent.at(cid) )
				continue;  // only add Wall force to molecules of component that should be considered

			for(unsigned int si=0; si<mi->numLJcenters(); ++si)
			{
				double y, y3, y9, ry, ryRel;
				const std::array<double,3> arrSite = mi->ljcenter_d_abs(si);
				const double* posSite = arrSite.data();
				ry = posSite[1];
				ryRel = (ry > _yOff) ? (ry - (_yOff+_dWidthHalf) ) : (ry - (_yOff-_dWidthHalf) );
				y = abs(ryRel);
				if(y < _yc){
					y3 = y * y * y;
					y9 = y3 * y3 * y3;

					double sig9_wi;
					sig9_wi = _sig3_wi[cid] * _sig3_wi[cid] * _sig3_wi[cid];
					f[1] = 4.0 * M_PI * _rhoW * _eps_wi[cid] * _sig3_wi[cid] * (sig9_wi / 5.0 / y9 - _sig3_wi[cid] / 2.0 / y3) / ryRel;
					_uPot_9_3[cid] += 4.0 * M_PI * _rhoW * _eps_wi[cid] * _sig3_wi[cid] * (sig9_wi / 45.0 / y9 - _sig3_wi[cid] / 6.0 / y3) - _uShift_9_3[cid];

					mi->Fljcenteradd(si, f);
	//				global_log->info() << "id=" << (*i).id() << ", ry=" << ry << ", ryRel=" << ryRel << ", f[1]=" << f[1] << endl;
				}
			}
		}
	}

	double u_pot;
	u_pot = _uPot_9_3[0] + domain -> getLocalUpotCompSpecific();
	domain->setLocalUpotCompSpecific(u_pot);
	for(unsigned cid = 0; cid < _nc; cid++) {
		_uPot_9_3[cid] = 0.0;
	}
} // end method calcTSLJ_9_3(...)

void Wall::calcTSLJ_10_4( ParticleContainer* partContainer, Domain* domain) {
	double regionLowCorner[3], regionHighCorner[3];

	/*! LJ-10-4 potential applied in y-direction */
	if(partContainer->getBoundingBoxMin(1)){ // if linked cell within the potential range (inside the potential's cutoff)
		for(unsigned d = 0; d < 3; d++){
			regionLowCorner[d] = partContainer->getBoundingBoxMin(d);
			regionHighCorner[d] = partContainer->getBoundingBoxMax(d);
		}

		//perform a check if the region is contained by the particleContainer???
		if(partContainer->isRegionInBoundingBox(regionLowCorner, regionHighCorner)){
			#if defined (_OPENMP)
			#pragma omp parallel shared(regionLowCorner, regionHighCorner)
			#endif
			{
				RegionParticleIterator begin = partContainer->iterateRegionBegin(regionLowCorner, regionHighCorner);
				RegionParticleIterator end = partContainer->iterateRegionEnd();

				for(RegionParticleIterator i = begin; i != end; ++i){
					//! so far for 1CLJ only, several 1CLJ-components possible
					double y, y2, y4, y5, y10, y11;
					unsigned cid = (*i).componentid();
					y = (*i).r(1) - _yOff;
					if(y < _yc){
						y2 = y * y;
						y4 = y2 * y2;
						y5 = y4 * y;
						y10 = y5 * y5;
						y11 = y10 * y;
						double f[3];
						for(unsigned d = 0; d < 3; d++) {
							f[d] = 0.0;
						}

						double sig2_wi = _sig2_wi[cid];
						double sig4_wi = _sig2_wi[cid] * _sig2_wi[cid];
						double sig5_wi = sig4_wi * _sig_wi[cid];
						double sig10_wi = sig5_wi * sig5_wi;
						double bracket = y + 0.61 * Delta;
						double bracket3 = bracket * bracket * bracket;
						double term1 = sig10_wi / y10;
						double term2 = sig4_wi / y4;
						double term3 = sig4_wi / (3 * Delta * bracket3);
						double preFactor = 2*M_PI*_eps_wi[cid]*_rhoW*sig2_wi*Delta;
						_uPot_10_4[cid] += preFactor * (2 / 5 * term1 - term2 - term3) - _uShift_10_4[cid];
						f[1] = preFactor * (4 * sig10_wi / y11 - 4 * sig4_wi / y5 - term3 * 3 / bracket);
						f[0] = 0.;
						f[2] = 0.;
						(*i).Fljcenteradd(0, f);
					}
				}
			}
		}
	}

	double u_pot;
	u_pot = _uPot_10_4[0] + domain -> getLocalUpotCompSpecific();
	domain->setLocalUpotCompSpecific(u_pot);
	for(unsigned cid = 0; cid < _nc; cid++) {
		_uPot_10_4[cid] = 0.0;
	}
} // end method calcTSLJ_10_4(...)


#include "Domain.h"
#include "longRange/Planar.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "plugins/NEMD/DistControl.h"

#include <vector>
#include <cmath>
#include <algorithm>

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"
#include "utils/FileUtils.h"


Planar::Planar(double /*cutoffT*/, double cutoffLJ, Domain* domain, DomainDecompBase* domainDecomposition,
		ParticleContainer* particleContainer, unsigned slabs, Simulation* simulation)
{
	Log::global_log->info() << "Long Range Correction for planar interfaces is used" << std::endl;
	_nStartWritingProfiles = 1;
	_nWriteFreqProfiles = 10;
	_nStopWritingProfiles = 100;
	cutoff = cutoffLJ;
	_domain = domain;
	_domainDecomposition = domainDecomposition;
	_particleContainer = particleContainer;
	_slabs = slabs;
	frequency=10;
	numComp=0;
}

Planar::~Planar() {
}

void Planar::init()
{
	//_smooth=true; // Deactivate this for transient simulations!
	_smooth=false;  //true; <-- only applicable to static density profiles

	std::vector<Component>&  components = *_simulation.getEnsemble()->getComponents();
	numComp=components.size();
	resizeExactly(numLJ, numComp);
	resizeExactly(numDipole, numComp);
	numLJSum=0;
	numDipoleSum = 0;
	resizeExactly(numLJSum2, numComp);
	resizeExactly(numDipoleSum2, numComp);
	for (unsigned i =0; i< numComp; i++){
		numLJSum2[i]=0;
		numDipoleSum2[i]=0;
	}
	for (unsigned i =0; i< numComp; i++){
		numLJ[i]=components[i].numLJcenters();
		if (components[i].numDipoles() > 0 || components[i].numCharges() != 0){
			numDipole[i]=1;
		} else{
			numDipole[i]=0;
		}
		numLJSum+=numLJ[i];
		numDipoleSum+=numDipole[i];
		for (unsigned j=i+1; j< numComp; j++){
			numLJSum2[j]+=numLJ[i];
			numDipoleSum2[j]+=numDipole[i];
		}
	}
	resizeExactly(muSquare, numDipoleSum);
	resizeExactly(uLJ, _slabs*numLJSum);
	resizeExactly(vNLJ, _slabs*numLJSum);
	resizeExactly(vTLJ, _slabs*numLJSum);
	resizeExactly(fLJ, _slabs*numLJSum);
	resizeExactly(rho_g, _slabs*numLJSum);
	resizeExactly(rho_l, _slabs*numLJSum);
	resizeExactly(fDipole, _slabs*numDipoleSum);
	resizeExactly(uDipole, _slabs*numDipoleSum);
	resizeExactly(vNDipole, _slabs*numDipoleSum);
	resizeExactly(vTDipole, _slabs*numDipoleSum);
	resizeExactly(rhoDipole, _slabs*numDipoleSum);
	resizeExactly(rhoDipoleL, _slabs*numDipoleSum);
	resizeExactly(eLong, numLJSum);

	unsigned counter=0;
	for (unsigned i =0; i< numComp; i++){		// Determination of the elongation of the Lennard-Jones sites
		for (unsigned j=0; j< components[i].numLJcenters(); j++){
			const LJcenter& ljcenteri = static_cast<const LJcenter&>(components[i].ljcenter(j));
			double dX[3];
			dX[0]=ljcenteri.rx();
			dX[1]=ljcenteri.ry();
			dX[2]=ljcenteri.rz();
			for (unsigned d=0; d<3; d++){
				dX[d]*=dX[d];
			}
			eLong[counter]=sqrt(dX[0]+dX[1]+dX[2]);
			counter++;
		}
	}
	for (unsigned i=0; i< _slabs*numLJSum; i++){
		rho_g[i]=0;
	}
	for (unsigned i=0; i < _slabs*numDipoleSum; i++){
		rhoDipole[i]=0;
	}

	unsigned dpcount=0;
	for (unsigned i=0; i<numComp; i++){
	//	for (unsigned j=0; j<numDipole[i]; j++){
		if (numDipole[i] != 0){
			double chargeBalance[3];
			for (unsigned d=0; d < 3; d++) chargeBalance[d] = 0;
			for (unsigned si=0; si < components[i].numCharges(); si++){
				double tq = components[i].charge(si).q();
				for (unsigned d=0; d < 3; d++) chargeBalance[d] += tq * components[i].charge(si).r()[d];
			}
			for (unsigned si=0; si < components[i].numDipoles(); si++){
				double tmy = components[i].dipole(si).absMy();
				for (unsigned d=0; d < 3; d++) chargeBalance[d] += tmy * components[i].dipole(si).e()[d];
			}

			// const Dipole& dipolej = static_cast<const Dipole&> (components[i].dipole(j));
			double my2 = 0;
			for (unsigned d=0; d < 3; d++)	my2 += chargeBalance[d]*chargeBalance[d];
			muSquare[dpcount]=my2;//dipolej.absMy()*dipolej.absMy();
			dpcount++;
#ifndef NDEBUG
			std::cout << dpcount << "\tmu: " << muSquare[dpcount-1] << std::endl;
#endif
		}
	}

	ymax=_domain->getGlobalLength(1);
	boxlength[0]=_domain->getGlobalLength(0);
	boxlength[2]=_domain->getGlobalLength(2);
	cutoff_slabs=cutoff*_slabs/ymax; // Number of slabs to cutoff
	delta=ymax/_slabs;
	V=ymax*_domain->getGlobalLength(0)*_domain->getGlobalLength(2); // Volume

	sint=_slabs;
	simstep = 0;

	temp=_domain->getTargetTemperature(0);

	if (_region.refPosID[0]+_region.refPosID[1] > 0) {
		_subject = getSubject();
		if (nullptr != _subject) {
			this->update(_subject);
			_subject->registerObserver(this);
			Log::global_log->info() << "Long Range Correction: Subject registered" << std::endl;
		} else {
			Log::global_log->error() << "Long Range Correction: Initialization of plugin DistControl is needed before! Program exit..." << std::endl;
			Simulation::exit(-1);
		}
	}
}

void Planar::readXML(XMLfileUnits& xmlconfig)
{
	Log::global_log->info() << "[Long Range Correction] Reading xml config" << std::endl;

	xmlconfig.getNodeValue("slabs", _slabs);
	xmlconfig.getNodeValue("smooth", _smooth);
	xmlconfig.getNodeValue("frequency", frequency);

	std::string strVal;
	_region.refPosID[0] = 0; _region.refPosID[1] = 0;

	// In some cases, like for density gradients in the bulk, the planar LRC may be wrong, since it was mainly developed for interfaces.
	// Therefore, a region can be specified wherein the LRC for the force is applied.
	// In order to still get correct state values in the bulks, the LRC for the virial and potential energy is also applied outside of the region.
	if (xmlconfig.getNodeValue("region/left", _region.refPos[0]) && xmlconfig.getNodeValue("region/right", strVal)) {
		// accept "box" as input
		_region.refPos[1] = (strVal == "box") ? _domain->getGlobalLength(1) : atof(strVal.c_str() );
		_region.actPos[0] = _region.refPos[0];
		_region.actPos[1] = _region.refPos[1];
		// read reference coords IDs
		xmlconfig.getNodeValue("region/left@refcoordsID", _region.refPosID[0] );
		xmlconfig.getNodeValue("region/right@refcoordsID", _region.refPosID[1] );
	} else {
		_region.refPos[0] = _region.actPos[0] = 0.0;
		_region.refPos[1] = _region.actPos[1] = _domain->getGlobalLength(1);
	}

	Log::global_log->info() << "Long Range Correction: using " << _slabs << " slabs for profiles to calculate LRC." << std::endl;
	Log::global_log->info() << "Long Range Correction: sampling profiles every " << frequency << "th simstep." << std::endl;
	Log::global_log->info() << "Long Range Correction: profiles are smoothed (averaged over time): " << std::boolalpha << _smooth << std::endl;
	Log::global_log->info() << "Long Range Correction: force corrections are applied to particles within yRef = " << _region.refPos[0] << " (refID: " << _region.refPosID[0] << ") and " << _region.refPos[1] << " (refID: " << _region.refPosID[1] << ")" << std::endl;
	Log::global_log->info() << "Long Range Correction: pot. energy and virial corrections are applied within whole domain" << std::endl;

	// write control
	bool bRet1 = xmlconfig.getNodeValue("writecontrol/start", _nStartWritingProfiles);
	bool bRet2 = xmlconfig.getNodeValue("writecontrol/frequency", _nWriteFreqProfiles);
	bool bRet3 = xmlconfig.getNodeValue("writecontrol/stop", _nStopWritingProfiles);
	if(_nWriteFreqProfiles < 1)
	{
		Log::global_log->error() << "Long Range Correction: Write frequency < 1! Programm exit ..." << std::endl;
		Simulation::exit(-1);
	}
	if(_nStopWritingProfiles <= _nStartWritingProfiles)
	{
		Log::global_log->error() << "Long Range Correction: Writing profiles 'stop' <= 'start'! Programm exit ..." << std::endl;
		Simulation::exit(-1);
	}
	bool bInputIsValid = (bRet1 && bRet2 && bRet3);
	if(true == bInputIsValid)
	{
		Log::global_log->info() << "Long Range Correction->writecontrol: Start writing profiles at simstep: " << _nStartWritingProfiles << std::endl;
		Log::global_log->info() << "Long Range Correction->writecontrol: Writing profiles with frequency: " << _nWriteFreqProfiles << std::endl;
		Log::global_log->info() << "Long Range Correction->writecontrol: Stop writing profiles at simstep: " << _nStopWritingProfiles << std::endl;
	}
	else
	{
		Log::global_log->error() << "Long Range Correction: Write control parameters not valid! Programm exit ..." << std::endl;
		Simulation::exit(-1);
	}
}

void Planar::calculateLongRange() {

	if (_smooth){
		const double delta_inv = 1.0 / delta;
		const double slabsPerV = _slabs / V;

		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		for(auto tempMol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid(); ++tempMol){
			unsigned cid=tempMol->componentid();

			for (unsigned i=0; i<numLJ[cid]; i++){
				int loc=(tempMol->ljcenter_d_abs(i)[1]) * delta_inv;
				if (loc < 0){
					loc=loc+_slabs;
				}
				else if (loc > sint-1){
					loc=loc-_slabs;
				}
				const int index = loc + _slabs * (i + numLJSum2[cid]);

				#if defined(_OPENMP)
				#pragma omp atomic
				#endif
				rho_g[index] += slabsPerV;
			}
			if (numDipole[cid] != 0){
				int loc=tempMol->r(1) * delta_inv;

				const int index = loc + _slabs * (numDipoleSum2[cid]);

				#if defined(_OPENMP)
				#pragma omp atomic
				#endif
				rhoDipole[index] += slabsPerV;
			}
		}
	}
	if (simstep % frequency == 0){	// The Density Profile is only calculated once in 10 simulation steps

		std::fill(rho_l.begin(), rho_l.end(), 0.0);
		std::fill(uLJ.begin(), uLJ.end(), 0.0);
		std::fill(vNLJ.begin(), vNLJ.end(), 0.0);
		std::fill(vTLJ.begin(), vTLJ.end(), 0.0);
		std::fill(fLJ.begin(), fLJ.end(), 0.0);

		std::fill(fDipole.begin(), fDipole.end(), 0.0);
		std::fill(uDipole.begin(), uDipole.end(), 0.0);
		std::fill(vNDipole.begin(), vNDipole.end(), 0.0);
		std::fill(vTDipole.begin(), vTDipole.end(), 0.0);
		std::fill(rhoDipoleL.begin(), rhoDipoleL.end(), 0.0);

		// Calculation of the density profile for s slabs
		if (!_smooth){
			const double delta_inv = 1.0 / delta;
			const double slabsPerV = _slabs / V;

			#if defined(_OPENMP)
			#pragma omp parallel
			#endif
			for(auto tempMol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid(); ++tempMol){
				unsigned cid=tempMol->componentid();

				for (unsigned i=0; i<numLJ[cid]; i++){
					int loc=(tempMol->ljcenter_d_abs(i)[1]) * delta_inv;
					if (loc < 0){
						loc=loc+_slabs;
					}
					else if (loc > sint-1){
						loc=loc-_slabs;
					}
					const int index = loc+_slabs*(i+numLJSum2[cid]);

					#if defined(_OPENMP)
					#pragma omp atomic
					#endif
					rho_l[index] += slabsPerV;
				}
				if (numDipole[cid] != 0){
					int loc=tempMol->r(1) * delta_inv;
					const int index = loc+_slabs*numDipoleSum2[cid];

					#if defined(_OPENMP)
					#pragma omp atomic
					#endif
					rhoDipoleL[index] += slabsPerV;
				}
			}
		}
		else{
			for (unsigned i=0; i<_slabs*numLJSum; i++){
				rho_l[i]=rho_g[i]/(simstep+1);
			}
			for (unsigned i=0; i<_slabs*numDipoleSum; i++){
				rhoDipoleL[i]=rhoDipole[i]/(simstep+1);
			}
		}

		// Distribution of the Density Profile to every node
		_domainDecomposition->collCommInit(_slabs*(numLJSum+numDipoleSum));
		for (unsigned i=0; i < _slabs*numLJSum; i++) {
			_domainDecomposition->collCommAppendDouble(rho_l[i]);
		}
		for (unsigned i=0; i< _slabs*numDipoleSum; i++){
			_domainDecomposition->collCommAppendDouble(rhoDipoleL[i]);
		}

		_domainDecomposition->collCommAllreduceSum();
		for (unsigned i=0; i < _slabs*numLJSum; i++) {
			rho_l[i] = _domainDecomposition->collCommGetDouble();
		}
		for (unsigned i=0; i< _slabs*numDipoleSum; i++){
			rhoDipoleL[i] = _domainDecomposition->collCommGetDouble();
		}

		_domainDecomposition->collCommFinalize();

		for (unsigned ci = 0; ci < numComp; ++ci){
			for (unsigned cj = 0; cj < numComp; ++cj){
				ParaStrm& params = _domain->getComp2Params()(ci,cj);
				params.reset_read();
				for (unsigned si = 0; si < numLJ[ci]; ++si) { // Long Range Correction for Lennard-Jones sites
					for (unsigned sj = 0; sj < numLJ[cj]; ++sj) {
						double eps24;
						double sig2;
						double shift6;
						double eps;
						params >> eps24;
						params >> sig2;
						params >> shift6;
						sig2=sqrt(sig2);
						eps=eps24/24;
						if (eLong[numLJSum2[ci]+si] ==0 && eLong[numLJSum2[cj]+sj] == 0){
							centerCenter(sig2,eps,ci,cj,si,sj);
						}
						else if (eLong[numLJSum2[ci]+si] ==0 || eLong[numLJSum2[cj]+sj] == 0){
							centerSite(sig2,eps,ci,cj,si,sj);
						}
						else{
							siteSite(sig2,eps,ci,cj,si,sj);
						}
					}
				}

				for (unsigned si=0; si< numDipole[ci]; si++){	//Long Range Correction for Dipoles
					for (unsigned sj=0; sj< numDipole[cj]; sj++){
						 dipoleDipole(ci,cj,si,sj);
					}
				}
			}
	 	}

		// Distribution of the Force, Energy and Virial to every Node
		_domainDecomposition->collCommInit(_slabs*(4*numLJSum+4*numDipoleSum));
		for (unsigned i=0; i<_slabs*numLJSum; i++){
			_domainDecomposition->collCommAppendDouble(uLJ[i]);
			_domainDecomposition->collCommAppendDouble(vNLJ[i]);
			_domainDecomposition->collCommAppendDouble(vTLJ[i]);
			_domainDecomposition->collCommAppendDouble(fLJ[i]);
		}
		for (unsigned i=0; i<_slabs*numDipoleSum; i++){
			_domainDecomposition->collCommAppendDouble(uDipole[i]);
			_domainDecomposition->collCommAppendDouble(fDipole[i]);
			_domainDecomposition->collCommAppendDouble(vNDipole[i]);
			_domainDecomposition->collCommAppendDouble(vTDipole[i]);
		}
		_domainDecomposition->collCommAllreduceSum();
		for (unsigned i=0; i<_slabs*numLJSum; i++){
			uLJ[i]=_domainDecomposition->collCommGetDouble();
			vNLJ[i]=_domainDecomposition->collCommGetDouble();
			vTLJ[i]=_domainDecomposition->collCommGetDouble();
			fLJ[i]=_domainDecomposition->collCommGetDouble();
		}
		for (unsigned i=0; i<_slabs*numDipoleSum; i++){
		      uDipole[i]=_domainDecomposition->collCommGetDouble();
		      fDipole[i]=_domainDecomposition->collCommGetDouble();
		      vNDipole[i]=_domainDecomposition->collCommGetDouble();
		      vTDipole[i]=_domainDecomposition->collCommGetDouble();
		}
		_domainDecomposition->collCommFinalize();
	}

	// Adding the Force to the Molecules; this is done in every timestep
	const double delta_inv = 1.0 / delta;

	double Upot_c=0.0;
	double Virial_c=0.0; // Correction used for the Pressure Calculation

#ifndef NDEBUG
	std::cout << "Long Range Correction: Using region boundaries: " << _region.actPos[0] << " and " << _region.actPos[1] << std::endl;
#endif

	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:Upot_c, Virial_c)
	#endif
	for (auto tempMol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid(); ++tempMol) {

		unsigned cid = tempMol->componentid();

		for (unsigned i=0; i<numLJ[cid]; i++){
			int loc=(tempMol->ljcenter_d_abs(i)[1]) * delta_inv;
			if (loc < 0){
				loc=loc+_slabs;
			}
			else if (loc > sint-1){
				loc=loc-_slabs;
			}
			double Fa[3]={0.0, 0.0, 0.0};
			const int index = loc + i * _slabs + _slabs * numLJSum2[cid];
			Fa[1] = fLJ[index];
			Upot_c += uLJ[index];
			Virial_c += 2 * vTLJ[index] + vNLJ[index];
			double Via[3] = {0.0, 0.0, 0.0};
			Via[0] = vTLJ[index];
			Via[1] = vNLJ[index];
			Via[2] = vTLJ[index];
			if ((tempMol->r(1) >= _region.actPos[0]) && (tempMol->r(1) <= _region.actPos[1])) {
				tempMol->Fljcenteradd(i, Fa);
			}
			tempMol->Viadd(Via);
//			tempMol->Uadd(uLJ[loc+i*s+_slabs*numLJSum2[cid]]);      // Storing potential energy onto the molecules is currently not implemented!
		}
		if (numDipole[cid] != 0){
			int loc = tempMol->r(1) * delta_inv;
			double Fa[3] = {0.0, 0.0, 0.0};
			const int index = loc + _slabs * numDipoleSum2[cid];
			Fa[1] = fDipole[index];
			Upot_c += uDipole[index];
			Virial_c += 2 * vTDipole[index] + vNDipole[index];
			double Via[3] = {0.0, 0.0, 0.0};
			Via[0] = vTDipole[index];
			Via[1] = vNDipole[index];
			Via[2] = vTDipole[index];
			if ((tempMol->r(1) >= _region.actPos[0]) && (tempMol->r(1) <= _region.actPos[1])) {
				tempMol->Fadd(Fa); // Force is stored on the center of mass of the molecule!
			}
			tempMol->Viadd(Via);
//			tempMol->Uadd(uDipole[loc+i*_slabs+_slabs*numDipoleSum2[cid]]); // Storing potential energy onto the molecules is currently not implemented!
		}
	}

	// Summation of the correction terms
	auto collComm = makeCollCommObjAllreduceAdd(_domainDecomposition->getCommunicator(), Upot_c, Virial_c);
	collComm.communicate();
	collComm.get(Upot_c, Virial_c);

	// Setting the Energy and Virial correction
	_domain->setUpotCorr(Upot_c);
	_domain->setVirialCorr(Virial_c);

	simstep++;
}


void Planar::centerCenter(double sig, double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj){
	double rc=sig/cutoff;
	double rc2=rc*rc;
	const double rc2_inv = 1.0 / rc2;
	double rc6=rc2*rc2*rc2;
	double rc12=rc6*rc6;
	double r,r2,r6,r12;
	double rhoI,rhoJ;
	double termU = 4*3.1416*delta*eps*sig*sig;
	double termF = 8*3.1416*delta*eps*sig;
	double termVN = 4*3.1416*delta*eps*sig*sig;
	double termVT = 2*3.1416*delta*eps*sig*sig;
	const unsigned slabsHalf = _slabs / 2;
	for (unsigned i=_domainDecomposition->getRank(); i<slabsHalf; i+=_domainDecomposition->getNumProcs()){
		const int index_i = i+si*_slabs+_slabs*numLJSum2[ci];

		rhoI=rho_l[index_i];
		vTLJ[index_i] += termVT*rhoI*(6*rc12*0.2-3*rc6*0.5)*rc2_inv;
		uLJ[index_i] += termU*rhoI*(rc12*0.2-rc6*0.5)*rc2_inv;
		for (unsigned j=i+1; j<i+slabsHalf; j++){
			const int index_j = j+sj*_slabs+_slabs*numLJSum2[cj];

			r=sig/((j-i)*delta);
			rhoJ=rho_l[index_j];
			if (j> i+cutoff_slabs){
				r2=r*r;
				r6=r2*r2*r2;
				r12=r6*r6;
				vTLJ[index_i]+=termVT*rhoJ*(r12*0.2-r6*0.5)/r2;
				vTLJ[index_j]+=termVT*rhoI*(r12*0.2-r6*0.5)/r2;
			}
			else{
				r2=rc2;
				r6=rc6;
				r12=rc12;
				vTLJ[index_j]+=termVT*rhoI*(r12*0.2*(6/r2-5/(r*r))-r6*0.5*(3/r2-2/(r*r)));
				vTLJ[index_i]+=termVT*rhoJ*(r12*0.2*(6/r2-5/(r*r))-r6*0.5*(3/r2-2/(r*r)));
			}
			uLJ[index_i] += termU*rhoJ*(r12*0.2-r6*0.5)/r2;
			uLJ[index_j] += termU*rhoI*(r12*0.2-r6*0.5)/r2;
			vNLJ[index_i] += termVN*rhoJ*(r12-r6)/(r*r);
			vNLJ[index_j] += termVN*rhoI*(r12-r6)/(r*r);
			fLJ[index_i] += -termF*rhoJ*(r12-r6)/r;
			fLJ[index_j] += termF*rhoI*(r12-r6)/r;
		}
		// Calculation of the Periodic boundary
		for (unsigned j=slabsHalf+i; j<_slabs; j++){
			const int index_j = j+sj*_slabs+_slabs*numLJSum2[cj];

			r=sig/((_slabs-j+i)*delta);
			rhoJ=rho_l[index_j];
			if (j <_slabs-cutoff_slabs+i){
				r2=r*r;
				r6=r2*r2*r2;
				r12=r6*r6;
				vTLJ[index_i] += termVT*rhoJ*(r12*0.2-r6*0.5)/r2;
				vTLJ[index_j] += termVT*rhoI*(r12*0.2-r6*0.5)/r2;
			}
			else{
				r2=rc2;
				r6=rc6;
				r12=rc12;
				vTLJ[index_j] += termVT*rhoI*(r12*0.2*(6/r2-5/(r*r))-r6*0.5*(3/r2-2/(r*r)));
				vTLJ[index_i] += termVT*rhoJ*(r12*0.2*(6/r2-5/(r*r))-r6*0.5*(3/r2-2/(r*r)));
			}
			uLJ[index_i] += termU*rhoJ*(r12*0.2-r6*0.5)/r2;
			uLJ[index_j] += termU*rhoI*(r12*0.2-r6*0.5)/r2;
			vNLJ[index_i] += termVN*rhoJ*(r12-r6)/(r*r);
			vNLJ[index_j] += termVN*rhoI*(r12-r6)/(r*r);
			fLJ[index_i] += termF*rhoJ*(r12-r6)/r;
			fLJ[index_j] += -termF*rhoI*(r12-r6)/r;
		}
	}

	// Calculation of the Forces on the slabs of the right hand side
	for (unsigned i=slabsHalf+_domainDecomposition->getRank(); i<_slabs; i+=_domainDecomposition->getNumProcs()){
		const int index_i = i+si*_slabs+_slabs*numLJSum2[ci];
		const int index_i_cj = i+si*_slabs+_slabs*numLJSum2[cj]; // TODO: check - this is really supposed to be a mix of the i and j variables?

		rhoI=rho_l[index_i];
		vTLJ[index_i_cj] += termVT*rhoI*(6*rc12*0.2-3*rc6*0.5)*rc2_inv;
		uLJ[index_i_cj] += termU*rhoI*(rc12*0.2-rc6*0.5)*rc2_inv;
		for (unsigned j=i+1; j<_slabs; j++) {
			const int index_j = j+sj*_slabs+_slabs*numLJSum2[cj];

			r=sig/((j-i)*delta);
			rhoJ=rho_l[index_j];
			if (j> i+cutoff_slabs){
				r2=r*r;
				r6=r2*r2*r2;
				r12=r6*r6;
				vTLJ[index_i] += termVT*rhoJ*(r12*0.2-r6*0.5)/r2;
				vTLJ[index_j] += termVT*rhoI*(r12*0.2-r6*0.5)/r2;
			}
			else{
				r2=rc2;
				r6=rc6;
				r12=rc12;
				vTLJ[index_j] += termVT*rhoI*(r12*0.2*(6/r2-5/(r*r))-r6*0.5*(3/r2-2/(r*r)));
				vTLJ[index_i] += termVT*rhoJ*(r12*0.2*(6/r2-5/(r*r))-r6*0.5*(3/r2-2/(r*r)));
			}
			uLJ[index_i] += termU*rhoJ*(r12*0.2-r6*0.5)/r2;
			uLJ[index_j] += termU*rhoI*(r12*0.2-r6*0.5)/r2;
			vNLJ[index_i] += termVN*rhoJ*(r12-r6)/(r*r);
			vNLJ[index_j] += termVN*rhoI*(r12-r6)/(r*r);
			fLJ[index_i] += -termF*rhoJ*(r12-r6)/r;
			fLJ[index_j] += termF*rhoI*(r12-r6)/r;
		}
	}
}

void Planar::centerSite(double sig, double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj){
	double sig2=sig*sig;
	double sig3=sig2*sig;
	double t = eLong[numLJSum2[ci]+si] + eLong[numLJSum2[cj]+sj]; // one of them is equal to zero.
	double rcPt=sig/(cutoff+t);
	double rcPt3=rcPt*rcPt*rcPt;
	double rcPt4=rcPt3*rcPt;
	double rcPt9=rcPt3*rcPt3*rcPt3;
	double rcPt10=rcPt9*rcPt;
	double rcMt=sig/(cutoff-t);
	double rcMt3=rcMt*rcMt*rcMt;
	double rcMt4=rcMt3*rcMt;
	double rcMt9=rcMt3*rcMt3*rcMt3;
	double rcMt10=rcMt9*rcMt;
	double rc=cutoff;
	double rc2=rc*rc;
	double termURC=-2*3.1416*eps*delta*sig3/(3*t)*((rcPt9-rcMt9)/15-(rcPt3-rcMt3)/2);
	double termFRC=-2*3.1416*eps*delta*sig2/(t*rc)*((rcPt10-rcMt10)/5-(rcPt4-rcMt4)/2);
	double termVNRC=termFRC/2;
	double termVTRC1=-3.1416*eps*delta*sig2/(2*t*rc)*((rcPt10-rcMt10)/5-(rcPt4-rcMt4)/2);
	double termVTRC2=termURC/2;
	double r,r2;
	double rhoI,rhoJ;
	for (unsigned i=_domainDecomposition->getRank(); i<_slabs/2; i+=_domainDecomposition->getNumProcs()){
		const int index_i = i+si*_slabs+_slabs*numLJSum2[ci];

		rhoI=rho_l[index_i];

		const int index_i_sj_cj = i+sj*_slabs+_slabs*numLJSum2[cj];
		rhoJ=rho_l[index_i_sj_cj]; // TODO: i, sj, cj? appears a few more times in other functions

		vTLJ[index_i]+=rhoJ*(termVTRC1*rc2+termVTRC2);
		uLJ[index_i]+=rhoJ*termURC;
		for (unsigned j=i+1; j<i+_slabs/2; j++){
			const int index_j = j+sj*_slabs+_slabs*numLJSum2[cj];

			r=(j-i)*delta; // xi in Werth2014
			r2=r*r;
			rhoJ=rho_l[index_j];
			if (j> i+cutoff_slabs){
				double rPt=sig/(r+t);
				double rPt3=rPt*rPt*rPt;
				double rPt4=rPt3*rPt;
				double rPt9=rPt3*rPt3*rPt3;
				double rPt10=rPt9*rPt;
				double rMt=sig/(r-t);
				double rMt3=rMt*rMt*rMt;
				double rMt4=rMt3*rMt;
				double rMt9=rMt3*rMt3*rMt3;
				double rMt10=rMt9*rMt;
				double termU=-2*3.1416*eps*delta*sig3/(3*t)*((rPt9-rMt9)/15-(rPt3-rMt3)/2);
				double termF=-2*3.1416*eps*delta*sig2/(t*r)*((rPt10-rMt10)/5-(rPt4-rMt4)/2);
				double termVN=termF/2;
				double termVT2=termU/2;
				uLJ[index_i]+=rhoJ*termU;
				uLJ[index_j]+=rhoI*termU;
				vNLJ[index_i]+=rhoJ*termVN*r2;
				vNLJ[index_j]+=rhoI*termVN*r2;
				vTLJ[index_j]+=rhoI*termVT2;
				vTLJ[index_i]+=rhoJ*termVT2;
				fLJ[index_i]+=-rhoJ*termF*r;
				fLJ[index_j]+=rhoI*termF*r;
			}
			else{
				uLJ[index_i]+=rhoJ*termURC;
				uLJ[index_j]+=rhoI*termURC;
				vNLJ[index_i]+=rhoJ*termVNRC*r2;
				vNLJ[index_j]+=rhoI*termVNRC*r2;
				vTLJ[index_j]+=rhoI*(termVTRC1*(rc2-r2)+termVTRC2);
				vTLJ[index_i]+=rhoJ*(termVTRC1*(rc2-r2)+termVTRC2);
				fLJ[index_i]+=-rhoJ*termFRC*r;
				fLJ[index_j]+=rhoI*termFRC*r;
			}
		}
		// Calculation of the Periodic boundary
		for (unsigned j=_slabs/2+i; j<_slabs; j++) {
			const int index_j = j+sj*_slabs+_slabs*numLJSum2[cj];

			r=(_slabs-j+i)*delta;
			r2=r*r;
			rhoJ=rho_l[index_j];
			if (j <_slabs-cutoff_slabs+i){
				double rPt=sig/(r+t);
				double rPt3=rPt*rPt*rPt;
				double rPt4=rPt3*rPt;
				double rPt9=rPt3*rPt3*rPt3;
				double rPt10=rPt9*rPt;
				double rMt=sig/(r-t);
				double rMt3=rMt*rMt*rMt;
				double rMt4=rMt3*rMt;
				double rMt9=rMt3*rMt3*rMt3;
				double rMt10=rMt9*rMt;
				double termU=-2*3.1416*eps*delta*sig3/(3*t)*((rPt9-rMt9)/15-(rPt3-rMt3)/2);
				double termF=-2*3.1416*eps*delta*sig2/(t*r)*((rPt10-rMt10)/5-(rPt4-rMt4)/2);
				double termVN=termF/2;
				double termVT2=termU/2;
				uLJ[index_i]+=rhoJ*termU;
				uLJ[index_j]+=rhoI*termU;
				vNLJ[index_i]+=rhoJ*termVN*r2;
				vNLJ[index_j]+=rhoI*termVN*r2;
				vTLJ[index_i]+=rhoJ*termVT2;
				vTLJ[index_j]+=rhoI*termVT2;
				fLJ[index_i]+=rhoJ*termF*r;
				fLJ[index_j]+=-rhoI*termF*r;
			}
			else{
				uLJ[index_i]+=rhoJ*termURC;
				uLJ[index_j]+=rhoI*termURC;
				vNLJ[index_i]+=rhoJ*termVNRC*r2;
				vNLJ[index_j]+=rhoI*termVNRC*r2;
				vTLJ[index_i]+=rhoJ*(termVTRC1*(rc2-r2)+termVTRC2);
				vTLJ[index_j]+=rhoI*(termVTRC1*(rc2-r2)+termVTRC2);
				fLJ[index_i]+=rhoJ*termFRC*r;
				fLJ[index_j]+=-rhoI*termFRC*r;
			}
		}
	}

	// Calculation of the Forces on the slabs of the right hand side
	for (unsigned i=_slabs/2+_domainDecomposition->getRank(); i<_slabs; i+=_domainDecomposition->getNumProcs()){
		const int index_i = i+si*_slabs+_slabs*numLJSum2[ci];

		rhoI=rho_l[index_i];

		const int index_i_sj_cj = i+sj*_slabs+_slabs*numLJSum2[cj];
		rhoJ=rho_l[index_i_sj_cj];

		vTLJ[index_i]+=rhoJ*(termVTRC1*rc2+termVTRC2);
		uLJ[index_i]+=rhoJ*termURC;
		for (unsigned j=i+1; j<_slabs; j++){
			const int index_j = j+sj*_slabs+_slabs*numLJSum2[cj];

			r=(j-i)*delta;
			r2=r*r;
			rhoJ=rho_l[index_j];
			if (j> i+cutoff_slabs){
				double rPt=sig/(r+t);
				double rPt3=rPt*rPt*rPt;
				double rPt4=rPt3*rPt;
				double rPt9=rPt3*rPt3*rPt3;
				double rPt10=rPt9*rPt;
				double rMt=sig/(r-t);
				double rMt3=rMt*rMt*rMt;
				double rMt4=rMt3*rMt;
				double rMt9=rMt3*rMt3*rMt3;
				double rMt10=rMt9*rMt;
				double termU=-2*3.1416*eps*delta*sig3/(3*t)*((rPt9-rMt9)/15-(rPt3-rMt3)/2);
				double termF=-2*3.1416*eps*delta*sig2/(t*r)*((rPt10-rMt10)/5-(rPt4-rMt4)/2);
				double termVN=termF/2;
				double termVT2=termU/2;
				uLJ[index_i]+=rhoJ*termU;
				uLJ[index_j]+=rhoI*termU;
				vNLJ[index_i]+=rhoJ*termVN*r2;
				vNLJ[index_j]+=rhoI*termVN*r2;
				vTLJ[index_j]+=rhoI*termVT2;
				vTLJ[index_i]+=rhoJ*termVT2;
				fLJ[index_i]+=-rhoJ*termF*r;
				fLJ[index_j]+=rhoI*termF*r;
			}
			else{
				uLJ[index_i]+=rhoJ*termURC;
				uLJ[index_j]+=rhoI*termURC;
				vNLJ[index_i]+=rhoJ*termVNRC*r2;
				vNLJ[index_j]+=rhoI*termVNRC*r2;
				vTLJ[index_j]+=rhoI*(termVTRC1*(rc2-r2)+termVTRC2);
				vTLJ[index_i]+=rhoJ*(termVTRC1*(rc2-r2)+termVTRC2);
				fLJ[index_i]+=-rhoJ*termFRC*r;
				fLJ[index_j]+=rhoI*termFRC*r;
			}
		}
	}
}

void Planar::siteSite(double sig, double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj){
	double sig2=sig*sig;
	double sig3=sig2*sig;
	double sig4=sig2*sig2;
	double t1 = eLong[numLJSum2[ci]+si];
	double t2 = eLong[numLJSum2[cj]+sj];
	double tP = t1 + t2; // tau+
	double tM = t1 - t2; // tau-
	double rcPtP=sig/(cutoff+tP);
	double rcPtP2=rcPtP*rcPtP;
	double rcPtP3=rcPtP2*rcPtP;
	double rcPtP8=rcPtP2*rcPtP3*rcPtP3;
	double rcPtP9=rcPtP8*rcPtP;
	double rcPtM=sig/(cutoff+tM);
	double rcPtM2=rcPtM*rcPtM;
	double rcPtM3=rcPtM2*rcPtM;
	double rcPtM8=rcPtM2*rcPtM3*rcPtM3;
	double rcPtM9=rcPtM8*rcPtM;
	double rcMtP=sig/(cutoff-tP);
	double rcMtP2=rcMtP*rcMtP;
	double rcMtP3=rcMtP2*rcMtP;
	double rcMtP8=rcMtP2*rcMtP3*rcMtP3;
	double rcMtP9=rcMtP8*rcMtP;
	double rcMtM=sig/(cutoff-tM);
	double rcMtM2=rcMtM*rcMtM;
	double rcMtM3=rcMtM2*rcMtM;
	double rcMtM8=rcMtM2*rcMtM3*rcMtM3;
	double rcMtM9=rcMtM8*rcMtM;
	double rc=cutoff;
	double rc2=rc*rc;
	double termURC=3.1416*eps*delta*sig4/(12*t1*t2)*((rcPtP8-rcPtM8-rcMtM8+rcMtP8)/30-(rcPtP2-rcPtM2-rcMtM2+rcMtP2));
	double termFRC=3.1416*eps*delta*sig3/(3*t1*t2*rc)*((rcPtP9-rcPtM9-rcMtM9+rcMtP9)/15-(rcPtP3-rcPtM3-rcMtM3+rcMtP3)/2);
	double termVNRC=termFRC/2;
	double termVTRC1=3.1416*eps*delta*sig3/(4*t1*t2*rc)*((rcPtP9-rcPtM9-rcMtM9+rcMtP9)/45-(rcPtP3-rcPtM3-rcMtM3+rcMtP3)/6);
	double termVTRC2=termURC/2;
	double r,r2;
	double rhoI,rhoJ;
	for (unsigned i=_domainDecomposition->getRank(); i<_slabs/2; i+=_domainDecomposition->getNumProcs()) {
		const int index_i = i+si*_slabs+_slabs*numLJSum2[ci];

		rhoI=rho_l[index_i];

		const int index_i_sj_cj = i+sj*_slabs+_slabs*numLJSum2[cj];
		rhoJ=rho_l[index_i_sj_cj];

		vTLJ[index_i]+=rhoJ*(termVTRC1*rc2+termVTRC2);
		uLJ[index_i]+=rhoJ*termURC;
		for (unsigned j=i+1; j<i+_slabs/2; j++){
			const int index_j = j+sj*_slabs+_slabs*numLJSum2[cj];

			r=(j-i)*delta;
			r2=r*r;
			rhoJ=rho_l[index_j];
			if (j> i+cutoff_slabs){
				double rPtP=sig/(r+tP);
				double rPtP2=rPtP*rPtP;
				double rPtP3=rPtP2*rPtP;
				double rPtP8=rPtP2*rPtP3*rPtP3;
				double rPtP9=rPtP8*rPtP;
				double rPtM=sig/(r+tM);
				double rPtM2=rPtM*rPtM;
				double rPtM3=rPtM2*rPtM;
				double rPtM8=rPtM2*rPtM3*rPtM3;
				double rPtM9=rPtM8*rPtM;
				double rMtP=sig/(r-tP);
				double rMtP2=rMtP*rMtP;
				double rMtP3=rMtP2*rMtP;
				double rMtP8=rMtP2*rMtP3*rMtP3;
				double rMtP9=rMtP8*rMtP;
				double rMtM=sig/(r-tM);
				double rMtM2=rMtM*rMtM;
				double rMtM3=rMtM2*rMtM;
				double rMtM8=rMtM2*rMtM3*rMtM3;
				double rMtM9=rMtM8*rMtM;
				double termU=3.1416*eps*delta*sig4/(12*t1*t2)*((rPtP8-rPtM8-rMtM8+rMtP8)/30-(rPtP2-rPtM2-rMtM2+rMtP2));
				double termF=3.1416*eps*delta*sig3/(3*t1*t2*r)*((rPtP9-rPtM9-rMtM9+rMtP9)/15-(rPtP3-rPtM3-rMtM3+rMtP3)/2);
				double termVN=termF/2;
				double termVT2=termU/2;
				uLJ[index_i]+=rhoJ*termU;
				uLJ[index_j]+=rhoI*termU;
				vNLJ[index_i]+=rhoJ*termVN*r2;
				vNLJ[index_j]+=rhoI*termVN*r2;
				vTLJ[index_j]+=rhoI*termVT2;
				vTLJ[index_i]+=rhoJ*termVT2;
				fLJ[index_i]+=-rhoJ*termF*r;
				fLJ[index_j]+=rhoI*termF*r;
			}
			else{
				uLJ[index_i]+=rhoJ*termURC;
				uLJ[index_j]+=rhoI*termURC;
				vNLJ[index_i]+=rhoJ*termVNRC*r2;
				vNLJ[index_j]+=rhoI*termVNRC*r2;
				vTLJ[index_j]+=rhoI*(termVTRC1*(rc2-r2)+termVTRC2);
				vTLJ[index_i]+=rhoJ*(termVTRC1*(rc2-r2)+termVTRC2);
				fLJ[index_i]+=-rhoJ*termFRC*r;
				fLJ[index_j]+=rhoI*termFRC*r;
			}
		}
		// Calculation of the Periodic boundary
		for (unsigned j=_slabs/2+i; j<_slabs; j++) {
			const int index_j = j+sj*_slabs+_slabs*numLJSum2[cj];

			r=(_slabs-j+i)*delta;
			r2=r*r;
			rhoJ=rho_l[index_j];
			if (j <_slabs-cutoff_slabs+i){
				double rPtP=sig/(r+tP);
				double rPtP2=rPtP*rPtP;
				double rPtP3=rPtP2*rPtP;
				double rPtP8=rPtP2*rPtP3*rPtP3;
				double rPtP9=rPtP8*rPtP;
				double rPtM=sig/(r+tM);
				double rPtM2=rPtM*rPtM;
				double rPtM3=rPtM2*rPtM;
				double rPtM8=rPtM2*rPtM3*rPtM3;
				double rPtM9=rPtM8*rPtM;
				double rMtP=sig/(r-tP);
				double rMtP2=rMtP*rMtP;
				double rMtP3=rMtP2*rMtP;
				double rMtP8=rMtP2*rMtP3*rMtP3;
				double rMtP9=rMtP8*rMtP;
				double rMtM=sig/(r-tM);
				double rMtM2=rMtM*rMtM;
				double rMtM3=rMtM2*rMtM;
				double rMtM8=rMtM2*rMtM3*rMtM3;
				double rMtM9=rMtM8*rMtM;
				double termU=3.1416*eps*delta*sig4/(12*t1*t2)*((rPtP8-rPtM8-rMtM8+rMtP8)/30-(rPtP2-rPtM2-rMtM2+rMtP2));
				double termF=3.1416*eps*delta*sig3/(3*t1*t2*r)*((rPtP9-rPtM9-rMtM9+rMtP9)/15-(rPtP3-rPtM3-rMtM3+rMtP3)/2);
				double termVN=termF/2;
				double termVT2=termU/2;
				uLJ[index_i]+=rhoJ*termU;
				uLJ[index_j]+=rhoI*termU;
				vNLJ[index_i]+=rhoJ*termVN*r2;
				vNLJ[index_j]+=rhoI*termVN*r2;
				vTLJ[index_j]+=rhoI*termVT2;
				vTLJ[index_i]+=rhoJ*termVT2;
				fLJ[index_i]+=rhoJ*termF*r;
				fLJ[index_j]+=-rhoI*termF*r;
			}
			else{
				uLJ[index_i]+=rhoJ*termURC;
				uLJ[index_j]+=rhoI*termURC;
				vNLJ[index_i]+=rhoJ*termVNRC*r2;
				vNLJ[index_j]+=rhoI*termVNRC*r2;
				vTLJ[index_j]+=rhoI*(termVTRC1*(rc2-r2)+termVTRC2);
				vTLJ[index_i]+=rhoJ*(termVTRC1*(rc2-r2)+termVTRC2);
				fLJ[index_i]+=rhoJ*termFRC*r;
				fLJ[index_j]+=-rhoI*termFRC*r;
			}
		}
	}

	// Calculation of the Forces on the slabs of the right hand side
	for (unsigned i=_slabs/2+_domainDecomposition->getRank(); i<_slabs; i+=_domainDecomposition->getNumProcs()) {
		const int index_i = i+si*_slabs+_slabs*numLJSum2[ci];

		rhoI=rho_l[index_i];

		const int index_i_sj_cj = i+sj*_slabs+_slabs*numLJSum2[cj];
		rhoJ=rho_l[index_i_sj_cj];

		vTLJ[index_i]+=rhoJ*(termVTRC1*rc2+termVTRC2);
		uLJ[index_i]+=rhoJ*termURC;
		for (unsigned j=i+1; j<_slabs; j++){
			const int index_j = j+sj*_slabs+_slabs*numLJSum2[cj];

			r=(j-i)*delta;
			r2=r*r;
			rhoJ=rho_l[index_j];
			if (j> i+cutoff_slabs){
				double rPtP=sig/(r+tP);
				double rPtP2=rPtP*rPtP;
				double rPtP3=rPtP2*rPtP;
				double rPtP8=rPtP2*rPtP3*rPtP3;
				double rPtP9=rPtP8*rPtP;
				double rPtM=sig/(r+tM);
				double rPtM2=rPtM*rPtM;
				double rPtM3=rPtM2*rPtM;
				double rPtM8=rPtM2*rPtM3*rPtM3;
				double rPtM9=rPtM8*rPtM;
				double rMtP=sig/(r-tP);
				double rMtP2=rMtP*rMtP;
				double rMtP3=rMtP2*rMtP;
				double rMtP8=rMtP2*rMtP3*rMtP3;
				double rMtP9=rMtP8*rMtP;
				double rMtM=sig/(r-tM);
				double rMtM2=rMtM*rMtM;
				double rMtM3=rMtM2*rMtM;
				double rMtM8=rMtM2*rMtM3*rMtM3;
				double rMtM9=rMtM8*rMtM;
				double termU=3.1416*eps*delta*sig4/(12*t1*t2)*((rPtP8-rPtM8-rMtM8+rMtP8)/30-(rPtP2-rPtM2-rMtM2+rMtP2));
				double termF=3.1416*eps*delta*sig3/(3*t1*t2*r)*((rPtP9-rPtM9-rMtM9+rMtP9)/15-(rPtP3-rPtM3-rMtM3+rMtP3)/2);
				double termVN=termF/2;
				double termVT2=termU/2;
				uLJ[index_i]+=rhoJ*termU;
				uLJ[index_j]+=rhoI*termU;
				vNLJ[index_i]+=rhoJ*termVN*r2;
				vNLJ[index_j]+=rhoI*termVN*r2;
				vTLJ[index_j]+=rhoI*termVT2;
				vTLJ[index_i]+=rhoJ*termVT2;
				fLJ[index_i]+=-rhoJ*termF*r;
				fLJ[index_j]+=rhoI*termF*r;
			}
			else{
				uLJ[index_i]+=rhoJ*termURC;
				uLJ[index_j]+=rhoI*termURC;
				vNLJ[index_i]+=rhoJ*termVNRC*r2;
				vNLJ[index_j]+=rhoI*termVNRC*r2;
				vTLJ[index_j]+=rhoI*(termVTRC1*(rc2-r2)+termVTRC2);
				vTLJ[index_i]+=rhoJ*(termVTRC1*(rc2-r2)+termVTRC2);
				fLJ[index_i]+=-rhoJ*termFRC*r;
				fLJ[index_j]+=rhoI*termFRC*r;
			}
		}
	}
}

void Planar::dipoleDipole(unsigned ci,unsigned cj,unsigned si,unsigned sj){
	double rc=cutoff;
	double rc2=rc*rc;
	double rc4=rc2*rc2;
	double rc6=rc4*rc2;
	double termU = 3.1416/4*muSquare[ci]*muSquare[cj]*delta / (3*temp);
	double termF = 3.1416 * muSquare[ci]*muSquare[cj]*delta / (3*temp);
	double termVN= 3.1416/2*muSquare[ci]*muSquare[cj]*delta / (3*temp);
	double termVT= termU;
	for (unsigned i=_domainDecomposition->getRank(); i<_slabs/2; i+=_domainDecomposition->getNumProcs()){
		const int index_i = i+si*_slabs+_slabs*numDipoleSum2[ci];

		double rhoI = rhoDipoleL[index_i];
		for (unsigned j=i+1; j<i+_slabs/2; j++){
			const int index_j = j+sj*_slabs+_slabs*numDipoleSum2[cj];

			double rhoJ = rhoDipoleL[index_j];
			double r=(j-i)*delta;
			double r2,r4,r6;
			if (j> i+cutoff_slabs){
			r2=r*r;
			r4=r2*r2;
			r6=r4*r2;
			}
			else{
			  r2=rc2;
			  r4=rc4;
			  r6=rc6;
			}
			fDipole[index_i] += termF*rhoJ/r6 * r;
			fDipole[index_j] -= termF*rhoI/r6 * r;
			uDipole[index_i] -= termU*rhoJ/r4;
			uDipole[index_j] -= termU*rhoI/r4;
			vNDipole[index_i]-= termVN*rhoJ/r6 *r*r;
			vNDipole[index_j]-= termVN*rhoI/r6 *r*r;
			vTDipole[index_i]-= termVT*rhoJ/r6 *(1.5*r2 - r*r);
			vTDipole[index_j]-= termVT*rhoI/r6 *(1.5*r2 - r*r);

		}
		// Calculation of the Periodic boundary
		for (unsigned j=_slabs/2+i; j<_slabs; j++){
			const int index_j = j+sj*_slabs+_slabs*numDipoleSum2[cj];

			double rhoJ = rhoDipoleL[index_j];
			double r=(_slabs-j+i)*delta;
			double r2,r4,r6;
			if (j <_slabs-cutoff_slabs+i){
			r2=r*r;
			r4=r2*r2;
			r6=r4*r2;
			}
			else{
			  r2=rc2;
			  r4=rc4;
			  r6=rc6;
			}
			fDipole[index_i] -= termF*rhoJ/r6 * r;
			fDipole[index_j] += termF*rhoI/r6 * r;
			uDipole[index_i] -= termU*rhoJ/r4;
			uDipole[index_j] -= termU*rhoI/r4;
			vNDipole[index_i]-= termVN*rhoJ/r6 *r*r;
			vNDipole[index_j]-= termVN*rhoI/r6 *r*r;
			vTDipole[index_i]-= termVT*rhoJ/r6 *(1.5*r2 - r*r);
			vTDipole[index_j]-= termVT*rhoI/r6 *(1.5*r2 - r*r);
		}
	}

	// Calculation of the Forces on the slabs of the right hand side
	for (unsigned i=_slabs/2+_domainDecomposition->getRank(); i<_slabs; i+=_domainDecomposition->getNumProcs()) {
		const int index_i = i+si*_slabs+_slabs*numDipoleSum2[ci];

		double rhoI = rhoDipoleL[index_i];
		for (unsigned j=i+1; j<_slabs; j++){
			const int index_j = j+sj*_slabs+_slabs*numDipoleSum2[cj];

			double rhoJ = rhoDipoleL[index_j];
			double r=(j-i)*delta;
			double r2,r4,r6;
			if (j> i+cutoff_slabs){
			r2=r*r;
			r4=r2*r2;
			r6=r4*r2;
			}
			else{
			  r2=rc2;
			  r4=rc4;
			  r6=rc6;
			}
			fDipole[index_i] += termF*rhoJ/r6 * r;
			fDipole[index_j] -= termF*rhoI/r6 * r;
			uDipole[index_i] -= termU*rhoJ/r4;
			uDipole[index_j] -= termU*rhoI/r4;
			vNDipole[index_i]-= termVN*rhoJ/r6 *r*r;
			vNDipole[index_j]-= termVN*rhoI/r6 *r*r;
			vTDipole[index_i]-= termVT*rhoJ/r6 *(1.5*r2 - r*r);
			vTDipole[index_j]-= termVT*rhoI/r6 *(1.5*r2 - r*r);
		}
	}
}

double Planar::lrcLJ(Molecule* mol){
	double potentialEnergy = 0.;
	unsigned cid=mol->componentid();
	for (unsigned i=0; i<numLJ[cid]; i++){
		int loc=(mol->r(1)+mol->ljcenter_d(i)[1])/delta;
		if (loc < 0){
			loc=loc+_slabs;
		}
		else if (loc > sint-1){
			loc=loc-_slabs;
		}
		potentialEnergy += uLJ[loc+i*_slabs+_slabs*numLJSum2[cid]];
	}
	for (unsigned i=0;i<numDipole[cid]; i++){
		int loc=(mol->r(1)+mol->dipole_d(i)[1])/delta;
		if (loc < 0){
			loc=loc+_slabs;
		}
		else if (loc > sint-1){
			loc=loc-_slabs;
		}
		potentialEnergy += uDipole[loc+i*_slabs+_slabs*numDipoleSum2[cid]];
	}
	return potentialEnergy;
}

void Planar::directDensityProfile(){
	_smooth = false;
}

void Planar::writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstepI)
{
	if( 0 != simstepI % _nWriteFreqProfiles || simstepI < _nStartWritingProfiles || simstepI > _nStopWritingProfiles)
		return;

	// writing .dat-files
	std::stringstream filenamestream;
	filenamestream << "LRC" << "_TS" << fill_width('0', 9) << simstepI << ".dat";

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	// int numprocs = domainDecomp->getNumProcs();
	if (rank!= 0)
		return;
#endif

	std::stringstream outputstream;
	double lengthPerSlab = _domain->getGlobalLength(1)/_slabs;

	// header
	outputstream << "                     pos";
	for(uint32_t si=0; si<numLJSum; ++si)
	{
		outputstream << "            rho_l_LRC[" << si << "]";
		outputstream << "                F_LRC[" << si << "]";
		outputstream << "                u_LRC[" << si << "]";
	}
	outputstream << std::endl;

	// data
	for(uint32_t pi=0; pi<_slabs; ++pi)
	{
		double pos = lengthPerSlab*(pi+0.5);
		outputstream << std::setw(24) << pos;
		for(uint32_t si=0; si<numLJSum; ++si)
		{
			outputstream << std::setw(24) << FORMAT_SCI_MAX_DIGITS << rho_l[_slabs*si+pi];
			if ((pos >= _region.actPos[0]) && (pos <= _region.actPos[1])) {
				outputstream << std::setw(24) << FORMAT_SCI_MAX_DIGITS << fLJ[_slabs*si+pi];
			} else {
				outputstream << std::setw(24) << FORMAT_SCI_MAX_DIGITS << 0.0;
			}
			outputstream << std::setw(24) << FORMAT_SCI_MAX_DIGITS << uLJ[_slabs*si+pi];
		}
		outputstream << std::endl;
	} // loop: pos

	// open file for writing
	// scalar
	std::ofstream fileout(filenamestream.str().c_str(), std::ios::out);
	fileout << outputstream.str();
	fileout.close();
}

SubjectBase* Planar::getSubject() {
	SubjectBase* subject = nullptr;
	std::list<PluginBase*>& plugins = *(global_simulation->getPluginList() );
	for (auto&& pit:plugins) {
		std::string name = pit->getPluginName();
		if(name == "DistControl") {
			subject = dynamic_cast<SubjectBase*>(pit);
		}
	}
	return subject;
}

void Planar::update(SubjectBase* subject) {
	auto* distControl = dynamic_cast<DistControl*>(subject);
	double dMidpointLeft, dMidpointRight, origin;
	dMidpointLeft = dMidpointRight = 0.;
	if(nullptr != distControl) {
		dMidpointLeft = distControl->GetInterfaceMidLeft();
		dMidpointRight = distControl->GetInterfaceMidRight();
#ifndef NDEBUG
		std::cout << "Long Range Correction: Received value for left interface position:  " << dMidpointLeft << std::endl;
		std::cout << "Long Range Correction: Received value for right interface position: " << dMidpointRight << std::endl;
#endif
	}

	for (unsigned short i=0; i<2; i++) {

		switch(_region.refPosID[i]) {
		case 0:
			origin = 0.;
			break;
		case 1:
			origin = dMidpointLeft;
			break;
		case 2:
			origin = dMidpointRight;
			break;
		default:
			origin = 0.;
		}
		_region.actPos[i] = origin + _region.refPos[i];
	}
}

/*
 * SphericalSampling.cpp
 *
 *  Created on: Aug 2024
 *      Author: JakNiem
 */


#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "longRange/Spherical.h"
#include <cmath>



using Log::global_log;


Spherical::Spherical(double /*cutoffT*/, double cutoffLJ, Domain* domain, DomainDecompBase* domainDecomposition,
					 ParticleContainer* particleContainer, Simulation* /*simulation*/):
        _cutoffLJ{cutoffLJ},
        _domain{domain}, 
        _domainDecomposition{domainDecomposition},
        _particleContainer{particleContainer}
    {
#ifdef DEBUG
    global_log->info() << "[SphericalLRC]: Constructor called."<< std::endl;
#endif
	global_log->info() << "Long Range Correction for spherical interfaces is used. " << std::endl;
	global_log->info() << "[CONSTRUCTOR STILL UNDER CONSTRCUTION]. " << std::endl;
    }

Spherical::~Spherical(){}
#ifdef DEBUG
    global_log->info() << "[SphericalLRC]: Spherical::~Spherical(){} called."<< std::endl;
#endif



void Spherical::readXML(XMLfileUnits& xmlconfig)
{
#ifdef DEBUG
    global_log->info() << "[SphericalLRC]: Spherical::readXML(...) called "<< std::endl;
#endif
	global_log->info() << "[Long Range Correction] reading XML for paramters of SphericalLRC." << std::endl;
	xmlconfig.getNodeValue("shells", _nShells);
	xmlconfig.getNodeValue("temperature", _T);
	xmlconfig.getNodeValue("bubble", _isBubble);   // default: false
	xmlconfig.getNodeValue("calculationFreq", _calcFreq);
	xmlconfig.getNodeValue("writeFreq", _writeFreq);
	xmlconfig.getNodeValue("disableLRC", _disableLRC);
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);

	_lenVector = _nShells+1;
	
    // making sure calcFreq and writeFreq are valid:
    if(_calcFreq < 1){
        _calcFreq = 100;
        global_log->warning() << "[Long Range Correction] Warning: _calcFreq invalid. setting _calcFreq = " << _calcFreq << std::endl; 
    }
    if(_writeFreq < 1){
        _writeFreq = _calcFreq * 10;
        global_log->warning() << "[Long Range Correction] Warning: _writeFreq invalid. setting _writeFreq = " << _writeFreq << std::endl; 
    }
	// making sure that writeFreq is a multiple of calcFreq:
	_writeFreq = (_writeFreq / _calcFreq) * _calcFreq;
	if (_writeFreq == 0) {
		_writeFreq = _calcFreq;
	}

	global_log->info() << "[Long Range Correction] Using " << _nShells << " shells for profiles to calculate LRC." << std::endl;
	if (_isBubble) {
		global_log->info() << "[Long Range Correction] System contains a bubble. COMalignerBubble plugin is required." << std::endl;
	} else {
		global_log->info() << "[Long Range Correction] System contains a droplet. COMaligner plugin is required."<< std::endl;
	}

/* 
	// init --> nötig? wird doch von MarDyn als nächstes gemacht! (vllt mal probieren, ob / wie oft gecalled wird)
	this->init();
 */
}



void Spherical::init()
{
	global_log->info() << "[Long Range Correction] Initializing. Is this function called once, twice or not at all?!" << std::endl;


    resizeVectors();
    resetVectors();


	_globalNumMols = _domain->getglobalNumMolecules(true, _particleContainer, _domainDecomposition);

	_globalBoxLength[0] = _domain->getGlobalLength(0);
	_globalBoxLength[1] = _domain->getGlobalLength(1);
	_globalBoxLength[2] = _domain->getGlobalLength(2);

	_globalCenter[0] = _globalBoxLength[0]/2.;
	_globalCenter[1] = _globalBoxLength[1]/2.;
	_globalCenter[2] = _globalBoxLength[2]/2.;

    _distMax = *(std::min_element(_globalBoxLength.begin(), _globalBoxLength.end())) * 0.5;
    // _shellWidth = _distMax/_nShells;
}


void Spherical::calculateLongRange()
{
#ifdef DEBUG
    global_log->info() << "[SphericalLRC]: Spherical::calculateLongRange() called "<< std::endl;
#endif

	int rank = _domainDecomposition->getRank();
	uint64_t simstep = _simulation.getSimulationStep();


//////////////////////////////////////////////////////////////
//					density calculations					//
//////////////////////////////////////////////////////////////


	//instantanes Dichteprofil:
    CommVar<std::vector<double>> rho_step;
    CommVar<std::vector<unsigned long>> numMolecules_step;
	resizeCommVarVector(rho_step); //muss der auch =0 gesetzt werden?
	resizeCommVarVector(numMolecules_step); //muss der auch =0 gesetzt werden?

	for (auto mol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); mol.isValid();
		 ++mol) {
		unsigned long molID = mol->getID();
		
		//distance to center: ksi
		double distCenter_x = mol->r(0) - _globalCenter[0];
		double distCenter_y = mol->r(1) - _globalCenter[1];
		double distCenter_z = mol->r(2) - _globalCenter[2];
		double distCenter2 = distCenter_x*distCenter_x + distCenter_y*distCenter_y + distCenter_z*distCenter_z;
		double distCenter = std::sqrt(distCenter2);
		
		//shell index of molecule:
		const unsigned int k = std::min(_nShells, static_cast<unsigned int>((distCenter/_distMax)*_nShells));  // Shell index of molecule 
		// shellidOfMolecule[mol] = k;


		numMolecules_step.local[k] ++;

	}
	_numsamples++;

	// Gather quantities:
#ifdef ENABLE_MPI
    MPI_Allreduce(numMolecules_step.local.data(), numMolecules_step.global.data(), _lenVector, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    for (unsigned long i = 0; i < _lenVector; i++) {
        numMolecules_step.global[i] = numMolecules_step.local[i];
    }
	//TODO: würde das hier nicht auch ohne for-loop gehen? (einfach numMolecules_step.global = numMolecules_step.local;) 
#endif

	// adding numMolecules_step.global to _numMolecules_accum
	std::transform(_numMolecules_accum.begin(), _numMolecules_accum.end(), numMolecules_step.global.begin(),
               _numMolecules_accum.begin(), std::plus<unsigned long>());
	
	// calculation of density from _numMolecules_accum. TODO: rewrite using stl!
	_shellLowerBound[0] = 0.;
	for(int i = 0; i < _lenVector-1; i++){
		_shellLowerBound[i+1] = _distMax*(i+1.)/(_nShells);
		// global_log->info() << "_shellLowerBound["<<i+1<<"]="<<_shellLowerBound[i+1] <<std::endl;
		_shellVolume[i] = 4./3. * M_PI * ( std::pow(_shellLowerBound[i+1],3) - std::pow(_shellLowerBound[i],3) );
		// global_log->info() << "_shellVolume["<<i<<"]="<<_shellVolume[i] <<std::endl;
		_density_avg[i] = _numMolecules_accum[i] / (_numsamples * _shellVolume[i]);
		// global_log->info() << "_density_avg["<<i<<"]="<<_density_avg[i] <<std::endl;
	}
	_shellVolume[_lenVector-1] = _globalBoxLength[0]*_globalBoxLength[1]*_globalBoxLength[2] - 4./3.*M_PI*std::pow(_shellLowerBound[_lenVector-1],3);
	_density_avg[_lenVector-1] = _numMolecules_accum[_lenVector-1] / (_numsamples * _shellVolume[_lenVector-1]);


	if ((simstep) % _calcFreq == 0) {  // 1000
		calculateBulkDensities();
		calculateTanhProfile(); // NEEDS DEBUGGING
	}
	
	//////////////////////////////////////////////////////////////////
	//					correction of U, F, ViN, ViT				//
	//////////////////////////////////////////////////////////////////


	//4 lines to initialize one vector? ugly ...
	CommVar<std::vector<double>> uCorrection_step; //by Shell
	resizeCommVarVector(uCorrection_step);
	std::fill(uCorrection_step.global.begin(), uCorrection_step.global.end(), 0.); 
	std::fill(uCorrection_step.local.begin(), uCorrection_step.local.end(), 0.); 

	//TODO: move to object? (only calc once):
	double rCutoff = 5.;	//TODO: find correct rc
	double rm = 16.; 		//TODO: find value for rm (hardcoded is quatsch)

	double epsilon = 1.; 	//TODO: find correct epsilon
	double sigma = 1.;		//TODO: find correct sigma
	double sigma6 = std::pow(epsilon, 6); 


	
	double generalPrefactor =  2.*M_PI*epsilon*sigma6;

	for (auto mol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); mol.isValid(); ++mol) {
		//MOLECULE PROPERTIES:
		unsigned long molID = mol->getID();
		//distance to center ("ksi" in Nitze2021):
		double distCenter_x = mol->r(0) - _globalCenter[0];
		double distCenter_y = mol->r(1) - _globalCenter[1];
		double distCenter_z = mol->r(2) - _globalCenter[2];
		double distCenter2 = distCenter_x*distCenter_x + distCenter_y*distCenter_y + distCenter_z*distCenter_z;
		double distCenter = std::sqrt(distCenter2);
		double distCenterInv = 1./(distCenter);
		double distCenterInv3 = std::pow(distCenterInv,3);

		//MOLECULE CORRECTION TERMS:
		double fCorrectionMol = 0.;
		double virNCorrectionMol = 0.;
		double virTCorrectionMol = 0.;

		double prefactorU = generalPrefactor;
		double prefactorF = 2.*generalPrefactor * distCenterInv3;
		double prefactorVirN = generalPrefactor * distCenterInv3; 
		double prefactorVirT = 4.*generalPrefactor * distCenterInv;

		//LOOP OVER SHELLS:
		for(int k = 0; k < _nShells+1; k++){ //using k as shell index for consitency with Nitzke2021, equation (14)
			//SHELL PROPERTIES:
			double Rk = _shellLowerBound[k]; //TODO: this is not the definition of Rk that Nitzke2021 used in their implementation!
			double Rk2 = Rk*Rk;
			double shellWidth = _distMax/_nShells; //TODO: calculation is performed for each molecule over again. perhaps rather save in Spherical::-object

			//MOLECULE<->SHELL INTERACTION RADIUS BOUNDS:
			double rLowerBound = rCutoff < std::abs(distCenter-Rk) ? distCenter-Rk : rCutoff;
			double rLowerBound2 = rLowerBound*rLowerBound;
			double rLowerBoundInv2 = std::pow(rLowerBound, -2);
			double rLowerBoundInv4 =  rLowerBoundInv2 * rLowerBoundInv2;
			double rLowerBoundInv6 =  rLowerBoundInv2 * rLowerBoundInv4;
			double rLowerBoundInv8 =  rLowerBoundInv4 * rLowerBoundInv4;
			double rLowerBoundInv10 = rLowerBoundInv4 * rLowerBoundInv6;
			double rLowerBoundInv12 = rLowerBoundInv6 * rLowerBoundInv6;

			double rUpperBound = std::min(distCenter+Rk, rm);
			double rUpperBound2 = rUpperBound*rUpperBound;
			double rUpperBoundInv2 = std::pow(rUpperBound, -2);
			double rUpperBoundInv4 = rUpperBoundInv2 * rUpperBoundInv2;
			double rUpperBoundInv6 = rUpperBoundInv2 * rUpperBoundInv4;
			double rUpperBoundInv8 = rUpperBoundInv4 * rUpperBoundInv4;
			double rUpperBoundInv10 = rUpperBoundInv4 * rUpperBoundInv6;
			double rUpperBoundInv12 = rUpperBoundInv6 * rUpperBoundInv6;

			//U CORRECTION (CONTRIBUTION OF CURRENT MOL TO TO SHELL K):
			uCorrection_step.local[k] += prefactorU * _virtual_density[k]*Rk*shellWidth * (
										.4 * sigma6 * (std::pow((distCenter + Rk), -10) - std::pow(rLowerBound, -10) ) 
										- (std::pow((distCenter + Rk), -4)  - std::pow(rLowerBound, -4))
										) / distCenter;
			//F CORRECTION (CONTRIBUTION OF SHELL K TO CURRENT MOL):
			fCorrectionMol += _virtual_density[k]*Rk*shellWidth *(
								sigma6 * ((1.2*rUpperBound2 + distCenter2 - Rk2)*rUpperBoundInv12 - (1.2*rLowerBound2 + distCenter2 - Rk2)*rLowerBoundInv12)	
								- 		 ((1.5*rUpperBound2 + distCenter2 - Rk2)*rUpperBoundInv6  - (1.5*rLowerBound2 + distCenter2 - Rk2)*rLowerBoundInv6)
								);

			//VIRIAL CORRECTION (CONTRIBUTION OF SHELL K TO CURRENT MOL):
			virNCorrectionMol += _virtual_density[k]*Rk*shellWidth *(
								  (1.5 * sigma6 * (rUpperBoundInv8 - rLowerBoundInv8) - 3. * (rUpperBoundInv2 - rLowerBoundInv2))
								+ 2.*(distCenter2- Rk2) * (1.2 * sigma6*(rUpperBoundInv10 - rLowerBoundInv10) - 1.5*(rUpperBoundInv4 - rLowerBoundInv4))
								+ std::pow((distCenter2- Rk2), 2) * (sigma6*(rUpperBoundInv12 - rLowerBoundInv12) - (rUpperBoundInv6 - rLowerBoundInv6))
								);
			virTCorrectionMol += _virtual_density[k]*Rk*shellWidth *(
								 1.2 * sigma6*(rUpperBoundInv10-rLowerBoundInv10) 
								-1.5 * (rUpperBoundInv4 - rLowerBoundInv4)
								);

		}
		fCorrectionMol *= prefactorF;
		double fCorrectionMolKartesian[3] = {distCenter_x * fCorrectionMol,
											 distCenter_y * fCorrectionMol,
											 distCenter_z * fCorrectionMol};
		virNCorrectionMol *= prefactorVirN;
		virTCorrectionMol = virTCorrectionMol*prefactorVirT - virNCorrectionMol;

		mol->Fadd(fCorrectionMolKartesian);
		mol->ViNadd(virNCorrectionMol);
		mol->ViTadd(virTCorrectionMol);
	}

	// Gather quantities:
#ifdef ENABLE_MPI
	MPI_Allreduce(uCorrection_step.local.data(), uCorrection_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
	for (unsigned long i = 0; i < _lenVector; i++) {
		uCorrection_step.global[i] = uCorrection_step.local[i];
	}
	//TODO: würde das hier nicht auch ohne for-loop gehen? (einfach numMolecules_step.global = numMolecules_step.local;) 
#endif


#ifdef DEBUG
    global_log->info() << "[SphericalLRC]: uCorrection_step for simstep "<< simstep << ":\n";
	for(int k = 0; k<_nShells+1; k++){
    	global_log->info() << "shell "<< k << " : " << uCorrection_step.global[k]<<"\n";
	}
	global_log->info() << std::flush;
#endif



}

void Spherical::writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep)
{
#ifdef DEBUG
    global_log->info() << "[SphericalLRC]: Spherical::writeProfiles(...) called "<< std::endl;
#endif

}





void Spherical::calculateBulkDensities() //the procedure used here is not good. I'll keep it for now, but 
{
	// NOTE: procedure to calculate the instantaneous inside and outside densities is risky here, because it makes
	// a priori asumptions about where to find the phase boundary! Better Ideas are welcome.
	_bulkBoundaries.inside_from = _nShells / 10; //10; 
	_bulkBoundaries.inside_to = _nShells / 5 + 5; // 35;
	_bulkBoundaries.outside_from = _nShells - 2 - _nShells / 10 - 5;  // 72;
	_bulkBoundaries.outside_to = _nShells - 1;
	
	// systemcenter is at half of domain size (densities etc are averaged spherically around systemcenter)
	double inside_from_coord = (_bulkBoundaries.inside_from*_globalCenter[0])/_nShells;
	double inside_to_coord = (_bulkBoundaries.inside_to*_globalCenter[0])/_nShells;
	double outside_from_coord = (_bulkBoundaries.outside_from*_globalCenter[0])/_nShells;
	double outside_to_coord = (_bulkBoundaries.outside_to*_globalCenter[0])/_nShells;

	global_log->warning() << "[Long Range Correction] Calculating rho_inside from " << inside_from_coord << " (Shell " << _bulkBoundaries.inside_from << ")"
																	<< " to " << inside_to_coord << " (Shell " << _bulkBoundaries.inside_to << ")"
														<< " and rho_out from " << outside_from_coord << " (Shell " << _bulkBoundaries.outside_from << ")"
																	<< " to " << outside_to_coord << " (Shell " << _bulkBoundaries.outside_to << ")" << std::endl;


	_rho_in = 0.;
	for (unsigned int i = _bulkBoundaries.inside_from; i < _bulkBoundaries.inside_to; i++) {
		_rho_in += _numMolecules_accum[i]/(_numsamples*_shellVolume[i]);
	}
	_rho_in /= (_bulkBoundaries.inside_to - _bulkBoundaries.inside_from);

	_rho_out = 0.;
	for (unsigned int i = _bulkBoundaries.outside_from; i < _bulkBoundaries.outside_to; i++) {
		_rho_out += _numMolecules_accum[i]/(_numsamples*_shellVolume[i]);
	}
	_rho_out /= (_bulkBoundaries.outside_to - _bulkBoundaries.outside_from);

	global_log->info() << "[Long Range Correction] Averaged rho: rho_in = " << _rho_in << " rho_out = " << _rho_out << std::endl;
}

void Spherical::calculateTanhProfile()
{
	
	//TODO: rework
	// info: uses (all time) average density

	// info: uses _shellLowerBound as radius. Implementation by Nitzke used upper(?) buond shell radius!

	// D0 with 1090 (Baidakov et al.)
	double Dmin = 0.0;
	double Dmax = 0.0;

	if (!_isBubble) { //droplet
		double r10 = _rho_out + 0.1 * (_rho_in - _rho_out);
		double r90 = _rho_out + 0.9 * (_rho_in - _rho_out);
		for (unsigned int i = 1; i < (_nShells - 10); i++) {  // TODO:WHY IS THERE A HARDCODED 10 HERE?!
			if (_density_avg[i] > r90) {
				Dmin = _shellLowerBound[i];
			}
		}
		for (unsigned int i = 1; i < (_nShells - 10); i++) {  // some value/limitation of the search space could be needed for the largest shells in case of droplet (low density)
			// This could be written much better with a "backwards" loop and/or a break condition
			unsigned int j = _nShells - i;
			if (_density_avg[j] < r10) {
				Dmax = _shellLowerBound[j];
			}
		}
	} else { //bubble
		double r10 = _rho_in + 0.1 * (_rho_out - _rho_in);
		double r90 = _rho_in + 0.9 * (_rho_out - _rho_in);
		for (unsigned int i = 1; i < (_nShells - 10); i++) {  // some value/limitation of the search space could be needed for the smallest shells in case of bubble (low density)
			if (_density_avg[i] < r10) {
				Dmin = _shellLowerBound[i];
			}
		}
		for (unsigned int i = 1; i < (_nShells - 10); i++) {
			unsigned int j = _nShells - i;
			if (_density_avg[j] > r90) {
				Dmax = _shellLowerBound[j];
			}
		}
	}

	double D0 = (Dmax - Dmin);
	double R0 = Dmin + 0.5 * D0;

	// density profile (TODO:replace forloop with stl-hack)
	for (unsigned int i = 0; i < _nShells; i++) {

		double a =  0.5 * (_rho_in + _rho_out);
		double b = 0.5 * (_rho_in - _rho_out);
		double c = 2 * (_shellLowerBound[i+1] - R0) / D0;
		_density_avg_fitted[i] =a - b * std::tanh(c);
		_virtual_density[i] = _density_avg_fitted[i] - _rho_out; 
	}
	_density_avg_fitted[_nShells] =_density_avg_fitted[_nShells-1];
	_virtual_density[_nShells] = _density_avg_fitted[_nShells] - _rho_out; 

	//TODO: vielleicht müssen kleinere werte hier gleich 0 gesetzt werden:
	/* 	// shift for correction, since outer bulk only needs homogeneous correction
		if (droplet) {
			for (unsigned int i = 0; i < NShells; i++) {
				if (rhoShellsT[i] >= 1.02 * rho_out) {
					rhoShellsT[i] -= rho_out;
				} else {
					rhoShellsT[i] = 0.0;
				}
			}
		} else {
			for (unsigned int i = 0; i < NShells; i++) {
				rhoShellsT[i] -= rho_out;
				// if (rhoShellsT[i] <= 0.998 * rho_out) {
				// 	rhoShellsT[i] -= rho_out;
				// } else {
				// 	rhoShellsT[i] = 0.0;
				// }
			}
		} */


}



//////////////////////////////////////////////////////////////


// Resize vectors
void Spherical::resizeVectors() {
    _shellVolume.resize(_lenVector);
    _shellLowerBound.resize(_lenVector);
    _numMolecules_accum.resize(_lenVector);
    _density_avg.resize(_lenVector);
    _density_avg_fitted.resize(_lenVector);
    _virtual_density.resize(_lenVector);
	//TODO: alle vektoren hier vorhanden?

}

// Fill vectors with zeros
void Spherical::resetVectors() {
    std::fill(_shellVolume.begin(), _shellVolume.end(), 0.); // korrekte 0 einfüllen: 0.0f, 0ul, ...
    std::fill(_shellLowerBound.begin(), _shellLowerBound.end(), 0.); // korrekte 0 einfüllen: 0.0f, 0ul, ...
    std::fill(_numMolecules_accum.begin(), _numMolecules_accum.end(), 0ul); // korrekte 0 einfüllen: 0.0f, 0ul, ...
    std::fill(_density_avg.begin(), _density_avg.end(), 0.); // korrekte 0 einfüllen: 0.0f, 0ul, ...
    std::fill(_density_avg_fitted.begin(), _density_avg_fitted.end(), 0.); // korrekte 0 einfüllen: 0.0f, 0ul, ...
    std::fill(_virtual_density.begin(), _virtual_density.end(), 0.); // korrekte 0 einfüllen: 0.0f, 0ul, ...
	//TODO: alle vektoren hier vorhanden?
}
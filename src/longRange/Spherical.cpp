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
#include "utils/FileUtils.h"



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

    initVectors();

	_globalNumMols = _domain->getglobalNumMolecules(true, _particleContainer, _domainDecomposition);

	_globalBoxLength[0] = _domain->getGlobalLength(0);
	_globalBoxLength[1] = _domain->getGlobalLength(1);
	_globalBoxLength[2] = _domain->getGlobalLength(2);

	_globalCenter[0] = _globalBoxLength[0]/2.;
	_globalCenter[1] = _globalBoxLength[1]/2.;
	_globalCenter[2] = _globalBoxLength[2]/2.;

    _distMax = *(std::min_element(_globalBoxLength.begin(), _globalBoxLength.end())) * 0.5;
    // _shellWidth = _distMax/_nShells;


	// TODO: rewrite using stl!
	_shellLowerBound[0] = 0.;
	for(int i = 0; i < _lenVector-1; i++){
		_shellLowerBound[i+1] = _distMax*(i+1.)/(_nShells); //TODO: move calculatoins of shellLowrBound and shellVolume to init()
		// global_log->info() << "_shellLowerBound["<<i+1<<"]="<<_shellLowerBound[i+1] <<std::endl;
		_shellVolume[i] = 4./3. * M_PI * ( std::pow(_shellLowerBound[i+1],3) - std::pow(_shellLowerBound[i],3) );
		// global_log->info() << "_shellVolume["<<i<<"]="<<_shellVolume[i] <<std::endl;
	}
	_shellVolume[_lenVector-1] = _globalBoxLength[0]*_globalBoxLength[1]*_globalBoxLength[2] - 4./3.*M_PI*std::pow(_shellLowerBound[_lenVector-1],3);
}


void Spherical::calculateLongRange()
{
#ifdef DEBUG
    global_log->info() << "[SphericalLRC]: Spherical::calculateLongRange() called "<< std::endl;
#endif

	int rank = _domainDecomposition->getRank();
	uint64_t simstep = _simulation.getSimulationStep();


//////////////////////////////////////////////////////////////
//					density calculations					// //TODO: JEDEN ZEITSCHRITT, ODER NUR IN CALCFREQ??
//////////////////////////////////////////////////////////////
	

	//instantanes Dichteprofil:
    CommVar<std::vector<unsigned long>> numMolecules_step;
	initCommVarVector(numMolecules_step); 

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
/* #ifdef ENABLE_MPI //TODO: MPI
    MPI_Allreduce(numMolecules_step.local.data(), numMolecules_step.global.data(), _lenVector, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
 */    for (unsigned long i = 0; i < _lenVector; i++) {
        numMolecules_step.global[i] = numMolecules_step.local[i];
    }
	//TODO: würde das hier nicht auch ohne for-loop gehen? (einfach numMolecules_step.global = numMolecules_step.local;) 
/* endif */

	// adding numMolecules_step.global to _numMolecules_accum
	std::transform(_numMolecules_accum.begin(), _numMolecules_accum.end(), numMolecules_step.global.begin(),
               _numMolecules_accum.begin(), std::plus<unsigned long>());
				


	//////////////////////////////////////////////////////////////////
	//			CORRECTION	CALCULATION (only in calcfreq)			//
	//////////////////////////////////////////////////////////////////
	if ((simstep) % _calcFreq == 0) {  // wird dann auch bei simstep==0 ausgeführt. ist das weise?
			// calculation_2_getShellPropertiesByLoopingOverAllMolecules();


		// CALCULATION OF DENSITY FROM _NUMMOLECULES_ACCUM. TODO: rewrite using stl! OR: use openmp
		for(int i = 0; i < _lenVector-1; i++){
			_density_avg[i] = _numMolecules_accum[i] / (_numsamples * _shellVolume[i]);
			// global_log->info() << "_density_avg["<<i<<"]="<<_density_avg[i] <<std::endl;
		}
		_density_avg[_lenVector-1] = _numMolecules_accum[_lenVector-1] / (_numsamples * _shellVolume[_lenVector-1]);

		calculateBulkDensities();
		calculateVirtualDensity(); 



		// calculationV1_fromPaper_withIsabelsPrefactorAndRlow_repParticles();
		calculationV3_isabelsMethod();
	}   



	//////////////////////////////////////////////////////////////////
	//				APPLY CORRECTION (each step)					//
	//////////////////////////////////////////////////////////////////
	if(simstep  > _calcFreq) //dont apply before calculated
	{
		for (auto mol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); mol.isValid(); ++mol) {
			//MOLECULE PROPERTIES:
			unsigned long molID = mol->getID();
			double distCenter_x = mol->r(0) - _globalCenter[0];
			double distCenter_y = mol->r(1) - _globalCenter[1];
			double distCenter_z = mol->r(2) - _globalCenter[2];
			double distCenter2 = distCenter_x*distCenter_x + distCenter_y*distCenter_y + distCenter_z*distCenter_z;
			double distCenter = std::sqrt(distCenter2);
			double distCenterInv = 1./distCenter;
			const unsigned int shellIDofMol = std::min(_nShells, static_cast<unsigned int>((distCenter/_distMax)*_nShells));  // Shell index of molecule 

			//APPLY CORRECTIONS:
			double FcorrMolApprox[3] = {.0,.0,.0};
			FcorrMolApprox[0] = _FcorrectionShell[shellIDofMol] * distCenter_x*distCenterInv; 
			FcorrMolApprox[1] = _FcorrectionShell[shellIDofMol] * distCenter_y*distCenterInv;
			FcorrMolApprox[2] = _FcorrectionShell[shellIDofMol] * distCenter_z*distCenterInv;
			

			// global_log->info() << "molid : shellid = "<< molID <<" :  "<< shellIDofMol << ":\n";
			// 	global_log->info() << "FcorrMolApprox["<<0<<"] = "<< _FcorrectionShell[shellIDofMol] << " * "<< distCenter_x << ":\n";
			// 	global_log->info() << "FcorrMolApprox["<<1<<"] = "<< _FcorrectionShell[shellIDofMol] << " * "<< distCenter_y << ":\n";
			// 	global_log->info() << "FcorrMolApprox["<<2<<"] = "<< _FcorrectionShell[shellIDofMol] << " * "<< distCenter_z << ":\n";


			mardyn_assert(!std::isnan(FcorrMolApprox[0])); // catches NaN
			mardyn_assert(!std::isnan(FcorrMolApprox[1])); // catches NaN
			mardyn_assert(!std::isnan(FcorrMolApprox[2])); // catches NaN



			// mol->Uadd(_UcorrectionShell[shellIDofMol]);
			mol->Fadd(FcorrMolApprox);
			mol->ViNadd(_VirNcorrectionShell[shellIDofMol]);
			mol->ViTadd(_VirTcorrectionShell[shellIDofMol]);
		}

		double UCorrSum_global = 0.;	//TODO (sum of all virial corrections??)
		double VirCorrSum_global = 0.;	//TODO (sum of all energy corrections??)
		_domain->setUpotCorr(UCorrSum_global);   //muss das in jedem step passieren? wenn die terme konstant sind, dann könnte ichs einfach nur bei "calculateLongrange" machen.	
		_domain->setVirialCorr(VirCorrSum_global);	
	}


/* #ifdef DEBUG
    global_log->info() << "[SphericalLRC]: UcorrectionSum_step for simstep "<< simstep << ":\n";
	for(int k = 0; k<_nShells+1; k++){
    	global_log->info() << "shell "<< k << " : " << UcorrectionSum_step.global[k]<<"\n";
	}
	global_log->info() << std::flush;
#endif */
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

void Spherical::calculateVirtualDensity()
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
			for (unsigned int i = 0; i < _nShells; i++) {
				if (rhoShellsT[i] >= 1.02 * rho_out) {
					rhoShellsT[i] -= rho_out;
				} else {
					rhoShellsT[i] = 0.0;
				}
			}
		} else {
			for (unsigned int i = 0; i < _nShells; i++) {
				rhoShellsT[i] -= rho_out;
				// if (rhoShellsT[i] <= 0.998 * rho_out) {
				// 	rhoShellsT[i] -= rho_out;
				// } else {
				// 	rhoShellsT[i] = 0.0;
				// }
			}
		} */


}








void Spherical::writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep)
{
#ifdef DEBUG
    global_log->info() << "[SphericalLRC]: Spherical::writeProfiles(...) called. This function currently doesn't do anything."<< std::endl;
#endif
	if(simstep == 0){return;}

	//SAMPLE MOLECULE DATA ON EACH TIMESTEP
	CommVar<std::vector<unsigned long>> numMolecules_step;
	CommVar<std::vector<double>> virN_step;
	CommVar<std::vector<double>> virT_step;
	CommVar<std::vector<double>> virX_step;
	CommVar<std::vector<double>> virY_step;
	CommVar<std::vector<double>> virZ_step;
	CommVar<std::vector<double>> velocityN_step;

	initCommVarVector(numMolecules_step);
	initCommVarVector(virN_step);
	initCommVarVector(virT_step);
	initCommVarVector(virX_step);
	initCommVarVector(virY_step);
	initCommVarVector(virZ_step);
	initCommVarVector(velocityN_step);

	for (auto pit = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
		const double distCenter_x = pit->r(0) - _globalCenter[0];
		const double distCenter_y = pit->r(1) - _globalCenter[1];
		const double distCenter_z = pit->r(2) - _globalCenter[2];
		const double distCenterSquared = std::pow(distCenter_x,2) + std::pow(distCenter_y,2) + std::pow(distCenter_z,2);
		const double distCenter = std::sqrt(distCenterSquared);
		// const double distMaxSquared = _distMax * _distMax;

		const unsigned int index = std::min(_nShells, static_cast<unsigned int>((distCenter/_distMax)*_nShells));  // Shell index of molecule 

		numMolecules_step.local[index] ++;

		const double u_x = pit->v(0);
		const double u_y = pit->v(1);
		const double u_z = pit->v(2);
		const double vi_x = pit->Vi(0);
		const double vi_y = pit->Vi(1);
		const double vi_z = pit->Vi(2);
		const double vi_n = pit->ViN();
		const double vi_t = pit->ViT();

		velocityN_step.local[index] = (distCenter_x*u_x + distCenter_y*u_y + distCenter_z*u_z)/distCenter;  // Radial 

		virN_step.local[index] += vi_n;
		virT_step.local[index] += vi_t;
		virX_step.local[index] += vi_x;
		virY_step.local[index] += vi_y;
		virZ_step.local[index] += vi_z;
	}

// Gather quantities. Note: MPI_Reduce instead of MPI_Allreduce! Therefore, only root has correct values
/* #ifdef ENABLE_MPI //TODO: MPI
	MPI_Reduce(virN_step.local.data(), virN_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(virT_step.local.data(), virT_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(virX_step.local.data(), virX_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(virY_step.local.data(), virY_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(virZ_step.local.data(), virZ_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(velocityN_step.local.data(), velocityN_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(numMolecules_step.local.data(), numMolecules_step.global.data(), _lenVector, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else */
	for (unsigned long i = 0; i < _lenVector; i++) {
		virN_step.global[i] = virN_step.local[i];
		virT_step.global[i] = virT_step.local[i];
		virX_step.global[i] = virX_step.local[i];
		virY_step.global[i] = virY_step.local[i];
		virZ_step.global[i] = virZ_step.local[i];
		velocityN_step.global[i] = velocityN_step.local[i];
		numMolecules_step.global[i] = numMolecules_step.local[i];
	}
/* #endif */

	//ADD TO ACCUM VECTORS
	for (int i = 0; i < _lenVector; i++) //TODO: use stl-algorithm instead of loop!
	{
		_VirN_accum[i] +=virN_step.global[i];
		_VirT_accum[i] +=virT_step.global[i];
		_VirX_accum[i] +=virX_step.global[i];
		_VirY_accum[i] +=virY_step.global[i];
		_VirZ_accum[i] +=virZ_step.global[i];
		_velocityN_accum[i] += velocityN_step.global[i];
		_numMolecules_accum_output[i] += numMolecules_step.global[i];
	}
	


	if(simstep%_writeFreq == 0){
        if (domainDecomp->getRank() == 0) {
            // Write output file
            std::stringstream ss;
            ss << std::setw(9) << std::setfill('0') << simstep;
            const std::string fname = _outputPrefix+ss.str()+".dat";
            std::ofstream ofs;
            ofs.open(fname, std::ios::out);
            ofs << std::setw(24) << "radius_(lower)"     	// Bin position (radius)
                << std::setw(24) << "density_step"   // Average number of molecules in bin per step
                << std::setw(24) << "_density_avg"       	// 
                << std::setw(24) << "_density_avg_fitted"   // 
                << std::setw(24) << "_virtual_density"      // 
               
			    << std::setw(24) << "_FcorrectionShell"        	// 
                << std::setw(24) << "_UcorrectionShell"        	// 
                << std::setw(24) << "_VirNcorrectionShell"        	// 
                << std::setw(24) << "_VirTcorrectionShell"        	// 

                << std::setw(24) << "p_(sphVir)"  // Pressure (using spherical vir)
                << std::setw(24) << "p_(sphVir)_direct"  // Pressure (using spherical vir)
                << std::setw(24) << "p_(xyzVir)_direct"   // Pressure (using xyz vir)
                << std::setw(24) << "p_n"        // Pressure in normal direction; 
                << std::setw(24) << "p_t"        // Pressure in tangantial direction; 
                << std::setw(24) << "p_n_direct"        // Pressure in normal direction; 
                << std::setw(24) << "p_t_direct"        // Pressure in tangantial direction; 
                << std::setw(24) << "p_n_corr"        // Pressure in normal direction; 
                << std::setw(24) << "p_t_corr"        // Pressure in tangantial direction; 

                << std::setw(24) << "VirN"          // virial -- for testing
                << std::setw(24) << "VirT"          // virial -- for testing
                << std::setw(24) << "VirN_direct"          // virial -- for testing
                << std::setw(24) << "VirT_direct"          // virial -- for testing
                << std::setw(24) << "VirN_corr"          // virial -- for testing
                << std::setw(24) << "VirT_corr"          // virial -- for testing
               
			    << std::setw(24) << "VirX_direct"          // virial -- for testing
                << std::setw(24) << "VirY_direct"          // virial -- for testing
                << std::setw(24) << "VirZ_direct"          // virial -- for testing

                // << std::setw(24) << "T"             // Temperature, assuming that there is no drift
                // << std::setw(24) << "T_driftcorr"   // Temperature without drift (i.e. "real" temperature)
                // << std::setw(24) << "v_r"        // Drift velocity in radial direction

                << std::setw(24) << "numMols_step"   // Average number of molecules in bin per step
                << std::setw(24) << "_numMolecules_accum"   // 
                << std::setw(24) << "_numMolsaccum_output"   // 
                << std::setw(24) << "_globalNumMols(const)"        	// 
                << std::setw(24) << "shellVolume"        	// 
				;
                // << std::setw(24) << "T_n"        // Temperature in radial direction
                // << std::setw(24) << "T_t"        // Temperature in tangantial direction
                // << std::setw(24) << "numSamples";    // Number of samples (<= _writeFrequency)
            ofs << std::endl;

            //this implementation is probaobly bad, but meant to be only temporary (#todo):
            std::vector<double> shell_centralRadius;
            shell_centralRadius.resize(_lenVector);
            for(int i = 1; i < _lenVector; i++){
                shell_centralRadius[i-1] = (_shellLowerBound[i-1] + _shellLowerBound[i])/2.;
            }
            shell_centralRadius[_lenVector-1] = (_shellLowerBound[_lenVector-1] + (_distMax/_nShells)); // not very precise
            // \end bad implementation

            for (unsigned long i = 0; i < _lenVector; i++) {
                // unsigned long numSamples {0ul};
                // double numMolsPerStep {std::nan("0")}; // Not an int as particles change bin during simulation
                // double T {std::nan("0")};
                // double T_driftcorr {std::nan("0")};
                // double ekin {std::nan("0")};
                // double p_xyz {std::nan("0")};
                // double p_sph {std::nan("0")};
                // double v_r {std::nan("0")};
                // double vir_n {std::nan("0")};
                // double vir_t {std::nan("0")};
                // double vir_n_direct {std::nan("0")};
                // double vir_t_direct {std::nan("0")};
                // double vir_n_corr {std::nan("0")};
                // double vir_t_corr {std::nan("0")};
                // double density_step {std::nan("0")};
                 
				// double v_drift_squared = v_r*v_r;  // is this reasonable at all?
				double T = _T; //TODO: hier die gesamttemperatur zu nehmen ist sehr bedenklich! korrekt wäre die temperatur, die tatsächlich in der sehll herrscht!
				// T           = (2*_ekin_accum[i]) / _doftotal_accum[i];
				// T_driftcorr = (2*_ekin_accum[i] - v_drift_squared*_mass_accum[i]) / _doftotal_accum[i];
				// ekin        = _ekin_accum[i] / numMols_accum;

				double density_step = numMolecules_step.global[i]/_shellVolume[i];
				double density_avg = _density_avg[i];

				double vir_n = _VirN_accum[i]/_numMolecules_accum_output[i];
				double vir_t = _VirT_accum[i]/_numMolecules_accum_output[i];
				double vir_n_corr = _VirNcorrectionShell[i];
				double vir_t_corr = _VirTcorrectionShell[i];
				double vir_n_direct = vir_n - vir_n_corr;
				double vir_t_direct = vir_t - vir_t_corr;

				double vir_x_direct = _VirX_accum[i]/_numMolecules_accum_output[i];
				double vir_y_direct = _VirY_accum[i]/_numMolecules_accum_output[i];
				double vir_z_direct = _VirZ_accum[i]/_numMolecules_accum_output[i];

				double p_n         = density_avg * (T + vir_n); 
				double p_t         = density_avg * (T + vir_t);
				double p_n_direct  = density_avg * (T + vir_n_direct);
				double p_t_direct  = density_avg * (T + vir_t_direct);
				double p_n_corr 	= density_avg * (T + vir_n_corr);
				double p_t_corr 	= density_avg * (T + vir_t_corr);

				double p_sph       = density_avg * ( T + (vir_n + 2.*vir_t)/(3.) );
				double p_sph_corr  = density_avg * ( T + (vir_n_corr + 2.*vir_t_corr)/(3.) );
				double p_sph_direct= density_avg * ( T + (vir_n_direct + 2.*vir_t_direct)/(3.) );
				double p_xyz_direct= density_avg * ( (_VirX_accum[i]+_VirY_accum[i]+_VirZ_accum[i])/(3.0*_numMolecules_accum_output[i]) + T);


                         
                //FOR TESTING ONNLY (could lead to crashes (or undef. bhv?))
                const double numMols_accum = static_cast<double>(_numMolecules_accum[i]);
                // \END for testing only

                ofs << FORMAT_SCI_MAX_DIGITS << _shellLowerBound[i]  // Radius bin
                    << FORMAT_SCI_MAX_DIGITS << density_step
                    << FORMAT_SCI_MAX_DIGITS << _density_avg[i]
                    << FORMAT_SCI_MAX_DIGITS << _density_avg_fitted[i]
                    << FORMAT_SCI_MAX_DIGITS << _virtual_density[i]

                    << FORMAT_SCI_MAX_DIGITS << _FcorrectionShell[i]
                    << FORMAT_SCI_MAX_DIGITS << _UcorrectionShell[i]
                    << FORMAT_SCI_MAX_DIGITS << _VirNcorrectionShell[i]
                    << FORMAT_SCI_MAX_DIGITS << _VirTcorrectionShell[i]

                    << FORMAT_SCI_MAX_DIGITS << p_sph
                    << FORMAT_SCI_MAX_DIGITS << p_sph_direct
                    << FORMAT_SCI_MAX_DIGITS << p_xyz_direct
                    << FORMAT_SCI_MAX_DIGITS << p_n
                    << FORMAT_SCI_MAX_DIGITS << p_t
                    << FORMAT_SCI_MAX_DIGITS << p_n_direct
                    << FORMAT_SCI_MAX_DIGITS << p_t_direct
                    << FORMAT_SCI_MAX_DIGITS << p_n_corr
                    << FORMAT_SCI_MAX_DIGITS << p_t_corr

                    << FORMAT_SCI_MAX_DIGITS << vir_n
                    << FORMAT_SCI_MAX_DIGITS << vir_t
                    << FORMAT_SCI_MAX_DIGITS << vir_n_direct
                    << FORMAT_SCI_MAX_DIGITS << vir_t_direct
                    << FORMAT_SCI_MAX_DIGITS << vir_n_corr
                    << FORMAT_SCI_MAX_DIGITS << vir_t_corr

                    << FORMAT_SCI_MAX_DIGITS << vir_x_direct
                    << FORMAT_SCI_MAX_DIGITS << vir_y_direct
                    << FORMAT_SCI_MAX_DIGITS << vir_z_direct

                    // << FORMAT_SCI_MAX_DIGITS << T
                    // << FORMAT_SCI_MAX_DIGITS << T_driftcorr
                    // << FORMAT_SCI_MAX_DIGITS << v_r

                    << FORMAT_SCI_MAX_DIGITS << numMolecules_step.global[i]
                    << FORMAT_SCI_MAX_DIGITS << _numMolecules_accum[i]
                    << FORMAT_SCI_MAX_DIGITS << _numMolecules_accum_output[i]
                    << FORMAT_SCI_MAX_DIGITS << _globalNumMols
                    << FORMAT_SCI_MAX_DIGITS << _shellVolume[i]
                    << std::endl;
            }
            ofs.close();
        }

        // Reset vectors to zero
		std::fill(_VirN_accum.begin(), _VirN_accum.end(), 0.); 				 
		std::fill(_VirT_accum.begin(), _VirT_accum.end(), 0.); 				 
		std::fill(_VirX_accum.begin(), _VirX_accum.end(), 0.); 				 
		std::fill(_VirY_accum.begin(), _VirY_accum.end(), 0.); 				 
		std::fill(_VirZ_accum.begin(), _VirZ_accum.end(), 0.); 				 
		std::fill(_velocityN_accum.begin(), _velocityN_accum.end(), 0.); 				 
		std::fill(_numMolecules_accum_output.begin(), _numMolecules_accum_output.end(), 0.); 				 








		
		
		
		
		
		


	}
}









void Spherical::calculationV1_fromPaper_withIsabelsPrefactorAndRlow_repParticles(){
/* 
	Implementation according to paper, using representative Particles to calculate Shellwise corrections.
	chagnes:
		-- used additiona prefactor *distanceFromCenter/4. that is only found in isabels code, not in the paper
		-- used rlow from isabels code, not from paper.
	>> incorrect results. droplets get destroyed, FCorr seems to be broken.
 */

	//TODO: move to object? (only calc once):
	double rCutoff = _cutoffLJ;	//TODO: find correct rc
	double rm = 16.; 		//TODO: find value for rm (hardcoded is quatsch)

	double epsilon = 1.; 	//TODO: find correct epsilon
	double sigma = 1.;		//TODO: find correct sigma
	double sigma6 = std::pow(sigma, 6); 

	double generalPrefactor =  2.*M_PI*epsilon*sigma6;




	/* method: 
		*		to determine the correction that applies to a particle in shell k,
		*		we use one representative particle per shell.
		*		the force / energy / viral that this particle experiences from the actual (not representative) density distribution
		*		is the correction, that will be applied to each particle of that shell in each time step.
		*		 
		*		(note: the method used by the legacy code followed the same approach, only using the actual particles per shell instead of representative ones.
		*		 this produced many unecessary calculations, lead to problems if there was a shell without an actual particle in it.
		*		 furthermore, the implemention was false (confusing influence a shell has on a particle with the influence a particle has on a shell),
		*		 leading to false results.)
		*/
	// LOOP OVER REPRESENTATIVE MOLECULES:
	for (int i = 0; i < _nShells+1; i++) { // i is the index of the shell that the representative mol lives in
		double _shellWidth = _distMax/_nShells;
		double distCenter_x = _shellWidth*(i+.5); // postition of representative molecule
		double distCenter_y = 0.;
		double distCenter_z = 0.;
		double distCenter2 = distCenter_x*distCenter_x + distCenter_y*distCenter_y + distCenter_z*distCenter_z;
		double distCenter = std::sqrt(distCenter2);
		double distCenterInv = 1./(distCenter);
		double distCenterInv2 = distCenterInv*distCenterInv;
		double distCenterInv3 = distCenterInv2*distCenterInv;

		//MOLECULE CORRECTION BY MOL:
		double uCorrectionMol    {0.};
		double fCorrectionMol    {0.};
		double virNCorrectionMol {0.};
		double virTCorrectionMol {0.};

			//mir ist unklar, warum jeweils *(distCenter/4.) gerechnet werden sollte.  entsprciht meiner ansicht nach nicht dem paper, scheint aber nur so zu funktioinieren, und isabel hats auch gemacht (zumindest bei fCorr, die anderen müsste man nochmal checken).
		double prefactorU 	 = (distCenter/4.) * generalPrefactor						;
		double prefactorF 	 = (distCenter/4.) * 2.*generalPrefactor * distCenterInv2 	;
		double prefactorVirN = (distCenter/4.) * generalPrefactor * distCenterInv3 		; 
		double prefactorVirT = (distCenter/4.) * 4.*generalPrefactor * distCenterInv 	;

		//LOOP OVER SHELLS:
		for(int k = 0; k < _nShells; k++){ //using k as shell index for consitency with Nitzke2021, equation (14)
			//SHELL PROPERTIES:
			double Rk = _shellLowerBound[k+1]; //TODO: this is the definition of Rk that Nitzke2021 used in their implementation
			double Rk2 = Rk*Rk;
			double shellWidth = _distMax/_nShells; //TODO: calculation is performed for each molecule over again. perhaps rather save in Spherical::-object

			//MOLECULE<->SHELL INTERACTION RADIUS BOUNDS:
			// double rLowerBound = rCutoff < std::abs(distCenter-Rk) ? distCenter-Rk : rCutoff;
			// double rLowerBound = std::max(std::abs(distCenter-Rk), rCutoff); //THIS SEEMS TO BE THE ROOT OF ALL EVIL
			
			
			
			double drShells = _distMax/(_nShells+1);  //WHY: warum +1? und müsste ich dann hier sogar +2 rechnen? (weil _nShells bei mir nur die shells zählt, und bei isabel die außenhülle einberechnet wird?)
			double deltaShells = rCutoff / drShells;
			double realk = distCenter / drShells;
			double lowerS = std::min(static_cast<double>(std::floor((realk - deltaShells))), static_cast<double>(_nShells + 1));
			double interS = std::max(static_cast<double>(std::ceil(abs(realk - deltaShells))), 1.0);
			double upperS = std::min(static_cast<double>(std::ceil(realk + deltaShells)), static_cast<double>(_nShells + 2));

			std::cout << "rep mol i " << i << ", k = "<< k<< ", lowerS = "<< lowerS<< std::endl;
			double rLowerBound;
			if(k < (lowerS + 1)){
				std::cout << "is   I" << std::endl;
				rLowerBound = distCenter - _shellLowerBound[k];
			} else if(k >= interS && k < upperS){
				std::cout << "is  II" << std::endl;
				rLowerBound = rCutoff; 
			} else if(k >= upperS && k < (_nShells + 2)){
				std::cout << "is III" << std::endl;
				rLowerBound = _shellLowerBound[k] - distCenter;
			} else {
				std::cout << "continue for mol at ksi = "<< distCenter<< ", k = "<< k<<", lowerS = " << lowerS<< std::endl;
				continue;
			}
			if (rLowerBound > rm){ continue; } // according to isabels code

			double rLowerBound2 = rLowerBound*rLowerBound;
			double rLowerBoundInv2 = std::pow(rLowerBound, -2);
			double rLowerBoundInv4 =  rLowerBoundInv2 * rLowerBoundInv2;
			double rLowerBoundInv6 =  rLowerBoundInv2 * rLowerBoundInv4;
			double rLowerBoundInv8 =  rLowerBoundInv4 * rLowerBoundInv4;
			double rLowerBoundInv10 = rLowerBoundInv4 * rLowerBoundInv6;
			double rLowerBoundInv12 = rLowerBoundInv6 * rLowerBoundInv6;

			double rUpperBound = std::min(distCenter+Rk, rm);
			// double rUpperBound = rm < std::abs(distCenter+Rk) ? rm : std::abs(distCenter+Rk);

			double rUpperBound2 = rUpperBound*rUpperBound;
			double rUpperBoundInv2 = std::pow(rUpperBound, -2);
			double rUpperBoundInv4 = rUpperBoundInv2 * rUpperBoundInv2;
			double rUpperBoundInv6 = rUpperBoundInv2 * rUpperBoundInv4;
			double rUpperBoundInv8 = rUpperBoundInv4 * rUpperBoundInv4;
			double rUpperBoundInv10 = rUpperBoundInv4 * rUpperBoundInv6;
			double rUpperBoundInv12 = rUpperBoundInv6 * rUpperBoundInv6;

			//U CORRECTION (CONTRIBUTION OF CURRENT MOL TO TO SHELL K):
			double term1=
								_virtual_density[k]*Rk*shellWidth 
								* (
								.4 * sigma6 * (std::pow((distCenter + Rk), -10) - std::pow(rLowerBound, -10) ) 
								- (std::pow((distCenter + Rk), -4)  - std::pow(rLowerBound, -4))
								) / distCenter;
			uCorrectionMol += term1;

			//F CORRECTION (CONTRIBUTION OF SHELL K TO CURRENT MOL):
			double term2=
								_virtual_density[k]*Rk*shellWidth 
								*(
								sigma6 * ((1.2*rUpperBound2 + distCenter2 - Rk2)*rUpperBoundInv12 - (1.2*rLowerBound2 + distCenter2 - Rk2)*rLowerBoundInv12)	
								- 		 ((1.5*rUpperBound2 + distCenter2 - Rk2)*rUpperBoundInv6  - (1.5*rLowerBound2 + distCenter2 - Rk2)*rLowerBoundInv6)
								);
			fCorrectionMol += term2;	

			//VIRIAL CORRECTION (CONTRIBUTION OF SHELL K TO CURRENT MOL):
			double term3=
									_virtual_density[k]*Rk*shellWidth 
									*(
									(1.5 * sigma6 * (rUpperBoundInv8 - rLowerBoundInv8) - 3. * (rUpperBoundInv2 - rLowerBoundInv2))
									+ 2.*(distCenter2- Rk2) * (1.2 * sigma6*(rUpperBoundInv10 - rLowerBoundInv10) - 1.5*(rUpperBoundInv4 - rLowerBoundInv4))
									+ std::pow((distCenter2- Rk2), 2) * (sigma6*(rUpperBoundInv12 - rLowerBoundInv12) - (rUpperBoundInv6 - rLowerBoundInv6))
									);
			virNCorrectionMol +=term3;


			double term4=
									_virtual_density[k]*Rk*shellWidth 
									*(
									1.2 * sigma6*(rUpperBoundInv10-rLowerBoundInv10) 
									-1.5 * (rUpperBoundInv4 - rLowerBoundInv4)
									); //hier fehlt noch die subtraktion des VirNCorrectionMol Terms (wird nach dem loop gemacht)
			virTCorrectionMol += 	term4;

		}
		

		uCorrectionMol 		*= prefactorU;
		fCorrectionMol 		*= prefactorF;
		virNCorrectionMol 	*= prefactorVirN;
		virTCorrectionMol 	= prefactorVirT * virTCorrectionMol - virNCorrectionMol;

		_UcorrectionShell[i]	= uCorrectionMol; 	
		_FcorrectionShell[i]	= fCorrectionMol;	
		_VirNcorrectionShell[i]	= .5* virNCorrectionMol; 
		_VirTcorrectionShell[i]	= .5* virTCorrectionMol;

		// //FOR TESTING!!: 
		// _VirNcorrectionShell[i] = fCorrectionMol; 
		// _VirTcorrectionShell[i] = uCorrectionMol; 

	}
}













void Spherical::calculationV2_isabelsCode___abandoned(){

/* 
	Reverse engineered isabels code / made it fit into the new Spherical.cpp. 
	abandoned in favor of V3, which yields same results, but better (less ragged & more efficient due to representative particles)
 */

    std::vector<double> UShells_Mean_local;
    std::vector<double> FShells_Mean_local;
    std::vector<double> PNShells_Mean_local;
    std::vector<double> PTShells_Mean_local;
	UShells_Mean_local.resize(_lenVector);
	FShells_Mean_local.resize(_lenVector);
	PNShells_Mean_local.resize(_lenVector);
	PTShells_Mean_local.resize(_lenVector);
	std::fill(UShells_Mean_local.begin(), UShells_Mean_local.end(), 0.);
	std::fill(FShells_Mean_local.begin(), FShells_Mean_local.end(), 0.);
	std::fill(PNShells_Mean_local.begin(), PNShells_Mean_local.end(), 0.);
	std::fill(PTShells_Mean_local.begin(), PTShells_Mean_local.end(), 0.);

	 
	
	

	//TODO: move to object? (only calc once):
	const double rCutoff = _cutoffLJ;	//TODO: find correct rc
	const double rm = 16.; 		//TODO: find value for rm (hardcoded is quatsch)  // == rcmax?

	const double epsilon = 1.; 	//TODO: find correct epsilon
	const double sigma = 1.;		//TODO: find correct sigma
	const double sigma6 = std::pow(sigma, 6); 

	const double generalPrefactor =  2.*M_PI*epsilon*sigma6;


	const int calcFreq = _calcFreq; // redundand, just a quickfix
	const uint64_t simstep = _simulation.getSimulationStep();


	if ((simstep) % _calcFreq == 0) {

		// U Correction of homogeneous system  for one component
		double UCORR = _rho_out*(8./3.)*M_PI*(1./(3.*std::pow(rCutoff,9))-1./std::pow(rCutoff,3));
    	double PCORR = _rho_out*(16./3.)*M_PI*(2./(3.*std::pow(rCutoff,9))-1./std::pow(rCutoff,3));
		// global_log->info() << "[Long Range Correction] Homogeneous term: rho_out = " << _rho_out << " UpotConstKorrLJ = " << UpotConstKorrLJ << " ; VirialConstKorrLJ = " << VirialConstKorrLJ << std::endl;
		global_log->info() << "[Long Range Correction] Alt. homog. term: rho_out = " << _rho_out << " UpotConstKorrLJ = " << UCORR      << " ; VirialConstKorrLJ = " << PCORR << std::endl;







		// Korrektur je Schale
		// std::fill(_FcorrectionShell.begin(), _FcorrectionShell.end(), 0.);		 
		// std::fill(_UcorrectionShell.begin(), _UcorrectionShell.end(), 0.);		 
		// std::fill(_VirNcorrectionShell.begin(), _VirNcorrectionShell.end(), 0.);		 
		// std::fill(_VirTcorrectionShell.begin(), _VirTcorrectionShell.end(), 0.);		





		double rlow, rlow2, rlowInv, rlowInv2, rdash2, rdashInv, UCorrTemp, rdash, rdashInv2, FCorrTemp, PNCorrTemp,
			PTCorrTemp;
		double ksi2 = 0.0;

		for (auto mol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); mol.isValid(); ++mol) {
				double distCenter_x = mol->r(0) - _globalCenter[0];
				double distCenter_y = mol->r(1) - _globalCenter[1];
				double distCenter_z = mol->r(2) - _globalCenter[2];
				double distCenter2 = distCenter_x*distCenter_x + distCenter_y*distCenter_y + distCenter_z*distCenter_z;
				double distCenter = std::sqrt(distCenter2);

				double drShells = _distMax/(_nShells+1);  //WHY: warum +1? und müsste ich dann hier sogar +2 rechnen? (weil _nShells bei mir nur die shells zählt, und bei isabel die außenhülle einberechnet wird?)
				double deltaShells = rCutoff / drShells;
				
				double realk = distCenter / drShells;
				unsigned long k = std::round(realk);

				double lowerS = std::min(static_cast<double>(std::floor((realk - deltaShells))), static_cast<double>(_nShells + 1));
				double interS = std::max(static_cast<double>(std::ceil(abs(realk - deltaShells))), 1.0);
				double upperS = std::min(static_cast<double>(std::ceil(realk + deltaShells)), static_cast<double>(_nShells + 2));



			unsigned long molID = mol->getID();

			double UCorrShells = 0.0;
			double FCorrShells = 0.0;
			double PNCorrShells = 0.0;
			double PTCorrShells = 0.0;


			//stuff from multi-site-generalization
			unsigned int ci = 0;
			unsigned int cj = 0;
			ParaStrm& params = _domain->getComp2Params()(ci, cj);
			params.reset_read();
			unsigned si = 0;
			unsigned sj = 0; 

			double eps24;
			double sigma;
			double sig2;
			double shift6;
			double eps;
			params >> eps24;
			params >> sig2;
			params >> shift6;
			// sigma = sqrt(sig2);
			// eps = eps24 / 24;
			double tau1 = 0.;
			double tau2 = 0.;
			// double tau1 = 0.5;
			// double tau2 = 0;
			double sigma6 = sig2 * sig2 * sig2;
			double factorU = -M_PI * drShells * eps24 * sigma6 / (6. * distCenter);
			double factorF = -factorU / distCenter;
			double factorP = 0.5 * factorF / distCenter;

			// RShells[i] == _shellLowerBound[i+1] !!
			const double rcmax = rm;
			std::cout << " ---- MOLECULE " << molID << ": lowerS = " << lowerS <<", interS = " << interS << ", upperS = " <<upperS << std::endl;
			for (unsigned long j = 1; j < (_lenVector + 1); j++) {
				double shellLowerBound2 = _shellLowerBound[j] * _shellLowerBound[j];

                /* int sumOfTrue = 0;
				if(j<(lowerS + 1)){
                    rlow = distCenter - _shellLowerBound[j];
                    sumOfTrue++;
                } 
				if(j >= interS && j < upperS){
                    rlow = rCutoff; 
                    sumOfTrue++;
                }
				if(j >= upperS && j < (_nShells + 2)){
                    rlow = _shellLowerBound[j] - distCenter;
                    sumOfTrue++;
                }

                if(sumOfTrue < 1){ continue;}
                if(sumOfTrue > 1){ std::cout << "error in sphericalLRC. more than one condition true. exiting." << std::cout; _simulation.exit(0)}
				 */

				if (_density_avg_fitted[j - 1] != 0.0) {
                    if(j<(lowerS + 1)){
                        rlow = distCenter - _shellLowerBound[j];
                    } else if(j >= interS && j < upperS){
                        rlow = rCutoff; 
                    } else if(j >= upperS && j < (_nShells + 2)){
                        rlow = _shellLowerBound[j] - distCenter;
                    } else {
                        continue;
                    }
					rlowInv = 1. / rlow;
					rlowInv2 = rlowInv * rlowInv;
					rdashInv = 1. / (_shellLowerBound[j] + distCenter);
					

					UCorrTemp = sigma6 * 0.2 * (pow(rdashInv, 10) - pow(rlowInv, 10)) -
								0.5 * (pow(rdashInv, 4) - pow(rlowInv, 4));
					if (rlow < rcmax) {
						rdash = std::min(rcmax, (_shellLowerBound[j] + distCenter));
						rdashInv = 1. / rdash;
						rdashInv2 = rdashInv * rdashInv;
						FCorrTemp =
							sigma6 * ((6. / 5. * rdash * rdash + ksi2 - shellLowerBound2) *
											pow(rdashInv2, 6) -
										(6. / 5. * rlow * rlow + ksi2 - shellLowerBound2) *
											pow(rlowInv2, 6)) -
							(1.5 * rdash * rdash + ksi2 - shellLowerBound2) * pow(rdashInv2, 3) +
							(1.5 * rlow * rlow + ksi2 - shellLowerBound2) * pow(rlowInv2, 3);
						FCorrShells =
							FCorrShells + FCorrTemp * factorF * _density_avg_fitted[j - 1] * _shellLowerBound[j];
						PNCorrTemp =  // 1.5*sigma6 * (pow(rdashInv2,4) - pow(rlowInv2,4))
							-3. * (rdashInv2 - rlowInv2) +
							2. * (ksi2 - shellLowerBound2) *
								(  // 6/5*sigma6 * (pow(rdashInv2,5) - pow(rlowInv2,5))
									-1.5 * (pow(rdashInv2, 2) - pow(rlowInv2, 2))) +
							pow((ksi2 - shellLowerBound2), 2) *
								(  // sigma6 * ( pow(rdashInv2,6) - pow(rlowInv2,6))
									-(pow(rdashInv2, 3) - pow(rlowInv2, 3)));
						PNCorrShells = PNCorrShells +
										PNCorrTemp * factorP * _density_avg_fitted[j - 1] * _shellLowerBound[j];
						PTCorrTemp = 6. / 5. * sigma6 * (pow(rdashInv2, 5) - pow(rlowInv2, 5)) -
										1.5 * (pow(rdashInv2, 2) - pow(rlowInv2, 2));
						PTCorrShells = PTCorrShells + 4. * ksi2 * PTCorrTemp * factorP *
															_density_avg_fitted[j - 1] * _shellLowerBound[j];
					}
					UCorrShells =
						UCorrShells + UCorrTemp * factorU * _density_avg_fitted[j - 1] * _shellLowerBound[j];
				}
			}

			int shellOfThisMol;

				if (k >= _nShells) {
					shellOfThisMol = _nShells;
				} else if (k == 0) {
					shellOfThisMol = k;
				} else {
					shellOfThisMol = k - 1;
				}


			PTCorrShells -= PNCorrShells;
			UShells_Mean_local[shellOfThisMol] += UCorrShells;
			FShells_Mean_local[shellOfThisMol] += FCorrShells;
			PNShells_Mean_local[shellOfThisMol] += PNCorrShells;
			PTShells_Mean_local[shellOfThisMol] += PTCorrShells;
			// std::cout << shellOfThisMol<<std::endl;
		}

		// Distribution of Shell Corrections to every node
/* #ifdef ENABLE_MPI //TODO: MPI
	MPI_Allreduce(numMolecules_step.local.data(), numMolecules_step.global.data(), _lenVector, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
	*/    

		

		for (unsigned long i = 0; i < _lenVector; i++) {
			if (_density_avg[i] != 0.0) { //WHY: is this the correct density? probably not ...
			//WTF: and if density_avg = 0 ? then we just leave it ... ok cool.
				UShells_Mean_local[i]  /= (_density_avg[i] * _shellVolume[i]);
				FShells_Mean_local[i]  /= (_density_avg[i] * _shellVolume[i]);
				PNShells_Mean_local[i] /= (_density_avg[i] * _shellVolume[i]);
				PTShells_Mean_local[i] /= (_density_avg[i] * _shellVolume[i]);
			}
			
			
			// _UcorrectionShell[i] = UShells_Mean_local[i];	
			// _FcorrectionShell[i] = FShells_Mean_local[i];	
			// TESTING:
			_UcorrectionShell[i] = FShells_Mean_local[i];//FOR TESTING PTShells_Mean_local[i];//FALSE! just for output ...
		}
	//TODO: würde das hier nicht auch ohne for-loop gehen? (einfach numMolecules_step.global = numMolecules_step.local;) 
/* endif */
	}
}


void Spherical::calculationV3_isabelsMethod(){
/* 
	-- reverse engineered from isabels code
	-- uses representative particle method
	-- main difference to V1: calculation of rlow
	-- includes prefactor distanceFromCenter/4., that seems to not be in paper.
 */

    std::vector<double> UShells_Mean_local;
    std::vector<double> FShells_Mean_local;
    std::vector<double> PNShells_Mean_local;
    std::vector<double> PTShells_Mean_local;
	UShells_Mean_local.resize(_lenVector);
	FShells_Mean_local.resize(_lenVector);
	PNShells_Mean_local.resize(_lenVector);
	PTShells_Mean_local.resize(_lenVector);
	std::fill(UShells_Mean_local.begin(), UShells_Mean_local.end(), 0.);
	std::fill(FShells_Mean_local.begin(), FShells_Mean_local.end(), 0.);
	std::fill(PNShells_Mean_local.begin(), PNShells_Mean_local.end(), 0.);
	std::fill(PTShells_Mean_local.begin(), PTShells_Mean_local.end(), 0.);

	 
	
	

	//TODO: move to object? (only calc once):
	const double rCutoff = _cutoffLJ;	//TODO: find correct rc
	const double rm = 16.; 		//TODO: find value for rm (hardcoded is quatsch)  // == rcmax?

	const double epsilon = 1.; 	//TODO: find correct epsilon
	const double sigma = 1.;		//TODO: find correct sigma
	const double sigma6 = std::pow(sigma, 6); 

	const double generalPrefactor =  2.*M_PI*epsilon*sigma6;


	const int calcFreq = _calcFreq; // redundand, just a quickfix
	const uint64_t simstep = _simulation.getSimulationStep();


	if ((simstep) % _calcFreq == 0) {

		// U Correction of homogeneous system  for one component
		double UCORR = _rho_out*(8./3.)*M_PI*(1./(3.*std::pow(rCutoff,9))-1./std::pow(rCutoff,3));
    	double PCORR = _rho_out*(16./3.)*M_PI*(2./(3.*std::pow(rCutoff,9))-1./std::pow(rCutoff,3));
		// global_log->info() << "[Long Range Correction] Homogeneous term: rho_out = " << _rho_out << " UpotConstKorrLJ = " << UpotConstKorrLJ << " ; VirialConstKorrLJ = " << VirialConstKorrLJ << std::endl;
		global_log->info() << "[Long Range Correction] Alt. homog. term: rho_out = " << _rho_out << " UpotConstKorrLJ = " << UCORR      << " ; VirialConstKorrLJ = " << PCORR << std::endl;







		// Korrektur je Schale
		// std::fill(_FcorrectionShell.begin(), _FcorrectionShell.end(), 0.);		 
		// std::fill(_UcorrectionShell.begin(), _UcorrectionShell.end(), 0.);		 
		// std::fill(_VirNcorrectionShell.begin(), _VirNcorrectionShell.end(), 0.);		 
		// std::fill(_VirTcorrectionShell.begin(), _VirTcorrectionShell.end(), 0.);		





		double ksi2 = 0.0;
		
		double drShells = _distMax/(_nShells+1);  //WHY: warum +1? und müsste ich dann hier sogar +2 rechnen? (weil _nShells bei mir nur die shells zählt, und bei isabel die außenhülle einberechnet wird?)

		for (int i = 0; i < _lenVector; i++) {
			double distCenter_x = _shellLowerBound[i] + drShells/2.;
			double distCenter_y = 0;
			double distCenter_z = 0;
			double distCenter2 = distCenter_x*distCenter_x + distCenter_y*distCenter_y + distCenter_z*distCenter_z;
			double distCenter = std::sqrt(distCenter2);

			double deltaShells = rCutoff / drShells;
			
			double realk = distCenter / drShells;
			unsigned long k = std::round(realk);

			double lowerS = std::min(static_cast<double>(std::floor((realk - deltaShells))), static_cast<double>(_nShells + 1));
			double interS = std::max(static_cast<double>(std::ceil(abs(realk - deltaShells))), 1.0);
			double upperS = std::min(static_cast<double>(std::ceil(realk + deltaShells)), static_cast<double>(_nShells + 2));
			
			std::cout << " ---- REP. MOLECULE " << i << ": lowerS = " << lowerS <<", interS = " << interS << ", upperS = " <<upperS << std::endl;


			//stuff from multi-site-generalization
			unsigned int ci = 0;
			unsigned int cj = 0;
			ParaStrm& params = _domain->getComp2Params()(ci, cj);
			params.reset_read();
			unsigned si = 0;
			unsigned sj = 0; 

			double eps24;
			double sigma;
			double sig2;
			double shift6;
			double eps;
			params >> eps24;
			params >> sig2;
			params >> shift6;
			// sigma = sqrt(sig2);
			// eps = eps24 / 24;
			double tau1 = 0.;
			double tau2 = 0.;
			// double tau1 = 0.5;
			// double tau2 = 0;
			double sigma6 = sig2 * sig2 * sig2;
			double factorU = -M_PI * drShells * eps24 * sigma6 / (6. * distCenter);
			double factorF = -factorU / distCenter;
			// double factorP = 0.5 * factorF / distCenter;
			double factorVirN = 0.5 * factorF / distCenter;
			double factorVirT = 0.5 * factorF * distCenter;

			double UCorrShells = 0.0;
			double FCorrShells = 0.0;
			double VirNCorrShells = 0.0;
			double VirTCorrShells = 0.0;
			// RShells[i] == _shellLowerBound[i+1] !!
			const double rcmax = rm;
			for (unsigned long j = 1; j < (_lenVector + 1); j++) {
				double shellLowerBound2 = _shellLowerBound[j] * _shellLowerBound[j];

				double rlow;
				if (_density_avg_fitted[j - 1] != 0.0) {
                    if(j<(lowerS + 1)){
                        rlow = distCenter - _shellLowerBound[j];
                    } else if(j >= interS && j < upperS){
                        rlow = rCutoff; 
                    } else if(j >= upperS && j < (_nShells + 2)){
                        rlow = _shellLowerBound[j] - distCenter;
                    } else {
                        continue;
                    }
					double rlowInv = 1. / rlow;
					double rlowInv2 = rlowInv * rlowInv;
					double rlowInv12 = pow(rlowInv2, 6);

					double rdashInv = 1. / (_shellLowerBound[j] + distCenter);

					double UCorrTemp = sigma6 * 0.2 * (pow(rdashInv, 10) - pow(rlowInv, 10)) -
								0.5 * (pow(rdashInv, 4) - pow(rlowInv, 4));
					if (rlow < rcmax) {
						double rdash = std::min(rcmax, (_shellLowerBound[j] + distCenter));
						rdashInv = 1. / rdash;
						double rdashInv2 = rdashInv * rdashInv;
						double FCorrTemp =
							sigma6 * ((6. / 5. * rdash * rdash + ksi2 - shellLowerBound2) *
											pow(rdashInv2, 6) -
										(6. / 5. * rlow * rlow + ksi2 - shellLowerBound2) *
											pow(rlowInv2, 6)) -
							(1.5 * rdash * rdash + ksi2 - shellLowerBound2) * pow(rdashInv2, 3) +
							(1.5 * rlow * rlow + ksi2 - shellLowerBound2) * pow(rlowInv2, 3);
						FCorrShells =
							FCorrShells + FCorrTemp * factorF * _density_avg_fitted[j - 1] * _shellLowerBound[j];
						// double PNCorrTemp =  // 1.5*sigma6 * (pow(rdashInv2,4) - pow(rlowInv2,4))
						// 	-3. * (rdashInv2 - rlowInv2) +
						// 	2. * (ksi2 - shellLowerBound2) *
						// 		(  // 6/5*sigma6 * (pow(rdashInv2,5) - pow(rlowInv2,5))
						// 			-1.5 * (pow(rdashInv2, 2) - pow(rlowInv2, 2))) +
						// 	pow((ksi2 - shellLowerBound2), 2) *
						// 		(  // sigma6 * ( pow(rdashInv2,6) - pow(rlowInv2,6))
						// 			-(pow(rdashInv2, 3) - pow(rlowInv2, 3)));
						double VirNCorrTemp =  // 1.5*sigma6 * (pow(rdashInv2,4) - pow(rlowInv2,4))
							_virtual_density[j]*_shellLowerBound[j+1]*drShells 
									*(
									(1.5 * sigma6 * (pow(rdash, -8) - pow(rlow, -8)) - 3. * (rdashInv2 - rlowInv2))
									+ 2.*(distCenter2- shellLowerBound2) * (1.2 * sigma6*(pow(rdash, -10) - pow(rlow, -10)) - 1.5*(pow(rdash, -4) - pow(rlow, -4)))
									+ std::pow((distCenter2- shellLowerBound2), 2) * (sigma6*(pow(rdash, -12) - rlowInv12) - (pow(rdash, -6) - pow(rlow, -6)))
									);
						VirNCorrShells += 
										VirNCorrTemp * factorVirN * _density_avg_fitted[j - 1] * _shellLowerBound[j];
						// double PTCorrTemp = 6. / 5. * sigma6 * (pow(rdashInv2, 5) - pow(rlowInv2, 5)) -
										// 1.5 * (pow(rdashInv2, 2) - pow(rlowInv2, 2));
						double VirTCorrTemp = 
									_virtual_density[j]*_shellLowerBound[j+1]*drShells 
									*(
									1.2 * sigma6*(pow(rdash, -10) -pow(rlow, -10)) 
									-1.5 * (pow(rdash, -4) - pow(rlow, -4))
									);
						VirTCorrShells += 4. * ksi2 * VirTCorrTemp * factorVirT *_density_avg_fitted[j - 1] * _shellLowerBound[j];
					}
					UCorrShells =
						UCorrShells + UCorrTemp * factorU * _density_avg_fitted[j - 1] * _shellLowerBound[j];
				}
			}

			VirTCorrShells -= VirNCorrShells;
			_UcorrectionShell[i] = UCorrShells;
			_FcorrectionShell[i] = FCorrShells;
			_VirNcorrectionShell[i] = VirNCorrShells;
			_VirTcorrectionShell[i] = VirTCorrShells;
		}

		// Distribution of Shell Corrections to every node
/* #ifdef ENABLE_MPI //TODO: MPI
	MPI_Allreduce(numMolecules_step.local.data(), numMolecules_step.global.data(), _lenVector, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
	*/    

		/* 

		for (unsigned long i = 0; i < _lenVector; i++) {
			// if (_density_avg[i] != 0.0) { //WHY: is this the correct density? probably not ...
			// //WTF: and if density_avg = 0 ? then we just leave it ... ok cool.
			// 	UShells_Mean_local[i]  /= (_density_avg[i] * _shellVolume[i]);
			// 	FShells_Mean_local[i]  /= (_density_avg[i] * _shellVolume[i]);
			// 	PNShells_Mean_local[i] /= (_density_avg[i] * _shellVolume[i]);
			// 	PTShells_Mean_local[i] /= (_density_avg[i] * _shellVolume[i]);
			// }
			
			
			_UcorrectionShell[i] = UShells_Mean_local[i];	
			_FcorrectionShell[i] = FShells_Mean_local[i];	
			// _VirNcorrectionShell[i] = PNShells_Mean_local[i];//FALSE! just for output ...
			// _VirTcorrectionShell[i] = PTShells_Mean_local[i];//FALSE! just for output ...
		} */
	//TODO: würde das hier nicht auch ohne for-loop gehen? (einfach numMolecules_step.global = numMolecules_step.local;) 
/* endif */
	}
}











// Functions for homogeneous corrections

double Spherical::TICCu(int n, double rcutoff, double sigma2) {
	return -(pow(rcutoff, (2 * n + 3))) / (pow(sigma2, n) * (2 * n + 3));
}

double Spherical::TICSu(int n, double rcutoff, double sigma2, double tau) {
	return -(pow((rcutoff + tau), (2 * n + 3)) - pow((rcutoff - tau), (2 * n + 3))) * rcutoff /
			   (pow(4. * sigma2, n) * tau * (n + 1) * (2 * n + 3)) +
		   (pow((rcutoff + tau), (2 * n + 4)) - pow((rcutoff - tau), (2 * n + 4))) /
			   (4. * pow(sigma2, n) * tau * (n + 1) * (2 * n + 3) * (2 * n + 4));
}

double Spherical::TISSu(int n, double rcutoff, double sigma2, double tau1, double tau2) {
	double tauPlus = tau1 + tau2;
	double tauMinus = tau1 - tau2;
	return -(pow((rcutoff + tauPlus), (2 * n + 4)) - pow((rcutoff + tauMinus), (2 * n + 4)) -
			 pow((rcutoff - tauMinus), (2 * n + 4)) + pow((rcutoff - tauPlus), (2 * n + 4))) *
			   rcutoff / (8. * pow(sigma2, n) * tau1 * tau2 * (n + 1) * (2 * n + 3) * (2 * n + 4)) +
		   (pow((rcutoff + tauPlus), (2 * n + 5)) - pow((rcutoff + tauMinus), (2 * n + 5)) -
			pow((rcutoff - tauMinus), (2 * n + 5)) + pow((rcutoff - tauPlus), (2 * n + 5))) /
			   (8. * pow(sigma2, n) * tau1 * tau2 * (n + 1) * (2 * n + 3) * (2 * n + 4) * (2 * n + 5));
}

double Spherical::TICCp(int n, double rcutoff, double sigma2) { return 2 * n * TICCu(n, rcutoff, sigma2); }

double Spherical::TICSp(int n, double rcutoff, double sigma2, double tau) {
	return -(pow((rcutoff + tau), (2 * n + 2)) - pow((rcutoff - tau), (2 * n + 2))) * pow(rcutoff, (2)) /
			   (4. * pow(sigma2, n) * tau * (n + 1)) -
		   3. * TICSu(n, rcutoff, sigma2, tau);
}

double Spherical::TISSp(int n, double rcutoff, double sigma2, double tau1, double tau2) {
	double tauPlus = tau1 + tau2;
	double tauMinus = tau1 - tau2;
	return -(pow((rcutoff + tauPlus), (2 * n + 3)) - pow((rcutoff + tauMinus), (2 * n + 3)) -
			 pow((rcutoff - tauMinus), (2 * n + 3)) + pow((rcutoff - tauPlus), (2 * n + 3))) *
			   pow(rcutoff, 2) / (8. * pow(sigma2, n) * tau1 * tau2 * (n + 1) * (2 * n + 3)) -
		   3. * TISSu(n, rcutoff, sigma2, tau1, tau2);
}












//////////////////////////////////////////////////////////////

	/* 	global_log->info() 	<< "[molid:fX:fY:fZ:virN:virT:shellID] = "
					<< molID 
					<< " : "
					<< FcorrMolApprox[0]
					<< " : "
					<< FcorrMolApprox[1]
					<< " : "
					<< FcorrMolApprox[2]
					<< " : "
					<< _VirNcorrectionShell[shellIDofMol]
					<< " : "
					<< _VirTcorrectionShell[shellIDofMol]
					<< " : "
					<< shellIDofMol
					<< std::endl; */


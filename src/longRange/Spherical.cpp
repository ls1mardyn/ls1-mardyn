
#include "Domain.h"
#include "longRange/Spherical.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"

#include <vector>
#include <cmath>
#include <algorithm>

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"
#include "utils/FileUtils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using Log::global_log;

Spherical::Spherical(double /*cutoffT*/, double cutoffLJ, Domain* domain, DomainDecompBase* domainDecomposition,
		ParticleContainer* particleContainer, Simulation* simulation)
{
	global_log->info() << "Long Range Correction for spherical interfaces is used" << endl;
	 // default values
	rc = cutoffLJ;
	_domain = domain;
	_domainDecomposition = domainDecomposition;
	_particleContainer = particleContainer;
	NShells = 60;
	NSMean = 50;
	numComp = 1;
	globalNumMols = 8560;
	rcmax = 8;
}

Spherical::~Spherical() {
}

void Spherical::init()
{
	global_log->info() << "[Spherical Long Range Correction] Initializing" << endl;

	vector<Component>&  components = *_simulation.getEnsemble()->getComponents();
	numComp = components.size(); // Noch nicht unterstützt
	resizeExactly(numLJ, numComp);

	numLJSum=0;
	resizeExactly(numLJSum2, numComp);

	for (unsigned i =0; i< numComp; i++){
		numLJSum2[i]=0;
	}
	for (unsigned i =0; i< numComp; i++) {
		numLJ[i]=components[i].numLJcenters();
		numLJSum+=numLJ[i];
		for (unsigned j=i+1; j< numComp; j++){
			numLJSum2[j]+=numLJ[i];
		}
	}

	resizeExactly(ksi, globalNumMols);
	resizeExactly(FcorrX, globalNumMols); // requires refreshIDs=1 in simulation options
	resizeExactly(FcorrY, globalNumMols); // requires refreshIDs=1 in simulation options
	resizeExactly(FcorrZ, globalNumMols); // requires refreshIDs=1 in simulation options
	resizeExactly(FcorrX_global, globalNumMols); // requires refreshIDs=1 in simulation options
	resizeExactly(FcorrY_global, globalNumMols); // requires refreshIDs=1 in simulation options
	resizeExactly(FcorrZ_global, globalNumMols); // requires refreshIDs=1 in simulation options
	resizeExactly(lowerS, globalNumMols); // requires refreshIDs=1 in simulation options
	resizeExactly(interS, globalNumMols); // requires refreshIDs=1 in simulation options
	resizeExactly(upperS, globalNumMols); // requires refreshIDs=1 in simulation options

	std::fill(ksi.begin(), ksi.end(), 0.0);
	std::fill(FcorrX.begin(), FcorrX.end(), 0.0);
	std::fill(FcorrY.begin(), FcorrY.end(), 0.0);
	std::fill(FcorrZ.begin(), FcorrZ.end(), 0.0);
	std::fill(FcorrX_global.begin(), FcorrX_global.end(), 0.0);
	std::fill(FcorrY_global.begin(), FcorrY_global.end(), 0.0);
	std::fill(FcorrZ_global.begin(), FcorrZ_global.end(), 0.0);
	std::fill(lowerS.begin(), lowerS.end(), 0.0);
	std::fill(interS.begin(), interS.end(), 0.0);
	std::fill(upperS.begin(), upperS.end(), 0.0);

	resizeExactly(RShells, NShells);
	resizeExactly(RShells2, NShells);
	resizeExactly(VShells, NShells);
	resizeExactly(rhoShellsTemp, NShells);
	resizeExactly(rhoShellsMean, NShells*NSMean);
	resizeExactly(rhoShellsAvg, NShells);
	resizeExactly(rhoShells, NShells);
	resizeExactly(rhoShells_global, NShells);
	resizeExactly(rhoShellsT, NShells);
	std::fill(RShells.begin(), RShells.end(), 0.0);
	std::fill(RShells2.begin(), RShells2.end(), 0.0);
	std::fill(VShells.begin(), VShells.end(), 0.0);
	std::fill(rhoShellsTemp.begin(), rhoShellsTemp.end(), 0.0);
	std::fill(rhoShellsMean.begin(), rhoShellsMean.end(), 0.0);
	std::fill(rhoShellsAvg.begin(), rhoShellsAvg.end(), 0.0);
	std::fill(rhoShells.begin(), rhoShells.end(), 0.0);
	std::fill(rhoShells_global.begin(), rhoShells_global.end(), 0.0);
	std::fill(rhoShellsT.begin(), rhoShellsT.end(), 0.0);

	boxlength[0]=_domain->getGlobalLength(0);
	boxlength[1]=_domain->getGlobalLength(0);
	boxlength[2]=_domain->getGlobalLength(2);
	for(unsigned short d=0;d<3;++d) { systemcenter[d] = 0.5*boxlength[d]; }
	
	_drShells = 0.5 * std::min({boxlength[0], boxlength[1], boxlength[2]}) / (NShells + 1);

	_deltaShells =  rc / _drShells;

	double drShells05 = 0.5*_drShells;
	double maxShells = NShells + _deltaShells;

	for (unsigned int i=0; i< NShells-1; i++){
		RShells[i] = (i+1) * _drShells;
		RShells2[i] = RShells[i]*RShells[i];
		VShells[i] = (4./3.)*M_PI * ( pow((RShells[i]+drShells05),3) - pow((RShells[i]-drShells05),3) );
	}
	RShells[NShells-1] = NShells * _drShells;
	RShells2[NShells-1] = RShells[NShells-1] * RShells[NShells-1];
	VShells[NShells-1] = (boxlength[0]*boxlength[1]*boxlength[2]) - (4./3.)*M_PI * pow((RShells[NShells-1]-drShells05),3);
}

void Spherical::readXML(XMLfileUnits& xmlconfig)
{
	xmlconfig.getNodeValue("shells", NShells);
	//xmlconfig.getNodeValue("frequency", frequency);
	global_log->info() << "Long Range Correction: using " << NShells << " slabs for profiles to calculate LRC." << endl;
	
	// init
	this->init();
}

void Spherical::calculateLongRange(){

	global_log->info() << "[Spherical Long Range Correction] Correcting" << endl;

	uint64_t simstep = _simulation.getSimulationStep();

	for (auto tempMol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid(); ++tempMol) {
		unsigned long molID = tempMol->getID();
		for(unsigned short d=0;d<3;++d) {
			FcorrX[molID] = systemcenter[0] - tempMol->r(0);
			FcorrY[molID] = systemcenter[1] - tempMol->r(1);
			FcorrZ[molID] = systemcenter[2] - tempMol->r(2);
		}
		ksi[molID] = std::sqrt((tempMol->r(0)-systemcenter[0])*(tempMol->r(0)-systemcenter[0]) 
							+ (tempMol->r(1)-systemcenter[1])*(tempMol->r(1)-systemcenter[1]) 
							+ (tempMol->r(2)-systemcenter[2])*(tempMol->r(2)-systemcenter[2]));
		ksi[molID] = sqrt(ksi[molID]);
		double realk = ksi[molID]/_drShells;
		lowerS[molID] = std::min( static_cast<double>(std::floor((realk-_deltaShells))), static_cast<double>(NShells));
		interS[molID] = std::max( static_cast<double>(std::ceil(abs(realk-_deltaShells))), 1.0 );
		upperS[molID] = std::max( static_cast<double>(std::ceil(realk+_deltaShells)), static_cast<double>(NShells+1));
		unsigned long k = std::round( realk );
		if (k >= NShells) {
			rhoShellsTemp[NShells] += 1.;
		} else {
			rhoShellsTemp[k] += 1.;
		}
	}
	for (unsigned int i=0; i< NShells; i++){
		rhoShellsTemp[i] = rhoShellsTemp[i] / VShells[i];
		rhoShellsAvg[i] = rhoShellsAvg[i] + rhoShellsTemp[i];
	}

	unsigned long MeanIndex = (static_cast<int>(std::floor( simstep/1000.0 ))) % NSMean;
	if ((simstep % 1000) == 0) {
		std::fill(rhoShellsMean.begin()+NShells*MeanIndex, rhoShellsMean.begin()+NShells*(MeanIndex+1), 0.0);
	}
	for (unsigned int i=0; i<NShells; i++) { rhoShellsMean[NShells*MeanIndex+i] += rhoShellsTemp[i]; }

	// mittlere Dichte
	if ( simstep % 10000 == 0) {
		std::fill(rhoShells.begin(), rhoShells.end(), 0.0);
		std::fill(rhoShells_global.begin(), rhoShells_global.end(), 0.0);

		for (unsigned int j=0; j<NShells; j++) {
			for (unsigned int i=0; i<NSMean; i++) {
				rhoShells[j] += 0.001/NSMean*rhoShellsMean[i*NShells+j];
			}
		}

		// Distribution of the Density Profile to every node
		_domainDecomposition->collCommInit(NShells);
		for (unsigned i=0; i < NShells; i++) {
			_domainDecomposition->collCommAppendDouble(rhoShells[i]);
		}
		_domainDecomposition->collCommAllreduceSum();
		for (unsigned i=0; i < NShells; i++) {
			rhoShells_global[i] = _domainDecomposition->collCommGetDouble();
		}
		_domainDecomposition->collCommFinalize();

		// rhol and rhov 
		double rhol = 0.0;
		for (unsigned int i=4; i<14; i++) {  // Werte für alle Schleifen zur rhol rhov Berechnung überprüfen, habe ich geändert wegen *0.1!
			rhol += 0.1*rhoShells[i];
		}

		for (unsigned int i=15; i<NShells; i++) {
			if ( abs(rhol-rhoShells[i]) < 0.1*rhol ) {
				rhol = ((i-1)*rhol + rhoShells[i])/i;
			}
		}

		double rhov = 0.0;
		for (unsigned int i=NShells-11; i<(NShells-1); i++) {
			rhov += 0.1*rhoShells[i];
		}

		for (unsigned int i=1; i<(NShells-20); i++) {
			unsigned int j = NShells-10-i;
			if ( abs(rhov-rhoShells[j]) < 0.2*rhov ) {
				rhov = ((10+i-1)*rhov + rhoShells[j])/(10+i);
			}
		}

		// D0 with 1090 (Baidakov et al.)
		double r10 = rhov + 0.1*(rhol-rhov);
		double r90 = rhov + 0.9*(rhol-rhov);
		double Dmin = 0.0;
		double Dmax = 0.0;
		for (unsigned int i=15; i<NShells; i++) {
			if ( rhoShells[i] > r90 ) {
				Dmin = RShells[i];
			}
		}
		for (unsigned int i=1; i<(NShells-20); i++) {
			unsigned int j = NShells-10-i;
			if ( rhoShells[j] < r10 ) {
				Dmax = RShells[j];
			}
		}

		double D0 = (Dmax-Dmin);

		// R0 
		double R0 = Dmin + 0.5*D0;

		// density profile
		for (unsigned int i=0; i<NShells; i++) {
			rhoShellsT[i] = RhoP(RShells[i], rhov, rhol, D0, R0);
		}

		// OPEN(30,FILE=TRIM(Dateiname)//"_tanh.csv",POSITION="APPEND",STATUS="OLD")
		// WRITE(30,"(F12.8,";",F12.8,";",F12.8,";",F12.8)") rhov, rhol, D0, R0
		// CLOSE(30)

		// shift for Force Correction
		for (unsigned int i=0; i<NShells; i++) {
			if ( rhoShellsT[i] >= 1.02*rhov) {
				rhoShellsT[i] -= rhov;
				double bnd = RShells[i];
			} else {
				rhoShellsT[i] = 0.0;
			}
		}

		// OPEN(39,FILE=TRIM(Dateiname)//"_T.csv",POSITION="APPEND",STATUS="OLD")
		// DO i=1, NShells, 1  
		// WRITE(39,"(F12.8,";",F12.8)") ksi(i), rhoShellsT(i,1) 
		// } 
		// CLOSE(39)

		// U Correction of homogeneous system
		// TODO: multi component
		// rhovtot = 0.0
		// DO i=1,AnzKomp,1
		//   rhovtot = rhovtot + rhov(i)
		// }
		double UpotKorrLJ = 0.0;
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
						//double tau1 = eLong[numLJSum2[ci]+si];  // ÜBERPRÜFEN
						//double tau2 = eLong[numLJSum2[cj]+sj];
						double tau1 = 0.0;
						double tau2 = 0.0;
						// double tau1 = SQRT(DOT_PRODUCT( rSiteBody(:,Si,i),rSiteBody(:,Si,i) ));   // ???? epsilon, sigma2, rsitebody, etc
						// double tau2 = SQRT(DOT_PRODUCT( rSiteBody(:,Sj,j),rSiteBody(:,Sj,j) ));
						if ( (tau1 == 0.0) && (tau2 == 0.0) ) {
						UpotKorrLJ = UpotKorrLJ
							+ rhov * (eps24/6.0)
							* (  TICCu(-6,rc,sig2)
								- TICCu(-3,rc,sig2) );
						// VirialKorrLJ = VirialKorrLJ&
						//     * (eps24/6.0)&
						//     * (  TICCp(-6,rc,sig2)&
						//        - TICCp(-3,rc,sig2) )
						}
						if ( (tau1 == 0.0) || (tau2 == 0.0) ) {
						double tau = max(tau1,tau2);   // Filtert den Wert Null raus
						UpotKorrLJ = UpotKorrLJ
							* (eps24/6.0)
							* (  TICSu(-6,rc,sig2,tau)
								- TICSu(-3,rc,sig2,tau) );
						// VirialKorrLJ = VirialKorrLJ&
						//     * (eps24/6.0)&
						//     * (  TICSp(-6,rc,sig2,tau)&
						//        - TICSp(-3,rc,sig2,tau) )
						}
						if ( (tau1 != 0.0) && (tau2 != 0.0) ) {
						UpotKorrLJ = UpotKorrLJ
							* (eps24/6.0)
							* (  TISSu(-6,rc,sig2,tau1,tau2)
								- TISSu(-3,rc,sig2,tau1,tau2) );
						// VirialKorrLJ = VirialKorrLJ&
						//     * (eps24/6.0)&
						//     * (  TISSp(-6,rc,sig2,tau1,tau2)&
						//        - TISSp(-3,rc,sig2,tau1,tau2) )
						}
					}
				}
				UpotKorrLJ *= 2.0*M_PI;
			}
		}
	}
	
	// Berechnung der Korrekturen in jedem Zeitschritt?
	double rlow, rlowInv, rlowInv2, rdashInv, UCorrTemp, rdash, rdashInv2, FCorrTemp, PNCorrTemp, PTCorrTemp;
	double UCorrSum = 0.0;

	for (auto tempMol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid(); ++tempMol) {
		unsigned long molID = tempMol->getID();
		double UCorrShells = 0.0;
		double FCorrShells = 0.0;
		double PNCorrShells = 0.0;
		double PTCorrShells = 0.0;
		double ksi2 = ksi[molID]*ksi[molID];
		for (unsigned ci = 0; ci < numComp; ++ci){		
			for (unsigned cj = 0; cj < numComp; ++cj){	
				ParaStrm& params = _domain->getComp2Params()(ci,cj);
				params.reset_read();
				for (unsigned si = 0; si < numLJ[ci]; ++si) { // Long Range Correction for Lennard-Jones sites
					for (unsigned sj = 0; sj < numLJ[cj]; ++sj) {
						double eps24;
						double sigma;
						double sig2;
						double shift6;
						double eps;
						params >> eps24;
						params >> sig2;
						params >> shift6;
						sigma=sqrt(sig2);
						eps=eps24/24;
						//double tau1 = eLong[numLJSum2[ci]+si];  // ÜBERPRÜFEN
						//double tau2 = eLong[numLJSum2[cj]+sj];
						double tau1 = 0.0;
						double tau2 = 0.0;
						double sigma6 = sig2*sig2*sig2;
						double factorU = -M_PI*_drShells*eps24*sigma6/(6*ksi[molID]);
						double factorF = -factorU/ksi[molID];
						double factorP = 0.5*factorF/ksi[molID];
						if ( (tau1 == 0.0) && (tau2 == 0.0) ) { // Center-Center
							for (unsigned long j=1; j<(lowerS[molID]); j++) { // Loop over Shells with smaller Radius
								if (rhoShellsT[j] != 0.0) {
									rlow = ksi[molID] - RShells[j];
									rlowInv = 1/rlow;
									rlowInv2 = rlowInv * rlowInv;
									rdashInv = 1/(RShells[j] + ksi[molID]);
									UCorrTemp = sigma6*0.2*(pow(rdashInv,10) - pow(rlowInv,10) )
												- 0.5*(pow(rdashInv,4) - pow(rlowInv,4) );
									if ( rlow < rcmax) {
										rdash = min( rcmax,(RShells[j] + ksi[molID]) );
										rdashInv = 1/rdash;
										rdashInv2 = rdashInv * rdashInv;
										FCorrTemp = sigma6
												* ( (6/5*rdash*rdash + ksi2 - RShells2[j]) * pow(rdashInv2,6)
													- (6/5*rlow*rlow + ksi2 - RShells2[j]) * pow(rlowInv2,6) )
													- (1.5*rdash*rdash + ksi2 - RShells2[j]) * pow(rdashInv2,3)
													+ (1.5*rlow*rlow + ksi2 - RShells2[j]) * pow(rlowInv2,3);
										FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
										PNCorrTemp = 1.5*sigma6 * (pow(rdashInv2,4) - pow(rlowInv2,4))
											- 3 * (rdashInv2 - rlowInv2)
											+ 2*(ksi2-RShells2[j]) * ( 6/5*sigma6 * (pow(rdashInv2,5) - pow(rlowInv2,5))
											- 1.5 * (pow(rdashInv2,2) - pow(rlowInv2,2)) )
											+ pow((ksi2-RShells2[j]),2) * ( sigma6 * ( pow(rdashInv2,6) - pow(rlowInv2,6)) 
											- (pow(rdashInv2,3) - pow(rlowInv2,3)) );
										PNCorrShells = PNCorrShells + PNCorrTemp*factorP*rhoShellsT[j]*RShells[j];
										PTCorrTemp = 6/5*sigma6 * (pow(rdashInv2,5) - pow(rlowInv2,5))
												- 1.5 * (pow(rdashInv2,2) - pow(rlowInv2,2)) ;
										PTCorrShells = PTCorrShells + 2*ksi2*PTCorrTemp*factorP*rhoShellsT[j]*RShells[j];
									}
									UCorrShells = UCorrShells + UCorrTemp*factorU*rhoShellsT[j]*RShells[j];
								}
							}
							for (unsigned long j=upperS[molID]; j<NShells; j++) {
								if (rhoShellsT[j] != 0.0) {
									rlow = RShells[j] - ksi[molID];
									rlowInv = 1/rlow;
									rlowInv2 = rlowInv * rlowInv;
									rdashInv = 1/(RShells[j] + ksi[molID]);
									UCorrTemp = sigma6*0.2*(pow(rdashInv,10) - pow(rlowInv,10) )
												- 0.5*(pow(rdashInv,4) - pow(rlowInv,4) );
									if ( rlow < rcmax) {
										rdash = min( rcmax,(RShells[j] + ksi[molID]) );
										rdashInv = 1/rdash;
										rdashInv2 = rdashInv * rdashInv;
										FCorrTemp = sigma * ( (6/5*rdash*rdash + ksi2 - RShells2[j]) * pow(rdashInv2,6)
													- (6/5*rlow*rlow + ksi2 - RShells2[j]) * pow(rlowInv2,6) )
													- (1.5*rdash*rdash + ksi2 - RShells2[j]) * pow(rdashInv2,3)
													+ (1.5*rlow*rlow + ksi2 - RShells2[j]) * pow(rlowInv2,3) ;
										FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
										PNCorrTemp = 1.5*sigma6 * (pow(rdashInv2,4) - pow(rlowInv2,4))
													- 3 * (rdashInv2 - rlowInv2)
													+ 2*(ksi2-RShells2[j]) * ( 6/5*sigma6 * (pow(rdashInv2,5) - pow(rlowInv2,5))
													- 1.5 * (pow(rdashInv2,2) - pow(rlowInv2,2)) ) 
													+ pow((ksi2-RShells2[j]),2) * ( sigma6 * ( pow(rdashInv2,6) - pow(rlowInv2,6))
													- (pow(rdashInv2,3) - pow(rlowInv2,3)) );
										PNCorrShells = PNCorrShells + PNCorrTemp*factorP*rhoShellsT[j]*RShells[j];
										PTCorrTemp = 6/5*sigma6 * (pow(rdashInv2,5) - pow(rlowInv2,5))
													- 1.5 * (pow(rdashInv2,2) - pow(rlowInv2,2));
										PTCorrShells = PTCorrShells + 2*ksi2*PTCorrTemp*factorP*rhoShellsT[j]*RShells[j];
									}
									UCorrShells = UCorrShells + UCorrTemp*factorU*rhoShellsT[j]*RShells[j];
								}
							}
							for (unsigned long j=interS[molID]; j<upperS[molID]; j++) { // Loop over partly contributing Shells, Indizes prüfen!
								if (rhoShellsT[j] != 0.0) {
									rlow = rc;
									rlowInv = 1/rlow;
									rlowInv2 = rlowInv * rlowInv;
									rdashInv = 1/(RShells[j] + ksi[molID]);
									UCorrTemp = sigma6*0.2*(pow(rdashInv,10) - pow(rlowInv,10) )
												- 0.5*(pow(rdashInv,4) - pow(rlowInv,4) );
									rdash = min( rcmax,(RShells[j] + ksi[molID]) );
									rdashInv = 1/rdash;
									rdashInv2 = rdashInv * rdashInv;
									FCorrTemp = sigma6 * ( (6/5*rdash*rdash + ksi2 - RShells2[j]) * pow(rdashInv2,6)
												- (6/5*rlow*rlow + ksi2 - RShells2[j]) * pow(rlowInv2,6) )
												- (1.5*rdash*rdash + ksi2 - RShells2[j]) * pow(rdashInv2,3)
												+ (1.5*rlow*rlow + ksi2 - RShells2[j]) * pow(rlowInv2,3) ;
									FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
									PNCorrTemp = 1.5*sigma6 * (pow(rdashInv2,4) - pow(rlowInv2,4))
												- 3 * (rdashInv2 - rlowInv2) 
												+ 2*(ksi2-RShells2[j]) * ( 6/5*sigma6 * (pow(rdashInv2,5) - pow(rlowInv2,5))
												- 1.5 * (pow(rdashInv2,2) - pow(rlowInv2,2)) ) 
												+ pow((ksi2-RShells2[j]),2) * ( sigma6 * ( pow(rdashInv2,6) - pow(rlowInv2,6)) 
												- (pow(rdashInv2,3) - pow(rlowInv2,3)) );
									PNCorrShells = PNCorrShells + PNCorrTemp*factorP*rhoShellsT[j]*RShells[j];
									PTCorrTemp = 6/5*sigma6 * (pow(rdashInv2,5) - pow(rlowInv2,5))
												- 1.5 * (pow(rdashInv2,2) - pow(rlowInv2,2)) ;
									PTCorrShells = PTCorrShells + 2*ksi2*PTCorrTemp*factorP*rhoShellsT[j]*RShells[j];
									UCorrShells = UCorrShells + UCorrTemp*factorU*rhoShellsT[j]*RShells[j];
								}
							}
						}
						if ( (tau1 == 0.0) || (tau2 == 0.0) ) { // Center-Site
						double tau = max(tau1,tau2);   // Filtert den Wert Null raus
						for (unsigned long j=1; j<(lowerS[molID]); j++) {
							if (rhoShellsT[j] != 0.0) {
							rlow = ksi[molID] - RShells[j];
							UCorrTemp = SICSu(-6,RShells[j] + ksi[molID],tau) * sigma6
										- SICSu(-6,rlow,tau) * sigma6
										- SICSu(-3,RShells[j] + ksi[molID],tau)
										+ SICSu(-3,rlow,tau);
							if ( rlow < rcmax) {
								rdash = min( rcmax,(RShells[j] + ksi[molID]) );
								FCorrTemp = SICSu(-6,rdash,tau) * sigma6
										- SICSu(-6,rlow,tau) * sigma6
										- SICSu(-3,rdash,tau)
										+ SICSu(-3,rlow,tau);
								FCorrShells = FCorrShells - 2*FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
								FCorrTemp = (rdash*rdash+ksi[molID]*ksi[molID]-RShells2[j])/rdash
										* ( -11*sigma6*SICSf(-13,rdash,tau)
											+ 5*SICSf(-7,rdash,tau) ) 
										- (rlow*rlow+ksi[molID]*ksi[molID]-RShells2[j])/rlow
										* ( -11*sigma6*SICSf(-13,rlow,tau)
											+ 5*SICSf(-7,rlow,tau) ) ;
								FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT[j]*RShells[j];

							}
							UCorrShells = UCorrShells - 2*UCorrTemp*factorU*rhoShellsT[j]*RShells[j];
							}
						}
						for (unsigned long j=upperS[molID]; j<NShells; j++) {
							if (rhoShellsT[j] != 0.0) {
							rlow = RShells[j] - ksi[molID];
							UCorrTemp = SICSu(-6,RShells[j] + ksi[molID],tau) * sigma6
										- SICSu(-6,rlow,tau) * sigma6
										- SICSu(-3,RShells[j] + ksi[molID],tau)
										+ SICSu(-3,rlow,tau);
							if ( rlow < rcmax) {
								rdash = min( rcmax,(RShells[j] + ksi[molID]) );
								FCorrTemp = SICSu(-6,rdash,tau) * sigma6
										- SICSu(-6,rlow,tau) * sigma6
										- SICSu(-3,rdash,tau)
										+ SICSu(-3,rlow,tau);
								FCorrShells = FCorrShells - 2*FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
								FCorrTemp = (rdash*rdash+ksi[molID]*ksi[molID]-RShells2[j])/rdash
										* ( -11*sigma6*SICSf(-13,rdash,tau)
											+ 5*SICSf(-7,rdash,tau) ) 
										- (rlow*rlow+ksi[molID]*ksi[molID]-RShells2[j])/rlow
										* ( -11*sigma6*SICSf(-13,rlow,tau)
											+ 5*SICSf(-7,rlow,tau) ) ;
								FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
							}
							UCorrShells = UCorrShells - 2*UCorrTemp*factorU*rhoShellsT[j]*RShells[j];
							}
						}
						for (unsigned long j=interS[molID]; j<upperS[molID]; j++) {
							if (rhoShellsT[j] != 0.0) {
							rlow = rc;
							UCorrTemp = SICSu(-6,RShells[j] + ksi[molID],tau) * sigma6
										- SICSu(-6,rlow,tau) * sigma6
										- SICSu(-3,RShells[j] + ksi[molID],tau)
										+ SICSu(-3,rlow,tau);
							rdash = min( rcmax,(RShells[j] + ksi[molID]) );
							FCorrTemp = SICSu(-6,rdash,tau) * sigma6
										- SICSu(-6,rlow,tau) * sigma6
										- SICSu(-3,rdash,tau)
										+ SICSu(-3,rlow,tau);
							FCorrShells = FCorrShells - 2*FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
							FCorrTemp = (rdash*rdash+ksi[molID]*ksi[molID]-RShells2[j])/rdash
										* ( -11*sigma6*SICSf(-13,rdash,tau)
											+ 5*SICSf(-7,rdash,tau) ) 
										- (rlow*rlow+ksi[molID]*ksi[molID]-RShells2[j])/rlow
										* ( -11*sigma6*SICSf(-13,rlow,tau)
											+ 5*SICSf(-7,rlow,tau) ) ;
							FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
							UCorrShells = UCorrShells - 2*UCorrTemp*factorU*rhoShellsT[j]*RShells[j];
							}
						}
						}
						if ( (tau1 == 0.0) && (tau2 == 0.0) ) { // Site-Site
						for (unsigned long j=1; j<(lowerS[molID]); j++) {
							if (rhoShellsT[j] != 0.0) {
							rlow = ksi[molID] - RShells[j];
							UCorrTemp = SISSu(-6,RShells[j] + ksi[molID],tau1,tau2) * sigma6
										- SISSu(-6,rlow,tau1,tau2) * sigma6
										- SISSu(-3,RShells[j] + ksi[molID],tau1,tau2)
										+ SISSu(-3,rlow,tau1,tau2);
							if ( rlow < rcmax) {
								rdash = min( rcmax,(RShells[j] + ksi[molID]) );
								FCorrTemp = SISSu(-6,rdash,tau1,tau2) * sigma6
										- SISSu(-6,rlow,tau1,tau2) * sigma6
										- SISSu(-3,rdash,tau1,tau2)
										+ SISSu(-3,rlow,tau1,tau2);
								FCorrShells = FCorrShells + 2*FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
								FCorrTemp = (rdash*rdash+ksi[molID]*ksi[molID]-RShells2[j])/rdash
										* ( -11*sigma6*SISSf(-13,rdash,tau1,tau2)
											+ 5*SISSf(-7,rdash,tau1,tau2) ) 
										- (rlow*rlow+ksi[molID]*ksi[molID]-RShells2[j])/rlow
										* ( -11*sigma6*SISSf(-13,rlow,tau1,tau2)
											+ 5*SISSf(-7,rlow,tau1,tau2) ) ;
								FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
							}
							UCorrShells = UCorrShells - 2*UCorrTemp*factorU*rhoShellsT[j]*RShells[j];
							}
						}
						for (unsigned long j=upperS[molID]; j<NShells; j++) {
							if (rhoShellsT[j] != 0.0) {
							rlow = RShells[j] - ksi[molID];
							UCorrTemp = SISSu(-6,RShells[j] + ksi[molID],tau1,tau2) * sigma6
										- SISSu(-6,rlow,tau1,tau2) * sigma6
										- SISSu(-3,RShells[j] + ksi[molID],tau1,tau2)
										+ SISSu(-3,rlow,tau1,tau2);
							if ( rlow < rcmax) {
								rdash = min( rcmax,(RShells[j] + ksi[molID]) );
								FCorrTemp = SISSu(-6,rdash,tau1,tau2) * sigma6
										- SISSu(-6,rlow,tau1,tau2) * sigma6
										- SISSu(-3,rdash,tau1,tau2)
										+ SISSu(-3,rlow,tau1,tau2);
								FCorrShells = FCorrShells + 2*FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
								FCorrTemp = (rdash*rdash+ksi[molID]*ksi[molID]-RShells2[j])/rdash
										* ( -11*sigma6*SISSf(-13,rdash,tau1,tau2)
											+ 5*SISSf(-7,rdash,tau1,tau2) ) 
										- (rlow*rlow+ksi[molID]*ksi[molID]-RShells2[j])/rlow
										* ( -11*sigma6*SISSf(-13,rlow,tau1,tau2)
											+ 5*SISSf(-7,rlow,tau1,tau2) ) ;
								FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
							}
							UCorrShells = UCorrShells - 2*UCorrTemp*factorU*rhoShellsT[j]*RShells[j];
							}
						}
						for (unsigned long j=interS[molID]; j<upperS[molID]; j++) {
							if (rhoShellsT[j] != 0.0) {
							rlow = rc;
							UCorrTemp = SISSu(-6,RShells[j] + ksi[molID],tau1,tau2) * sigma6
										- SISSu(-6,rlow,tau1,tau2) * sigma6
										- SISSu(-3,RShells[j] + ksi[molID],tau1,tau2)
										+ SISSu(-3,rlow,tau1,tau2);
							rdash = min( rcmax,(RShells[j] + ksi[molID]) );
							FCorrTemp = SISSu(-6,rdash,tau1,tau2) * sigma6
										- SISSu(-6,rlow,tau1,tau2) * sigma6
										- SISSu(-3,rdash,tau1,tau2)
										+ SISSu(-3,rlow,tau1,tau2);
							FCorrShells = FCorrShells + 2*FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
							FCorrTemp = (rdash*rdash+ksi[molID]*ksi[molID]-RShells2[j])/rdash
										* ( -11*sigma6*SISSf(-13,rdash,tau1,tau2)
											+ 5*SISSf(-7,rdash,tau1,tau2) ) 
										- (rlow*rlow+ksi[molID]*ksi[molID]-RShells2[j])/rlow
										* ( -11*sigma6*SISSf(-13,rlow,tau1,tau2)
											+ 5*SISSf(-7,rlow,tau1,tau2) ) ;
							FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT[j]*RShells[j];
							UCorrShells = UCorrShells - 2*UCorrTemp*factorU*rhoShellsT[j]*RShells[j];
							}
						}
						}
					}
				}
				}
				UCorrSum = UCorrSum + UCorrShells;
				for(unsigned short d=0;d<3;++d) {
					FcorrX[molID] = FcorrX[molID]/ksi[molID]*FCorrShells;
					FcorrY[molID] = FcorrY[molID]/ksi[molID]*FCorrShells;
					FcorrZ[molID] = FcorrZ[molID]/ksi[molID]*FCorrShells;
				}
				PTCorrShells -= PNCorrShells;  // Was passiert hiermit? Wird bei nächstem Mol. wieder auf 0 gesetzt und ist kein Array.
				if (molID == 1888) {
					std::cout << "Local rank " << _domainDecomposition->getRank() << " ---->" << PNCorrShells << " " << PTCorrShells << std::endl;
				}

				// OPEN(29,FILE=TRIM(Dateiname)//"_U.csv",POSITION="APPEND",STATUS="OLD")
				// WRITE(29,"(F12.8,";",F12.8)") ksi[molID], UCorrShells
				// CLOSE(29)
				// OPEN(49,FILE=TRIM(Dateiname)//"_PN.csv",POSITION="APPEND",STATUS="OLD")
				// WRITE(49,"(F12.8,";",F12.8)") ksi[molID], PNCorrShells
				// CLOSE(49)
				// OPEN(59,FILE=TRIM(Dateiname)//"_PT.csv",POSITION="APPEND",STATUS="OLD")
				// WRITE(59,"(F12.8,";",F12.8)") ksi[molID], PTCorrShells
				// CLOSE(59)


		
		}
	}

	double PTCorrShells=1.0; // DUMMY

	// Distribution of the Force, Energy and Virial to every Node
	_domainDecomposition->collCommInit(3*globalNumMols);
	for (unsigned i=0; i<globalNumMols; i++){
		_domainDecomposition->collCommAppendDouble(FcorrX[i]);
		_domainDecomposition->collCommAppendDouble(FcorrY[i]);
		_domainDecomposition->collCommAppendDouble(FcorrZ[i]);
	}
	_domainDecomposition->collCommAllreduceSum();
	for (unsigned i=0; i<globalNumMols; i++){
		FcorrX_global[i]=_domainDecomposition->collCommGetDouble();
		FcorrY_global[i]=_domainDecomposition->collCommGetDouble();
		FcorrZ_global[i]=_domainDecomposition->collCommGetDouble();
	}
	_domainDecomposition->collCommFinalize();

	_domainDecomposition->collCommInit(2);
	_domainDecomposition->collCommAppendDouble(UCorrSum);
	_domainDecomposition->collCommAppendDouble(PTCorrShells);
	_domainDecomposition->collCommAllreduceSum();
	double UCorrSum_global = _domainDecomposition->collCommGetDouble();
	double PTCorrShells_global = _domainDecomposition->collCommGetDouble();
	_domainDecomposition->collCommFinalize();

	for (auto tempMol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid(); ++tempMol) {
		unsigned long molID = tempMol->getID();
		double Fa[3]={0.0, 0.0, 0.0};
		Fa[0] = FcorrX_global[molID];
		Fa[1] = FcorrY_global[molID];
		Fa[2] = FcorrZ_global[molID];
		//tempMol->Fljcenteradd(i, Fa);  // Ich weiß nicht, was i ist...
		////tempMol->Viadd(??);
		////tempMol->Uadd(??);
	}

	// Setting the Energy and Virial correction
	//_domain->setUpotCorr(UCorrSum_global);
	//_domain->setVirialCorr(PTCorrShells_global); //????

}


void Spherical::centerCenter(double sig, double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj) {
	double a = 1;
}

void Spherical::centerSite(double sig, double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj){
	// double sig2=sig*sig;
	// double sig3=sig2*sig;
	// double t = eLong[numLJSum2[ci]+si] + eLong[numLJSum2[cj]+sj]; // one of them is equal to zero.
	global_log->error() << "LongRangeCorrection: Center-Site correction not yet supported. Program exit ..." << endl;
    Simulation::exit(-1);
}

void Spherical::siteSite(double sig, double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj){
	// double sig2=sig*sig;
	// double sig3=sig2*sig;
	// double sig4=sig2*sig2;
	// double t1 = eLong[numLJSum2[ci]+si];
	// double t2 = eLong[numLJSum2[cj]+sj];
	// double tP = t1 + t2; // tau+ 
	// double tM = t1 - t2; // tau-

	global_log->error() << "LongRangeCorrection: Site-Site correction not yet supported. Program exit ..." << endl;
    Simulation::exit(-1);
}

void Spherical::writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstepI)
{
	double a = 1;
}



double Spherical::RhoP(double r, double rhov, double rhol, double D0, double R0) {

  return 0.5*(rhol+rhov) - 0.5*(rhol-rhov)*tanh(2*(r-R0)/D0);

}


double Spherical::SICCu(int n, double r) {

  return -( pow(r,(2*n + 2)) ) /( (n+1) );

}

double Spherical::SICSu(int n, double r, double tau) {

  return ( pow((r+tau),(2*n+3)) - pow((r-tau),(2*n+3)) ) / ( 4 * tau * (n+1) * (2*n+3) );

}

double Spherical::SISSu(int n, double r, double tau1, double tau2) {

  double tauPlus = tau1+tau2;
  double tauMinus = tau1-tau2;

  return (pow((r+tauPlus),(2*n+4)) - pow((r+tauMinus),(2*n+4))
  			- pow((r-tauMinus),(2*n+4)) + pow((r-tauPlus),(2*n+4)) )
			/ ( 8*tau1*tau2*(n+1)*(2*n+3)*(2*n+4) );
}


double Spherical::SICSf(int n, double r, double tau) {

  return ( pow((r+tau),(n+3)) - pow((r-tau),(n+3)) )
          / ( 4 * tau * (0.5*n+1) * (n+3) );

}

double Spherical::SISSf(int n, double r, double tau1, double tau2) {

  double tauPlus = tau1+tau2;
  double tauMinus = tau1-tau2;

  return ( pow((r+tauPlus),(n+4)) - pow((r+tauMinus),(n+4))
             - pow((r-tauMinus),(n+4)) + pow((r-tauPlus),(n+4)) ) 
            /( 8*tau1*tau2*(0.5*n+1)*(n+3)*(n+4) );

}

double Spherical::TICCu(int n, double rcutoff, double sigma2) {

  return -( pow(rcutoff,(2*n+3)) )
            /( pow(sigma2,(n)*(2*n+3)) );
}


double Spherical::TICSu(int n, double rcutoff, double sigma2, double tau) {

  return -( pow((rcutoff+tau),(2*n+3)) - pow((rcutoff-tau),(2*n+3)) ) * rcutoff  
            /( pow(4*sigma2,(n))*tau*(n+1)*(2*n+3) )
          +  ( pow((rcutoff+tau),(2*n+4)) - pow((rcutoff-tau),(2*n+4)) )
            /( 4*pow(sigma2,n)*tau*(n+1)*(2*n+3)*(2*n+4) );
}


double Spherical::TISSu(int n, double rcutoff, double sigma2, double tau1, double tau2) {

  double tauPlus = tau1+tau2;
  double tauMinus = tau1-tau2;
  return -(   pow((rcutoff+tauPlus),(2*n+4)) - pow((rcutoff+tauMinus),(2*n+4))
               - pow((rcutoff-tauMinus),(2*n+4)) + pow((rcutoff-tauPlus),(2*n+4)) ) * rcutoff
            /( 8*pow(sigma2,(n))*tau1*tau2*(n+1)*(2*n+3)*(2*n+4) )
          +  (   pow((rcutoff+tauPlus),(2*n+5)) - pow((rcutoff+tauMinus),(2*n+5))
               - pow((rcutoff-tauMinus),(2*n+5)) + pow((rcutoff-tauPlus),(2*n+5)) )
            /( 8*pow(sigma2,(n))*tau1*tau2*(n+1)*(2*n+3)*(2*n+4)*(2*n+5) );
}


double Spherical::TICCp(int n, double rcutoff, double sigma2) {

  return 2*n * TICCu(n,rcutoff,sigma2);
}


double Spherical::TICSp(int n, double rcutoff, double sigma2, double tau) {

  return -( pow((rcutoff+tau),(2*n+2)) - pow((rcutoff-tau),(2*n+2))) * pow(rcutoff,(2))
            /( 4*pow(sigma2,(n))*tau*(n+1) )
          - 3*TICSu(n,rcutoff,sigma2,tau);
}


double Spherical::TISSp(int n, double rcutoff, double sigma2, double tau1, double tau2) {

  double tauPlus = tau1+tau2;
  double tauMinus = tau1-tau2;
  return -(   pow((rcutoff+tauPlus),(2*n+3)) - pow((rcutoff+tauMinus),(2*n+3))
               - pow((rcutoff-tauMinus),(2*n+3)) + pow((rcutoff-tauPlus),(2*n+3)) ) * pow(rcutoff,2)
            /( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3) )
          - 3*TISSu(n,rcutoff,sigma2,tau1,tau2);
}

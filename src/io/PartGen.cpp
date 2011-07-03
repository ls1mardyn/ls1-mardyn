#include "PartGen.h"

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "ensemble/GrandCanonical.h"

#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;

// initialize constants
const double PartGen::_atomicMassDim(1.660539E-27);
const double PartGen::_avogadroDim(6.02214E+23);
const double PartGen::_boltzmannDim(1.38065E-23);
const double PartGen::_cocoDim(2.307076e-28);
const double PartGen::_debyeDim(2.0819435e-11);

PartGen::PartGen(){
}

//! @brief set the phase space file name
void PartGen::setPhaseSpaceFile(string filename){
	
}

//! @brief set the phase space header file name (can be identical to the
//         phase space file
void PartGen::setPhaseSpaceHeaderFile(string filename){
	srand(time(0));

	ifstream inpfStream;
	inpfStream.open(filename.c_str());
	
	ignoreLines(inpfStream,2);
	
	// Read in reference values
	readRefValues(inpfStream);

	// Read in state point
	readStatePoint(inpfStream );

	// Algorithm
	readAlgorithm(inpfStream, _simBoxRatio);

	// calculate Size of the Simulation Box
	setVectorValues(_simBoxLength, 3, 1.0);
	_simBoxLength[0] = pow( (_numberOfMolecules * pow(_simBoxRatio[0],2)) /
	                        (_rho*_simBoxRatio[1]*_simBoxRatio[2]), (1.0/3.0));
	_simBoxLength[1] = pow( (_numberOfMolecules * pow(_simBoxRatio[1],2)) /
	                        (_rho*_simBoxRatio[0]*_simBoxRatio[2]), (1.0/3.0));
	_simBoxLength[2] = pow( (_numberOfMolecules * pow(_simBoxRatio[2],2)) /
	                        (_rho*_simBoxRatio[0]*_simBoxRatio[1]), (1.0/3.0));

	// Visualisierung
	// NO VISUALISATION PROCESSED
	ignoreLines(inpfStream, 7);

	// Molekuelbeschreibung
	int MaxAnzSites = getIntParamValue(inpfStream);
	int MaxAnzCharges = getIntParamValue(inpfStream);
	int MaxAnzDipole = getIntParamValue(inpfStream);
	int MaxAnzQuadrupole = getIntParamValue(inpfStream);
	int MaxAnzTersoff = getIntParamValue(inpfStream);

	ignoreLines(inpfStream,1);
	_epsilonRF = getDoubleParamValue(inpfStream);
	ignoreLines(inpfStream,3);
	

	// Mixing rules
	setMatrixValues(_etaLB, _numberOfComponents, _numberOfComponents, 1.0);
	setMatrixValues(_xiLB, _numberOfComponents, _numberOfComponents, 1.0);
	for(int j=0; j<_numberOfComponents-1; j++){
		for(int i=j+1; i<_numberOfComponents; i++){
			getDoubleParamValues(inpfStream, _etaLB[i][j], _xiLB[i][j]);
			_etaLB[j][i] = _etaLB[i][j]; // symmetrisch
			_xiLB[j][i] = _xiLB[i][j]; // symmetrisch
		}
	}

	// 
	setVectorValues(_numSites, _numberOfComponents, 0);
	setVectorValues(_numDipoles, _numberOfComponents, 0);
	setVectorValues(_numQuadrupoles, _numberOfComponents, 0);
	setMatrixValues(_iBodyDummy, _numberOfComponents, 3, 0.0);
	set3DMatrixValues(_rSiteBody, _numberOfComponents, MaxAnzSites, 3, 0.0);
	set3DMatrixValues(_rDipoleBody, _numberOfComponents, MaxAnzDipole, 3, 0.0);
	set3DMatrixValues(_rQuadrupoleBody, _numberOfComponents, MaxAnzQuadrupole, 3, 0.0);
	setMatrixValues(_epsilonSite, _numberOfComponents, MaxAnzSites, 0.0);
	setMatrixValues(_mSite, _numberOfComponents, MaxAnzSites, 0.0);
	setMatrixValues(_sigmaSite, _numberOfComponents, MaxAnzSites, 0.0);
	set3DMatrixValues(_eMyBody, _numberOfComponents, MaxAnzDipole, 3, 0.0);
	setMatrixValues(_absMy, _numberOfComponents, MaxAnzDipole, 0.0);
	set3DMatrixValues(_eQBody, _numberOfComponents, MaxAnzQuadrupole, 3, 0.0);
	setMatrixValues(_absQ, _numberOfComponents, MaxAnzDipole, 0.0);
	
	setVectorValues(_numCharges, _numberOfComponents, 0);
	setVectorValues(_numTersoff, _numberOfComponents, 0);
	set3DMatrixValues(_rChargeBody, _numberOfComponents, MaxAnzCharges, 3, 0.0);
	set3DMatrixValues(_rTersoffBody, _numberOfComponents, MaxAnzTersoff, 3, 0.0);
	setMatrixValues(_mCharge, _numberOfComponents, MaxAnzCharges, 0);
	setMatrixValues(_qCharge, _numberOfComponents, MaxAnzCharges, 0);
	setMatrixValues(_mTersoff, _numberOfComponents, MaxAnzTersoff, 0);
	setMatrixValues(_ATersoff, _numberOfComponents, MaxAnzTersoff, 0);
	setMatrixValues(_BTersoff, _numberOfComponents, MaxAnzTersoff, 0);
	setMatrixValues(_lambdaTersoff, _numberOfComponents, MaxAnzTersoff, 0);
	setMatrixValues(_muTersoff, _numberOfComponents, MaxAnzTersoff, 0);
	setMatrixValues(_RTersoff, _numberOfComponents, MaxAnzTersoff, 0);
	setMatrixValues(_STersoff, _numberOfComponents, MaxAnzTersoff, 0);
	setMatrixValues(_cTersoff, _numberOfComponents, MaxAnzTersoff, 0);
	setMatrixValues(_dTersoff, _numberOfComponents, MaxAnzTersoff, 0);
	setMatrixValues(_hTersoff, _numberOfComponents, MaxAnzTersoff, 0);
	setMatrixValues(_nTersoff, _numberOfComponents, MaxAnzTersoff, 0);
	setMatrixValues(_betaTersoff, _numberOfComponents, MaxAnzTersoff, 0);

	// Loop over all Components
	for(int comp=0; comp<_numberOfComponents; comp++)
	{
		 ignoreLines(inpfStream,3);
		_numSites[comp] = getIntParamValue(inpfStream);
		_numCharges[comp] = getIntParamValue(inpfStream);
		_numDipoles[comp] = getIntParamValue(inpfStream);
		_numQuadrupoles[comp] = getIntParamValue(inpfStream);
		_numTersoff[comp] = getIntParamValue(inpfStream);
		if( (_numSites[comp] > MaxAnzSites)    || (_numCharges[comp] > MaxAnzCharges) ||
		    (_numDipoles[comp] > MaxAnzDipole) || (_numQuadrupoles[comp] > MaxAnzQuadrupole)
	                                   || (_numTersoff[comp] > MaxAnzTersoff) )
		{
			cout << "Exceeded maximal number of interaction sites." << endl;
			exit(1);
		}
		ignoreLines(inpfStream,1);
		for(int site=0; site<_numSites[comp]; site++)
		{
			getDoubleParamValues(inpfStream, _rSiteBody[comp][site][0], 
			                                 _rSiteBody[comp][site][1], 
			                                 _rSiteBody[comp][site][2]);
		}
		for(int cg=0; cg < this->_numCharges[comp]; cg++)
		{
			 getDoubleParamValues(inpfStream, _rChargeBody[comp][cg][0], 
			                                  _rChargeBody[comp][cg][1], 
			                                  _rChargeBody[comp][cg][2]);
		}
		for(int dip=0; dip<_numDipoles[comp]; dip++)
		{
			getDoubleParamValues(inpfStream, _rDipoleBody[comp][dip][0], 
			                                 _rDipoleBody[comp][dip][1], 
			                                 _rDipoleBody[comp][dip][2]);
		}
		for(int quad=0; quad<_numQuadrupoles[comp]; quad++)
		{
			getDoubleParamValues(inpfStream, _rQuadrupoleBody[comp][quad][0], 
			                                 _rQuadrupoleBody[comp][quad][1], 
			                                 _rQuadrupoleBody[comp][quad][2]);
		}
		for(int tter=0; tter < this->_numTersoff[comp]; tter++)
		{
			 getDoubleParamValues(inpfStream, _rChargeBody[comp][tter][0], 
			                                  _rChargeBody[comp][tter][1], 
			                                  _rChargeBody[comp][tter][2]);
		}

		for(int site=0; site<_numSites[comp]; site++)
		{
			getDoubleParamValues(inpfStream, _epsilonSite[comp][site], _mSite[comp][site], _sigmaSite[comp][site]);
		}
		for(int cg = 0; cg < this->_numCharges[comp]; cg++)
		{
			 getDoubleParamValues( inpfStream, _mCharge[comp][cg],
			                 _qCharge[comp][cg] );
		}
		for(int dip=0; dip<_numDipoles[comp]; dip++){
			getDoubleParamValues(inpfStream, _eMyBody[comp][dip][0], 
			                                 _eMyBody[comp][dip][1], 
			                                 _eMyBody[comp][dip][2], 
			                                 _absMy[comp][dip]);
			double dr2 = dotprod(_eMyBody[comp][dip], _eMyBody[comp][dip]);
			if(dr2 != 0.0){
				double dr = sqrt(dr2);
				vecmult(_eMyBody[comp][dip], (1/dr));
			}
		}
		for(int quad=0; quad<_numQuadrupoles[comp]; quad++){
			getDoubleParamValues(inpfStream, _eQBody[comp][quad][0], 
			                                 _eQBody[comp][quad][1], 
			                                 _eQBody[comp][quad][2], 
			                                 _absQ[comp][quad]);
			double dr2 = dotprod(_eQBody[comp][quad], _eQBody[comp][quad]);
			if(dr2 != 0.0){
				double dr = sqrt(dr2);
				vecmult(_eMyBody[comp][quad], (1/dr));
			}
		}    
		for(int tter = 0; tter < this->_numTersoff[comp]; tter++)
		{
			 getDoubleParamValues( inpfStream, _mTersoff[comp][tter],
			                 _ATersoff[comp][tter],
			                 _BTersoff[comp][tter],
			                 _lambdaTersoff[comp][tter] );
			 getDoubleParamValues( inpfStream, _muTersoff[comp][tter],
			                 _RTersoff[comp][tter],
			                 _STersoff[comp][tter],
			                 _cTersoff[comp][tter] );
			 getDoubleParamValues( inpfStream, _dTersoff[comp][tter],
			                 _hTersoff[comp][tter],
			                 _nTersoff[comp][tter],
			                 _betaTersoff[comp][tter] );
		}
		ignoreLines(inpfStream, 1);
		getDoubleParamValues( inpfStream, _iBodyDummy[comp][0], 
		                                  _iBodyDummy[comp][1], 
		                                  _iBodyDummy[comp][2] );
		for(int dim=0; dim<3; dim++){
			_iBodyDummy[comp][dim] *= _atomicMassDim*1.0E-20;
		}
	}
	// VISUALISIERUNG WIRD NICHT DURCHGEFUEHRT
	inpfStream.close();

	// reduce all Molecule Parameters
	//
	// welchen Sinn ergeben diese Koeffizienten (1e-10 etc.)?
	// vielleicht war die Idee diese: in readRefValues werden
	// eps_ref, m_ref und sig_ref kuenstlich in das Einheitensystem
	// J, kg, m uebertragen. Da aber die Potentialparameter alle in
	// K, u, A gegeben sind, muss diese Umwandlung jetzt bei jeder
	// Verwendung dieser Variablen nochmal rueckwaerts durchgefuehrt
	// werden.
	//
	matmult(_iBodyDummy, 1 / (_mRefDim*pow(_sigmaRefDim,2)));
	mat3dmult(_rSiteBody, 1.0E-10 / _sigmaRefDim);
	mat3dmult(_rChargeBody, 1.0E-10 / _sigmaRefDim);
	mat3dmult(_rDipoleBody, 1.0E-10 / _sigmaRefDim);
	mat3dmult(_rQuadrupoleBody, 1.0E-10 / _sigmaRefDim);
	mat3dmult(_rTersoffBody, 1.0E-10 / _sigmaRefDim);
	matmult(_epsilonSite, _boltzmannDim / _epsilonRefDim);
	matmult(_mSite, _atomicMassDim / _mRefDim);
	matmult(_sigmaSite, 1.0E-10 / _sigmaRefDim);
	/*
	 * These lines do not appear to make a lot of sense:
	 *
	matmult(_absMy, 1 / (1.0E+24 * sqrt(10.0*_epsilonRefDim*(pow(_sigmaRefDim,3.0))) ));
	matmult(_absQ, 1 / (1.0E+34 * sqrt(10.0*_epsilonRefDim*(pow(_sigmaRefDim,5.0))) ));
	 */

	/*
	 * This version assumes that parameters are given in e, D, and DA.
	 */
	double qRefDim = sqrt(_sigmaRefDim * _epsilonRefDim / _cocoDim);  // [e]
	double muRefDim = qRefDim * _sigmaRefDim / _debyeDim;  // [D]
	double qdrRefDim = muRefDim * (1.0e+10 * _sigmaRefDim);  // [DA]
	matmult(_absMy, 1.0 / muRefDim);
	matmult(_absQ, 1.0 / qdrRefDim);
	
	matmult(_mCharge, _atomicMassDim / _mRefDim);
	matmult(_qCharge, 1.0 / qRefDim);
	matmult(_mTersoff, _atomicMassDim / _mRefDim);
	matmult(_ATersoff, _boltzmannDim / _epsilonRefDim);  // given in units of K
	matmult(_BTersoff, _boltzmannDim / _epsilonRefDim);  // given in units of K
	matmult(_lambdaTersoff, 1.0e+10 * _sigmaRefDim);  // given in units of 1/A
	matmult(_muTersoff, 1.0e+10 * _sigmaRefDim);  // given in units of 1/A
	matmult(_RTersoff, 1.0E-10 / _sigmaRefDim);
	matmult(_STersoff, 1.0E-10 / _sigmaRefDim);
 
	// transform all components to a principle axis system
	principleAxisTransform();

	// calculate the total mass of each component
	setVectorValues(_massOfComps,_numberOfComponents,0.0);
	for(int comp=0; comp<_numberOfComponents; comp++){ 
		for(int site=0; site<_numSites[comp]; site++){
			_massOfComps[comp] += _mSite[comp][site];
		}
	}    
}
	
//! @brief read the phase space components and header information
//! the LJ cutoff radius must be passed using readAlgorithm
void PartGen::readPhaseSpaceHeader(
	 Domain* domain, double timestep
){
	vector<Component>& dcomponents = domain->getComponents();
	domain->setCurrentTime(0);
	domain->setGlobalTemperature(_temperature);
	domain->setGlobalLength(0, _simBoxLength[0]);
	domain->setGlobalLength(1, _simBoxLength[1]);
	domain->setGlobalLength(2, _simBoxLength[2]);
	dcomponents.resize(_numberOfComponents);

	for(int comp=0; comp<_numberOfComponents; comp++) {
		dcomponents[comp].setID(comp);
		unsigned int numljcenters   = getNumSites(comp);
		unsigned int numcharges     = getNumCharges(comp);
		unsigned int numdipoles     = getNumDipoles(comp);
		unsigned int numquadrupoles = getNumQuadrupoles(comp);
		unsigned int numtersoff     = getNumTersoff(comp);
		for(unsigned int site=0; site<numljcenters; site++)
		{
			dcomponents[comp].addLJcenter(getSitePos(comp, site, 0),
			                              getSitePos(comp, site, 1),
			                              getSitePos(comp, site, 2),
			                              getMSite(comp, site),
			                              getEpsilonSite(comp, site),
			                              getSigmaSite(comp, site),
				    this->_LJCutoffRadius,
				    getShiftSite(comp, site));
		}
		for(unsigned int cg = 0; cg < numcharges; cg++)
		{
			 dcomponents[comp].addCharge( getChargePos(comp, cg, 0),
			                              getChargePos(comp, cg, 1),
			                              getChargePos(comp, cg, 2),
				    getChargeMass(comp, cg),
				    getCharge(comp, cg) );
		}
		for(unsigned int dip=0; dip<numdipoles; dip++)
		{
			dcomponents[comp].addDipole(getDipolePos(comp, dip, 0),
			                            getDipolePos(comp, dip, 1),
			                            getDipolePos(comp, dip, 2),
			                            getEMyBody(comp, dip, 0),
			                            getEMyBody(comp, dip, 1),
			                            getEMyBody(comp, dip, 2),
			                            getAbsMy(comp, dip));
		}
		for(unsigned int quad=0; quad<numquadrupoles; quad++)
		{
			dcomponents[comp].addQuadrupole(getQuadrupolePos(comp, quad, 0),
			                                getQuadrupolePos(comp, quad, 1),
			                                getQuadrupolePos(comp, quad, 2),
			                                getEQBody(comp, quad, 0),
			                                getEQBody(comp, quad, 1),
			                                getEQBody(comp, quad, 2),
			                                getAbsQ(comp, quad));
		}
		for(unsigned int tter = 0; tter < numtersoff; tter++)
		{
			 dcomponents[comp].addTersoff( getTersoffPos(comp, tter, 0),
			                               getTersoffPos(comp, tter, 1),
			                               getTersoffPos(comp, tter, 2),
				     getTersoffMass(comp, tter),
				     getTersoffA(comp, tter),
				     getTersoffB(comp, tter),
				     getTersoff_lambda(comp, tter),
				     getTersoff_mu(comp, tter),
				     getTersoffR(comp, tter),
				     getTersoffS(comp, tter),
				     getTersoff_c(comp, tter),
				     getTersoff_d(comp, tter),
				     getTersoff_h(comp, tter),
				     getTersoff_n(comp, tter),
				     getTersoff_beta(comp, tter) );
		}
		double IDummy1,IDummy2,IDummy3;
		IDummy1 = getIDummy(comp, 0);
		IDummy2 = getIDummy(comp, 1);
		IDummy3 = getIDummy(comp, 2);
		if(IDummy1>0.) dcomponents[comp].setI11(IDummy1);
		if(IDummy2>0.) dcomponents[comp].setI22(IDummy2);
		if(IDummy3>0.) dcomponents[comp].setI33(IDummy3);
	}
	vector<double>& dmixcoeff = domain->getmixcoeff();
	dmixcoeff.clear();
	for(int comp1=0; comp1<_numberOfComponents-1; comp1++){
		for(int comp2=comp1+1; comp2<_numberOfComponents; comp2++){
			dmixcoeff.push_back(getXi(comp1, comp2));
			dmixcoeff.push_back(getEta(comp1, comp2));
		}
	}
	domain->setepsilonRF(_epsilonRF);
}
	
//! @brief read the actual phase space information
unsigned long PartGen::readPhaseSpace(
	 ParticleContainer* particleContainer,
	 list<ChemicalPotential>* lmu, Domain* domain,
	 DomainDecompBase* domainDecomp
) {
	vector<double> bBoxMin;
	vector<double> bBoxMax;
	//_globalLength[0] = partGen->getLSimBox(0);
	//_globalLength[1] = partGen->getLSimBox(1);
	//_globalLength[2] = partGen->getLSimBox(2);

	bBoxMin.resize(3);
	bBoxMax.resize(3);
	for (int i=0; i<3; i++) {
		bBoxMin[i] = domainDecomp->getBoundingBoxMin(i, domain);
		bBoxMax[i] = domainDecomp->getBoundingBoxMax(i, domain);
	}
	//partGen->createHomogeneousDist(particleContainer, bBoxMin, bBoxMax, this, domainDecomp);
	createClusters(particleContainer, bBoxMin, bBoxMax, domain, domainDecomp);

	vector<Component>& dcomponents = domain->getComponents();
	vector<unsigned long> partsPerComp;
	partsPerComp.resize(_numberOfComponents);
	domain->setglobalNumMolecules(domainDecomp->countMolecules(particleContainer, partsPerComp));
	for(unsigned int i=0; i<partsPerComp.size(); i++) {
		dcomponents[i].setNumMolecules(partsPerComp[i]);
		domain->setglobalRotDOF(partsPerComp[i]*dcomponents[i].getRotationalDegreesOfFreedom());
	}
	domain->setglobalRho(domain->getglobalNumMolecules()/(_simBoxLength[0]*_simBoxLength[1]*_simBoxLength[2]));
	
	unsigned long maxid = 0;
	for( Molecule* pp = particleContainer->begin();
	     pp != particleContainer->end();
	     pp = particleContainer->next() )
	{
	    if(pp->id() > maxid) maxid = pp->id();
	    std::list<ChemicalPotential>::iterator cpit;
	    bool brk = true;
	    for(cpit = lmu->begin(); cpit != lmu->end(); cpit++)
	    {
	       if(!cpit->hasSample())
	 {
	    if(pp->componentid() == cpit->getComponentID())
	    {
		       cpit->storeMolecule(*pp);
	    }
	    else brk = false;
	 }
			}
			if(brk) break;
	 }
	 return maxid;
}





void PartGen::setClusterFile(double gasDensity, double fluidDensity, double volPercOfFluid, string clusterFileName){
	_clusterFile = clusterFileName;
	_gasDensity = gasDensity;
	_fluidDensity = fluidDensity;

	double averageRho =   fluidDensity*volPercOfFluid/100.0 
		+ gasDensity*(1-volPercOfFluid/100.0);
	_simBoxLength[0] = pow( (_numberOfMolecules * pow(_simBoxRatio[0],2)) /
	      (averageRho*_simBoxRatio[1]*_simBoxRatio[2]), (1.0/3.0));
	_simBoxLength[1] = pow( (_numberOfMolecules * pow(_simBoxRatio[1],2)) /
	      (averageRho*_simBoxRatio[0]*_simBoxRatio[2]), (1.0/3.0));
	_simBoxLength[2] = pow( (_numberOfMolecules * pow(_simBoxRatio[2],2)) /
	      (averageRho*_simBoxRatio[0]*_simBoxRatio[1]), (1.0/3.0));

}

void PartGen::createHomogeneousDist(ParticleContainer* particleContainer, vector<double> &bBoxMin, vector<double> &bBoxMax, Domain* domain, DomainDecompBase* domainDecomp){


	vector<int> globalFccCells; 
	vector<int> localFccCellsMin; 
	vector<int> localFccCellsMax; 
	vector<double> fccCellLength; 
	
	setVectorValues(globalFccCells, 3, 0);
	setVectorValues(localFccCellsMin, 3, 0);
	setVectorValues(localFccCellsMax, 3, 0);
	setVectorValues(fccCellLength, 3, 0);

	for(int dim=0; dim<3; dim++){
		globalFccCells[dim] = (int) ceil( pow( _rho/4.0, 1.0/3.0 )*_simBoxLength[dim]);
		fccCellLength[dim] = _simBoxLength[dim]/globalFccCells[dim];
		localFccCellsMin[dim] = (int) floor(bBoxMin[dim]/fccCellLength[dim]);
		localFccCellsMax[dim] = (int) floor(bBoxMax[dim]/fccCellLength[dim]);
	}

	vector<vector<double> > fccOffsets;
	setMatrixValues(fccOffsets, 4, 3, 0.0);
	fccOffsets[1][0] = 0.5*fccCellLength[0];
	fccOffsets[1][1] = 0.5*fccCellLength[1];
	fccOffsets[2][0] = 0.5*fccCellLength[0];
	fccOffsets[2][2] = 0.5*fccCellLength[2];
	fccOffsets[3][1] = 0.5*fccCellLength[1];
	fccOffsets[3][2] = 0.5*fccCellLength[2];

	vector<double> r_;
	setVectorValues(r_, 3, 0.0);

	int molCount = 0;
	for(int fcc=0; fcc<4; fcc++){
		for(int iz=localFccCellsMin[2]; iz<=localFccCellsMax[2]; iz++){
			for(int iy=localFccCellsMin[1]; iy<=localFccCellsMax[1]; iy++){
				for(int ix=localFccCellsMin[0]; ix<=localFccCellsMax[0]; ix++){
					// Position
					r_[0] = ix*(fccCellLength[0]) + fccOffsets[fcc][0];
					r_[1] = iy*(fccCellLength[1]) + fccOffsets[fcc][1];
					r_[2] = iz*(fccCellLength[2]) + fccOffsets[fcc][2];

					if(not domainDecomp->procOwnsPos(r_[0], r_[1], r_[2], domain)){
						// if the position is not in the domain of this proc,
						// the molecule must not be created
						continue;
					}
					addParticle(molCount, r_[0], r_[1], r_[2], particleContainer, 
					            domain, domainDecomp);
					molCount++;
				}
			}
		}
	}
}

void PartGen::createClusters(ParticleContainer* particleContainer, vector<double> &bBoxMin, vector<double> &bBoxMax, Domain* domain, DomainDecompBase* domainDecomp){

	readLocalClusters(domain, domainDecomp);
	
	vector<int> globalFccCells; 
	vector<int> clusterFccCellsMin; 
	vector<int> clusterFccCellsMax; 
	vector<double> fccCellLength; 
	
	setVectorValues(globalFccCells, 3, 0);
	setVectorValues(clusterFccCellsMin, 3, 0);
	setVectorValues(clusterFccCellsMax, 3, 0);
	setVectorValues(fccCellLength, 3, 0);

	// fluid properties
	for(int dim=0; dim<3; dim++){
		globalFccCells[dim] = (int) ceil(pow( _fluidDensity/4., 1./3.)*_simBoxLength[dim]);
		fccCellLength[dim] = _simBoxLength[dim]/globalFccCells[dim];
	}
	vector<vector<double> > fccOffsets;
	setMatrixValues(fccOffsets, 4, 3, 0.0);
	fccOffsets[1][0] = 0.5*fccCellLength[0];
	fccOffsets[1][1] = 0.5*fccCellLength[1];
	fccOffsets[2][0] = 0.5*fccCellLength[0];
	fccOffsets[2][2] = 0.5*fccCellLength[2];
	fccOffsets[3][1] = 0.5*fccCellLength[1];
	fccOffsets[3][2] = 0.5*fccCellLength[2];

	//cout << "INITIAL DISTANCE BETWEEN PARTICLES: " << fccCellLength[0] << endl;

	double securityOffset = 0.5*fccCellLength[0];

	vector<double> clusterPos;
	clusterPos.resize(3);
	double radius;

	vector<double> r_;
	setVectorValues(r_, 3, 0.0);

	int molCount = 0;
	for(unsigned int cluster=0; cluster<_localClusters.size(); cluster++){
		radius = _localClusters[cluster][3];
		for(int dim=0; dim<3; dim++){
			clusterPos[dim] = _localClusters[cluster][dim];
			clusterFccCellsMin[dim] = (int) floor((clusterPos[dim]-radius)/fccCellLength[dim]);
			clusterFccCellsMax[dim] = (int) floor((clusterPos[dim]+radius)/fccCellLength[dim]);
		}
		
 
		for(int fcc=0; fcc<4; fcc++){
			for(int iz=clusterFccCellsMin[2]; iz<=clusterFccCellsMax[2]; iz++){
				for(int iy=clusterFccCellsMin[1]; iy<=clusterFccCellsMax[1]; iy++){
					for(int ix=clusterFccCellsMin[0]; ix<=clusterFccCellsMax[0]; ix++){
						// Position
						r_[0] = ix*(fccCellLength[0]) + fccOffsets[fcc][0];
						r_[1] = iy*(fccCellLength[1]) + fccOffsets[fcc][1];
						r_[2] = iz*(fccCellLength[2]) + fccOffsets[fcc][2];
		
		        if(sqrt(pow(r_[0] - clusterPos[0],2.0) +
		                pow(r_[1] - clusterPos[1],2.0) +
		                pow(r_[2] - clusterPos[2],2.0)) > radius){
		          // molecule is too far away from the center of the cluster
		          continue;
		        }
		        if(not domainDecomp->procOwnsPos(r_[0], r_[1], r_[2], domain)){
		          // if the position is not in the domain of this proc,
		          // the molecule must not be created
		          continue;
		        }
		        if(belongsToPreviousCluster(r_[0], r_[1], r_[2], cluster)){
		          // some other cluster already created this molecule
		          continue;
		        }

						addParticle(molCount, r_[0], r_[1], r_[2], particleContainer, 
						            domain, domainDecomp);
						molCount++;
					}
				}
			}
		}
	}

	// gas properties
	vector<int> fccCellsMin; 
	vector<int> fccCellsMax; 
	setVectorValues(fccCellsMin, 3, 0);
	setVectorValues(fccCellsMax, 3, 0);
	
	for(int dim=0; dim<3; dim++){
		globalFccCells[dim] = (int) ceil(pow( _gasDensity/4., 1./3.)*_simBoxLength[dim]);
		fccCellLength[dim] = _simBoxLength[dim]/globalFccCells[dim];
		fccCellsMin[dim] = (int) floor(bBoxMin[dim]/fccCellLength[dim]);
		fccCellsMax[dim] = (int) floor(bBoxMax[dim]/fccCellLength[dim]);
	}

	fccOffsets[1][0] = 0.5*fccCellLength[0];
	fccOffsets[1][1] = 0.5*fccCellLength[1];
	fccOffsets[2][0] = 0.5*fccCellLength[0];
	fccOffsets[2][2] = 0.5*fccCellLength[2];
	fccOffsets[3][1] = 0.5*fccCellLength[1];
	fccOffsets[3][2] = 0.5*fccCellLength[2];

	// GAS
	for(int fcc=0; fcc<4; fcc++){
		for(int iz=fccCellsMin[2]; iz<=fccCellsMax[2]; iz++){
			for(int iy=fccCellsMin[1]; iy<=fccCellsMax[1]; iy++){
				for(int ix=fccCellsMin[0]; ix<=fccCellsMax[0]; ix++){
					// Position
					r_[0] = ix*(fccCellLength[0]) + fccOffsets[fcc][0];
					r_[1] = iy*(fccCellLength[1]) + fccOffsets[fcc][1];
					r_[2] = iz*(fccCellLength[2]) + fccOffsets[fcc][2];
		
		      if(not domainDecomp->procOwnsPos(r_[0], r_[1], r_[2], domain)){
		        // if the position is not in the domain of this proc,
		        // the molecule must not be created
		        continue;
		      }
		      if(closeToAnyCluster(r_[0], r_[1], r_[2], securityOffset)){
		        // some other cluster already created this molecule
		        continue;
		      }
	
	        addParticle(molCount, r_[0], r_[1], r_[2], particleContainer, 
	        domain, domainDecomp);
	        molCount++;
	      }
	    }
		}
	}
}

void PartGen::printConfig(){
	cout << "##########################################" << endl;
	cout << "### Result of PartGen Const. ###" << endl;
	cout << "##########################################" << endl;

	cout << "_epsilonRefDim     : " << _epsilonRefDim << endl;
	cout << "_mRefDim           : " << _mRefDim << endl; 
	cout << "_sigmaRefDim       : " << _sigmaRefDim << endl; 
	cout << "_numberOfComponents: " << _numberOfComponents << endl; 
	for(int i=0; i<_numberOfComponents; i++){
		cout << "_numMolsPerComp["<<i<<"] : " << _numMolsPerComp[i] << endl; 
	}
	
	cout << "_rho               : " << _rho << endl; 
	cout << "_temperature       : " << _temperature << endl; 
	cout << "_simBoxLength[0]   : " << _simBoxLength[0] << endl; 
	cout << "_simBoxLength[1]   : " << _simBoxLength[1] << endl; 
	cout << "_simBoxLength[2]   : " << _simBoxLength[2] << endl; 

}

double PartGen::getTemperature(){
	return _temperature;
}
double PartGen::getLSimBox(int dim){
	return _simBoxLength[dim];
}
double PartGen::getEpsilonRF(){
	return _epsilonRF;
}
int PartGen::getNumComps(){
	return _numberOfComponents;
}

int PartGen::getNumSites(int comp){
	return _numSites[comp];
}
int PartGen::getNumDipoles(int comp){
	return _numDipoles[comp];
}
int PartGen::getNumQuadrupoles(int comp){
	return _numQuadrupoles[comp];
}

int PartGen::getNumCharges(int comp){
	return _numCharges[comp];
}
int PartGen::getNumTersoff(int comp){
	return _numTersoff[comp];
}

double PartGen::getSitePos(int comp, int site, int dim){
	return _rSiteBody[comp][site][dim];
}
double PartGen::getDipolePos(int comp, int site, int dim){
	return _rDipoleBody[comp][site][dim];
}
double PartGen::getQuadrupolePos(int comp, int site, int dim){
	return _rQuadrupoleBody[comp][site][dim];
}

double PartGen::getChargePos(int comp, int site, int dim)
{
	 return this->_rChargeBody[comp][site][dim];
}
double PartGen::getTersoffPos(int comp, int site, int dim)
{
	 return this->_rTersoffBody[comp][site][dim];
}

double PartGen::getEpsilonSite(int comp, int site){
	return _epsilonSite[comp][site];
}
bool PartGen::getShiftSite(int comp, int site){
	return _shiftSite[comp][site];
}
double PartGen::getMSite(int comp, int site){
	return _mSite[comp][site];
}
double PartGen::getSigmaSite(int comp, int site){
	return _sigmaSite[comp][site];
}
double PartGen::getEMyBody(int comp, int site, int dim){
	return _eMyBody[comp][site][dim];
}
double PartGen::getAbsMy(int comp, int site){
	return _absMy[comp][site];
}
double PartGen::getEQBody(int comp, int site, int dim){
	return _eQBody[comp][site][dim];
}
double PartGen::getAbsQ(int comp, int site){
	return _absQ[comp][site];
}
double PartGen::getIDummy(int comp, int dim){
	return _iBodyDummy[comp][dim];
}
double PartGen::getEta(int comp1, int comp2){
	return _etaLB[comp1][comp2];
}
double PartGen::getXi(int comp1, int comp2){
	return _xiLB[comp1][comp2];
}

double PartGen::getChargeMass(int comp, int site)
{
	 return this->_mCharge[comp][site];
}
double PartGen::getCharge(int comp, int site)
{
	 return this->_qCharge[comp][site];
}
double PartGen::getTersoffMass(int comp, int site)
{
	 return this->_mTersoff[comp][site];
}
double PartGen::getTersoffA(int comp, int site)
{
	 return this->_ATersoff[comp][site];
}
double PartGen::getTersoffB(int comp, int site)
{
	 return this->_BTersoff[comp][site];
}
double PartGen::getTersoff_lambda(int comp, int site)
{
	 return this->_lambdaTersoff[comp][site];
}
double PartGen::getTersoff_mu(int comp, int site)
{
	 return this->_muTersoff[comp][site];
}
double PartGen::getTersoffR(int comp, int site)
{
	 return this->_RTersoff[comp][site];
}
double PartGen::getTersoffS(int comp, int site)
{
	 return this->_STersoff[comp][site];
}
double PartGen::getTersoff_c(int comp, int site)
{
	 return this->_cTersoff[comp][site];
}
double PartGen::getTersoff_d(int comp, int site)
{
	 return this->_dTersoff[comp][site];
}
double PartGen::getTersoff_h(int comp, int site)
{
	 return this->_hTersoff[comp][site];
}
double PartGen::getTersoff_n(int comp, int site)
{
	 return this->_nTersoff[comp][site];
}
double PartGen::getTersoff_beta(int comp, int site)
{
	 return this->_betaTersoff[comp][site];
}

void PartGen::addParticle(int id, double x, double y, double z, ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp){

	// Component-ID
	int Comp_ = randCompID();
	vector<double> v_;
	vector<double> q_;
	vector<double> D_;

	setVectorValues(v_, 3, 0.0);
	setVectorValues(q_, 4, 0.0);
	setVectorValues(D_, 3, 0.0);

	// Velocity
	for(int dim=0; dim<3; dim++){
		v_[dim]=randdouble(-0.5, 0.5);
	}
	// Velocity Correction
	double vCorr = sqrt( 3.0*_temperature / (dotprod(v_,v_)*_massOfComps[Comp_]) );
	vecmult(v_, vCorr);

	// Orientation
	//  vector<vector<double> > quaternions;
	//  setMatrixValues(quaternions, numberOfMolecules, 4, 0.0);
	for(int q=0; q<4; q++){
		q_[q] = randchoice(-0.5, 0.5);
	}

	// Rotation
	vector<double> wBody;
	setVectorValues(wBody, 3, 0.0);
	vector<vector<double> > AT;
	setMatrixValues(AT, 3, 3, 0.0);
	double sumIw2 = 0;
	for(int dim=0; dim<3; dim++){
		wBody[dim]=randdouble(-.5,.5);
		sumIw2 += _iBody[Comp_][dim]*pow(wBody[dim],2);
	}

	if(sumIw2==0.0){
		vecmult(wBody, 0.0);
	}
	else{
		double wKorr = sqrt(_degreesOfFreedom[Comp_]*_temperature/sumIw2 );
		vecmult(wBody, wKorr);
	}

	AT[0][0] = q_[0]*q_[0] + q_[1]*q_[1] - q_[2]*q_[2] - q_[3]*q_[3];
	AT[0][1] = 2.0 * ( q_[1]*q_[2] - q_[0]*q_[3] );
	AT[0][2] = 2.0 * ( q_[1]*q_[3] + q_[0]*q_[2] );
	AT[1][0] = 2.0 * ( q_[1]*q_[2] + q_[0]*q_[3] );
	AT[1][1] = q_[0]*q_[0] - q_[1]*q_[1] + q_[2]*q_[2] - q_[3]*q_[3];
	AT[1][2] = 2.0 * ( q_[2]*q_[3] - q_[0]*q_[1] );
	AT[2][0] = 2.0 * ( q_[1]*q_[3] - q_[0]*q_[2] );
	AT[2][1] = 2.0 * ( q_[2]*q_[3] + q_[0]*q_[1] );
	AT[2][2] = q_[0]*q_[0] - q_[1]*q_[1] - q_[2]*q_[2] + q_[3]*q_[3];

	for(int ix=0; ix<3; ix++){
		for(int iy=0; iy<3; iy++){
			D_[ix] += AT[ix][iy]*(_iBody[Comp_][iy]*wBody[iy]);
		}
	}    
	Molecule m1 = Molecule(id, Comp_, 
	     x, y, z,
	     v_[0], v_[1], v_[2],
	     q_[0], q_[1], q_[2], q_[3],
	     D_[0], D_[1], D_[2], &domain->getComponents());
	particleContainer->addParticle(m1);
}

void PartGen::readLocalClusters(Domain* domain, DomainDecompBase* domainDecomp){

	//cout.precision(2);
	ifstream clusterStream;
	clusterStream.open(_clusterFile.c_str());
	vector<double> shiftedSphere; 
	vector<double> sphere; 
	double distanceToDomain;
	sphere.resize(4); // x y z r
	shiftedSphere.resize(4); // x y z r
	while(clusterStream) {
		clusterStream >> sphere[0] >> sphere[1] >> sphere[2] >> sphere[3];
		sphere[0]*=_simBoxLength[0];
		sphere[1]*=_simBoxLength[1];
		sphere[2]*=_simBoxLength[2];
		sphere[3]*=_simBoxLength[0]; // radius
		shiftedSphere[3] = sphere[3];

		for(int ix=-1; ix <=1; ix++){ 
			for(int iy=-1; iy <=1; iy++){ 
				for(int iz=-1; iz <=1; iz++){ 
					shiftedSphere[0] = sphere[0] + ix*_simBoxLength[0];
					shiftedSphere[1] = sphere[1] + iy*_simBoxLength[1];
					shiftedSphere[2] = sphere[2] + iz*_simBoxLength[2];
					distanceToDomain = domainDecomp->guaranteedDistance(shiftedSphere[0],
					          shiftedSphere[1], shiftedSphere[2], domain);
					if(distanceToDomain <= shiftedSphere[3]){ // sphere[3] = radius
						_localClusters.push_back(shiftedSphere);
						if(domainDecomp->getRank()==0){
							//cout << "Added cluster: " << shiftedSphere[0] << " / " << shiftedSphere[1] << " / " << shiftedSphere[2] << " , r: " << shiftedSphere[3] << endl;
						}
					}
				}
			}
		}
	}
	clusterStream.close();
}





bool PartGen::belongsToPreviousCluster(double x, double y, double z, int clusterid){
	bool belongsToOther = false;
	for(int i=0; i<clusterid; i++){
		if(sqrt(pow(x - _localClusters[i][0],2.0) +
			pow(y - _localClusters[i][1],2.0) +
			pow(z - _localClusters[i][2],2.0)) <= _localClusters[i][3]){
			belongsToOther = true;
			break;
		}
	}
	return belongsToOther;
}

bool PartGen::closeToAnyCluster(double x, double y, double z, double offset){
	bool belongsToOther = false;
	for(unsigned int i=0; i<_localClusters.size(); i++){
		if(sqrt(pow(x - _localClusters[i][0],2.0) +
			pow(y - _localClusters[i][1],2.0) +
			pow(z - _localClusters[i][2],2.0)) <= _localClusters[i][3] + offset){
			belongsToOther = true;
			break;
		}
	}
	return belongsToOther;
}




void PartGen::ignoreLines(ifstream &inpfStream, int numLines){
	for(int i=0; i<numLines; i++){
		inpfStream.ignore(1024, '\n');
	}
}

int PartGen::getIntParamValue(ifstream &inpfStream){
	string dummy; 
	int value;
	dummy = "dummy";
	while(dummy != "="){
		inpfStream >> dummy;
	}
	inpfStream >> value;
	inpfStream.ignore(1024, '\n');
	return value;
}

double PartGen::getDoubleParamValue(ifstream &inpfStream){
	string dummy; 
	double value;
	dummy = "dummy";
	while(dummy != "="){
		inpfStream >> dummy;
	}
	inpfStream >> value;
	inpfStream.ignore(1024, '\n');
	return value;
}

void PartGen::removePrefix(ifstream &inpfStream){
	string dummy; 
	dummy = "dummy";
	while(dummy != "="){
		inpfStream >> dummy;
	}
}

void PartGen::getDoubleParamValues(ifstream &inpfStream, double &val1, double &val2){
	removePrefix(inpfStream);
	inpfStream >> val1 >> val2;
	inpfStream.ignore(1024, '\n');
}
void PartGen::getDoubleParamValues(ifstream &inpfStream, double &val1, double &val2, double &val3){
	removePrefix(inpfStream);
	inpfStream >> val1 >> val2 >> val3;
	inpfStream.ignore(1024, '\n');
}
void PartGen::getDoubleParamValues(ifstream &inpfStream, double &val1, double &val2, double &val3, double &val4){
	removePrefix(inpfStream);
	inpfStream >> val1 >> val2 >> val3 >> val4;
	inpfStream.ignore(1024, '\n');
}

void PartGen::readRefValues(ifstream &inpfStream){
	ignoreLines(inpfStream,2);
	_epsilonRefDim = getDoubleParamValue(inpfStream);
	_mRefDim = getDoubleParamValue(inpfStream);
	_sigmaRefDim = getDoubleParamValue(inpfStream);

	_epsilonRefDim *= _boltzmannDim;     // not reduced
	_mRefDim *= _atomicMassDim;    // not reduced
	_sigmaRefDim *= 1.0E-10;        // not reduced
	
	/*
	 * ACHTUNG: in der eingelesenen Datei
	 * haben eps_ref, m_ref und sig_ref
	 * die Einheiten K, u und A.
	 *
	 * Hier werden sie (eigentlich unnoetigerweise)
	 * umgestellt auf J, kg und m.
	 */
	
	ignoreLines(inpfStream, 1);
}

void PartGen::readStatePoint(ifstream &inpfStream){
	ignoreLines(inpfStream,3);
	_numberOfComponents = getIntParamValue(inpfStream);
	_numMolsPerComp.resize(_numberOfComponents);
	_numberOfMolecules = 0;
	for(int i=0; i<_numberOfComponents; i++){
		_numMolsPerComp[i] = getIntParamValue(inpfStream);
		_numberOfMolecules += _numMolsPerComp[i];
	}
	_rho = getDoubleParamValue(inpfStream);               // not reduced
	_rho *= 1.0E+03*pow(_sigmaRefDim,3)*_avogadroDim;       // reduced
	_temperature = getDoubleParamValue(inpfStream);       // not reduced;
	_temperature *= _boltzmannDim/_epsilonRefDim;           // reduced
	ignoreLines(inpfStream, 1);
}

void PartGen::readAlgorithm(ifstream &inpfStream, 
							vector<double> &simBoxRatio){
	ignoreLines(inpfStream,2);
	getDoubleParamValue(inpfStream);   // dt, not reduced
	//_dt *= 1.0E-15 * sqrt(_epsilonRefDim/_mRefDim)/_sigmaRefDim; // reduced
	this->_LJCutoffRadius = getDoubleParamValue(inpfStream);  // rc, not reduced
	//_rc *= 1.0E-10 / _sigmaRefDim;          // reduced
	// THIS VERSION DOES NOT PROCESS EquiSchritte AND SimSchritte
	ignoreLines(inpfStream,2);
	simBoxRatio.resize(3);
	for(int i=0; i<3; i++){
		simBoxRatio[i] = getDoubleParamValue(inpfStream);
	}
	// THIS VERSION DOES NOT PROCESS Output
	ignoreLines(inpfStream, 2);
}

void PartGen::set3DMatrixValues(vector<vector<vector<double> > > &matrix, int numx, int numy, int numz, double value){
	matrix.resize(numx);
	for(int i=0; i<numx; i++){
		matrix[i].resize(numy);
		for(int j=0; j<numy; j++){
			matrix[i][j].resize(numz);
			for(int k=0; k<numz; k++){
	matrix[i][j][k] = value;
	    }
		}
	}
}

void PartGen::setMatrixValues(vector<vector<double> > &matrix, int numrows, int numcols, double value){
	matrix.resize(numrows);
	for(int i=0; i<numrows; i++){
		matrix[i].resize(numcols);
		for(int j=0; j<numcols; j++){
			matrix[i][j] = value;
		}
	}
}

void PartGen::setVectorValues(vector<int> &vect, int num, int value){
	vect.resize(num);
	for(int i=0; i<num; i++){
		vect[i] = value;
	}
}

void PartGen::setVectorValues(vector<double> &vect, int num, double value){
	vect.resize(num);
	for(int i=0; i<num; i++){
		vect[i] = value;
	}
}

double PartGen::dotprod(vector<double> &v1, vector<double> &v2){
	if(v1.size() != v2.size()){
		cout << "SIZE OF TWO VECTOR FOR dotprod NOT EQUAL!" << endl;
		exit(1);
	}
	double erg = 0;
	for(unsigned int i=0; i<v1.size(); i++){
		erg += v1[i]*v2[i];
	}
	return erg;
}

void PartGen::vecmult(vector<double> &v, double value){
	for(unsigned int i=0; i<v.size(); i++){
		v[i] *= value;
	}
}

double PartGen::max(double a, double b, double c){
	double max;
	if(a>b) max = a;
	else max = b;
	if(c>max) max=c;
	return max;
}

void PartGen::vecadd(vector<double> &v1, vector<double> &v2){
	if(v1.size() != v2.size()){
		cout << "SIZE OF TWO VECTOR FOR vecadd NOT EQUAL!" << endl;
		exit(1);
	}
	for(unsigned int i=0; i<v1.size(); i++){
		v1[i] += v2[i];
	}
}

void PartGen::matmult(vector<vector<double> > &m, double value){
	for(unsigned int i=0; i<m.size(); i++){
		for(unsigned int j=0; j<m[i].size(); j++){
			m[i][j] *= value;
		}
	}
}

vector<vector<double> > PartGen::matmult(vector<vector<double> > &m1, vector<vector<double> > &m2){
	vector<vector<double> > erg;
	erg.resize(3);
	for(int i=0; i<3; i++){
		erg[i].resize(3);
	}
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			erg[i][j] = 0.0;
			for(int k=0; k<3; k++){
	erg[i][j] += m1[i][k] * m2[k][j];
	    }
		}
	}
	return erg;
}

void PartGen::matmult(vector<double> &v, vector<vector<double> > &m){
	vector<double> erg;
	erg.resize(3);
	for(int i=0; i<3; i++){
		erg[i] = 0.0;
		for(int j=0; j<3; j++){
			erg[i] += v[j] * m[j][i];
		}
	}
	for(int i=0; i<3; i++){
		v[i] = erg[i];
	}
}


void PartGen::mat3dmult(vector<vector<vector<double> > > &m, double value){
	for(unsigned int i=0; i<m.size(); i++){
		for(unsigned int j=0; j<m[i].size(); j++){
			for(unsigned int k=0; k<m[i][j].size(); k++){
	m[i][j][k] *= value;
	    }
		}
	}
}


void PartGen::solveLGS(vector<vector<double> > &A, vector<double> &x){
	x.resize(3);
	
	double faktor = 0;
	for(int j = 0; j<3; j++){
		// Pivot-Suche: groesstes Element in der aktuellen Spalte (unterhald der Diagonalen) ermitteln
		int indexpivot = j;
		double pivot = fabs(A[j][j]);
		for(int i=j+1; i<3; i++){
			if(fabs(A[i][j])>pivot){
	indexpivot = i;
	pivot = fabs(A[i][j]);
	    }
		}
		// Wenn das Pivotelement nicht in der aktuellen Zeile ist, Zeilen tauschen 
		if(indexpivot!=j){
			double temp;
			//  Zeilen in der Matirx A tauschen 
			for(int i=j; i<3; i++){
	temp = A[j][i];
	A[j][i] = A[indexpivot][i];
	A[indexpivot][i] = temp;
	    }
	    //          /* zugehoerige Elemente im Vektor b*/
	    //         temp = this.b.getElement(j);
	    //         this.b.setElement(j,this.b.getElement(indexpivot));
	    //         this.b.setElement(indexpivot,temp);
		}
		/* Eliminationsschritt */
		for(int i=j+1; i<3; i++){
			faktor = A[i][j] / A[j][j];
			for(int k=j+1; k<3; k++){
	A[i][k] = A[i][k] - faktor*A[j][k];
	    }
	    //this.b.setElement(i,this.b.getElement(i)-faktor*this.b.getElement(j));
		}
	}

	/* Substitution */
	for(int i=2; i>=0; i--){
		//x[i] = b[i];
		x[i] = 0;
		for(int j =i+1; j<3; j++){
			x[i] = x[i] - A[i][j]*x[j];
		}
		if(fabs(A[i][i]) > 1E-12){
			x[i] = x[i]/A[i][i];
		}
		else{
			x[i] = 1.0;
		}
	}
	double norm = sqrt(pow(x[0],2) + pow(x[1],2) + pow(x[2],2));
	for(int i=0; i<3; i++){
		x[i] /= norm;
	}
}

void PartGen::getEigenvecs(vector<vector<double> > &m, vector<vector<double> > &eigenvecs){
	if(m.size() != 3){
		cout << "WRONG MATRIX SIZE FOR eigenvecs!" << endl;
		exit(1);
	}
	eigenvecs.resize(3);
	for(int i=0; i<3; i++){
		eigenvecs[i].resize(3);
	}

	if(fabs(m[0][1]) < 1E-13 && fabs(m[0][2]) < 1E-13 && fabs(m[1][0]) < 1E-13 &&
		 fabs(m[1][2]) < 1E-13 && fabs(m[2][0]) < 1E-13 && fabs(m[2][1]) < 1E-13 ){
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
	if(i==j){
		eigenvecs[i][j] = 1.0;
	}
	else{
		eigenvecs[i][j] = 0.0;
	}
	    }
		}
		return;
	}

	// determinant of matrix:
	// | x    r    s |        | x-l  r   s  |
	// | r    y    t |  ==>   |  r  y-l  t  |
	// | s    t    z |        |  s   t  z-l |
	// ==> characteristic polynomial:
	// -1*l^3 + (x+y+z)*l^2 + (r^2+s^2+t^2-xy-xz-yz)*l + (xyz+2rst-r^2z-s^2y-t^2x)
	// = a*l^3 + b*l^2 + c*l + d 
	double a = -1; 
	double b = m[0][0] + m[1][1] + m[2][2]; 
	double c = pow(m[0][1],2) + pow(m[0][2],2) + pow(m[1][2],2)
		- m[0][0]*m[1][1] - m[0][0]*m[2][2] - m[1][1]*m[2][2];
	double d = m[0][0]*m[1][1]*m[2][2] + 2*m[0][1]*m[0][2]*m[1][2]
		- pow(m[0][1],2)*m[2][2] - pow(m[0][2],2)*m[1][1] - pow(m[1][2],2)*m[0][0];
	
	double p = 3*a*c-b*b;
	double q = 27*a*a*d - 9*a*b*c + 2*b*b*b;
	
	if(q*q+4*p*p*p >= 0){
		cout << "COULD NOT DETERMINE EIGENVALUES!" << endl;
		exit(1);
	}
	double phi = acos(-q/(2*sqrt(-p*p*p)));
	double y1 =  2*sqrt(-p)*cos(phi/3);
	double y2 = -2*sqrt(-p)*cos(phi/3+acos(-1.)/3);
	double y3 = -2*sqrt(-p)*cos(phi/3-acos(-1.)/3);

	vector<double> eigenvalues;
	eigenvalues.resize(3);
	eigenvalues[0] = (y1-b)/(3*a);
	eigenvalues[1] = (y2-b)/(3*a);
	eigenvalues[2] = (y3-b)/(3*a);
	
	for(int eigenval=0; eigenval<3; eigenval++){

		vector<vector<double> > tmat;
		vector<double> eigenvec;
		tmat.resize(3);
		eigenvec.resize(3);
		for(int i=0; i<3; i++){
			tmat[i].resize(3);
			for(int j=0; j<3; j++){
	tmat[i][j] = m[i][j];
	    }
	    tmat[i][i] -= eigenvalues[eigenval];
		}
		    
		solveLGS(tmat, eigenvec);
		// EACH ROW IS A EIGENVECTOR! (usally, columns are eigenvectors)
		for(int i=0; i<3; i++){
			eigenvecs[eigenval][i] = eigenvec[i];
		}
	}
}

void PartGen::principleAxisTransform(){
	setVectorValues(_degreesOfFreedom, _numberOfComponents, 0);
	setMatrixValues(_iBody, _numberOfComponents, 3, 0.0);
	setMatrixValues(_invIBody, _numberOfComponents, 3, 0.0);
	
	for(int comp=0; comp<_numberOfComponents; comp++){
		vector<double> dr;
		setVectorValues(dr, 3, 0.0);
		// Massenmittelpunkt berechnen
		double massSum = 0.0;
		for(int site=0; site<_numSites[comp];site++){
			massSum += _mSite[comp][site];
			dr[0] -= _rSiteBody[comp][site][0] * _mSite[comp][site];
			dr[1] -= _rSiteBody[comp][site][1] * _mSite[comp][site];
			dr[2] -= _rSiteBody[comp][site][2] * _mSite[comp][site];
		}
		vecmult(dr, 1.0 / massSum);

		for(int site=0; site<_numSites[comp]; site++){
			vecadd(_rSiteBody[comp][site], dr);
		}
		for(int dip=0; dip<_numDipoles[comp]; dip++){
			vecadd(_rDipoleBody[comp][dip], dr);
		}
		for(int quad=0; quad<_numQuadrupoles[comp]; quad++){
			vecadd(_rQuadrupoleBody[comp][quad], dr);
		}

		vector<vector<double> > A;
		createTraegheitsMatrix(A, _rSiteBody[comp], _mSite[comp], _numSites[comp]);

		vector<vector<double> > vecs;
		getEigenvecs(A, vecs);

		//     if(comp == 3){
		//       vecs[0][0] = -0.94747869;
		//       vecs[0][1] = 0.0;
		//       vecs[0][2] = 0.31981892;
		//       vecs[1][0] = -0.3198189;
		//       vecs[1][1] = 0.0;
		//       vecs[1][2] = -0.94747869;
		//       vecs[2][0] = 0.0;
		//       vecs[2][1] = 1.0;
		//       vecs[2][2] = 0.0;
		//     }

		
		for(int site=0; site<_numSites[comp]; site++){
			matmult(_rSiteBody[comp][site], vecs);
		}
		for(int dip=0; dip< _numDipoles[comp]; dip++){
			matmult(_rDipoleBody[comp][dip],vecs);
			matmult(_eMyBody[comp][dip],vecs);
		}    
		for(int quad=0; quad< _numQuadrupoles[comp]; quad++){
			matmult(_rQuadrupoleBody[comp][quad],vecs);
			matmult(_eQBody[comp][quad],vecs);
		}

		createTraegheitsMatrix(A, _rSiteBody[comp], _mSite[comp], _numSites[comp]);
		_degreesOfFreedom[comp] = 3;
		for(int direction=0; direction<3; direction++){
			if(_iBodyDummy[comp][direction] <=0.0){
	if(A[direction][direction] <= max(A[0][0],A[1][1],A[2][2])*100.0*1.0E-16){
		_iBody[comp][direction] = 0.0;
		_invIBody[comp][direction] = 0.0;
		_degreesOfFreedom[comp] -= 1;
	}
	else{
		_iBody[comp][direction] = A[direction][direction];
		_invIBody[comp][direction] = 1.0 / _iBody[comp][direction];
	}
	    }
	    else{
	_iBody[comp][direction] = _iBodyDummy[comp][direction];
	_invIBody[comp][direction] = 1.0 / _iBody[comp][direction];
	    }
		}
	}
}

void PartGen::createTraegheitsMatrix(vector<vector<double> > &matrix, vector<vector<double> > &sitesPos, vector<double> &masses, int numsites){
	setMatrixValues(matrix, 3, 3, 0.0);
	
	for(int site=0; site<numsites; site++){
		// Haupttrdgheitsmomente
		matrix[0][0] += (pow(sitesPos[site][1],2) + pow(sitesPos[site][2],2)) * masses[site];
		matrix[1][1] += (pow(sitesPos[site][0],2) + pow(sitesPos[site][2],2)) * masses[site];
		matrix[2][2] += (pow(sitesPos[site][0],2) + pow(sitesPos[site][1],2)) * masses[site];
		// Deviationsmomente
		matrix[0][1] += - (sitesPos[site][0] * sitesPos[site][1]) * masses[site];
		matrix[0][2] += - (sitesPos[site][0] * sitesPos[site][2]) * masses[site];
		matrix[1][2] += - (sitesPos[site][1] * sitesPos[site][2]) * masses[site];
		// Symmetrie
		matrix[1][0] = matrix[0][1];
		matrix[2][0] = matrix[0][2];
		matrix[2][1] = matrix[1][2];
	}
}

int PartGen::randint(int a, int b){
	if(b<a){
		cout << "ERROR in randint(int a, int b)" << endl;
		exit(1);
	}
	return a + ((int) ((double (b+1-a))*rand()/RAND_MAX));
}

double PartGen::randchoice(double a, double b){
	if(rand() > RAND_MAX/2) return a;
	else return b;
}

double PartGen::randdouble(double a, double b){
	return a+rand()*(b-a)/(RAND_MAX);
}

int PartGen::randCompID(){
	int randNum = randint(0, _numberOfMolecules-1);
	int tempNumMols = 0;
	for(int i=0; i<_numberOfComponents; i++){
		tempNumMols += _numMolsPerComp[i];
		if (randNum < tempNumMols){
			return i;
		}
	}
	cerr << "Error in randCompID, to many molecules?"  << endl;
	exit(1);
	return -1; /* Silence warnings about missing return statement */
}



// int main(int argc, char *argv[]){

//   // 7) Alle noch undefinierten Felder mit ALLOCATE erstellen
//   //***************************************************************************
//   setMatrixValues(_r_, _numberOfMolecules, 3, 0.0);
//   setMatrixValues(_v_, _numberOfMolecules, 3, 0.0);
//   setMatrixValues(_D_, _numberOfMolecules, 3, 0.0);
//   setVectorValues(_Komp_, _numberOfMolecules,0);


//   // StartKonfiguration.f90
//   //***************************************************************************
//   //**  1) Anfangspositionen (Translation: kfz-Gitter;  Rotation: beliebig)  **
//   //**  2) Anfangsgeschwindigkeiten, -drall per Zufallsgenerator             **
//   //**     Gesamtimpuls zu Null bestimmen (Translation/Rotation)             **
//   //**  3) Erstellen der Anfangs-Liste                                       **
//   //**  4) Nullter-Schritt: Krdfte/Momente = 0                               **
//   //***************************************************************************
//   //
//   //random.seed()
//   //
//   //! 1) Anfangspositionen (Translation: kfz-Gitter;  Rotation: beliebig)
//   //!**************************************************************************
//   //! Translation

//   int molid;
//   // zufdllige Komponenten-Molek|l Zuordnung
//   for(int comp=0; comp<numberOfComponents-1; comp++){
//     for(int mol=0; mol<numMolsPerComp[comp]; mol++){
//       molid=randint(0,numberOfMolecules-1);
//       while(_Komp_[molid]!=0){
//   molid=randint(0,numberOfMolecules-1);   
//       }
//       _Komp_[molid] = comp+1;
//     }
//   }
//   for(int mol=0; mol<numberOfMolecules; mol++){
//     if(_Komp_[mol]==0){
//       _Komp_[mol] = numberOfComponents;
//     }
//   }

//   // Rotation


//   //cout << "Quaternion: " << quaternions[0][0] << endl;
	

//   // 2) Anfangsgeschwindigkeiten, -drall per Zufallsgenerator
//   //    Gesamtimpuls zu Null bestimmen (Translation/Rotation)
//   //*****************************************************************************

//   vector<double> dv;
//   setVectorValues(dv, 3, 0.0);
//   for(int mol=0; mol<numberOfMolecules; mol++){
//      // Impuls = 0 setzen
//      for(int dim=0; dim<3;dim++){
//        dv[dim] += _v_[mol][dim];
//      }
//   }
//   // Impuls = 0 setzen
//   for(int dim=0; dim<3; dim++){
//     dv[dim] = -dv[dim]/numberOfMolecules;
//   }
//   for(int mol=0; mol<numberOfMolecules; mol++){
//     for(int dim=0; dim<3; dim++){
//       _v_[mol][dim] += dv[dim];
//     }
//   }

//   // Nachkorrektur der Temperatur, da Verdnderung durch Impuls=0
//   double sumMv2 = 0.0;
//   for(int mol=0; mol<numberOfMolecules; mol++){
//     sumMv2 += dotprod(_v_[mol],_v_[mol])*massOfComps[_Komp_[mol]-1];
//   }
//   double vKorr = sqrt( 3.0*numberOfMolecules*temperature/sumMv2 );
//   matmult(_v_, vKorr);
//   double betaTrans = 1.0;

//   // Rotation




//   // Impuls = 0 setzen
//   int x = 0;
//   vector<double> dD;
//   setVectorValues(dD, 3, 0.0);
//   for(int mol=0; mol<numberOfMolecules; mol++){
//     if(_degreesOfFreedom[_Komp_[mol]-1]!=0){
//       for(int dim=0; dim<3; dim++){
//   dD[dim] -= _D_[mol][dim];
//       }
//       x += 1;
//     }
//   }
//   if(x==0){ //f|r "nur" 1CLJ Null im Nenner vermeiden!
//     for(int dim=0; dim<3; dim++){
//       dD[dim] = 0.0;
//     }
//   }
//   else{
//     for(int dim=0; dim<3; dim++){
//       dD[dim] = dD[dim] / ((double) x);
//     }
//   }
//   for(int mol=0; mol<numberOfMolecules; mol++){
//     if(_degreesOfFreedom[_Komp_[mol]-1]!=0){   // 1CLJ nicht korrigieren
//       for(int dim=0; dim<3; dim++){
//   _D_[mol][dim] += dD[dim];
//       }
//     }
//   }

//   // Nachkorrektur der Temperatur, da Verdnderung durch Impuls=0
//   double sumIw2 = 0.0;
//   vector<vector<double> > A;
//   setMatrixValues(A, 3, 3, 0.0);
//   for(int mol=0; mol<numberOfMolecules; mol++){
//     vector<double> q = quaternions[mol];
//     A[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
//     A[0][1] = 2.0 * ( q[1]*q[2] + q[0]*q[3] );
//     A[0][2] = 2.0 * ( q[1]*q[3] - q[0]*q[2] );
//     A[1][0] = 2.0 * ( q[1]*q[2] - q[0]*q[3] );
//     A[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
//     A[1][2] = 2.0 * ( q[2]*q[3] + q[0]*q[1] );
//     A[2][0] = 2.0 * ( q[1]*q[3] + q[0]*q[2] );
//     A[2][1] = 2.0 * ( q[2]*q[3] - q[0]*q[1] );
//     A[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
//     vector<double> v1;
//     setVectorValues(v1, 3, 0.0);
//     for(int ix=0; ix<3; ix++){
//       for(int iy=0; iy<3; iy++){
//   v1[ix] += A[ix][iy]*_D_[mol][iy];
//       }
//     }
//     double wBody2; // WHAT IS THIS?
//     wBody2 = dotprod(v1, _invIBody[_Komp_[mol]-1]);
//     for(int dim=0; dim<3; dim++){
//       sumIw2 += _iBody[_Komp_[mol]-1][dim]*pow(wBody2, 2);
//     }
//   }

	


//   double betaRot;
//   if(sumIw2==0.0){  // damit bei "nur" 1CLJ keine Null im Nenner!
//     betaRot = 1.0;
//   }
//   else{
//     int numFHGs = 0;
//     for(int comp=0; comp<numberOfComponents; comp++){
//       numFHGs += _degreesOfFreedom[comp]*numMolsPerComp[comp];
//     }
//     betaRot = sqrt( numFHGs*temperature/sumIw2 );
//   }
//   matmult(_D_, betaRot);


//}


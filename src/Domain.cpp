#include "Domain.h"

#include "datastructures/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include <cmath>

#include <iostream>
#include <sstream>
#include <string>

using namespace std;

// used by Domain::init_Corr() -----------
//! @todo WHAT DOES THIS METHOD DO?
double TICCu(int n,double rc,double sigma2);
//! @todo WHAT DOES THIS METHOD DO?
double TICSu(int n,double rc,double sigma2,double tau);
//! @todo WHAT DOES THIS METHOD DO?
double TISSu(int n,double rc,double sigma2,double tau1,double tau2);
//! @todo WHAT DOES THIS METHOD DO?
double TICCv(int n,double rc,double sigma2);
//! @todo WHAT DOES THIS METHOD DO?
double TICSv(int n,double rc,double sigma2,double tau);
//! @todo WHAT DOES THIS METHOD DO?
double TISSv(int n,double rc,double sigma2,double tau1,double tau2);
//----------------------------------------


Domain::Domain(int rank){
  _localRank = rank;
  _localUpot = 0;
  _localVirial = 0;   
  _globalUpot = 0;
  _globalVirial = 0; 
  _globalRho = 0;
  _globalRotDOF = 0;
  _globalLength[0] = 0;
  _globalLength[1] = 0;
  _globalLength[2] = 0;
  _globalBetaTrans = 1.0;
  _globalBetaRot = 1.0;
  _globalTemperature = 1.0;
  _localSummv2 = 0.0;
  _localSumIw2 = 0.0; 
  _currentTime = 0.0;
}

void Domain::setLocalUpot(double Upot) {_localUpot = Upot;}

double Domain::getLocalUpot() const {return _localUpot; }

void Domain::setLocalVirial(double Virial) {_localVirial = Virial;}

double Domain::getLocalVirial() const {return _localVirial; }

double Domain::getGlobalBetaTrans() const { return _globalBetaTrans; }

double Domain::getGlobalBetaRot() const { return _globalBetaRot; }

double Domain::getGlobalTemperature() const { return _globalTemperature; }

void Domain::setGlobalTemperature(double temp) { _globalTemperature = temp; }

vector<double> & Domain::getmixcoeff() { return _mixcoeff; }

double Domain::getepsilonRF() const { return _epsilonRF; }

void Domain::setepsilonRF(double erf) { _epsilonRF = erf; }

unsigned long Domain::getglobalNumMolecules() const { return _globalNumMolecules; }

void Domain::setglobalNumMolecules(unsigned long glnummol) { _globalNumMolecules = glnummol; }

double Domain::getGlobalPressure() const { return _globalTemperature*_globalRho+_globalRho*getAverageGlobalVirial()/3.; }

double Domain::getAverageGlobalVirial() const { return _globalVirial/_globalNumMolecules; }

double Domain::getAverageGlobalUpot() const { return _globalUpot/_globalNumMolecules; }

void Domain::setLocalSummv2(double summv2){ _localSummv2 = summv2; }
    
void Domain::setLocalSumIw2(double sumIw2){ _localSumIw2 = sumIw2; } 

int Domain::getlocalRank(){ return _localRank;}

unsigned long Domain::getinpversion(){ return _inpversion;}

void Domain::setinpversion(unsigned long inpv){ _inpversion = inpv;}

double Domain::getglobalRho(){ return _globalRho;}

void Domain::setglobalRho(double grho){ _globalRho = grho;}

unsigned long Domain::getglobalRotDOF(){ return _globalRotDOF;}

void Domain::setglobalRotDOF(unsigned long grotdof){ _globalRotDOF = grotdof;}

double Domain::getCurrentTime(){ return _currentTime;}

void Domain::setCurrentTime(double curtime){ _currentTime = curtime;}

void Domain::advanceTime(double timestep){ _currentTime += timestep;}

vector<Component>& Domain::getComponents(){
  return _components; 
}

void Domain::addComponent(Component component){
  _components.push_back(component);
}

Comp2Param& Domain::getComp2Params(){
  return _comp2params; 
}

double Domain::getGlobalLength(int index) const {
  return _globalLength[index];
}

void Domain::setGlobalLength(int index, double length) {
  _globalLength[index] = length;
}

void Domain::calculateGlobalValues(DomainDecompBase* domainDecomp,
     ParticleContainer* particleContainer){
  // number of molecules on the local process. After the reduce operation
  // num_molecules will contain the global number of molecules
  unsigned long numMolecules = particleContainer->getNumberOfParticles();
  double Upot = _localUpot;
  double Virial = _localVirial;
  double summv2 = _localSummv2;
  double sumIw2 = _localSumIw2;

  // To calculate Upot, Ukin and Pressure, intermediate values from all      
  // processes are needed. Here the         
  // intermediate values of all processes are summed up so that the root    
  // process can calculate the final values. to be able to calculate all     
  // values at this point, the calculation of the intermediate value sum_v2  
  // had to be moved from Thermostat to upd_postF and the final calculations  
  // of m_Ukin, m_Upot and Pressure had to be moved from Thermostat / upd_F  
  // to this point           
  domainDecomp->collCommInit(5);
  domainDecomp->collCommAppendDouble(Upot);
  domainDecomp->collCommAppendDouble(Virial);
  domainDecomp->collCommAppendDouble(summv2);
  domainDecomp->collCommAppendDouble(sumIw2);
  domainDecomp->collCommAppendUnsLong(numMolecules);
  domainDecomp->collCommAllreduceSum();
  Upot = domainDecomp->collCommGetDouble();
  Virial = domainDecomp->collCommGetDouble();
  summv2 = domainDecomp->collCommGetDouble();
  sumIw2 = domainDecomp->collCommGetDouble();
  numMolecules = domainDecomp->collCommGetUnsLong();
  domainDecomp->collCommFinalize();

  // Process 0 has to add the dipol correction:
  // m_UpotCorr and m_VirialCorr already contain constant (internal) dipol correction
  _globalUpot = Upot + _UpotCorr;
  _globalVirial = Virial + _VirialCorr;

  //! @todo comment
  _globalBetaTrans=sqrt(3.*numMolecules*_globalTemperature/summv2);
  _globalBetaRot=1.0;
  //! @todo comment
  if(sumIw2!=0.) _globalBetaRot=sqrt(_globalRotDOF*_globalTemperature/sumIw2);    
}

void Domain::calculateVelocitySums(ParticleContainer* partCont){
  Molecule* tempMolecule;
  _localSummv2 = 0.0;
  _localSumIw2 = 0.0;
  
  
  for(tempMolecule = partCont->begin(); tempMolecule != partCont->end(); tempMolecule = partCont->next()){
    tempMolecule->calculate_mv2_Iw2(_localSummv2, _localSumIw2);
  }
}


void Domain::writeCheckpoint(string filename, ParticleContainer* particleContainer,
                             DomainDecompBase* domainDecomp) const {
                     
  ofstream checkpointfilestream(filename.c_str());
  checkpointfilestream << "MDProject\t20070111"<< endl;
  checkpointfilestream << "currentTime\t" << _currentTime << endl;
  checkpointfilestream << "Temperature\t" << _globalTemperature << endl;
  checkpointfilestream << "Length\t" << _globalLength[0] << "\t" << _globalLength[1] << "\t" << _globalLength[2] << endl;
  checkpointfilestream << "NumberOfComponents\t" << _components.size() << endl;
  for(vector<Component>::const_iterator pos=_components.begin();pos!=_components.end();++pos){
    pos->write(checkpointfilestream);
  }
  unsigned int numperline=_components.size();
  unsigned int iout=0;
  for(vector<double>::const_iterator pos=_mixcoeff.begin();pos!=_mixcoeff.end();++pos){
    checkpointfilestream << *pos;
    iout++;
    // 2 parameters (xi and eta)
    if(iout/2>=numperline) {
      checkpointfilestream << endl;
      iout=0;
      --numperline;
    }
    else if(!(iout%2)) {
      checkpointfilestream << "\t";
    }
    else {
      checkpointfilestream << " ";
    }
  }

  checkpointfilestream << _epsilonRF << endl;
  checkpointfilestream << "NumberOfMolecules\t" << _globalNumMolecules << endl;

  checkpointfilestream << "MoleculeFormat\t" << "\tICRVQD" << endl;
  checkpointfilestream.close();
  
  domainDecomp->writeMoleculesToFile(filename, particleContainer);


}

void Domain::initParameterStreams(double cutoffRadius){
  _comp2params.initialize(_components,_mixcoeff,_epsilonRF, cutoffRadius); 
}

void Domain::initFarFieldCorr(double cutoffRadius) {
  double UpotCorrLJ=0.;
  double VirialCorrLJ=0.;
  double MySelbstTerm=0.;
  unsigned int numcomp=_components.size();
  unsigned long nummolecules=0;
  for(unsigned int i=0;i<numcomp;++i) {
    Component& ci=_components[i];
    nummolecules+=ci.numMolecules();
    unsigned int numljcentersi=ci.numLJcenters();
    unsigned int numdipolesi=ci.numDipoles();
    for(unsigned int j=0;j<numcomp;++j) {
      Component& cj=_components[j];
      unsigned int numljcentersj=cj.numLJcenters();
      ParaStrm& params=_comp2params(i,j);
      params.reset_read();
      // LJ centers
      for(unsigned int si=0;si<numljcentersi;++si) {
        double xi=ci.ljcenter(si).rx();
        double yi=ci.ljcenter(si).ry();
        double zi=ci.ljcenter(si).rz();
        double tau1=sqrt(xi*xi+yi*yi+zi*zi);
        for(unsigned int sj=0;sj<numljcentersj;++sj) {
          double xj=cj.ljcenter(sj).rx();
          double yj=cj.ljcenter(sj).ry();
          double zj=cj.ljcenter(sj).rz();
          double tau2=sqrt(xj*xj+yj*yj+zj*zj);
          if(tau1+tau2>=cutoffRadius){
            cerr << "Error calculating cutoff corrections, rc too small" << endl;
          }
          double eps24;
          params >> eps24;
          double sig2;
          params >> sig2;
          double fac=double(ci.numMolecules())*double(cj.numMolecules())*eps24;
          if(tau1==0. && tau2==0.){
            UpotCorrLJ+=fac*(TICCu(-6,cutoffRadius,sig2)-TICCu(-3,cutoffRadius,sig2));
            VirialCorrLJ+=fac*(TICCv(-6,cutoffRadius,sig2)-TICCv(-3,cutoffRadius,sig2));
          }
          else if(tau1!=0. && tau2!=0.) {
            UpotCorrLJ+=fac*(TISSu(-6,cutoffRadius,sig2,tau1,tau2)-TISSu(-3,cutoffRadius,sig2,tau1,tau2));
            VirialCorrLJ+=fac*(TISSv(-6,cutoffRadius,sig2,tau1,tau2)-TISSv(-3,cutoffRadius,sig2,tau1,tau2));
          }
          else {
            if(tau2==0.) {
              tau2=tau1;
            }
            UpotCorrLJ+=fac*(TICSu(-6,cutoffRadius,sig2,tau2)-TICSu(-3,cutoffRadius,sig2,tau2));
            VirialCorrLJ+=fac*(TICSv(-6,cutoffRadius,sig2,tau2)-TICSv(-3,cutoffRadius,sig2,tau2));
          }
        }
      }
      // Dipoles
      if(i==j){
        double summy2=0.;
        for(unsigned int si=0;si<numdipolesi;++si){
          // i==j
          for(unsigned int sj=0;sj<numdipolesi;++sj){
            double my2;
            params >> my2;
            summy2+=my2;
          }
        }
        // ci==cj
        MySelbstTerm+=summy2*ci.numMolecules();
      }
    }
  }

  double fac=M_PI*_globalRho/(3.*_globalNumMolecules);
  UpotCorrLJ*=fac;
  VirialCorrLJ*=-fac;
        
  double epsRFInvrc3=2.*(_epsilonRF-1.)/((cutoffRadius*cutoffRadius*cutoffRadius)*(2.*_epsilonRF+1.));
  MySelbstTerm*=-0.5*epsRFInvrc3;

  _UpotCorr=UpotCorrLJ+MySelbstTerm;
  _VirialCorr=VirialCorrLJ+3.*MySelbstTerm;

}

// Helper functions ================================================================================
// used by Domain::init_Corr()

double TICCu(int n,double rc,double sigma2) {
  return -pow(rc,2*n+3) / (pow(sigma2,n)*(2*n+3));
}

double TICSu(int n,double rc,double sigma2,double tau){
  return -( pow(rc+tau,2*n+3) - pow(rc-tau,2*n+3) ) * rc / ( 4*pow(sigma2,n)*tau*(n+1)*(2*n+3) ) +  ( pow(rc+tau,2*n+4) - pow(rc-tau,2*n+4) ) / ( 4*pow(sigma2,n)*tau*(n+1)*(2*n+3)*(2*n+4) );
}
    
double TISSu(int n,double rc,double sigma2,double tau1,double tau2){
  double tauMinus,tauPlus;
  tauPlus = tau1+tau2;
  tauMinus = tau1-tau2;
  return -(   pow(rc+tauPlus,2*n+4) - pow(rc+tauMinus,2*n+4) - pow(rc-tauMinus,2*n+4) + pow(rc-tauPlus,2*n+4) ) * rc / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3)*(2*n+4) ) +  (   pow(rc+tauPlus,2*n+5) - pow(rc+tauMinus,2*n+5) - pow(rc-tauMinus,2*n+5) + pow(rc-tauPlus,2*n+5) ) / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3)*(2*n+4)*(2*n+5) );
}

double TICCv(int n,double rc,double sigma2){
  return 2*n * TICCu(n,rc,sigma2);
}

double TICSv(int n,double rc,double sigma2,double tau){
  return -( pow(rc+tau,2*n+2) - pow(rc-tau,2*n+2) ) * rc*rc / ( 4*pow(sigma2,n)*tau*(n+1) ) - 3*TICSu(n,rc,sigma2,tau);
}

double TISSv(int n,double rc,double sigma2,double tau1,double tau2){
  double tauMinus,tauPlus;
  tauPlus = tau1+tau2;
  tauMinus = tau1-tau2;
  return -(   pow(rc+tauPlus,2*n+3) - pow(rc+tauMinus,2*n+3) - pow(rc-tauMinus,2*n+3) + pow(rc-tauPlus,2*n+3) ) * rc*rc / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3) ) - 3*TISSu(n,rc,sigma2,tau1,tau2);
}



#include "thermostats/VelocityScalingThermostat.h"

#include "Domain.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "utils/Logger.h"

#include <iomanip>

using namespace std;
using Log::global_log;

VelocityScalingThermostat::VelocityScalingThermostat() : _globalBetaTrans(1), _globalBetaRot(1), _globalVelocity(NULL), _componentwise(false) {
	_globalVelocity = new double[3];
	_molDirVelThreshold = 50;
}

VelocityScalingThermostat::~VelocityScalingThermostat() {
	delete[] _globalVelocity;
}

void VelocityScalingThermostat::setBetaTrans(int componentId, double beta) {
	_componentBetaTrans[componentId] = beta;
}

void VelocityScalingThermostat::setBetaRot(int componentId, double beta) {
	_componentBetaRot[componentId] = beta;
}

void VelocityScalingThermostat::setGlobalVelocity(double v[3]) {
	for(int d = 0; d < 3; d++) {
		_globalVelocity[d] = v[d];
	}
}
void VelocityScalingThermostat::setVelocity(int componentId, double v[3]) {
	if( _componentVelocity.find(componentId) == _componentVelocity.end() ) {
		_componentVelocity[componentId] = new double[3];
	}
	for(int d = 0; d < 3; d++) {
		_componentVelocity[componentId][d] = v[d];
	}
}

// after the scaling factor is calculated in domain->calculateGlobalValues() the velocity is scaled  here (1D OR 3D)
// additionally the kinetic part of the stresses, the pressure and the heatflux is calculated here
void VelocityScalingThermostat::apply(ParticleContainer *moleculeContainer, DomainDecompBase *dode) {
	Molecule *molecule;
	if(_componentwise ) {
		for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()) {
			int thermostatId;
			double betaTrans = _globalBetaTrans;
			double betaRot = _globalBetaRot;
			unsigned cid = molecule->componentid();
			double alphaTrans = _globalAlphaTrans;
			thermostatId = _simulation.getDomain()->getThermostat(cid);
			if(_simulation.getDomain()->isThermostatLayer() == true || _simulation.getDomain()->isThermostatWallLayer() == true){
			  thermostatId = _simulation.getDomain()->moleculeInLayer(molecule->r(0), molecule->r(1), molecule->r(2), molecule->componentid());
			}
			betaTrans = _componentBetaTrans[thermostatId];
			betaRot   = _componentBetaRot[thermostatId];
			alphaTrans = _componentAlphaTrans[thermostatId];
			map<int, int> Dim = _simulation.getDomain()->getDim();
			map<int, double*>::iterator vIter;
			if(thermostatId != -10){
			  if((vIter = _componentVelocity.find(thermostatId)) != _componentVelocity.end()) {
			      	if(_simulation.getDomain()->isScaling1Dim(thermostatId) && _simulation.getDomain()->getAlphaTransCorrection(thermostatId) == false){
                                        molecule->scale_v_1Dim(alphaTrans, Dim[thermostatId]);
				}
				else{
                                        molecule->scale_v(betaTrans);
				}
			  }
			  else {
				double *v = vIter->second;
				molecule->vsub(v[0], v[1], v[2]);
				molecule->scale_v(betaTrans);
				molecule->vadd(v[0], v[1], v[2]);
			  }
			  molecule->scale_D(betaRot);
			}
		}
	}else{
		double betaTrans = _globalBetaTrans;
		double betaRot = _globalBetaRot;
		global_log->debug() << "Beta rot: " << betaRot << endl;
		global_log->debug() << "Beta trans: " << betaTrans << endl;
		for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()) {
			if(_globalVelocity == NULL) {
				molecule->scale_v(betaTrans);
			}
			else {
				molecule->vsub(_globalVelocity[0], _globalVelocity[1], _globalVelocity[2]);
				molecule->scale_v(betaTrans);
				molecule->vadd(_globalVelocity[0], _globalVelocity[1], _globalVelocity[2]);
			}
			molecule->scale_D(betaRot);
		}
	}
	
	for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()) {	
		// virial calculation (kinetic and force part) for stresses in solids
		if(_simulation.getDomain()->isStressCalculating(molecule->componentid()) && !(_simulation.getDomain()->getSimstep() % _simulation.getDomain()->getStressRecordTimeStep())){
		  for(int d = 0; d < 3; d++){
		    long double v_d = molecule->v(d) - molecule->getDirectedVelocityStress(d);
		    for(int e = 0; e < 3; e++){
		      long double v_e = molecule->v(e) - molecule->getDirectedVelocityStress(e);
		      long double kinTrans = v_d * v_e * molecule->mass();
		      molecule->addVirialKin(d, e, kinTrans);
		    }
		  }
		}
		if(_simulation.getDomain()->isBulkPressure(molecule->componentid())){
		  if((molecule->r(0) >= _simulation.getDomain()->getBulkBoundary(0)-_simulation.getLJCutoff() && molecule->r(0) <= _simulation.getDomain()->getBulkBoundary(1)+_simulation.getLJCutoff() && molecule->r(1) >= _simulation.getDomain()->getBulkBoundary(2)-_simulation.getLJCutoff() && molecule->r(1) <= _simulation.getDomain()->getBulkBoundary(3)+_simulation.getLJCutoff())){
		    for (unsigned short d = 0; d < 3; ++d){
			long double v_d = molecule->v(d) - molecule->getDirectedVelocityStress(d);
			long double kinTrans = v_d * v_d * molecule->mass();
			molecule->addPressureKin(d, kinTrans);
		    }
		  }
		}
		if(_simulation.getDomain()->isBarostat(molecule->componentid()) && ((_simulation.getDomain()->getSimstep() >= _simulation.getBarostatTimeInit() && _simulation.getDomain()->getSimstep() <= _simulation.getBarostatTimeEnd()))){
		  if((molecule->r(0) >= _simulation.getDomain()->getControl_bottom(0) && molecule->r(0) <= _simulation.getDomain()->getControl_top(0) && molecule->r(1) >= _simulation.getDomain()->getControl_bottom(1) && molecule->r(1) <= _simulation.getDomain()->getControl_top(1))){
		    for (unsigned short d = 0; d < 3; ++d){
			long double v_d = molecule->v(d) - molecule->getDirectedVelocity(d);
			long double kinTrans = v_d * v_d * molecule->mass();
			molecule->addPressureKin_barostat(d, kinTrans);
		    }
		  }
		}
		if(_simulation.getDomain()->isConfinement(molecule->componentid())){
		  if(molecule->r(0) >= _simulation.getDomain()->getConfinementEdge(0) && molecule->r(0) <= _simulation.getDomain()->getConfinementEdge(1) && molecule->r(1) >= _simulation.getDomain()->get_confinementMidPoint(3)-_simulation.getDomain()->getConfinementEdge(5)-0.755 && molecule->r(1) <= _simulation.getDomain()->get_confinementMidPoint(1)+_simulation.getDomain()->getConfinementEdge(5)+0.755 && !(_simulation.getDomain()->getSimstep() % _simulation.getDomain()->getConfinementRecordTimeStep())){
		    for (unsigned short d = 0; d < 3; ++d){
			long double v_d = molecule->v(d) - molecule->getDirectedVelocityConfinement(d);
			long double kinTrans = v_d * v_d * molecule->mass();
			molecule->addPressureKinConfinement(d, kinTrans);
			for(int e = 0; e < 3; e++){
			  long double v_e = molecule->v(e) - molecule->getDirectedVelocityConfinement(e);
			  long double kinTrans = v_d * v_e * molecule->mass();
			  molecule->addVirialKinConfinement(d, e, kinTrans);
			}
		    }
		  }
		}
	}
}

// calculation of the directed velocities: moving average; size of control volumes: 1*1*z_max OR other
void VelocityScalingThermostat::calculateDirectedVelocities(ParticleContainer *moleculeContainer, DomainDecompBase *dode) {
   double xMax = _simulation.getDomain()->getGlobalLength(0);
   double yMax = _simulation.getDomain()->getGlobalLength(1);
   double roundX, roundY;
   unsigned CID, xunStress_tot = 0, yunStress_tot = 0, zunStress_tot = 0, xunSlab_tot = 0, yunSlab_tot = 0, zunSlab_tot = 0, xunConf_tot = 0, yunConf_tot = 0, xun = 0, yun = 0, zun = 0;
   unsigned long xunConf = 0, yunConf = 0, xunStress = 0, yunStress = 0, zunStress = 0, xunSlab = 0, yunSlab = 0, zunSlab = 0;
   double deltaX_stress = 0.0, deltaY_stress = 0.0, deltaZ_stress = 0.0, deltaX_slab = 0.0, deltaY_slab = 0.0, deltaZ_slab = 0.0, deltaX_conf = 0.0, deltaY_conf = 0.0;
   
   string moved ("moved");
   string fixed ("fixed");
   unsigned cid_moved = _simulation.getDomain()->getCidMovement(moved, _simulation.getDomain()->getNumberOfComponents()) - 1;
   unsigned cid_fixed = _simulation.getDomain()->getCidMovement(fixed, _simulation.getDomain()->getNumberOfComponents()) - 1;
   
   // ------- PREPARATION --------
   if(_simulation.isRecordingSlabProfile()){
    xunSlab_tot = _simulation.getDomain()->getUniversalNProfileUnitsSlab(0);
    yunSlab_tot = _simulation.getDomain()->getUniversalNProfileUnitsSlab(1);
    zunSlab_tot = _simulation.getDomain()->getUniversalNProfileUnitsSlab(2);
    deltaX_slab = 1/_simulation.getDomain()->getUniversalInvProfileUnitSlab(0);
    deltaY_slab = 1/_simulation.getDomain()->getUniversalInvProfileUnitSlab(1);
    deltaZ_slab = 1/_simulation.getDomain()->getUniversalInvProfileUnitSlab(2);
   }
   if(_simulation.isRecordingStressProfile()){
    xunStress_tot = _simulation.getDomain()->getUniversalNProfileUnits_Stress(0);
    yunStress_tot = _simulation.getDomain()->getUniversalNProfileUnits_Stress(1);
    zunStress_tot = _simulation.getDomain()->getUniversalNProfileUnits_Stress(2);
    deltaX_stress = 1/_simulation.getDomain()->getUniversalInvProfileUnit_Stress(0);
    deltaY_stress = 1/_simulation.getDomain()->getUniversalInvProfileUnit_Stress(1);
    deltaZ_stress = 1/_simulation.getDomain()->getUniversalInvProfileUnit_Stress(2);
   }
   if(_simulation.isRecordingConfinementProfile()){
    xunConf_tot = _simulation.getDomain()->getUniversalNProfileUnitsConfinement(0);
    yunConf_tot = _simulation.getDomain()->getUniversalNProfileUnitsConfinement(1);
    deltaX_conf = 1/_simulation.getDomain()->getUniversalInvProfileUnitConfinement(0);
    deltaY_conf = 1/_simulation.getDomain()->getUniversalInvProfileUnitConfinement(1);
   }

   // ------- INITIALIZATION --------
    for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
     _countedMoleculesAll[cid] = 0;
     _universalCountedMoleculesAll[cid];
     for (unsigned d = 0; d < 3; d++){
	_directedVelocityAll[d][cid] = 0.0;
	_universalDirectedVelocityAll[d][cid] = 0.0;
     }
     for (double xCoord = 0.5; xCoord < xMax; xCoord++)
       for (double yCoord = 0.5; yCoord < yMax; yCoord++){
	  _countedMolecules[cid][xCoord][yCoord] = 0;
	  _universalCountedMolecules[cid][xCoord][yCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocity[d][cid][xCoord][yCoord] = 0.0;
	    _universalDirectedVelocity[d][cid][xCoord][yCoord] = 0.0;
	    if(!dode->getRank() && this->_universalPriorVelSums[d][cid][xCoord][yCoord].size() == 0){
	      this->_averagedVelSum[d][cid][xCoord][yCoord] = 0.0;
	      this->_averagedN[cid][xCoord][yCoord] = 0;
	    }
	  }
       }
    }
   
  if(!_simulation.isMovingAverageSigma2D()){
   // same for SlabProfile, StressProfile and Confinement summed up over x timesteps
   for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
    if(_simulation.isRecordingSlabProfile() && (_simulation.getSimulationStep() == 1 || (_simulation.getSimulationStep() > _simulation.getInitStatistics()))){
     for (unsigned xCoord = 0; xCoord < xunSlab_tot; xCoord++){
       for (unsigned yCoord = 0; yCoord < yunSlab_tot; yCoord++){
	 for (unsigned zCoord = 0; zCoord < zunSlab_tot; zCoord++){
	  _countedMoleculesSlab[cid][xCoord][yCoord][zCoord] = 0;
	  _universalCountedMoleculesSlab[cid][xCoord][yCoord][zCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocitySlab[d][cid][xCoord][yCoord][zCoord] = 0.0;
	    _universalDirectedVelocitySlab[d][cid][xCoord][yCoord][zCoord] = 0.0;
	    if(!dode->getRank() && this->_universalPriorVelSumsSlab[d][cid][xCoord][yCoord][zCoord].size() == 0){
	      this->_averagedVelSumSlab[d][cid][xCoord][yCoord][zCoord] = 0.0;
	      this->_averagedNSlab[cid][xCoord][yCoord][zCoord] = 0;
	    }
	  }
	 }
       }
     }
    }
    if(_simulation.isRecordingStressProfile() && (_simulation.getSimulationStep() == 1 || (_simulation.getSimulationStep() > _simulation.getInitStatistics()))){
     for (unsigned xCoord = 0; xCoord < xunStress_tot; xCoord++){
       for (unsigned yCoord = 0; yCoord < yunStress_tot; yCoord++){
	 for (unsigned zCoord = 0; zCoord < zunStress_tot; zCoord++){
	  _countedMoleculesStress[cid][xCoord][yCoord][zCoord] = 0;
	  _universalCountedMoleculesStress[cid][xCoord][yCoord][zCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocityStress[d][cid][xCoord][yCoord][zCoord] = 0.0;
	    _universalDirectedVelocityStress[d][cid][xCoord][yCoord][zCoord] = 0.0;
	    if(!dode->getRank() && this->_universalPriorVelSumsStress[d][cid][xCoord][yCoord][zCoord].size() == 0){
	      this->_averagedVelSumStress[d][cid][xCoord][yCoord][zCoord] = 0.0;
	      this->_averagedNStress[cid][xCoord][yCoord][zCoord] = 0;
	    }
	  }
	 }
       }
     }
    }
    if(_simulation.isRecordingConfinementProfile() && (_simulation.getSimulationStep() == 1 || (_simulation.getSimulationStep() > _simulation.getInitStatistics()))){
     for (unsigned xCoord = 0; xCoord < xunConf_tot; xCoord++){
       for (unsigned yCoord = 0; yCoord < yunConf_tot; yCoord++){
	  _countedMoleculesConfinement[cid][xCoord][yCoord] = 0;
	  _universalCountedMoleculesConfinement[cid][xCoord][yCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocityConfinement[d][cid][xCoord][yCoord] = 0.0;
	    _universalDirectedVelocityConfinement[d][cid][xCoord][yCoord] = 0.0;
	    if(!dode->getRank() && this->_universalPriorVelSumsConfinement[d][cid][xCoord][yCoord].size() == 0){
	      this->_averagedVelSumConfinement[d][cid][xCoord][yCoord] = 0.0;
	      this->_averagedNConfinement[cid][xCoord][yCoord] = 0;
	    }
	  }
       }
     }
    }
   }
  }
   
   // ------- DATA RECORDING --------
   Molecule *molecule;
   
   for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
      roundX = floor(molecule->r(0)) + 0.5;
      roundY = floor(molecule->r(1)) + 0.5;
      CID = molecule->componentid();
      if(CID == cid_moved || CID == cid_fixed){
	 _countedMoleculesAll[CID]++;
	for (unsigned d = 0; d < 3; d++)
		_directedVelocityAll[d][CID] += molecule->v(d);     
      }
      else{
	_countedMolecules[CID][roundX][roundY]++;
	for (unsigned d = 0; d < 3; d++)
		_directedVelocity[d][CID][roundX][roundY] += molecule->v(d);
      }
     if(!_simulation.isMovingAverageSigma2D()){ 
      if((CID != cid_moved && CID != cid_fixed) && _simulation.getSimulationStep() >= _simulation.getInitStatistics()){
       if(_simulation.isRecordingSlabProfile()){
	xun = floor(molecule->r(0) / deltaX_slab);
	yun = floor(molecule->r(1) / deltaY_slab);
	zun = floor(molecule->r(2) / deltaZ_slab);
	_countedMoleculesSlab[CID][xun][yun][zun]++;
	for (unsigned d = 0; d < 3; d++)
	  _directedVelocitySlab[d][CID][xun][yun][zun] += molecule->v(d);
       }
       if(_simulation.isRecordingStressProfile()){
	xun = floor(molecule->r(0) / deltaX_stress);
	yun = floor(molecule->r(1) / deltaY_stress);
	zun = floor(molecule->r(2) / deltaZ_stress);
	_countedMoleculesStress[CID][xun][yun][zun]++;
	for (unsigned d = 0; d < 3; d++)
	  _directedVelocityStress[d][CID][xun][yun][zun] += molecule->v(d);
       }
       if(_simulation.isRecordingConfinementProfile()){
	xun = floor((molecule->r(0)-_simulation.getDomain()->getConfinementEdge(0)) / deltaX_conf);
	yun = floor((molecule->r(1)-(_simulation.getDomain()->get_confinementMidPoint(3)-_simulation.getDomain()->getConfinementEdge(5))) / deltaY_conf);
	_countedMoleculesConfinement[CID][xun][yun]++;
	for (unsigned d = 0; d < 3; d++)
	  _directedVelocityConfinement[d][CID][xun][yun] += molecule->v(d);
       }
      }
     }
   }
   
   // ------- PARALLELIZATION --------
   //Parallelization for cid_moved and cid_fixed
   if(_simulation.getSimulationStep() >= 1){
     for (unsigned cid = 0; cid < _countedMoleculesAll.size(); cid++){
	dode->collCommInit(4);
	dode->collCommAppendUnsLong(_countedMoleculesAll[cid]);
	for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocityAll[d][cid]);
	dode->collCommAllreduceSum();
	_universalCountedMoleculesAll[cid] = dode->collCommGetUnsLong();
	for(int d = 0; d < 3; d++)
		_universalDirectedVelocityAll[d][cid] = dode->collCommGetDouble();
	dode->collCommFinalize();
     }
   }
	
   if(_simulation.getSimulationStep() >= 1){
     for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      dode->collCommInit(4*round(xMax)*round(yMax));
	for(double x = 0.5; x < xMax; x++)
	  for(double y = 0.5; y < yMax; y++){
	    dode->collCommAppendUnsLong(_countedMolecules[cid][x][y]);
	    for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocity[d][cid][x][y]);
	  }
      dode->collCommAllreduceSum();
	for(double x = 0.5; x < xMax; x++)
	  for(double y = 0.5; y < yMax; y++){
	    _universalCountedMolecules[cid][x][y] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++){
		_universalDirectedVelocity[d][cid][x][y] = dode->collCommGetDouble();
	    }
	  }
      dode->collCommFinalize();

	if(!dode->getRank())
	{
	  for(double x = 0.5; x < xMax; x++)
	      for(double y = 0.5; y < yMax; y++){
		this->_universalPriorN[cid][x][y].push_back(this->_universalCountedMolecules[cid][x][y]);
		this->_averagedN[cid][x][y] += this->_universalCountedMolecules[cid][x][y];
		for(int d = 0; d < 3; d++){
		  this->_universalPriorVelSums[d][cid][x][y].push_back(this->_universalDirectedVelocity[d][cid][x][y]);
		  this->_averagedVelSum[d][cid][x][y] += this->_universalDirectedVelocity[d][cid][x][y];
		  if(this->_universalCountedMolecules[cid][x][y] > _molDirVelThreshold)
		    this->_currentDirectedVel[d][cid][x][y] = this->_universalDirectedVelocity[d][cid][x][y]/this->_universalCountedMolecules[cid][x][y];
		  else
		    this->_currentDirectedVel[d][cid][x][y] = 0.0;
		}
	      }
	}
	
	if(!dode->getRank())
	{
	  for(double x = 0.5; x < xMax; x++)
	      for(double y = 0.5; y < yMax; y++){
		for(int d = 0; d < 3; d++){
                  if(this->_averagedN[cid][x][y] > _molDirVelThreshold)
                    this->_universalDirectedVelocity[d][cid][x][y] = this->_averagedVelSum[d][cid][x][y]/this->_averagedN[cid][x][y];
                  else
                    this->_universalDirectedVelocity[d][cid][x][y] = 0.0;  
		}
		if(this->_universalPriorVelSums[0][cid][x][y].size() == _simulation.getDirectedVelocityTime()){
		  for(int d = 0; d < 3; d++){
		    this->_averagedVelSum[d][cid][x][y] -= this->_universalPriorVelSums[d][cid][x][y].front();
		    this->_universalPriorVelSums[d][cid][x][y].pop_front();
		  }
		  this->_averagedN[cid][x][y] -= this->_universalPriorN[cid][x][y].front();
		  this->_universalPriorN[cid][x][y].pop_front();
		}
	      }
	}

	dode->collCommInit(3*round(xMax)*round(yMax));
	for(double x = 0.5; x < xMax; x++)
	      for(double y = 0.5; y < yMax; y++)
		for(int d = 0; d < 3; d++){
		  if(_simulation.isShearRate() && _simulation.getDomain()->getShearWidth() == 0.0)
			dode->collCommAppendDouble(2*this->_currentDirectedVel[d][cid][x][y] - this->_universalDirectedVelocity[d][cid][x][y]);
		  else
			dode->collCommAppendDouble(this->_universalDirectedVelocity[d][cid][x][y]);
		}
	dode->collCommBroadcast();
	for(double x = 0.5; x < xMax; x++)
	      for(double y = 0.5; y < yMax; y++)
		for(int d = 0; d < 3; d++)
		  this->_universalDirectedVelocity[d][cid][x][y] = dode->collCommGetDouble();
	dode->collCommFinalize();
     }
   }
  if(!_simulation.isMovingAverageSigma2D()){   
   if(_simulation.getSimulationStep() > _simulation.getInitStatistics()){
    //Parallelization for SlabProfile, StressProfile and Confinement
    if(_simulation.isRecordingSlabProfile()){
    for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      dode->collCommInit(4 * xunSlab_tot * yunSlab_tot * zunSlab_tot);
	for(unsigned x = 0; x < xunSlab_tot; x++)
	  for(unsigned y = 0; y < yunSlab_tot; y++)
	   for(unsigned z = 0; z < zunSlab_tot; z++){ 
	    dode->collCommAppendUnsLong(_countedMoleculesSlab[cid][x][y][z]);
	    for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocitySlab[d][cid][x][y][z]);
	  }
      dode->collCommAllreduceSum();
	for(unsigned x = 0; x < xunSlab_tot; x++)
	  for(unsigned y = 0; y < yunSlab_tot; y++)
	   for(unsigned z = 0; z < zunSlab_tot; z++){ 
	    _universalCountedMoleculesSlab[cid][x][y][z] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++)
	      _universalDirectedVelocitySlab[d][cid][x][y][z] = dode->collCommGetDouble();
	   }
      dode->collCommFinalize();
      
      if(!dode->getRank())
	{
	  for(unsigned x = 0; x < xunSlab_tot; x++)
	    for(unsigned y = 0; y < yunSlab_tot; y++)
	      for(unsigned z = 0; z < zunSlab_tot; z++){ 
		this->_universalPriorNSlab[cid][x][y][z].push_back(this->_universalCountedMoleculesSlab[cid][x][y][z]);
		this->_averagedNSlab[cid][x][y][z] += this->_universalCountedMoleculesSlab[cid][x][y][z];
		for(int d = 0; d < 3; d++){
		  this->_universalPriorVelSumsSlab[d][cid][x][y][z].push_back(this->_universalDirectedVelocitySlab[d][cid][x][y][z]);
		  this->_averagedVelSumSlab[d][cid][x][y][z] += this->_universalDirectedVelocitySlab[d][cid][x][y][z];
		  if(this->_universalCountedMoleculesSlab[cid][x][y][z] > _molDirVelThreshold)
		    this->_currentDirectedVelSlab[d][cid][x][y][z] = this->_universalDirectedVelocitySlab[d][cid][x][y][z]/this->_universalCountedMoleculesSlab[cid][x][y][z];
		  else
		    this->_currentDirectedVelSlab[d][cid][x][y][z] = 0.0;
		}
	      }
	}

	if(!dode->getRank())
	{
	  for(unsigned x = 0; x < xunSlab_tot; x++)
	    for(unsigned y = 0; y < yunSlab_tot; y++)
	      for(unsigned z = 0; z < zunSlab_tot; z++){
		for(int d = 0; d < 3; d++){
                  if(this->_averagedNSlab[cid][x][y][z]  > _molDirVelThreshold)
                    this->_universalDirectedVelocitySlab[d][cid][x][y][z] = this->_averagedVelSumSlab[d][cid][x][y][z]/this->_averagedNSlab[cid][x][y][z];
                  else
                    this->_universalDirectedVelocitySlab[d][cid][x][y][z] = 0.0;    
		}
		if(this->_universalPriorVelSumsSlab[0][cid][x][y][z].size() == _simulation.getDirectedVelocityTime()){
		  for(int d = 0; d < 3; d++){
		    this->_averagedVelSumSlab[d][cid][x][y][z] -= this->_universalPriorVelSumsSlab[d][cid][x][y][z].front();
		    this->_universalPriorVelSumsSlab[d][cid][x][y][z].pop_front();
		  }
		  this->_averagedNSlab[cid][x][y][z] -= this->_universalPriorNSlab[cid][x][y][z].front();
		  this->_universalPriorNSlab[cid][x][y][z].pop_front();
		}
	      }
	}

	dode->collCommInit(3 * xunSlab_tot * yunSlab_tot * zunSlab_tot);
	for(unsigned x = 0; x < xunSlab_tot; x++)
	    for(unsigned y = 0; y < yunSlab_tot; y++)
	      for(unsigned z = 0; z < zunSlab_tot; z++)
		for(int d = 0; d < 3; d++){
		  if(_simulation.isShearRate() && _simulation.getDomain()->getShearWidth() == 0.0)		
			dode->collCommAppendDouble(2*this->_currentDirectedVelSlab[d][cid][x][y][z] - this->_universalDirectedVelocitySlab[d][cid][x][y][z]);
		  else
			dode->collCommAppendDouble(this->_universalDirectedVelocitySlab[d][cid][x][y][z]);
		}
	dode->collCommBroadcast();
	for(unsigned x = 0; x < xunSlab_tot; x++)
	    for(unsigned y = 0; y < yunSlab_tot; y++)
	      for(unsigned z = 0; z < zunSlab_tot; z++)
		for(int d = 0; d < 3; d++)
		  this->_universalDirectedVelocitySlab[d][cid][x][y][z] = dode->collCommGetDouble();
	dode->collCommFinalize();
    }
    }
    if(_simulation.isRecordingStressProfile()){
    for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      dode->collCommInit(4 * xunStress_tot * yunStress_tot * zunStress_tot);
	for(unsigned x = 0; x < xunStress_tot; x++)
	  for(unsigned y = 0; y < yunStress_tot; y++)
	   for(unsigned z = 0; z < zunStress_tot; z++){ 
	    dode->collCommAppendUnsLong(_countedMoleculesStress[cid][x][y][z]);
	    for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocityStress[d][cid][x][y][z]);
	  }
      dode->collCommAllreduceSum();
	for(unsigned x = 0; x < xunStress_tot; x++)
	  for(unsigned y = 0; y < yunStress_tot; y++)
	   for(unsigned z = 0; z < zunStress_tot; z++){
	    _universalCountedMoleculesStress[cid][x][y][z] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++)
		_universalDirectedVelocityStress[d][cid][x][y][z] = dode->collCommGetDouble();
	  }
      dode->collCommFinalize();
      
      if(!dode->getRank())
	{
	  for(unsigned x = 0; x < xunStress_tot; x++)
	    for(unsigned y = 0; y < yunStress_tot; y++)
	      for(unsigned z = 0; z < zunStress_tot; z++){ 
		this->_universalPriorNStress[cid][x][y][z].push_back(this->_universalCountedMoleculesStress[cid][x][y][z]);
		this->_averagedNStress[cid][x][y][z] += this->_universalCountedMoleculesStress[cid][x][y][z];
		for(int d = 0; d < 3; d++){
		  this->_universalPriorVelSumsStress[d][cid][x][y][z].push_back(this->_universalDirectedVelocityStress[d][cid][x][y][z]);
		  this->_averagedVelSumStress[d][cid][x][y][z] += this->_universalDirectedVelocityStress[d][cid][x][y][z];
		  if(this->_universalCountedMoleculesStress[cid][x][y][z] > _molDirVelThreshold)
		    this->_currentDirectedVelStress[d][cid][x][y][z] = this->_universalDirectedVelocityStress[d][cid][x][y][z]/this->_universalCountedMoleculesStress[cid][x][y][z];
		  else
		    this->_currentDirectedVelStress[d][cid][x][y][z] = 0.0;
		}
	      }
	}

	if(!dode->getRank())
	{
	  for(unsigned x = 0; x < xunStress_tot; x++)
	    for(unsigned y = 0; y < yunStress_tot; y++)
	      for(unsigned z = 0; z < zunStress_tot; z++){
		for(int d = 0; d < 3; d++){
                  if(this->_averagedNStress[cid][x][y][z]  > _molDirVelThreshold)
                    this->_universalDirectedVelocityStress[d][cid][x][y][z] = this->_averagedVelSumStress[d][cid][x][y][z]/this->_averagedNStress[cid][x][y][z];
                  else
                    this->_universalDirectedVelocityStress[d][cid][x][y][z] = 0.0;  
		}
		if(this->_universalPriorVelSumsStress[0][cid][x][y][z].size() == _simulation.getDirectedVelocityTime()){
		  for(int d = 0; d < 3; d++){
		    this->_averagedVelSumStress[d][cid][x][y][z] -= this->_universalPriorVelSumsStress[d][cid][x][y][z].front();
		    this->_universalPriorVelSumsStress[d][cid][x][y][z].pop_front();
		  }
		  this->_averagedNStress[cid][x][y][z] -= this->_universalPriorNStress[cid][x][y][z].front();
		  this->_universalPriorNStress[cid][x][y][z].pop_front();
		}
	      }
	}

	dode->collCommInit(3 * xunStress_tot * yunStress_tot * zunStress_tot);
	for(unsigned x = 0; x < xunStress_tot; x++)
	  for(unsigned y = 0; y < yunStress_tot; y++)
	   for(unsigned z = 0; z < zunStress_tot; z++)
		for(int d = 0; d < 3; d++){
		  if(_simulation.isShearRate() && _simulation.getDomain()->getShearWidth() == 0.0)		
			dode->collCommAppendDouble(2*this->_currentDirectedVelStress[d][cid][x][y][z] - this->_universalDirectedVelocityStress[d][cid][x][y][z]);
		  else
			dode->collCommAppendDouble(this->_universalDirectedVelocityStress[d][cid][x][y][z]);
		}
	dode->collCommBroadcast();
	for(unsigned x = 0; x < xunStress_tot; x++)
	  for(unsigned y = 0; y < yunStress_tot; y++)
	   for(unsigned z = 0; z < zunStress_tot; z++)
		for(int d = 0; d < 3; d++)
		  this->_universalDirectedVelocityStress[d][cid][x][y][z] = dode->collCommGetDouble();
	dode->collCommFinalize();
    }
    }
    if(_simulation.isRecordingConfinementProfile()){
    for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      dode->collCommInit(4 * xunConf_tot * yunConf_tot);
	for(unsigned x = 0; x < xunConf_tot; x++)
	  for(unsigned y = 0; y < yunConf_tot; y++){
	    dode->collCommAppendUnsLong(_countedMoleculesConfinement[cid][x][y]);
	    for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocityConfinement[d][cid][x][y]);
	  }
      dode->collCommAllreduceSum();
	for(unsigned x = 0; x < xunConf_tot; x++)
	  for(unsigned y = 0; y < yunConf_tot; y++){
	    _universalCountedMoleculesConfinement[cid][x][y] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++){
		_universalDirectedVelocityConfinement[d][cid][x][y] = dode->collCommGetDouble();
	    }
	  }
      dode->collCommFinalize();
      
      if(!dode->getRank())
	{
	  for(unsigned x = 0; x < xunConf_tot; x++)
	    for(unsigned y = 0; y < yunConf_tot; y++){
		this->_universalPriorNConfinement[cid][x][y].push_back(this->_universalCountedMoleculesConfinement[cid][x][y]);
		this->_averagedNConfinement[cid][x][y] += this->_universalCountedMoleculesConfinement[cid][x][y];
		for(int d = 0; d < 3; d++){
		  this->_universalPriorVelSumsConfinement[d][cid][x][y].push_back(this->_universalDirectedVelocityConfinement[d][cid][x][y]);
		  this->_averagedVelSumConfinement[d][cid][x][y] += this->_universalDirectedVelocityConfinement[d][cid][x][y];
		  if(this->_universalCountedMoleculesConfinement[cid][x][y] > _molDirVelThreshold)
		    this->_currentDirectedVelConfinement[d][cid][x][y] = this->_universalDirectedVelocityConfinement[d][cid][x][y]/this->_universalCountedMoleculesConfinement[cid][x][y];
		  else
		    this->_currentDirectedVelConfinement[d][cid][x][y] = 0.0;
		}
	      }
	}

	if(!dode->getRank())
	{
	  for(unsigned x = 0; x < xunConf_tot; x++)
	    for(unsigned y = 0; y < yunConf_tot; y++){
	      for(int d = 0; d < 3; d++){
                  if(this->_averagedNConfinement[cid][x][y]  > _molDirVelThreshold)
                    this->_universalDirectedVelocityConfinement[d][cid][x][y] = this->_averagedVelSumConfinement[d][cid][x][y]/this->_averagedNConfinement[cid][x][y];
                  else
                    this->_universalDirectedVelocityConfinement[d][cid][x][y] = 0.0;
	      }
		if(this->_universalPriorVelSumsConfinement[0][cid][x][y].size() == _simulation.getDirectedVelocityTime()){
		  for(int d = 0; d < 3; d++){
		    this->_averagedVelSumConfinement[d][cid][x][y] -= this->_universalPriorVelSumsConfinement[d][cid][x][y].front();
		    this->_universalPriorVelSumsConfinement[d][cid][x][y].pop_front();
		  }
		  this->_averagedNConfinement[cid][x][y] -= this->_universalPriorNConfinement[cid][x][y].front();
		  this->_universalPriorNConfinement[cid][x][y].pop_front();
		  
		}
	      }
	}

	dode->collCommInit(3 * xunConf_tot * yunConf_tot);
	for(unsigned x = 0; x < xunConf_tot; x++)
	  for(unsigned y = 0; y < yunConf_tot; y++)
		for(int d = 0; d < 3; d++){
		//Original
		  if(_simulation.isShearRate() && _simulation.getDomain()->getShearWidth() == 0.0)
			dode->collCommAppendDouble(2*this->_currentDirectedVelConfinement[d][cid][x][y] - this->_universalDirectedVelocityConfinement[d][cid][x][y]);
		  else
			dode->collCommAppendDouble(this->_universalDirectedVelocityConfinement[d][cid][x][y]);  	
		}
	dode->collCommBroadcast();
	for(unsigned x = 0; x < xunConf_tot; x++)
	  for(unsigned y = 0; y < yunConf_tot; y++)
		for(int d = 0; d < 3; d++)
		  this->_universalDirectedVelocityConfinement[d][cid][x][y] = dode->collCommGetDouble();
	dode->collCommFinalize();
    }
    }
   }
  }

   // ------- DATA ASSIGNMENT --------	
   // attachment of the _universalDirectedVelocity to each molecule in dependence on the control volume it's part of
    if(_simulation.getSimulationStep() >= _simulation.getInitStatistics()){
    for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
     unsigned cid = molecule->componentid();
     if(cid == cid_moved || cid == cid_fixed){
	for (int d=0; d<3; d++){
		double vAverage = 0.0;
		if(_universalCountedMoleculesAll[cid] > _molDirVelThreshold)
			vAverage = _universalDirectedVelocityAll[d][cid]/_universalCountedMoleculesAll[cid];
		molecule->setDirectedVelocity(d, vAverage); 
		molecule->setDirectedVelocitySlab(d, vAverage);
		molecule->setDirectedVelocityStress(d, vAverage); 
		molecule->setDirectedVelocityConfinement(d, vAverage);
	}
     }else{
	double x = floor(molecule->r(0)) + 0.5;
	double y = floor(molecule->r(1)) + 0.5;
	for (int d=0; d<3; d++){
		molecule->setDirectedVelocity(d, _universalDirectedVelocity[d][cid][x][y]);
		if(_simulation.isMovingAverageSigma2D()){
			if(_simulation.isRecordingSlabProfile())	
				molecule->setDirectedVelocitySlab(d, _universalDirectedVelocity[d][cid][x][y]); 
			if(_simulation.isRecordingStressProfile())
				molecule->setDirectedVelocityStress(d, _universalDirectedVelocity[d][cid][x][y]); 
			if(_simulation.isRecordingConfinementProfile())
				molecule->setDirectedVelocityConfinement(d, _universalDirectedVelocity[d][cid][x][y]);
		}
	}
       if(!_simulation.isMovingAverageSigma2D()){
	if(_simulation.isRecordingSlabProfile()){
		xunSlab = floor(molecule->r(0) / deltaX_slab);
		yunSlab = floor(molecule->r(1) / deltaY_slab);
		zunSlab = floor(molecule->r(2) / deltaZ_slab);
		for (int d=0; d<3; d++)
		molecule->setDirectedVelocitySlab(d, _universalDirectedVelocitySlab[d][cid][xunSlab][yunSlab][zunSlab]); 
	}
	if(_simulation.isRecordingStressProfile()){
		xunStress = floor(molecule->r(0) / deltaX_stress);
		yunStress = floor(molecule->r(1) / deltaY_stress);
		zunStress = floor(molecule->r(2) / deltaZ_stress);
		for (int d=0; d<3; d++)
		molecule->setDirectedVelocityStress(d, _universalDirectedVelocityStress[d][cid][xunStress][yunStress][zunStress]); 
	}
	if(_simulation.isRecordingConfinementProfile()){
		xunConf = floor((molecule->r(0)-_simulation.getDomain()->getConfinementEdge(0)) / deltaX_conf);
		yunConf = floor((molecule->r(1)-(_simulation.getDomain()->get_confinementMidPoint(3)-_simulation.getDomain()->getConfinementEdge(5))) / deltaY_conf);
		if (xunConf >= 0 && xunConf < xunConf_tot && yunConf >= 0 && yunConf < yunConf_tot){
			for (int d=0; d<3; d++){
			molecule->setDirectedVelocityConfinement(d, _universalDirectedVelocityConfinement[d][cid][xunConf][yunConf]);
			}
		}
	}
       }
     }
    }
   }
}

// calculation of the directed velocities: average of neighbour list; size of control volumes: 1*1*z_max OR other --> suitable for more efficient and less parallelized calculation
void VelocityScalingThermostat::calculateDirectedVelocitiesNeighbourList(ParticleContainer *moleculeContainer, DomainDecompBase *dode) {
   double roundX, roundY;
   unsigned CID, xunConf_tot = 0, yunConf_tot = 0, xun = 0, yun = 0, zun = 0;
   unsigned xunConf = 0, yunConf = 0, xunStress = 0, yunStress = 0, zunStress = 0, xunSlab = 0, yunSlab = 0, zunSlab = 0;
   double deltaX_stress = 0.0, deltaY_stress = 0.0, deltaZ_stress = 0.0, deltaX_slab = 0.0, deltaY_slab = 0.0, deltaZ_slab = 0.0, deltaX_conf = 0.0, deltaY_conf = 0.0;
   
   string moved ("moved");
   string fixed ("fixed");
   unsigned cid_moved = _simulation.getDomain()->getCidMovement(moved, _simulation.getDomain()->getNumberOfComponents()) - 1;
   unsigned cid_fixed = _simulation.getDomain()->getCidMovement(fixed, _simulation.getDomain()->getNumberOfComponents()) - 1;
   
   std::map<unsigned,  std::map<double, std::map<double, bool> > > list;
   std::map<unsigned,  std::map<unsigned, std::map<unsigned, std::map<unsigned, bool> > > > list_Slab;
   std::map<unsigned,  std::map<unsigned, std::map<unsigned, std::map<unsigned, bool> > > > list_Stress;
   std::map<unsigned,  std::map<unsigned, std::map<unsigned, bool> > > list_Confinement;
   
   // ------- PREPARATION --------
   // data transfer for determination of the control volumes
   if(_simulation.isRecordingSlabProfile()){
    deltaX_slab = 1/_simulation.getDomain()->getUniversalInvProfileUnitSlab(0);
    deltaY_slab = 1/_simulation.getDomain()->getUniversalInvProfileUnitSlab(1);
    deltaZ_slab = 1/_simulation.getDomain()->getUniversalInvProfileUnitSlab(2);
   }
   if(_simulation.isRecordingStressProfile()){
    deltaX_stress = 1/_simulation.getDomain()->getUniversalInvProfileUnit_Stress(0);
    deltaY_stress = 1/_simulation.getDomain()->getUniversalInvProfileUnit_Stress(1);
    deltaZ_stress = 1/_simulation.getDomain()->getUniversalInvProfileUnit_Stress(2);
   }
   if(_simulation.isRecordingConfinementProfile()){
    xunConf_tot = _simulation.getDomain()->getUniversalNProfileUnitsConfinement(0);
    yunConf_tot = _simulation.getDomain()->getUniversalNProfileUnitsConfinement(1);
    deltaX_conf = 1/_simulation.getDomain()->getUniversalInvProfileUnitConfinement(0);
    deltaY_conf = 1/_simulation.getDomain()->getUniversalInvProfileUnitConfinement(1);
   }
	
   // ------- INITIALIZATION --------	
   // initialize molecule counting; number of molecules of a certain component-id in each control volume
   Molecule *molecule;
   for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
      roundX = floor(molecule->r(0)) + 0.5;
      roundY = floor(molecule->r(1)) + 0.5;
      CID = molecule->componentid();
      
      list[CID][roundX][roundY] = true;
      
      _universalCountedMoleculesAll[CID] = 0;
      for (unsigned d = 0; d < 3; d++){
	_universalDirectedVelocityAll[d][CID] = 0.0;
      }
      _universalCountedMolecules[CID][roundX][roundY] = 0;
      for (unsigned d = 0; d < 3; d++){
	    _universalDirectedVelocity[d][CID][roundX][roundY] = 0.0;
	    if(this->_universalPriorVelSums[d][CID][roundX][roundY].size() == 0){
	      this->_averagedVelSum[d][CID][roundX][roundY] = 0.0;
	      this->_averagedN[CID][roundX][roundY] = 0;
	    }
      }
            
      if((CID != cid_moved && CID != cid_fixed) && _simulation.getSimulationStep() >= _simulation.getInitStatistics()){
       if(_simulation.isRecordingSlabProfile()){
	xun = floor(molecule->r(0) / deltaX_slab);
	yun = floor(molecule->r(1) / deltaY_slab);
	zun = floor(molecule->r(2) / deltaZ_slab);
	
	list_Slab[CID][xun][yun][zun] = true;
	
	_universalCountedMoleculesSlab[CID][xun][yun][zun] = 0;
	for (unsigned d = 0; d < 3; d++){
	    _universalDirectedVelocitySlab[d][CID][xun][yun][zun] = 0.0;
	    if(this->_universalPriorVelSumsSlab[d][CID][xun][yun][zun].size() == 0){
	      this->_averagedVelSumSlab[d][CID][xun][yun][zun] = 0.0;
	      this->_averagedNSlab[CID][xun][yun][zun] = 0;
	    }
	}
       }
       if(_simulation.isRecordingStressProfile()){
	xun = floor(molecule->r(0) / deltaX_stress);
	yun = floor(molecule->r(1) / deltaY_stress);
	zun = floor(molecule->r(2) / deltaZ_stress);
	
	list_Stress[CID][xun][yun][zun] = true;
	
	_universalCountedMoleculesStress[CID][xun][yun][zun] = 0;
	for (unsigned d = 0; d < 3; d++){
	    _universalDirectedVelocityStress[d][CID][xun][yun][zun] = 0.0;
	    if(this->_universalPriorVelSumsStress[d][CID][xun][yun][zun].size() == 0){
	      this->_averagedVelSumStress[d][CID][xun][yun][zun] = 0.0;
	      this->_averagedNStress[CID][xun][yun][zun] = 0;
	    }
	}
       }
       if(_simulation.isRecordingConfinementProfile()){
	xun = floor((molecule->r(0)-_simulation.getDomain()->getConfinementEdge(0)) / deltaX_conf);
	yun = floor((molecule->r(1)-(_simulation.getDomain()->get_confinementMidPoint(3)-_simulation.getDomain()->getConfinementEdge(5))) / deltaY_conf);
	
	list_Confinement[CID][xun][yun] = true;
	
	_universalCountedMoleculesConfinement[CID][xun][yun] = 0;
	for (unsigned d = 0; d < 3; d++){
	    _universalDirectedVelocityConfinement[d][CID][xun][yun] = 0.0;
	    if(this->_universalPriorVelSumsConfinement[d][CID][xun][yun].size() == 0){
	      this->_averagedVelSumConfinement[d][CID][xun][yun] = 0.0;
	      this->_averagedNConfinement[CID][xun][yun] = 0;
	    }
	}
       }
      }
   }
   
 // ------- DATA RECORDING --------  
 if(_simulation.getSimulationStep() >= 1){
   for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
      roundX = floor(molecule->r(0)) + 0.5;
      roundY = floor(molecule->r(1)) + 0.5;
      CID = molecule->componentid();
      if(CID == cid_moved || CID == cid_fixed){
	 _universalCountedMoleculesAll[CID]++;
	for (unsigned d = 0; d < 3; d++)
		_universalDirectedVelocityAll[d][CID] += molecule->v(d);     
      }
      else{
	_universalCountedMolecules[CID][roundX][roundY]++;
	for (unsigned d = 0; d < 3; d++)
		_universalDirectedVelocity[d][CID][roundX][roundY] += molecule->v(d);
      }
      
      if((CID != cid_moved && CID != cid_fixed) && _simulation.getSimulationStep() >= _simulation.getInitStatistics()){
       if(_simulation.isRecordingSlabProfile()){
	xun = floor(molecule->r(0) / deltaX_slab);
	yun = floor(molecule->r(1) / deltaY_slab);
	zun = floor(molecule->r(2) / deltaZ_slab);
	_universalCountedMoleculesSlab[CID][xun][yun][zun]++;
	for (unsigned d = 0; d < 3; d++)
	  _universalDirectedVelocitySlab[d][CID][xun][yun][zun] += molecule->v(d);
       }
       if(_simulation.isRecordingStressProfile()){
	xun = floor(molecule->r(0) / deltaX_stress);
	yun = floor(molecule->r(1) / deltaY_stress);
	zun = floor(molecule->r(2) / deltaZ_stress);
	_universalCountedMoleculesStress[CID][xun][yun][zun]++;
	for (unsigned d = 0; d < 3; d++)
	  _universalDirectedVelocityStress[d][CID][xun][yun][zun] += molecule->v(d);
       }
       if(_simulation.isRecordingConfinementProfile()){
	xun = floor((molecule->r(0)-_simulation.getDomain()->getConfinementEdge(0)) / deltaX_conf);
	yun = floor((molecule->r(1)-(_simulation.getDomain()->get_confinementMidPoint(3)-_simulation.getDomain()->getConfinementEdge(5))) / deltaY_conf);
	_universalCountedMoleculesConfinement[CID][xun][yun]++;
	for (unsigned d = 0; d < 3; d++)
	  _universalDirectedVelocityConfinement[d][CID][xun][yun] += molecule->v(d);
       }
      }
   }
   for(std::map<unsigned,  std::map<double, std::map<double, bool> > >::iterator out = list.begin(); out != list.end(); ++out )
	for(std::map<double, map<double, bool> >::iterator in = out->second.begin(); in != out->second.end(); ++in )
	  for(std::map<double, bool>::iterator inner = in->second.begin(); inner != in->second.end(); ++inner )
	    if(list[out->first][in->first][inner->first] == true){	

		this->_universalPriorN[out->first][in->first][inner->first].push_back(this->_universalCountedMolecules[out->first][in->first][inner->first]);
		this->_averagedN[out->first][in->first][inner->first] += this->_universalCountedMolecules[out->first][in->first][inner->first];
		for(int d = 0; d < 3; d++){
			this->_universalPriorVelSums[d][out->first][in->first][inner->first].push_back(this->_universalDirectedVelocity[d][out->first][in->first][inner->first]);
			this->_averagedVelSum[d][out->first][in->first][inner->first] += this->_universalDirectedVelocity[d][out->first][in->first][inner->first];
			if(this->_universalCountedMolecules[out->first][in->first][inner->first] != 0)
				this->_currentDirectedVel[d][out->first][in->first][inner->first] = this->_universalDirectedVelocity[d][out->first][in->first][inner->first]/this->_universalCountedMolecules[out->first][in->first][inner->first];
			else
				this->_currentDirectedVel[d][out->first][in->first][inner->first] = 0.0;
		}
		for(int d = 0; d < 3; d++){
			if(this->_averagedN[out->first][in->first][inner->first] != 0)
				this->_universalDirectedVelocity[d][out->first][in->first][inner->first] = this->_averagedVelSum[d][out->first][in->first][inner->first]/this->_averagedN[out->first][in->first][inner->first];
			else
				this->_universalDirectedVelocity[d][out->first][in->first][inner->first] = 0.0;  
		}
		if(this->_universalPriorVelSums[0][out->first][in->first][inner->first].size() == _simulation.getDirectedVelocityTime()){
			for(int d = 0; d < 3; d++){
				this->_averagedVelSum[d][out->first][in->first][inner->first] -= this->_universalPriorVelSums[d][out->first][in->first][inner->first].front();
				this->_universalPriorVelSums[d][out->first][in->first][inner->first].pop_front();
			}
			this->_averagedN[out->first][in->first][inner->first] -= this->_universalPriorN[out->first][in->first][inner->first].front();
			this->_universalPriorN[out->first][in->first][inner->first].pop_front();
		}
		
		for(int d = 0; d < 3; d++){
			if(_simulation.isShearRate() && _simulation.getDomain()->getShearWidth() == 0.0)
				_universalDirectedVelocity[d][out->first][in->first][inner->first] = 2*this->_currentDirectedVel[d][out->first][in->first][inner->first] - this->_universalDirectedVelocity[d][out->first][in->first][inner->first];
			else
				_universalDirectedVelocity[d][out->first][in->first][inner->first] = _universalDirectedVelocity[d][out->first][in->first][inner->first];
		}
	    }
     
   if(_simulation.getSimulationStep() > _simulation.getInitStatistics()){
	//Parallelization for SlabProfile, StressProfile and Confinement
	if(_simulation.isRecordingSlabProfile()){
		for(std::map<unsigned,  std::map<unsigned, std::map<unsigned, std::map<unsigned, bool> > > >::iterator out = list_Slab.begin(); out != list_Slab.end(); ++out )
			for(std::map<unsigned, std::map<unsigned, std::map<unsigned, bool> > >::iterator mid = out->second.begin(); mid != out->second.end(); ++mid )
				for(std::map<unsigned, std::map<unsigned, bool> >::iterator in = mid->second.begin(); in != mid->second.end(); ++in )
				  for(std::map<unsigned, bool>::iterator inner = in->second.begin(); inner != in->second.end(); ++inner )
				    if(list_Slab[out->first][mid->first][in->first][inner->first] == true){
		
					this->_universalPriorNSlab[out->first][mid->first][in->first][inner->first].push_back(this->_universalCountedMoleculesSlab[out->first][mid->first][in->first][inner->first]);
					this->_averagedNSlab[out->first][mid->first][in->first][inner->first] += this->_universalCountedMoleculesSlab[out->first][mid->first][in->first][inner->first];
					for(int d = 0; d < 3; d++){
						this->_universalPriorVelSumsSlab[d][out->first][mid->first][in->first][inner->first].push_back(this->_universalDirectedVelocitySlab[d][out->first][mid->first][in->first][inner->first]);
						this->_averagedVelSumSlab[d][out->first][mid->first][in->first][inner->first] += this->_universalDirectedVelocitySlab[d][out->first][mid->first][in->first][inner->first];
						if(this->_universalCountedMoleculesSlab[out->first][mid->first][in->first][inner->first] != 0)
							this->_currentDirectedVelSlab[d][out->first][mid->first][in->first][inner->first] = this->_universalDirectedVelocitySlab[d][out->first][mid->first][in->first][inner->first]/this->_universalCountedMoleculesSlab[out->first][mid->first][in->first][inner->first];
						else
							this->_currentDirectedVelSlab[d][out->first][mid->first][in->first][inner->first] = 0.0;
					}

					for(int d = 0; d < 3; d++){
						if(this->_averagedNSlab[out->first][mid->first][in->first][inner->first] != 0)
							this->_universalDirectedVelocitySlab[d][out->first][mid->first][in->first][inner->first] = this->_averagedVelSumSlab[d][out->first][mid->first][in->first][inner->first]/this->_averagedNSlab[out->first][mid->first][in->first][inner->first];
						else
							this->_universalDirectedVelocitySlab[d][out->first][mid->first][in->first][inner->first] = 0.0;    
					}
					if(this->_universalPriorVelSumsSlab[0][out->first][mid->first][in->first][inner->first].size() == _simulation.getDirectedVelocityTime()){
						for(int d = 0; d < 3; d++){
							this->_averagedVelSumSlab[d][out->first][mid->first][in->first][inner->first] -= this->_universalPriorVelSumsSlab[d][out->first][mid->first][in->first][inner->first].front();
							this->_universalPriorVelSumsSlab[d][out->first][mid->first][in->first][inner->first].pop_front();
						}
						this->_averagedNSlab[out->first][mid->first][in->first][inner->first] -= this->_universalPriorNSlab[out->first][mid->first][in->first][inner->first].front();
						this->_universalPriorNSlab[out->first][mid->first][in->first][inner->first].pop_front();
					}
					
					for(int d = 0; d < 3; d++){
						if(_simulation.isShearRate() && _simulation.getDomain()->getShearWidth() == 0.0)		
							_universalDirectedVelocitySlab[d][out->first][mid->first][in->first][inner->first] = 2*this->_currentDirectedVelSlab[d][out->first][mid->first][in->first][inner->first] - this->_universalDirectedVelocitySlab[d][out->first][mid->first][in->first][inner->first];
						else
							_universalDirectedVelocitySlab[d][out->first][mid->first][in->first][inner->first] = _universalDirectedVelocitySlab[d][out->first][mid->first][in->first][inner->first];	
					}
				    }
	}
	if(_simulation.isRecordingStressProfile()){
		for(std::map<unsigned,  std::map<unsigned, std::map<unsigned, std::map<unsigned, bool> > > >::iterator out = list_Stress.begin(); out != list_Stress.end(); ++out )
			for(std::map<unsigned, std::map<unsigned, std::map<unsigned, bool> > >::iterator mid = out->second.begin(); mid != out->second.end(); ++mid )
				for(std::map<unsigned, std::map<unsigned, bool> >::iterator in = mid->second.begin(); in != mid->second.end(); ++in )
				  for(std::map<unsigned, bool>::iterator inner = in->second.begin(); inner != in->second.end(); ++inner )
				    if(list_Stress[out->first][mid->first][in->first][inner->first] == true){
					    
					this->_universalPriorNStress[out->first][mid->first][in->first][inner->first].push_back(this->_universalCountedMoleculesStress[out->first][mid->first][in->first][inner->first]);
					this->_averagedNStress[out->first][mid->first][in->first][inner->first] += this->_universalCountedMoleculesStress[out->first][mid->first][in->first][inner->first];
					for(int d = 0; d < 3; d++){
						this->_universalPriorVelSumsStress[d][out->first][mid->first][in->first][inner->first].push_back(this->_universalDirectedVelocityStress[d][out->first][mid->first][in->first][inner->first]);
						this->_averagedVelSumStress[d][out->first][mid->first][in->first][inner->first] += this->_universalDirectedVelocityStress[d][out->first][mid->first][in->first][inner->first];
						if(this->_universalCountedMoleculesStress[out->first][mid->first][in->first][inner->first] != 0)
							this->_currentDirectedVelStress[d][out->first][mid->first][in->first][inner->first] = this->_universalDirectedVelocityStress[d][out->first][mid->first][in->first][inner->first]/this->_universalCountedMoleculesStress[out->first][mid->first][in->first][inner->first];
						else
							this->_currentDirectedVelStress[d][out->first][mid->first][in->first][inner->first] = 0.0;
					}

					for(int d = 0; d < 3; d++){
						if(this->_averagedNStress[out->first][mid->first][in->first][inner->first] != 0)
							this->_universalDirectedVelocityStress[d][out->first][mid->first][in->first][inner->first] = this->_averagedVelSumStress[d][out->first][mid->first][in->first][inner->first]/this->_averagedNStress[out->first][mid->first][in->first][inner->first];
						else
							this->_universalDirectedVelocityStress[d][out->first][mid->first][in->first][inner->first] = 0.0;  
					}
					if(this->_universalPriorVelSumsStress[0][out->first][mid->first][in->first][inner->first].size() == _simulation.getDirectedVelocityTime()){
						for(int d = 0; d < 3; d++){
							this->_averagedVelSumStress[d][out->first][mid->first][in->first][inner->first] -= this->_universalPriorVelSumsStress[d][out->first][mid->first][in->first][inner->first].front();
							this->_universalPriorVelSumsStress[d][out->first][mid->first][in->first][inner->first].pop_front();
						}
						this->_averagedNStress[out->first][mid->first][in->first][inner->first] -= this->_universalPriorNStress[out->first][mid->first][in->first][inner->first].front();
						this->_universalPriorNStress[out->first][mid->first][in->first][inner->first].pop_front();
					}
					
					for(int d = 0; d < 3; d++){
						if(_simulation.isShearRate() && _simulation.getDomain()->getShearWidth() == 0.0)		
							_universalDirectedVelocityStress[d][out->first][mid->first][in->first][inner->first] = 2*this->_currentDirectedVelStress[d][out->first][mid->first][in->first][inner->first] - this->_universalDirectedVelocityStress[d][out->first][mid->first][in->first][inner->first];
						else
							_universalDirectedVelocityStress[d][out->first][mid->first][in->first][inner->first] = this->_universalDirectedVelocityStress[d][out->first][mid->first][in->first][inner->first];
					}
				    }
	}
	if(_simulation.isRecordingConfinementProfile()){
		for(std::map<unsigned,  std::map<unsigned, std::map<unsigned, bool> > >::iterator out = list_Confinement.begin(); out != list_Confinement.end(); ++out )
			for(std::map<unsigned, std::map<unsigned, bool> >::iterator in = out->second.begin(); in != out->second.end(); ++in )
			  for(std::map<unsigned, bool>::iterator inner = in->second.begin(); inner != in->second.end(); ++inner )
			    if(list_Confinement[out->first][in->first][inner->first] == true){
	    
				this->_universalPriorNConfinement[out->first][in->first][inner->first].push_back(this->_universalCountedMoleculesConfinement[out->first][in->first][inner->first]);
				this->_averagedNConfinement[out->first][in->first][inner->first] += this->_universalCountedMoleculesConfinement[out->first][in->first][inner->first];
				for(int d = 0; d < 3; d++){
					this->_universalPriorVelSumsConfinement[d][out->first][in->first][inner->first].push_back(this->_universalDirectedVelocityConfinement[d][out->first][in->first][inner->first]);
					this->_averagedVelSumConfinement[d][out->first][in->first][inner->first] += this->_universalDirectedVelocityConfinement[d][out->first][in->first][inner->first];
					if(this->_universalCountedMoleculesConfinement[out->first][in->first][inner->first] != 0)
						this->_currentDirectedVelConfinement[d][out->first][in->first][inner->first] = this->_universalDirectedVelocityConfinement[d][out->first][in->first][inner->first]/this->_universalCountedMoleculesConfinement[out->first][in->first][inner->first];
					else
						this->_currentDirectedVelConfinement[d][out->first][in->first][inner->first] = 0.0;
					}

				for(int d = 0; d < 3; d++){
					if(this->_averagedNConfinement[out->first][in->first][inner->first] != 0)
						this->_universalDirectedVelocityConfinement[d][out->first][in->first][inner->first] = this->_averagedVelSumConfinement[d][out->first][in->first][inner->first]/this->_averagedNConfinement[out->first][in->first][inner->first];
					else
						this->_universalDirectedVelocityConfinement[d][out->first][in->first][inner->first] = 0.0;
				}
				if(this->_universalPriorVelSumsConfinement[0][out->first][in->first][inner->first].size() == _simulation.getDirectedVelocityTime()){
					for(int d = 0; d < 3; d++){
						this->_averagedVelSumConfinement[d][out->first][in->first][inner->first] -= this->_universalPriorVelSumsConfinement[d][out->first][in->first][inner->first].front();
						this->_universalPriorVelSumsConfinement[d][out->first][in->first][inner->first].pop_front();
					}
					this->_averagedNConfinement[out->first][in->first][inner->first] -= this->_universalPriorNConfinement[out->first][in->first][inner->first].front();
					this->_universalPriorNConfinement[out->first][in->first][inner->first].pop_front();
				}

				for(int d = 0; d < 3; d++){
					if(_simulation.isShearRate() && _simulation.getDomain()->getShearWidth() == 0.0)
						_universalDirectedVelocityConfinement[d][out->first][in->first][inner->first] = 2*this->_currentDirectedVelConfinement[d][out->first][in->first][inner->first] - this->_universalDirectedVelocityConfinement[d][out->first][in->first][inner->first];
					else
						_universalDirectedVelocityConfinement[d][out->first][in->first][inner->first] = _universalDirectedVelocityConfinement[d][out->first][in->first][inner->first];  
				}
			  }
	}
   }
 }
 
 // ------- DATA ASSIGNMENT --------
 // attachment of the _universalDirectedVelocity to each molecule in dependence on the control volume part of
 if(_simulation.getSimulationStep() >= _simulation.getInitStatistics()){
    for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
     unsigned cid = molecule->componentid();
     if(cid == cid_moved || cid == cid_fixed){
	for (int d=0; d<3; d++){
		double vAverage = 0.0;
		if(_universalCountedMoleculesAll[cid] > _molDirVelThreshold)
			vAverage = _universalDirectedVelocityAll[d][cid]/_universalCountedMoleculesAll[cid];
		molecule->setDirectedVelocity(d, vAverage); 
		molecule->setDirectedVelocitySlab(d, vAverage);
		molecule->setDirectedVelocityStress(d, vAverage); 
		molecule->setDirectedVelocityConfinement(d, vAverage);
	}
     }else{
	double x = floor(molecule->r(0)) + 0.5;
	double y = floor(molecule->r(1)) + 0.5;
	for (int d=0; d<3; d++)
		molecule->setDirectedVelocity(d, _universalDirectedVelocity[d][cid][x][y]); 
     
	if(_simulation.isRecordingSlabProfile()){
		xunSlab = floor(molecule->r(0) / deltaX_slab);
		yunSlab = floor(molecule->r(1) / deltaY_slab);
		zunSlab = floor(molecule->r(2) / deltaZ_slab);
		for (int d=0; d<3; d++)
			molecule->setDirectedVelocitySlab(d, _universalDirectedVelocitySlab[d][cid][xunSlab][yunSlab][zunSlab]); 
	}
	if(_simulation.isRecordingStressProfile()){
		xunStress = floor(molecule->r(0) / deltaX_stress);
		yunStress = floor(molecule->r(1) / deltaY_stress);
		zunStress = floor(molecule->r(2) / deltaZ_stress);
		for (int d=0; d<3; d++)
			molecule->setDirectedVelocityStress(d, _universalDirectedVelocityStress[d][cid][xunStress][yunStress][zunStress]); 
	}
	if(_simulation.isRecordingConfinementProfile()){
		xunConf = floor((molecule->r(0)-_simulation.getDomain()->getConfinementEdge(0)) / deltaX_conf);
		yunConf = floor((molecule->r(1)-(_simulation.getDomain()->get_confinementMidPoint(3)-_simulation.getDomain()->getConfinementEdge(5))) / deltaY_conf);
		if (xunConf >= 0 && xunConf < xunConf_tot && yunConf >= 0 && yunConf < yunConf_tot){
			for (int d=0; d<3; d++)
				molecule->setDirectedVelocityConfinement(d, _universalDirectedVelocityConfinement[d][cid][xunConf][yunConf]);
		}
	}
     }
    }
   }
}

// calculation of the directed velocities: simple average; size of control volumes: 1*1*z_max OR other
void VelocityScalingThermostat::calculateDirectedVelocitiesSimple(ParticleContainer *moleculeContainer, DomainDecompBase *dode) {
   double xMax = _simulation.getDomain()->getGlobalLength(0);
   double yMax = _simulation.getDomain()->getGlobalLength(1);
   double roundX, roundY;
   unsigned CID, xunStress_tot = 0, yunStress_tot = 0, zunStress_tot = 0, xunSlab_tot = 0, yunSlab_tot = 0, zunSlab_tot = 0, xunConf_tot = 0, yunConf_tot = 0, xun = 0, yun = 0, zun = 0;
   unsigned long xunConf = 0, yunConf = 0, xunStress = 0, yunStress = 0, zunStress = 0, xunSlab = 0, yunSlab = 0, zunSlab = 0;
   double deltaX_stress = 0.0, deltaY_stress = 0.0, deltaZ_stress = 0.0, deltaX_slab = 0.0, deltaY_slab = 0.0, deltaZ_slab = 0.0, deltaX_conf = 0.0, deltaY_conf = 0.0;
   
   string moved ("moved");
   string fixed ("fixed");
   unsigned cid_moved = _simulation.getDomain()->getCidMovement(moved, _simulation.getDomain()->getNumberOfComponents()) - 1;
   unsigned cid_fixed = _simulation.getDomain()->getCidMovement(fixed, _simulation.getDomain()->getNumberOfComponents()) - 1;

   // ------- PREPARATION --------
   // data transfer for determination of the control volumes
   if(_simulation.isRecordingSlabProfile()){
    xunSlab_tot = _simulation.getDomain()->getUniversalNProfileUnitsSlab(0);
    yunSlab_tot = _simulation.getDomain()->getUniversalNProfileUnitsSlab(1);
    zunSlab_tot = _simulation.getDomain()->getUniversalNProfileUnitsSlab(2);
    deltaX_slab = 1/_simulation.getDomain()->getUniversalInvProfileUnitSlab(0);
    deltaY_slab = 1/_simulation.getDomain()->getUniversalInvProfileUnitSlab(1);
    deltaZ_slab = 1/_simulation.getDomain()->getUniversalInvProfileUnitSlab(2);
   }
   if(_simulation.isRecordingStressProfile()){
    xunStress_tot = _simulation.getDomain()->getUniversalNProfileUnits_Stress(0);
    yunStress_tot = _simulation.getDomain()->getUniversalNProfileUnits_Stress(1);
    zunStress_tot = _simulation.getDomain()->getUniversalNProfileUnits_Stress(2);
    deltaX_stress = 1/_simulation.getDomain()->getUniversalInvProfileUnit_Stress(0);
    deltaY_stress = 1/_simulation.getDomain()->getUniversalInvProfileUnit_Stress(1);
    deltaZ_stress = 1/_simulation.getDomain()->getUniversalInvProfileUnit_Stress(2);
   }
   if(_simulation.isRecordingConfinementProfile()){
    xunConf_tot = _simulation.getDomain()->getUniversalNProfileUnitsConfinement(0);
    yunConf_tot = _simulation.getDomain()->getUniversalNProfileUnitsConfinement(1);
    deltaX_conf = 1/_simulation.getDomain()->getUniversalInvProfileUnitConfinement(0);
    deltaY_conf = 1/_simulation.getDomain()->getUniversalInvProfileUnitConfinement(1);
   }

   // ------- INITIALIZATION --------
   // initialize molecule counting; number of molecules of a certain component-id in each control volume
   if(((_simulation.getSimulationStep()-1)%_simulation.getDirectedVelocityTime() == 0)){
    for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      _countedMoleculesAll[cid] = 0;
      for (unsigned d = 0; d < 3; d++){
	_directedVelocityAll[d][cid] = 0.0;
      }
     for (double xCoord = 0.5; xCoord < xMax; xCoord++)
       for (double yCoord = 0.5; yCoord < yMax; yCoord++){
	  _countedMolecules[cid][xCoord][yCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocity[d][cid][xCoord][yCoord] = 0.0;
	  }
       }
    }
   }
  if(!_simulation.isSimpleAverageSigma2D()){
   // same for SlabProfile, StressProfile and Confinement summed up over x timesteps
   for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
    if(_simulation.isRecordingSlabProfile() && (_simulation.getSimulationStep() == 1+_simulation.getSimulationStart() || (_simulation.getSimulationStep() > _simulation.getInitStatistics() && (_simulation.getSimulationStep()-1)%_simulation.getDirectedVelocityTime() == 0))){
     for (unsigned xCoord = 0; xCoord < xunSlab_tot; xCoord++){
       for (unsigned yCoord = 0; yCoord < yunSlab_tot; yCoord++){
	 for (unsigned zCoord = 0; zCoord < zunSlab_tot; zCoord++){
	  _countedMoleculesSlab[cid][xCoord][yCoord][zCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocitySlab[d][cid][xCoord][yCoord][zCoord] = 0.0;
	    _universalDirectedVelocitySlab[d][cid][xCoord][yCoord][zCoord] = 0.0;
	  }
	 }
       }
     }
    }
    if(_simulation.isRecordingStressProfile() && (_simulation.getSimulationStep() == 1+_simulation.getSimulationStart() || (_simulation.getSimulationStep() > _simulation.getInitStatistics() && (_simulation.getSimulationStep()-1)%_simulation.getDirectedVelocityTime() == 0))){
     for (unsigned xCoord = 0; xCoord < xunStress_tot; xCoord++){
       for (unsigned yCoord = 0; yCoord < yunStress_tot; yCoord++){
	 for (unsigned zCoord = 0; zCoord < zunStress_tot; zCoord++){
	  _countedMoleculesStress[cid][xCoord][yCoord][zCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocityStress[d][cid][xCoord][yCoord][zCoord] = 0.0;
	    _universalDirectedVelocityStress[d][cid][xCoord][yCoord][zCoord] = 0.0;
	  }
	 }
       }
     }
    }
    if(_simulation.isRecordingConfinementProfile() && (_simulation.getSimulationStep() == 1+_simulation.getSimulationStart() || (_simulation.getSimulationStep() > _simulation.getInitStatistics() && (_simulation.getSimulationStep()-1)%_simulation.getDirectedVelocityTime() == 0))){
     for (unsigned xCoord = 0; xCoord < xunConf_tot; xCoord++){
       for (unsigned yCoord = 0; yCoord < yunConf_tot; yCoord++){
	  _countedMoleculesConfinement[cid][xCoord][yCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocityConfinement[d][cid][xCoord][yCoord] = 0.0;
	    _universalDirectedVelocityConfinement[d][cid][xCoord][yCoord] = 0.0;
	  }
       }
     }
    }
   }
  }
   
   Molecule *molecule;
   // ------- DATA RECORDING --------

   for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
      roundX = floor(molecule->r(0)) + 0.5;
      roundY = floor(molecule->r(1)) + 0.5;
      CID = molecule->componentid();
      
      if(CID == cid_moved || CID == cid_fixed){
	 _countedMoleculesAll[CID]++;
	for (unsigned d = 0; d < 3; d++)
		_directedVelocityAll[d][CID] += molecule->v(d);     
      }
      else{
	_countedMolecules[CID][roundX][roundY]++;
	for (unsigned d = 0; d < 3; d++)
		_directedVelocity[d][CID][roundX][roundY] += molecule->v(d);
      }
      
     if(!_simulation.isSimpleAverageSigma2D()){ 
      if((CID != cid_moved && CID != cid_fixed) && _simulation.getSimulationStep() >= _simulation.getInitStatistics()){
       if(_simulation.isRecordingSlabProfile()){
	xun = floor(molecule->r(0) / deltaX_slab);
	yun = floor(molecule->r(1) / deltaY_slab);
	zun = floor(molecule->r(2) / deltaZ_slab);
	_countedMoleculesSlab[CID][xun][yun][zun]++;
	for (unsigned d = 0; d < 3; d++)
	  _directedVelocitySlab[d][CID][xun][yun][zun] += molecule->v(d);
       }
       if(_simulation.isRecordingStressProfile()){
	xun = floor(molecule->r(0) / deltaX_stress);
	yun = floor(molecule->r(1) / deltaY_stress);
	zun = floor(molecule->r(2) / deltaZ_stress);
	_countedMoleculesStress[CID][xun][yun][zun]++;
	for (unsigned d = 0; d < 3; d++)
	  _directedVelocityStress[d][CID][xun][yun][zun] += molecule->v(d);
       }
       if(_simulation.isRecordingConfinementProfile()){
	xun = floor((molecule->r(0)-_simulation.getDomain()->getConfinementEdge(0)) / deltaX_conf);
	yun = floor((molecule->r(1)-(_simulation.getDomain()->get_confinementMidPoint(3))) / deltaY_conf);
	_countedMoleculesConfinement[CID][xun][yun]++;
	for (unsigned d = 0; d < 3; d++)
	  _directedVelocityConfinement[d][CID][xun][yun] += molecule->v(d);
       }
      }
     }
   }
   
   // ------- PARALLELIZATION --------
   // for cid_fixed and cid_moved
   if(_simulation.getSimulationStep() >= 1+_simulation.getSimulationStart()  && _simulation.getSimulationStep()%_simulation.getDirectedVelocityTime() == 0){
     for (unsigned cid = 0; cid < _countedMoleculesAll.size(); cid++){
	dode->collCommInit(4);
	dode->collCommAppendUnsLong(_countedMoleculesAll[cid]);
	for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocityAll[d][cid]);
	dode->collCommAllreduceSum();
	_universalCountedMoleculesAll[cid] = dode->collCommGetUnsLong();
	for(int d = 0; d < 3; d++)
		_universalDirectedVelocityAll[d][cid] = dode->collCommGetDouble();
	dode->collCommFinalize();
     }
   }
   
   // carried out every X-th time step (X = _simulation.getDirectedVelocityTime())
   if(_simulation.getSimulationStep() >= 1+_simulation.getSimulationStart() && _simulation.getSimulationStep()%_simulation.getDirectedVelocityTime() == 0){
     for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      dode->collCommInit(4*round(xMax)*round(yMax));
	for(double x = 0.5; x < xMax; x++)
	  for(double y = 0.5; y < yMax; y++){
	    dode->collCommAppendUnsLong(_countedMolecules[cid][x][y]);
	    for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocity[d][cid][x][y]);
	  }
      dode->collCommAllreduceSum();
	for(double x = 0.5; x < xMax; x++)
	  for(double y = 0.5; y < yMax; y++){
	    _universalCountedMolecules[cid][x][y] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++){
		_universalDirectedVelocity[d][cid][x][y] = dode->collCommGetDouble();
	    }
	  }
      dode->collCommFinalize();
      for(double x = 0.5; x < xMax; x++)
	  for(double y = 0.5; y < yMax; y++){
	    for(int d = 0; d < 3; d++){
		if (_universalCountedMolecules[cid][x][y] > _molDirVelThreshold)
		  _universalDirectedVelocity[d][cid][x][y] = _universalDirectedVelocity[d][cid][x][y]/_universalCountedMolecules[cid][x][y];
		else
		  _universalDirectedVelocity[d][cid][x][y] = 0.0;
	    }
	  }
     }
   }
  if(!_simulation.isSimpleAverageSigma2D()){   
   if(_simulation.getSimulationStep() > _simulation.getInitStatistics() && _simulation.getSimulationStep() > _simulation.getSimulationStart() && _simulation.getSimulationStep()%_simulation.getDirectedVelocityTime() == 0){
    //Parallelization for SlabProfile, StressProfile and Confinement
    if(_simulation.isRecordingSlabProfile()){
    for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      dode->collCommInit(4 * xunSlab_tot * yunSlab_tot * zunSlab_tot);
	for(unsigned x = 0; x < xunSlab_tot; x++)
	  for(unsigned y = 0; y < yunSlab_tot; y++)
	   for(unsigned z = 0; z < zunSlab_tot; z++){ 
	    dode->collCommAppendUnsLong(_countedMoleculesSlab[cid][x][y][z]);
	    for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocitySlab[d][cid][x][y][z]);
	  }
      dode->collCommAllreduceSum();
	for(unsigned x = 0; x < xunSlab_tot; x++)
	  for(unsigned y = 0; y < yunSlab_tot; y++)
	   for(unsigned z = 0; z < zunSlab_tot; z++){ 
	    _universalCountedMoleculesSlab[cid][x][y][z] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++){
		_universalDirectedVelocitySlab[d][cid][x][y][z] = dode->collCommGetDouble();
	    }
	  }
      dode->collCommFinalize();
      for(unsigned x = 0; x < xunSlab_tot; x++)
	  for(unsigned y = 0; y < yunSlab_tot; y++)
	   for(unsigned z = 0; z < zunSlab_tot; z++){ 
	    for(int d = 0; d < 3; d++){
		if (_universalCountedMoleculesSlab[cid][x][y][z] > _molDirVelThreshold)
		  _universalDirectedVelocitySlab[d][cid][x][y][z] = _universalDirectedVelocitySlab[d][cid][x][y][z]/_universalCountedMoleculesSlab[cid][x][y][z];
		else
		  _universalDirectedVelocitySlab[d][cid][x][y][z] = 0.0;
	    }
	  }
    }
    }
    if(_simulation.isRecordingStressProfile()){
    for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      dode->collCommInit(4 * xunStress_tot * yunStress_tot * zunStress_tot);
	for(unsigned x = 0; x < xunStress_tot; x++)
	  for(unsigned y = 0; y < yunStress_tot; y++)
	   for(unsigned z = 0; z < zunStress_tot; z++){ 
	    dode->collCommAppendUnsLong(_countedMoleculesStress[cid][x][y][z]);
	    for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocityStress[d][cid][x][y][z]);
	  }
      dode->collCommAllreduceSum();
	for(unsigned x = 0; x < xunStress_tot; x++)
	  for(unsigned y = 0; y < yunStress_tot; y++)
	   for(unsigned z = 0; z < zunStress_tot; z++){
	    _universalCountedMoleculesStress[cid][x][y][z] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++){
		_universalDirectedVelocityStress[d][cid][x][y][z] = dode->collCommGetDouble();
	    }
	  }
      dode->collCommFinalize();
      for(unsigned x = 0; x < xunStress_tot; x++)
	  for(unsigned y = 0; y < yunStress_tot; y++)
	   for(unsigned z = 0; z < zunStress_tot; z++){
	    for(int d = 0; d < 3; d++){
		if (_universalCountedMoleculesStress[cid][x][y][z] > _molDirVelThreshold)
		  _universalDirectedVelocityStress[d][cid][x][y][z] = _universalDirectedVelocityStress[d][cid][x][y][z]/_universalCountedMoleculesStress[cid][x][y][z];
		else
		  _universalDirectedVelocityStress[d][cid][x][y][z] = 0.0;
	    }
	  }
    }
    }
    if(_simulation.isRecordingConfinementProfile()){
    for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      dode->collCommInit(4 * xunConf_tot * yunConf_tot);
	for(unsigned x = 0; x < xunConf_tot; x++)
	  for(unsigned y = 0; y < yunConf_tot; y++){
	    dode->collCommAppendUnsLong(_countedMoleculesConfinement[cid][x][y]);
	    for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocityConfinement[d][cid][x][y]);
	  }
      dode->collCommAllreduceSum();
	for(unsigned x = 0; x < xunConf_tot; x++)
	  for(unsigned y = 0; y < yunConf_tot; y++){
	    _universalCountedMoleculesConfinement[cid][x][y] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++){
		_universalDirectedVelocityConfinement[d][cid][x][y] = dode->collCommGetDouble();
	    }
	  }
      dode->collCommFinalize();
      for(unsigned x = 0; x < xunConf_tot; x++)
	  for(unsigned y = 0; y < yunConf_tot; y++){
	    for(int d = 0; d < 3; d++){
		if (_universalCountedMoleculesConfinement[cid][x][y] > _molDirVelThreshold)
		  _universalDirectedVelocityConfinement[d][cid][x][y] = _universalDirectedVelocityConfinement[d][cid][x][y]/_universalCountedMoleculesConfinement[cid][x][y];
		else
		  _universalDirectedVelocityConfinement[d][cid][x][y] = 0.0;
	    }
	  }
    }
    }
   }
  }
  
   // ------- DATA ASSIGNMENT --------
   // attachment of the _universalDirectedVelocity to each molecule in dependence on the control volume part of
   // carried out every time step since molecules move over time ---> very important!!!
    if(_simulation.getSimulationStep() >= _simulation.getInitStatistics()+_simulation.getDirectedVelocityTime() && _simulation.getSimulationStep() >= _simulation.getSimulationStart()+_simulation.getDirectedVelocityTime()){ 
     for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
      unsigned cid = molecule->componentid();
      double x = floor(molecule->r(0)) + 0.5;
      double y = floor(molecule->r(1)) + 0.5;
      if(cid == cid_moved || cid == cid_fixed){
	for (int d=0; d<3; d++){
		double vAverage = 0.0;
		if(_universalCountedMoleculesAll[cid] > _molDirVelThreshold)
			vAverage = _universalDirectedVelocityAll[d][cid]/_universalCountedMoleculesAll[cid];
		molecule->setDirectedVelocity(d, vAverage); 
		molecule->setDirectedVelocitySlab(d, vAverage);
		molecule->setDirectedVelocityStress(d, vAverage); 
		molecule->setDirectedVelocityConfinement(d, vAverage);
	}
      }else{
       for (int d=0; d<3; d++){
	molecule->setDirectedVelocity(d, _universalDirectedVelocity[d][cid][x][y]); 
	if(_simulation.isSimpleAverageSigma2D()){
	    if(_simulation.isRecordingSlabProfile())	
		molecule->setDirectedVelocitySlab(d, _universalDirectedVelocity[d][cid][x][y]); 
	    if(_simulation.isRecordingStressProfile())
		molecule->setDirectedVelocityStress(d, _universalDirectedVelocity[d][cid][x][y]); 
	    if(_simulation.isRecordingConfinementProfile())
		molecule->setDirectedVelocityConfinement(d, _universalDirectedVelocity[d][cid][x][y]);
	}
       }
       if(!_simulation.isSimpleAverageSigma2D()){
        if(_simulation.isRecordingSlabProfile()){
	 xunSlab = floor(molecule->r(0) / deltaX_slab);
	 yunSlab = floor(molecule->r(1) / deltaY_slab);
	 zunSlab = floor(molecule->r(2) / deltaZ_slab);
	 for (int d=0; d<3; d++)
	  molecule->setDirectedVelocitySlab(d, _universalDirectedVelocitySlab[d][cid][xunSlab][yunSlab][zunSlab]); 
        }
        if(_simulation.isRecordingStressProfile()){
	 xunStress = floor(molecule->r(0) / deltaX_stress);
	 yunStress = floor(molecule->r(1) / deltaY_stress);
	 zunStress = floor(molecule->r(2) / deltaZ_stress);
	 for (int d=0; d<3; d++)
	  molecule->setDirectedVelocityStress(d, _universalDirectedVelocityStress[d][cid][xunStress][yunStress][zunStress]); 
        }
        if(_simulation.isRecordingConfinementProfile()){
	 xunConf = floor((molecule->r(0)-_simulation.getDomain()->getConfinementEdge(0)) / deltaX_conf);
	 yunConf = floor((molecule->r(1)-(_simulation.getDomain()->get_confinementMidPoint(3))) / deltaY_conf);
	 if (xunConf >= 0 && xunConf < xunConf_tot && yunConf >= 0 && yunConf < yunConf_tot){
	  for (int d=0; d<3; d++){
	    molecule->setDirectedVelocityConfinement(d, _universalDirectedVelocityConfinement[d][cid][xunConf][yunConf]);
	  }
	 }
	}
       }
      }
     }
    }
}

// calculation of the directed velocities: simple average; size of control volumes: 1*1*1
void VelocityScalingThermostat::calculateDirectedVelocitiesSimpleSigma3D(ParticleContainer *moleculeContainer, DomainDecompBase *dode) {
   double xMax = _simulation.getDomain()->getGlobalLength(0);
   double yMax = _simulation.getDomain()->getGlobalLength(1);
   double zMax = _simulation.getDomain()->getGlobalLength(2);
   double roundX, roundY, roundZ;
   unsigned CID;
   
   string moved ("moved");
   string fixed ("fixed");
   unsigned cid_moved = _simulation.getDomain()->getCidMovement(moved, _simulation.getDomain()->getNumberOfComponents()) - 1;
   unsigned cid_fixed = _simulation.getDomain()->getCidMovement(fixed, _simulation.getDomain()->getNumberOfComponents()) - 1;

   // ------- INITIALIZATION --------
   // initialize molecule counting; number of molecules of a certain component-id in each control volume
   if(((_simulation.getSimulationStep()-1)%_simulation.getDirectedVelocityTime() == 0)){
    for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
     _countedMoleculesAll[cid] = 0;
     for (unsigned d = 0; d < 3; d++){
	_directedVelocityAll[d][cid] = 0.0;
     }	    
     for (double xCoord = 0.5; xCoord < xMax; xCoord++)
       for (double yCoord = 0.5; yCoord < yMax; yCoord++)
        for (double zCoord = 0.5; zCoord < zMax; zCoord++){
	  _countedMolecules3D[cid][xCoord][yCoord][zCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocity3D[d][cid][xCoord][yCoord][zCoord] = 0.0;
	  }
       }
    }
   }
   
   Molecule *molecule;
   // ------- DATA RECORDING --------

   for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
      roundX = floor(molecule->r(0)) + 0.5;
      roundY = floor(molecule->r(1)) + 0.5;
      roundZ = floor(molecule->r(2)) + 0.5;
      CID = molecule->componentid();
      if(CID == cid_moved || CID == cid_fixed){
	 _countedMoleculesAll[CID]++;
	for (unsigned d = 0; d < 3; d++)
		_directedVelocityAll[d][CID] += molecule->v(d);     
      }else{
	_countedMolecules3D[CID][roundX][roundY][roundZ]++;
	for (unsigned d = 0; d < 3; d++)
		_directedVelocity3D[d][CID][roundX][roundY][roundZ] += molecule->v(d);
      }
   }
     
   // ------- PARALLELIZATION --------
   // for cid_fixed and cid_moved
   if(_simulation.getSimulationStep() >= 1+_simulation.getSimulationStart()  && _simulation.getSimulationStep()%_simulation.getDirectedVelocityTime() == 0){
     for (unsigned cid = 0; cid < _countedMoleculesAll.size(); cid++){
	dode->collCommInit(4);
	dode->collCommAppendUnsLong(_countedMoleculesAll[cid]);
	for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocityAll[d][cid]);
	dode->collCommAllreduceSum();
	_universalCountedMoleculesAll[cid] = dode->collCommGetUnsLong();
	for(int d = 0; d < 3; d++)
		_universalDirectedVelocityAll[d][cid] = dode->collCommGetDouble();
	dode->collCommFinalize();
     }
   }
   
   // carried out every X-th time step (X = _simulation.getDirectedVelocityTime())
   if(_simulation.getSimulationStep() >= 1+_simulation.getSimulationStart() && _simulation.getSimulationStep()%_simulation.getDirectedVelocityTime() == 0){
     for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      dode->collCommInit(4*round(xMax)*round(yMax)*round(zMax));
	for(double x = 0.5; x < xMax; x++)
	  for(double y = 0.5; y < yMax; y++)
	   for(double z = 0.5; z < zMax; z++){
	    dode->collCommAppendUnsLong(_countedMolecules3D[cid][x][y][z]);
	    for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocity3D[d][cid][x][y][z]);
	   }
      dode->collCommAllreduceSum();
	for(double x = 0.5; x < xMax; x++)
	  for(double y = 0.5; y < yMax; y++)
	   for(double z = 0.5; z < zMax; z++){
	    _universalCountedMolecules3D[cid][x][y][z] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++)
		_universalDirectedVelocity3D[d][cid][x][y][z] = dode->collCommGetDouble();
	   }
      dode->collCommFinalize();
      for(double x = 0.5; x < xMax; x++)
	  for(double y = 0.5; y < yMax; y++)
	   for(double z = 0.5; z < zMax; z++){
	    for(int d = 0; d < 3; d++){
		if (_universalCountedMolecules3D[cid][x][y][z] > _molDirVelThreshold)
		  _universalDirectedVelocity3D[d][cid][x][y][z] = _universalDirectedVelocity3D[d][cid][x][y][z]/_universalCountedMolecules3D[cid][x][y][z];
		else
		  _universalDirectedVelocity3D[d][cid][x][y][z] = 0.0;
	    }
	   }
     }
   }
   
   // ------- DATA ASSIGNMENT --------
   // attachment of the _universalDirectedVelocity to each molecule in dependence on the control volume part of
   // carried out every time step since molecules move over time ---> very important!!!
    if(_simulation.getSimulationStep() >= _simulation.getInitStatistics()+_simulation.getDirectedVelocityTime() && _simulation.getSimulationStep() >= _simulation.getSimulationStart()+_simulation.getDirectedVelocityTime()){ 
     for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
      unsigned cid = molecule->componentid();
      double x = floor(molecule->r(0)) + 0.5;
      double y = floor(molecule->r(1)) + 0.5;
      double z = floor(molecule->r(2)) + 0.5;
      if(cid == cid_moved || cid == cid_fixed){
	for (int d=0; d<3; d++){
		double vAverage = 0.0;
		if(_universalCountedMoleculesAll[cid] > _molDirVelThreshold)
			vAverage = _universalDirectedVelocityAll[d][cid]/_universalCountedMoleculesAll[cid];
		molecule->setDirectedVelocity(d, vAverage); 
		molecule->setDirectedVelocitySlab(d, vAverage);
		molecule->setDirectedVelocityStress(d, vAverage); 
		molecule->setDirectedVelocityConfinement(d, vAverage);
	}
      }else{
       for (int d=0; d<3; d++){
	molecule->setDirectedVelocity(d, _universalDirectedVelocity3D[d][cid][x][y][z]); 
	if(_simulation.isRecordingSlabProfile())	
		molecule->setDirectedVelocitySlab(d, _universalDirectedVelocity3D[d][cid][x][y][z]); 
	if(_simulation.isRecordingStressProfile())
		molecule->setDirectedVelocityStress(d, _universalDirectedVelocity3D[d][cid][x][y][z]); 
	if(_simulation.isRecordingConfinementProfile())
		molecule->setDirectedVelocityConfinement(d, _universalDirectedVelocity3D[d][cid][x][y][z]);
      }
      }
    }
   }
}

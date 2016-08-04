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

void VelocityScalingThermostat::apply(ParticleContainer *moleculeContainer, DomainDecompBase *dode) {
	Molecule *molecule;
	
	if(_componentwise ) {
		for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()) {
			int thermostatId;
			double betaTrans = _globalBetaTrans;
			double betaRot = _globalBetaRot;
			int cid = molecule->componentid();
			double alphaTrans = _globalAlphaTrans;
			thermostatId = _simulation.getDomain()->getThermostat(cid);
			betaTrans = _componentBetaTrans[thermostatId];
			betaRot   = _componentBetaRot[thermostatId];
			alphaTrans = _componentAlphaTrans[thermostatId];
			map<int, int> Dim = _simulation.getDomain()->getDim();
			map<int, double*>::iterator vIter;
			if( (vIter = _componentVelocity.find(thermostatId)) != _componentVelocity.end()) {
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
		
		// Calculate directed velocities in control volumes of dimension [1 * 1 * Z_max]
		if(_simulation.getBoolDirVel() == true)
		  calculateDirectedVelocities(moleculeContainer, dode);
		
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
	else {
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
}

// Calculates directed velocities for control volumes of dimension [1 * 1 * Z_max] averaged for each timestep
void VelocityScalingThermostat::calculateDirectedVelocities(ParticleContainer *moleculeContainer, DomainDecompBase *dode) {
   double xMax = _simulation.getDomain()->getGlobalLength(0);
   double yMax = _simulation.getDomain()->getGlobalLength(1);
   double roundX, roundY;
   unsigned CID, xunStress_tot = 0, yunStress_tot = 0, zunStress_tot = 0, xunSlab_tot = 0, yunSlab_tot = 0, zunSlab_tot = 0, xunConf_tot = 0, yunConf_tot = 0, xun = 0, yun = 0, zun = 0;
   unsigned long xunConf, yunConf = 0, xunStress = 0, yunStress = 0, zunStress = 0, xunSlab = 0, yunSlab = 0, zunSlab = 0;
   double deltaX_stress = 0.0, deltaY_stress = 0.0, deltaZ_stress = 0.0, deltaX_slab = 0.0, deltaY_slab = 0.0, deltaZ_slab = 0.0, deltaX_conf = 0.0, deltaY_conf = 0.0;
   
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

   // initialize molecule counting; number of molecules of a certain component-id in each control volume
   if(((_simulation.getSimulationStep()-1)%_simulation.getDirectedVelocityTime() == 0)){
    for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
     for (double xCoord = 0.5; xCoord < xMax; xCoord++)
       for (double yCoord = 0.5; yCoord < yMax; yCoord++){
	  _countedMolecules[cid][xCoord][yCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocity[cid][xCoord][yCoord][d] = 0.0;
	    _universalDirectedVelocity[cid][xCoord][yCoord][d] = 0.0;
	  }
       }
    }
   }
   // same for SlabProfile, StressProfile and Confinement summed up over x timesteps
   for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
    if(_simulation.isRecordingSlabProfile() && (_simulation.getSimulationStep() == 1 || (_simulation.getSimulationStep() > _simulation.getInitStatistics() && (_simulation.getSimulationStep()-1)%_simulation.getDirectedVelocityTime() == 0))){
     for (unsigned xCoord = 0; xCoord < xunSlab_tot; xCoord++){
       for (unsigned yCoord = 0; yCoord < yunSlab_tot; yCoord++){
	 for (unsigned zCoord = 0; zCoord < zunSlab_tot; zCoord++){
	  _countedMoleculesSlab[cid][xCoord][yCoord][zCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocitySlab[cid][xCoord][yCoord][zCoord][d] = 0.0;
	    _universalDirectedVelocitySlab[cid][xCoord][yCoord][zCoord][d] = 0.0;
	  }
	 }
       }
     }
    }
    if(_simulation.isRecordingStressProfile() && (_simulation.getSimulationStep() == 1 || (_simulation.getSimulationStep() > _simulation.getInitStatistics() && (_simulation.getSimulationStep()-1)%_simulation.getDirectedVelocityTime() == 0))){
     for (unsigned xCoord = 0; xCoord < xunStress_tot; xCoord++){
       for (unsigned yCoord = 0; yCoord < yunStress_tot; yCoord++){
	 for (unsigned zCoord = 0; zCoord < zunStress_tot; zCoord++){
	  _countedMoleculesStress[cid][xCoord][yCoord][zCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocityStress[cid][xCoord][yCoord][zCoord][d] = 0.0;
	    _universalDirectedVelocityStress[cid][xCoord][yCoord][zCoord][d] = 0.0;
	  }
	 }
       }
     }
    }
    if(_simulation.isRecordingConfinementProfile() && (_simulation.getSimulationStep() == 1 || (_simulation.getSimulationStep() > _simulation.getInitStatistics() && (_simulation.getSimulationStep()-1)%_simulation.getDirectedVelocityTime() == 0))){
     for (unsigned xCoord = 0; xCoord < xunConf_tot; xCoord++){
       for (unsigned yCoord = 0; yCoord < yunConf_tot; yCoord++){
	  _countedMoleculesConfinement[cid][xCoord][yCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocityConfinement[cid][xCoord][yCoord][d] = 0.0;
	    _universalDirectedVelocityConfinement[cid][xCoord][yCoord][d] = 0.0;
	  }
       }
     }
    }
   }
   
   Molecule *molecule;

   for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
      roundX = floor(molecule->r(0)) + 0.5;
      roundY = floor(molecule->r(1)) + 0.5;
      CID = molecule->componentid();
      _countedMolecules[CID][roundX][roundY]++;
      for (unsigned d = 0; d < 3; d++)
	_directedVelocity[CID][roundX][roundY][d] += molecule->v(d);
      
      if(_simulation.getSimulationStep() >= _simulation.getInitStatistics()){
       if(_simulation.isRecordingSlabProfile()){
	xun = floor(molecule->r(0) / deltaX_slab);
	yun = floor(molecule->r(1) / deltaY_slab);
	zun = floor(molecule->r(2) / deltaZ_slab);
	_countedMoleculesSlab[CID][xun][yun][zun]++;
	for (unsigned d = 0; d < 3; d++)
	  _directedVelocitySlab[CID][xun][yun][zun][d] += molecule->v(d);
       }
       if(_simulation.isRecordingStressProfile()){
	xun = floor(molecule->r(0) / deltaX_stress);
	yun = floor(molecule->r(1) / deltaY_stress);
	zun = floor(molecule->r(2) / deltaZ_stress);
	_countedMoleculesStress[CID][xun][yun][zun]++;
	for (unsigned d = 0; d < 3; d++)
	  _directedVelocityStress[CID][xun][yun][zun][d] += molecule->v(d);
       }
       if(_simulation.isRecordingConfinementProfile()){
	xun = floor((molecule->r(0)-_simulation.getDomain()->getConfinementEdge(0)) / deltaX_conf);
	yun = floor((molecule->r(1)-(_simulation.getDomain()->get_confinementMidPoint(3))) / deltaY_conf);
	_countedMoleculesConfinement[CID][xun][yun]++;
	for (unsigned d = 0; d < 3; d++)
	  _directedVelocityConfinement[CID][xun][yun][d] += molecule->v(d);
       }
      }
   }
     
   //Parallelization
   if(_simulation.getSimulationStep() >= 1 && _simulation.getSimulationStep()%_simulation.getDirectedVelocityTime() == 0){
     for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      dode->collCommInit(4*round(xMax)*round(yMax));
	for(double x = 0.5; x < xMax; x++)
	  for(double y = 0.5; y < yMax; y++){
	    dode->collCommAppendUnsLong(_countedMolecules[cid][x][y]);
	    for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocity[cid][x][y][d]);
	  }
      dode->collCommAllreduceSum();
	for(double x = 0.5; x < xMax; x++)
	  for(double y = 0.5; y < yMax; y++){
	    _universalCountedMolecules[cid][x][y] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++){
		_universalDirectedVelocity[cid][x][y][d] = dode->collCommGetDouble();
		if (_universalCountedMolecules[cid][x][y] != 0)
		  _universalDirectedVelocity[cid][x][y][d] = _universalDirectedVelocity[cid][x][y][d]/_universalCountedMolecules[cid][x][y];
		else
		  _universalDirectedVelocity[cid][x][y][d] = 0.0;
		_universalDirectedVelocity[cid][x][y][d] = _universalDirectedVelocity[cid][x][y][d]/_simulation.getDirectedVelocityTime();
	    }
	  }
      dode->collCommFinalize();
     }
   }
     
   if(_simulation.getSimulationStep() > _simulation.getInitStatistics() && _simulation.getSimulationStep()%_simulation.getDirectedVelocityTime() == 0){
    //Parallelization for SlabProfile, StressProfile and Confinement
    if(_simulation.isRecordingSlabProfile()){
    for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      dode->collCommInit(4 * xunSlab_tot * yunSlab_tot * zunSlab_tot);
	for(unsigned x = 0; x < xunSlab_tot; x++)
	  for(unsigned y = 0; y < yunSlab_tot; y++)
	   for(unsigned z = 0; z < zunSlab_tot; z++){ 
	    dode->collCommAppendUnsLong(_countedMoleculesSlab[cid][x][y][z]);
	    for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocitySlab[cid][x][y][z][d]);
	  }
      dode->collCommAllreduceSum();
	for(unsigned x = 0; x < xunSlab_tot; x++)
	  for(unsigned y = 0; y < yunSlab_tot; y++)
	   for(unsigned z = 0; z < zunSlab_tot; z++){ 
	    _universalCountedMoleculesSlab[cid][x][y][z] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++){
		_universalDirectedVelocitySlab[cid][x][y][z][d] = dode->collCommGetDouble();
		if (_universalCountedMoleculesSlab[cid][x][y][z] != 0)
		  _universalDirectedVelocitySlab[cid][x][y][z][d] = _universalDirectedVelocitySlab[cid][x][y][z][d]/_universalCountedMoleculesSlab[cid][x][y][z];
		else
		  _universalDirectedVelocitySlab[cid][x][y][z][d] = 0.0;
		_universalDirectedVelocitySlab[cid][x][y][z][d] = _universalDirectedVelocitySlab[cid][x][y][z][d]/_simulation.getDirectedVelocityTime();
	    }
	  }
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
		dode->collCommAppendDouble(_directedVelocityStress[cid][x][y][z][d]);
	  }
      dode->collCommAllreduceSum();
	for(unsigned x = 0; x < xunStress_tot; x++)
	  for(unsigned y = 0; y < yunStress_tot; y++)
	   for(unsigned z = 0; z < zunStress_tot; z++){
	    _universalCountedMoleculesStress[cid][x][y][z] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++){
		_universalDirectedVelocityStress[cid][x][y][z][d] = dode->collCommGetDouble();
		if (_universalCountedMoleculesStress[cid][x][y][z] != 0)
		  _universalDirectedVelocityStress[cid][x][y][z][d] = _universalDirectedVelocityStress[cid][x][y][z][d]/_universalCountedMoleculesStress[cid][x][y][z];
		else
		  _universalDirectedVelocityStress[cid][x][y][z][d] = 0.0;
		_universalDirectedVelocityStress[cid][x][y][z][d] = _universalDirectedVelocityStress[cid][x][y][z][d]/_simulation.getDirectedVelocityTime();
	    }
	  }
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
		dode->collCommAppendDouble(_directedVelocityConfinement[cid][x][y][d]);
	  }
      dode->collCommAllreduceSum();
	for(unsigned x = 0; x < xunConf_tot; x++)
	  for(unsigned y = 0; y < yunConf_tot; y++){
	    _universalCountedMoleculesConfinement[cid][x][y] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++){
		_universalDirectedVelocityConfinement[cid][x][y][d] = dode->collCommGetDouble();
		if (_universalCountedMoleculesConfinement[cid][x][y] != 0)
		  _universalDirectedVelocityConfinement[cid][x][y][d] = _universalDirectedVelocityConfinement[cid][x][y][d]/_universalCountedMoleculesConfinement[cid][x][y];
		else
		  _universalDirectedVelocityConfinement[cid][x][y][d] = 0.0;
		_universalDirectedVelocityConfinement[cid][x][y][d] = _universalDirectedVelocityConfinement[cid][x][y][d]/_simulation.getDirectedVelocityTime();
	    }
	  }
      dode->collCommFinalize();
    }
    }
   }
   
   
   // attachment of the _universalDirectedVelocity to each molecule in dependence on the control volume part of
    if(_simulation.getSimulationStep() >= _simulation.getInitStatistics() && _simulation.getSimulationStep()%_simulation.getDirectedVelocityTime() == 0){ 
    for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
     unsigned cid = molecule->componentid();
     CID = molecule->componentid();
     if(_simulation.getSimulationStep() >= 1 && _simulation.getSimulationStep()%_simulation.getDirectedVelocityTime() == 0){
      double x = floor(molecule->r(0)) + 0.5;
      double y = floor(molecule->r(1)) + 0.5;
      
      for (int d=0; d<3; d++){
	molecule->setDirectedVelocity(d, _universalDirectedVelocity[cid][x][y][d]); 
      }
     }
     
      if(_simulation.isRecordingSlabProfile()){
	xunSlab = floor(molecule->r(0) / deltaX_slab);
	yunSlab = floor(molecule->r(1) / deltaY_slab);
	zunSlab = floor(molecule->r(2) / deltaZ_slab);
	for (int d=0; d<3; d++)
	  molecule->setDirectedVelocitySlab(d, _universalDirectedVelocitySlab[cid][xunSlab][yunSlab][zunSlab][d]); 
      }
      if(_simulation.isRecordingStressProfile()){
	xunStress = floor(molecule->r(0) / deltaX_stress);
	yunStress = floor(molecule->r(1) / deltaY_stress);
	zunStress = floor(molecule->r(2) / deltaZ_stress);
	for (int d=0; d<3; d++)
	  molecule->setDirectedVelocityStress(d, _universalDirectedVelocityStress[cid][xunStress][yunStress][zunStress][d]); 
      }
      if(_simulation.isRecordingConfinementProfile()){
	xunConf = floor((molecule->r(0)-_simulation.getDomain()->getConfinementEdge(0)) / deltaX_conf);
	yunConf = floor((molecule->r(1)-(_simulation.getDomain()->get_confinementMidPoint(3))) / deltaY_conf);
	if (xunConf >= 0 && xunConf < xunConf_tot && yunConf >= 0 && yunConf < yunConf_tot){
	  for (int d=0; d<3; d++){
	    molecule->setDirectedVelocityConfinement(d, _universalDirectedVelocityConfinement[cid][xunConf][yunConf][d]);
	  }
	}
      }
     }
   }
}

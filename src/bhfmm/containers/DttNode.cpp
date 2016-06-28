#include "DttNode.h"
#include "Simulation.h"
#include "Domain.h"
#include "UniformPseudoParticleContainer.h"
#include "utils/Logger.h"
#include "bhfmm/utils/RotationParameterLookUp.h"
#include <stdlib.h>

static const bool debug = false;

dtt::DttNode::DttNode(ParticleCell& particles,int threshold,double ctr[3],double domLen[3],
		int order ,int depth, bool srcOnly) :
		_mpCell(order), _occ(true),_leafParticles(), _threshold(threshold), _order(order),
		_splitable(false), _depth(depth), _srcOnly(srcOnly){
	for(unsigned int i=0; i<3; i++){
		_ctr[i] = ctr[i];
		_domLen[i] = domLen[i];			
	}
	if(particles.getMoleculeCount()==0){
		_occ = false;
		return;
	}
	double scalar = _domLen[0]*_domLen[0]+_domLen[1]*_domLen[1]+_domLen[2]*_domLen[2];
	double radius = sqrt(scalar)*0.5;
	double c[3];
	for(int k=0; k<3; k++){
		c[k]=_ctr[k];
	}		
	
	bhfmm::Vector3<double> v3Ctr(c);
	if(!_srcOnly){
		_mpCell.local.setCenter(v3Ctr);
		_mpCell.local.setRadius(radius);
	}		
	_mpCell.multipole.setCenter(v3Ctr);
	_mpCell.multipole.setRadius(radius);
	
	int pCount = particles.getMoleculeCount();
	if(threshold == 0){ pCount *= -1; }

	if(_depth <= 0 && _threshold >= pCount){
		_splitable = false;

		int currentParticleCount = particles.getMoleculeCount();

		// loop over all particles in the cell
		for (int i = 0; i < currentParticleCount; i++) {
			Molecule& molecule1 = particles.moleculesAt(i);
			_leafParticles.addParticle(&molecule1);
		} // current particle closed

	}else{ 
		_splitable = true;

		std::vector<ParticleCell> particleContainer(8,ParticleCell());
		divideParticles(particles, particleContainer);
		
		double child_ctr[3];		
		double child_domLen[3]; 
		child_domLen[0] = _domLen[0]/2;
		child_domLen[1] = _domLen[1]/2;
		child_domLen[2] = _domLen[2]/2;
				
		for(int i=0; i<8; i++){
			child_ctr[0] = i%2 == 0 ? _ctr[0] - child_domLen[0]/2 : _ctr[0] + child_domLen[0]/2;
			child_ctr[1] = i%4 > 1 	? _ctr[1] - child_domLen[1]/2 : _ctr[1] + child_domLen[1]/2;
			child_ctr[2] = i > 3 ? _ctr[2] - child_domLen[2]/2 : _ctr[2] + child_domLen[2]/2;
	
			_children.push_back( new dtt::DttNode(particleContainer[i],_threshold,
				child_ctr,child_domLen,_order,_depth-1,_srcOnly) );
		}
	}
}

bool dtt::DttNode::upwardPass(){
	if(!_occ){ return false; }

	if(!_splitable){
		int currentParticleCount = _leafParticles.getMoleculeCount();

		int Occupied = 0;

		// loop over all particles in the cell
		for (int i = 0; i < currentParticleCount; i++) {
			++Occupied;
			Molecule& molecule1 = _leafParticles.moleculesAt(i);
			
			int ni= molecule1.numCharges();

			for(int j=0; j<ni; j++){
				const std::array<double,3> dii = molecule1.charge_d(j);
				const Charge& chargei=static_cast<const Charge&> (molecule1.component()->charge(j));
				double dr[3];

				for(int k=0; k<3; k++){
					dr[k]=molecule1.r(k)+dii[k];
				}	// for k closed

				bhfmm::Vector3<double> site_pos_vec3(dr);
				_mpCell.multipole.addSource(site_pos_vec3, chargei.q());				
			}// for j closed
		} // current particle closed

		_mpCell.occ = 1;
		return true;
	}

	for(int i = 0; i < 8; i++){
		if(_children[i]->_occ){
			_children[i]->upwardPass();
			_mpCell.multipole.addMultipoleParticle(_children[i]->_mpCell.multipole);
		}
	}	
	return true;
}

void dtt::DttNode::downwardPass(){
	if(_splitable){
		for(unsigned int i = 0; i < 8 ; i++){
			if(!_children[i]->_occ) continue;
			_mpCell.local.actOnLocalParticle(_children[i]->_mpCell.local);
			_children[i]->downwardPass();
		}
	}else{
		int currentParticleCount = _leafParticles.getMoleculeCount();
		bhfmm::SolidHarmonicsExpansion leLocal(_order);
		double u = 0;
		double uSum = 0.0;
		double f[3] = {0.0, 0.0, 0.0};
		bhfmm::Vector3<double>f_vec3;
		double virialSum=0.0;
		double P_xxSum=0.0;
		double P_yySum=0.0;
		double P_zzSum=0.0;

		for (int i = 0; i < currentParticleCount; i++) {
			Molecule& molecule1 = _leafParticles.moleculesAt(i);
			int ni = molecule1.numCharges();

			for(int j=0; j<ni; j++){
				const std::array<double,3> dii = molecule1.charge_d(j);
				const Charge& chargei=static_cast<const Charge&> (molecule1.component()->charge(j));
				bhfmm::Vector3<double> dr;

				for(int k=0; k<3; k++){
					dr[k]=molecule1.r(k)+dii[k];
				}// for k closed

				_mpCell.local.actOnTarget(dr,chargei.q(),u,f_vec3);
				f[0] = f_vec3[0];
				f[1] = f_vec3[1];
				f[2] = f_vec3[2];
	
				double virial = 0.0;
				for(int l=0; l<3; l++){
					virial +=-f[l]*dr[l];
				}
				P_xxSum +=0.5*-f[0]*dr[0];
				P_yySum +=0.5*-f[1]*dr[1];
				P_zzSum +=0.5*-f[2]*dr[2];
				molecule1.Fchargeadd(j, f);
				uSum +=0.5*chargei.q()*u;
				virialSum +=0.5*virial;
			}
		}	
	
//		_domain->addLocalUpot(uSum);
//		_domain->addLocalVirial(virialSum);
//		_domain->addLocalP_xx(P_xxSum);
//		_domain->addLocalP_yy(P_yySum);
//		_domain->addLocalP_zz(P_zzSum);
	}
}

std::vector<ParticleCell> dtt::DttNode::getLeafParticleCells(){
	std::vector<ParticleCell> retval(0);
	if(!_splitable){
		retval.push_back(_leafParticles);	
	}else{
		std::vector<ParticleCell> lower;
		for(unsigned int i = 0; i < 8; i++){
			if(!_children[i]-> _occ) continue;
			lower = _children[i]->getLeafParticleCells();
			retval.insert(retval.end(),lower.begin(),lower.end());
		}
	}	
	return retval;
}

void dtt::DttNode::p2p(bhfmm::VectorizedChargeP2PCellProcessor * v_c_p2p_c_p){
	// TODO
//	_leafParticles.convertAoSToSoACharge();
	v_c_p2p_c_p->preprocessCell(_leafParticles);
	v_c_p2p_c_p->processCell(_leafParticles);
	v_c_p2p_c_p->postprocessCell(_leafParticles);
//	_leafParticles.convertSoAToAoSCharge();
}
		
void dtt::DttNode::p2p(std::vector<ParticleCell> leafParticlesFar, 
	bhfmm::VectorizedChargeP2PCellProcessor * v_c_p2p_c_p, bhfmm::Vector3<double> shift){
	std::vector<double> _shift;
	for(int i=0; i<3; i++){
		_shift.push_back(shift[i]);
  	} 
	if(!_splitable){
		//TODO only convert once
//		_leafParticles.convertAoSToSoACharge();
		v_c_p2p_c_p->preprocessCell(_leafParticles);
		for(unsigned int i = 0; i < leafParticlesFar.size(); i++){
			//TODO only convert once!! 
//			leafParticlesFar[i].convertAoSToSoACharge(_shift);
			v_c_p2p_c_p->preprocessCell(leafParticlesFar[i]);
			v_c_p2p_c_p->processCellPair(_leafParticles,leafParticlesFar[i]);
			v_c_p2p_c_p->postprocessCell(leafParticlesFar[i]);
			//leafParticlesFar[i].convertSoAToAoSCharge(_shift);
		}
		v_c_p2p_c_p->postprocessCell(_leafParticles);
//		_leafParticles.convertSoAToAoSCharge();
	}else{
		for(unsigned int i = 0; i < 8; i++){
			_children[i]->p2p(leafParticlesFar, v_c_p2p_c_p,shift);
		}
	}
}

void dtt::DttNode::m2l( bhfmm::SHMultipoleParticle& multipole, bhfmm::Vector3<double> periodicShift){
	_mpCell.local.addMultipoleParticle(multipole, periodicShift);
}

void dtt::DttNode::divideParticles( ParticleCell& particles, std::vector<ParticleCell>& cell_container ){
  	int child;
	int currentParticleCount = particles.getMoleculeCount();
	// loop over all particles in the cell
	for (int i = 0; i < currentParticleCount; i++) {
		Molecule& molecule1 = particles.moleculesAt(i);
		int c[3];
		c[0] = molecule1.r(0) < _ctr[0] ? 0 : 1;
		c[1] = molecule1.r(1) < _ctr[1] ? 1 : 0;
		c[2] = molecule1.r(2) < _ctr[2] ? 1 : 0;
	
		child = c[0] + 2*c[1] + 4*c[2];
		cell_container[child].addParticle(&molecule1);
	}
}	

int dtt::DttNode::getMaxDepth(){
	int max = 0;
	int act = 0;
	if(!_splitable) return 0;
	for(unsigned int i=0; i<8; i++){
		if(_children[i]->_occ)
			act = _children[i]->getMaxDepth();
		if(max < act){
			max = act;
		}
	}
	return max+1;
}

void dtt::DttNode::printSplitable(bool print){
	if(!print){ return; }
	if(!_splitable){
		if(!_occ)std::cout << "EMPTY";
		return;
	}
	for(int i = 0; i < 8; i++){
		if(!_children[i]->_occ){
			std::cout << "EMPTY";	 
			continue;
		}
		std::cout<<i<<"[";
		_children[i]->printSplitable(print);
		std::cout<<"]";	
	}
}

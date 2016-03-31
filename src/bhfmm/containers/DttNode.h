#ifndef DTTNODE_H_
#define DTTNODE_H_
#include "particleContainer/ParticleCell.h"
#include "PseudoParticleContainer.h"
#include <vector>

#include "bhfmm/cellProcessors/VectorizedChargeP2PCellProcessor.h"

class DttNodeTest;

namespace dtt{
	class DttNode;
}

class dtt::DttNode{
	friend class ::DttNodeTest;

	public:
		DttNode(int o):_mpCell(o){};
		DttNode(ParticleCell particles,int threshold,double ctr[3],double domLen[3], 
			int order,int depth=0, bool srcOnly=false);
		
		~DttNode() { for(unsigned int i=0;i<_children.size();i++){delete _children[i];} }			
		
		double _ctr[3], _domLen[3];
		bhfmm::MpCell _mpCell;
		bool _occ;	
		ParticleCell _leafParticles;
		
		bool get_children(std::vector<DttNode*> & ch)
			{ ch = _children; return _splitable; };
		
		bool upwardPass();
		void downwardPass();
		void p2p(bhfmm::VectorizedChargeP2PCellProcessor * v_c_p2p_c_p);
		void p2p(std::vector<ParticleCell> leafParticlesFar,
			bhfmm::VectorizedChargeP2PCellProcessor * v_c_p2p_c_p, bhfmm::Vector3<double> shift);
		void m2l(bhfmm::SHMultipoleParticle& multipole, bhfmm::Vector3<double> periodicShift);
		
		std::vector<ParticleCell> getLeafParticleCells();	
		int getMaxDepth();
		void printSplitable(bool print);
	
	private:
		double _threshold;
		int _order;
    		bool _splitable;
		std::vector<DttNode*> _children; 
		int _depth;
		bool _srcOnly;
		void initTree(ParticleCell particles);		
    void divideParticles(ParticleCell particles, std::vector<ParticleCell>& cellContainer);
};

#endif /* DTTNODE_H_ */

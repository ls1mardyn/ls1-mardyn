#ifndef DTTNODE_H_
#define DTTNODE_H_

#include "particleContainer/ParticleCell.h"
#include "PseudoParticleContainer.h"
#include "bhfmm/cellProcessors/VectorizedChargeP2PCellProcessor.h"
#include "bhfmm/utils/Vector3.h"

#include <vector>
#include <cassert>

class DttNodeTest;

namespace bhfmm {
class DttNode;
}

class bhfmm::DttNode {
	friend class ::DttNodeTest;

public:
	DttNode(int o) :
			_mpCell(o) {
	}

	DttNode(const std::vector<Molecule *>& particles, int threshold, Vector3<double> ctr,
			Vector3<double> domLen, int order, int depth = 0, bool srcOnly = false);

	~DttNode() {
		for (unsigned int i = 0; i < _children.size(); i++) {
			delete _children[i];
		}
	}

	bool get_children(std::vector<DttNode*> & ch) {
		ch = _children;
		return not _isLeafNode;
	}

	void upwardPass();
	void downwardPass();
	void p2p(VectorizedChargeP2PCellProcessor * v_c_p2p_c_p);
	void p2p(std::vector<ParticleCell> leafParticlesFar,
			VectorizedChargeP2PCellProcessor * v_c_p2p_c_p,
			Vector3<double> shift);
	void m2l(const SHMultipoleParticle& multipole,
			Vector3<double> periodicShift);

	std::vector<ParticleCell> getLeafParticleCells();
	int getMaxDepth() const;
	void printSplitable(bool print) const;

	bool isEmpty() const {
		return _mpCell.occ == 0;
	}

	bool isOccupied() const {
		return not isEmpty();
	}

	Vector3<double> getCenter() const {
		return _ctr;
	}
	Vector3<double> getSize() const {
		return _domLen;
	}
	double getSize(int d) const {
		assert(d < 2 and d >= 0);
		return _domLen[d];
	}
	MpCell& getMpCell() {
		return _mpCell;
	}

private:
	Vector3<double> _ctr, _domLen;
	MpCell _mpCell;
	ParticleCell _leafParticles;

	double _threshold;
	int _order;
	bool _isLeafNode;
	std::vector<DttNode*> _children;
	int _depth;
	bool _srcOnly;
	//void initTree(ParticleCell particles);
	void divideParticles(const std::vector<Molecule *>& particles,
			std::array<std::vector<Molecule *>, 8>& cell_container) const;
};

#endif /* DTTNODE_H_ */

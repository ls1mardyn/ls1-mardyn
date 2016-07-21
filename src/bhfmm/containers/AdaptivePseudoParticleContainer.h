#ifndef ADAPTIVEPSEUDOPARTICLECONTAINER_H_
#define ADAPTIVEPSEUDOPARTICLECONTAINER_H_

#include "PseudoParticleContainer.h"
#include "DttNode.h"
#include <vector>
#include <cmath>
#include <math.h>
#include <stdlib.h>

#include "bhfmm/cellProcessors/VectorizedChargeP2PCellProcessor.h"

namespace bhfmm {

typedef struct {
	dtt::DttNode *target, *source;
	bhfmm::Vector3<double> shift;
} TargetSourceTupel;

class AdaptivePseudoParticleContainer: public PseudoParticleContainer {
public:
	AdaptivePseudoParticleContainer(double domainLength[3], int threshold,
			int orderOfExpansions, bool periodic = true) :
			PseudoParticleContainer(orderOfExpansions), _periodicBC(periodic), _threshold(
					threshold), root(0), halo_node(0), _subdivisionFactor(0) {
		assert(_threshold > 0);
		for (int i = 0; i < 3; i++) {
			_domainLength[i] = domainLength[i];
		}
	}

	AdaptivePseudoParticleContainer(double domainLength[3],
			int orderOfExpansions, double cellLength[3], int subdivisionFactor,
			bool periodic) :
			PseudoParticleContainer(orderOfExpansions), _periodicBC(periodic), _threshold(
					0), root(0), _subdivisionFactor(subdivisionFactor) {
		for (int i = 0; i < 3; i++) {
			_domainLength[i] = domainLength[i];
			_cellLength[i] = cellLength[i];
		}
	}

	~AdaptivePseudoParticleContainer() {
		delete root;
	}
	;

	void clear();
	void build(ParticleContainer* pc);
	void upwardPass(P2MCellProcessor * cp);
	void horizontalPass(VectorizedChargeP2PCellProcessor * cp);
	void downwardPass(L2PCellProcessor *cp);

	void processMultipole(ParticleCell& /*cell*/) {
	}
	void processFarField(ParticleCell& /*cell*/) {
	}
	void processTree() {
	}
	void printTimers() {
	}
	;

private:
	bool _periodicBC;

	VectorizedChargeP2PCellProcessor * vc_p2p_cp;
	ParticleCell _Cells;
	int _threshold;
	std::vector<TargetSourceTupel> stack;
	dtt::DttNode *root, *halo_node;
	double _domainLength[3];
	double _cellLength[3];
	int _subdivisionFactor;

	void buildHaloTrees();
	void work_on_stack();
	void push_on_stack(dtt::DttNode * trg, dtt::DttNode * src,
			Vector3<double> shift, int split);
	void print_stack_op(int op, dtt::DttNode * t = NULL,
			dtt::DttNode * s = NULL);

};
//AdaptivePseudoParticleContainer

}//nemspace bhfmm

#endif //ADAPTIVEPSEUDOPARTICLECONTAINER_H_

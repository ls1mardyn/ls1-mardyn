#ifndef ADAPTIVEPSEUDOPARTICLECONTAINER_H_
#define ADAPTIVEPSEUDOPARTICLECONTAINER_H_

#include "bhfmm/utils/Vector3.h"
#include "PseudoParticleContainer.h"
#include "DttNode.h"
#include "bhfmm/cellProcessors/VectorizedChargeP2PCellProcessor.h"

#include <vector>
#include <cmath>
#include <math.h>
#include <stdlib.h>



namespace bhfmm {

typedef struct {
	DttNode *target, *source;
	Vector3<double> shift;
} TargetSourceTupel;

class AdaptivePseudoParticleContainer: public PseudoParticleContainer {
public:
	AdaptivePseudoParticleContainer(double domainLength[3], int threshold,
			int orderOfExpansions, bool periodic = true) :
			PseudoParticleContainer(orderOfExpansions), _periodicBC(periodic), _threshold(
					threshold), root(0), halo_node(0), _domainLength(
					domainLength), _subdivisionFactor(0) {
		assert(_threshold > 0);
	}

	AdaptivePseudoParticleContainer(double domainLength[3],
			int orderOfExpansions, double cellLength[3], int subdivisionFactor,
			bool periodic) :
			PseudoParticleContainer(orderOfExpansions), _periodicBC(periodic), _threshold(
					0), root(0), _domainLength(domainLength), _cellLength(
					cellLength), _subdivisionFactor(subdivisionFactor) {
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
	std::vector<Molecule *> _particles;
	int _threshold;
	std::vector<TargetSourceTupel> stack;
	DttNode *root, *halo_node;
	Vector3<double> _domainLength, _cellLength;
	int _subdivisionFactor;

	void buildHaloTrees();
	void work_on_stack();
	void push_on_stack(DttNode * trg, DttNode * src,
			Vector3<double> shift, int split);
	void print_stack_op(int op, DttNode * t = NULL,
			DttNode * s = NULL);

};
//AdaptivePseudoParticleContainer

}//nemspace bhfmm

#endif //ADAPTIVEPSEUDOPARTICLECONTAINER_H_

#include "AdaptivePseudoParticleContainer.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"

namespace bhfmm {

const int split_target = 0, split_source = 1, split_internal = 3, NO_WORK = 5,
		M2L = 6, P2P = 7, PUSH = 8, POP = 9;

const double epsilon = 0.001;

const bool debug = false;

void AdaptivePseudoParticleContainer::clear() {
	//TODO FastMultipoleMethod.cpp calls this before computation
	//_Cells.removeAllParticles();
	//delete root;
	//delete halo_node;
}

void AdaptivePseudoParticleContainer::build(ParticleContainer* pc) {
	double lowBound[3] = { 0.0, 0.0, 0.0 };
	double highBound[3] = { _domainLength[0], _domainLength[1], _domainLength[2]};

	ParticleIterator tM;
	for (tM = pc->iterator(); tM.hasNext(); tM.next()) {
		if (tM->inBox(lowBound, highBound)) {
			_particles.push_back(&(*tM));
		}
	}

	Vector3<double> ctr = _domainLength * 0.5;

	if (_threshold == 0) {
		int depth = log2(
				(_domainLength[0] / _cellLength[0]) * _subdivisionFactor);
		root = new DttNode(_particles, _threshold, ctr, _domainLength, _maxOrd,
				depth);
	} else {
		root = new DttNode(_particles, _threshold, ctr, _domainLength, _maxOrd,
				0);
	}

	TargetSourceTupel tst;
	tst.source = root;
	tst.target = root;
	bhfmm::Vector3<double> shift;
	for (int i = 0; i < 3; i++) {
		shift[i] = 0.0;
	}
	tst.shift = shift;
	stack.push_back(tst);

	if (_periodicBC) {
		buildHaloTrees();
	}
}

void AdaptivePseudoParticleContainer::upwardPass(P2MCellProcessor* /*cp*/) {
	root->upwardPass();
}

void AdaptivePseudoParticleContainer::horizontalPass(
		VectorizedChargeP2PCellProcessor* cp) {
	// M2l and P2P
	vc_p2p_cp = cp;
	while (stack.size() > 0) {
		work_on_stack();
	}
}

void AdaptivePseudoParticleContainer::downwardPass(L2PCellProcessor* /*cp*/) {
	// L2L and L2P
	root->downwardPass();
}

void AdaptivePseudoParticleContainer::buildHaloTrees() {
	int root_depth = root->getMaxDepth();
	double ctr[3] = { _domainLength[0] / 2, _domainLength[1] / 2,
			_domainLength[2] / 2 };

	bhfmm::Vector3<double> shift;
	halo_node = new DttNode(_particles, _threshold, ctr, _domainLength,
			_maxOrd, root_depth, true);
	TargetSourceTupel tst;
	tst.target = root;
	tst.source = halo_node;

	for (int x = -1; x <= 1; x++) {
		shift[0] = _domainLength[0] * (double) x;
		for (int y = -1; y <= 1; y++) {
			shift[1] = _domainLength[1] * (double) y;
			for (int z = -1; z <= 1; z++) {
				if (x == 0 && y == 0 && z == 0)
					continue;

				shift[2] = _domainLength[2] * (double) z;
				tst.shift = shift;
				stack.push_back(tst);
			}
		}
	}
}

void AdaptivePseudoParticleContainer::print_stack_op(int op,
		DttNode * /*t*/, DttNode * /*s*/) {
	if (!debug) {
		return;
	}
	switch (op) {
	case NO_WORK:
		std::cout << "		--NO WORK--" << std::endl;
		break;
	case M2L:
		std::cout << "		--M2L--";
		break;
	case P2P:
		std::cout << "		--P2P--" << std::endl;
		break;
	case PUSH:
		std::cout << "    --PUSH--" << std::endl;
		break;
	case POP:
		std::cout << "    --POP--" << std::endl;
		break;
	}
}

void AdaptivePseudoParticleContainer::push_on_stack(DttNode * trg,
		DttNode * src, Vector3<double> shift, int split = split_target) {
	std::vector<DttNode*> *src_chs = new std::vector<DttNode*>();
	std::vector<DttNode*> *trg_chs = new std::vector<DttNode*>();

	TargetSourceTupel tst;
	tst.shift = shift;
	switch (split) {
	case split_target:
		if (trg->get_children(*trg_chs)) {
			tst.source = src;
			for (unsigned int i = 0; i < trg_chs->size(); i++) {
				tst.target = (*trg_chs)[i];
				if (tst.target->isEmpty()) {
					print_stack_op(NO_WORK, tst.target, tst.source);
					continue;
				}
				print_stack_op(PUSH, tst.target, tst.source);
				stack.push_back(tst);
			}
		} else {
			print_stack_op(P2P);
			trg->p2p(src->getLeafParticleCells(), vc_p2p_cp, shift);
		}
		break;
	case split_source:
		if (src->get_children(*src_chs)) {
			tst.target = trg;
			for (unsigned int i = 0; i < src_chs->size(); i++) {
				tst.source = (*src_chs)[i];
				if (tst.source->isEmpty()) {
					print_stack_op(NO_WORK, tst.target, tst.source);
					continue;
				}
				print_stack_op(PUSH, tst.target, tst.source);
				stack.push_back(tst);
			}
		} else {
			print_stack_op(P2P);
			trg->p2p(src->getLeafParticleCells(), vc_p2p_cp, shift);
		}
		break;
	case split_internal:
		if (trg->get_children(*trg_chs)) {
			mardyn_assert(src->get_children(*src_chs));
			for (unsigned int i = 0; i < trg_chs->size(); i++) {
				for (unsigned int j = i; j < trg_chs->size(); j++) {
					tst.target = (*trg_chs)[i];
					tst.source = (*trg_chs)[j];
					if (tst.target->isEmpty() or tst.source->isEmpty()) {
						print_stack_op(NO_WORK, tst.target, tst.source);
						continue;
					}
					print_stack_op(PUSH, tst.target, tst.source);
					stack.push_back(tst);
				}
			}
		} else {
			print_stack_op(P2P);
			trg->p2p(vc_p2p_cp);
		}
		break;
	default:
		break;
	}
	delete src_chs;
	delete trg_chs;
}

void AdaptivePseudoParticleContainer::work_on_stack() {
	TargetSourceTupel tst = stack.back();
	stack.pop_back();
	print_stack_op(POP, tst.target, tst.source);

	double comp_r = tst.target->getMpCell().local.getRadius()
			- tst.source->getMpCell().multipole.getRadius();

	//TODO: hardcoded cells of equal size
	bool sameSize = std::abs(comp_r) <= epsilon;

	double dist = (tst.target->getCenter() - (tst.source->getCenter() + tst.shift)).L2Norm();

	if (dist > epsilon) {
//		TODO: assuming cubic cells
		if (((tst.target->getSize(0) + tst.source->getSize(0)) / dist < 1 + epsilon) and sameSize) {
			print_stack_op(M2L);
			tst.target->m2l(tst.source->getMpCell().multipole, tst.shift);
			if (tst.shift[0] == 0.0 && tst.shift[1] == 0.0
					&& tst.shift[2] == 0.0) {
				tst.source->m2l(tst.target->getMpCell().multipole, tst.shift);
			}
		} else if (comp_r < -epsilon) {
			push_on_stack(tst.target, tst.source, tst.shift, split_source);
		} else {
			push_on_stack(tst.target, tst.source, tst.shift, split_target);
		}
	} else {
		push_on_stack(tst.target, tst.source, tst.shift, split_internal);
	}
}

} //namespace bhfmm


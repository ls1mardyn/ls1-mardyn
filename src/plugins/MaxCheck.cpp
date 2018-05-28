/*
 * MaxCheck.cpp
 *
 *  Created on: 28.05.2018
 *      Author: mheinen
 */

#include "MaxCheck.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;

MaxCheck::MaxCheck()
{
}

MaxCheck::~MaxCheck()
{
}

void MaxCheck::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	global_log->debug() << "MaxCheck enabled" << std::endl;
}

void MaxCheck::readXML(XMLfileUnits& xmlconfig) {

	// targets
	uint32_t numTargets = 0;
	XMLfile::Query query = xmlconfig.query("targets/target");
	numTargets = query.card();
	global_log->info() << "MaxCheck: Number of component targets: " << numTargets << endl;
	if(numTargets < 1) {
		global_log->warning() << "MaxCheck: No target parameters specified. Program exit ..." << endl;
		Simulation::exit(-1);
	}
	string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator nodeIter;
	for( nodeIter = query.begin(); nodeIter; nodeIter++ ) {
		xmlconfig.changecurrentnode(nodeIter);

		uint32_t cid_ub;
		MaxVals mv;
		xmlconfig.getNodeValue("@cid", cid_ub);

		mv.F = mv.F2 = 0.;
		mv.v = mv.v2 = 0.;

		xmlconfig.getNodeValue("Fmax", mv.F);
		global_log->info() << "MaxCheck: Fmax(cid="<<cid_ub<<"): " << mv.F << endl;
		mv.F2 = mv.F * mv.F;

		xmlconfig.getNodeValue("vmax", mv.v);
		global_log->info() << "MaxCheck: vmax(cid="<<cid_ub<<"): " << mv.v << endl;
		mv.v2 = mv.v * mv.v;

		_maxVals[cid_ub] = mv;
	}  // loop over 'targets/target' nodes

}

void MaxCheck::afterForces(
        ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
        unsigned long simstep) {
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const ParticleIterator begin = particleContainer->iterator();

		uint64_t id;
		uint32_t cid_ub;
		double r[3];
		double F[3];
		double v[3];
		MaxVals absVals;

		for (ParticleIterator it = begin; it.hasNext(); it.next())
		{
			id=it->id();
			cid_ub=it->componentid()+1;
			for(uint8_t d=0; d<3; ++d) {
				r[d]=it->r(d);
				F[d]=it->F(d);
				v[d]=it->v(d);
			}

			// calc abs vals
			absVals.F2 = this->calcSquaredVectorLength(F);
			absVals.v2 = this->calcSquaredVectorLength(v);

			// mark for deletion
			MaxVals &mv = _maxVals[cid_ub];
			if(mv.F > 0. && absVals.F2 > mv.F2)
				it.deleteCurrentParticle();
//				_deletions.push_back(&(*it));

			if(mv.v > 0. && absVals.v2 > mv.v2)
				it.deleteCurrentParticle();
//				_deletions.push_back(&(*it));
		}

//		// perform deletions
//		for(auto it:_deletions) {
//			particleContainer->deleteMolecule(*it, true);
//			_deletions.pop_back();
//		}

	} // end pragma omp parallel
}


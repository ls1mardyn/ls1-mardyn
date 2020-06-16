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

MaxCheck::MaxCheck() {
}

MaxCheck::~MaxCheck() {
}

void MaxCheck::init(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {
	global_log->debug() << "MaxCheck enabled" << std::endl;
}

void MaxCheck::readXML(XMLfileUnits& xmlconfig) {

	// Timestep control
	_control.start = 0;
	_control.freq = 1;
	_control.stop = 10000;
	xmlconfig.getNodeValue("control/start", _control.start);
	xmlconfig.getNodeValue("control/frequency", _control.freq);
	xmlconfig.getNodeValue("control/stop", _control.stop);
	global_log->info() << "MaxCheck is acting start:freq:stop = "
			<< _control.start << ":" << _control.freq << ":" << _control.stop
			<< endl;

	// yrange
	Domain* domain = global_simulation->getDomain();
	_yrange.min = 0.;
	_yrange.max = domain->getGlobalLength(1);
	xmlconfig.getNodeValue("yrange/min", _yrange.min);
	xmlconfig.getNodeValue("yrange/max", _yrange.max);

	// targets
	uint32_t numTargets = 0;
	XMLfile::Query query = xmlconfig.query("targets/target");
	numTargets = query.card();
	global_log->info() << "MaxCheck: Number of component targets: "
			<< numTargets << endl;
	if (numTargets < 1) {
		global_log->warning()
				<< "MaxCheck: No target parameters specified. Program exit ..."
				<< endl;
		Simulation::exit(-1);
	}
	string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator nodeIter;
	for (nodeIter = query.begin(); nodeIter; nodeIter++) {
		xmlconfig.changecurrentnode(nodeIter);

		uint32_t cid_ub;
		MaxVals mv;
		xmlconfig.getNodeValue("@cid", cid_ub);

		mv.F = mv.F2 = 0.;
		mv.v = mv.v2 = 0.;
		mv.method = MCM_UNKNOWN;

		xmlconfig.getNodeValue("@method", mv.method);
		global_log->info() << "MaxCheck: Method(cid=" << cid_ub << "): "
				<< mv.method << endl;

		xmlconfig.getNodeValue("Fmax", mv.F);
		global_log->info() << "MaxCheck: Fmax(cid=" << cid_ub << "): " << mv.F
				<< endl;
		mv.F2 = mv.F * mv.F;

		xmlconfig.getNodeValue("vmax", mv.v);
		global_log->info() << "MaxCheck: vmax(cid=" << cid_ub << "): " << mv.v
				<< endl;
		mv.v2 = mv.v * mv.v;

		_maxVals[cid_ub] = mv;
	}  // loop over 'targets/target' nodes

}

void MaxCheck::siteWiseForces(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep)
{
	if (simstep < _control.start || simstep > _control.stop
			|| simstep % _control.freq != 0)
		return;
	this->checkMaxVals(particleContainer, domainDecomp, simstep);
}

void MaxCheck::afterForces(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep < _control.start || simstep > _control.stop
			|| simstep % _control.freq != 0)
		return;
	this->checkMaxVals(particleContainer, domainDecomp, simstep);
}

void MaxCheck::endStep(ParticleContainer *particleContainer,
		DomainDecompBase *domainDecomp, Domain *domain, unsigned long simstep
		) {
	if (simstep < _control.start || simstep > _control.stop
			|| simstep % _control.freq != 0)
		return;
	this->checkMaxVals(particleContainer, domainDecomp, simstep);
}

void MaxCheck::checkMaxVals(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep) {

#if defined(_OPENMP)
#pragma omp parallel
#endif
	{

		uint64_t id;
		uint32_t cid_ub;
		double r[3];
		double F[3];
		double v[3];
		MaxVals absVals;

		for (auto it = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
			id = it->getID();
			cid_ub = it->componentid() + 1;
			for (uint8_t d = 0; d < 3; ++d) {
				r[d] = it->r(d);
				F[d] = it->F(d);
				v[d] = it->v(d);
			}

			if(r[1] < _yrange.min || r[1] > _yrange.max)
				continue;

			// calc abs vals
			absVals.F2 = this->calcSquaredVectorLength(F);
			absVals.v2 = this->calcSquaredVectorLength(v);

			MaxVals &mv = _maxVals[cid_ub];

			if (MCM_LIMIT_TO_MAX_VALUE == mv.method) {
				if (mv.F > 0. && absVals.F2 > mv.F2) {
					double Fabs = sqrt(absVals.F2);
					double scale = mv.F / Fabs;
					it->scale_F(scale);
				}

				if (mv.v > 0. && absVals.v2 > mv.v2) {
					double vabs = sqrt(absVals.v2);
					double scale = mv.v / vabs;
					it->scale_v(scale);
				}
			} else if (MCM_LIMIT_TO_MAX_VALUE_OVERLAPS == mv.method) {
				if (mv.F > 0. && absVals.F2 > mv.F2) {
					double Fabs = sqrt(absVals.F2);
					double scale = mv.F / Fabs;
					it->scale_F(scale);

					if (mv.v > 0. && absVals.v2 > mv.v2) {
						double vabs = sqrt(absVals.v2);
						scale = mv.v / vabs;
						it->scale_v(scale);
					}
				}
			} else if (MCM_DELETE_PARTICLES == mv.method) {
				if ( (mv.F > 0. && absVals.F2 > mv.F2) || (mv.v > 0. && absVals.v2 > mv.v2) )
				    particleContainer->deleteMolecule(it, false);
				//				_deletions.push_back(&(*it));
			}
		}

		//		// perform deletions
		//		for(auto it:_deletions) {
		//			particleContainer->deleteMolecule(*it, true);
		//			_deletions.pop_back();
		//		}

	} // end pragma omp parallel
}

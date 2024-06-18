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
#include <array>


MaxCheck::MaxCheck() {
}

MaxCheck::~MaxCheck() {
}

void MaxCheck::init(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {
	Log::global_log->debug() << "MaxCheck enabled" << std::endl;
}

void MaxCheck::readXML(XMLfileUnits& xmlconfig) {

	// Timestep control
	_control.start = 0;
	_control.freq = 1;
	_control.stop = 10000;
	xmlconfig.getNodeValue("control/start", _control.start);
	xmlconfig.getNodeValue("control/frequency", _control.freq);
	xmlconfig.getNodeValue("control/stop", _control.stop);
	Log::global_log->info() << "[MaxCheck] Active period: start:freq:stop = "
			<< _control.start << ":" << _control.freq << ":" << _control.stop
			<< std::endl;

	// range
	Domain* domain = global_simulation->getDomain();
	_range.inclusive = true;
	_range.xmin = 0.;
	_range.xmax = domain->getGlobalLength(0);
	_range.ymin = 0.;
	_range.ymax = domain->getGlobalLength(1);
	_range.zmin = 0.;
	_range.zmax = domain->getGlobalLength(2);
	xmlconfig.getNodeValue("range/inclusive", _range.inclusive);
	xmlconfig.getNodeValue("range/xmin", _range.xmin);
	xmlconfig.getNodeValue("range/xmax", _range.xmax);
	xmlconfig.getNodeValue("range/ymin", _range.ymin);
	xmlconfig.getNodeValue("range/ymax", _range.ymax);
	xmlconfig.getNodeValue("range/zmin", _range.zmin);
	xmlconfig.getNodeValue("range/zmax", _range.zmax);

	// Warn if old config style is used as it has big influence on the simulation
	if (not xmlconfig.query("yrange").empty()) {
		Log::global_log->warning() << "[MaxCheck] <yrange> is deprecated! Use <range> and <ymin>/<ymax> instead." << std::endl;
	}

	Log::global_log->info() << "[MaxCheck] Apply " << ((_range.inclusive) ? "inside" : "outside") << " range:"
					   << " x = " << _range.xmin << " - " << _range.xmax << " ;"
					   << " y = " << _range.ymin << " - " << _range.ymax << " ;"
					   << " z = " << _range.zmin << " - " << _range.zmax << std::endl;

	// targets
	uint32_t numTargets = 0;
	XMLfile::Query query = xmlconfig.query("targets/target");
	numTargets = query.card();
	Log::global_log->info() << "[MaxCheck] Number of component targets: " << numTargets << std::endl;
	if (numTargets < 1) {
		Log::global_log->warning() << "[MaxCheck] No target parameters specified. Program exit ..." << std::endl;
		Simulation::exit(-1);
	}
	std::string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator nodeIter;
	for (nodeIter = query.begin(); nodeIter; nodeIter++) {
		xmlconfig.changecurrentnode(nodeIter);

		uint32_t cid_ub;
		MaxVals mv;
		xmlconfig.getNodeValue("@cid", cid_ub);

		// These values are later checked if >=0 --> do nothing if quantity not specified
		mv.F = mv.F2 = -1.;
		mv.v = mv.v2 = -1.;
		mv.M = mv.M2 = -1.;
		mv.L = mv.L2 = -1.;
		mv.method = MCM_UNKNOWN;

		xmlconfig.getNodeValue("@method", mv.method);
		Log::global_log->info() << "[MaxCheck] Method(cid=" << cid_ub << "): " << mv.method << std::endl;

		xmlconfig.getNodeValue("Fmax", mv.F);
		if (mv.F < 0.) {
			Log::global_log->info() << "[MaxCheck] Fmax(cid=" << cid_ub << ") not set" << std::endl;
		} else {
			Log::global_log->info() << "[MaxCheck] Fmax(cid=" << cid_ub << "): " << mv.F << std::endl;
			mv.F2 = mv.F * mv.F;
		}

		xmlconfig.getNodeValue("vmax", mv.v);
		if (mv.v < 0.) {
			Log::global_log->info() << "[MaxCheck] vmax(cid=" << cid_ub << ") not set" << std::endl;
		} else {
			Log::global_log->info() << "[MaxCheck] vmax(cid=" << cid_ub << "): " << mv.v << std::endl;
			mv.v2 = mv.v * mv.v;
		}

		xmlconfig.getNodeValue("Mmax", mv.M);
		if (mv.M < 0.) {
			Log::global_log->info() << "[MaxCheck] Mmax(cid=" << cid_ub << ") not set" << std::endl;
		} else {
			Log::global_log->info() << "[MaxCheck] Mmax(cid=" << cid_ub << "): " << mv.M << std::endl;
			mv.M2 = mv.M * mv.M;
		}

		xmlconfig.getNodeValue("Lmax", mv.L);
		if (mv.L < 0.) {
			Log::global_log->info() << "[MaxCheck] Lmax(cid=" << cid_ub << ") not set" << std::endl;
		} else {
			Log::global_log->info() << "[MaxCheck] Lmax(cid=" << cid_ub << "): " << mv.L << std::endl;
			mv.L2 = mv.L * mv.L;
		}

		_maxVals[cid_ub] = mv;
	}  // loop over 'targets/target' nodes

}

void MaxCheck::beforeEventNewTimestep(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep
	)
{
	if (simstep < _control.start || simstep > _control.stop || simstep % _control.freq != 0) {
		return;
	}
	this->checkMaxVals(particleContainer, domainDecomp, simstep);
}

void MaxCheck::afterForces(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep
)
{
	if (simstep < _control.start || simstep > _control.stop || simstep % _control.freq != 0) {
		return;
	}
	this->checkMaxVals(particleContainer, domainDecomp, simstep);
}

void MaxCheck::checkMaxVals(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep) {

#if defined(_OPENMP)
#pragma omp parallel
#endif
	{
		uint32_t cid_ub;
		std::array<double,3> r;
		std::array<double,3> F;
		std::array<double,3> v;
		std::array<double,3> M;
		std::array<double,3> L;
		MaxVals absVals;

		for (auto it = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
			cid_ub = it->componentid() + 1;
			for (uint8_t d = 0; d < 3; ++d) {
				r[d] = it->r(d);
				F[d] = it->F(d);
				v[d] = it->v(d);
				M[d] = it->M(d);
				L[d] = it->D(d);
			}

			// check range
			bool isInside = this->moleculeInsideRange(r);
			if (_range.inclusive and not isInside) {
				// If the range is inclusive, we skip if the particle is not in the range.
				continue;
			}
			if (not _range.inclusive and isInside) {
				// If the range is exclusive, we skip if the particle is in the range.
 				continue;
			}

			// calc abs vals
			absVals.F2 = this->calcSquaredVectorLength(F);
			absVals.v2 = this->calcSquaredVectorLength(v);
			absVals.M2 = this->calcSquaredVectorLength(M);
			absVals.L2 = this->calcSquaredVectorLength(L);

			MaxVals &mv = _maxVals[cid_ub];

			if (MCM_LIMIT_TO_MAX_VALUE == mv.method) {
				if (mv.F >= 0. && absVals.F2 > mv.F2) {
					double Fabs = sqrt(absVals.F2);
					double scale = mv.F / Fabs;
					it->scale_F(scale);
				}

				if (mv.v >= 0. && absVals.v2 > mv.v2) {
					double vabs = sqrt(absVals.v2);
					double scale = mv.v / vabs;
					it->scale_v(scale);
				}

				if (mv.M >= 0. && absVals.M2 > mv.M2) {
					double Mabs = sqrt(absVals.M2);
					double scale = mv.M / Mabs;
					it->scale_M(scale);
				}

				if (mv.L >= 0. && absVals.L2 > mv.L2) {
					double Labs = sqrt(absVals.L2);
					double scale = mv.L / Labs;
					it->scale_D(scale);
				}
			} else if (MCM_LIMIT_TO_MAX_VALUE_OVERLAPS == mv.method) {
				if (mv.F >= 0. && absVals.F2 > mv.F2) {
					double Fabs = sqrt(absVals.F2);
					double scale = mv.F / Fabs;
					it->scale_F(scale);

					if (mv.v >= 0. && absVals.v2 > mv.v2) {
						double vabs = sqrt(absVals.v2);
						scale = mv.v / vabs;
						it->scale_v(scale);
					}
				}
				if (mv.M >= 0. && absVals.M2 > mv.M2) {
					double Mabs = sqrt(absVals.M2);
					double scale = mv.M / Mabs;
					it->scale_M(scale);

					if (mv.L >= 0. && absVals.L2 > mv.L2) {
						double Labs = sqrt(absVals.L2);
						scale = mv.L / Labs;
						it->scale_D(scale);
					}
				}
			} else if (MCM_DELETE_PARTICLES == mv.method) {
				if ( (mv.F >= 0. && absVals.F2 > mv.F2) ||
					 (mv.v >= 0. && absVals.v2 > mv.v2) ||
					 (mv.M >= 0. && absVals.M2 > mv.M2) ||
					 (mv.L >= 0. && absVals.L2 > mv.L2) ) {
				    particleContainer->deleteMolecule(it, false);
				}
			}
		}
	} // end pragma omp parallel
}

bool MaxCheck::moleculeInsideRange(std::array<double,3>& r)
{
	return r[0] >= _range.xmin && r[0] < _range.xmax && r[1] >= _range.ymin && r[1] < _range.ymax && r[2] >= _range.zmin && r[2] < _range.zmax;
}

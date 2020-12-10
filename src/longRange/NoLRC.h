#ifndef NOLRC_H__
#define NOLRC_H__

#include "Domain.h"
#include "LongRangeCorrection.h"

#include "utils/Logger.h"
using Log::global_log;

class Simulation;
class Domain;
class LongRangeCorrection;

class NoLRC: public LongRangeCorrection{

public:
	NoLRC(double cutoffRadius, double cutoffRadiusLJ,  Domain* domain, Simulation* simulation) {
        global_log->info() << "No long range correction is used: UpotCorr = VirialCorr = 0" << std::endl;
        _domain = domain;
    };
	virtual ~NoLRC() {}

	virtual void init() {}
	virtual void readXML(XMLfileUnits& xmlconfig) {}
	virtual void calculateLongRange() {
        _domain->setUpotCorr(0.);
        _domain->setVirialCorr(0.);
      };
	virtual void writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {}

private:
    Domain* _domain;
};

#endif /* __NOLRC_H__ */

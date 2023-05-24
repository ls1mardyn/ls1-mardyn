#ifndef NOLRC_H__
#define NOLRC_H__

#include "Domain.h"
#include "LongRangeCorrection.h"

#include "utils/Logger.h"

class Simulation;
class Domain;
class LongRangeCorrection;

/**
 * Dummy long-range correction class which provides no correction.
 */
class NoLRC: public LongRangeCorrection{

public:
	NoLRC(double /* cutoffRadius */, double /* cutoffRadiusLJ */,  Domain* domain, Simulation* /* simulation */) {
        _domain = domain;
    };
	virtual ~NoLRC() {}

	virtual void init() { global_log->info() << "No long range correction is used: UpotCorr = VirialCorr = 0" << std::endl; }
	virtual void readXML(XMLfileUnits& /* xmlconfig */) {}
	virtual void calculateLongRange() {
        _domain->setUpotCorr(0.);
        _domain->setVirialCorr(0.);
      };
	virtual void writeProfiles(DomainDecompBase* /* domainDecomp */, Domain* /* domain */, unsigned long /* simstep */) {}

private:
    Domain* _domain;
};

#endif /* __NOLRC_H__ */

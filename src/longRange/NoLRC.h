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
	~NoLRC() override = default;

	virtual void init() { Log::global_log->info() << "No long range correction is used: UpotCorr = VirialCorr = 0" << std::endl; }
	virtual void readXML(XMLfileUnits& /* xmlconfig */) {}
	virtual void calculateLongRange() {
        _domain->setUpotCorr(0.);
        _domain->setVirialCorr(0.);
      };
	void writeProfiles(DomainDecompBase* /* domainDecomp */, Domain* /* domain */, unsigned long /* simstep */) {}

    // Get potential energy correction per molecule
    double getUpotCorr(Molecule* /* mol */) override { return 0.0; };

private:
    Domain* _domain;
};

#endif /* __NOLRC_H__ */

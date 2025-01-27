#include "io/EnergyLogWriter.h"

#include <ostream>

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/xmlfileUnits.h"


void EnergyLogWriter::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
	Log::global_log->info() << "Init global energy log." << std::endl;

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	// int numprocs = domainDecomp->getNumProcs();
	if (rank!= 0)
		return;
#endif

	std::stringstream outputstream;
	outputstream.write(reinterpret_cast<const char*>(&_writeFrequency), 8);

	std::ofstream fileout(_outputFilename, std::ios::out | std::ios::binary);
	fileout << outputstream.str();
	fileout.close();
}

void EnergyLogWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "Write frequency: " << _writeFrequency << std::endl;
	_outputFilename = "global_energy.log";
	xmlconfig.getNodeValue("outputfilename", _outputFilename);
	Log::global_log->info() << "Output filename: " << _outputFilename << std::endl;
}

void EnergyLogWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                              unsigned long simstep) {

	if( 0 != (simstep % _writeFrequency) ) {
		return;
	}

	unsigned long nNumMolsGlobalEnergyLocal = 0ul;
	double UkinLocal = 0.;
	double UkinTransLocal = 0.;
	double UkinRotLocal = 0.;
	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:nNumMolsGlobalEnergyLocal,UkinLocal,UkinTransLocal,UkinRotLocal)
	#endif
	{
		for (auto moleculeIter = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
			 moleculeIter.isValid(); ++moleculeIter) {
			nNumMolsGlobalEnergyLocal++;
			UkinLocal += moleculeIter->U_kin();
			UkinTransLocal += moleculeIter->U_trans();
			UkinRotLocal += moleculeIter->U_rot();
		}
	}

	// calculate global values
#ifdef ENABLE_PERSISTENT
	auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), nNumMolsGlobalEnergyLocal, UkinLocal, UkinTransLocal, UkinRotLocal);
	collComm.persistent();
	unsigned long nNumMolsGlobalEnergyGlobal;
	double UkinGlobal;
	double UkinTransGlobal;
	double UkinRotGlobal;
	collComm.get(nNumMolsGlobalEnergyGlobal, UkinGlobal, UkinTransGlobal, UkinRotGlobal);
#else
	domainDecomp->collCommInit(4);
	domainDecomp->collCommAppendUnsLong(nNumMolsGlobalEnergyLocal);
	domainDecomp->collCommAppendDouble(UkinLocal);
	domainDecomp->collCommAppendDouble(UkinTransLocal);
	domainDecomp->collCommAppendDouble(UkinRotLocal);
	domainDecomp->collCommAllreduceSum();
	unsigned long nNumMolsGlobalEnergyGlobal = domainDecomp->collCommGetUnsLong();
	double UkinGlobal = domainDecomp->collCommGetDouble();
	double UkinTransGlobal = domainDecomp->collCommGetDouble();
	double UkinRotGlobal = domainDecomp->collCommGetDouble();
	domainDecomp->collCommFinalize();
#endif

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	if (rank!= 0)
		return;
#endif

	const double globalUpot = domain->getGlobalUpot();
	const double globalT = domain->getGlobalCurrentTemperature();
	const double globalPressure = domain->getGlobalPressure() ;

	std::stringstream outputstream;
	outputstream.write(reinterpret_cast<const char*>(&nNumMolsGlobalEnergyGlobal), 8);
	outputstream.write(reinterpret_cast<const char*>(&globalUpot), 8);
	outputstream.write(reinterpret_cast<const char*>(&UkinGlobal), 8);
	outputstream.write(reinterpret_cast<const char*>(&UkinTransGlobal), 8);
	outputstream.write(reinterpret_cast<const char*>(&UkinRotGlobal), 8);
	outputstream.write(reinterpret_cast<const char*>(&globalT), 8);
	outputstream.write(reinterpret_cast<const char*>(&globalPressure), 8);

	std::ofstream fileout(_outputFilename, std::ios::app | std::ios::binary);
	fileout << outputstream.str();
	fileout.close();
}

void EnergyLogWriter::finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {}
